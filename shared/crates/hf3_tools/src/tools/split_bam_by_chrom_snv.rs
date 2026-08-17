//! Split input name-sorted BAM file(s) into temporary per-chromosome BAM files 
//! based on the chromosome of each alignment. Unlike for SVs, only reads with
//! single ~end-to-end alignments are printed and used for SNV calling.
//! 
//! Output only includes usable on-target duplex PacBioStrand reads, as defined 
//! by the presence of the minimap2 `cs:Z:` and HiFiRe3 `dd:Z:` tags, with the 
//! requested minimum number of sequencing passes per insert to ensure duplex
//! basecalling accuracy.
//! 
//! Along the way, also extract the two source strands from the original,
//! unmerged PacBio unaligned bams basecalled --by-strand. QNAMEs kept from the
//! first pass through the aligned duplex-read bam are used to filter and assign  
//! unaligned strand-level reads to their appropriate chromosome file.
//! 
//! Support multiple input BAM files for multi-sample variant calling.

// imports
use std::error::Error;
use std::str::from_utf8_unchecked;
use rustc_hash::FxHashMap;
use rust_htslib::bam::{
    Reader, Read, Writer, 
    Record as BamRecord, Header, HeaderView, 
    Format, record::Aux
};
use rust_htslib::tpool::ThreadPool;
use mdi::pub_key_constants;
use mdi::workflow::{Workflow, Config, Counters};
use mdi::OutputFile;
use genomex::genome::{Chroms, TargetRegions};
use genomex::bam::tags as bam_tags;
use genomex::bam::qual::median_qual_aln;
use genomex::bam::cigar::{get_clip_left, get_clip_right};
use crate::formats::hf3_tags::*;
use crate::snvs::check_pacbio_strand;

// constants
const TOOL: &str = "split_by_chrom_snv";
pub_key_constants!(
    // from environment variables
    N_CPU
    IS_COMPOSITE_GENOME
    MIN_AVG_BASE_QUAL
    MIN_N_PASSES
    NAME_BAM_FILES // harvested upstream from the input directories
    HIFI_BAM_FILES // set by user to the same files as for `basecall PacBio`
    FAIL_BAM_FILES
    INDEX_FILE_PREFIX_WRK
    SNV_SAMPLES_FILE
    // counter keys
    N_ALNS
    N_USABLE_ALNS
    N_BASES
    N_ALNS_BY_GENOME
    N_BASES_BY_GENOME
    N_ALNS_BY_SAMPLE
    N_BASES_BY_SAMPLE
);
const ML_TAG: &[u8] = BASE_MODS.as_bytes();
const MM_TAG: &[u8] = BASE_MOD_PROBS.as_bytes();
const IP_TAG: &[u8] = INTER_PULSE_DURATION.as_bytes();
const PW_TAG: &[u8] = PULSE_WIDTH.as_bytes();
const CS_TAG: &[u8] = DIFFERENCE_STRING.as_bytes();
const DD_TAG: &[u8] = STRAND_DIFFERENCES.as_bytes();
const SB_TAG: &[u8] = SAMPLE_BIT.as_bytes(); // set here
const MIN_MAPQ: u8  = 50; // TODO: expose as options?
const MAX_CLIP: u32 = 25;

// main function called by xxx_tools main()
pub fn main() -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_u32_env(&[N_CPU]);
    cfg.set_bool_env(&[IS_COMPOSITE_GENOME]);
    cfg.set_u8_env(&[MIN_AVG_BASE_QUAL]);
    cfg.set_f64_env(&[MIN_N_PASSES]);
    cfg.set_string_env(&[
        NAME_BAM_FILES, HIFI_BAM_FILES, FAIL_BAM_FILES,
        INDEX_FILE_PREFIX_WRK, SNV_SAMPLES_FILE
    ]);
    let min_avg_base_qual = *cfg.get_u8(MIN_AVG_BASE_QUAL);
    let min_n_passes      = *cfg.get_f64(MIN_N_PASSES) as f32;

    // validate we are working with the expected read data type
    check_pacbio_strand(TOOL, &mut cfg)?;

    // initialize counters
    let mut ctrs = Counters::new(TOOL, &[
        (N_ALNS,        "alignments processed"),
        (N_USABLE_ALNS, "usable on-target alignments in output"),
        (N_BASES,       "reference bases in usable on-target alignments"),
    ]);
    ctrs.add_keyed_counters(&[
        (N_ALNS_BY_GENOME,  "on-target alignments by genome"),
        (N_BASES_BY_GENOME, "reference bases in on-target alignments by genome"),
        (N_ALNS_BY_SAMPLE,  "on-target alignments by sample"),
        (N_BASES_BY_SAMPLE, "reference bases in on-target alignments by sample"),
    ]);

    // initialize the tool
    let mut w = Workflow::new(TOOL, cfg, ctrs);
    w.log.initializing();

    // collect the working chromosomes
    let chroms = Chroms::new(&mut w.cfg);
    let targets = TargetRegions::from_env(&mut w, false);
    let on_target_chroms = targets.get_region_chroms(&chroms);
    let is_composite = *w.cfg.get_bool(IS_COMPOSITE_GENOME);
    let chrom_file_prefix = w.cfg.get_string(INDEX_FILE_PREFIX_WRK);

    // use a thread pool for BAM reading and writing
    let tpool = ThreadPool::new(w.cfg.get_u32(N_CPU) - 1)?;

    // initialize the output BAM writers, one per target chromosome shared over all samples
    let name_bam_files = w.cfg.get_string(NAME_BAM_FILES);
    let name_bam_paths = name_bam_files.split(',').collect::<Vec<&str>>();
    let name_bam = Reader::from_path(name_bam_paths[0])?;
    let header_view = name_bam.header().clone(); // for TID lookups
    let header = Header::from_template(&header_view); // shared header for each output BAM writer
    let mut duplex_writers: FxHashMap<u32, Writer> = FxHashMap::default(); // TID -> file writer named by our padded chrom index
    let mut strand_writers: FxHashMap<u32, Writer> = FxHashMap::default();
    for (chrom, chrom_index) in on_target_chroms.iter() {
        let tid = header_view
            .tid(chrom.as_bytes())
            .expect(format!("{} not found in BAM header", chrom).as_str()) as u32;
        let chrom_index_padded = format!("{:02}", chrom_index);
        let mut duplex_writer = Writer::from_path(
            format!("{}.chr{}.duplex.bam", chrom_file_prefix, chrom_index_padded),
            &header,
            Format::Bam
        ).expect(&format!("Failed to create duplex BAM writer for chrom {}", chrom));
        duplex_writer.set_thread_pool(&tpool)?;
        duplex_writers.insert(tid, duplex_writer);
        let mut strand_writer = Writer::from_path(
            format!("{}.chr{}.strand.bam", chrom_file_prefix, chrom_index_padded),
            &header,
            Format::Bam
        ).expect(&format!("Failed to create strand BAM writer for chrom {}", chrom));
        strand_writer.set_thread_pool(&tpool)?;
        strand_writers.insert(tid, strand_writer);
    }
    drop(name_bam);

    // run through multiple duplex name BAM files to support single and multi-sample analyses
    w.log.print("streaming duplex BAM records");
    let samples_file = w.cfg.get_string(SNV_SAMPLES_FILE);
    let header = vec!["sample_bit", "sample_name"];
    let mut samples_file = OutputFile::open_file(&samples_file, b'\t', Some(&header)); 
    let mut sample_bit: u32 = 1;
    let mut qname_to_chrom: FxHashMap<String, u32> = FxHashMap::default();
    for name_bam_path in name_bam_paths {
        let bam_file_name = name_bam_path.split('/').last().unwrap();
        let sample_name = bam_file_name.split('.').nth(0).unwrap();
        samples_file.write_record(vec![&sample_bit.to_string(), sample_name]);
        let mut name_bam = Reader::from_path(name_bam_path)?;
        name_bam.set_thread_pool(&tpool)?;

        // process input BAM records
        eprintln!("    {}", sample_name);
        let mut aln = BamRecord::new();
        while let Some(result) = name_bam.read(&mut aln) {
            match result {
                Ok(_)  => print_duplex_aln(
                    &mut aln, 
                    &header_view, 
                    is_composite, 
                    min_avg_base_qual,
                    min_n_passes,
                    &mut duplex_writers, 
                    &mut w.ctrs, 
                    sample_bit, 
                    sample_name,
                    &mut qname_to_chrom,
                )?,
                Err(_) => panic!("duplex BAM parsing failed")
            }
        }
        sample_bit <<= 1; 
    }   

    // use the collected QNAMEs to also store the original--by-strand reads
    // for all kept duplex reads for use in final error correction
    w.log.print("streaming by-strand BAM records");
    print_strand_reads(
        &w.cfg, HIFI_BAM_FILES, &tpool, &qname_to_chrom, &mut strand_writers
    )?;
    print_strand_reads(
        &w.cfg, FAIL_BAM_FILES, &tpool, &qname_to_chrom, &mut strand_writers
    )?;

    // report counter values
    w.ctrs.print_grouped(&[
        &[N_ALNS, N_USABLE_ALNS, N_BASES],
        &[N_ALNS_BY_GENOME],
        &[N_BASES_BY_GENOME],
        &[N_ALNS_BY_SAMPLE],
        &[N_BASES_BY_SAMPLE],
    ]);
    Ok(())
}

// print and count on-target alignments
fn print_duplex_aln(
    aln:          &mut BamRecord, 
    header_view:  &HeaderView,
    is_composite: bool,
    min_avg_base_qual: u8,
    min_n_passes:      f32,
    duplex_writers: &mut FxHashMap<u32, Writer>,
    ctrs:           &mut Counters,
    sample_bit:     u32,
    sample_name:    &str,
    qname_to_chrom: &mut FxHashMap<String, u32>,
) -> Result<(), Box<dyn Error>> {
    ctrs.increment(N_ALNS);

    // filters below are roughly ordered by likelihood to hit and computational cost

    // skip low confidence alignments, including unmapped reads
    if aln.mapq() < MIN_MAPQ || 

    // skip non-duplex reads used for SV calling but not SNV calling  
        !aln.aux(DD_TAG).is_ok() ||

    // require a minimum number of PacBio passes in duplex reads
        bam_tags::get_tag_f32_default(aln, PACBIO_EFF_COVERAGE, 0.0) < min_n_passes || 
    
    // rare mapped reads lack a minimap2 cs tag for unknown reasons
        !aln.aux(CS_TAG).is_ok() ||

    // require single ~end-to-end alignment on reads, which excludes:
    //    - SV-containing reads
    //    - reads with overly large end clips relative to their outermost RE sites
        get_clip_left( aln) > MAX_CLIP ||
        get_clip_right(aln) > MAX_CLIP ||

    // skip reads with low average base quality 
    // TODO: skip reads with too many no-call N bases? 
        median_qual_aln(aln) < min_avg_base_qual
    {
        return Ok(()); 
    }

    // skip reads in untargeted samples that map to other than nuclear chromosomes
    let tid = aln.tid() as u32;
    if let Some(writer) = duplex_writers.get_mut(&tid){
        ctrs.increment(N_USABLE_ALNS);
        ctrs.increment_keyed(N_ALNS_BY_SAMPLE, sample_name);

        // add the sample bit tag for multi-sample comparison
        aln.push_aux(SB_TAG, Aux::U32(sample_bit)).unwrap();

        // commit on-target reads to temporary BAM files
        writer.write(aln)?; 

        // store committed QNAMEs for printing to --by-strand files
        let duplex_qname = unsafe { from_utf8_unchecked(aln.qname()) };
        qname_to_chrom.insert(duplex_qname.to_string(), tid);

        // increment counters
        let cigar_view = aln.cigar();
        let n_bases = cigar_view.end_pos() as usize - aln.pos() as usize;
        ctrs.add_to(N_BASES, n_bases); 
        ctrs.add_to_keyed(N_BASES_BY_SAMPLE, sample_name, n_bases);
        if is_composite {
            let chrom = unsafe{ from_utf8_unchecked(header_view.tid2name(tid))}; 
            let genome_name = chrom.split_once('_').unwrap().1; // e.g., chr1_hs1
            ctrs.increment_keyed(N_ALNS_BY_GENOME, genome_name);
            ctrs.add_to_keyed(N_BASES_BY_GENOME,  genome_name, n_bases);
        }                
    }
    Ok(())
}

// print --by-strand reads for the same printed duplex reads
fn print_strand_reads(
    cfg:       &Config,
    bam_files: &str,
    tpool:     &ThreadPool,
    qname_to_chrom: &FxHashMap<String, u32>,
    strand_writers: &mut FxHashMap<u32, Writer>,
) -> Result<(), Box<dyn Error>> {
    let bam_files = cfg.get_string(bam_files);
    let bam_paths = bam_files.split(',').collect::<Vec<&str>>();
    let mut aln = BamRecord::new();
    for bam_path in bam_paths {
        let bam_file_name = bam_path.split('/').last().unwrap();
        eprintln!("    {}", bam_file_name);
        let mut bam = Reader::from_path(bam_path)?;
        bam.set_thread_pool(tpool)?;
        while let Some(aln_result) = bam.read(&mut aln) {
            match aln_result {
                Ok(_)  => {
                    // m21026_251212_225317/90968331/ccs/fwd   4       *       0       255
                    // m21026_251212_225317/94245975/ccs/rev   4       *       0       255
                    // m21026_251212_225317/105452886/ccs      16      chr3_hs1        124859916       60
                    let strand_qname = unsafe { from_utf8_unchecked(aln.qname()) };
                    let parts: Vec<&str> = strand_qname.split('/').collect();
                    if parts.len() < 2 { continue }
                    let duplex_qname = format!("{}/{}/ccs", parts[0], parts[1]);   
                    if let Some(tid) = qname_to_chrom.get(&duplex_qname){
                        let _ = aln.remove_aux(ML_TAG); // debulk for caching
                        let _ = aln.remove_aux(MM_TAG);
                        let _ = aln.remove_aux(IP_TAG);
                        let _ = aln.remove_aux(PW_TAG);
                        strand_writers.get_mut(&tid).unwrap().write(&aln)?; 
                    }
                },
                Err(_) => panic!("strand BAM parsing failed")
            }
        }
    }
    Ok(())

}
