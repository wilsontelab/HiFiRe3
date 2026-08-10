//! Process reads with first alignments on a specific chromosome, provided as a 
//! message on a channel.

// imports
use std::error::Error;
use rustc_hash::{FxHashMap, FxHashSet};
use crossbeam::channel::{Receiver, Sender};
use minimap2::{Aligner as Minimap2};
use rust_htslib::bam::{Reader, Read, Record as BamRecord};
use mdi::pub_key_constants;
use mdi::workflow::Config;
use crate::snvs::*;

// constants
pub_key_constants!(
    // from environment variables
    MIN_FRAGMENT_READS
    MIN_HOMOZYGOUS_READS
    INDEX_FILE_PREFIX_WRK
    GENOME_REPEAT_MASKER_BED
    GENOME_SIMPLE_REPEAT_BED
);

// process chromosomes received on the channel
pub fn process_chrom(
    tool:     &SnvAnalysisTool,
    rx_chrom: Receiver<(String, u8)>,
    tx_data:  Sender<SnvChromWorkerData>,
) -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_usize_env(&[MIN_FRAGMENT_READS, MIN_HOMOZYGOUS_READS]);
    cfg.set_string_env(&[INDEX_FILE_PREFIX_WRK, GENOME_REPEAT_MASKER_BED, GENOME_SIMPLE_REPEAT_BED]);
    let chrom_file_prefix = cfg.get_string(INDEX_FILE_PREFIX_WRK); // created by split_bam_by_chrom
    let rmsk_simple_repeats_bed = cfg.get_string(GENOME_REPEAT_MASKER_BED);
    let trf_simple_repeats_bed = cfg.get_string(GENOME_SIMPLE_REPEAT_BED);

    // process chromosomes received on the channel
    for (chrom_name, chrom_index) in rx_chrom.iter() {
        let chrom_index_padded = format!("{:02}", chrom_index);
        let chrom_bam_path  = format!("{}.chr{}.bam", chrom_file_prefix, &chrom_index_padded);
        eprintln!("    {}", chrom_name);

        // open the input BAM file
        // all reads are on-target and have first alignment on chrom
        let mut chrom_bam = Reader::from_path(&chrom_bam_path)?;

        // assemble the chromosome worker tool
        let mut worker = SnvChromWorker {
            chrom: chrom_name.clone(),
            chrom_index,
            chrom_tid: tool.fa.fai().tid(&chrom_name).expect("Failed to get chrom TID"),
            simple_repeats: SimpleRepeats::new(
                tool, &chrom_name, rmsk_simple_repeats_bed, trf_simple_repeats_bed
            ),
            min_fragment_reads:   *cfg.get_usize(MIN_FRAGMENT_READS),
            min_homozygous_reads: *cfg.get_usize(MIN_HOMOZYGOUS_READS),
            minimap2:       Minimap2::builder().map_hifi().with_cigar(),
            frag_vars:      FragmentVariants::new(),
            encoding:       AlignmentEncoding::new(), // read encoding for visualization
            tracking_variants: Vec::with_capacity(1024),
            seq0_bases:     Vec::with_capacity(MAX_EXPECTED_READ_LEN), // used with cs_map for consensus calling
            cs_map:         Vec::with_capacity(MAX_EXPECTED_READ_LEN),
            hap_vs_ref:     Vec::with_capacity(256), // consensus encoding for visualization
            hap_vars:       FxHashMap::default(),
            hap_votes:      FxHashMap::default(),
            var_tgt_pos0:   None,
            tgt_bases:      String::with_capacity(128),
            alt_bases:      String::with_capacity(128),
            alt_qual:       Vec::with_capacity(128),
            allowed:        true,
            cs_op:          ':',
            op_val:         String::with_capacity(128),
            variant_tally:       VariantsTally::new(),
            variant_reads_tally: VariantReadsTally::new(),
            // debug: ReFragment { start0: 112213159, end1: 112220205 },
            // show_debug: false,
        };
        worker.hap_vars.insert(Haplotype::Haplotype1, FxHashSet::default());
        worker.hap_vars.insert(Haplotype::Haplotype2, FxHashSet::default());
        worker.hap_vars.insert(Haplotype::Homozygous, FxHashSet::default());
        worker.hap_votes.insert(Haplotype::Haplotype1, 0);
        worker.hap_votes.insert(Haplotype::Haplotype2, 0);

        // process alignment records one at a time, add to growing RE fragment collections
        let mut aln = BamRecord::new();
        let mut chrom_aln_count:      usize = 0;
        let mut chrom_aln_count_used: usize = 0;
        let mut fragment_reads = FragmentReads::new();
        while let Some(result) = chrom_bam.read(&mut aln) {
            match result {
                Ok(_)  => {
                    chrom_aln_count += 1;
                    chrom_aln_count_used += fragment_reads.insert(&aln);
                },
                Err(_) => panic!("BAM parsing failed")
            }
        }

        // post-process read groups by re-aligning reads to fragment consensus(es)
        let mut haplotype_consensuses = HaplotypeConsensuses::new();
        let mut reads_on_reference = FragmentHaplotypes::new();
        let mut reads_on_haplotype = FragmentHaplotypes::new();
        worker.analyze_reads(
            tool, fragment_reads, 
            &mut haplotype_consensuses,
            &mut reads_on_reference,
            &mut reads_on_haplotype,
        );

        // finish processing and writing pileup and variants
        let variants_file_path = format!(
            "{}.chr{}.snv_indel.variants.txt.bgz", 
            chrom_file_prefix, &chrom_index_padded
        );
        let variant_reads_file_path = format!(
            "{}.chr{}.snv_indel.variant_reads.txt.bgz", 
            chrom_file_prefix, &chrom_index_padded
        );
        let reads_on_reference_path = format!(
            "{}.chr{}.fragments.on_reference.bed.bgz", 
            chrom_file_prefix, &chrom_index_padded
        );
        let reads_on_haplotype_path = format!(
            "{}.chr{}.fragments.on_haplotype.bed.bgz", 
            chrom_file_prefix, &chrom_index_padded
        );
        let variant_metadata = VariantsTally::write_sorted(
            tool, &mut worker, &mut haplotype_consensuses,
            variants_file_path
        );
        let variant_reads_metadata = VariantReadsTally::write_sorted(
            tool, &mut worker,
            variant_reads_file_path
        );
        let reads_on_reference_metadata = FragmentHaplotypes::write_sorted(
            tool, &mut worker, &mut haplotype_consensuses,
            &reads_on_reference, 
            reads_on_reference_path
        );
        let reads_on_haplotype_metadata = FragmentHaplotypes::write_sorted(
            tool, &mut worker, &mut haplotype_consensuses,
            &reads_on_haplotype, 
            reads_on_haplotype_path
        );

        // send error corrected metadata to main thread
        tx_data.send(SnvChromWorkerData::TotalAlnCount(chrom_aln_count))?;
        tx_data.send(SnvChromWorkerData::UsableAlnCount((chrom_name.clone(), chrom_aln_count_used)))?;
        tx_data.send(SnvChromWorkerData::VariantMetadata(variant_metadata))?;
        tx_data.send(SnvChromWorkerData::VariantReadsMetadata(variant_reads_metadata))?;
        tx_data.send(SnvChromWorkerData::ReadsOnReferenceMetadata(reads_on_reference_metadata))?;
        tx_data.send(SnvChromWorkerData::ReadsOnHaplotypeMetadata(reads_on_haplotype_metadata))?;
    }
    Ok(())
}
