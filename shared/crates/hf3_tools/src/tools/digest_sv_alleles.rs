//! Annotate SVs in a Sniffles2 or CuteSV VCF with the size of the RE fragments 
//! that contain SV junctions for a requested series of REs, as well as the 
//! distance from the junction to the nearest RE site. Steps are to:
//!     - scan SV records from the VCF
//!     - reject Sniffles2 symbolic INV records (INV junctions are in BND records)
//!     - reject SVs that are too close to chrom ends
//!     - reject SVs that are too small
//!     - collect all unique QNAMEs that supported any passing junction
//!     - collect SEQ for those reads from the BAM file that generated the VCF
//!     - for each junction, assemble its reference sequence, i.e., scaffold, from VCF data
//!     - use mini_consensus to assemble a read consensus bounding the junction, i.e., a haplotype
//!     - locate the junction in the consensus
//!     - in silico digest the consensus with each enzyme and report:
//!         - fragment sizes in INFO/<RE_NAME>SIZE
//!         - nearest RE site distances in INFO/<RE_NAME>DIST

// imports
use std::error::Error;
use std::str::from_utf8_unchecked;
use crossbeam::channel::{bounded, Receiver, Sender};
use rustc_hash::FxHashMap;
use rust_htslib::bcf;
use rust_htslib::bcf::Read as VcfRead;
use rust_htslib::bam;
use rust_htslib::bam::Read as BamRead;
use rust_htslib::faidx;
use minimap2::{Aligner as Minimap2, Mapping};
use mdi::pub_key_constants;
use mdi::workflow::{Workflow, Config, Counters};
use mini_consensus::{Resolver, ResolverConfig};
use crate::benchmark::*;

// constants for environment variable, config, and counter keys, etc.
const TOOL: &str = "digest_sv_alleles";
pub_key_constants!(
    // from environment variables
    N_CPU
    GENOME_FASTA
    BENCHMARK_BAM
    TRUTHSET_VCF
    TRUTHSET_ALGORITHM
    RESTRICTION_ENZYMES
    MIN_SV_LEN
    // counter keys
    N_INPUT_SV_RECORDS
    N_OUTPUT_JXN_RECORDS
    N_UNIQUE_QNAMES
);
const CHANNEL_CAPACITY: usize = 100;

/// Tool for processing VCF SV calls.
struct DigestionTool {
    is_sniffles: bool,
    n_cpu: usize,
    min_sv_len: u32,
    jxns: Vec<bcf::Record>,
    scfs: Vec<Vec<u8>>,
    seqs: FxHashMap<Vec<u8>, Vec<u8>>,
    res:  Vec<ReSite>,
}

// function called by hf3_tools main()
pub fn main() -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_string_env(&[
        GENOME_FASTA, BENCHMARK_BAM, 
        TRUTHSET_VCF, TRUTHSET_ALGORITHM,
        RESTRICTION_ENZYMES
    ]);
    cfg.set_usize_env( &[N_CPU, MIN_SV_LEN]);

    // initialize counters
    let ctrs = Counters::new(TOOL, &[
        (N_INPUT_SV_RECORDS,   "number of input SV records in VCF"),
        (N_OUTPUT_JXN_RECORDS, "number of output junction records after size and other limited filtering"),
        (N_UNIQUE_QNAMES,      "number of unique reads used to assemble junction allele consensuses"),
    ]);

    // initialize the tool
    let mut w = Workflow::new(TOOL, cfg, ctrs);
    w.log.initializing();

    // create the digestion tool
    let mut tool = DigestionTool {
        is_sniffles: w.cfg.get_string(TRUTHSET_ALGORITHM) == "sniffles",
        n_cpu:      *w.cfg.get_usize(N_CPU),
        min_sv_len: *w.cfg.get_usize(MIN_SV_LEN) as u32,
        jxns: Vec::with_capacity(10_000),
        scfs: Vec::with_capacity(10_000),
        seqs: FxHashMap::default(),
        res: w.cfg.get_string(RESTRICTION_ENZYMES).split(',').map(|re: &str|{
            let re: Vec<_> = re.split(':').collect();
            ReSite::from_re_name(re[0], Some(re[1]))
        }).collect(),
    };
    
    // execute the workflow
    let vcf_out = scan_vcf(&mut w, &mut tool)?;
    collect_read_seqs(&mut w, &mut tool)?;
    resolve_consensuses(&mut w, &mut tool, vcf_out)?;

    // print counts
    w.ctrs.print_grouped(&[
        &[N_INPUT_SV_RECORDS, N_OUTPUT_JXN_RECORDS, N_UNIQUE_QNAMES],
    ]);
    Ok(())
}

// -----------------------------------------------------------------------------
// initial VCF scan
// -----------------------------------------------------------------------------
/// Scan the input VCF file to collect SV calls and unique QNAMEs.
fn scan_vcf (w: &mut Workflow, tool: &mut DigestionTool) -> Result<bcf::Writer, Box<dyn Error>> {
    w.log.print("scanning input VCF");

    // initialize all readers
    let truthset_vcf = w.cfg.get_string(TRUTHSET_VCF).to_string();
    let mut vcf = bcf::Reader::from_path(truthset_vcf)?;
    let header = vcf.header().to_owned();
    let genome_fasta = w.cfg.get_string(GENOME_FASTA).to_string();
    let faidx = faidx::Reader::from_path(genome_fasta)?;
    let mut vcf_out = get_bcf_writer(&tool, bcf::Header::from_template(&header))?;

    // process SV records one at a time
    // ignore any non-SV variants, although they are not expected from benchmark VCF files
    for result in vcf.records() {
        let mut jxn = result?;
        let Some(sv_type) = get_vcf_sv_type(&jxn)? else { continue };
        let sv_type_str = unsafe{ from_utf8_unchecked(&sv_type) };
        w.ctrs.increment(N_INPUT_SV_RECORDS);

        // short-circuit Sniffles2 symbolic inversions, they aren't unique junctions 
        // (see notes in build_scaffold())
        if tool.is_sniffles && sv_type_str == "INV" { continue }
        let alt = get_vcf_alt(&jxn)?;

        // collect metadata on the left side of the SV; check position on chrom
        let (chrom_1, pos1_1) = get_vcf_left_node(&jxn)?;
        let chrom_len_1 = faidx.fetch_seq_len(&chrom_1) as u32;
        if pos1_1 <= SCAFFOLD_FLANK_LEN || pos1_1 + SCAFFOLD_FLANK_LEN >= chrom_len_1 { 
            continue
        }

        // collect metadata on the right side of the SV; check position on chrom
        let (chrom_2, pos1_2) = get_vcf_right_node(&jxn, &chrom_1, alt, &sv_type)?;
        let chrom_len_2 = faidx.fetch_seq_len(&chrom_2) as u32;
        if pos1_2 <= SCAFFOLD_FLANK_LEN || pos1_2 + SCAFFOLD_FLANK_LEN >= chrom_len_2 { 
            continue
        }

        // skip SVs that are too small
        let sv_len = if chrom_1 == chrom_2 {
            (pos1_2 as i32 - pos1_1 as i32).unsigned_abs()
        } else {
            999_999_999
        };
        let min_sv_len: u32 = match sv_type_str {
            "INV" | "BND" => SCAFFOLD_FLANK_LEN as u32, // ensure that scaffolds don't cross both INV junctions
            _ => tool.min_sv_len
        };
        if sv_len < min_sv_len { continue } 

        // collect information required to collect SEQs and build consensus for site analysis
        let Some(qnames) = get_vcf_qnames(&jxn)? else { continue }; 
        let downsampled_qnames = &qnames[0..qnames.len().min(MAX_CONSENSUS_SEQS)];
        for qname in downsampled_qnames {
            tool.seqs.insert(qname.clone(), Vec::new());
        } 

        // build the scaffold sequence for later use
        //  - throughout, truthset VCF does not provide information on bases inserted 
        //    at junctions, but these bases should resolve in the consensus
        //  - Sniffles2 INFO/STRAND is NOT information about the junction path, it  
        //    summarizes the orientation(s) of supporting reads    
        let mut scf: Vec<u8>;   
        let scf2: Vec<u8>;
        match sv_type_str {
            "DEL" => {
                // POS: 2707, REF: TCTgca...ccg, ALT: T, SVLEN=-259; END=2966; (2707 + 259 = 2966)
                // DEL is a single junction that is a fully resolved SV (sub)path
                scf  = get_ref_seq(&faidx, &chrom_1, pos1_1 - SCAFFOLD_FLANK_LEN + 1, pos1_1, false)?;
                scf2 = get_ref_seq(&faidx, &chrom_2, pos1_2 + 1, pos1_2 + SCAFFOLD_FLANK_LEN, false)?;  
            },
            "INS" => {
                // POS: 2566, REF: c, ALT: cCCTCACCCTCACCCTCACCCTCACCCTGACCCTGACCCTGACCCTCACCCTGACCCTCACCCTCACCCTCA, SVLEN=71; END=2566;
                // INS is a single junction that is a fully resolved SV (sub)path (but DEL+INS at same junction is ambiguous)
                scf  = get_ref_seq(&faidx, &chrom_1, pos1_1 - SCAFFOLD_FLANK_LEN + 1, pos1_1, false)?;
                scf2 = get_ref_seq(&faidx, &chrom_2, pos1_2 + 1, pos1_2 + SCAFFOLD_FLANK_LEN, false)?;
                scf.extend(alt[1..].to_ascii_uppercase());
            },
            "DUP" => {
                // POS: 3009031, REF: G, ALT: <DUP>, SVLEN=1394; END=3010425; (3009031 + 1394 = 3010425)
                // DUP is a single junction that is a fully resolved SV (sub)path
                scf  = get_ref_seq(&faidx, &chrom_2, pos1_2 - SCAFFOLD_FLANK_LEN + 1, pos1_2, false)?;
                scf2 = get_ref_seq(&faidx, &chrom_1, pos1_1 + 1, pos1_1 + SCAFFOLD_FLANK_LEN, false)?;
            },
            "INV" => {
                // INV encoding
                //  Sniffles2 
                //    A symbolic representation of the overall inversion, ignored. 
                //      chr10  91443587  Sniffles2.INV.B69S9  N  <INV>              SVTYPE=INV;SVLEN=2442;END=91446029;SUPPORT=45;
                //    Two matching BND records for the actual junctions; these are what we care about.
                //      chr10  91443587  Sniffles2.BND.CF3S9  N  N]chr10:91446029]  SVTYPE=BND;SUPPORT=25;
                //      chr10  91443588  Sniffles2.BND.CF6S9  N  [chr10:91446031[N  SVTYPE=BND;SUPPORT=22;
                //    A mysterious entry that seems to encode a direct continuity from base 91443587 to the next, should be ignored.
                //      chr10  91443587  Sniffles2.BND.CF2S9  N  N[chr10:91443588[  SVTYPE=BND;SUPPORT=20;
                // CuteSV
                //    Two INV records encoding the paired ++ and -- junctions, used as is.
                //      chr10  91443587  cuteSV.INV.41        c  <INV>              SVTYPE=INV;SVLEN=2442;END=91446029;RE=47;STRAND=++;
                //      chr10  91443588  cuteSV.INV.42        c  <INV>              SVTYPE=INV;SVLEN=2442;END=91446030;RE=46;STRAND=--;
                // in all, (POS + LEN = END)
                let strand = get_vcf_strand_bytes(
                    &jxn, 
                    &SvSource::CuteSV,
                    b"INV",
                    alt,
                )?;
                let strand_str = unsafe{ from_utf8_unchecked(&strand) };
                match strand_str {
                    "++" => {
                        scf  = get_ref_seq(&faidx, &chrom_1, pos1_1 - SCAFFOLD_FLANK_LEN + 1, pos1_1, false)?;
                        scf2 = get_ref_seq(&faidx, &chrom_2, pos1_2 - SCAFFOLD_FLANK_LEN + 1, pos1_2, true)?;
                    },
                    _ => {
                        scf  = get_ref_seq(&faidx, &chrom_1, pos1_1 + 1, pos1_1 + SCAFFOLD_FLANK_LEN, true)?;
                        scf2 = get_ref_seq(&faidx, &chrom_2, pos1_2 + 1, pos1_2 + SCAFFOLD_FLANK_LEN, false)?;
                    }
                }
            },
            "BND" => {
                // CHR1: chr5, POS1: 16396441, REF: C, ALT: [chr1:145169227[C, CHR2=chr1; (no SVLEN or END)
                // orientation encoded in ALT
                // BND is a single junction from either a fully or partially resolved SV (sub)path
                let bnd_type = BndType::from_alt(alt)?;
                match bnd_type {
                    BndType::AfterRight => {
                        scf  = get_ref_seq(&faidx,  &chrom_1, pos1_1 - SCAFFOLD_FLANK_LEN + 1, pos1_1, false)?;
                        scf2 = get_ref_seq(&faidx, &chrom_2, pos1_2, pos1_2 + SCAFFOLD_FLANK_LEN - 1, false)?;
                    },
                    BndType::AfterRcLeft => {
                        scf  = get_ref_seq(&faidx,  &chrom_1, pos1_1 - SCAFFOLD_FLANK_LEN + 1, pos1_1, false)?;
                        scf2 = get_ref_seq(&faidx, &chrom_2, pos1_2 - SCAFFOLD_FLANK_LEN + 1, pos1_2, true)?;
                    },
                    BndType::BeforeLeft => {
                        scf  = get_ref_seq(&faidx,  &chrom_2, pos1_2 - SCAFFOLD_FLANK_LEN + 1, pos1_2, false)?;
                        scf2 = get_ref_seq(&faidx, &chrom_1, pos1_1, pos1_1 + SCAFFOLD_FLANK_LEN - 1, false)?;
                    },
                    BndType::BeforeRcRight => {
                        scf  = get_ref_seq(&faidx,  &chrom_2, pos1_2, pos1_2 + SCAFFOLD_FLANK_LEN - 1, true)?;
                        scf2 = get_ref_seq(&faidx, &chrom_1, pos1_1, pos1_1 + SCAFFOLD_FLANK_LEN - 1, false)?;
                    },
                }
            },
            _ => panic!("Failed to match sv_type {}", sv_type_str)
        }
        scf.extend(scf2);
        tool.scfs.push(scf);

        // prepare jxn record for new tags in parallel thread following SEQ recovery
        vcf_out.translate(&mut jxn);
        tool.jxns.push(jxn);      
    };

    // return the input header for modifcation by the VCF writer
    w.ctrs.add_to(N_OUTPUT_JXN_RECORDS, tool.jxns.len());
    w.ctrs.add_to(N_UNIQUE_QNAMES, tool.seqs.len());
    Ok(vcf_out)
}

/// Establish a bcf::Writer updates to the required new tags:
fn get_bcf_writer(tool: &DigestionTool, mut header: bcf::Header) -> Result<bcf::Writer, Box<dyn Error>>{
    for re in &tool.res {
        let record = format!("##INFO=<ID={},Number=1,Type=Integer,Description=\"{} fragment size\">", re.size_tag_str, re.name);
        header.push_record(record.as_bytes());  
        let record = format!("##INFO=<ID={},Number=1,Type=Integer,Description=\"{} nearest site distance\">", re.dist_tag_str, re.name);
        header.push_record(record.as_bytes());  
    }
    let vcf = bcf::Writer::from_stdout(
        &header, 
        false, // compress output with BGZF
        bcf::Format::Vcf
    )?;
    Ok(vcf)
}

/// Fetch a reference sequence given chrom, 1-based start, and 1-based end.
fn get_ref_seq(
    faidx:  &faidx::Reader, 
    chrom:  &str, 
    start1: u32, 
    end1:   u32,
    rc:     bool,
) -> Result<Vec<u8>, Box<dyn Error>> {
    let seq = faidx.fetch_seq(
        chrom, 
        (start1 - 1) as usize, 
        (end1   - 1) as usize
    )?;
    let mut seq = seq.to_ascii_uppercase();
    if rc { 
        seq = seq.iter().rev().map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _    => b'N',
        }).collect(); 
    }
    Ok(seq)
}

// -----------------------------------------------------------------------------
// SEQ extraction
// -----------------------------------------------------------------------------
/// Fill the SEQs for the set of unique QNAMEs found in any SV read into the tool.
fn collect_read_seqs(w: &mut Workflow, tool: &mut DigestionTool) -> Result<(), Box<dyn Error>> {
    w.log.print("extracting read sequences");

    // initialize reading
    let benchmark_bam = w.cfg.get_string(BENCHMARK_BAM).to_string();
    let mut bam = bam::Reader::from_path(benchmark_bam)?;
    bam.set_threads(tool.n_cpu)?;
    let mut aln = bam::Record::new();

    // collect SV-containing read sequences on first encounter in BAM
    while let Some(result) = bam.read(&mut aln) {
        result?;
        match tool.seqs.get_mut(aln.qname()){
            Some(seq) => if seq.is_empty() {
                seq.extend(aln.seq().as_bytes()); // orientation does not matter for mini-consensus
            },
            None => continue
        }
    }
    Ok(())
}

// -----------------------------------------------------------------------------
// site finding in read consensuses
// -----------------------------------------------------------------------------
/// Resolve consensus sequences of SV reads on their scaffolds and perform RE
/// site finding.
fn resolve_consensuses(
    w: &mut Workflow,
    tool: &DigestionTool,
    mut vcf_out: bcf::Writer,
) -> Result<(), Box<dyn Error>> {
    w.log.print("resolving consensuses and performing RE site analysis");

    // create channels for parallel processing
    let (tx_jxn_data, rx_jxn_data) =    
        bounded::<(bcf::Record, usize, &DigestionTool)>(CHANNEL_CAPACITY);
    let (tx_jxn_mod, rx_jxn_mod) =
        bounded::<bcf::Record>(CHANNEL_CAPACITY);
    let n_worker_threads = (*w.cfg.get_usize(N_CPU))
        .saturating_sub(1)
        .max(1);

    // spawn worker threads
    crossbeam::scope(|scope| -> Result<(), Box<dyn Error>> {
        
        // Workers consume junctions and return modified records.
        for worker_i0 in 0..n_worker_threads {
            let rx_jxn_data = rx_jxn_data.clone();
            let tx_jxn_mod = tx_jxn_mod.clone();
            scope.spawn(move |_| {
                process_jxn_tags(worker_i0, rx_jxn_data, tx_jxn_mod)
                    .expect("Junction processing failed");
            });
        }
        drop(tx_jxn_mod);
        drop(rx_jxn_data);

        // Produce inputs concurrently with output collection. Clone each
        // record so tool.jxns remains available to the processing code.
        scope.spawn(move |_| {
            for (jxn_i0, jxn) in tool.jxns.iter().enumerate() {
                tx_jxn_data
                    .send((jxn.clone(), jxn_i0, tool))
                    .expect("Junction workers disconnected");
            }
            drop(tx_jxn_data);
        });

        // Write results on this thread. If writing fails, keep draining the
        // channel so workers cannot deadlock on a full output channel.
        for jxn in rx_jxn_mod.iter() {
            vcf_out.write(&jxn)?
        }
        Ok(())
    }).expect("Crossbeam scope panicked")
}

/// Process site tags from a single junction in parallel. Modify the owned 
/// record and send it back with new tags added.
fn process_jxn_tags(
    worker_i0: usize,
    rx_jxn_data: Receiver<(bcf::Record, usize, &DigestionTool)>,
    tx_jxn_mod:  Sender<bcf::Record>,
) -> Result<(), Box<dyn Error>> {

    // one consensus resolver per thread
    let (mut resolver, mut poa_pool) = Resolver::with_capacity(
        ResolverConfig::default(),
        100, 
        (SCAFFOLD_FLANK_LEN * 2 + 1000) as usize,
    );

    // process one kept junction at a time
    for (mut jxn, jxn_i0, tool) in rx_jxn_data.iter() {
        let scf = &tool.scfs[jxn_i0];

        // initialize the SV record for site analysis
        let Some(qnames) = get_vcf_qnames(&jxn)? else { continue };
        let downsampled_qnames = &qnames[0..qnames.len().min(MAX_CONSENSUS_SEQS)];

        // build the SV allele haplotype consensus sequence
        let seqs: Vec<_> = downsampled_qnames.iter().map(|qname|{
            tool.seqs.get(qname).unwrap().to_vec()
        }).collect();
        resolver.set_scaffold_with_aligner(scf);
        resolver.par_set_seqs(seqs, |_, _| true);
        let haplotype = resolver.get_consensus(&mut poa_pool);
        let haplotype_str = unsafe{ from_utf8_unchecked(&haplotype) };

        // locate the junction in the consensus and characterize de novo base insertions
        let minimap2 = Minimap2::builder()
            .map_hifi()
            .with_cigar() // required if CS tag is desired
            .with_seq(scf)
            .expect("Failed to build minimap2 hap_on_scf index");
        let hap_on_scf = minimap2.map(
            &haplotype, true, false, 
            None, None, None
        )?;
        let mut hap_on_scf_map = vec![None::<usize>; scf.len()];
        for mapping in &hap_on_scf {
            if !mapping.is_primary { break }
            process_cs_tag(&mut hap_on_scf_map, mapping);
        }
        let mut jxn_pos0 = None::<usize>;
        if let (
            Some(mut left_bkpt0), 
            Some(mut right_bkpt0)
        ) = (
            hap_on_scf_map[..SCAFFOLD_FLANK_LEN as usize].iter().rev().find_map(|&x| x),
            hap_on_scf_map[scf.len() - SCAFFOLD_FLANK_LEN as usize..].iter().find_map(|&x| x)
        ) {
            if left_bkpt0 > right_bkpt0 {
                (left_bkpt0, right_bkpt0) = (right_bkpt0, left_bkpt0)
            }
            jxn_pos0 = Some(left_bkpt0 + (right_bkpt0 - left_bkpt0) / 2);
        }

        // get fragment metadata for each RE
        for re in &tool.res {
            let (frag_size, min_site_dist) = if let Some(jxn_pos0) = jxn_pos0 {
                get_site_data(haplotype_str, &re.site, jxn_pos0)
            } else {
                (-1, -1)
            };
            jxn.push_info_integer(&re.size_tag, &[frag_size])?;
            jxn.push_info_integer(&re.dist_tag, &[min_site_dist])?;
        }

        // write the final VCF record
        tx_jxn_mod.send(jxn)?;
    };

    eprintln!("Worker {} done.", worker_i0 + 1);
    Ok(())
}

/// Process a cs:Z:tag to help fill hap_on_scf_map from one primary mapping.
fn process_cs_tag(
    hap_on_scf_map: &mut Vec<Option<usize>>,
    mapping: &Mapping,
) { 
    let Some(aln) = &mapping.alignment else { return };
    let Some(cs) = &aln.cs else { return };
    let mut hap_pos0 = mapping.query_start as usize;
    let mut scf_pos0 = mapping.target_start as usize;
    let mut chars = cs.chars();
    let mut cs_op = chars.next().unwrap();
    let mut op_val = String::new();
    while let Some(char) = chars.next() {
        if char.is_alphanumeric() {
            op_val.push(char);
        } else {
            handle_cs_op(hap_on_scf_map, cs_op, &op_val, &mut hap_pos0, &mut scf_pos0);
            cs_op = char;
            op_val.clear();
        }
    }
    handle_cs_op(hap_on_scf_map, cs_op, &op_val, &mut hap_pos0, &mut scf_pos0);
}

/// Process one cs:Z:tag operation to add haplotype position to hap_on_scf_map. 
fn handle_cs_op(
    hap_on_scf_map: &mut Vec<Option<usize>>,
    cs_op : char,
    op_val: &str,
    hap_pos0: &mut usize, 
    scf_pos0: &mut usize,
) {
    match cs_op {
        // :[0-9]+   Identical sequence length
        ':' => {
            let op_len = op_val.parse::<usize>().unwrap();
            let mut hap_pos: std::ops::Range<_> = *hap_pos0..*hap_pos0 + op_len;
            hap_on_scf_map[*scf_pos0..*scf_pos0 + op_len].fill_with(|| hap_pos.next());
            *hap_pos0 += op_len;
            *scf_pos0 += op_len;
        },
        // *[acgtn][acgtn]   Substitution: target to query
        '*' => {
            //     S
            // rrrrRrrrr
            // qqqqQqqqq
            //     A
            hap_on_scf_map[*scf_pos0] = Some(*hap_pos0);
            *hap_pos0 += 1;
            *scf_pos0 += 1;
        },
        // +[acgtn]+   Insertion to the target
        '+' => {
            //    *INI     insertions may have heteroduplex bases within homoduplex query run
            // rrrr   Rrrr
            // qqqqQqqqqqq
            //    aA Aa
            let op_len = op_val.len();
            *hap_pos0 += op_len;
        },
        // -[acgtn]+   Deletion from the target
        '-' => {
            //     DDD
            // rrrrRrrrrrr
            // qqqq   Qqqq
            //   aA   Aa
            let op_len = op_val.len();
            hap_on_scf_map[*scf_pos0..*scf_pos0 + op_len].fill(Some(*hap_pos0 - 1));
            *scf_pos0 += op_len;
        },
        _   => panic!("Unexpected operation in cs tag: {}", cs_op),
    }
}

/// Find the RE sites to the left and to the right of the ~center of the junction
/// in the haplotype consensus; return (fragment size, nearest site distance).
fn get_site_data(haplotype_str: &str, site: &str, jxn_pos0: usize) -> (i32, i32) {
    // -------------------------|---------------------------
    //   ######               ######            ######                     
    let half_site_len = site.len() / 2;

    let boundary = (jxn_pos0 + half_site_len).min(haplotype_str.len());
    let left_site_pos = haplotype_str[..boundary]
        .rfind(site)
        .map(|i0| i0 + half_site_len)
        .unwrap_or(0);

    let boundary = jxn_pos0.saturating_sub(half_site_len);
    let mut right_site_pos = haplotype_str[boundary..]
        .find(site)
        .map(|i0| boundary + i0 + half_site_len)
        .unwrap_or(haplotype_str.len());
    if right_site_pos < left_site_pos { right_site_pos = left_site_pos }

    let frag_size = right_site_pos - left_site_pos;
    let left_site_dist  = jxn_pos0.saturating_sub(left_site_pos);
    let right_site_dist = right_site_pos.saturating_sub(jxn_pos0);

    (frag_size as i32, left_site_dist.min(right_site_dist) as i32)
}
