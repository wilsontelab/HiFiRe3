//! Merge the two strands of PacBioStrand CCS reads into a single sequence.
// input:
//     "fail" and "hifi" CCS reads as FAIL_BAM_FILE, HIFI_BAM_FILE
//      called with the by-strand option, and optionally with kinetics tags
// output:
//     READ_BAM_FILE, with:
//         - one read per original ZMW, either:
//              - merging both strands for error corrected SNV and SV detection
//              - passing a single usable strand as is for SV detection only
//              - suppressing reads with no usable strands
//     PACBIO_BASECALL_KINETICS file with kinetics values stratified by 3-base context

// modules
mod channel_reader;
mod channel_worker;
mod channel_kinetics;

// dependencies
use std::error::Error;
use crossbeam::channel::bounded;
use minimap2::{Aligner as Minimap2};
use faimm::IndexedFasta;
use rust_htslib::bam::{Reader, Read, Writer, Record as BamRecord, record::Aux, Header, Format};
use mdi::pub_key_constants;
use mdi::workflow::{Config, Counters};
use crate::formats::hf3_tags::*;
use channel_kinetics::kinetics::{StrandMetadata, FrameCounts};
use channel_reader::usable::{UnusableReason, StrandSource};

/// BufferedStrand holds the needed data for a single PacBio strand.
struct BufferedStrand {
    pub source:   StrandSource, 
    pub usable:   bool,
    pub reason:   UnusableReason,
    pub ff:       u8, // same names as PacBio tags
    pub ec:       f32,
    pub seq:      String,
    pub qual:     Vec<u8>, 
    pub ip:       Option<Vec<u8>>, // inter-pulse durations
    pub pw:       Option<Vec<u8>>, // pulse widths
}

// MergePair holds the data needed to merge two strands.
struct StrandPair {
    sources: StrandSource,   // strand pairs might come from one file or one strand from each of fail and hifi files
    qname:   Vec<u8>,
    ff:      u8,
    this:    BufferedStrand, // ultimately becomes the forward reference strand whenever possible
    prev:    BufferedStrand, // the partner strand encountered previously, before encountering this strand (unless reordered by strand)
}

/// KineticsInstance supports kinetics data collection in the side-effect channel.
struct KineticsInstance {
    strand_metadata: StrandMetadata,
    frame_counts:    FrameCounts,
}

/// MergeResult carries the outcome of attempting to merge two strands.
struct MergeResult {
    sources: StrandSource, // strand pairs might come from one file or one strand from each of fail and hifi files
    reason:  &'static str, // string encoding of one or two UnusableReasons
    outcome: &'static str, // string encoding of the strand merging outcome
    qname:   Vec<u8>,
    seq:     String,
    qual:    Vec<u8>,
    ff:      u8,  // bitwise OR of the two strand's fail flags (or just this strand if single)
    ec:      f32, // summed coverage of the two strands (or just one strand if reporting single strand as the read)
    dt:      Option<u8>,
    dd:      Option<String>,
    sk:      Option<Vec<u16>>,
}

// constants for environment variable, config, and counter keys, etc.
const TOOL: &str = "basecall_pacbio";
pub_key_constants!(
    // from environment variables
    N_CPU
    HIFI_BAM_FILE // used in main() for header copying
    READ_BAM_FILE // the output BAM file for merged reads
    GENOME_FASTA
    MINIMAP2_INDEX_WRK
    // counter keys
    N_READS
    N_READS_BY_REASON  // strand usability reasons
    FAIL_BY_REASON
    HIFI_BY_REASON
    BOTH_BY_REASON
    N_READS_BY_OUTCOME // strand merging outcomes
    FAIL_BY_OUTCOME
    HIFI_BY_OUTCOME
    BOTH_BY_OUTCOME
    N_ORPHAN_STRANDS
    // merge outcomes
    UNUSABLE_READ
    UNMAPPED_READ
    ORPHAN_STRAND
    // correction to reference types
    CORRECTED_TO_REFERENCE
);
const CHANNEL_CAPACITY: usize = 16; // channel buffer size
const MIN_READ_LEN: usize = 250;   // TODO: expose read length limits as options?
const MAX_READ_LEN: usize = 25000; // MAX_READ_LEN impacts the size of the Smith-Waterman buffer

// main basecall pacbio function called by hf3_tools main()
pub fn main() -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_usize_env(&[N_CPU]);
    cfg.set_string_env(&[HIFI_BAM_FILE, READ_BAM_FILE, GENOME_FASTA, MINIMAP2_INDEX_WRK]);

    // initialize counters
    let mut ctrs = Counters::new(TOOL, &[
        (N_READS,            "number of reads with two opposing strands"),
        (N_ORPHAN_STRANDS,   "number of orphan strands with a QNAME only encountered once in the input files"),
    ]);
    ctrs.add_keyed_counters(&[
        (N_READS_BY_REASON,  "number of reads by usability reason"),
        (FAIL_BY_REASON,     "number of reads with both strands from fail.bam by usability reason"),
        (HIFI_BY_REASON,     "number of reads with both strands from hifi.bam by usability reason"),
        (BOTH_BY_REASON,     "number of reads with one strand from each bam file by usability reason"),
        (N_READS_BY_OUTCOME, "number of reads by strand merging outcome"),
        (FAIL_BY_OUTCOME,    "number of reads with both strands from fail.bam by strand merging outcome"),
        (HIFI_BY_OUTCOME,    "number of reads with both strands from hifi.bam by strand merging outcome"),
        (BOTH_BY_OUTCOME,    "number of reads with one strand from each bam file by strand merging outcome"),
    ]);
    let mut ctrs2 = Counters::new(TOOL, &[]);
    ctrs2.add_keyed_counters(&[
        (CORRECTED_TO_REFERENCE, "number of bases corrected to reference by strand base pattern"),
    ]);

    // instantiate shared minimap2 aligner
    let mm2_index = cfg.get_string(MINIMAP2_INDEX_WRK).to_string();
    let mut minimap2 = Minimap2::builder()
        .map_hifi()
        // .with_cigar() // required if CS tag is desired
        .with_index(mm2_index, None)
        .expect("Failed to build minimap2 index");
    minimap2.mapopt.bw      = 500;
    minimap2.mapopt.bw_long = 3300; // same as the PacBio platform alignment default

    // instantiate the FAI sequence fetcher
    let genome_fasta = cfg.get_string(GENOME_FASTA).to_string();
    let fa = IndexedFasta::from_file(&genome_fasta)
        .expect("Error opening genome FASTA file");

    // read the header from HIFI_BAM_FILE to use for the output
    let hifi_path = cfg.get_string(HIFI_BAM_FILE).to_string();
    let hifi_bam = Reader::from_path(&hifi_path)?;
    let header = Header::from_template(hifi_bam.header());

    // create output BAM header and writer
    let output_path = cfg.get_string(READ_BAM_FILE).to_string();
    let mut bam_writer = Writer::from_path(&output_path, &header, Format::Bam)?;

    // create channels for fan-in/fan-out parallel processing, and kinetics side-effect
    let (tx_strand_pair, rx_strand_pair) 
        = bounded::<StrandPair>(CHANNEL_CAPACITY);
    let (tx_kinetics, rx_kinetics) 
        = bounded::<KineticsInstance>(CHANNEL_CAPACITY);
    let (tx_corr_to_ref, rx_corr_to_ref) 
        = bounded::<String>(CHANNEL_CAPACITY);
    let (tx_merge_result, rx_merge_result) 
        = bounded::<MergeResult>(CHANNEL_CAPACITY);

    // spawn threads
    crossbeam::scope(|scope| {

        // reader: read strands from input BAMs and dispatch for merging
        let tx_mr = tx_merge_result.clone();
        scope.spawn(move |_| {
            channel_reader::stream_bam_files(
                tx_strand_pair, 
                tx_mr, // for short-circuting merge worker
            ).unwrap();
        });

        // workers: merge strand pairs from reader in parallel
        // create various types of merge results and kinetics side-effects
        let n_cpu = *cfg.get_usize(N_CPU);
        let n_workers: usize = n_cpu.max(4) - 3;// reserve one thread for reader, kinetics, and output
        for _ in 0..n_workers {
            let rx_strand_pair = rx_strand_pair.clone();
            let tx_kinetics = tx_kinetics.clone();
            let tx_corr_to_ref = tx_corr_to_ref.clone();
            let tx_merge_result = tx_merge_result.clone();
            scope.spawn(|_| {
                channel_worker::merge_strand_pairs(
                    rx_strand_pair, 
                    tx_kinetics,
                    tx_corr_to_ref,
                    tx_merge_result,
                    &minimap2,
                    &fa,
                ).unwrap();
            });
        }
        drop(tx_kinetics);
        drop(tx_corr_to_ref);
        drop(tx_merge_result); 

        // kinetics side-effect receiver channel
        scope.spawn(move |_| {
            channel_kinetics::collect_kinetics(
                rx_kinetics
            ).unwrap();
        });

        // corrected to reference receiver channel
        scope.spawn(|_| {
            for base_pattern in rx_corr_to_ref {
                ctrs2.increment_keyed(CORRECTED_TO_REFERENCE, &base_pattern);
            }
        });

        // process merge results as they arrive in the main thread
        for merge_result in rx_merge_result {

            // update counters
            ctrs.increment(N_READS);
            ctrs.increment_keyed(N_READS_BY_REASON, merge_result.reason);
            ctrs.increment_keyed(merge_result.sources.to_reason_key(),  merge_result.reason);
            ctrs.increment_keyed(N_READS_BY_OUTCOME, merge_result.outcome);
            ctrs.increment_keyed(merge_result.sources.to_outcome_key(), merge_result.outcome);

            // write merged read if at least one usable strand
            if merge_result.outcome != UNUSABLE_READ &&
               merge_result.outcome != UNMAPPED_READ &&
               merge_result.outcome != ORPHAN_STRAND {
                let mut read = BamRecord::new();
                read.set(
                    &merge_result.qname, 
                    None, 
                    merge_result.seq.as_bytes(), 
                    &merge_result.qual
                );
                read.push_aux(
                    PACBIO_FAIL.as_bytes(), 
                    Aux::U8(merge_result.ff)
                ).unwrap();
                read.push_aux(
                    PACBIO_EFF_COVERAGE.as_bytes(), 
                    Aux::Float(merge_result.ec)
                ).unwrap();
                if let Some(dt) = merge_result.dt {
                    read.push_aux( // include dt even if 0 as a flag that duplex basecalling was performed
                        STRAND_DIFFERENCE_TYPES.as_bytes(),
                        Aux::U8(dt)
                    ).unwrap();
                }
                if let Some(dd) = merge_result.dd {  
                    read.push_aux(
                        STRAND_DIFFERENCES.as_bytes(), 
                        Aux::String(&dd)
                    ).unwrap();
                }
                if let Some(sk) = merge_result.sk { 
                    if sk.len() > 0 {
                        read.push_aux(
                            SUBSTITUTION_KINETICS.as_bytes(), 
                            Aux::ArrayU16(sk.as_slice().into())
                        ).unwrap();
                    }
                }
                bam_writer.write(&read).unwrap();
            } else if merge_result.outcome == ORPHAN_STRAND {
                ctrs.increment(N_ORPHAN_STRANDS);
            }
        }
    }).expect("Crossbeam scope panicked");

    // print counts
    ctrs.print_grouped(&[
        &[N_READS, N_ORPHAN_STRANDS],
        &[N_READS_BY_REASON],
        &[FAIL_BY_REASON],
        &[HIFI_BY_REASON],
        &[BOTH_BY_REASON],
        &[N_READS_BY_OUTCOME],
        &[FAIL_BY_OUTCOME],
        &[HIFI_BY_OUTCOME],
        &[BOTH_BY_OUTCOME],
    ]);
    ctrs2.print_grouped(&[
        &[CORRECTED_TO_REFERENCE],
    ]);
    Ok(())
}
