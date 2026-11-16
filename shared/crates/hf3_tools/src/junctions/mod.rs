//! Module to handle SV junction analysis according to expectations
//! and encoding of HiFiRe3 libraries. 
//
// READ PATHS - a set of observed nodes over all alignments per QNAME, plus a projection (n_nodes = n_alns + 1 + 1 projection)
//         >---------------------->                   end-to-end read path,     canonical orientation, no junctions
//         >---------------->~~~~~>                   projected read path,      canonical orientation, no junctions
//         <----------------------<                   end-to-end read path, non-canonical orientation, no junctions
//         <~~~~~~~<--------------<                   projected read path,  non-canonical orientation, no junctions
//         >----------------------o------------>      partial digestion path,   canonical orientation, no junctions
//  >------X---------------------->                   chimeric read path, false junction
//         >-------------|--------------|-----------> SV read path, canonical orientation, two junctions
//         >--------------|-------------|-----------> similar SV read path, one exact junction match, one within match tolerance
//         <--------------|-------------|-----------< similar SV read path, non-canonical orientation
//         >--------------|----->~~~~>                partial SV read path, projected, one junction
//  >------X--------------|-------------|-----------> chimeric SV read path, one false junction, two true junctions

// JUNCTIONS - two ordered breakpoint nodes, with metadata (offset, orientation, count, QNAMEs, etc.)
//  ---|-----  two junctions with counts before fuzzy matching, canonical orientation
//  ----|----  collapsed to one junction after fuzzy matching
//  ----X----  chimeric junctions persist with metadata for rejecting them

// ALIGNMENTS - a unique reference span with nodes at each end, orientation, index in read
//  >---------------------->   various alignment segments
//  >---------------->~~~~~>   each with 3 nodes: end5_observed, end3_observed, and end3_projected
//  <----------------------<   end5_observed and end3_projected taken as site_pos when EXPECTING_ENDPOINT_RE_SITES
//  <~~~~~~~<--------------<
//  >-------------|~~~~~>
//  <--------------|
//  |-------------|
//  |----------->
//  >------X
//  X--------------|

// modules 
mod flag;
mod junction;
mod read_path;
mod alignment;
mod grouping;
mod deduplication;
mod samples;
mod vcf;

// re-exports
pub use flag::*;
pub use junction::*;
pub use read_path::*;
pub use alignment::*;
pub use grouping::*;
pub use samples::*;

// dependencies
use crossbeam::channel::Sender;
use rust_htslib::bam::HeaderView;
use genomex::genome::{Chroms, TargetRegions, Genes, Exclusions};

/// The JunctionAnalysisTool collects structs and methods for SV analysis
/// at the genome level after all chromosomes have been processed.
pub struct JunctionAnalysisTool {

    // global configuration parameters
    pub n_cpu:     u32,

    // chromosomes and regions
    pub chroms:     Chroms,
    pub targets:    TargetRegions,
    pub genes:      Genes,
    pub exclusions: Exclusions,

    // grouping parameters
    pub group_breakpoint_distance: usize,
    pub group_stem_distance:       u32,

    // deduplication parameters
    pub is_ont: bool,
    pub deduplicate_reads: bool,

    // reporting parameters
    pub min_n_observed: u16,
    pub min_sv_mapq:    u8,

    // final junction output files
    pub final_jxns_file_1: String,
    pub final_jxns_file_2: String,
}

/// The JunctionChromWorker tool collects structs and methods 
/// for SV analysis while processing a single chromosome.
pub struct JunctionChromWorker <'a, 'b>{

    // global configuration parameters
    pub expecting_endpoint_re_sites: bool,

    // chromosome parsing
    pub chroms: &'a Chroms, // HiFiRe3 ordered chroms
    pub header: HeaderView, // BAM header contigs

    // structures for aggregating over all reads
    pub first_alns: Vec<AlignmentSegment>,

    // result channel
    pub tx_data: &'b Sender<ChromWorkerData>,
}

// ChromWorkerData enum allows difference types of data to be
// transmitted to the main thread for aggregation.
pub enum ChromWorkerData {
    ReadPath(SvReadPath),
    DistalAln(AlignmentSegment),
    Junction((OrderedJunction, JunctionInstance)),
    ChromReadCount((String, usize)),
}
