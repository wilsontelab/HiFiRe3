//! Module to handle SNV/indel analysis according to expectations and encoding 
//! of HiFiRe3 libraries. 

// modules 
pub mod fragment;
pub mod dd_tag;
pub mod cs_tag;
pub mod simple_repeat;
pub mod variant;
pub mod strand_data;
pub mod analyze_reads;
pub mod haplotype;
pub mod encoding;

// re-exports
pub use fragment::*;
pub use dd_tag::*;
pub use simple_repeat::*;
pub use variant::*;
pub use strand_data::*;
pub use haplotype::*;
pub use encoding::*;

// imports
use std::error::Error;
use rustc_hash::{FxHashMap, FxHashSet};
use minimap2::{Aligner as Minimap2, PresetSet};
use faimm::IndexedFasta;
use serde::Serialize;
use mdi::pub_key_constants;
use mdi::workflow::Config;
use genomex::genome::{TargetRegions, Exclusions};
use crate::snvs::analyze_reads::poa::*;
use crate::tools::type_aliases::*;

// constants
pub_key_constants!(
    // from environment variables
    SEQUENCING_PLATFORM
    LIBRARY_TYPE
);
pub const MIN_HAPLOTYPE_READS: u16 = 2;  // a haplotype must have >=2 matching reads to be called heterozygous
pub const MIN_SNV_INDEL_QUAL: u8 = 27;
pub const MAX_EXPECTED_READ_LEN:  usize = 10000; // use for allocating recycled objects

/// Ensure that PacBio SNV analysis is performed on a library from the 
/// PacBioStrand sequencing platform.
pub fn check_pacbio_strand(tool: &str, cfg: &mut Config) -> Result<(), Box<dyn Error>> {
    cfg.set_string_env(&[SEQUENCING_PLATFORM, LIBRARY_TYPE]);
    let sequencing_platform = cfg.get_string(SEQUENCING_PLATFORM);
    let library_type        = cfg.get_string(LIBRARY_TYPE);
    if sequencing_platform != "PacBioStrand" || 
       library_type        != "HiFi" {
        return Err(format!(
            "{} requires PacBioStrand HiFi reads; found {} {} reads", 
            tool, sequencing_platform, library_type).into()
        );
    }
    Ok(())
}

/// SnvAnalysisTool collects tools for SNV analysis shared with all chromosome 
/// workers.
pub struct SnvAnalysisTool {

    // global configuration parameters
    pub n_cpu: u32,

    // chromosomes and regions
    pub targets:    TargetRegions,
    pub exclusions: Exclusions,
    pub fa:         IndexedFasta,
}

/// SnvChromWorker collects tools for SNV analysis that are specific to 
/// processing of a single specific chromosome in parallel.
pub struct SnvChromWorker{

    // chromosome parsing
    pub chrom:       ChromName,
    pub chrom_index: ChromIndex1,
    pub chrom_tid:   usize,
    pub simple_repeats: SimpleRepeats,

    // option parameters
    pub min_fragment_reads:   usize,
    pub min_homozygous_reads: usize,

    // processing tools
    pub poa:            Poa,
    pub minimap2:       Minimap2<PresetSet>,
    pub frag_vars:      FragmentVariants,
    pub encoding:       AlignmentEncoding,  // read encoding for visualization
    pub tracking_variants: Vec<Variant>, // potentially heterozygous variants
    pub seq0_bases:     Vec<String>, // used with cs_map for consensus calling
    pub cs_map:         Vec<FxHashMap<String, u8>>,
    pub hap_vs_ref:     Vec<String>, // consensus encoding for visualization
    pub hap_vars:       FxHashMap<Haplotype, FxHashSet<Variant>>,
    pub hap_votes:      FxHashMap<Haplotype, usize>,
    pub str_matches:    Vec<bool>,
    pub var_tgt_pos0:   Option<SeqPos0>, 
    pub tgt_bases:      UppercaseACGTN,
    pub alt_bases:      UppercaseACGTN,
    pub min_qual:       u8,
    pub allowed:        bool,
    pub cs_op:          char,
    pub op_val:         String,
    // pub debug:      ReFragment,
    // pub show_debug: bool,

    // metadata aggregation
    pub variant_tally:       VariantsTally,
    pub variant_reads_tally: VariantReadsTally,
}

/// SnvChromWorkerData allows difference types of metadata to be trasmitted to 
/// the main thread for aggregation over the entire input.
pub enum SnvChromWorkerData {
    TotalAlnCount(usize),
    UsableAlnCount((ChromName, usize)),
    VariantMetadata(VariantMetadata),
    VariantReadsMetadata(VariantReadsMetadata),
    ReadsOnReferenceMetadata(FragmentHaplotypeMetadata),
    ReadsOnHaplotypeMetadata(FragmentHaplotypeMetadata),
}
