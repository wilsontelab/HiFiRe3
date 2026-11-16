//! SV/junction matching during benchmarking.

// imports
// use std::error::Error;
use std::str::from_utf8_unchecked;
use serde::Serialize;
use rust_htslib::bcf;
use genomex::sam::junction::JunctionType;
use crate::junctions::FinalJunction;
use crate::tools::compare_sv_calls::ComparisonTool;
use super::*;

// -----------------------------------------------------------------------------
// Handle mixtures of composite and non-composite genomes
// -----------------------------------------------------------------------------
/// Convert a chrom name into a primary genome canonical chrom name.
pub fn simplify_chrom_name(
    chrom: &str, 
    is_composite_genome: bool
) -> Option<String> {
    let chrom = if is_composite_genome {
        let parts: Vec<_> = chrom.split("_").collect();
        let last_i0 = parts.len() - 1;
        if parts[last_i0] == "dm6" { return None } // drop fly composite alignments
        parts[0..last_i0].join("_") // drop the genome suffix
    } else {
        chrom.to_string()
    };
    if chrom.contains("_") { 
        None // drop non-canonical chroms 
    } else { 
        Some(chrom) 
    }
}

/// Check whether an SV breakpoint is on target.
pub fn sv_is_on_target(
    chrom_1: &str,
    pos1_1:  u32,
    chrom_2: &str,
    pos1_2:  u32,
    tool: &ComparisonTool
) -> bool {
    tool.targets.get_pos_region(chrom_1, pos1_1).2 == INSIDE_REGION ||
    tool.targets.get_pos_region(chrom_2, pos1_2).2 == INSIDE_REGION
}

/// Check whether an SV breakpoint is excluded.
pub fn sv_is_excluded(
    chrom_1: &str,
    pos1_1:  u32,
    chrom_2: &str,
    pos1_2:  u32,
    tool: &ComparisonTool
) -> bool {
    tool.primary_genome_exclusions.pos_in_region(&chrom_1, pos1_1) ||
    tool.primary_genome_exclusions.pos_in_region(&chrom_2, pos1_2)
}

// -----------------------------------------------------------------------------
// SV sources and types for HF3 vs. VCF comparison
// -----------------------------------------------------------------------------
/// Enum for tracking SV/junction sources.
#[repr(u8)]
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize)]
pub enum SvSource {
    #[serde(rename = "HF3")]
    HiFiRe3 = 1,
    #[serde(rename = "SNF")]
    Sniffles = 2,
    #[serde(rename = "CSV")]
    CuteSV = 3,
}

/// Enum for normalizing SV/junction types.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Serialize)]
#[serde(rename_all = "UPPERCASE")]
pub enum SvType {
    Del,
    Dup,
    Inv,
    Trans,
    Ins,
}
impl SvType {
    pub fn from_vcf(jxn: &bcf::Record, chrom_1: &str, alt: &[u8]) ->  (Self, Vec<u8>) {
        let sv_type_vcf = get_vcf_sv_type_str(jxn);
        let sv_type = match sv_type_vcf.as_str() {
            "DEL" => SvType::Del,
            "DUP" => SvType::Dup,
            "INV" => SvType::Inv,
            "INS" => SvType::Ins,
            "BND" => {
                let bnd_alt = unpack_bnd_alt(alt).expect("Failed to unpack BND ALT."); 
                if chrom_1 == bnd_alt.chrom {
                    match &bnd_alt.bnd_type {
                        BndType::AfterRight    => SvType::Del,
                        BndType::AfterRcLeft   => SvType::Inv,
                        BndType::BeforeLeft    => SvType::Dup,
                        BndType::BeforeRcRight => SvType::Inv,
                    }
                } else {
                    SvType::Trans
                }
            },
            _ => panic!("Unknown VCF SVTYPE: {}", sv_type_vcf),
        };
        (sv_type, sv_type_vcf.as_bytes().to_vec())
    }
    pub fn from_final_junction(item: &FinalJunction) -> Self {
        match item.jxn_type {
            x if x == JunctionType::Deletion as u8 => SvType::Del,
            x if x == JunctionType::Duplication as u8 => SvType::Dup,
            x if x == JunctionType::Inversion as u8 => SvType::Inv,
            x if x == JunctionType::Translocation as u8 => SvType::Trans,
            _ => panic!("Proper final junction type not expected.")
        }
    }
}

// -----------------------------------------------------------------------------
// SV junction calls
// -----------------------------------------------------------------------------
/// Table for aggregating SV/junction calls.
#[derive(Serialize)]
pub struct JunctionCall {
    pub jxn_id:  u32,

    pub sv_type: SvType,
    pub strands: String,

    pub source:  SvSource,

    pub chrom_1: String,
    pub pos1_1:  u32,
    pub chrom_2: String,
    pub pos1_2:  u32, // VCF as is, HiFiRe3 adjusted to Sniffles 
    pub pos1_2_csv:  u32, // HiFiRe3 adjusted to CuteSV

    pub sv_size:     u32,
    pub frag_size:   u32,
    pub stem_length: u32,

    pub offset:  i16,
    pub jxn_seq: String,

    pub n_observed: u16,
    pub mapq: u8,
    pub divergence: f32,

    pub jxn_fail: u8,
    pub aln_fail: u8,
    pub all_chimeric: u8,

    pub mosaic: bool,
    pub vaf: f32,
    
    pub qnames: String,
}

/// Aggergated algorithm outcomeos for grouped calls.
#[repr(u8)]
pub enum AlgorithmOutcomes{
    HiFiRe3  = 1,
    Sniffles = 2,
    CuteSV   = 4,
}

/// Table for aggregating SV/junction calls.
#[derive(Serialize)]
pub struct JunctionCallCollapsed {
    pub jxn_id:  u32,

    pub sv_type: SvType,
    pub strands: String,

    pub n_algorithms: u8,
    pub algorithms: u8,
    pub hf3:  u16, // n_observed
    pub snf:  u16,
    pub csv:  u16,

    pub chrom_1: String,
    pub pos1_1:  u32,
    pub chrom_2: String,
    pub pos1_2:  u32, // VCF as is, HiFiRe3 adjusted to Sniffles 
    pub pos1_2_csv:  u32, // HiFiRe3 adjusted to CuteSV

    pub sv_size:       u32,
    pub frag_size_hf3: u32,
    pub frag_size_vcf: u32,

    pub offset:  i16,
    pub jxn_seq: String,

    pub mapq: u8,
    pub divergence: f32,

    pub jxn_fail: u8,
    pub aln_fail: u8,
    pub all_chimeric: u8,

    pub mosaic: bool,
    pub vaf: f32,
    
    pub qnames: String,
}

// -----------------------------------------------------------------------------
//sorted QNAME lists
// -----------------------------------------------------------------------------
/// Assemble a sorted qnames concatenation from a comma-delimited String (HF3).
pub fn get_sorted_qnames_str(qnames: &str) -> String {
    let mut qnames: Vec<_> = qnames
        .split(",")
        .map(str::trim)
        .filter(|&s| !s.is_empty())
        .collect();
    sort_and_concatenate_qnames(&mut qnames[1..]) // since formatted as `,xxx,...`
}
/// Assemble a sorted qnames concatenation from a Vec<Vec<u8>> (VCF).
pub fn get_sorted_qnames_bytes(qnames: Vec<Vec<u8>>) -> String {
    let mut qnames: Vec<_> = qnames.iter().map(|x|{
        unsafe{ from_utf8_unchecked(x) }
    }).collect();
    sort_and_concatenate_qnames(&mut qnames)
}
fn sort_and_concatenate_qnames(qnames: &mut [&str]) -> String {
    qnames.sort_unstable();
    qnames.join(",")
}
