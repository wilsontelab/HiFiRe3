//! Final junction handling during benchmarking.

// imports
// use std::error::Error;
use crate::tools::compare_sv_calls::ComparisonTool;
use super::*;

/// Convert a HiFiRe3 chrom index into a non-composite canonical chrom name.
pub fn hf3_chrom_index_to_name(
    chrom_index1: u8, 
    tool: &ComparisonTool, 
) -> Option<String> {
    simplify_chrom_name(
        &tool.chroms.rev_index[&chrom_index1], 
        tool.is_composite_genome
    )
}
