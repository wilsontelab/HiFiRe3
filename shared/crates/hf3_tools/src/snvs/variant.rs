/// Support for calling and counting specific variants from error-corrected reads.

// imports
use serde::{Serialize, Serializer};
use super::*;

/// A Variant encodes a specific SNV or indel, or a series of operations,
/// observed beginning at a single reference position on a known chromosome
/// or on a single resolved haplotype.
/// 
/// A Variant allows any number of reference bases to be replaced by any number 
/// of non-reference bases, so it is equally capable of representing 
/// substitutions, insertions, deletions, and complex indels. 
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Debug)]
pub struct Variant {
    // leftmost reference position when n_tgt_bases > 0, or the position preceding an insertion
    pub ref_pos0: ChromPos0,
    // the coordinate on target (either reference or haplotype) that matches ref_pos0
    pub tgt_pos0: SeqPos0,
    // for substitutions and deletions, the expected bases replaced by alt_bases         
    pub tgt_bases: Option<UppercaseACGTN>,
    // for substitutions and insertions, the bases replacing the expected bases    
    pub alt_bases: Option<UppercaseACGTN>,
    // whether the variant in an indel
    #[serde(serialize_with = "serialize_indel")]
    pub is_indel: bool,
    pub re_fragment: ReFragment,
    // whether tgt_pos0 is relative to a reference chromosome or a haplotype consensus
    #[serde(serialize_with = "serialize_haplotype")]
    pub haplotype: Haplotype,
}
/// Helper function to serialize is_indel as u8.
pub fn serialize_indel<S: Serializer>(
    b: &bool, 
    serializer: S
) -> Result<S::Ok, S::Error>{
    serializer.serialize_u8(*b as u8)
}
impl Variant {
    /// Create a new Variant instance with the specified fields.
    pub fn new(
        ref_pos0:    ChromPos0, 
        tgt_pos0:    SeqPos0, 
        tgt_bases:   &str,
        alt_bases:   &str,
        re_fragment: &ReFragment,
        haplotype:   Haplotype,
    ) -> Self {
        Variant {
            ref_pos0,
            tgt_pos0,
            tgt_bases: if tgt_bases.is_empty() {
                None 
            } else { 
                Some(tgt_bases.to_string()) 
            },
            alt_bases: if alt_bases.is_empty() {
                None 
            } else { 
                Some(alt_bases.to_string()) 
            },
            is_indel: tgt_bases.len() as usize != alt_bases.len(),
            re_fragment: *re_fragment,
            haplotype,
        }
    }

    /// Get the signed difference in ref vs. alt length.
    pub fn alt_minus_ref(&self) -> i32 {
        self.alt_bases.as_ref().map_or(0, |alt| alt.len() as i32) - 
        self.tgt_bases.as_ref().map_or(0, |alt| alt.len() as i32)
    }

    /// Pack a variant into a string representation for printing to 
    /// variant_reads file.
    pub fn to_string(
        &self, 
        tgt_start0: SeqPos0, 
        min_qual:   PhredQual
    ) -> String {
        let tgt_bases = self.tgt_bases
            .as_deref()
            .unwrap_or("-");
        let alt_bases = self.alt_bases
            .as_deref()
            .unwrap_or("-");
        format!(
            "{}:{}:{}:{}", 
            tgt_start0 + self.tgt_pos0,
            tgt_bases,
            alt_bases,
            min_qual
        )
    }
}

/// A VariantLocation records just the position and type of variant in a 
/// ReFragment. It is used to flag whether subclonal variant matches a clonal
/// variant of the same type called at the same position, which can flag the 
/// subclonal variant as untrustworthy.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct VariantLocation {
    pub ref_pos0:    SeqPos0,
    pub is_indel:    bool,
    pub re_fragment: ReFragment, // not haplotype here, we seek to compare across haplotypes
}
