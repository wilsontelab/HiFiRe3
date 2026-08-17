//! Handling of PacBio initial three-strand SNV/indel error correction tags and
//! the associated information carried in the dt, dd, and sk tags.
//! 
//! The logic describe below is essentially the same as used for final subclonal
//! variant calling, except that the POA approach used there forces all 
//! non-homoduplex variants to the haplotype consensus, i.e., no N bases are
//! needed since those three-strand consensus sequences aren't retained further. 
//! 
//! Initial three-strand error correction compares two PacBio strand consensuses 
//! (this, prev) to each other and to the reference genome (ref) to determine 
//! the final basecalling output where:
//! - homoduplex bases with Watson-Crick complementary strands are:
//!     - committed as sequenced, regardless of the reference (mis)match
//!     - committed with the higher base quality of the two strands
//!     - allowed to call clonal SNV and indel variants downstream
//! - heteroduplex bases where one strand matches reference are:
//!    - committed as the reference base
//!    - committed with the base quality of the strand that matched reference
//!    - not relevant to variant calling as they were error-corrected to reference
//!    - tracked with kinetics data for analyzing the reason for strand differences
//! - heteroduplex bases where neither strand matches reference, or there is no reference, are:
//!    - committed as one or more N bases in SEQ
//!        - for unresolved heteroduplex substitutions, a single N is reported
//!        - for unresolved heteroduplex indels, the N-track length matches the longer strand
//!             e.g., NN is reported if one strand reported two bases where the other reported none
//!    - committed with base quality 0 over all reported N bases
//!    - not allowed to call clonal SNV and indel variants downstream
//!
//! Two (this, prev) or three (this, prev, ref) strands are compared to determine 
//! the basecalling result. Coercion to reference and N/! bases are how duplex 
//! error correction is manifest, whereas homoduplex bases persist with the 
//! ability to call clonal variants (sub clonal variant calling return the 
//! the unmerged duplex strands).
//! 
//!       homoduplex            heteroduplex  
//!                        indel      substitution  
//! this  1  1  1  1     1  1  1  1    1  1  1  1     
//! prev  1  1  1  1     2  2  2  2    2  2  2  2  
//! ref   ?  1  2  ?     ?  1  2  3    ?  1  2  3  
//! read  1  1  1  1     N  1  2  N    N  1  2  N  
//! qual  X  X  X  X     !  1  2  !    !  1  2  !  X=strand max, 1|2=reported strand, !=0  
//! dd:Z  :  =  *  ^     !  +  -  #    ?  >  <  &  

// imports
use super::*;

// strand_merger outcome flag bits
pub const PERFECT_MATCH: u8                = 0;   // the strand sequences were exactly the same, and matched the reference perfectly
pub const HAS_STRAND_CLIP: u8              = 1;   // bases on prev_strand were clipped when aligning to this_strand
pub const HAS_REF_UNALIGNED: u8            = 2;   // bases on this_strand were not aligned to the reference genome, but others were
pub const HAS_HOMODUPLEX_INDEL: u8         = 4;   // the strands agreed on an indel relative to reference
pub const HAS_STRAND_INDEL: u8             = 8;   // the strands have an indel between them, regardless of its match to reference
pub const HAS_STRAND_INDEL_NEITHER_REF: u8 = 16;  // the strands have an indel between them that does not match reference on either strand
pub const HAS_HOMODUPLEX_SUBS: u8          = 32;  // the strands agreed on a base substitution relative to reference
pub const HAS_STRAND_SUBS: u8              = 64;  // the strands have a base substitution between them, regardless of its match to reference
pub const HAS_STRAND_SUBS_NEITHER_REF: u8  = 128; // the strands have a base substitution between them that does not match reference on either strand

// base to use in merged output SEQ when strands differ and neither strand matches the reference
// homoduplex strands always print bases as sequenced regardless of reference match
// heteroduplex strands print reference bases when one strand matches the reference
pub const SEQ_MASKED_BASE: char = 'N'; 

// DD tag operations
//   prev_on_this clip operations
pub const PREV_CLIP_OP: &str = "~";  // marks end clips of prev_on_this where prev.seq did not match this.seq bases

//   homoduplex operations
pub const HOMODUP_UNKNOWN:  &str = ":"; // homoduplex bases that could not be validated against a good reference alignment (none in the read were)
pub const HOMODUP_REF:      &str = "="; // homoduplex bases that DID     match the reference
pub const HOMODUP_NOT_REF:  &str = "*"; // homoduplex bases that DID NOT match the reference
pub const HOMODUP_NOT_ALN:  &str = "^"; // homoduplex bases that did not align to reference (although others in the read did)

//   heteroduplex indel operations
pub const HETERODUP_INDEL_UNKNOWN:     &str = "!"; // heteroduplex indels that could not be validated against a good reference alignment
pub const HETERODUP_INS_VS_REF:        &str = "+"; // heteroduplex indels that DID match reference on at least one strand
pub const HETERODUP_DEL_VS_REF:        &str = "-"; //    INS|DEL, i.e, +|- identifies the change on the non-reference strand relative to reference
pub const HETERODUP_INDEL_NEITHER_REF: &str = "#"; // heteroduplex indels that DID NOT match reference on either strand

//   heteroduplex base substitution operations
pub const HETERODUP_SUBS_UNKNOWN:      &str = "?"; // heteroduplex substitutions that could not be validated against a good reference alignment
pub const HETERODUP_SUBS_THIS_REF:     &str = ">"; // this.seq base (listed first  in the op value) matched the reference base
pub const HETERODUP_SUBS_PREV_REF:     &str = "<"; // prev.seq base (listed second in the op value) matched the reference base
pub const HETERODUP_SUBS_NEITHER_REF:  &str = "&"; // heteroduplex substitutions that DID NOT match reference on either strand 

/// The types of values found in a `dd:Z:` tag mask.
#[repr(u8)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum DdMaskType {
    Homoduplex = 0,
    EndClipped = 1,
    UnresolvedHeteroduplex = 2,
    CorrectedToReference = 3,
}
impl SnvChromWorker {

    /// Convert a `dd:Z:` tag into a `Vec<DdMaskType>` indicating whether each 
    /// read position was either:
    ///     - error corrected to reference 
    ///     - reported as an N base at an unresolved heteroduplex
    /// 
    /// during basecalling. Homoduplex bases can always call variants. 
    /// EndClipped and UnresolvedHeteroduplex bases will never call variants  
    /// since they were reported as N. CorrectedToReference bases will be 
    /// reported in lower case to guard against the rare situation where the 
    /// reference base was the incorrect choice, i.e., when heteroduplex bases 
    /// were encountered at a true SNP. 
    /// 
    /// For reverse strand alignments, the read mask is reversed to match the read 
    /// SEQ order in the BamRecord.
    pub fn get_dd_mask(
        read: &ReadInstance,
    ) -> Vec<DdMaskType> {
        let read_len = read.seq_bytes.len();
        let mut mask: Vec<DdMaskType> = vec![DdMaskType::Homoduplex; read_len];
        let mut offset0: usize = 0;
        let mut chars = read.dd.chars();
        let mut op = chars.next().unwrap();
        let mut val: String = String::with_capacity(128);
        while let Some(char) = chars.next() {
            if char.is_alphanumeric() {
                val.push(char);
            } else {
                Self::add_to_dd_mask(&mut mask, &mut offset0, op, &val);
                op = char;
                val.clear();
            }
        }
        Self::add_to_dd_mask(&mut mask, &mut offset0, op, &val);
        if read.is_reverse { mask.reverse(); }
        mask
    }

    /// Add one dd tag operation to the read mask.
    fn add_to_dd_mask(
        mask:    &mut Vec<DdMaskType>, 
        offset0: &mut usize,
        op:      char, 
        val:     &str, 
    ) {
        match op {
            // prev_on_this clip operations
            //      two-strand validation of a reference variant is impossible
            '~' => Self::set_dd_mask(mask, offset0, DdMaskType::EndClipped, val.parse::<usize>().unwrap()),
            // homoduplex operations
            //      always allowed to call variants, but only * operations are expected to do so
            //      as alignment outcomes will presumably continue to be the same
            ':' => Self::set_dd_mask(mask, offset0, DdMaskType::Homoduplex, val.parse::<usize>().unwrap()),
            '=' => Self::set_dd_mask(mask, offset0, DdMaskType::Homoduplex, val.parse::<usize>().unwrap()),
            '*' => Self::set_dd_mask(mask, offset0, DdMaskType::Homoduplex, 1), // always come one read base at a time
            '^' => Self::set_dd_mask(mask, offset0, DdMaskType::Homoduplex, val.parse::<usize>().unwrap()),
            // heteroduplex indel operations
            //       ! and # never allowed to call variants since they weren't validated by both read strands
            //       + and - are never expected to call variants as they were error-corrected to reference
            '!' => Self::set_dd_mask(mask, offset0, DdMaskType::UnresolvedHeteroduplex, val.len()), // unknown read and op have same length
            '+' => Self::set_dd_mask(mask, offset0, DdMaskType::CorrectedToReference, 0),      // heterodup insertions relative to ref not included in read
            '-' => Self::set_dd_mask(mask, offset0, DdMaskType::CorrectedToReference,   val.len()), // whereas deletions were committed as ref bases
            '#' => Self::set_dd_mask(mask, offset0, DdMaskType::UnresolvedHeteroduplex, val.len()),
            // heteroduplex base substitution operations
            //       ? and & never allowed to call variants since they weren't validated by both read strands
            //       > and < are never expected to call variants as they were error-corrected to reference
            '?' => Self::set_dd_mask(mask, offset0, DdMaskType::UnresolvedHeteroduplex, 1),
            '>' => Self::set_dd_mask(mask, offset0, DdMaskType::CorrectedToReference,   1),
            '<' => Self::set_dd_mask(mask, offset0, DdMaskType::CorrectedToReference,   1),
            '&' => Self::set_dd_mask(mask, offset0, DdMaskType::UnresolvedHeteroduplex, 1),
            _   => panic!("Unexpected operation in dd tag: {}", op),
        }
    }

    /// Update a block of contiguous read positions in the mask to false as 
    /// needed and increment the position offset.
    fn set_dd_mask(
        mask:    &mut Vec<DdMaskType>,
        offset0: &mut usize, 
        mask_type: DdMaskType, 
        len: usize
    ) {
        if mask_type != DdMaskType::Homoduplex {
            if len > 0 {
                for i in *offset0..(*offset0 + len) {
                    mask[i] = mask_type;
                }
            } else if *offset0 > 0 {
                // heteroduplex insertions relative to ref were corrected to ref
                // thus have no corresponding query base, record mask on previous query base
                mask[*offset0 - 1] = mask_type;
            }
        }
        *offset0 += len;
    }
}
