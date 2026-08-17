//! Support for collecting read data organized by their source RE fragments.

// imports
use std::str::from_utf8_unchecked;
use rustc_hash::FxHashMap;
use rust_htslib::bam::Record as BamRecord;
use serde::Serialize;
use genomex::bam::tags as bam_tags;
use crate::formats::hf3_tags::*;
use super::*;

/// A unique read span on a known chromosome corresponding to a RE fragment. For 
/// RE-based PacBio sequencing only a relatively limited number of unique read 
/// spans are expected.
/// 
/// The SEQ of each read assigned to a specific ReFragment is flush out to the 
/// exact start0 and end1 of the ReFragment, although the read's alignment may
/// start up to three bases inside of that position due to clipping allowed to 
/// account for RFLP SNPs. 
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Debug)]
pub struct ReFragment {
    pub start0: ChromPos0, // BED half-open coordinates (not site positions)
    pub end1:   ChromPos1, // start0 and end1 are oriented to ref top strand
}
impl ReFragment {

    /// Return the standardized ReFragment that matches a newly encountered read.
    pub fn from_aln(aln: &BamRecord) -> Option<Self> {
        let cigar = aln.cigar();
        let left_clip  = cigar.leading_softclips();
        let right_clip = cigar.trailing_softclips();
        if left_clip  > 3 ||
           right_clip > 3 {
            return None;
        }
        let re_fragment = ReFragment {
            start0: (aln.pos() - left_clip) as ChromPos0,
            end1: (cigar.end_pos() + right_clip) as ChromPos1
        };
        Some(re_fragment)
    }
}

/// A ReadInstance collects only the essential information from each encountered 
/// read for later re-analysis.
pub struct ReadInstance {
    pub sample_bit: SampleBit,
    pub qname:      QName,
    pub is_reverse: bool,
    pub aln_score:  u32,
    pub seq_bytes:  Vec<BaseByteACGTN>, // needed to find informative duplex bases
    pub cs: String, // as created by minimap2 during alignment
    pub dd: String, // as created by hf3_tools pre-alignment
    pub qry_pos0:   SeqPos0,   // where cs aln starts on SEQ, i.e., AFTER re-orientation to ref top strand
    pub aln_start0: ChromPos0, // where cs aln starts on ref top strand, not necessarily at site pos
}
impl ReadInstance {

    /// Create a new ReadInstance from its BamRecord.
    pub fn from_aln(aln: &BamRecord) -> Self {
        ReadInstance {
            sample_bit: bam_tags::get_tag_u32(aln, SAMPLE_BIT),
            qname:      unsafe { from_utf8_unchecked(aln.qname()).to_string() },
            is_reverse: aln.is_reverse(),
            aln_score:  bam_tags::get_tag_u32_default(aln, ALN_SCORE, 0), 

            // SAM SEQ is reference top-strand oriented
            seq_bytes:  aln.seq().as_bytes(), 

            // minimap2 cs tag is reference top-strand oriented
            cs: bam_tags::get_tag_str(aln, DIFFERENCE_STRING),  

            // dd tag is NOT reference top-strand oriented when is_reverse (done by get_dd_mask)
            dd: bam_tags::get_tag_str(aln, STRAND_DIFFERENCES),

            // qry_pos0, ref_start0, AND ref_end1 are reference top-strand oriented
            qry_pos0:   aln.cigar().leading_softclips() as u32, // thus, the 3' clip when is_reverse
            aln_start0: aln.pos() as ChromPos0
        }
    }
}

/// FragmentReads collects encountered ReadInstances for a given ReFragment.
/// One FragmentReads object is instantiated by each SnvChromWorker.
pub struct FragmentReads{
    pub instances: FxHashMap<ReFragment, Vec<ReadInstance>>
}
impl FragmentReads {

    /// Create a new FragmentReads HashMap.
    pub fn new() -> Self{
        let mut instances = FxHashMap::default();
        instances.reserve(4096); // about 16 Mb of ReFragments to start out
        Self{instances}
    }
    
    /// Add a ReadInstance to the FragmentReads HashMap. Reject reads that don't
    /// cleanly end within 3 bp of their nominated RE sites.
    pub fn insert(&mut self, aln: &BamRecord) -> usize {
        let Some(re_fragment) = ReFragment::from_aln(aln) 
            else { return 0; };
        let read_instance= ReadInstance::from_aln(aln);
        self.instances
            .entry(re_fragment)
            .or_insert_with(|| Vec::with_capacity(8))
            .push(read_instance);  
        1
    }
}

/// A ReadMapEntry carries bits of information a specific Variant in a specific
/// ReadInstance, including whether it was observed there and at what quality.
#[derive(Clone, Copy)]
pub struct ReadMapEntry {
    has_var: bool, // immutable record of whether a read reported a specific variant
    pub is_informative: bool, // false if the variant had N bases or bases error-corrected to reference
    pub min_qual: PhredQual,
}
impl ReadMapEntry{
    /// Create a new empty ReadMapEntry.
    pub fn new() -> Self{
        Self { 
            has_var: false, 
            is_informative: true,
            min_qual: 0 
        }
    }
    /// Get the immutable `has_var` value of a ReadMapEntry.
    pub fn has_var(&self) -> bool {
        self.has_var
    }
}

/// ReadMap collects information of the specific ReadInstances that reported a 
/// given Variant. Allocation is one-time fixed.
pub struct ReadMap {
    pub n_matching_reads: u16,
    pub n_informative: u16,
    pub zyg_int: u8, // from 0==heterozygous(0.5) to 100=fully homozygous
    pub read_map: Vec<ReadMapEntry>,
}
impl ReadMap {
    /// Create a new ReadMap.
    pub fn new(n_reads: usize) -> Self{
        Self{
            n_matching_reads: 0,
            n_informative: 0, // n_matching_reads + reads that could have called the variant
            zyg_int: 0,
            read_map: vec![ReadMapEntry::new(); n_reads],
        }
    } 
}

/// FragmentVariants collects encountered ReferenceVariants for a given 
/// ReFragment over all FragmentReads. One FragmentVariants object is 
/// instantiated per SnvChromWorker that is reset as needed per ReFragment.
pub struct FragmentVariants {
    pub n_reads: usize, // total reads assigned to the ReFragment
    pub variant_map: FxHashMap<Variant, ReadMap>,
}
impl FragmentVariants {

    /// Create a new FragmentVariants map.
    pub fn new() -> Self{
        let mut variant_map = FxHashMap::default();
        variant_map.reserve(128);
        Self{
            n_reads: 0,
            variant_map,
        }
    }

    /// Reset a FragmentVariants map to initialize collection of a new set of 
    /// ReFragment variants.
    pub fn reset(&mut self, n_reads: usize){
        self.n_reads = n_reads;
        self.variant_map.clear();
    }

    /// Add one Variant from a ReadInstance to its RefFragment's 
    /// FragmentVariants map.
    pub fn insert(
        &mut self, 
        variant:  Variant,
        read_i:   ReadIndex,
        min_qual: PhredQual,
    ) {
        let vmap = self.variant_map
            .entry(variant.clone())
            .or_insert_with(|| ReadMap::new(self.n_reads));
        vmap.n_matching_reads += 1;
        vmap.read_map[read_i] = ReadMapEntry{
            has_var:        true,
            is_informative: true,
            min_qual,
        };
    }
}
