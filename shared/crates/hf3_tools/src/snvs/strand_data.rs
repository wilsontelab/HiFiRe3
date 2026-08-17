//! Support for loading the original strand data, prior to duplex merging.
//! These original strands are used for haplotype consensus building and
//! subclonal variant calling.

// imports
use rustc_hash::FxHashMap;
use rust_htslib::bam::Record as BamRecord;
use super::*;

// constants
const QUAL_FLANK_BASES: usize = 3; // variant min_qual considers this many flanking bases

/// Required information for a given strand of a duplex read.
pub struct SourceStrand {
    pub seq:  Vec<BaseByteACGT>,
    pub qual: Vec<PhredQual>
}
impl SourceStrand {
    /// Create a new empty SourceStrand of a duplex read. This is replaced
    /// later in place from a BamRecord, so needs no capacity.
    pub fn new() -> Self {
        Self {
            seq:  Vec::new(),
            qual: Vec::new()
        }
    }

    /// Get the minimum base quality in the region surrounding a variant.
    pub fn get_min_qual(&self, mut start0: usize, mut end1: usize) -> u8 {
        start0 = start0.saturating_sub(QUAL_FLANK_BASES);
        end1 = (end1 + QUAL_FLANK_BASES).min(self.qual.len());
        *self.qual[start0..end1].iter().min().unwrap()
    }
}

/// A collection of SourceStrand pairs for all observed ReFragments.
pub struct SourceStrands {
    pub by_read: FxHashMap<QName, (SourceStrand, SourceStrand)>,
} 
impl SourceStrands {

    /// Create a new empty SourceStrands object for both strands of an as yet
    /// undetermined set of working reads. Reserve the minimal capacity they
    /// will are expected to require.
    pub fn new(worker: &SnvChromWorker, re_fragments: &Vec<ReFragment>) -> Self{
        let mut source_strands = Self {
            by_read: FxHashMap::default(),
        };
        source_strands.by_read.reserve(
            re_fragments.len() * worker.min_fragment_reads * 2
        );
        source_strands
    }

    /// Add a source strand sequence harvested from a --by-strand read files.
    /// Ignore reads not pre-filled into the by_read HashMap. Strand order 
    /// in the tuple of sequences is arbitrary; one strand sequence will
    /// be the reverse complement of the haplotype, which one is not yet known.
    pub fn insert(
        &mut self, 
        duplex_qname: QName,
        strand_aln: &BamRecord,
    ){
        if let Some(strands) = self.by_read.get_mut(&duplex_qname){
            if strands.0.seq.is_empty() {
                strands.0.seq  = strand_aln.seq().as_bytes();
                strands.0.qual = strand_aln.qual().to_vec();
            } else {
                strands.1.seq  = strand_aln.seq().as_bytes();
                strands.1.qual = strand_aln.qual().to_vec();
            }
        }
    }
}
