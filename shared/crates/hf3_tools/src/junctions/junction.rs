//! Structure and methods for aggregating and describing detected junctions.

// dependencies
use std::mem::take;
use genomex::sam::junction::JunctionType;
use rustc_hash::FxHashMap;
use rust_htslib::bam::record::Record as BamRecord;
use rayon::prelude::*;
use mdi::OutputCsv;
use genomex::bam::{tags, cigar};
use genomex::sam::{junction::Junction, SamRecord};
use crate::formats::hf3_tags::*;
use crate::inserts::ReadLevelMetadata;
use super::{JunctionAnalysisTool, JxnFailureFlag};

// constants
const INSTANCE_CAPACITY: usize = 100; // pre-allocate space for junction instances
const NULL_STEM_LENGTH: i32 = 999999 as i32; // when STEM_LENGTH3 is not available

/// An OrderedJunction describes a single detected junction in the canonical 
/// orientation, by orienting a deserialized genomex::Junction JUNCTION tag 
/// using this module's extract_junction().
/// 
/// Values are the same for all instances of a given junction and used to 
/// key junction aggregation, where junction "sameness" is based on the 
/// breakpoint nodes and the number of bases in the junction microhomology 
/// or insertion, i.e, its alignment offset. Thus, an OrderedJunction 
/// considers the number of bases in the junction path, but allows for 
/// base differences resulting from basecalling  errors (but not indels, 
/// which are handled during fuzzy matching).
#[derive(PartialEq, Eq, Hash, Clone, Copy)]
pub struct OrderedJunction {
    /* ------------------------------------------- */
    // two ordered breakpoint nodes and an alignment offset define a junction
    // each of these values is needed for fuzzy matching and grouping of junctions
    pub node_aln5_end3: isize, // +|- chrom_index << 29 | junction node pos1
    pub node_aln3_end5: isize, // node_aln3_end5 = 5' end of the 3'-most alignment in the canonical orientation
    pub offset:         i16,   // (+)offset = insertion, (-)offset = microhomology overlap
    /* ------------------------------------------- */
    // junction types and sizes derive deterministically from nodes 
    // so are the same for all JunctionInstances
    pub jxn_type: u8,
    strands:      u8,
    sv_size:      u32, // modal size selected when OrderedJunctions are aggregated by fuzzy matching
    // other node-level properties that aren't part of a Junction are added later
}

/// A JunctionInstance is one observation of a given OrderedJunction. 
/// These values derive from one read's alignments as parsed by this 
/// module's extract_junction().
/// 
/// Usually there is one instance of a given OrderedJunction per read, 
/// but multiple instances could occur if the read has multiple 
/// identical alignments, e.g., in tandem repeat regions. Fields `aln5_i``
/// and `qry_pos1_aln5_end3` distinguish multiple instances of the same
/// OrderedJunction in one read by recording the position of the
/// junction within the original read.
pub struct JunctionInstance {
    /* ------------------------------------------- */
    // junction properties that collapse to a single best value after fuzzy matching
    /* ------------------------------------------- */
    jxn_seq: String, // bases in the canonical orientation; not used for junction grouping
    /* ------------------------------------------- */
    // values maintained as lists of values per junction instance
    /* ------------------------------------------- */
    // values that are the same for all read alignments (even those not flanking the junction)
    // the same read_data may be present twice if the read has multiple identical junctions
    // different JunctionInstances from the same read will either all be duplicates or not
    pub read_data: ReadLevelMetadata, // qname, insert_size, node1, node2, channel, n_jxns, duplicate
    /* ------------------------------------------- */
    // values that place (multiple) occurrences of the (same) junction in one read
    // these values will differ for two instances within the same read_data
    // these offsets reflect the original, not the reoriented, read sequence
    aln5_i:             u8,  // index of the 5' alignment in the read's alignments
    qry_pos1_aln5_end3: u32, // number of bases from read start to the 3' end of the 5' alignment
    /* ------------------------------------------- */
    // values that combine properties of the two flanking alignments
    jxn_orientation:  u8, // 0|1 bool if this junction instance was flipped relative the original read orientation
    jxn_failure_flag: u8, // JXN_FAILURE_FLAG tag value
    aln_failure_flag: u8, // bitwise OR of the ALN_FAILURE_FLAG tag values on each side of the junction
    /* ------------------------------------------- */
    // values that select one value from the two flanking alignments
    min_stem_length: u32, // BEST (STEM_LENGTH5.abs(), STEM_LENGTH3.abs()) of a flanking alignment (only one has to pass)
    min_mapq:        u8,  // WORST mapping quality of a flanking alignment (a junction is only as good as its worst flank)
    max_divergence:  f32, // WORST DIVERGENCE tag value of a flanking alignment
}

/// A JunctionInstances collection holds individual observations of a 
/// given OrderedJunction. These values might differ between instances 
/// of the same junction, including jxn_seq, which depends on the 
/// sequenced insert bases, not just the breakpoint nodes and offset.
#[derive(Clone)]
pub struct JunctionInstances {
    // see JunctionInstance for extended field descriptions
    /* ------------------------------------------- */
    pub jxn_seqs: Vec<String>, // bases in the canonical orientation; not used for junction grouping
    /* ------------------------------------------- */
    pub read_data: Vec<ReadLevelMetadata>, // qname, insert_size, node1, node2, channel, n_jxns, duplicate
    /* ------------------------------------------- */
    aln5_is: Vec<u8>,              // index of the 5' alignment in the read's alignments
    qry_pos1_aln5_end3s: Vec<u32>, // number of bases from read start to the 3' end of the 5' alignment
    /* ------------------------------------------- */
    pub jxn_orientations:  Vec<u8>, // 0|1 bool if this junction instance was flipped relative the original read orientation
    jxn_failure_flags: Vec<u8>,     // JXN_FAILURE_FLAG tag value
    aln_failure_flags: Vec<u8>,     // bitwise OR of the ALN_FAILURE_FLAG tag value on each side of the junction
    /* ------------------------------------------- */
    min_stem_lengths: Vec<u32>, // BEST (STEM_LENGTH5.abs(), STEM_LENGTH3.abs()) of a flanking alignment (only one has to pass)
    min_mapqs:        Vec<u8>,  // WORST mapping quality of a flanking alignment (a junction is only as good as its worst flank)
    max_divergences:  Vec<f32>, // WORST DIVERGENCE tag value of a flanking alignment
}
impl JunctionInstances {
    /// Create a new empty set of JunctionInstances. 
    pub fn new() -> Self {
        JunctionInstances {
            /* ------------------------------------------- */
            jxn_seqs:            Vec::with_capacity(INSTANCE_CAPACITY), // pre-allocate more than usually needed
            /* ------------------------------------------- */
            read_data:           Vec::with_capacity(INSTANCE_CAPACITY),
            /* ------------------------------------------- */
            aln5_is:             Vec::with_capacity(INSTANCE_CAPACITY),
            qry_pos1_aln5_end3s: Vec::with_capacity(INSTANCE_CAPACITY),
            /* ------------------------------------------- */
            jxn_orientations:    Vec::with_capacity(INSTANCE_CAPACITY),
            jxn_failure_flags:   Vec::with_capacity(INSTANCE_CAPACITY),
            aln_failure_flags:   Vec::with_capacity(INSTANCE_CAPACITY),
            /* ------------------------------------------- */
            min_stem_lengths:    Vec::with_capacity(INSTANCE_CAPACITY),
            min_mapqs:           Vec::with_capacity(INSTANCE_CAPACITY),
            max_divergences:     Vec::with_capacity(INSTANCE_CAPACITY),
        }
    }
    /// Add a JunctionInstance to the collection.
    pub fn add_instance(&mut self, instance: JunctionInstance) {
        /* ------------------------------------------- */
        self.jxn_seqs           .push(instance.jxn_seq); // pass owned values
        /* ------------------------------------------- */
        self.read_data          .push(instance.read_data);
        /* ------------------------------------------- */
        self.aln5_is            .push(instance.aln5_i);
        self.qry_pos1_aln5_end3s.push(instance.qry_pos1_aln5_end3);
        /* ------------------------------------------- */
        self.jxn_orientations   .push(instance.jxn_orientation);
        self.jxn_failure_flags  .push(instance.jxn_failure_flag);
        self.aln_failure_flags  .push(instance.aln_failure_flag);
        /* ------------------------------------------- */
        self.min_stem_lengths   .push(instance.min_stem_length);
        self.min_mapqs          .push(instance.min_mapq);
        self.max_divergences    .push(instance.max_divergence);
    } 
    /// Add multiple JunctionInstances to the collection.
    pub fn extend(&mut self, other: JunctionInstances) {
        /* ------------------------------------------- */
        self.jxn_seqs           .extend(other.jxn_seqs); // take ownership of Vecs
        /* ------------------------------------------- */
        self.read_data          .extend(other.read_data);
        /* ------------------------------------------- */
        self.aln5_is            .extend(other.aln5_is);
        self.qry_pos1_aln5_end3s.extend(other.qry_pos1_aln5_end3s);
        /* ------------------------------------------- */
        self.jxn_orientations   .extend(other.jxn_orientations);
        self.jxn_failure_flags  .extend(other.jxn_failure_flags);
        self.aln_failure_flags  .extend(other.aln_failure_flags);
        /* ------------------------------------------- */
        self.min_stem_lengths   .extend(other.min_stem_lengths);
        self.min_mapqs          .extend(other.min_mapqs);
        self.max_divergences    .extend(other.max_divergences);
    }

    /// Get the best, most frequent jxn_seq from a set of JunctionInstances,
    /// preferring to use non-duplicate reads if possible.
    fn get_best_jxn_seq(&mut self, best_offset_abs: usize) -> String {
        if self.jxn_seqs.len() == 1 {
            return take(&mut self.jxn_seqs[0]); // a singleton junction, safe to take as is
        }
        let mut counts = self.jxn_seqs.iter()
            .zip(&self.read_data)
            .filter(|(seq, rd)| seq.len() == best_offset_abs && !rd.is_duplicate)
            .map(|(seq, _)| seq)
            .fold(FxHashMap::default(), |mut acc, seq| {
                *acc.entry(seq).or_insert(0_u16) += 1;
                acc
            });
        if counts.is_empty() { // handle edge case where all reads with best_offset_abs are duplicates
            counts = self.jxn_seqs.iter()
                .filter(|seq| seq.len() == best_offset_abs)
                .fold(FxHashMap::default(), |mut acc, seq| {
                    *acc.entry(seq).or_insert(0_u16) += 1;
                    acc
                });
        }
        if counts.is_empty() { // handle extreme edge case where no reads have the best_offset_abs length (?)
            counts = self.jxn_seqs.iter()
                .fold(FxHashMap::default(), |mut acc, seq| {
                    *acc.entry(seq).or_insert(0_u16) += 1;
                    acc
                });
        }
        counts.into_iter() // guaranteed to have at least one entry here
            .max_by_key(|&(_, count)| count)
            .map(|(seq, _)| seq.clone())
            .unwrap()
    }

    /// Get the number of junction instances,
    /// optionally excluding those marked as duplicates.
    fn n_instances(&self, allow_dups: bool) -> u16 {
        if !allow_dups {
            self.read_data.iter().filter(|rd| !rd.is_duplicate).count() as u16
        } else {
            self.read_data.len() as u16
        }
    }

    /// Get the number of unique reads that generated all junction instances,
    /// optionally excluding those marked as duplicates.
    fn n_unique_reads(&self, allow_dups: bool) -> u16 {
        let iter = self.read_data.iter();
        let mut qnames: FxHashMap<&str, usize> = FxHashMap::default();
        if !allow_dups {
            iter.filter_map(|rd| if rd.is_duplicate { None } else { Some(&rd.qname) } )
                .for_each(|qname| { *qnames.entry(qname).or_insert(0) += 1; });
        } else {
            iter.map(|rd| &rd.qname)
                .for_each(|qname| { *qnames.entry(qname).or_insert(0) += 1; });
        };
        qnames.len() as u16
    }
}

/// Extract the OrderedJunction and JunctionInstance from the two alignments flanking a junction,
/// where aln5 already carries the JUNCTION tag that aggregates junction metadata.
pub fn extract_junction(
    aln5_i:    u8,
    aln5:      &BamRecord,
    aln3:      &BamRecord,
    read_data: &ReadLevelMetadata,
    qlen:      u32,
) -> (OrderedJunction, JunctionInstance) { 
    // the 5' alignment carries the encoded junction metadata
    let jxn_tag = tags::get_tag_str(aln5, JUNCTION); // tag must exist if this function is called
    let jxn = Junction::deserialize(&jxn_tag).orient_junction();
    // only consider meaningful 3' stem lengths
    let stem_length3 = tags::get_tag_i32_default(aln3, STEM_LENGTH3, NULL_STEM_LENGTH).unsigned_abs(); 
    let stem_length3_wrk = if stem_length3 > 1 { stem_length3 } else { NULL_STEM_LENGTH as u32 }; 
    (
        OrderedJunction {
            /* ------------------------------------------- */
            node_aln5_end3: jxn.node_aln5_end3,
            node_aln3_end5: jxn.node_aln3_end5,
            offset:         jxn.offset,
            /* ------------------------------------------- */
            jxn_type:       jxn.jxn_type as u8,
            strands:        jxn.strands  as u8,
            sv_size:        jxn.sv_size  as u32,
        },
        JunctionInstance {
            /* ------------------------------------------- */
            jxn_seq:   jxn.jxn_seq,
            /* ------------------------------------------- */
            read_data: read_data.clone(),
            /* ------------------------------------------- */
            aln5_i:             aln5_i,
            qry_pos1_aln5_end3: cigar::get_query_end1(&aln5, qlen),
            /* ------------------------------------------- */
            jxn_orientation:  jxn.was_flipped as u8,
            jxn_failure_flag: tags::get_tag_u8_default(aln5, JXN_FAILURE_FLAG, 0),
            aln_failure_flag: (
                tags::get_tag_u8_default(aln5, ALN_FAILURE_FLAG, 0) | 
                tags::get_tag_u8_default(aln3, ALN_FAILURE_FLAG, 0)
            ),
            /* ------------------------------------------- */
            min_stem_length:  tags::get_tag_i32_default(aln5, STEM_LENGTH5, 1).unsigned_abs().min(stem_length3_wrk),
            min_mapq:         aln5.mapq().min(aln3.mapq()),
            max_divergence:   tags::get_tag_f32_default(aln5, DIVERGENCE, 0.0).max(
                              tags::get_tag_f32_default(aln3, DIVERGENCE, 0.0)),
        }
    )
}

/// A FinalJunction describes one fully grouped, deduplicated and 
/// resolved junction call. Some fields are initially printed 
/// with null values to be updated later.
#[derive(serde::Serialize, serde::Deserialize, Clone, Debug)]
pub struct FinalJunction {
    /* ------------------------------------------- */
    // two best ordered breakpoint nodes and an alignment offset define a junction
    #[serde(rename = "#chrom_index1_1")]
    pub chrom_index1_1:  u8,  // node1, broken into its constituent parts for chrom-pos sorting
    pub ref_pos1_1:      u32,
    pub strand_index0_1:     u8,
    pub chrom_index1_2:  u8,  // node2
    pub ref_pos1_2:      u32,
    pub strand_index0_2:     u8,
    pub offset:  i16,    // most frequent junction offset in a fuzzy-matched group
    pub jxn_seq:     String, // most frequent junction bases among instances with the best junction offset
    /* ------------------------------------------- */
    // junction properties that derive deterministically from nodes
    pub jxn_type: u8,
    strands:      u8,
    pub sv_size:      u32,
    /* ------------------------------------------- */
    // values maintained as lists of values per junction instance
    /* ------------------------------------------- */
    // values expanded and updated from read_data (ReadLevelMetadata); includes all instances, even duplicates
    pub qnames:        String, // comma-delimited of read QNAMEs that sequenced this junction; may have duplicates if a read had multiple identical junctions
    pub insert_sizes:  String, // comma-delimited insert sizes of qnames
    outer_node1s:  String, // comma-delimited outer node1 values of qnames
    outer_node2s:  String, // comma-delimited outer node2 values of qnames
    channels:      String, // comma-delimited ONT channel numbers of qnames; 0 if not applicable
    n_jxns:        String, // comma-delimited number of junctions in qnames
    is_duplicates: String, // comma-delimited integer bools if each QNAME was marked as a duplicate read; always duplicated or not for identical qnames
    is_duplexes:   String, // comma-delimited integer bools if duplex basecalling was performed on each read
    /* ------------------------------------------- */
    // values that place (multiple) occurrences of the (same) junction in one read
    aln5_is:             String, // comma-delimited index of the 5' alignment in a QNAME's alignments
    qry_pos1_aln5_end3s: String, // comma-delimited number of bases from read start to the 3' end of the 5' alignment
    /* ------------------------------------------- */
    // values that combine properties of the two flanking alignments
    jxn_orientations:  String, // comma-delimited 0|1 bool if this junction instance was flipped relative the original read orientation
    jxn_failure_flags: String, // comma-delimited JXN_FAILURE_FLAG tag value
    aln_failure_flags: String, // comma-delimited bitwise OR of the ALN_FAILURE_FLAG tag value on each side of the junction
    /* ------------------------------------------- */
    // values that select one value from the two flanking alignments
    min_stem_lengths:  String, // comma-delimited BEST (STEM_LENGTH5.abs(), STEM_LENGTH3.abs()) of a flanking alignment (only one has to pass)
    min_mapqs:         String, // comma-delimited WORST mapping quality of a flanking alignment (a junction is only as good as its worst flank)
    max_divergences:   String, // comma-delimited WORST DIVERGENCE tag value of a flanking alignment
    /* ------------------------------------------- */
    // counts of supporting instances and reads
    n_instances:       u16, // number of junction instances that sequenced this FinalJunction (the largest of all the counts)
    n_reads:           u16, // number of unique reads that sequenced this FinalJunction (can be < n_instances)
    pub n_instances_dedup: u16, // number of deduplicated instances after read-level deduplication
    n_reads_dedup:     u16, // number of unique reads after read-level deduplication (the smallest of all the counts)
    /* ------------------------------------------- */
    // values that reflect aggregated quality metrics over all deduplicated read instances
    has_multi_jxn_read:      u8,  // integer bool; true if any QNAME had multiple junctions of any type
    has_multi_instance_read: u8,  // integer bool; true if any QNAME had multiple instances of this junction
    has_duplex_read:         u8,  // integer bool; true if any read with this junction had duplex basecalling performed
    is_bidirectional:        u8,  // integer bool; true if junction instances were observed in both orientations
    pub jxn_failure_flag:        u8,  // bitwise OR of all instance jxn_failure_flags
    pub aln_failure_flag:        u8,  // bitwise OR of all instance aln_failure_flags
    pub any_was_chimeric:    u8,  // integer bool; true if any jxn_failure_flag had a chimeric bit set
    pub any_was_not_chimeric:    u8,  // integer bool; true if any jxn_failure_flag did not have a chimeric bit set
    pub min_stem_length:         u32, // best stem length over all instances; thus, whether any instance was a high-quality detection
    pub max_min_mapq:            u8,  // best min_mapq over all instances
    pub min_max_divergence:      f32, // best max_divergence over all instances
    /* ------------------------------------------- */
    // additional junction properties that derive deterministically from nodes
    pub is_intergenomic: u8,     // integer bool; true if the junction connects different genomes of a composite reference
    target_1:       String, // name of the target region containing breakpoint 1; * if untargeted or not on target
    target_dist_1:  i32,    // signed target distance as ref_pos1_1 - target_1 center; 0 if untargeted or not on target
    target_2:       String, // name of target region containing breakpoint 2; * if untargeted or not on target
    target_dist_2:  i32,    // signed target distance as ref_pos1_2 - target_2 center; 0 if untargeted or not on target
    genes_1:        String, // comma-delimited gene(s) overlapping breakpoint 1; * if not in a gene
    gene_dists_1:   String, // comma-delimited signed gene distances as ref_pos1_1 - gene_1 center; 0 if not in a gene
    genes_2:        String, // comma-delimited gene(s) overlapping breakpoint 2; * if not in a gene
    gene_dists_2:   String, // comma-delimited signed gene distances as ref_pos1_2 - gene_2 center; 0 if not in a gene
    is_excluded_1:  u8,     // integer bool; true if breakpoint 1 overlaps an excluded region
    is_excluded_2:  u8,     // integer bool; true if breakpoint 2 overlaps an excluded region
    /* ------------------------------------------- */
    // additional read-level properties that derive deterministically from the global config
    pub sample_bits:  u32, // bitwise OR of the sample bits of all samples that sequenced this SV
    n_samples:    u8,  // number of unique samples that sequenced this SV
    /* ------------------------------------------- */
    // additional junction properties added later, initialized to zero
    pub bkpt_coverage_1: u16, // total number of any type of alignment that crossed breakpoint 1
    pub bkpt_coverage_2: u16, // total number of any type of alignment that crossed breakpoint 2
}
impl FinalJunction {

    /// Convert a single best OrderedJunction and its JunctionInstances into a FinalJunction.
    /// Not all `instances` necessarily match exactly to `best_jxn` if they were fuzzy grouped.
    pub fn from_best_junction(
        best_jxn:      &OrderedJunction,
        mut instances: JunctionInstances, // possibly a concatenation of JunctionInstances from multiple fuzzy-matched OrderedJunctions
        tool:          &JunctionAnalysisTool,
        ont_cache:     &mut FxHashMap<u32, FxHashMap<u8, Vec<usize>>>, // pre-allocated deduplication caches
        other_cache:   &mut FxHashMap<(isize, isize), Vec<usize>>,
    ) -> Self {

        // unpack the breakpoint nodes from the best_jxn
        let (chrom_1, chrom_index1_1, ref_pos1_1, is_reverse_1) = 
            SamRecord::unpack_signed_node(best_jxn.node_aln5_end3, &tool.chroms);
        let (chrom_2, chrom_index1_2, ref_pos1_2, is_reverse_2) = 
            SamRecord::unpack_signed_node(best_jxn.node_aln3_end5, &tool.chroms);

        // determine the best offset length for jxn_seq extraction
        // this unsigned value is the length of either an insertion or a microhomology span
        let best_offset_abs = best_jxn.offset.unsigned_abs() as usize;

        // mark duplicates as needed, preferring reads that match the best_jxn offset when possible
        let n_instances = instances.n_instances(true);
        if n_instances > 1 {
            if tool.is_ont {
                Self::mark_duplicates_by_channel(&mut instances, ont_cache);
            } else if tool.deduplicate_reads {
                Self::mark_duplicates_by_span(   &mut instances, other_cache, best_offset_abs);
            }
        }

        // get the most frequently observed (non-duplicate) value; must match the best_jxn offset
        let best_jxn_seq = instances.get_best_jxn_seq(best_offset_abs);

        // get intersections with genome regions
        let (target_1, target_dist_1) = tool.targets.get_pos_region_unpadded(&chrom_1, ref_pos1_1);
        let (target_2, target_dist_2) = tool.targets.get_pos_region_unpadded(&chrom_2, ref_pos1_2);
        let (genes_1,  gene_dists_1)  = tool.genes  .get_pos_region_unpadded(&chrom_1, ref_pos1_1);
        let (genes_2,  gene_dists_2)  = tool.genes  .get_pos_region_unpadded(&chrom_2, ref_pos1_2);

        // construct and return the FinalJunction
        let rd = &instances.read_data;
        let sample_bits = rd.iter().fold(0_u32, |acc, rd| acc | rd.sample_bit);
        FinalJunction {
            /* ------------------------------------------- */
            chrom_index1_1,
            ref_pos1_1,
            strand_index0_1: is_reverse_1 as u8,
            chrom_index1_2,
            ref_pos1_2,
            strand_index0_2: is_reverse_2 as u8,
            offset:          best_jxn.offset, // NOT best_offset_abs, we want the signed value here
            jxn_seq:         best_jxn_seq,
            /* ------------------------------------------- */
            jxn_type:        best_jxn.jxn_type,
            strands:         best_jxn.strands,
            sv_size:         best_jxn.sv_size,
            /* ------------------------------------------- */
            qnames:        Self::join_with_flanked_commas(
                              &rd.iter().map(|rd| rd.qname.clone()).collect::<Vec<_>>(),
                              100),
            insert_sizes:  Self::join_with_flanked_commas(
                              &rd.iter().map(|rd| rd.insert_size).collect::<Vec<_>>(), 
                              6),
            outer_node1s:  Self::join_with_flanked_commas(
                              &rd.iter().map(|rd| rd.node1).collect::<Vec<_>>(), 
                              14),
            outer_node2s:  Self::join_with_flanked_commas(
                              &rd.iter().map(|rd| rd.node2).collect::<Vec<_>>(), 
                              14),
            channels:      Self::join_with_flanked_commas(
                              &rd.iter().map(|rd| rd.channel).collect::<Vec<_>>(),  
                              6),
            n_jxns:        Self::join_with_flanked_commas(
                              &rd.iter().map(|rd| rd.n_jxns).collect::<Vec<_>>(),     
                              2),
            is_duplicates: Self::join_with_flanked_commas(
                              &rd.iter().map(|rd| rd.is_duplicate as u8).collect::<Vec<_>>(), 
                              1),
            is_duplexes:   Self::join_with_flanked_commas(
                              &rd.iter().map(|rd| rd.is_duplex as u8).collect::<Vec<_>>(), 
                              1),
            /* ------------------------------------------- */
            aln5_is:             Self::join_with_flanked_commas(&instances.aln5_is, 2),
            qry_pos1_aln5_end3s: Self::join_with_flanked_commas(&instances.qry_pos1_aln5_end3s, 6),
            /* ------------------------------------------- */
            jxn_orientations:  Self::join_with_flanked_commas(&instances.jxn_orientations, 1),
            jxn_failure_flags: Self::join_with_flanked_commas(&instances.jxn_failure_flags, 3),
            aln_failure_flags: Self::join_with_flanked_commas(&instances.aln_failure_flags, 3),
            /* ------------------------------------------- */
            min_stem_lengths:  Self::join_with_flanked_commas(&instances.min_stem_lengths, 5),
            min_mapqs:         Self::join_with_flanked_commas(&instances.min_mapqs, 2),
            max_divergences:   Self::join_with_flanked_commas(&instances.max_divergences, 10),
            /* ------------------------------------------- */
            n_instances,
            n_reads:           instances.n_unique_reads(true),
            n_instances_dedup: instances.n_instances(   false),
            n_reads_dedup:     instances.n_unique_reads(false),
            /* ------------------------------------------- */
            has_multi_jxn_read: rd.iter().any(|rd| rd.n_jxns > 1) as u8,
            has_multi_instance_read: {
                let mut qname_counts: FxHashMap<&str, usize> = FxHashMap::default();
                for rd in rd {
                    *qname_counts.entry(&rd.qname).or_insert(0) += 1;
                }
                qname_counts.values().any(|&count| count > 1) as u8
            },
            has_duplex_read:      rd.iter().any(|rd| rd.is_duplex) as u8,
            is_bidirectional:    (instances.jxn_orientations.iter().any(|&ori| ori == 0) &&
                                  instances.jxn_orientations.iter().any(|&ori| ori == 1)) as u8,
            jxn_failure_flag:     instances.jxn_failure_flags.iter().fold(0u8, |acc, &flag| acc | flag),
            aln_failure_flag:     instances.aln_failure_flags.iter().fold(0u8, |acc, &flag| acc | flag),
            any_was_chimeric:     instances.jxn_failure_flags.iter().any(|&flag|  JxnFailureFlag::is_chimeric_jxn_u8(flag)) as u8,
            any_was_not_chimeric: instances.jxn_failure_flags.iter().any(|&flag| !JxnFailureFlag::is_chimeric_jxn_u8(flag)) as u8,
            min_stem_length:     *instances.min_stem_lengths.iter().min().unwrap(),
            max_min_mapq:        *instances.min_mapqs.iter().max().unwrap(),
            min_max_divergence:  *instances.max_divergences.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap(),
            /* ------------------------------------------- */
            is_intergenomic: !tool.chroms.is_same_genome_suffix(&chrom_1, &chrom_2) as u8,
            target_1:       target_1.names[0].clone(), // targets are non-overlapping, only ever a single value expected
            target_dist_1:  target_dist_1[0],
            target_2:       target_2.names[0].clone(),
            target_dist_2:  target_dist_2[0],
            genes_1:        genes_1.names.join(","), // genes may overlap each other, so join with commas
            gene_dists_1:   gene_dists_1.iter().map(|&d| d.to_string()).collect::<Vec<_>>().join(","),
            genes_2:        genes_2.names.join(","),
            gene_dists_2:   gene_dists_2.iter().map(|&d| d.to_string()).collect::<Vec<_>>().join(","),
            is_excluded_1:  tool.exclusions.pos_in_region(&chrom_1, ref_pos1_1) as u8,
            is_excluded_2:  tool.exclusions.pos_in_region(&chrom_2, ref_pos1_2) as u8,
            /* ------------------------------------------- */
            // sample_names: Self::join_with_flanked_commas(&vec![tool.data_name.clone()], tool.data_name.len()),
            sample_bits: sample_bits,
            n_samples:   sample_bits.count_ones() as u8,
            /* ------------------------------------------- */
            bkpt_coverage_1: 0, // added downstream after alignment segments are merge sorted
            bkpt_coverage_2: 0,
        }
    }

    /// Convert two or more merging FinalJunctions into a single output FinalJunction.
    pub fn from_merged_junctions(jxns: &[FinalJunction]) -> Self {
        let sample_bits = jxns.iter().fold(0_u32, |acc, jxn| acc | jxn.sample_bits); // previously bit-shifted by source
        let best_jxn  = jxns.iter().max_by_key(|jxn| jxn.n_instances_dedup).unwrap();
        FinalJunction {
            /* ------------------------------------------- */
            chrom_index1_1:  best_jxn.chrom_index1_1,
            ref_pos1_1:      best_jxn.ref_pos1_1,
            strand_index0_1: best_jxn.strand_index0_1,
            chrom_index1_2:  best_jxn.chrom_index1_2,
            ref_pos1_2:      best_jxn.ref_pos1_2,
            strand_index0_2: best_jxn.strand_index0_2,
            offset:          best_jxn.offset,
            jxn_seq:         best_jxn.jxn_seq.clone(),
            /* ------------------------------------------- */
            jxn_type:        best_jxn.jxn_type,
            strands:         best_jxn.strands,
            sv_size:         best_jxn.sv_size,
            /* ------------------------------------------- */
            qnames:        Self::join_with_flanked_commas(
                              &jxns.iter().map(|jxn| jxn.qnames[1..jxn.qnames.len() - 1].to_string()).collect::<Vec<_>>(),
                              1000),
            insert_sizes:  Self::join_with_flanked_commas(
                              &jxns.iter().map(|jxn| jxn.insert_sizes[1..jxn.insert_sizes.len() - 1].to_string()).collect::<Vec<_>>(), 
                              60),
            outer_node1s:  Self::join_with_flanked_commas(
                              &jxns.iter().map(|jxn| jxn.outer_node1s[1..jxn.outer_node1s.len() - 1].to_string()).collect::<Vec<_>>(), 
                              140),
            outer_node2s:  Self::join_with_flanked_commas(
                              &jxns.iter().map(|jxn| jxn.outer_node2s[1..jxn.outer_node2s.len() - 1].to_string()).collect::<Vec<_>>(), 
                              140),
            channels:      Self::join_with_flanked_commas(
                              &jxns.iter().map(|jxn| jxn.channels[1..jxn.channels.len() - 1].to_string()).collect::<Vec<_>>(),  
                              60),
            n_jxns:        Self::join_with_flanked_commas(
                              &jxns.iter().map(|jxn| jxn.n_jxns[1..jxn.n_jxns.len() - 1].to_string()).collect::<Vec<_>>(),     
                              20),
            is_duplicates: Self::join_with_flanked_commas(
                              &jxns.iter().map(|jxn| jxn.is_duplicates[1..jxn.is_duplicates.len() - 1].to_string()).collect::<Vec<_>>(), 
                              10),
            is_duplexes:   Self::join_with_flanked_commas(
                              &jxns.iter().map(|jxn| jxn.is_duplexes[1..jxn.is_duplexes.len() - 1].to_string()).collect::<Vec<_>>(), 
                              10),
            /* ------------------------------------------- */
            aln5_is:            Self::join_with_flanked_commas(
                                &jxns.iter().map(|jxn| jxn.aln5_is[1..jxn.aln5_is.len() - 1].to_string()).collect::<Vec<_>>(), 
                              10),
            qry_pos1_aln5_end3s: Self::join_with_flanked_commas(
                                &jxns.iter().map(|jxn| jxn.qry_pos1_aln5_end3s[1..jxn.qry_pos1_aln5_end3s.len() - 1].to_string()).collect::<Vec<_>>(), 
                              10),
            /* ------------------------------------------- */
            jxn_orientations:   Self::join_with_flanked_commas(
                              &jxns.iter().map(|jxn| jxn.jxn_orientations[1..jxn.jxn_orientations.len() - 1].to_string()).collect::<Vec<_>>(), 
                              10),
            jxn_failure_flags:  Self::join_with_flanked_commas(
                              &jxns.iter().map(|jxn| jxn.jxn_failure_flags[1..jxn.jxn_failure_flags.len() - 1].to_string()).collect::<Vec<_>>(), 
                              10),
            aln_failure_flags:  Self::join_with_flanked_commas(
                              &jxns.iter().map(|jxn| jxn.aln_failure_flags[1..jxn.aln_failure_flags.len() - 1].to_string()).collect::<Vec<_>>(), 
                              10),
            /* ------------------------------------------- */
            min_stem_lengths:  Self::join_with_flanked_commas(
                              &jxns.iter().map(|jxn| jxn.min_stem_lengths[1..jxn.min_stem_lengths.len() - 1].to_string()).collect::<Vec<_>>(), 
                              1),
            min_mapqs:         Self::join_with_flanked_commas(
                              &jxns.iter().map(|jxn| jxn.min_mapqs[1..jxn.min_mapqs.len() - 1].to_string()).collect::<Vec<_>>(), 
                              1),
            max_divergences:   Self::join_with_flanked_commas(
                              &jxns.iter().map(|jxn| jxn.max_divergences[1..jxn.max_divergences.len() - 1].to_string()).collect::<Vec<_>>(), 
                              1),
            /* ------------------------------------------- */
            n_instances:       jxns.iter().map(|jxn| jxn.n_instances).sum(),
            n_reads:           jxns.iter().map(|jxn| jxn.n_reads).sum(),
            n_instances_dedup: jxns.iter().map(|jxn| jxn.n_instances_dedup).sum(),
            n_reads_dedup:     jxns.iter().map(|jxn| jxn.n_reads_dedup).sum(),
            /* ------------------------------------------- */
            has_multi_jxn_read:      jxns.iter().any(|jxn| jxn.has_multi_jxn_read == 1) as u8,
            has_multi_instance_read: jxns.iter().any(|jxn| jxn.has_multi_instance_read == 1) as u8,
            has_duplex_read:         jxns.iter().any(|jxn| jxn.has_duplex_read == 1) as u8,
            is_bidirectional:        jxns.iter().any(|jxn| jxn.is_bidirectional == 1) as u8,
            jxn_failure_flag:        jxns.iter().fold(0u8, |acc, jxn| acc | jxn.jxn_failure_flag),
            aln_failure_flag:        jxns.iter().fold(0u8, |acc, jxn| acc | jxn.aln_failure_flag),
            any_was_chimeric:        jxns.iter().any(|jxn| jxn.any_was_chimeric == 1) as u8,
            any_was_not_chimeric:    jxns.iter().any(|jxn| jxn.any_was_not_chimeric == 1) as u8,
            min_stem_length:     jxns.iter().map(|jxn| jxn.min_stem_length).min().unwrap(),
            max_min_mapq:        jxns.iter().map(|jxn| jxn.max_min_mapq).max().unwrap(),
            min_max_divergence:  jxns.iter().map(|jxn| (jxn.min_max_divergence * 10000.0) as u32).min().unwrap() as f32 / 10000.0,
            /* ------------------------------------------- */
            is_intergenomic: best_jxn.is_intergenomic,
            target_1:        best_jxn.target_1.clone(), // targets are non-overlapping, only ever a single value expected
            target_dist_1:   best_jxn.target_dist_1,
            target_2:        best_jxn.target_2.clone(),
            target_dist_2:   best_jxn.target_dist_2,
            genes_1:         best_jxn.genes_1.clone(), // genes may overlap each other, so join with commas
            gene_dists_1:    best_jxn.gene_dists_1.clone(),
            genes_2:         best_jxn.genes_2.clone(),
            gene_dists_2:    best_jxn.gene_dists_2.clone(),
            is_excluded_1:   best_jxn.is_excluded_1,
            is_excluded_2:   best_jxn.is_excluded_2,
            /* ------------------------------------------- */
            // sample_names: Self::join_with_flanked_commas(&vec![tool.data_name.clone()], tool.data_name.len()),
            sample_bits: sample_bits,
            n_samples:   sample_bits.count_ones() as u8,
            /* ------------------------------------------- */
            bkpt_coverage_1: jxns.iter().map(|jxn| jxn.bkpt_coverage_1).sum(),
            bkpt_coverage_2: jxns.iter().map(|jxn| jxn.bkpt_coverage_2).sum(),
        }
    }

    /// Join values to strings with beginning and trailing comma separators.
    fn join_with_flanked_commas<T: std::fmt::Display>(d: &[T], unit_size: usize) -> String {
        use std::fmt::Write; // scope to avoid conflict with std::io::Write
        let capacity = d.len() * (unit_size + 1) + 1;
        let mut acc = String::with_capacity(capacity);
        acc.push(',');
        for item in d { write!(&mut acc, "{},", item).unwrap(); }
        acc
    }

    /// Write all final junctions to two sorted files: one sorted by breakpoint 
    /// node1, the other by breakpoint node2. Note that the nodes are not reordered,
    /// the same rows are just sorted differently.
    pub fn write_sorted(
        mut jxns: Vec<Self>,
        tool: &JunctionAnalysisTool,
    ){
        // print sorted by breakpoint 1
        jxns.par_sort_unstable_by(|a, b| 
                  (a.chrom_index1_1, a.ref_pos1_1, a.strand_index0_1)
            .cmp(&(b.chrom_index1_1, b.ref_pos1_1, b.strand_index0_1))
        );
        let writer = OutputCsv::open(&tool.final_jxns_file_1, Some(tool.n_cpu));
        writer.serialize_all(&jxns);

        // print sorted by breakpoint 2
        jxns.par_sort_unstable_by(|a, b| 
                  (a.chrom_index1_2, a.ref_pos1_2, a.strand_index0_2)
            .cmp(&(b.chrom_index1_2, b.ref_pos1_2, b.strand_index0_2))
        );
        let writer = OutputCsv::open(&tool.final_jxns_file_2, Some(tool.n_cpu));
        writer.serialize_all(&jxns);
    }

    /// Get a string representation of the junction offset type.
    pub fn get_offset_type(&self) -> &'static str {
             if self.offset > 0 { "insertion" } 
        else if self.offset < 0 { "microhomology" } 
                           else { "blunt"}
    }

    /// Get a string representation of the junction chimeric status.
    pub fn get_chimericity(&self) -> &'static str {
        if self.any_was_chimeric == 0 { "had no chimeric" } 
                                 else { "had >=1 chimeric" }
    }

    /// Get a string representation of the junction genome pairing.
    pub fn get_genomicity(&self) -> &'static str {
        if self.is_intergenomic == 0 { "single genome" } 
                                else { "intergenomic" }
    }

    /// Get a string representation of the junction coverage.
    pub fn get_coverage_type(&self) -> &'static str {
        if self.n_instances_dedup == 1 { "singleton" } 
         else if self.n_instances == 2 { "two instances" }
                                  else { ">=3 instances" }
    }

    /// Get the median insert size
    pub fn get_median_insert_size(&self) -> u32 {
        let mut insert_sizes: Vec<_> = self.insert_sizes
            .split(',')
            .map(str::trim)
            .filter(|&s| !s.is_empty())
            .filter_map(|s| {
                if let Ok(num) = s.parse::<i32>() {
                    Some(num.unsigned_abs()) // drop the sign, which reflects size selection filtering
                } else {
                    eprintln!("bad insert size: {}", s);
                    None
                }
            })
            .collect();
        if insert_sizes.len() == 0 { return 0 }
        insert_sizes.sort_unstable();
        let n = insert_sizes.len();
        if n % 2 == 0 {
            ((insert_sizes[n / 2 - 1] as f64 + insert_sizes[n / 2] as f64) / 2.0) as u32
        } else {
            insert_sizes[n / 2]
        }
    }  

    /// Convert strand indices to "+-", "-+", "++", "--".
    pub fn get_strand_encoding(&self) -> String {
        let strand1 = match self.strand_index0_1 {
            0 => "+",
            _ => "-", // thus, as would be used in paired end reads
        };
        let strand2 = match self.strand_index0_2 {
            0 => "-",
            _ => "+",
        };
        format!("{}{}", strand1, strand2)
    } 

    /// Adjust previously extracted HiFiRe3 junction metadata to 
    /// Sniffles-compatible VCF encoding. 
    /// 
    /// Resulting corrected positions and `sv_size` were empirically determined
    /// to adjust most HiFiRe3 junctions to Sniffles VCF coordinates and are
    /// consistent with the behavior of Sniffles in working from the last
    /// retained, i.e., aligned, bases, which is broadly similar to HiFiRe3, 
    /// except that HiFiRe3 performs complete junction characterization.
    /// 
    /// CuteSV POS, SVLEN and END encoding are NOT the same. See 
    pub fn adjust_to_vcf(
        &self,
        chrom_1: &mut String,
        pos1_1:  &mut u32,
        chrom_2: &mut String,
        pos1_2:  &mut u32,
        strands: &mut String, // as +- (DEL), -+ (DUP), ++ (INVF), -- (INVR)
    ) -> u32 {
        // HiFiRe3 and VCF both sort same-chrom junctions by position
        // but breakpoint assignments differ.       
        let sv_size: u32;
        match self.jxn_type {
            // DEL (+-)
            x if x == JunctionType::Deletion as u8 => {
                // -----MMM-----
                // -----MMM~~~~~~~MMM-----
                // H      1       2
                // S      1      2        H=HiFiRe3, S=Sniffles, empirically derived
                // -----|-----
                // -----~~~~~~~-----
                // H   1       2
                // S   1      2
                *pos1_2 -= 1;
                sv_size = *pos1_2 - *pos1_1; // thus does not count microhomology bases in the deletion
            },
            // DUP (-+)
            x if x == JunctionType::Duplication as u8 => {
                // -----MMM~~~~~~MMM~~~~~~MMM-----
                // -----MMM~~~~~~MMM-----
                // H    1          2
                // S   1           2
                // -----~~~~~~~|~~~~~~~-----
                // -----~~~~~~~-----
                // H    1     2
                // S   1      2
                *pos1_1 -= 1;
                sv_size = *pos1_2 - *pos1_1;
            },
            // INV-F (++)
            x if x == JunctionType::Inversion as u8 && self.strand_index0_1 == 0 => {
                // -----MMM7654321
                // -----MMM1234567mmm-----
                // H      1         2
                // S      1         2     
                // -----7654321
                // -----1234567-----
                // H   1      2
                // S   1      2
                sv_size = *pos1_2 - *pos1_1;
            },
            // INV-R (--)
            x if x == JunctionType::Inversion as u8 && self.strand_index0_1 == 1 => {
                //         7654321MMM-----
                // -----mmm1234567MMM-----
                // H    1         2
                // V   1         2
                //      7654321-----
                // -----1234567-----
                // H    1      2
                // V   1      2
                *pos1_2 += 1;
                sv_size = *pos1_2 - *pos1_1;
            },
            // TRANS - orient junctions using lexical-sorted chroms
            _ => {
                // TODO: adjust positions for microhomologies
                if chrom_1 > chrom_2 {
                    (*chrom_1, *pos1_1, *chrom_2, *pos1_2) = 
                    (chrom_2.clone(), pos1_2.clone(), chrom_1.clone(), pos1_1.clone());
                    *strands = match strands.as_str() {
                        "-+" => "+-".to_string(),
                        "+-" => "-+".to_string(),
                        _ => strands.clone() // INV types do not invert, they are the same regardless of node order
                    }
                }
                sv_size = 999_999_999;
            },
        } 
        sv_size
    }

    /// Adjust previously extracted HiFiRe3 junction metadata to CuteSV-compatible 
    /// VCF encoding. This is necessary for fuzzy matching of DEL events where
    /// HiFiRe3 and Sniffles2 `pos1_2` differ significantly from CuteSV.
    pub fn adjust_pos1_2_to_cutesv (
        &self,
        pos1_2: u32,
    ) -> u32 {  
        if self.jxn_type == JunctionType::Deletion as u8 {
            (pos1_2 as i32 - self.offset as i32) as u32
        } else {
            pos1_2
        }
    }
}
