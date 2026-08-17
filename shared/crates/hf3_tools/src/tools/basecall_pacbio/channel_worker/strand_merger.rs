//! Support for comparing and merging two PacBio strands into a single 
//! duplex consensus read to enforce duplex error correction. When
//! possible, performs three-strand error correction wherein bases
//! expected by a reference alignment are used to resolve strand
//! heteroduplex differences. See crate::snvs::tags for more detail.

// dependencies
use std::iter::repeat_n;
use rustc_hash::FxHashMap;
use crossbeam::channel::Sender;
use minimap2::{Aligner as Minimap2, Built, Strand};
use faimm::IndexedFasta;
use mdi::pub_key_constants;
use genomex::sequence::{rc_acgt_str, rc_acgtn_str, Aligner};
use crate::tools::basecall_pacbio::channel_kinetics::kinetics;
use crate::snvs::dd_tag::*;
use super::{StrandPair, KineticsInstance, MAX_READ_LEN};

// constants
pub_key_constants!(
    FAILED_TWO_STRAND_LEN_DELTA
    FAILED_TWO_STRAND_ALN_QUAL
);
// Smith-Waterman thresholds and operations
// TODO: expose MAX_SHIFT, MAX_SCORE_DIFF_PER_KB MIN_REF_ALN_LEN and/or MIN_MAPQ as options
const MAX_SHIFT: usize           = 50;   // allowed alignment shift from identity for fast alignment mode
const MAX_SCORE_DIFF_PER_KB: f64 = 50.0; // require a minimal match between strands
const SW_DEL_OP: &str            = "-";  // never written to SEQ or DD (except as the distinct HETERODUP_DEL_VS_REF op of the same char)
// ref_on_this minimap2 thresholds and operations
const MM_F_NO_PRINT_2ND: [u64; 1] = [16384]; // minimap2 flag to suppress printing/return of secondary alignments
const MM2_K: usize           = 19;  // minimizer k-mer length from map-hifi preset
const MM2_W: usize           = 19;  // minimizer window size  from map-hifi preset
const MAX_EXTENSION: usize   = MM2_K + MM2_W; // maximum alowed extension when determining how much of reference to extract for alignment to this.seq
const MIN_REF_ALN_LEN: usize = 100; // minimum ref_on_this alignment length to consider this.seq as a good reference match
const MIN_MAPQ: u32          = 50;  // minimum ref_on_this mapping quality  to consider this.seq as a good reference match
const REF_UNALIGNED_OP: &str = "x"; // marks positions in ref_on_this where this.seq did not align to reference (either end clip or internal between alns); never written to SEQ or DD
// DD tag thresholds and operations
pub const DD_CAPACITY: usize = 1000; // to minimize re-allocation while building DD tags

/// StrandMerger structure for merging two PacBio strands.
pub struct StrandMerger{
    aligner:             Aligner,
    pub has_ref_on_this: bool,
    ref_on_this:         Vec<String>,
    prev_on_this:        Vec<String>,
    prev_qual:           Vec<u8>,
    prev_offset:         usize,
    ident_unknown_len:   usize,
    ident_ref_len:       usize,
    dd_ops:              Vec<(&'static str, String)>,
    outcomes:            u8,
    is_dd_fusable_op:    FxHashMap<&'static str, bool>,
}
impl StrandMerger {
    /// Create a new StrandMerger structure.
    pub fn new() -> Self {
        StrandMerger { 
            aligner: Aligner::new_fast(
                MAX_READ_LEN + 1,
                MAX_READ_LEN + 1,
                MAX_SHIFT,
            ),
            has_ref_on_this:   false,
            ref_on_this:       Vec::with_capacity(MAX_READ_LEN),
            prev_on_this:      Vec::with_capacity(MAX_READ_LEN),
            prev_qual:         Vec::with_capacity(MAX_READ_LEN),
            prev_offset:       0,
            ident_unknown_len: 0,
            ident_ref_len:     0,
            dd_ops:            Vec::with_capacity(DD_CAPACITY),
            outcomes:          PERFECT_MATCH,
            is_dd_fusable_op: FxHashMap::from_iter([
                // end clip operation
                (PREV_CLIP_OP,     true),  // recorded as counts of runs of clipped bases
                // homoduplex operations
                (HOMODUP_UNKNOWN,  false), // recorded as the length of the identity run
                (HOMODUP_REF,      false), // recorded as the length of the identity run
                (HOMODUP_NOT_REF,  false), // recorded as one homoduplex base followed by the expected reference base(s), if any
                (HOMODUP_NOT_ALN,  true),  // recorded as the length of the identity run
                // heteroduplex indel operations
                (HETERODUP_INDEL_UNKNOWN,     true), // recorded as strings of indel bases on one strand relative to the other (and ref)
                (HETERODUP_INS_VS_REF,        true), // matched by an equal number N bases in output read unless one strand is ref
                (HETERODUP_DEL_VS_REF,        true),
                (HETERODUP_INDEL_NEITHER_REF, true),
                // heteroduplex substitution operations
                (HETERODUP_SUBS_UNKNOWN,     false), // recorded as the base context around the substitution, one substitution at a time
                (HETERODUP_SUBS_THIS_REF,    false), // matched by one N in output read unless one strand is ref
                (HETERODUP_SUBS_PREV_REF,    false),
                (HETERODUP_SUBS_NEITHER_REF, false),
            ]),
        }
    }

    /// Align this.seq to the reference genome to assemble ref_on_this map.
    /// Along the way, swap this and prev strands to make this == forward whenever possible.
    pub fn set_ref_on_this(
        &mut self,
        strand_pair: StrandPair,
        minimap2:    &Minimap2<Built>,
        fa:          &IndexedFasta,
        second_pass: bool, // true if this is the second call to this function after strand swapping
    ) -> Option<StrandPair> { // None if unmapped, otherwise a possibly swapped strand_pair

        // perform initial mapping of this.seq to the reference genome
        let mut this_maps = minimap2.map(
            strand_pair.this.seq.as_bytes(), 
            false, false,
            None, Some(&MM_F_NO_PRINT_2ND), None
        ).expect("Minimap2 mapping failure in strand merger");

        // short-circuit unmapped reads, they are dropped by caller as unusable
        // all other reads will return a StrandPair even if mappings had low MAPQ and ref_on_this is all REF_UNALIGNED_OP
        let n_maps = this_maps.len();
        if n_maps == 0 || 
           this_maps[0].target_id < 0 { // tid == -1 indicates unmapped
            return None;
        }

        // sort alignments to query order on this.seq
        if n_maps > 1 {
            this_maps.sort_by_key(|m| m.query_start);
        }

        // on first pass, re-orient this and prev so this is forward strand whenever possible
        if !second_pass && 
            this_maps[0         ].strand == Strand::Reverse &&  // nothing to do if this is already forward
            this_maps[n_maps - 1].strand == Strand::Reverse {   // nothing to do on RF combination, rc would also yield RF
            return self.set_ref_on_this(
                StrandPair {
                    sources: strand_pair.sources,
                    qname:   strand_pair.qname,
                    ff:      strand_pair.ff,
                    this:    strand_pair.prev, // flip this and prev so that this 5' end is forward
                    prev:    strand_pair.this,
                },
                minimap2, fa, true
            );
        }
        // recast this_to_ref alignment to create ref_on_this, in parallel to structure of prev_on_this
        let this_len = strand_pair.this.seq.len();
        self.has_ref_on_this = false;
        self.ref_on_this.clear();
        self.ref_on_this.resize(this_len, REF_UNALIGNED_OP.to_string());
        for map in this_maps {

            // ignore alignments that are too short, have poor mapping quality, or have large indels that foul fast-mode SW
            let this_start0 = map.query_start  as usize; // half open [this_start0, this_end1)
            let this_end1   = map.query_end    as usize;
            let ref_start0  = map.target_start as usize;
            let ref_end1    = map.target_end   as usize;
            let this_span   = this_end1.saturating_sub(this_start0);
            let ref_span    = ref_end1.saturating_sub(ref_start0);
            let delta_span  = (this_span as isize - ref_span as isize).unsigned_abs();
            if map.mapq   < MIN_MAPQ || 
               this_span  < MIN_REF_ALN_LEN || 
               ref_span   < MIN_REF_ALN_LEN ||
               delta_span > MAX_SHIFT {
                continue; 
            }

            // collect the reference sequence for alignment in the same orientation as this.seq (usually but not always forward strand)
            // extend the recovered bases at the 5' and 3' ends to deal with imprecise alignment boundaries in minimap2 fast mode
            let chrom = map.target_name
                .expect(&format!("Failed to extract chrom from minimap2 tid {}", map.target_id)); // TODO: is this the same as fa tid?
            let fai_tid = fa.fai().tid(&chrom)
                .expect(&format!("Failed to map minimap2 tid {}, chrom {}, to fasta tid", map.target_id, chrom));
            let chrom_size = fa.fai().size(fai_tid).unwrap();
            let extension5 = this_start0.min(MAX_EXTENSION);
            let extension3 = (this_len - this_end1).min(MAX_EXTENSION);
            let ref_seq = if map.strand == Strand::Reverse { // this.seq stays the same, we reverse-complement ref_seq to match
                let ref_seq = fa.view(
                    fai_tid, 
                    ref_start0.saturating_sub(extension3),
                    (ref_end1 + extension5).min(chrom_size) // extensions deal with imprecise alignment boundary in minimap2 fast mode
                )
                    .expect("Failed to extract ref_seq in strand_merger")
                    .to_string()
                    .to_ascii_uppercase(); // since GENOME_FASTA may have lower-case bases in repetitive regions
                rc_acgtn_str(&ref_seq)
            } else {
                fa.view(
                    fai_tid, 
                    ref_start0.saturating_sub(extension5), 
                    (ref_end1 + extension3).min(chrom_size)
                )
                    .expect("Failed to extract ref_seq in strand_merger")
                    .to_string()
                    .to_ascii_uppercase() // since GENOME_FASTA may have lower-case bases in repetitive regions
            };

            // use Smith-Waterman to align ref_seq to this.seq
            // use only the portion of this.seq expected to correspond to ref_seq as extended above
            let this_aln_start0 = this_start0.saturating_sub(extension5);
            let this_aln_end1 = (this_end1 + extension3).min(this_len);
            let this_tgt = &strand_pair.this.seq[this_aln_start0..this_aln_end1];
            let mut aln = self.aligner.align(&ref_seq, this_tgt, None, true);

            // reject poor ref_on_this alignments, e.g., due to indels larger than max_shift that foul alignment
            // panic on the unexpected outcome of failed alignment of ref_seq to this.seq
            if !aln.has_min_weighted_score(MAX_SCORE_DIFF_PER_KB, ||{
                panic!(
                    "Failed alignment of faimm ref_seq to this.seq in strand_merger; \n{}\n\n{}", 
                    ref_seq, this_tgt
                );
            }){ continue; }

            // add SW alignment to ref_on_this map
            // reads with >1 alignment may have a partial overwrite of junction (micro)homologies
            let aln_tgt_start0 = this_aln_start0 + aln.tgt_start0;
            let aln_tgt_end0   = this_aln_start0 + aln.tgt_end0;
            self.ref_on_this[aln_tgt_start0..=aln_tgt_end0].swap_with_slice(aln.qry_on_tgt.as_mut_slice());
            self.has_ref_on_this = true;
        }

        // return the possibly re-oriented StrandPair
        return Some(strand_pair);
    }

    /// Align rc(prev.seq) to this.seq to assemble prev_on_this map.
    pub fn set_prev_on_this(
        &mut self,
        strand_pair: &StrandPair,
    ) -> Option<&'static str> { // Option<failure_reason>

        // reject reads where strands are too different in length, should be a rare outcome
        let this = &strand_pair.this.seq;
        let prev = rc_acgt_str(&strand_pair.prev.seq);
        let this_len = this.len();
        let prev_len = prev.len();
        let delta_len  = (this_len as isize - prev_len as isize).unsigned_abs();
        if delta_len > MAX_SHIFT { return Some(FAILED_TWO_STRAND_LEN_DELTA); }

        // align rc(prev.seq) to this.seq, where this.seq is usually the forward reference strand
        let mut aln = self.aligner.align(&prev, this, None, true);

        // reject poor alignments between strands that might occur from failed strands or other poor ccs
        // panic on the unexpected outcome of failed alignment between strands
        if !aln.has_min_weighted_score(MAX_SCORE_DIFF_PER_KB, ||{
            panic!(
                "Failed alignment of prev.seq to this.seq in strand_merger; \n{}\n\n{}", 
                prev, this
            ) 
        }){ return Some(FAILED_TWO_STRAND_ALN_QUAL); }

        // add SW alignment to prev_on_this map, including terminal clips
        self.prev_on_this.clear();
        self.prev_on_this.resize(this_len, PREV_CLIP_OP.to_string()); // terminal clip operator
        self.prev_on_this[aln.tgt_start0..=aln.tgt_end0].swap_with_slice(aln.qry_on_tgt.as_mut_slice());
        self.prev_offset = aln.qry_start0;

        // set prev_qual to match the orientation of prev_on_this, but prev_qual is a simple Vec, not mapped to this
        self.prev_qual.clear();
        self.prev_qual.extend(strand_pair.prev.qual.iter().rev()); // reverse prev qual to match rc(prev.seq)

        // return the absence of prev_on_this failure
        None
    }

    /// Merge two PacBio strands into a single consensus read, masking
    /// heteroduplex bases where the strands disagree to N, unless
    /// one strand matches the reference genome in which case it is reported.
    /// Heteroduplex indels are omitted from the output SEQ if they
    /// cannot be verified by the reference genome since false
    /// insertions are more common than false deletions.
    /// 
    /// Create a difference tag similar to minimap2 cs tag to record the
    /// strand differences, and a kinetic parameters tag at positions
    /// with heteroduplex base substitutions.
    pub fn merge_strands(
        &mut self,
        strand_pair: &StrandPair,
        tx_kinetics: &Sender<KineticsInstance>,
        tx_corr_to_ref: &Sender<String>,
    ) -> (String, Vec<u8>, u8, String, Vec<u16>) {

        // initialize merger state
        self.ident_unknown_len = 0; 
        self.ident_ref_len     = 0; 
        self.dd_ops.clear();
        self.outcomes = PERFECT_MATCH;
        let this_len  = strand_pair.this.seq.len();
        let prev_len  = strand_pair.prev.seq.len();
        let mut seq  = String::with_capacity(this_len * 2);
        let mut qual = Vec::with_capacity(this_len * 2);
        let mut sk_tag: Vec<u16> = Vec::with_capacity(this_len / 4);

        // loop through this.seq on base at a time
        // there is one corresponding String entry in each of prev_on_this and ref_on_this
        for this_offset in 0..this_len {

            // collect the relevant bases and flags
            let this_base       = &strand_pair.this.seq[this_offset..=this_offset]; // always a single base
            let prev_bases   = &self.prev_on_this[this_offset]; // either one (un)matched base, - for a deletion, or multiple bases for an insertion
            let ref_bases    = &self.ref_on_this[ this_offset]; // as above
            let is_heteroduplex = this_base != prev_bases;
            let this_is_ref     = this_base == ref_bases;

            // homoduplex bases that agree between strands
            if !is_heteroduplex {
                seq.push_str(this_base); // commit duplex-validated sequenced bases as is regardless of ref match
                qual.push(strand_pair.this.qual[this_offset].max(self.prev_qual[self.prev_offset])); // use the higher base quality of the two strands for homoduplex bases
                self.prev_offset += 1;
                if self.has_ref_on_this {
                    if this_is_ref {
                        // only record kinetics benchmarks for bases fully validated by both strands and ref
                        kinetics::push(
                            &strand_pair.this, this_offset, 
                            is_heteroduplex, true, this_is_ref, 
                            tx_kinetics
                        );
                        self.ident_ref_len += 1;
                    } else {
                        if self.ident_ref_len > 0 {
                            self.dd_ops.push((HOMODUP_REF, self.ident_ref_len.to_string()));
                            self.ident_ref_len = 0;
                        }
                        if ref_bases == REF_UNALIGNED_OP {
                            self.dd_ops.push((HOMODUP_NOT_ALN, this_base.to_string()));
                        } else if ref_bases == SW_DEL_OP {
                            // thus, *T alone is a T base not present in reference, i.e., a read insertion
                            self.dd_ops.push((HOMODUP_NOT_REF, this_base.to_string())); 
                        } else {
                            // *TA, *TAC, etc. are a T that matches one or more reference bases
                            self.dd_ops.push((HOMODUP_NOT_REF, format!("{}{}", this_base, ref_bases))); 
                        }
                        self.outcomes |= if ref_bases == REF_UNALIGNED_OP{
                            HAS_REF_UNALIGNED
                        } else if ref_bases.len() > 1 || 
                                  ref_bases == SW_DEL_OP {
                            HAS_HOMODUPLEX_INDEL
                        } else {
                            HAS_HOMODUPLEX_SUBS
                        };
                    }
                } else {
                    self.ident_unknown_len += 1;
                }

            // heteroduplex bases that differ between strands
            // commit the strand that matches the reference, if any
            } else {

                // finish any previous homoduplex run
                if self.ident_ref_len > 0 {
                    self.dd_ops.push((HOMODUP_REF, self.ident_ref_len.to_string()));
                    self.ident_ref_len = 0;
                }
                if self.ident_unknown_len > 0 {
                    self.dd_ops.push((HOMODUP_UNKNOWN, self.ident_unknown_len.to_string()));
                    self.ident_unknown_len = 0;
                }

                // collect additional information on prev strand
                let prev_is_ref  = prev_bases == ref_bases;
                let n_prev_bases = prev_bases.len();

                // short-circuit end clips
                if prev_bases == PREV_CLIP_OP {
                    self.dd_ops.push((PREV_CLIP_OP, SEQ_MASKED_BASE.to_string()));
                    if this_is_ref {
                        seq.push_str(this_base);
                        qual.push(strand_pair.this.qual[this_offset]);
                    } else {
                        seq.push(SEQ_MASKED_BASE);
                        qual.push(0);
                    }
                    // prev_offset starts after a left clip
                    self.outcomes |= HAS_STRAND_CLIP;

                // insertion on prev relative to this
                } else if n_prev_bases > 1 {
                    // keep this_base that the insertion bases were recorded on as is, it is always a match
                    let indel_bases = prev_bases[0..n_prev_bases - 1].to_string();
                    if self.has_ref_on_this {
                        if this_is_ref {
                            self.dd_ops.push((HETERODUP_INS_VS_REF, indel_bases));
                            seq.push_str(this_base);
                            qual.push(strand_pair.this.qual[this_offset]);
                        } else if prev_is_ref {
                            self.dd_ops.push((HETERODUP_DEL_VS_REF, indel_bases));
                            seq.push_str(prev_bases); // the last base matches the one this_base after the indel
                            qual.extend(self.prev_qual[self.prev_offset..self.prev_offset + n_prev_bases].iter()); 
                        } else {
                            let masked_bases: String = repeat_n(SEQ_MASKED_BASE, indel_bases.len()).collect();
                            self.dd_ops.push((HETERODUP_INDEL_NEITHER_REF, indel_bases));
                            seq.push_str(&masked_bases);
                            for _ in 0..masked_bases.len() { qual.push(0); }
                            seq.push_str(this_base);
                            qual.push(strand_pair.this.qual[this_offset]);
                            self.outcomes |= HAS_STRAND_INDEL_NEITHER_REF;
                        }
                        self.ident_ref_len += 1; // just the one base after the indel
                    } else {
                        let masked_bases: String = repeat_n(SEQ_MASKED_BASE, indel_bases.len()).collect();
                        self.dd_ops.push((HETERODUP_INDEL_UNKNOWN, indel_bases));
                        seq.push_str(&masked_bases);
                        for _ in 0..masked_bases.len() { qual.push(0); }
                        seq.push_str(this_base);
                        qual.push(strand_pair.this.qual[this_offset]);
                        self.ident_unknown_len += 1;
                    }
                    self.prev_offset += n_prev_bases;
                    self.outcomes |= HAS_STRAND_INDEL;

                // deletion on prev relative to this
                } else if prev_bases == SW_DEL_OP { 
                    let indel_bases = this_base.to_string();
                    if self.has_ref_on_this {
                        if this_is_ref {
                            self.dd_ops.push((HETERODUP_DEL_VS_REF, indel_bases));
                            seq.push_str(this_base);
                            qual.push(strand_pair.this.qual[this_offset]);
                        } else if prev_is_ref {
                            self.dd_ops.push((HETERODUP_INS_VS_REF, indel_bases));
                        } else {
                            let masked_bases: String = repeat_n(SEQ_MASKED_BASE, indel_bases.len()).collect();
                            self.dd_ops.push((HETERODUP_INDEL_NEITHER_REF, indel_bases));
                            seq.push_str(&masked_bases);
                            for _ in 0..masked_bases.len() { qual.push(0); }
                            self.outcomes |= HAS_STRAND_INDEL_NEITHER_REF;
                        }
                    } else {
                            let masked_bases: String = repeat_n(SEQ_MASKED_BASE, indel_bases.len()).collect();
                            self.dd_ops.push((HETERODUP_INDEL_UNKNOWN, indel_bases));
                            seq.push_str(&masked_bases);
                            for _ in 0..masked_bases.len() { qual.push(0); }
                    }
                    // does not advance on prev_strand.seq
                    self.outcomes |= HAS_STRAND_INDEL;

                // base substitution between strands
                } else { 

                    // record the substitution base plus one-base context on either side
                    // note that qry3 is the reverse-complement of tgt3 to match prev_strand ip/pw orientation
                    // thus heteroduplex alignment of ATC (qry = rc(prev_strand.seq)) to ACC (tgt) 
                    // yields tag value *ACCGAT, which parses as:
                    //            ***  <<< the three frames of kinetics data used for assessment per strand
                    //           1 2 3
                    //     5' ---A C C--- tgt    = this_strand.seq
                    //     3' ---T A G--- qry_rc = prev_strand.seq; either strand might be reference
                    //        ---6 5 4---
                    //            ***
                    let ref_is_known = ref_bases != REF_UNALIGNED_OP;
                    let (this3, this_sk) = kinetics::push(
                        &strand_pair.this, this_offset, 
                        is_heteroduplex, ref_is_known, this_is_ref, 
                        tx_kinetics
                    );
                    let (prev3, prev_sk) = kinetics::push(
                        &strand_pair.prev, prev_len - self.prev_offset - 1, 
                        is_heteroduplex, ref_is_known, prev_is_ref, 
                        tx_kinetics
                    );
                    let dd_op = if self.has_ref_on_this {
                        if this_is_ref {
                            seq.push_str(this_base);
                            qual.push(strand_pair.this.qual[this_offset]);
                            let _ = tx_corr_to_ref.send(format!("this ref {} alt {}", this_base, prev_bases));
                            HETERODUP_SUBS_THIS_REF
                        } else if prev_is_ref {
                            seq.push_str(prev_bases);
                            qual.push(self.prev_qual[self.prev_offset]);
                            let _ = tx_corr_to_ref.send(format!("prev ref {} alt {}", prev_bases, this_base));
                            HETERODUP_SUBS_PREV_REF
                        } else {
                            seq.push(SEQ_MASKED_BASE);
                            qual.push(0);
                            self.outcomes |= HAS_STRAND_SUBS_NEITHER_REF;
                            HETERODUP_SUBS_NEITHER_REF
                        } 
                    } else {
                        seq.push(SEQ_MASKED_BASE);
                        qual.push(0);
                        HETERODUP_SUBS_UNKNOWN
                    };
                    self.dd_ops.push((dd_op, format!("{}{}", this3, prev3)));
                    sk_tag.extend(this_sk);
                    sk_tag.extend(prev_sk);
                    self.prev_offset += 1;
                    self.outcomes |= HAS_STRAND_SUBS;
                }
            }
        }

        // finish up any pending identical run
        if self.ident_unknown_len > 0 {
            self.dd_ops.push((HOMODUP_UNKNOWN, self.ident_unknown_len.to_string()));
        }
        if self.ident_ref_len > 0 {
            self.dd_ops.push((HOMODUP_REF, self.ident_ref_len.to_string()));
        }

        // collapse runs of the same op to create the final dd tag
        let mut dd_tag = String::with_capacity(DD_CAPACITY);
        let mut dd_iter = self.dd_ops.iter_mut().peekable();
        let mut wrk_op = dd_iter.next().unwrap();
        while let Some(next_op) = dd_iter.peek() {
            if wrk_op.0 == next_op.0 && 
               self.is_dd_fusable_op[wrk_op.0] {
                wrk_op.1.push_str(&next_op.1); // dd op fusing entails concatenation of op bases
                dd_iter.next();
            } else {
                Self::commit_dd_op(wrk_op, &mut dd_tag);
                wrk_op = dd_iter.next().unwrap();
            }
        }
        Self::commit_dd_op(wrk_op, &mut dd_tag);

        // return the results
        (seq, qual, self.outcomes, dd_tag, sk_tag)
    }

    /// Commit a dd tag op, including counting bases for certain ops.
    fn commit_dd_op(wrk_op: &mut (&'static str, String), dd_tag: &mut String) {
        if wrk_op.0 == PREV_CLIP_OP || // some base sequences are reduced to base counts in the dd tag
           wrk_op.0 == HOMODUP_NOT_ALN { 
            wrk_op.1 = wrk_op.1.len().to_string(); 
        }
        dd_tag.push_str(wrk_op.0);
        dd_tag.push_str(&wrk_op.1);
    }
}
