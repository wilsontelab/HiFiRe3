//! Support for aligning reads to the haplotype consensus for analyze_reads.
//! This script does the final, main work of subclonal SNV calling.

// imports
// use std::cmp::Ordering;
use std::str::from_utf8_unchecked;
use std::iter::repeat_n;
use minimap2::{Aligner as Minimap2, Built, Strand};
use crate::snvs::*;
use super::poa::*;
use super::*;

// constants
const ONE_THIRD:  f64 = 1.0 / 3.0;
const TWO_THIRDS: f64 = 2.0 / 3.0;
const POA_ANCHOR_LEN: usize = 5;
const POA_ANCHOR_SPAN: usize = POA_ANCHOR_LEN * 2;

impl SnvChromWorker {

    /// Assign one read to each haplotype consensus to resolve its haplotype.
    pub(super) fn assign_read_to_haplotype(
        &mut self,
        read_i:   ReadIndex,
    ) -> Option<Haplotype> {
        // if self.show_debug {
        //     eprintln!("");
        //     eprintln!("read_i {read_i}");
        // }

        // collect haplotype votes for this read by comparing read_on_ref 
        // tracking variants to ref_on_hap variants
        self.reset_hap_votes();
        for variant in &self.tracking_variants {
            let vmap = self.frag_vars.variant_map.get_mut(&variant).unwrap();
            let haplotype_matches = (
                self.hap_vars[&Haplotype::Haplotype1].contains(variant),
                self.hap_vars[&Haplotype::Haplotype2].contains(variant)
            );
            // if self.show_debug &&
            //    (vmap.read_map[read_i].has_var() || 
            //     vmap.read_map[read_i].is_informative) {
            //     eprintln!("variant {:?}", variant);
            //     eprintln!("haplotype_matches {:?}", haplotype_matches);
            // }
            if vmap.read_map[read_i].has_var() {
                match haplotype_matches {
                    m if m == (true, false) => {
                        *self.hap_votes.get_mut(&Haplotype::Haplotype1).unwrap() += 1;
                    },
                    m if m == (false, true) => {
                        *self.hap_votes.get_mut(&Haplotype::Haplotype2).unwrap() += 1;
                    },
                    // homozygous (true, true) and subclonal (false, false) don't vote
                    _ => {} 
                }
            } else if vmap.read_map[read_i].is_informative {
                match haplotype_matches {
                    m if m == (false, true) => {
                        *self.hap_votes.get_mut(&Haplotype::Haplotype1).unwrap() += 1;
                    },
                    m if m == (true, false) => {
                        *self.hap_votes.get_mut(&Haplotype::Haplotype2).unwrap() += 1;
                    },
                    _ => {}
                }
            }
            // non-informative variant no-calls don't vote
        }   
        // if self.show_debug {
        //     eprintln!("self.hap_votes {:?}", self.hap_votes);
        // }

        // read assignments are resolved by majority vote across informative 
        // variants when the bias is sufficient, otherwise they are discarded
        let n_hap1_votes = self.hap_votes[&Haplotype::Haplotype1] as f64;
        let n_hap2_votes = self.hap_votes[&Haplotype::Haplotype2] as f64;
        let frac_hap1 = n_hap1_votes / (n_hap1_votes + n_hap2_votes);
        if frac_hap1 >= TWO_THIRDS  { // 1 of 1, 2 of 3, 3 of 4, 4 of 5, 4 of 6, 5 of 7, 6 of 8, etc.
            Some(Haplotype::Haplotype1)
        } else if frac_hap1 <= ONE_THIRD {
            Some(Haplotype::Haplotype2)
        } else {
            // reads with no votes end here since frac_hap1 is NaN, which compares false
            // also ambiguous reads with 1 of 2, 2 of 4, 3 of 5, 5 of 8, etc.
            None 
        }
    }

    /// Build an error-corrected read consensus from its two constituent 
    /// strands and align it to its haplotype consensus to call subclonal
    /// variants relative to that consensus.
    pub(super) fn align_to_haplotype_consensus(
        &mut self,
        source_strands:     &mut SourceStrands,
        reads_on_haplotype: &mut FragmentHaplotypes,
        re_fragment:   &ReFragment,
        haplotype:     Haplotype,
        ref_pos0_map:  &mut Vec<ChromPos0>,
        str0_pos0_map: &mut Vec<usize>,
        str1_pos0_map: &mut Vec<usize>,
        reads:         &[ReadInstance],
        read_is:       Vec<ReadIndex>,
        hap_seq:       String,
        mm2_hap:       Minimap2<Built>,
    ) {
        let n_hap_bases = hap_seq.len();
        let n_reads = reads.len();
        let n_haplotype_reads = read_is.len();
        self.frag_vars.reset(n_haplotype_reads);

        // for each read assigned to the haplotype:
        //   - establish the POA consensus of its two strands on the haplotype consensus
        //   - align that read consensus to the haplotype consensus to call subclonal variants
        // by using the haplotype consensus as the POA seed, only homoduplex 
        // strand variants relative to the haplotype are retained in the error-
        // corrected read consensus
        for read_j in 0..n_haplotype_reads {
            let read_i = read_is[read_j];
            let read = &reads[read_i];
            if let Some(strands) = source_strands.by_read.get_mut(&read.qname){
                str0_pos0_map.clear();
                str1_pos0_map.clear();
                self.str_matches.clear();
                self.str_matches.extend(repeat_n(true, n_hap_bases));

                // use minimap2 to establish a minizer map of each strand on haplotype
                if !self.fill_strand_on_hap(
                    &mm2_hap,
                    &mut strands.0.seq, 
                    str0_pos0_map, 
                ) { continue; };
                if !self.fill_strand_on_hap(
                    &mm2_hap,
                    &mut strands.1.seq, 
                    str1_pos0_map, 
                ) { continue; };

                // use the minimizer map to efficiently step through haplotype
                // to create a read consensus by selective application of POA
                // to only local regions with candidate variants
                let read_consensus = self.do_minimizer_optimized_poa(
                    str0_pos0_map, str1_pos0_map,
                    &hap_seq, n_hap_bases, strands
                );

                // align the three-strand read consensus to the haplotype consensus
                let read_on_hap = &mm2_hap.map(
                    read_consensus.as_bytes(), 
                    true, 
                    false, 
                    None, 
                    Some(&MM_F_NO_PRINT_2ND), 
                    None
                ).expect("Minimap2 failed at read_on_hap")[0];
                let Some(aln) = &read_on_hap.alignment else { continue; };
                let Some(cs) = &aln.cs else { continue; }; // never expected to fail
                self.encoding.prepare_read_on_hap(
                    re_fragment, read, read_on_hap.target_start as usize
                );
                self.process_cs_tag(
                    reads_on_haplotype, re_fragment, haplotype, 
                    Some((read_on_hap.query_start, read_on_hap.target_start)), 
                    Some(cs), ref_pos0_map,
                    // not read_i; read_j indexes into frag_vars
                    read_j, read, 
                    &source_strands.by_read[&read.qname], 
                    &str0_pos0_map, &str1_pos0_map
                ); 
                reads_on_haplotype.insert_encoding(
                    re_fragment, haplotype, self.encoding.clone()
                ); 
            }
        }

        // call subclonal variants aggregated over all haplotype reads (might be none)
        for variant in self.frag_vars.variant_map.keys(){
            let vmap =  &self.frag_vars.variant_map[&variant];
            let read_js: Vec<ReadIndex> = vmap.read_map.iter()
                .enumerate()
                .filter_map(|(read_j, r)|{
                    if r.has_var() { Some(read_j) } else { None }
                }).collect();
            let max_min_qual = read_js.iter()
                .map(|read_j| {
                    let read_i = read_is[*read_j];
                    let min_qual= vmap.read_map[*read_j].min_qual;
                    self.variant_reads_tally.add_subclonal_variant(
                        &reads[read_i], re_fragment, &haplotype, 
                        variant, min_qual, 
                        n_haplotype_reads, n_reads
                    );
                    min_qual
                })
                .max()
                .unwrap_or_default();
            self.variant_tally.add_subclonal(
                &variant, reads, &read_is, &read_js, 
                max_min_qual
            );
        }
    } 

    /// Use minimap2 to establish a minizer-assisted map of each strand on 
    /// haplotype. Reverse-complement the one strand that is expected to need it  
    /// for later strand-resolved use in POA.
    fn fill_strand_on_hap(
        &mut self,
        mm2_hap:      &Minimap2<Built>,
        str_seq:      &mut Vec<u8>,
        str_pos0_map: &mut Vec<usize>,
    ) -> bool {
        let str_on_hap = &mut mm2_hap.map(
            str_seq, 
            true, 
            false, 
            None, 
            Some(&MM_F_NO_PRINT_2ND), 
            None
        ).expect("Minimap2 failed at str_on_hap")[0];
        let Some(aln) = &str_on_hap.alignment else { return false; };
        let Some(cs) = &aln.cs else { return false; }; // never expected to fail

        if str_on_hap.strand == Strand::Reverse {
            str_on_hap.query_start = str_seq.len() as i32 - str_on_hap.query_end;
            *str_seq = Poa::reverse_complement(str_seq);
        }

        self.process_cs_str_on_hap(
            str_pos0_map, 
            str_on_hap.query_start as usize, // leftmost, see above
            str_on_hap.target_start as usize, 
            cs
        );
        true
    }
    
    /// Parse a cs tag from a single read strand sequence aligned to the  
    /// haplotype consensus.
    fn process_cs_str_on_hap(
        &mut self,
        str_pos0_map: &mut Vec<usize>,
        mut str_pos0: usize,
        mut hap_pos0: usize,
        cs: &str,
    ) {
        // fill any left-side gaps in the alignment; cannot call variants
        if hap_pos0 > 0 {
            (0..hap_pos0).for_each(|hap_pos0| {
                str_pos0_map.push(0);
                self.str_matches[hap_pos0] = false;
            });
        }

        // process the cs tag span
        self.reset_cs_variant();
        let mut chars = cs.chars();
        self.cs_op = chars.next().unwrap();
        self.op_val.clear();
        while let Some(char) = chars.next() {
            if char.is_alphanumeric() {
                self.op_val.push(char);
            } else {
                self.handle_cs_str_on_hap(
                    str_pos0_map, &mut str_pos0, &mut hap_pos0
                );
                self.cs_op = char;
                self.op_val.clear();
            }
        }
        self.handle_cs_str_on_hap(
            str_pos0_map, &mut str_pos0, &mut hap_pos0
        );

        // fill any right-side gaps in the alignment; cannot call variants
        while hap_pos0 < self.str_matches.len() {
            str_pos0_map.push(0);
            self.str_matches[hap_pos0] = false;
            hap_pos0 += 1;
        }
    }
    fn handle_cs_str_on_hap (
        &mut self,
        str_pos0_map: &mut Vec<usize>,
        str_pos0:     &mut usize,
        hap_pos0:     &mut usize,
    ){
        match self.cs_op {
            ':' => { // :[0-9]+   Identical sequence length
                let len = self.op_val.parse::<usize>().unwrap();
                (0..len).for_each(|_| {
                    str_pos0_map.push(*str_pos0);
                    *str_pos0 += 1;
                    *hap_pos0 += 1;
                });
            },
            '*' => { // *[acgtn][acgtn]   Substitution: target to query
                //     S
                // rrrrRrrrr
                // qqqqQqqqq
                //     A
                self.str_matches[*hap_pos0] = false;
                str_pos0_map.push(*str_pos0);
                *str_pos0 += 1;
                *hap_pos0 += 1;
            },
            '+' => { // +[acgtn]+   Insertion to the target
                //    *III 
                // rrrr   Rrrr
                // qqqqQqqqqqq
                //    aA Aa
                let n_ins_bases = self.op_val.len();
                self.str_matches[*hap_pos0 - 1] = false; // force both flanking bases to POA
                self.str_matches[*hap_pos0] = false;
                *str_pos0 += n_ins_bases;
            },
            '-' => { // -[acgtn]+   Deletion from the target
                //     DDD
                // rrrrRrrrrrr
                // qqqq   Qqqq
                //   aA   Aa
                let n_del_bases = self.op_val.len();
                (0..n_del_bases).for_each(|_| {
                    str_pos0_map.push(*str_pos0 - 1);
                    self.str_matches[*hap_pos0] = false;
                    *hap_pos0 += 1;
                });
            },
            _ => panic!("Unexpected CS tag operation: {}", self.cs_op),
        }
    }

    /// Perform local POA while committing identical spans as is.
    fn do_minimizer_optimized_poa(
        &mut self,
        str0_pos0_map: &Vec<usize>,
        str1_pos0_map: &Vec<usize>,
        hap_seq:       &str,
        n_hap_bases:   usize,
        strands:       &(SourceStrand, SourceStrand)
    ) -> String {
        let mut read_consensus = String::new();

        let mut hap_pos0:    usize = 0; // leftmost pos0 of the next encountered chunk
        let mut left_start0: usize = 0; // leftmost pos0 of the uncommitted match span left of variant
        let mut left_end1:   usize = 0; // righmost pos1 of the uncommitted match span left of variant
        let mut is_variant:  bool = false; // if true, match_left is set and we have a variant needing POA

        self.str_matches.chunk_by(|a, b| a == b).for_each(|m|{
            let is_hap_match = m[0];
            let n_hap_pos = m.len();

            // in a span where both strands matched the haplotype
            if is_hap_match {

                // initialize the first (and possibly only) matching span
                if !is_variant {
                    left_start0 = hap_pos0;
                    left_end1 = hap_pos0 + n_hap_pos;

                // process a variant span by POA if sufficient anchors
                // too-short matching anchors are included in the POA span 
                } else if n_hap_pos >= POA_ANCHOR_SPAN {
                    let hap_start0 = left_end1.saturating_sub(POA_ANCHOR_LEN);
                    let hap_end1 = (hap_pos0 + POA_ANCHOR_LEN).min(n_hap_bases);
                    if hap_start0 > left_start0 {
                        read_consensus.push_str(&hap_seq[left_start0..hap_start0]);
                    }
                    self.poa.seed_new_graph(hap_seq[hap_start0..hap_end1].as_bytes());

                    let str_start0 = str0_pos0_map[hap_start0];
                    let str_end1 = str0_pos0_map[hap_end1 - 1] + 1;
                    self.poa.add_read(&strands.0.seq[str_start0..str_end1]);

                    let str_start0 = str1_pos0_map[hap_start0];
                    let str_end1 = str1_pos0_map[hap_end1 - 1] + 1;
                    self.poa.add_read(&strands.1.seq[str_start0..str_end1]);

                    let var_consensus = self.poa.get_heaviest_path();
                    let var_consensus = unsafe{ from_utf8_unchecked(&var_consensus) };
                    read_consensus.push_str(var_consensus);

                    left_start0 = hap_end1; // thus, not including the flanking bases committed with POA
                    left_end1 = hap_pos0 + n_hap_pos;
                    is_variant = false;
                } 

            // in a span where at least one strand differed from haplotype
            } else {

                // handle flanking variant gaps, commit as haplotype consensus
                if hap_pos0 + n_hap_pos == n_hap_bases {
                    // do nothing, handled below after for_each terminates
                } else if hap_pos0 == 0 {
                    read_consensus.push_str(&hap_seq[0..n_hap_pos]);
                
                // in read middle, flag that we have a variant span pending right anchor for POA
                } else {
                    is_variant = true;
                }
            }
            hap_pos0 += n_hap_pos;
        });
        if left_start0 < n_hap_bases {
            read_consensus.push_str(&hap_seq[left_start0..n_hap_bases]);
        }
        read_consensus
    }
}

        // // align each haplotype read to its consensus if not already done
        // for read_j in 0..n_haplotype_reads {
        //     let read_i = read_is[read_j];
        //     let read = &reads[read_i];
        //     if let Some(assignments) = &assignments {
        //         let assignment = &assignments[read_j];
        //         self.encoding.prepare_read_on_hap(
        //             re_fragment, read, assignment.target_start as usize
        //         );
        //         self.process_cs_tag(
        //             reads_on_haplotype, re_fragment, haplotype, 
        //             Some((assignment.query_start, assignment.target_start)), 
        //             Some(&assignment.cs), ref_pos0_map,
        //             // not read_i; read_j indexes into frag_vars
        //             read_j, read, &read_masks[read_i] 
        //         );  
        //     } else { // this path is used by fully homozygous fragments
        //         let mm2_hap = mm2_hap.as_ref().unwrap();
        //         let read_on_hap = &mm2_hap.map(
        //             &read.seq_bytes, 
        //             true, 
        //             false, 
        //             None, 
        //             Some(&MM_F_NO_PRINT_2ND), 
        //             None
        //         ).expect("Minimap2 failed at read_on_hap")[0];
        //         let Some(aln) = &read_on_hap.alignment else { return; };
        //         let Some(cs) = &aln.cs else { return; }; // never expected to fail
        //         self.encoding.prepare_read_on_hap(
        //             re_fragment, read, read_on_hap.target_start as usize
        //         );
        //         self.process_cs_tag(
        //             reads_on_haplotype, re_fragment, haplotype, 
        //             Some((read_on_hap.query_start, read_on_hap.target_start)), 
        //             Some(cs), ref_pos0_map,
        //             // not read_i; read_j indexes into frag_vars
        //             read_j, read, &read_masks[read_i] 
        //         );     
        //     }
        //     reads_on_haplotype.insert_encoding(
        //         re_fragment, haplotype, self.encoding.clone()
        //     ); 
        // }

// pub struct ReadAssignment {
//     pub haplotype:    Haplotype,
//     // pub cs:           String,
//     // pub query_start:  i32,
//     // pub target_start: i32,
// } 

    // /// Get the read assignment for an unambiguously assigned read haplotype.
    // fn get_read_haplotype(
    //     read:      &ReadInstance,
    //     mm2_hap:   &Minimap2<Built>,
    //     haplotype: Haplotype,
    // ) -> Option<ReadAssignment> {
    //     let read_on_hap = &mm2_hap.map(
    //         &read.seq_bytes, 
    //         true, // cs not used, but need full score calculation
    //         false, 
    //         None, 
    //         Some(&MM_F_NO_PRINT_2ND), 
    //         None
    //     ).expect("Minimap2 failed at read_on_hap")[0];
    //     let Some(aln) = &read_on_hap.alignment else { return None; };
    //     let Some(cs) = &aln.cs else { return None; }; // don't expect failures
    //     Some(ReadAssignment{
    //         haplotype, 
    //         cs: cs.clone(), 
    //         query_start:  read_on_hap.query_start,
    //         target_start: read_on_hap.target_start,
    //     })
    // }

    // /// Get the read assignment for read with conflicting haplotype votes
    // /// by alignment score. Deprecated: too prone to wrong conclusions.
    // fn _get_read_haplotype_by_score(
    //     &mut self,
    //     read:     &ReadInstance,
    //     mm2_hap1: &Minimap2<Built>,
    //     mm2_hap2: &Minimap2<Built>,
    // ) -> Option<ReadAssignment> {
    //     let read_on_hap1 = &mm2_hap1.map(
    //         &read.seq_bytes, 
    //         true, // cs not used, but need full score calculation
    //         false, 
    //         None, 
    //         Some(&MM_F_NO_PRINT_2ND), 
    //         None
    //     ).expect("Minimap2 failed at read_on_hap1")[0];
    //     let read_on_hap2 = &mm2_hap2.map(
    //         &read.seq_bytes, 
    //         true, // cs not used, but need full score calculation
    //         false, 
    //         None, 
    //         Some(&MM_F_NO_PRINT_2ND), 
    //         None
    //     ).expect("Minimap2 failed at read_on_hap2")[0];
    //     let Some(aln1) = &read_on_hap1.alignment else { return None; };
    //     let Some(aln2) = &read_on_hap2.alignment else { return None; };
    //     let Some(score1) = &aln1.alignment_score else { return None; };
    //     let Some(score2) = &aln2.alignment_score else { return None; };
    //     let Some(cs1) = &aln1.cs else { return None; }; // don't expect failures
    //     let Some(cs2) = &aln2.cs else { return None; };
    //     match score1.cmp(score2) {
    //         Ordering::Greater => Some(ReadAssignment{
    //             haplotype: Haplotype::Haplotype1, 
    //             cs: cs1.clone(), 
    //             query_start:  read_on_hap1.query_start,
    //             target_start: read_on_hap1.target_start,
    //         }),
    //         Ordering::Less => Some(ReadAssignment{
    //             haplotype: Haplotype::Haplotype2, 
    //             cs: cs2.clone(), 
    //             query_start:  read_on_hap2.query_start,
    //             target_start: read_on_hap2.target_start,
    //         }),
    //         // rare reads where we truly can't determine the haplotype are dropped
    //         Ordering::Equal => None,
    //     }

    // }
