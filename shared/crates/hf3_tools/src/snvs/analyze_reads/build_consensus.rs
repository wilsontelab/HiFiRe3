//! Support for haplotype consensus building for analyze_reads.

// imports
use std::iter::repeat_n;
use minimap2::{Aligner as Minimap2, Built};
use crate::snvs::*;
use super::*;

impl SnvChromWorker {

    /// Use unambiguous haplotype reads to build a haplotype consensus sequence.
    /// Like all of the contributing reads, the returned consensus always has 
    /// the same endpoints as the host ReFragment.
    /// 
    /// This function is called with two or more input reads, but never one.
    pub(super) fn build_haplotype_consensus(
        &mut self,
        haplotype_consensuses: &mut HaplotypeConsensuses,
        re_fragment:  &ReFragment,
        haplotype:    Haplotype,
        ref_pos0_map: &mut Vec<ChromPos0>,
        reads:        &[ReadInstance],
        read_is:      &[ReadIndex],
    ) -> Option<(String, Minimap2<Built>)> {

        // select the read to use as the index during consensus assembly
        // prefer the one with the highest initial reference alignment score
        let read0_i = *read_is.iter()
            .max_by(|&&a, &&b|{
                reads[a].aln_score.cmp(&reads[b].aln_score)
            }).unwrap();

        // initialize the index strand sequence as the consensus builder comparator
        // consensus assembly uses read0_i SEQ as found on top ref strand
        let seq0_bytes = &reads[read0_i].seq_bytes;
        let mm2_seq0 = self.minimap2.clone()
            .with_seq(seq0_bytes)
            .expect("Failed to initialize minimap2 in build_haplotype_consensus()");
        self.seq0_bases.clear();
        self.seq0_bases.extend(seq0_bytes.iter()
            .map(|&b| match b {
                b'A' => "A".to_string(),
                b'C' => "C".to_string(),
                b'G' => "G".to_string(),
                b'T' => "T".to_string(),
                _    => "N".to_string(),
            })
        );
        self.cs_map.clear();
        self.cs_map.extend(self.seq0_bases.iter()
            .map(|b| {
                let mut m = FxHashMap::default();
                m.insert(b.clone(), 1);
                m
            })
        );

        // align each remaining read to target to count alternative bases
        // as throughout, read SEQ is always ref top strand
        read_is.iter().for_each(|read_i| {
            if *read_i == read0_i { return; } // skip the read in use as index target
            let read_on_seq0 = &mm2_seq0.map(
                &reads[*read_i].seq_bytes, 
                true, 
                false, 
                None, 
                Some(&MM_F_NO_PRINT_2ND), 
                None
            ).expect("Minimap2 failed at read_on_seq0")[0];
            let Some(aln) = &read_on_seq0.alignment else { return; };
            let Some(cs) = &aln.cs else { return; };
            self.process_cs_read_on_index(
                read_on_seq0.target_start as usize, 
                cs
            );
        });

        // scan the matrix by base to establish the haplotype consensus
        let cs_map_len = self.cs_map.len();
        let mut consensus = String::with_capacity(cs_map_len + 100);
        for cs_map_pos0 in 0..cs_map_len {
            let bases = self.cs_map[cs_map_pos0].iter()
                .max_by(|a, b| a.1.cmp(&b.1))
                .map(|(bases, _)| bases)
                .unwrap();
            if bases != "-" {
                consensus.push_str(bases); // could be multiple bases at an insertion
            }
        }

        // align the fragment reference span to the haplotype consensus, i.e., ref_on_hap
        let mm2_hap = self.minimap2.clone()
            .with_seq(consensus.as_bytes())
            .expect("Failed to initialize minimap2 in build_haplotype_consensus()");
        let (ref_seq, _) = haplotype_consensuses.get(re_fragment, Haplotype::Unspecified);
        let ref_on_hap = &mm2_hap.map(
            ref_seq.as_bytes(), 
            true, 
            false, 
            None, 
            Some(&MM_F_NO_PRINT_2ND), 
            None
        ).expect("Minimap2 failed at ref_on_hap")[0];
        let Some(aln) = &ref_on_hap.alignment else { return None; };
        let Some(cs) = &aln.cs else { return None; }; // never expected to fail

        // create a map of reference positions per each haplotype consensus position
        // and a hap_vs_ref encoding to retain a memory of clonal variants relative to reference
        ref_pos0_map.clear();
        self.hap_vs_ref.clear();
        let hap_vs_ref = self.process_cs_ref_on_hap(
            ref_pos0_map, re_fragment, haplotype, 
            ref_on_hap.target_start, ref_on_hap.target_end, ref_on_hap.query_start, 
            consensus.len(), cs, 
            reads, read_is,
        );

        // cache the haplotype consensus for printing
        haplotype_consensuses.insert(
            re_fragment, haplotype, 
            consensus.clone(), Some(hap_vs_ref)
        );

        // return the built minimap2 for subsequent read alignment
        Some((consensus, mm2_hap))
    }

    /// Parse a cs tag from one read aligned to an index read during haplotype
    /// consensus building.
    fn process_cs_read_on_index(
        &mut self,
        mut tgt_pos0: usize,
        cs:           &str,
    ){
        let mut chars = cs.chars();
        self.cs_op = chars.next().unwrap();
        self.op_val.clear();
        while let Some(char) = chars.next() {
            if char.is_alphanumeric() {
                self.op_val.push(char);
            } else {
                self.handle_cs_read_on_index(&mut tgt_pos0);
                self.cs_op = char;
                self.op_val.clear();
            }
        }
        self.handle_cs_read_on_index(&mut tgt_pos0);
    }
    fn handle_cs_read_on_index (
        &mut self,
        tgt_pos0: &mut usize,
    ){
        match self.cs_op {
            ':' => { // :[0-9]+   Identical sequence length
                (0..self.op_val.parse::<usize>().unwrap()).for_each(|_| {
                    *self.cs_map[*tgt_pos0].get_mut(&self.seq0_bases[*tgt_pos0]).unwrap() += 1;
                    *tgt_pos0 += 1;
                });
            },
            '*' => { // *[acgtn][acgtn]   Substitution: target to query
                self.cs_map[*tgt_pos0]
                    .entry(self.op_val[1..=1].to_ascii_uppercase())
                    .and_modify(|n| *n += 1)
                    .or_insert(1);
                *tgt_pos0 += 1;
            },
            '+' => { // +[acgtn]+   Insertion to the target
                let alt = self.op_val.to_ascii_uppercase();
                self.cs_map[*tgt_pos0 - 1]
                    .entry(format!("{}{}", self.seq0_bases[*tgt_pos0 - 1], alt))
                    .and_modify(|n| *n += 1)
                    .or_insert(1);
            },
            '-' => { // -[acgtn]+   Deletion from the target
                (0..self.op_val.len()).for_each(|_| {
                    self.cs_map[*tgt_pos0]
                        .entry("-".to_string())
                        .and_modify(|n| *n += 1)
                        .or_insert(1);
                    *tgt_pos0 += 1;
                });
            },
            _ => panic!("Unexpected CS tag operation: {}", self.cs_op),
        }
    }

    /// Parse a cs tag from a reference sequence aligned to the haplotype 
    /// consensus.
    fn process_cs_ref_on_hap(
        &mut self,
        ref_pos0_map: &mut Vec<ChromPos0>,
        re_fragment:  &ReFragment,
        haplotype:    Haplotype,
        target_start: i32,
        target_end:   i32,
        query_start:  i32,
        n_hap_bases:  usize,
        cs:           &str,
        reads:        &[ReadInstance],
        read_is:      &[ReadIndex],
    ) -> String {
        self.reset_cs_variant();

        if target_start > 0 {
            let count = target_start as usize;
            ref_pos0_map.extend(repeat_n(re_fragment.start0, count));
            self.hap_vs_ref.extend(repeat_n("+".to_string(), count));
        }

        let mut ref_pos0 = re_fragment.start0 + query_start as u32;
        let mut chars = cs.chars();
        self.cs_op = chars.next().unwrap();
        self.op_val.clear();
        while let Some(char) = chars.next() {
            if char.is_alphanumeric() {
                self.op_val.push(char);
            } else {
                self.handle_cs_op_ref_on_hap(
                    ref_pos0_map, &mut ref_pos0,
                    reads, read_is, 
                    re_fragment, haplotype,
                );
                self.cs_op = char;
                self.op_val.clear();
            }
        }
        self.handle_cs_op_ref_on_hap(
            ref_pos0_map, &mut ref_pos0,
            reads, read_is, 
            re_fragment, haplotype,
        );

        if target_end < n_hap_bases as i32 - 1 {
            let count = n_hap_bases - 1 - target_end as usize;
            ref_pos0_map.extend(repeat_n(re_fragment.end1 - 1, count));
            self.hap_vs_ref.extend(repeat_n("+".to_string(), count));
        }

        self.hap_vs_ref.join("")
    }
    fn handle_cs_op_ref_on_hap (
        &mut self,
        ref_pos0_map: &mut Vec<ChromPos0>,
        ref_pos0:     &mut ChromPos0,
        reads:        &[ReadInstance],
        read_is:      &[ReadIndex],
        re_fragment:  &ReFragment,
        haplotype:    Haplotype,
    ){
        match self.cs_op {
            ':' => { // :[0-9]+   Identical sequence length
                if self.var_tgt_pos0.is_some() { // commit any preceding variant stretch
                    if self.allowed {
                        let ref_pos0 = self.var_tgt_pos0.unwrap();
                        let mut variant = Variant::new(
                            ref_pos0,
                            ref_pos0,
                            &self.tgt_bases,
                            &self.alt_bases,
                            re_fragment, haplotype
                        );
                        // if self.show_debug {
                        //     eprintln!("ref_on_hap variant {:?}", variant);
                        // }
                        self.variant_tally.add_clonal(
                            &variant, reads, read_is
                        );
                        variant.haplotype = Haplotype::Unspecified;
                        self.hap_vars.get_mut(&haplotype).unwrap().insert(variant);
                    }
                    self.reset_cs_variant();
                }
                let len = self.op_val.parse::<usize>().unwrap();
                (0..len).for_each(|_| {
                    ref_pos0_map.push(*ref_pos0);
                    *ref_pos0 += 1;
                });
                self.hap_vs_ref.push(format!("={}", len));
            },
            '*' => { // *[acgtn][acgtn]   Substitution: target to query
                //     S
                // rrrrRrrrr
                // qqqqQqqqq
                //     A
                let tgt = self.op_val[1..=1].to_ascii_uppercase(); 
                let alt = self.op_val[0..=0].to_ascii_uppercase(); // yes, 0 is consensus
                self.tgt_bases.push_str(&tgt);
                self.alt_bases.push_str(&alt);
                self.allowed &= alt != "N";
                self.allowed &= !self.simple_repeats.binary_search(*ref_pos0, 1);
                if self.var_tgt_pos0.is_none() { self.var_tgt_pos0 = Some(*ref_pos0); }
                ref_pos0_map.push(*ref_pos0);
                self.hap_vs_ref.push(alt);
                *ref_pos0    += 1;
            },
            '+' => { // +[acgtn]+   Insertion to the target
                //     III
                // rrrrRrrrrrr
                // qqqq   Qqqq
                //   aA   Aa
                let n_del_bases = self.op_val.len() as u32; // yes, del in consensus
                let tgt = self.op_val.to_ascii_uppercase();
                self.tgt_bases.push_str(&tgt);
                self.allowed &= !self.simple_repeats.binary_search(*ref_pos0, n_del_bases);
                if self.var_tgt_pos0.is_none() { self.var_tgt_pos0 = Some(*ref_pos0); }
                if let Some(prev) = self.hap_vs_ref.pop() {
                    if prev.starts_with("=") {
                        let mut len = prev[1..prev.len()]
                            .parse::<u32>()
                            .unwrap();
                        len -= 1;
                        if len > 0 {
                            self.hap_vs_ref.push(format!("={}", len));
                        }
                        self.hap_vs_ref.push("-".to_string());                          
                    }
                }
                *ref_pos0    += n_del_bases;
                // no action on alt_bases
            },
            '-' => { // -[acgtn]+   Deletion from the target
                //    *DDD   
                // rrrr   Rrrr
                // qqqqQqqqqqq
                //    aA Aa
                let n_ins_bases = self.op_val.len() as u32; // yes, ins in consensus
                let alt = self.op_val.to_ascii_uppercase();
                self.alt_bases.push_str(&alt);
                self.allowed &= !alt.contains("N");
                self.allowed &= !self.simple_repeats.binary_search(*ref_pos0 - 1, 2);
                if self.var_tgt_pos0.is_none() { self.var_tgt_pos0 = Some(*ref_pos0 - 1); }
                (0..n_ins_bases).for_each(|_| {
                    ref_pos0_map.push(*ref_pos0 - 1);
                    self.hap_vs_ref.push("+".to_string());
                });
                // no action on tgt_bases or ref_pos0
            },
            _ => panic!("Unexpected CS tag operation: {}", self.cs_op),
        }
    }

    // /// Use unambiguous haplotype reads to build a haplotype consensus sequence.
    // /// Like all of the contributing reads, the returned consensus always has 
    // /// the same endpoints as the host ReFragment.
    // /// 
    // /// This function is called with two or more input reads, but never one.
    // /// 
    // /// The consensus is built from original PacBio --by-strand strand consensuses.
    // pub(super) fn build_haplotype_consensus_poa(
    //     &mut self,
    //     source_strands: &SourceStrands,
    //     haplotype_consensuses: &mut HaplotypeConsensuses,
    //     re_fragment:  &ReFragment,
    //     haplotype:    Haplotype,
    //     ref_pos0_map: &mut Vec<ChromPos0>,
    //     reads:        &[ReadInstance],
    //     read_is:      &[ReadIndex],
    // ) -> Option<(Vec<u8>, Minimap2<Built>)> {

    //     // get the reference sequence of the fragment on the top strand
    //     let (ref_seq, _) = haplotype_consensuses.get(re_fragment, Haplotype::Unspecified);

    //     // use partial order alignment on original strand sequences to establish
    //     // the haplotype consensus
    //     let mut strands: Vec<&[u8]> = Vec::with_capacity(read_is.len() * 2);
    //     for read_i in read_is {
    //         let qname = &reads[*read_i].qname;
    //         if let Some(seqs) = source_strands.by_read.get(qname){
    //             if !seqs.1.is_empty() {
    //                 strands.push(&seqs.0);
    //                 strands.push(&seqs.1);                    
    //             }
    //         }   
    //     }
    //     if strands.len() < 4 { return None; }
    //     let haplotype_consensus = consensus(&strands, 0, &self.poa_config)
    //         .expect("failed to generate POA consensus");

    //     // align the fragment reference span to the haplotype consensus, i.e., ref_on_hap
    //     let mm2_hap = self.minimap2.clone()
    //         .with_seq(&haplotype_consensus.sequence)
    //         .expect("Failed to initialize minimap2 in build_haplotype_consensus()");
    //     let ref_on_haps = &mm2_hap.map(
    //         ref_seq.as_bytes(), 
    //         true, 
    //         false, 
    //         None, 
    //         Some(&MM_F_NO_PRINT_2ND), 
    //         None
    //     ).expect("Minimap2 failed at ref_on_hap");
    //     if ref_on_haps.len() == 0 { return None; }
    //     let ref_on_hap = &ref_on_haps[0];
    //     let Some(aln) = &ref_on_hap.alignment else { return None; };
    //     let Some(cs) = &aln.cs else { return None; }; // never expected to fail

    //     // create a map of reference positions per each haplotype consensus position
    //     // and a hap_vs_ref encoding to retain a memory of clonal variants relative to reference
    //     ref_pos0_map.clear();
    //     self.hap_vs_ref.clear();
    //     let hap_vs_ref = self.process_cs_ref_on_hap(
    //         ref_pos0_map, re_fragment, haplotype, 
    //         ref_on_hap.target_start, ref_on_hap.target_end, ref_on_hap.query_start, 
    //         haplotype_consensus.sequence.len(), cs, 
    //         reads, read_is,
    //     );

    //     // cache the haplotype consensus for printing
    //     haplotype_consensuses.insert(
    //         re_fragment, haplotype, 
    //         unsafe { from_utf8_unchecked(&haplotype_consensus.sequence).to_string() }, 
    //         Some(hap_vs_ref)
    //     );

    //     // return the built minimap2 for subsequent read alignment
    //     Some((haplotype_consensus.sequence, mm2_hap))
    // }
}
