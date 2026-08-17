//! Support for parsing multiple types of cs tags (read_on_ref, read_on_hap,
//! etc.) for identifying and calling variants.

// imports
use super::*;

impl SnvChromWorker {

    /// Process a cs:Z:tag to add a read to a growing fragment variant list. 
    pub fn process_cs_tag(
        &mut self,
        frag_haps:    &mut FragmentHaplotypes,
        re_fragment:  &ReFragment,
        haplotype:    Haplotype,
        mapping:      Option<(i32, i32)>, // minimap2 query_start, target_start
        cs_tag:       Option<&String>,
        ref_pos0_map: &mut Vec<ChromPos0>, // used different depending on tgt_is_hap bool
        read_i:       ReadIndex,
        read:         &ReadInstance,
        source_strands: &(SourceStrand, SourceStrand),
        str0_pos0_map:  &Vec<usize>,
        str1_pos0_map:  &Vec<usize>,
    ) { 
        let is_read_on_hap = haplotype != Haplotype::Unspecified;
        let (
            mut qry_pos0, // position on the query read on same strand as tgt_pos0
            mut tgt_pos0, // position on either chromosome or haplotype consensus
        ) = if is_read_on_hap {(
            mapping.unwrap().0 as u32, // properly oriented since all qry are ref top strand
            mapping.unwrap().1 as u32,
        )} else {( // read_on_ref
            read.qry_pos0,
            read.aln_start0,
        )};

        self.reset_cs_variant();

        let mut chars = if is_read_on_hap {
            cs_tag.unwrap().chars()
        } else {
            read.cs.chars()
        };
        self.cs_op = chars.next().unwrap();
        self.op_val.clear();
        
        while let Some(char) = chars.next() {
            if char.is_alphanumeric() {
                self.op_val.push(char);
            } else {
                self.handle_cs_op(
                    frag_haps, ref_pos0_map,
                    re_fragment, haplotype, is_read_on_hap, 
                    &mut qry_pos0, &mut tgt_pos0, read_i, 
                    source_strands, str0_pos0_map, str1_pos0_map
                );
                self.cs_op = char;
                self.op_val.clear();
            }
        }
        self.handle_cs_op(
            frag_haps, ref_pos0_map,
            re_fragment, haplotype, is_read_on_hap, 
            &mut qry_pos0, &mut tgt_pos0, read_i, 
            source_strands, str0_pos0_map, str1_pos0_map
        );
    }

    /// Process one cs:Z:tag operation to add to the growing fragment variant 
    /// list. 
    fn handle_cs_op(
        &mut self,
        frag_haps:    &mut FragmentHaplotypes,
        ref_pos0_map: &mut Vec<ChromPos0>,
        re_fragment:  &ReFragment,
        haplotype:    Haplotype,
        is_read_on_hap: bool,
        qry_pos0:     &mut SeqPos0, 
        tgt_pos0:     &mut SeqPos0,
        read_i:       ReadIndex,
        source_strands: &(SourceStrand, SourceStrand),
        str0_pos0_map:  &Vec<usize>,
        str1_pos0_map:  &Vec<usize>,
    ) {
        match self.cs_op {

            // :[0-9]+   Identical sequence length
            ':' => {
                if self.var_tgt_pos0.is_some() { // commit any preceding variant stretch
                    if self.allowed {
                        let tgt_pos0 = self.var_tgt_pos0.unwrap();
                        let ref_pos0 = if is_read_on_hap {
                            ref_pos0_map[tgt_pos0 as usize]
                        } else {
                            tgt_pos0
                        };
                        let variant = Variant::new(
                            ref_pos0,
                            tgt_pos0,
                            &self.tgt_bases,
                            &self.alt_bases,
                            re_fragment, haplotype
                        );
                        self.frag_vars.insert(variant.clone(), read_i, self.min_qual);
                        frag_haps.insert_variant(re_fragment, haplotype, variant);
                    }
                    self.reset_cs_variant();
                }
                let len = self.op_val.parse::<u32>().unwrap();
                self.encoding.add_identity(len);
                if is_read_on_hap {
                    *qry_pos0 += len;
                    *tgt_pos0 += len;
                } else {
                    for _ in 0..len {
                        ref_pos0_map.push(*qry_pos0);
                        *qry_pos0 += 1;
                        *tgt_pos0 += 1;
                    }
                }
            },

            // *[acgtn][acgtn]   Substitution: target to query
            '*' => {
                //     S
                // rrrrRrrrr
                // qqqqQqqqq
                //     A
                let tgt = self.op_val[0..=0].to_ascii_uppercase();
                let alt = self.op_val[1..=1].to_ascii_uppercase();
                self.tgt_bases.push_str(&tgt);
                self.alt_bases.push_str(&alt);
                self.allowed &= tgt != "N";
                self.allowed &= alt != "N";
                self.allowed &= if is_read_on_hap {
                    let ref_pos0 = ref_pos0_map[*tgt_pos0 as usize];
                    !self.simple_repeats.binary_search(ref_pos0, 1)
                } else {
                    ref_pos0_map.push(*qry_pos0);
                    !self.simple_repeats.binary_search(*tgt_pos0, 1)
                };
                let low_qual = if is_read_on_hap {
                    let i0 = str0_pos0_map[*tgt_pos0 as usize];
                    let mut min_qual = source_strands.0.get_min_qual(i0, i0 + 1);
                    let i0 = str1_pos0_map[*tgt_pos0 as usize];
                    min_qual = min_qual.min(source_strands.1.get_min_qual(i0, i0 + 1));
                    self.min_qual = self.min_qual.min(min_qual);
                    min_qual <= MIN_SNV_INDEL_QUAL
                } else {
                    false
                };
                self.encoding.add_substitution(&alt, self.allowed, low_qual);
                if self.var_tgt_pos0.is_none() { self.var_tgt_pos0 = Some(*tgt_pos0); }
                *qry_pos0 += 1;
                *tgt_pos0 += 1;
            },

            // +[acgtn]+   Insertion to the target
            '+' => {
                //    *INI     insertions may have heteroduplex bases within homoduplex query run
                // rrrr   Rrrr
                // qqqqQqqqqqq
                //    aA Aa
                let n_ins_bases = self.op_val.len();
                let alt = self.op_val.to_ascii_uppercase();
                self.alt_bases.push_str(&alt);
                self.allowed &= !alt.contains("N");
                self.allowed &= if is_read_on_hap {
                    let ref_pos0 = ref_pos0_map[*tgt_pos0 as usize - 1];
                    !self.simple_repeats.binary_search(ref_pos0, 2)
                } else {
                    !self.simple_repeats.binary_search(*tgt_pos0 - 1, 2)
                };
                let low_qual = if is_read_on_hap {
                    let start0 = str0_pos0_map[*tgt_pos0 as usize - 1];
                    let end1   = str0_pos0_map[*tgt_pos0 as usize] + 1;
                    let mut min_qual = source_strands.0.get_min_qual(start0, end1);
                    let start0 = str1_pos0_map[*tgt_pos0 as usize - 1];
                    let end1   = str1_pos0_map[*tgt_pos0 as usize] + 1;
                    min_qual = min_qual.min(source_strands.1.get_min_qual(start0, end1));
                    self.min_qual = self.min_qual.min(min_qual);
                    min_qual <= MIN_SNV_INDEL_QUAL
                } else {
                    false
                };
                self.encoding.add_insertion(self.allowed, low_qual);
                if self.var_tgt_pos0.is_none() { self.var_tgt_pos0 = Some(*tgt_pos0 - 1); }
                *qry_pos0 += n_ins_bases as u32;
                // no action on tgt_bases or tgt_pos0
            },

            // -[acgtn]+   Deletion from the target
            '-' => {
                //     DDD
                // rrrrRrrrrrr
                // qqqq   Qqqq
                //   aA   Aa
                // heteroduplex indels in read strands are always reported as N bases
                // so do not expect heteroduplex indels to lead to falsely missing bases
                let n_del_bases = self.op_val.len() as u32;
                let tgt = self.op_val.to_ascii_uppercase();
                self.tgt_bases.push_str(&tgt);
                self.allowed &= !tgt.contains("N");
                self.allowed &= if is_read_on_hap {
                    let ref_pos0 = ref_pos0_map[*tgt_pos0 as usize];
                    !self.simple_repeats.binary_search(ref_pos0, n_del_bases)
                } else {
                    for _ in 0..n_del_bases {
                        ref_pos0_map.push(*qry_pos0 - 1); // never used downstream
                    }
                    !self.simple_repeats.binary_search(*tgt_pos0, n_del_bases)
                };
                let low_qual = if is_read_on_hap {
                    let start0 = str0_pos0_map[*tgt_pos0 as usize];
                    let mut min_qual = source_strands.0.get_min_qual(start0, start0 + 2);
                    let start0 = str1_pos0_map[*tgt_pos0 as usize];
                    min_qual = min_qual.min(source_strands.1.get_min_qual(start0, start0 + 2));  
                    self.min_qual = self.min_qual.min(min_qual);
                    min_qual <= MIN_SNV_INDEL_QUAL
                } else {
                    false
                };
                self.encoding.add_deletion(n_del_bases, self.allowed, low_qual);
                if self.var_tgt_pos0.is_none() { self.var_tgt_pos0 = Some(*tgt_pos0); }
                *tgt_pos0    += n_del_bases;
                // no action on qry_pos0, alt_bases, and N check not applicable
            },
            _   => panic!("Unexpected operation in cs tag: {}", self.cs_op),
        }
    }

    // /// Convert a Smith-Waterman Alignment into the equivalent minimap2 cs tag.
    // /// TODO: move this to genomex crate.
    // pub fn get_cs_tag(
    //     // &self,
    //     aln: &Alignment,
    //     tgt: &str, // the target sequenced that generated the alignment
    // ) -> String { 
    //     //    M operations carry the query base in the array slot (could be a base mismatch)
    //     //    I operations carry the inserted base prepended to the NEXT target postion
    //     //    D operations carry "-" in place of the query base that was deleted relative to target
    //     let mut cs = String::with_capacity(256);
    //     let mut del_val: String = String::with_capacity(256);
    //     let mut identity_len = 0_usize;
    //     for tgt_i0 in aln.tgt_start0..=aln.tgt_end0 {
    //         let tgt_base = &tgt[tgt_i0..=tgt_i0];
    //         let aln_val = aln.qry_on_tgt[tgt_i0 - aln.tgt_start0].as_str();
    //         if tgt_base == aln_val {
    //             if del_val.len() > 0 {
    //                 cs.push_str(&format!("-{}", del_val.to_ascii_lowercase()));
    //                 del_val.clear();
    //             }
    //             identity_len += 1;
    //         } else {
    //             if identity_len > 0 { 
    //                 cs.push_str(&format!(":{identity_len}"));
    //                 identity_len = 0;
    //             }
    //             if aln_val == "-" {
    //                 del_val.push_str(&tgt_base);
    //             } else {
    //                 if del_val.len() > 0 {
    //                     cs.push_str(&format!("-{}", del_val.to_ascii_lowercase()));
    //                     del_val.clear();
    //                 }
    //                 if aln_val.len() > 1 {
    //                     let ins_bases = &aln_val[0..aln_val.len() - 1];
    //                     cs.push_str(&format!("+{}", ins_bases.to_ascii_lowercase()));
    //                     identity_len = 1;
    //                 } else {
    //                     cs.push_str(&format!("*{}{}", tgt_base.to_ascii_lowercase(), aln_val.to_ascii_lowercase()));
    //                 }
    //             }
    //         }

    //     }
    //     if identity_len > 0 { cs.push_str(&format!(":{identity_len}")) }
    //     cs
    // }
}
