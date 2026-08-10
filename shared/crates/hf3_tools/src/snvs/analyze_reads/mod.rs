//! Support for building haplotype consensuses from RE fragments and using them
//! to call rare variants without reference alignment bias.
//! 
//! Function `analyze_reads` in this script consumes the most time in 
//! `analyze SNVs` action execution.
//! 
//! Note that all sequences are reference top-strand oriented throughout
//! even if a rare basecalled-read was originally in the reverse orientiation.

// modules
mod build_consensus;
mod align_reads;

// imports
use super::*;
use align_reads::*;

// constants
const MIN_HAPLOTYPE_READS: u16 = 2;  // a haplotype must have >=2 matching reads to be called heterozygous
const MM_F_NO_PRINT_2ND: [u64; 1] = [16384]; // minimap2 flag to suppress return of secondary alignments

impl SnvChromWorker {

    /// Reset a FragmentAnalysisTool for each newly encountered ReFragment.
    pub fn reset_read_analysis(&mut self, n_reads: usize){
        self.frag_vars.reset(n_reads); 
        self.tracking_variants.clear();
        self.hap_vars.get_mut(&Haplotype::Haplotype1).unwrap().clear();
        self.hap_vars.get_mut(&Haplotype::Haplotype2).unwrap().clear();
        self.hap_vars.get_mut(&Haplotype::Homozygous).unwrap().clear();
    }

    /// Reset FragmentAnalysisTool fields used during cs tag variant processing.
    pub fn reset_cs_variant(&mut self){
        self.var_tgt_pos0 = None;
        self.tgt_bases.clear();
        self.alt_bases.clear();
        self.alt_qual.clear();
        self.allowed      = true; 
    }    

    /// Reset read haplotype voting.
    pub fn reset_hap_votes(&mut self){
        *self.hap_votes.get_mut(&Haplotype::Haplotype1).unwrap() = 0;
        *self.hap_votes.get_mut(&Haplotype::Haplotype2).unwrap() = 0;
    }    

    /// Build homozygous or heterozygous consensuses for RE fragment haplotypes,
    /// and use them to call clonal and subclonal variants.
    pub fn analyze_reads(
        &mut self,
        tool: &SnvAnalysisTool,
        fragment_reads: FragmentReads,
        haplotype_consensuses: &mut HaplotypeConsensuses,
        reads_on_reference:    &mut FragmentHaplotypes,
        reads_on_haplotype:    &mut FragmentHaplotypes,
    ){
        // allocate recycled objects used in fragment read parsing
        let mut ref_pos_maps: Vec<Vec<ChromPos0>> = Vec::with_capacity(128); // ref pos to read pos per read
        let mut read_masks: Vec<Vec<DdMaskType>> = Vec::with_capacity(128); // expanded dd:Z: tags per read
        let mut ref_pos0_map1: Vec<ChromPos0> = Vec::with_capacity(MAX_EXPECTED_READ_LEN); // consensus to ref pos maps
        let mut ref_pos0_map2: Vec<ChromPos0> = Vec::with_capacity(MAX_EXPECTED_READ_LEN);

        // process all observed putative ReFragments one at a time
        for re_fragment in fragment_reads.instances.keys() {
            // self.show_debug = *re_fragment == self.debug;
            // if !self.show_debug { continue; }

            // require a minimum ReFragment coverage to continue
            let reads = fragment_reads.instances.get(&re_fragment).unwrap();
            let n_reads = reads.len();
            if n_reads < self.min_fragment_reads { continue; }

            // cache the RE fragment reference sequence for use below and during visualization
            haplotype_consensuses.insert_reference(tool, self, re_fragment);

            // collect all variants over all fragment reads into a variant-by-read matrix
            self.reset_read_analysis(n_reads);
            ref_pos_maps.clear(); // maps start at the first aligned read base
            read_masks.clear();   // masks start at the first read base, even if clipped
            reads.iter().enumerate().for_each(|(read_i, read)|{
                self.encoding.prepare_read_on_ref(re_fragment, read);
                ref_pos_maps.push(Vec::new());
                read_masks.push(Self::get_dd_mask(read));
                self.process_cs_tag(
                    reads_on_reference, re_fragment, Haplotype::Unspecified, 
                    None, None, 
                    &mut ref_pos_maps[read_i],
                    read_i, read, &read_masks[read_i]
                );
                reads_on_reference.insert_encoding(
                    re_fragment, 
                    Haplotype::Unspecified,
                    self.encoding.clone()
                );
            });

            // determine which reads could have reported each variant, i.e., were informative
            for (variant, vmap) in &mut self.frag_vars.variant_map {
                // only substitutions and insertions can be declared uninformative
                // simple deletions had no read bases to report as N values
                if let Some(alt_bases) = &variant.alt_bases {
                    let n_alt_bases = alt_bases.len();
                    reads.iter().enumerate().for_each(|(read_i, read)|{
                        let map_pos0 = variant.ref_pos0.saturating_sub(read.aln_start0) as usize;
                        if let Some(read_pos0) = ref_pos_maps[read_i].get(map_pos0) {
                            let read_pos0 = if variant.tgt_bases.is_none() {
                                *read_pos0 + 1
                            } else {
                                *read_pos0
                            } as usize;
                            let max1 = (read_pos0 + n_alt_bases).min(read.seq_bytes.len());
                            if read.seq_bytes[read_pos0..max1].contains(&b'N') ||
                               read_masks[read_i][read_pos0..max1].contains(&DdMaskType::CorrectedToReference) {
                                vmap.read_map[read_i].is_informative = false;
                            } else {
                                vmap.n_informative += 1;
                            }
                        } else {
                            vmap.n_informative += 1;
                        }
                    });
                } else {
                    vmap.n_informative = n_reads as u16;
                }
            }

            // remove homozyogous variants as they do not help parse haplotypes
            // collect tracking variants, i.e., recurrent allowed SNVs or indels
            let variants: Vec<Variant> = self.frag_vars.variant_map.keys().cloned().collect();
            let mut has_tracking_snp = false;
            for variant in variants {
                let vmap = self.frag_vars.variant_map.get_mut(&variant).unwrap();
                if vmap.n_matching_reads >= vmap.n_informative { // no informative alleles lacked this variant
                    self.frag_vars.variant_map.remove(&variant);
                } else {
                    let vaf = vmap.n_matching_reads as f64 / vmap.n_informative as f64;
                    vmap.zyg_int = ((vaf - 0.5).abs() * 200.0) as u8;
                    if vmap.n_matching_reads >= MIN_HAPLOTYPE_READS {
                        has_tracking_snp |= !variant.is_indel;
                        self.tracking_variants.push(variant)
                    }
                }
            }

            // short-circuit to final processing of purely homozygous fragments
            // with no potential remaining heterozygous or subclonal tracking variants
            let n_tracking_variants = self.tracking_variants.len();
            if n_tracking_variants == 0 {
                // require a minimum read count to avoid false homozygosity
                if n_reads < self.min_homozygous_reads { continue; }
                let hap1_read_is: Vec<ReadIndex> = (0..n_reads).map(|i| i).collect();
                let mm2_hap1 = self.build_haplotype_consensus(
                    haplotype_consensuses, re_fragment, Haplotype::Homozygous, 
                    &mut ref_pos0_map1, 
                    &reads, &hap1_read_is
                );
                let Some(mm2_hap1) = mm2_hap1 else { continue; };
                self.align_to_haplotype_consensus(
                    reads_on_haplotype, re_fragment, Haplotype::Homozygous,
                    &mut ref_pos0_map1, 
                    reads, hap1_read_is, &read_masks,
                    Some(mm2_hap1), None
                );
                continue;
            }

            // restrict tracking variants to more reliable SNPs if any are available
            // but use indels for haplotype assignment if they are all that exists
            // noting that tracking variants are always allowed, i.e., not in simple repeats
            if has_tracking_snp {
                self.tracking_variants.retain(|variant| !variant.is_indel);
            }

            // search for the variant to use as the index for defining haplotype consensuses
            // prefer the variant with vaf nearest to 0.5
            // let index_var = self.tracking_variants.iter()
            //     .min_by_key(|&variant| {
            //         let vmap =  &self.frag_vars.variant_map[variant];
            //         vmap.zyg_int // smallest values are the closest to vaf==0.5
            //     })
            //     .unwrap();
            // let index_vmap = &self.frag_vars.variant_map[&index_var];

            // find the most likely zygosity among the tracking variants by majority vote
            // when two zygosities are equally frequent, choose the vaf nearest to 0.5 (zyg_int closest to zero)
            let mut tracking_zyg_ints: Vec<_> = self.tracking_variants.iter()
                .map(|variant|{
                    let vmap =  &self.frag_vars.variant_map[variant];
                    vmap.zyg_int // smallest values are the closest to vaf==0.5
                }).collect();
            tracking_zyg_ints.sort_unstable();
            let mut tracking_zyg_ints: Vec<(isize, u8)> = tracking_zyg_ints
                    .chunk_by(|a, b| a == b)
                    .map(|chunk| (-(chunk.len() as isize), chunk[0]))
                    .collect();
            tracking_zyg_ints.sort_unstable();

            // select a single index variant with the most frequently observed zygosity
            // it does not have to be perfect, just sufficient to build consensuses
            let index_zyg_int = tracking_zyg_ints[0].1;
            let index_variant = self.tracking_variants.iter()
                .find(|&variant| {
                    let vmap =  &self.frag_vars.variant_map[variant];
                    vmap.zyg_int == index_zyg_int
                })
                .unwrap();
            let index_vmap = &self.frag_vars.variant_map[&index_variant];
            // if self.show_debug {
            //     for variant in &self.tracking_variants {
            //         eprintln!("tracking_variant {:?}", variant);
            //         let vmap = self.frag_vars.variant_map.get(&variant).unwrap();
            //         eprintln!("n_matching_reads {:?}", vmap.n_matching_reads);
            //     }
            //     eprintln!("index_var {:?}", index_variant);
            //     eprintln!("index_vmap.n_matching_reads {}", index_vmap.n_matching_reads);
            // }

            // find the reads that contributed to each index haplotype
            // such reads must be informative to avoid false negative variant calls
            let mut hap1_read_is: Vec<ReadIndex> = (0..n_reads)
                .filter_map(|read_i|{
                    if index_vmap.read_map[read_i].has_var() &&
                       index_vmap.read_map[read_i].is_informative { 
                        Some(read_i) 
                    } else { None }
                }).collect();
            let mut hap2_read_is: Vec<ReadIndex> = (0..n_reads)
                .filter_map(|read_i|{
                    if !index_vmap.read_map[read_i].has_var() &&
                        index_vmap.read_map[read_i].is_informative { 
                        Some(read_i) 
                    } else { None }
                }).collect();
            if hap1_read_is.len() == 0 || 
               hap2_read_is.len() == 0 { 
                continue;
            }

            // build consensus sequences of each initial haplotype
            let mm2_hap1 = self.build_haplotype_consensus(
                haplotype_consensuses, re_fragment, Haplotype::Haplotype1, 
                &mut ref_pos0_map1, 
                &reads, &hap1_read_is
            );
            let mm2_hap2 = self.build_haplotype_consensus(
                haplotype_consensuses, re_fragment, Haplotype::Haplotype2, 
                &mut ref_pos0_map2, 
                &reads, &hap2_read_is
            );
            let Some(mm2_hap1) = mm2_hap1 else { continue; };
            let Some(mm2_hap2) = mm2_hap2 else { continue; };

            // assign all reads to their final haplotype by comparing 
            // read_on_ref variants to ref_on_hap variants
            hap1_read_is.clear();
            hap2_read_is.clear();
            let mut hap1_assignments: Vec<ReadAssignment> = Vec::new();
            let mut hap2_assignments: Vec<ReadAssignment> = Vec::new();
            (0..n_reads).for_each(|read_i| {
                if let Some(assignment) = self.assign_read_to_haplotype(
                    &reads[read_i], read_i, &mm2_hap1, &mm2_hap2
                ){
                    match assignment.haplotype {
                        Haplotype::Haplotype1 => {
                            hap1_read_is.push(read_i);
                            hap1_assignments.push(assignment);
                        },
                        Haplotype::Haplotype2 =>  {
                            hap2_read_is.push(read_i);
                            hap2_assignments.push(assignment);
                        },
                        _ => {}
                    }
                }
            });

            // align each read to its haplotype consensus to call subclonal variants
            if hap1_read_is.len() > 0 {
                self.align_to_haplotype_consensus(
                    reads_on_haplotype, re_fragment, Haplotype::Haplotype1,
                    &mut ref_pos0_map1, 
                    reads, hap1_read_is, &read_masks,
                    None, Some(hap1_assignments)
                );
            }
            if hap2_read_is.len() > 0 {
                self.align_to_haplotype_consensus(
                    reads_on_haplotype, re_fragment, Haplotype::Haplotype2,
                    &mut ref_pos0_map2, 
                    reads, hap2_read_is, &read_masks,
                    None, Some(hap2_assignments)
                );
            }
        }
    }

}
