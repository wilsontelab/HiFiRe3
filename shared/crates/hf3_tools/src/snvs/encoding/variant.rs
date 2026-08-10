//! Support for creating encoded variant-level files for downstream use.

// imports
use rustc_hash::{FxHashMap, FxHashSet};
use serde::{Serialize, Serializer};
use mdi::OutputCsv;
use crate::snvs::*;

/// Clonality lists the types of variant calls. Unlike encodings, a single
/// output file includes all variants calls.
#[derive(PartialEq, Eq, Clone, Copy, Serialize)]
#[repr(u8)]
pub enum Clonality {
    Clonal    = 1,
    Subclonal = 0,
}
/// Helper function to serialize VariantZygosity as u8.
pub fn serialize_clonality<S: Serializer>(
    c: &Clonality, 
    serializer: S
) -> Result<S::Ok, S::Error>{
    serializer.serialize_u8(*c as u8)
}

/// VariantInstances holds the read count of a specific Variant and the 
/// fragments, samples, and reads that contributed to the count.
pub struct VariantInstances {
    n_matching_reads:  u16, // set for subclonal only; always 0 for clonal variants
    n_haplotype_reads: u16,
    n_reads:           u16,
    sample_bits:       SampleBits,
    max_avg_qual:      PhredQual,
    clonal:    Clonality,
    qnames:    Vec<QName>,
}
impl VariantInstances {

    /// Create a new empty VariantInstances object.
    fn new() -> Self {
        VariantInstances {
            n_matching_reads:  0,
            n_haplotype_reads: 0,
            n_reads:           0,
            sample_bits:       0,
            max_avg_qual:      0,
            clonal:    Clonality::Subclonal,
            qnames:    Vec::new(), // auto-allocate since many variants will have few reads
        }
    }
}

/// VariantMetadata reports summary results of variant calling and counting,
/// excluding variants in simple repeats.
pub struct VariantMetadata {
    pub n_variants:       usize,
    pub n_substitutions:  usize,
    pub n_insertions:     usize,
    pub n_deletions:      usize,
    pub n_clonal:         usize,
    pub n_subclonal:      usize,
    pub variant_count:    usize,
    pub variant_coverage: usize,
}
impl VariantMetadata {
    fn new(n_variants: usize) -> Self {
        VariantMetadata {
            n_variants,
            n_substitutions:  0,
            n_insertions:     0,
            n_deletions:      0,
            n_clonal:         0,
            n_subclonal:      0,
            variant_count:    0,
            variant_coverage: 0,
        }
    }
}

/// A VariantRecord is a specific SNV or indel as written to file for 
/// downstream use. 
#[derive(Serialize)]
struct VariantRecord<'a> {
    chrom_index:       ChromIndex1,
    variant:           &'a Variant, // specific variant observed at this position
    n_repeat_bases:    u32,
    context_base_left: UppercaseACGTN,
    context_base_right:UppercaseACGTN,
    n_matching_reads:  u16,
    n_haplotype_reads: u16,
    n_reads:           u16,
    n_multivariant_reads: u16, // number of reads that had at least one other subclonal variant
    sample_bits:      SampleBits,
    n_samples:        u32,
    #[serde(serialize_with = "serialize_clonality")]
    clonal:           Clonality,
    matches_clonal:   u8,
    max_avg_qual:     PhredQual, // max_avg_qual set on subclonal
    qnames: CommaDelimited, // comma-delimited list of QNAMEs with this variant
}

/// VariantsTally aggregates accumulated VariantInstances per Variant.
pub struct VariantsTally {
    pub tally:  FxHashMap<Variant, VariantInstances>,
    pub clonal: FxHashSet<VariantLocation>
}
impl VariantsTally {

    /// Create a new empty VariantsTally object. On tally is instantiated per 
    /// SnvChromWorker.
    pub fn new() -> Self {
        let mut tally = FxHashMap::default();
        tally.reserve(1_048_576);
        let mut clonal = FxHashSet::default();
        clonal.reserve(1_048_576);
        Self{ tally, clonal }
    }

    /// Update the instances tally of a specific clonal Variant derived from
    /// alignment to reference.
    pub fn add_clonal(
        &mut self, 
        variant:   &Variant, 
        reads:     &[ReadInstance], // all ReFragment reads
        read_is:   &[ReadIndex],    // indices into reads for reads with the variant
    ) {
        let instances = self.tally
            .entry(variant.clone())
            .or_insert_with(VariantInstances::new);
        instances.n_haplotype_reads += read_is.len() as u16;
        instances.n_reads           += reads.len() as u16;
        instances.clonal    = Clonality::Clonal;
        for read_i in read_is {
            let read = &reads[*read_i];
            instances.sample_bits |= read.sample_bit;
            instances.qnames.push(read.qname.clone());
        }
        self.clonal.insert(VariantLocation { 
            ref_pos0:    variant.ref_pos0, 
            is_indel:    variant.is_indel, 
            re_fragment: variant.re_fragment, 
        });
    }

    /// Update the instances tally of a specific subclonal Variant derived from
    /// alignment to a haplotype consensus.
    pub fn add_subclonal(
        &mut self, 
        variant: &Variant, 
        reads:   &[ReadInstance], // all ReFragment reads
        read_is: &[ReadIndex],    // read indices for the haplotype being processed
        read_js: &[ReadIndex],    // indices into read_is for reads with the variant
        max_avg_qual: PhredQual,
    ) {
        let instances = self.tally
            .entry(variant.clone())
            .or_insert_with(VariantInstances::new);
        instances.n_matching_reads  += read_js.len() as u16;
        instances.n_haplotype_reads += read_is.len() as u16;
        instances.n_reads           += reads.len() as u16;
        instances.max_avg_qual = instances.max_avg_qual.max(max_avg_qual);
        instances.clonal       = Clonality::Subclonal;
        for read_j in read_js {
            let read = &reads[read_is[*read_j]];
            instances.sample_bits |= read.sample_bit;
            instances.qnames.push(read.qname.clone());
        }
    }

    /// Sort and write a set of VariantInstances to a temporary file for the
    /// working chromosome.
    pub fn write_sorted(
        tool:   &SnvAnalysisTool,
        worker: &mut SnvChromWorker,
        haplotype_consensuses: &mut HaplotypeConsensuses,
        file_path: String,
    ) -> VariantMetadata {
        let mut csv = OutputCsv::open_csv(
            &file_path, 
            b'\t', 
            false, 
            Some(tool.n_cpu),
        );
        let mut variants = worker.variant_tally.tally.keys()
            .filter_map(|v|{
                let excluded  =  tool.exclusions.pos_in_region(&worker.chrom, v.ref_pos0 + 1);
                let on_target = !tool.targets.has_data || 
                                       tool.targets.pos_in_region(&worker.chrom, v.ref_pos0 + 1);
                if !excluded && on_target { Some(v.clone()) } else { None }
            }).collect::<Vec<_>>();
        variants.sort_unstable();
        let mut md = VariantMetadata::new(variants.len());
        for variant in variants {
            let instances = &worker.variant_tally.tally[&variant];
            let (context_base_left, context_base_right) = if instances.clonal == Clonality::Clonal {
                ("NA".to_string(), "NA".to_string())
            } else {
                let (hap_seq, _) = haplotype_consensuses.cache
                    .get_mut(&(variant.re_fragment, variant.haplotype))
                    .expect("Failed to get haplotype seq from cache.");
                let n_tgt_bases = variant.tgt_bases.as_ref().map_or(0, |s| s.len());
                let (left_pos0, right_pos0) = if n_tgt_bases > 0 {
                    (
                        variant.tgt_pos0.saturating_sub(1) as usize,
                        variant.tgt_pos0 as usize + n_tgt_bases
                    )
                } else {
                    (
                        variant.tgt_pos0 as usize,
                        variant.tgt_pos0 as usize + 1                    
                    )
                };
                let right_pos0 = right_pos0.min(hap_seq.len() - 1);
                (
                    hap_seq[left_pos0..=left_pos0].to_string(),
                    hap_seq[right_pos0..=right_pos0].to_string()
                )
            };
            let record = VariantRecord {
                chrom_index:   worker.chrom_index,
                variant:       &variant,
                n_repeat_bases: worker.simple_repeats.get_n_repeat_bases(&variant.re_fragment),
                context_base_left,
                context_base_right,
                n_matching_reads:     instances.n_matching_reads,
                n_haplotype_reads:    instances.n_haplotype_reads,
                n_reads:              instances.n_reads,
                n_multivariant_reads: instances.qnames.iter().filter(|&qname|{
                    worker.variant_reads_tally.tally.get(qname)
                        .map_or(false, |read|{
                            read.n_variants > 1
                        })
                }).count() as u16,
                sample_bits:   instances.sample_bits,
                n_samples:     instances.sample_bits.count_ones(),
                clonal:        instances.clonal,
                matches_clonal: if instances.clonal == Clonality::Clonal { 0 } else {
                    worker.variant_tally.clonal.contains(&VariantLocation { 
                        ref_pos0:    variant.ref_pos0, 
                        is_indel:    variant.is_indel, 
                        re_fragment: variant.re_fragment, // on either haplotype
                    }) as u8
                },
                max_avg_qual:  instances.max_avg_qual,
                qnames:   instances.qnames.join(",")
            };
            csv.serialize(&record);

            let alt_minus_ref = variant.alt_minus_ref();
            if alt_minus_ref == 0 {
                md.n_substitutions += 1;
            } else if alt_minus_ref > 0 {
                md.n_insertions += 1;
            } else {
                md.n_deletions += 1;
            }
            match record.clonal {
                Clonality::Subclonal => {
                    md.n_subclonal += 1
                },
                Clonality::Clonal => {
                    md.n_clonal += 1
                }
            }
            md.variant_count    += record.n_matching_reads as usize;
            md.variant_coverage += record.n_reads as usize; 
        }
        md
    }
}
