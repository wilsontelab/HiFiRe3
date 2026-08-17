//! Support for creating a qname, i.e., read-level summary file across all 
//! samples. Only reads with at least one allowed subclonal variant are 
//! recorded, and only allowed subclonal variants are tallied.

// imports
use rustc_hash::FxHashMap;
use mdi::OutputCsv;
use crate::snvs::*;

/// VariantReadInstance records metadata about all allowed subclonal variants 
/// found in a single read with at least one such variant. 
pub struct VariantReadInstance {
    re_fragment:  ReFragment,
    haplotype:    Haplotype,
    sample_bit:   SampleBit,
    n_bases:           u32, // number of bases in read, as compared to re_fragment size
    n_haplotype_reads: u16, // number of reads in this read's haplotype
    n_reads:           u16, // number of fragment reads in both haplotypes
    pub n_variants:    u16,
    n_low_qual:   u16,
    n_snv:        u16, // single-nucleotide variant
    n_mnv:        u16, // equal-length multi-nucleotide variant
    n_del:        u16, // simple deletion relative to haplotype
    n_ins:        u16, // simple insertion relative to haplotype
    n_complex:    u16, // a base span replaced with variant span of unequal length
    variants:     Vec<String>,
}
impl VariantReadInstance {
    /// Create a new empty VariantReadInstance object. One instance is 
    /// instantiated per recorded read.
    pub fn new(
        read:        &ReadInstance,
        re_fragment: &ReFragment,
        haplotype:   &Haplotype,
    ) -> Self {
        Self{ 
            re_fragment:  *re_fragment,
            haplotype:    *haplotype,
            sample_bit:   read.sample_bit,
            n_bases:      read.seq_bytes.len() as u32,
            n_haplotype_reads:  0,
            n_reads:      0,
            n_variants:   0,
            n_low_qual:   0,
            n_snv:        0, // single-nucleotide variant
            n_mnv:        0, // equal-length multi-nucleotide variant
            n_del:        0, // simple deletion relative to haplotype
            n_ins:        0, // simple insertion relative to haplotype
            n_complex:    0, // a base span replaced with variant span of unequal length
            variants:     Vec::new(),
        }
    }
}

/// VariantReadsMetadata reports summary results of reads with subclonal variants.
pub struct VariantReadsMetadata {
    pub n_variant_reads: usize,
    pub n_indel_only:    usize,
    pub n_one_snv:       usize,
    pub n_two_snv:       usize,
    pub n_three_snv:     usize,
    pub n_four_snv:      usize,
    pub n_five_snv:      usize, // five or more SNVs, exlusive of indels
}
impl VariantReadsMetadata {
    fn new(n_variant_reads: usize) -> Self {
        VariantReadsMetadata {
            n_variant_reads,
            n_indel_only: 0,
            n_one_snv:    0,
            n_two_snv:    0,
            n_three_snv:  0,
            n_four_snv:   0,
            n_five_snv:   0,
        }
    }
}

/// A fully assembled set of metadata about a VariantReadInstance as printed to 
/// file. 
#[derive(Serialize)]
struct VariantReadRecord {
    qname:       QName, // functionally a BED file, even if qname is not the proper ref
    chrom_index: ChromIndex1,
    re_fragment: ReFragment,
    #[serde(serialize_with = "serialize_haplotype")]
    haplotype:   Haplotype,
    n_repeat_bases: u32,
    sample_bit:  SampleBit,
    n_bases:     u32, // number of bases in read, as compared to re_fragment size
    n_haplotype_reads: u16,
    n_reads:     u16,
    n_variants:  u16,
    n_low_qual:  u16,
    n_snv:       u16, // single-nucleotide variant
    n_mnv:       u16, // equal-length multi-nucleotide variant
    n_del:       u16, // simple deletion relative to haplotype
    n_ins:       u16, // simple insertion relative to haplotype
    n_complex:   u16, // a base span replaced with variant span of unequal length
    variants:    String,
}

/// VariantReadsTally aggregates accumulated Variants per VariantReadRecord.
pub struct VariantReadsTally {
    pub tally: FxHashMap<QName, VariantReadInstance>
}
impl VariantReadsTally {

    /// Create a new empty VariantReadsTally object. On tally is instantiated  
    /// per SnvChromWorker.
    pub fn new() -> Self {
        let mut tally = FxHashMap::default();
        tally.reserve(0x10000);
        Self{ tally }
    }

    /// Add a Variant to a VariantReadInstance instance in the tally.
    pub fn add_subclonal_variant(
        &mut self,
        read:        &ReadInstance,
        re_fragment: &ReFragment,
        haplotype:   &Haplotype,
        variant:     &Variant, // allowed subclonal variants only
        min_qual:    PhredQual,
        n_haplotype_reads: usize,
        n_reads:           usize,
    ){
        let instance = self.tally
            .entry(read.qname.clone())
            .or_insert_with(|| VariantReadInstance::new(read, re_fragment, haplotype));
        instance.n_haplotype_reads = n_haplotype_reads as u16;
        instance.n_reads           = n_reads as u16;
        instance.n_variants += 1;
        if min_qual < MIN_SNV_INDEL_QUAL {
            instance.n_low_qual += 1;
        }
        if variant.tgt_bases.is_none() {
            instance.n_ins += 1;
        } else if variant.alt_bases.is_none() {
            instance.n_del += 1;
        } else {   
            let n_tgt_bases = variant.tgt_bases
                .as_ref()
                .map(|s| s.len() as u32)
                .unwrap_or(0);
            let n_alt_bases = variant.alt_bases
                .as_ref()
                .map(|s| s.len() as u32)
                .unwrap_or(0);
            if n_tgt_bases != n_alt_bases {
                instance.n_complex += 1;
            } else if n_tgt_bases == 1 {
                instance.n_snv += 1;
            } else {
                instance.n_mnv += 1;
            }
        } 
        instance.variants.push(variant.to_string(
            re_fragment.start0, 
            min_qual
        ));
    }

    /// Sort and write a set of VariantReadRecords to a temporary file for the
    /// working chromosome.
    pub fn write_sorted(
        tool:   &SnvAnalysisTool,
        worker: &mut SnvChromWorker,
        file_path: String,
    ) -> VariantReadsMetadata {
        let mut csv = OutputCsv::open_csv(
            &file_path, 
            b'\t', 
            false, 
            Some(tool.n_cpu),
        );
        let tally = &worker.variant_reads_tally.tally;
        let mut md = VariantReadsMetadata::new(tally.len());
        for (qname, instance) in tally {
            let record = VariantReadRecord {
                qname:       qname.clone(),
                re_fragment: instance.re_fragment,
                haplotype:   instance.haplotype,
                n_repeat_bases: worker.simple_repeats.get_n_repeat_bases(&instance.re_fragment),
                chrom_index: worker.chrom_index,
                sample_bit:  instance.sample_bit,
                n_bases:     instance.n_bases, // number of bases in read, as compared to re_fragment size
                n_haplotype_reads: instance.n_haplotype_reads,
                n_reads:     instance.n_reads,
                n_variants:  instance.n_variants,
                n_low_qual:  instance.n_low_qual,
                n_snv:       instance.n_snv, // single-nucleotide variant
                n_mnv:       instance.n_mnv, // equal-length multi-nucleotide variant
                n_del:       instance.n_del, // simple deletion relative to haplotype
                n_ins:       instance.n_ins, // simple insertion relative to haplotype
                n_complex:   instance.n_complex, // a base span replaced with variant span of unequal length
                variants:    instance.variants.join(",")
            };
            csv.serialize(&record);

            match record.n_snv {
                0 => md.n_indel_only += 1,
                1 => md.n_one_snv    += 1,
                2 => md.n_two_snv    += 1,
                3 => md.n_three_snv  += 1,
                4 => md.n_four_snv   += 1,
                _ => md.n_five_snv   += 1,
            }
        }
        md
    }
}
