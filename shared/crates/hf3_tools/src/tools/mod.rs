
// modules

// general type aliases
pub mod type_aliases;

// ONT read processing
pub mod trim_ont;     // the main trimming tool for HiFiRe3 ONT basecalling
pub mod reformat_ont; // alternative to reformat legacy basecall and trim files

// PacBio by-strand read merging
pub mod basecall_pacbio;

// read alignment analysis
pub mod post_process_ultima; // replace the pre-aignment actions on non-Ultima reads
pub mod analyze_alignments;  // applied to the minimap2 alignment stream
pub mod analyze_inserts;     // applied to all alignments after site localization

// structural variant analysis
pub mod split_bam_by_chrom_sv; // create temporary chrom-level BAM files
pub mod analyze_svs;           // analyze chrom-level BAM files then aggregate SV data
pub mod merge_svs;             // combine final junctions files across different library types

// structural variant benchmarking
pub mod digest_sv_alleles;     // digest allelic consensus sequences to determine RE fragment sizes
pub mod compare_sv_calls;      // run SV benchmark comparison

// single-nucleotide variant/indel analysis
pub mod split_bam_by_chrom_snv; // create temporary chrom-level BAM files for SNV/indel analysis
pub mod analyze_snvs;           // analyze RE fragment haplotypes for clonal and subclonal variants