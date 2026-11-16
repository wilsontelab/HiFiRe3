//! Export a set of FinalJunctions to VCF in the same format as Sniffles2.

// dependencies
use super::FinalJunction;

impl FinalJunction {
    /// Export a set of FinalJunctions to VCF in the same format as Sniffles2.
    pub fn _export_to_vcf(mut _jxns: Vec<Self>){

    }
}


// $ zcat GIAB_PB_SNF_hs1.benchmark.sniffles.vcf.gz | more

// ##fileformat=VCFv4.2
// ##source=Sniffles2_2.8.1
// ##fileDate="2026/11/06 18:20:34"
// ##contig=<ID=chr1,length=248387328>
// ...
// ##contig=<ID=chr21,length=45090682>

// ##ALT=<ID=INS,Description="Insertion">
// ##ALT=<ID=DEL,Description="Deletion">
// ##ALT=<ID=DUP,Description="Duplication">
// ##ALT=<ID=INV,Description="Inversion">
// ##ALT=<ID=BND,Description="Breakend; Translocation">

// ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
// ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
// ##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
// ##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
// ##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase-block, zero if none or not phased">
// ##FORMAT=<ID=ID,Number=1,Type=String,Description="Individual sample SV ID for multi-sample output">

// ##FILTER=<ID=PASS,Description="All filters passed">
// ##FILTER=<ID=GT,Description="Genotype filter">
// ##FILTER=<ID=SUPPORT_MIN,Description="Minimum read support filter">
// ##FILTER=<ID=STDEV_POS,Description="SV Breakpoint standard deviation filter">
// ##FILTER=<ID=STDEV_LEN,Description="SV length standard deviation filter">
// ##FILTER=<ID=COV_MIN,Description="Minimum coverage filter">
// ##FILTER=<ID=COV_MIN_GT,Description="Minimum coverage filter (missing genotype)">
// ##FILTER=<ID=COV_CHANGE_DEL,Description="Coverage change filter for DEL">
// ##FILTER=<ID=COV_CHANGE_DUP,Description="Coverage change filter for DUP">
// ##FILTER=<ID=COV_CHANGE_INS,Description="Coverage change filter for INS">
// ##FILTER=<ID=COV_CHANGE_FRAC_US,Description="Coverage fractional change filter: upstream-start">
// ##FILTER=<ID=COV_CHANGE_FRAC_SC,Description="Coverage fractional change filter: start-center">
// ##FILTER=<ID=COV_CHANGE_FRAC_CE,Description="Coverage fractional change filter: center-end">
// ##FILTER=<ID=COV_CHANGE_FRAC_ED,Description="Coverage fractional change filter: end-downstream">
// ##FILTER=<ID=COV_VAR,Description="Coverage variance exceeded">
// ##FILTER=<ID=MOSAIC_VAF,Description="Mosaic variant allele fraction filter">
// ##FILTER=<ID=NOT_MOSAIC_VAF,Description="Variant allele fraction filter for non-mosaic">
// ##FILTER=<ID=ALN_NM,Description="Length adjusted mismatch filter">
// ##FILTER=<ID=STRAND_BND,Description="Strand support filter for BNDs">
// ##FILTER=<ID=STRAND,Description="Strand support filter for germline SVs">
// ##FILTER=<ID=STRAND_MOSAIC,Description="Strand support filter for mosaic SVs">
// ##FILTER=<ID=SVLEN_MIN,Description="SV length filter">
// ##FILTER=<ID=SVLEN_MIN_MOSAIC,Description="SV length filter for mosaic SVs (min)">
// ##FILTER=<ID=SVLEN_MAX_MOSAIC,Description="SV length filter for mosaic SVs (max)">
// ##FILTER=<ID=SINGLE_BREAK,Description="A single break point was detected but not classified as an SV.">
// ##FILTER=<ID=INLINE_SA,Description="INLINE/CIGAR-based SV is mostly supported by SA reads">
// ##FILTER=<ID=MOSAIC_SV_CLOSE_EDGE,Description="For mosaic SVs, the location is close to the end of the read (either end)">
// ##FILTER=<ID=GT_FAILED,Description="Sniffles was unable to genotype this call.">
// ##FILTER=<ID=UNRESOLVED,Description="Call is within an unresolved region.">

// ##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Structural variation with precise breakpoints">
// ##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Structural variation with imprecise breakpoints">
// ##INFO=<ID=MOSAIC,Number=0,Type=Flag,Description="Structural variation classified as putative mosaic">
// ##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variation">
// ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variation">
// ##INFO=<ID=CHR2,Number=1,Type=String,Description="Mate chromsome for BND SVs">
// ##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of reads supporting the structural variation">
// ##INFO=<ID=SUPPORT_INLINE,Number=1,Type=Integer,Description="Number of reads supporting an INS/DEL SV (non-split events only), or reads for a DUP SV reported as INS by the mapper.">
// ##INFO=<ID=SUPPORT_SA,Number=1,Type=Integer,Description="Number of reads supporting a DEL SV through supplementary alignments (split events)">
// ##INFO=<ID=SUPPORT_LONG,Number=1,Type=Integer,Description="Number of soft-clipped reads putatively supporting the long insertion SV">
// ##INFO=<ID=SUPPORT_UNSCALED,Number=1,Type=Float,Description="Base support for long insertions, before adjustment">
// ##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variation">
// ##INFO=<ID=STDEV_POS,Number=1,Type=Float,Description="Standard deviation of structural variation start position">
// ##INFO=<ID=STDEV_LEN,Number=1,Type=Float,Description="Standard deviation of structural variation length">
// ##INFO=<ID=COVERAGE,Number=.,Type=Float,Description="Coverages near upstream, start, center, end, downstream of structural variation">
// ##INFO=<ID=STRAND,Number=1,Type=String,Description="Strands of supporting reads for structural variant">
// ##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count, summed up over all samples">
// ##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="List of read support for all samples">
// ##INFO=<ID=CONSENSUS_SUPPORT,Number=1,Type=Integer,Description="Number of reads that support the generated insertion (INS) consensus sequence">
// ##INFO=<ID=RNAMES,Number=.,Type=String,Description="Names of supporting reads (if enabled with --output-rnames)">
// ##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant Allele Fraction">
// ##INFO=<ID=COVERAGE_VAR,Number=1,Type=Float,Description="Variance of coverage across large events">
// ##INFO=<ID=NM,Number=.,Type=Float,Description="Mean number of query alignment length adjusted mismatches of supporting reads">
// ##INFO=<ID=PHASE,Number=.,Type=String,Description="Phasing information derived from supporting reads, represented as list of: HAPLOTYPE,PHASESET,HAPLOTYPE_SUPPORT,PHASESET_SUPPO
// RT,HAPLOTYPE_FILTER,PHASESET_FILTER">
// ##INFO=<ID=LASM,Number=0,Type=Flag,Description="Local assembly used to detect the structural variant">

// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE
// chr1    2566    Sniffles2.INS.0S0       c       cCCTCACCCTCACCCTCACCCTCACCCTGACCCTGACCCTGACCCTCACCCTGACCCTCACCCTCACCCTCA        60      PASS    IMPRECISE;SVTYPE=INS;SVLEN=71;END
// =2566;SUPPORT=3;RNAMES=m64018_210505_030429/14092581/ccs,m64018_210503_204922/153747606/ccs,m64049_210424_024840/112986691/ccs;COVERAGE=5,5,5,5,22;STRAND=+-;NM=0.026;PHASE=0,NUL
// L,3,3,PASS,FAIL;STDEV_LEN=5.033;STDEV_POS=54.271;SUPPORT_LONG=0;SUPPORT_UNSCALED=3;VAF=0.600    GT:GQ:DR:DV:PS  0/1:11:2:3:.

// chr1    2707    Sniffles2.DEL.53DS0     TCTgcaggcgcagaggcgcggcgcgccggcgcaggcgcagagaggcgcggcgcgccggcgcaggcgcagagaggcgcggcgcgccggcgcaggcgcagagaggcgcggcgcgccggcgcaggcgcagagaggcgcgg
// cgcgccggcgcaggcgcagagaggcgcggcgcgccggcgcaggcgcagagaggcgcggcgcgccggcgcaggcgcagagaggcgcggcgcgccggcgcaggcgcagagaggcgcggcgcgccg     T       60      PASS    PRECISE;SVTYPE=DEL;SVLEN=
// -259;END=2966;SUPPORT=4;RNAMES=m64049_210425_091619/25821819/ccs,m64018_210503_204922/153747606/ccs,m64049_210424_024840/112986691/ccs,m64018_210505_030429/14092581/ccs;COVERAGE
// =5,5,16,21,22;STRAND=+-;NM=0.027;PHASE=0,NULL,4,4,PASS,FAIL;STDEV_LEN=0.000;STDEV_POS=0.000;SUPPORT_SA=0;VAF=0.286      GT:GQ:DR:DV:PS  0/1:12:10:4:.

// chr1    3227    Sniffles2.INS.2S0       g       gGGGGGCGGGGGGGAGGTGCGGCGCAGGCGCACAGAGGCGCGTCCCCTGGTGGGGCGGGGGGGAGGCGCGGCGCAGGCGCACAGAGGCGCGTCCCCTGGT    60      PASS    IMPRECISE
// ;SVTYPE=INS;SVLEN=99;END=3227;SUPPORT=5;RNAMES=m64049_210424_024840/23463407/ccs,m64049_210425_091619/57608274/ccs,m64018_210505_030429/42600080/ccs,m64049_210425_091619/1409029
// 40/ccs,m64109_200813_162416/101188340/ccs;COVERAGE=15,22,22,22,22;STRAND=+-;NM=0.011;PHASE=0,NULL,5,5,PASS,FAIL;STDEV_LEN=3.055;STDEV_POS=43.466;SUPPORT_LONG=0;SUPPORT_UNSCALED=
// 5;VAF=0.227     GT:GQ:DR:DV:PS  0/1:2:17:5:.

// chr1    3635    Sniffles2.DEL.53ES0     ccctcttgctcacagtgtagtggcagcacgcccgcctgctggcagctggggacactgccgga  c       60      PASS    PRECISE;SVTYPE=DEL;SVLEN=-61;END=3696;SUPPORT=10;
// RNAMES=m64049_210425_091619/46205079/ccs,m64049_210424_024840/23463407/ccs,m64049_210425_091619/13961553/ccs,m64109_200813_162416/91096314/ccs,m64109_200805_204709/36438866/ccs,
// m64109_200813_162416/101188340/ccs,m64109_200810_062248/42666785/ccs,m64018_210505_030429/93194285/ccs,m64018_210503_204922/180420907/ccs,m64049_210425_091619/49610952/ccs;COVER
// AGE=22,22,22,22,23;STRAND=+-;NM=0.013;PHASE=0,NULL,10,10,PASS,FAIL;STDEV_LEN=0.516;STDEV_POS=0.000;SUPPORT_SA=0;VAF=0.455       GT:GQ:DR:DV:PS  0/1:60:12:10:.

// chr1    4728    Sniffles2.DEL.53FS0     tcacacaactaatcactcaaaatcatccttacactaaaaatgctaaactatacaatttctagaagaaacaatagaagaaaagctatgtgcctttgggtttggtaatgaattttaacaaatgatacaaaaggttgata
// tacacaaaagaaatgacattgtggttttcttaatatttaaagtttatactctggaagacaccttgttaagagaacaaaaagacaagccacatattgaagaaaatatttgcaaaatacacatctgagaaagaatttgtcttcaaaatatataaaaaaatattaaaactaaacaataag
// ttaaacagcccaactaaaaatgca        t       60      PASS    PRECISE;SVTYPE=DEL;SVLEN=-337;END=5065;SUPPORT=7;RNAMES=m64049_210425_091619/57608274/ccs,m64018_210505_030429/42600080/c
// cs,m64049_210425_091619/140902940/ccs,m64049_210425_091619/74187246/ccs,m64049_210422_203419/80282725/ccs,m64018_210503_204922/174130135/ccs,m64049_210426_154417/37750720/ccs;CO
// VERAGE=23,27,29,30,32;STRAND=+-;NM=0.010;PHASE=0,NULL,7,7,PASS,FAIL;STDEV_LEN=0.000;STDEV_POS=0.000;SUPPORT_SA=0;VAF=0.241      GT:GQ:DR:DV:PS  0/1:8:22:7:.

// chr1    11669   Sniffles2.DEL.540S0     TGTCCCTCTGTCTCTGCCAACCAGTTAACCTGCTGCTTCCTGGAGGGAGACAGTCCCTCA    T       58      PASS    PRECISE;SVTYPE=DEL;SVLEN=-59;END=11728;SUPPORT=21
// ;RNAMES=m64109_200815_033514/147786132/ccs,m64049_210426_154417/7341946/ccs,m64018_210503_204922/180420907/ccs,m64049_210424_024840/23463407/ccs,m64049_210425_091619/49610952/cc
// s,m64109_200805_204709/82512543/ccs,m64049_210425_091619/46205079/ccs,m64109_200805_204709/36438866/ccs,m64109_200805_204709/8192129/ccs,m64049_210425_091619/65210563/ccs,m64109
// _200815_033514/172032787/ccs,m64109_200805_204709/95946793/ccs,m64109_200805_204709/66192270/ccs,m64049_210425_091619/13961553/ccs,m64109_200813_162416/101188340/ccs,m64109_2008
// 10_062248/147587175/ccs,m64109_200810_062248/42666785/ccs,m64049_210425_091619/128256494/ccs,m64109_200810_062248/158860825/ccs,m64018_210505_030429/93194285/ccs,m64109_200813_1
// 62416/91096314/ccs;COVERAGE=33,32,32,32,31;STRAND=+-;NM=0.012;PHASE=0,NULL,21,21,PASS,FAIL;STDEV_LEN=0.000;STDEV_POS=0.000;SUPPORT_SA=0;VAF=0.656       GT:GQ:DR:DV:PS  0/1:51:11
// :21:.

// chr1    128171818       Sniffles2.BND.C63S0     t       t]chr22:9623086]        60      PASS    PRECISE;SVTYPE=BND;SUPPORT=27;RNAMES=m64018_210503_204922/121375419/ccs,m64049_21
// 0422_203419/52625523/ccs,m64049_210425_091619/65471772/ccs,m64018_210505_030429/176818097/ccs,m64109_200805_204709/117966216/ccs,m64018_210503_204922/61802184/ccs,m64018_210505_
// 030429/138741616/ccs,m64109_200813_162416/163973170/ccs,m64018_210503_204922/87360845/ccs,m64018_210505_030429/19531319/ccs,m64109_200805_204709/175506474/ccs,m64109_200815_0335
// 14/115804484/ccs,m64049_210425_091619/89588653/ccs,m64109_200815_033514/56295514/ccs,m64018_210503_204922/173474684/ccs,m64018_210503_204922/17434108/ccs,m64049_210426_154417/11
// 4230140/ccs,m64109_200805_204709/35849486/ccs,m64049_210424_024840/47186850/ccs,m64049_210426_154417/166987221/ccs,m64018_210505_030429/46399896/ccs,m64018_210505_030429/2339886
// 9/ccs,m64049_210425_091619/52626862/ccs,m64049_210425_091619/108135756/ccs,m64018_210505_030429/121505875/ccs,m64049_210426_154417/175834531/ccs,m64109_200810_062248/106627249/c
// cs;COVERAGE=65,62,63,8,8;STRAND=+-;NM=254.222;CHR2=chr22;PHASE=0,NULL,27,27,PASS,FAIL;STDEV_POS=0.000;VAF=0.614 GT:GQ:DR:DV:PS  0/1:60:17:27:.
