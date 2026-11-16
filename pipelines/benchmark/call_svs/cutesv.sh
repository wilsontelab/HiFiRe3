# run CuteSV on library reads re-aligned using HiFiRe3 minimap2 parameters

export TRUTHSET_ALGORITHM=cutesv
export TRUTHSET_VCF=${CUTESV_VCF}

# suggested ONT parameters
if [ "${SEQUENCING_PLATFORM}" = "ONT" ]; then
    MAX_CLUSTER_BIAS_INS=100
    DIFF_RATIO_MERGING_INS=0.3
    MAX_CLUSTER_BIAS_DEL=100
    DIFF_RATIO_MERGING_DEL=0.3
# fallback to PacBio HiFi parameters
else 
    MAX_CLUSTER_BIAS_INS=1000
    DIFF_RATIO_MERGING_INS=0.9
    MAX_CLUSTER_BIAS_DEL=1000
    DIFF_RATIO_MERGING_DEL=0.5
fi

echo "calling SVs with CuteSV"
rm -f ${CUTESV_VCF_TMP}
rm -f ${CUTESV_VCF}
cuteSV \
    --threads ${N_CPU} \
    --report_readid \
    --min_mapq ${MIN_SV_MAPQ} \
    --min_support ${MIN_N_OBSERVED} \
    --min_size ${MIN_SV_LEN} \
    --max_size -1 \
    --max_cluster_bias_INS ${MAX_CLUSTER_BIAS_INS} \
    --diff_ratio_merging_INS ${DIFF_RATIO_MERGING_INS} \
    --max_cluster_bias_DEL ${MAX_CLUSTER_BIAS_DEL} \
    --diff_ratio_merging_DEL ${DIFF_RATIO_MERGING_DEL} \
    ${BENCHMARK_BAM} \
    ${GENOME_FASTA} \
    ${CUTESV_VCF_TMP} \
    ${TMP_DIR_WRK_SMALL} 2>&1 | grep -v "Finished"
checkPipe
echo "done"

echo "compressing output vcf"
bgzip -@ ${N_CPU} ${CUTESV_VCF_TMP}

echo
echo "determining sizes of RE fragments that contain called SVs"
${HF3_TOOLS_BIN} digest_sv_alleles > ${CUTESV_FRAGSIZES_VCF}
checkPipe
echo "done"

# $ cuteSV --help
# usage: cuteSV <sorted.bam> <reference.fa> <output.vcf> <work_dir>
#         Current version: v2.1.4
#         For PacBio CLR data:
#                 --max_cluster_bias_INS          100
#                 --diff_ratio_merging_INS        0.3
#                 --max_cluster_bias_DEL  200
#                 --diff_ratio_merging_DEL        0.5
#         For PacBio CCS(HIFI) data:
#                 --max_cluster_bias_INS          1000
#                 --diff_ratio_merging_INS        0.9
#                 --max_cluster_bias_DEL  1000
#                 --diff_ratio_merging_DEL        0.5
#         For ONT data:
#                 --max_cluster_bias_INS          100
#                 --diff_ratio_merging_INS        0.3
#                 --max_cluster_bias_DEL  100
#                 --diff_ratio_merging_DEL        0.3

# options:
#   -t, --threads THREADS
#                         Number of threads to use.[16]
#   -b, --batches BATCHES
#                         Batch of genome segmentation interval.[10000000]
#   -S, --sample SAMPLE   Sample name/id
#   --report_readid       Enable to report supporting read ids for each SV.
#   --ignore_sequence     Do not output sequences for SVs.

# Collection of SV signatures:
#   -p, --max_split_parts MAX_SPLIT_PARTS
#                         Maximum number of split segments a read may be aligned before it is ignored. All split segments are considered when using -1. (Recommand -1 when applying assembly-based
#                         alignment.)[7]
#   -q, --min_mapq MIN_MAPQ
#                         Minimum mapping quality value of alignment to be taken into account.[20]
#   -r, --min_read_len MIN_READ_LEN
#                         Ignores reads that only report alignments with not longer than bp.[500]
#   -md, --merge_del_threshold MERGE_DEL_THRESHOLD
#                         Maximum distance of deletion signals to be merged. In our paper, I used -md 500 to process HG002 real human sample data.[0]
#   -mi, --merge_ins_threshold MERGE_INS_THRESHOLD
#                         Maximum distance of insertion signals to be merged. In our paper, I used -mi 500 to process HG002 real human sample data.[100]
#   -include_bed INCLUDE_BED
#                         Optional given bed file. Only detect SVs in regions in the BED file. [NULL]

# Generation of SV clusters:
#   -s, --min_support MIN_SUPPORT
#                         Minimum number of reads that support a SV to be reported.[10]
#   -l, --min_size MIN_SIZE
#                         Minimum size of SV to be reported.[30]
#   -L, --max_size MAX_SIZE
#                         Maximum size of SV to be reported. All SVs are reported when using -1. [100000]
#   -sl, --min_siglength MIN_SIGLENGTH
#                         Minimum length of SV signal to be extracted.[10]

# Computing genotypes:
#   --genotype            Enable to generate genotypes.
#   --gt_round GT_ROUND   Maximum round of iteration for alignments searching if perform genotyping.[500]
#   --read_range READ_RANGE
#                         The interval range for counting reads distribution.[1000]

# Advanced:
#   --max_cluster_bias_INS MAX_CLUSTER_BIAS_INS
#                         Maximum distance to cluster read together for insertion.[100]
#   --diff_ratio_merging_INS DIFF_RATIO_MERGING_INS
#                         Do not merge breakpoints with basepair identity more than [0.3] for insertion.
#   --max_cluster_bias_DEL MAX_CLUSTER_BIAS_DEL
#                         Maximum distance to cluster read together for deletion.[200]
#   --diff_ratio_merging_DEL DIFF_RATIO_MERGING_DEL
#                         Do not merge breakpoints with basepair identity more than [0.5] for deletion.
#   --max_cluster_bias_INV MAX_CLUSTER_BIAS_INV
#                         Maximum distance to cluster read together for inversion.[500]
#   --max_cluster_bias_DUP MAX_CLUSTER_BIAS_DUP
#                         Maximum distance to cluster read together for duplication.[500]
#   --max_cluster_bias_TRA MAX_CLUSTER_BIAS_TRA
#                         Maximum distance to cluster read together for translocation.[50]
#   --diff_ratio_filtering_TRA DIFF_RATIO_FILTERING_TRA
#                         Filter breakpoints with basepair identity less than [0.6] for translocation.
#   --remain_reads_ratio REMAIN_READS_RATIO
#                         The ratio of reads remained in cluster. Set lower when the alignment data have high quality but recommand over 0.5.[1.0]