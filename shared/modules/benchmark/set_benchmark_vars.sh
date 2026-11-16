# action:
#     set environment variables for performing benchmarking
# requires:
#     source $MODULES_DIR/library/set_library_vars.sh
# usage:
#     source $MODULES_DIR/analyze/set_analysis_vars.sh

export TRUTHSET_NAME=`basename $TRUTHSET_DIR`

# analysis identifiers
export BENCHMARK_IDENTIFIER=${DATA_NAME}.benchmark
export BENCHMARK_PREFIX=${TASK_DIR}/${BENCHMARK_IDENTIFIER}

# SV benchmark output
export BENCHMARK_BAM=${BENCHMARK_PREFIX}.coordinate.bam
export SNIFFLES_VCF=${BENCHMARK_PREFIX}.sniffles.vcf.gz
export SNIFFLES_FRAGSIZES_VCF=${BENCHMARK_PREFIX}.sniffles.frag_sizes.vcf.gz
export CUTESV_VCF_TMP=${BENCHMARK_PREFIX}.cutesv.vcf
export CUTESV_VCF=${BENCHMARK_PREFIX}.cutesv.vcf.gz
export CUTESV_FRAGSIZES_VCF=${BENCHMARK_PREFIX}.cutesv.frag_sizes.vcf.gz
export BENCHMARK_ADMIXED_CALLS=${BENCHMARK_PREFIX}.${TRUTHSET_NAME}.admixed_calls.txt
export BENCHMARK_AGGREGATED_CALLS=${BENCHMARK_PREFIX}.${TRUTHSET_NAME}.aggregated_calls.txt
