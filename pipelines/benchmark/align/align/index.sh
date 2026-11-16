# actions:
#     index benchmark bam file

echo "indexing benchmark bam file"
samtools index ${BENCHMARK_BAM}
checkPipe

echo "done"
