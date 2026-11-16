# run Sniffles2 on library reads re-aligned using HiFiRe3 minimap2 parameters
# unfortunately, Sniffles2 does not have an option for integer min-support

export TRUTHSET_ALGORITHM=sniffles
export TRUTHSET_VCF=${SNIFFLES_VCF}

echo "calling SVs with Sniffles2"
rm -f ${SNIFFLES_VCF}
sniffles \
    --input ${BENCHMARK_BAM} \
    --vcf ${SNIFFLES_VCF} \
    --reference ${GENOME_FASTA} \
    --threads ${N_CPU} \
    --tmp-dir ${TMP_DIR_WRK_SMALL} \
    --minsvlen ${MIN_SV_LEN} \
    --mapq ${MIN_SV_MAPQ} \
    --max-del-seq-len 500 \
    --output-rnames \
    --allow-overwrite \
    --mosaic \
    --mosaic-include-germline 2>&1 | grep -v "sniffles.worker"
checkPipe
echo "done"

echo
echo "determining sizes of RE fragments that contain called SVs"
${HF3_TOOLS_BIN} digest_sv_alleles > ${SNIFFLES_FRAGSIZES_VCF}
checkPipe
echo "done"


# $ sniffles --help
# usage: sniffles --input SORTED_INPUT.bam [--vcf OUTPUT.vcf] [--snf MERGEABLE_OUTPUT.snf] [--threads 4] [--mosaic]

# Sniffles2: A fast structural variant (SV) caller for long-read sequencing data
#  Version 2.8.1

#  Usage example A - Call SVs for a single sample:
#        sniffles --input sorted_indexed_alignments.bam --vcf output.vcf
#        ... OR, with CRAM input and bgzipped+tabix indexed VCF output:
#          sniffles --input sample.cram --vcf output.vcf.gz
#     Usage example C - Determine genotypes for a set of known SVs (force calling):
#        sniffles --input sample.bam --genotype-vcf input_known_svs.vcf --vcf output_genotypes.vcf

# Common parameters:
#   -i, --input IN [IN ...]                                      For single-sample calling: A coordinate-sorted and indexed .bam/.cram (BAM/CRAM format)
#                                                                file containing aligned reads. - OR - For multi-sample calling: Multiple .snf files
#                                                                (generated before by running Sniffles2 for individual samples with --snf)
#   -v, --vcf OUT.vcf                                            VCF output filename to write the called and refined SVs to. If the given filename ends
#                                                                with .gz, the VCF file will be automatically bgzipped and a .tbi index built for it.
#                                                                (default: None)
#   --snf OUT.snf                                                Sniffles2 file (.snf) output filename to store candidates for later multi-sample
#                                                                calling (default: None)
#   --reference reference.fasta                                  (Optional) Reference sequence the reads were aligned against. To enable output of
#                                                                deletion SV sequences, this parameter must be set. (default: None)
#   --phase                                                      Determine phase for SV calls (requires the input alignments to be phased)
#   -t, --threads N                                              Number of parallel threads to use (speed-up for multi-core CPUs) (default: 4)
#   -c, --contig CONTIG                                          (Optional) Only process the specified contigs. May be given more than once. (default:
#                                                                None)
#   --regions REGIONS.bed                                        (Optional) Only process the specified regions. (default: None)
#   --tmp-dir TMP_DIR                                            (Optional) Directory where temporary files are written, must exist. If it doesn't,
#                                                                default path is used (default: )
#   --all-contigs                                                (Optional) Process all contigs in the input file including small ones (<1MB). (default:
#                                                                False)

# SV Filtering parameters:
#   --minsvlen N                                                 Minimum SV length (in bp). May be prefixed with tilde (e.g. ~50) to allow for slightly
#                                                                smaller SVs if strongly supported. (default: ~50)
#   --mapq N                                                     Alignments with mapping quality lower than this value will be ignored
#   --no-qc, --qc-output-all                                     Output all SV candidates, disregarding quality control steps. (default: False)

# Multi-Sample Calling / Combine parameters:
#   --combine-pctseq COMBINE_PCTSEQ                              Minimum alignment distance as percent of SV length to be merged. Set to 0 to disable
#                                                                alignments for merging. (default: 0.7)
#   --combine-max-inmemory-results COMBINE_MAX_INMEMORY_RESULTS  Maximum number of .snf input files to keep results in memory for. If the number of
#                                                                input files exceeds this value, --no-sort should be given as well to keep the output in
#                                                                a single file. If sorting is requested, the result will be sorted but calls encountered
#                                                                out of order will be written to extra files and NOT be included in the result.
#                                                                (default: 20)
#   --combine-population population.snf                          Name of a population SNF to enable population annotation. (default: None)

# SV Postprocessing, QC and output parameters:
#   --output-rnames                                              Output names of all supporting reads for each SV in the RNAMEs info field (default:
#                                                                False)
#   --no-sort                                                    Do not sort output VCF by genomic coordinates (may slightly improve performance)
#                                                                (default: False)
#   --max-del-seq-len N                                          Maximum deletion sequence length to be output. Deletion SVs longer than this value will
#                                                                be written to the output as symbolic SVs. (default: 50000)
#   --symbolic                                                   Output all SVs as symbolic, including insertions and deletions, instead of reporting
#                                                                nucleotide sequences. (default: False)
#   --allow-overwrite                                            Allow overwriting output files if already existing (default: False)

# Mosaic calling mode parameters:
#   --mosaic                                                     Set Sniffles run mode to detect rare, somatic and mosaic SVs (default: False)
#   --mosaic-af-min F                                            Minimum allele frequency for mosaic SVs to be output (default: 0.05)
#   --mosaic-include-germline                                    Report germline SVs as well in mosaic mode (default: False)

# Developer parameters:
#   --tandem-repeats IN.bed                                      (Optional) Input .bed file containing tandem repeat annotations for the reference
#                                                                genome. (default: None)
