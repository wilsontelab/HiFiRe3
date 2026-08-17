# utility functions for loading RE site data and related
# caller must call stopSpinner()
snvData_create <- "asNeeded"

# # get track data
# hf3_getPileup <- function(sourceId, coord, readType){
#     if(readType == "error_corrected"){
#         hf3_getTrackData_bgz(sourceId, "errorCorrectedPileupBgz", coord, use_chrom = FALSE, debug = FALSE)
#     } else {
#         hf3_getTrackData_bgz(sourceId, "allReadsPileupBgz", coord, use_chrom = FALSE, debug = FALSE)
#     }
# }
# hf3_getVariants <- function(sourceId, coord, readType){
#     if(readType == "error_corrected"){
#         hf3_getTrackData_bgz(sourceId, "errorCorrectedVariantsBgz", coord, use_chrom = FALSE, debug = FALSE)
#     } else {
#         hf3_getTrackData_bgz(sourceId, "allReadsVariantsBgz", coord, use_chrom = FALSE, debug = FALSE)
#     }
# }
# hf3_getEncodings <- function(sourceId, coord, readType){
#     if(readType == "error_corrected"){
#         hf3_getTrackData_bgz(sourceId, "errorCorrectedEncodingsBgz", coord, use_chrom = FALSE, debug = FALSE)
#     } else {
#         hf3_getTrackData_bgz(sourceId, "allReadsEncodingsBgz", coord, use_chrom = FALSE, debug = FALSE)
#     }
# }

# get cached data, full tables
hf3_cached_create <- "asNeeded"
hf3_getVariants_cached <- function(sourceId){
    fileType <-"variantsBgz"
    sessionCache$get(
        fileType, 
        key = sourceId, 
        permanent = TRUE,
        from = "ram",
        create = hf3_cached_create,
        createFn = function(...) {
            startSpinner(session, message = "loading variants")
            dataFilePath <- getSourceFilePath(sourceId, fileType)
            d <- fread(
                cmd = paste("zcat", dataFilePath),
                col.names  =  names(hf3_bgzColumns[[fileType]]), 
                colClasses = unname(hf3_bgzColumns[[fileType]])
            )
            d[, ":="(
                n_tgt_bases = nchar(tgt_bases),
                n_alt_bases = nchar(alt_bases)
            )]
            d[clonal == 1, vaf := n_haplotype_reads / n_reads]
            d[, ":="(
                is_snv = n_tgt_bases == 1 & n_tgt_bases == n_alt_bases,
                vaf_bin = floor(vaf * 51) / 51
            )]
            d
        }  
    )$value
}
hf3_getVariantReads_cached <- function(sourceId){
    fileType <-"variantReadsBgz"
    sessionCache$get(
        fileType, 
        key = sourceId, 
        permanent = TRUE,
        from = "ram",
        create = hf3_cached_create,
        createFn = function(...) {
            startSpinner(session, message = "loading variant reads")
            dataFilePath <- getSourceFilePath(sourceId, fileType)
            d <- fread(
                cmd = paste("zcat", dataFilePath),
                col.names  =  names(hf3_bgzColumns[[fileType]]), 
                colClasses = unname(hf3_bgzColumns[[fileType]])
            )
            setkey(d, qname)
            d
        }  
    )$value
}
hf3_getFragments_cached <- function(sourceId, reads_on){
    fileType <- if (reads_on == "reference") {
        "fragmentsOnReferenceBgz"
    } else {
        "fragmentsOnHaplotypeBgz"
    }
    sessionCache$get(
        fileType, 
        key = sourceId, 
        permanent = TRUE,
        from = "ram",
        create = hf3_cached_create,
        createFn = function(...) {
            startSpinner(session, message = paste("loading reads on", reads_on))
            dataFilePath <- getSourceFilePath(sourceId, fileType)
            d <- fread(
                cmd = paste("zcat", dataFilePath),
                col.names  =  names(hf3_bgzColumns[[fileType]]), 
                colClasses = unname(hf3_bgzColumns[[fileType]])
            )
            d[, ":="(
                n_bases = end1 - start0,
                n_bases_bin = floor((end1 - start0) / 250) * 250
            )]
            d
        }  
    )$value
}

hf3_cached_create2 <- "asNeeded"  
hf3_haplotype_cols <- c("chrom_index1", "start0", "end1", "haplotype")
# establish the list of variants to be kept on their own merit, keeping 
# high-quality, subclonal SNVs on:
#   - heterozygous reads (at this point, still including n_haplotype_reads==2)
#   - homozygous reads that are not multivariant
# all matching variants from all samples are present, one row per variant, 
# at this point, sample_bits may carry more than one sample if an SNV
# was detected in more than one sample and might still match a clonal SNV 
# position (these are filtered below)
hf3_getValidSubclonalSnvs <- function(sourceId){
    sessionCache$get(
        "validSubclonalSnvs", 
        key = sourceId, 
        permanent = TRUE,
        from = "ram",
        create = hf3_cached_create2,
        # create = "once",
        createFn = function(...) {
            startSpinner(session, message = "loading valid subclonal SNVs")

            # collect high-quality subclonal SNVs
            variants <- hf3_getVariants_cached(sourceId)
            subclonal_snvs <- variants[
                is_snv == TRUE & 
                clonal == 0 & # subclonal might have >1 valid read instance in one or more samples
                max_min_qual >= 27 # binned qual levels are 3,10,17,22,27,35,40
            ]
            subclonal_snvs[, n_valid_matching_reads := fcase(
                haplotype == 3, n_matching_reads - n_multivariant_reads, # homozyogous fragments
                default = n_matching_reads # heterozygous fragments, multivariant reads permitted
            )]  

            # message(paste(nrow(subclonal_snvs), " = number of high-quality subclonal SNVs"))
            print(subclonal_snvs[, .N, keyby = .(n_matching_reads)])
            print(subclonal_snvs[, .N, keyby = .(n_valid_matching_reads)])
            # print(subclonal_snvs[, .N, keyby = .(n_samples)])
            # print(subclonal_snvs[, .N, keyby = .(matches_clonal)])

            # reject homozygous multivariant to yield valid SNVs
            subclonal_snvs <- subclonal_snvs[n_valid_matching_reads > 0]

            message(paste(nrow(subclonal_snvs), " = number of valid subclonal SNVs"))
            print(subclonal_snvs[, .N, keyby = .(n_samples)])
            print(subclonal_snvs[, .N, keyby = .(matches_clonal)])

            subclonal_snvs
        }  
    )$value
}
# similarly establish the list of 1-base indels
hf3_getValidSubclonalIndels <- function(sourceId){
    sessionCache$get(
        "validSubclonalIndels", 
        key = sourceId, 
        permanent = TRUE,
        from = "ram",
        create = hf3_cached_create2,
        # create = "once",
        createFn = function(...) {
            startSpinner(session, message = "loading valid subclonal indels")

            # collect high-quality subclonal 1-base indels
            variants <- hf3_getVariants_cached(sourceId)
            subclonal_indels <- variants[
                n_tgt_bases + n_alt_bases == 1 & 
                clonal == 0 & # subclonal might have >1 valid read instance in one or more samples
                max_min_qual >= 27 # binned qual levels are 3,10,17,22,27,35,40
            ]  
            subclonal_indels[, n_valid_matching_reads := fcase(
                haplotype == 3, n_matching_reads - n_multivariant_reads, # homozyogous fragments
                default = n_matching_reads # heterozygous fragments, multivariant reads permitted
            )]  

            # message(paste(nrow(subclonal_indels), " = number of high-quality subclonal 1-base indels"))
            print(subclonal_indels[, .N, keyby = .(n_matching_reads)])
            print(subclonal_indels[, .N, keyby = .(n_valid_matching_reads)])
            # print(subclonal_indels[, .N, keyby = .(n_samples)])
            # print(subclonal_indels[, .N, keyby = .(matches_clonal)])

            # reject homozygous multivariant to yield valid SNVs
            subclonal_indels <- subclonal_indels[n_valid_matching_reads > 0]

            message(paste(nrow(subclonal_indels), " = number of valid subclonal 1-base indels"))
            print(subclonal_indels[, .N, keyby = .(n_samples)])
            print(subclonal_indels[, .N, keyby = .(matches_clonal)])

            subclonal_indels
        }  
    )$value
}
# establish the list of fragment-haplotypes to be kept on their own merit, rejecting:
#   - (heterozygous) fragment-haplotypes with only two reads
#   - entire fragments with too little base complexity, i.e., too many masked simple repeats
# this is a fragment/haplotype-level analysis and result; all sample data are 
# aggregated when making decisions about coverage sufficiency for determining 
# accurate haplotype consensuses, etc.
hf3_getValidHaplotypes <- function(sourceId){
    sessionCache$get(
        "validHaplotypes", 
        key = sourceId, 
        permanent = TRUE,
        from = "ram",
        create = hf3_cached_create2,
        createFn = function(...) {
            startSpinner(session, message = "loading valid haplotypes")

            # collect and characterize reads_on_hap fragment-haplotypes
            haplotypes <- hf3_getFragments_cached(sourceId, "haplotype")
            haplotypes[, n_bases := end1 - start0]
            haplotypes[, n_unmasked_bases := n_bases - n_repeat_bases]
            haplotypes[, ":="(
                n_valid_reads = fcase(
                    haplotype == 3, n_reads - n_multivariant_reads, # homozyogous fragments
                    default = n_reads # heterozygous fragments, multivariant reads permitted
                ),
                frac_unmasked = n_unmasked_bases / n_bases
            )]

            message(paste(nrow(haplotypes), " = number of fragment-haplotypes, unfiltered"))
            print(haplotypes[, .(
                .N, 
                n_reads          = sum(n_reads),
                n_valid_reads    = sum(n_valid_reads),
                n_bases          = sprintf("%.2e", sum(n_bases)), 
                n_unmasked_bases = sprintf("%.2e", sum(n_unmasked_bases))
            ), keyby = .(haplotype)])
            print(haplotypes[, .(
                n_reads          = sum(n_reads),
                n_valid_reads    = sum(n_valid_reads),
                n_bases          = sprintf("%.2e", sum(n_bases)), 
                n_unmasked_bases = sprintf("%.2e", sum(n_unmasked_bases))
            ), keyby = .(round(frac_unmasked * 20, 0) / 20)])

            # reject low complexity and low coverage heterozygous to yield valid haplotypes
            haplotypes <- haplotypes[
                n_valid_reads >= 3 &  # need three reads to establish a reliable haplotype consensus
                frac_unmasked >= 0.75 # exclude low complexity fragments; ~99.5% of fragments pass this filter
            ]

            message(paste(nrow(haplotypes), " = number of valid high-complexity fragment-haplotypes"))
            print(haplotypes[, .(
                .N, 
                n_reads          = sum(n_reads),
                n_valid_reads    = sum(n_valid_reads),
                n_bases          = sprintf("%.2e", sum(n_bases)), 
                n_unmasked_bases = sprintf("%.2e", sum(n_unmasked_bases))
            ), keyby = .(haplotype)])
            print(
                haplotypes[n_valid_reads <= 65, .N, keyby = .(haplotype, n_valid_reads)] %>% 
                dcast(n_valid_reads ~ haplotype, value.var = "N")
            )

            haplotypes
        }  
    )$value
}
# use the list of valid SNVs to further reject entire fragment-haplotypes, where
# untrustworthy fragment-haplotypes are rejected if they have 
#   - excessive subclonal SNVs (often across multiple samples and positions)
#   - excessive multivariant reads
# this is a haplotype-level analysis and result; all sample SNVs are aggregated 
# when using SNV counts to reject suspicious, misbehaving haplotypes and entire 
# haplotypes are kept or rejected
hf3_getTrustworthyHaplotypes <- function(sourceId, subclonal_snvs, haplotypes){
    sessionCache$get(
        "trustworthyHaplotypes", 
        key = sourceId, 
        permanent = TRUE,
        from = "ram",
        create = hf3_cached_create2,
        createFn = function(...) {
            startSpinner(session, message = "loading trustworthy haplotypes")

            # left-join valid haplotypes to valid SNVs
            haplotype_snvs <- merge(
                haplotypes[,     .SD, .SDcols = c(hf3_haplotype_cols, "n_multivariant_reads")], 
                subclonal_snvs[, .SD, .SDcols = c(hf3_haplotype_cols, "n_valid_matching_reads")], 
                by = hf3_haplotype_cols, 
                all.x = TRUE, 
                all.y = FALSE, 
                sort = TRUE
            )
            haplotype_filter <- haplotype_snvs[,
                .( 
                    n_snv_bases = sum(n_valid_matching_reads, na.rm = TRUE) 
                ),
                keyby = c(hf3_haplotype_cols, "n_multivariant_reads")
            ]

            print(
                haplotype_filter[, .N, keyby = .(haplotype, n_snv_bases)] %>% 
                dcast(n_snv_bases ~ haplotype, value.var = "N")
            )
            print(
                haplotype_filter[, .N, keyby = .(haplotype, n_multivariant_reads)] %>% 
                dcast(n_multivariant_reads ~ haplotype, value.var = "N")
            )

            # inner-join haplotypes and the filter based on total valid SNV count
            haplotypes$n_multivariant_reads <- NULL
            haplotypes <- merge(
                haplotype_filter[
                    n_snv_bases <= 5 & 
                    n_multivariant_reads <= 5
                ],
                haplotypes,
                by = hf3_haplotype_cols, 
                all.x = FALSE, 
                all.y = FALSE,
                sort = TRUE
            )

            message(paste(nrow(haplotypes), " = number of trustworthy fragment-haplotypes, <= 5 SNVs"))
            
            haplotypes
        }  
    )$value
}
# use trustworthy haplotypes to remove untrustworthy subclonal SNVs from the list
# also enforce:
#   - n_samples==1 as true subclonal SNVs will only be in one sample
#   - matches_clonal==0 (false), to further catch read-haplotype mismatches
hf3_getTrustworthySubclonalSnvs <- function(sourceId, subclonal_snvs, haplotypes){
    sessionCache$get(
        "trustworthySubclonalSnvs", 
        key = sourceId, 
        permanent = TRUE,
        from = "ram",
        create = hf3_cached_create2,
        # create = "once",
        createFn = function(...) {
            startSpinner(session, message = "loading trustworthy SNVs")

            # # report some statistics
            # tmp <- merge(
            #     haplotypes[, .SD, .SDcols = hf3_haplotype_cols],
            #     subclonal_snvs[
            #         matches_clonal == 0
            #     ], 
            #     by = hf3_haplotype_cols, 
            #     all.x = FALSE, 
            #     all.y = FALSE, 
            #     sort = TRUE
            # ) 
            # print(tmp[, .N, keyby = .(n_valid_matching_reads)])
            # print(tmp[, .N, keyby = .(n_samples)])

            # inner-join trustworthy haplotypes to valid SNVs to yield trustworthy SNVs
            subclonal_snvs <- merge(
                haplotypes[, .SD, .SDcols = hf3_haplotype_cols],
                subclonal_snvs[
                    n_samples == 1 & 
                    matches_clonal == 0
                ], 
                by = hf3_haplotype_cols, 
                all.x = FALSE, 
                all.y = FALSE, 
                sort = TRUE
            )
            subclonal_snvs[, sample := hf3_getSampleNames(sourceId, sample_bits, as_string = FALSE)]

            message(paste(nrow(subclonal_snvs), " = number of trustworthy single-sample subclonal SNVs"))
            print(subclonal_snvs[, .N, keyby = .(n_valid_matching_reads)])
            print(subclonal_snvs[n_matching_reads == 1, .(n_read_snvs = .N), by = .(qnames)][, .N, keyby = .(n_read_snvs)])

            # parse ref and alt bases into mutation types
            comp <- c("A" = "T", "C" = "G", "G" = "C", "T" = "A")
            subclonal_snvs[, mutation := fcase(
                tgt_bases %in% c("C", "T"), paste0(     tgt_bases,  ">",      alt_bases),
                tgt_bases %in% c("A", "G"), paste0(comp[tgt_bases], ">", comp[alt_bases])
            )]
            subclonal_snvs[, context := fcase(
                is.na(context_base_left),       NA_character_,
                grepl("N", context_base_left),  NA_character_,
                grepl("N", context_base_right), NA_character_,
                tgt_bases %in% c("C", "T"), paste0(
                    context_base_left,      
                    "[", mutation, "]",
                    context_base_right
                ),
                tgt_bases %in% c("A", "G"), paste0(
                    comp[context_base_right],
                    "[", mutation, "]",
                    comp[context_base_left]
                )
            )]

            # print(
            #     subclonal_snvs[, .N, keyby = .(tgt_bases, alt_bases)] %>% 
            #     dcast(tgt_bases ~ alt_bases, value.var = "N")
            # )
            print(
                subclonal_snvs[, .N, keyby = .(mutation, sample_bits)] %>% 
                dcast(mutation ~ sample_bits, value.var = "N")
            )
            print(
                subclonal_snvs[!is.na(context), .N, keyby = .(mutation, context, sample_bits)] %>% 
                dcast(mutation + context ~ sample_bits, value.var = "N")
            )

            subclonal_snvs
        }  
    )$value
}
# similarly get trustrworthy indels
hf3_getTrustworthySubclonalIndels <- function(sourceId, subclonal_indels, haplotypes){
    sessionCache$get(
        "trustworthySubclonalIndels", 
        key = sourceId, 
        permanent = TRUE,
        from = "ram",
        create = hf3_cached_create2,
        # create = "once",
        createFn = function(...) {
            startSpinner(session, message = "loading trustworthy indels")

            # inner-join trustworthy haplotypes to valid indels to yield trustworthy indels
            subclonal_indels <- merge(
                haplotypes[, .SD, .SDcols = hf3_haplotype_cols],
                subclonal_indels[
                    n_samples == 1 & 
                    matches_clonal == 0
                ], 
                by = hf3_haplotype_cols, 
                all.x = FALSE, 
                all.y = FALSE, 
                sort = TRUE
            )
            subclonal_indels[, sample := hf3_getSampleNames(sourceId, sample_bits, as_string = FALSE)]

            message(paste(nrow(subclonal_indels), " = number of trustworthy single-sample subclonal 1-base indels"))

            # parse ref and alt bases into mutation types
            comp <- c("A" = "T", "C" = "G", "G" = "C", "T" = "A")
            subclonal_indels[, mutation := fcase(
                n_tgt_bases == 1, paste0("-", fcase( # 1-base deletion
                    tgt_bases %in% c("C", "T"),      tgt_bases,
                    tgt_bases %in% c("A", "G"), comp[tgt_bases]
                )),
                default = paste0("+", fcase( # 1-base insertion
                    alt_bases %in% c("C", "T"),      alt_bases,
                    alt_bases %in% c("A", "G"), comp[alt_bases]
                ))
            )]
            subclonal_indels[, context := fcase(
                is.na(context_base_left),       NA_character_,
                grepl("N", context_base_left),  NA_character_,
                grepl("N", context_base_right), NA_character_,
                default = fcase(
                    n_tgt_bases == 1, fcase(
                        tgt_bases %in% c("C", "T"), paste0(
                            context_base_left,      
                            "[", mutation, "]",
                            context_base_right
                        ),
                        tgt_bases %in% c("A", "G"), paste0(
                            comp[context_base_right],
                            "[", mutation, "]",
                            comp[context_base_left]
                        )                        
                    ),
                    default = fcase(
                        alt_bases %in% c("C", "T"), paste0(
                            context_base_left,      
                            "[", mutation, "]",
                            context_base_right
                        ),
                        alt_bases %in% c("A", "G"), paste0(
                            comp[context_base_right],
                            "[", mutation, "]",
                            comp[context_base_left]
                        )                        
                    )
                )
            )]

            print(
                subclonal_indels[!is.na(context), .N, keyby = .(mutation, context, sample_bits)] %>% 
                dcast(mutation + context ~ sample_bits, value.var = "N")
            )

            subclonal_indels
        }  
    )$value
}
# expand trustworthy haplotypes to one read per row with metadata for further  
# filtering and grouping, to begin the sample-level analysis of haplotypes
# remove multivariant reads from homozygous haplotypes - it is possible that 
# this will leave zero fragments in a haplotype for a specific sample, which is OK
# note that all reads are present, even the majority invariant reads
hf3_getHaplotypeReads <- function(sourceId, haplotypes){
    sessionCache$get(
        "haplotypeReads", 
        key = sourceId, 
        permanent = TRUE,
        from = "ram",
        create = hf3_cached_create2,
        createFn = function(...) {
            startSpinner(session, message = "expanding haplotype reads")

            # expand haplotypes data to one row per source read
            haplotype_reads <- haplotypes[, 
                .(
                    sample_bit      = as.integer(strsplit(sample_bitss, ",")[[1]]),
                    # qname           =            strsplit(qnames, ",")[[1]],
                    n_read_variants = as.integer(strsplit(n_variantss, ",")[[1]])
                ),
                keyby = c(hf3_haplotype_cols, "n_unmasked_bases")
            ]

            message(paste(nrow(haplotype_reads), " = number of reads in trustworthy haplotypes"))
            print(
                haplotype_reads[n_read_variants <= 10, .N, keyby = .(n_read_variants, haplotype)] %>%
                dcast(n_read_variants ~ haplotype, value.var = "N")
            )

            #  remove multivariant homozogyous reads to match the valid SNVs filter above
            haplotype_reads <- haplotype_reads[
                haplotype != 3 | # all heterozygous haplotype reads are informative
                n_read_variants <= 1 # cannot trust multivariant homozogyous reads
            ]

            message(paste(nrow(haplotype_reads), " = number of valid reads in trustworthy haplotypes"))
            print(
                haplotype_reads[n_read_variants <= 10, .N, keyby = .(n_read_variants, haplotype)] %>%
                dcast(n_read_variants ~ haplotype, value.var = "N")
            )

            haplotype_reads
        }  
    )$value
}

# varTypeIs <- list(
#     snp     = 1, # single-nucleotide polymorphism, e.g. A>G
#     del1    = 2, # deletion of 1 base, e.g. A>-
#     ins1    = 3, # insertion of 1 base, e.g. ->A
#     indel0  = 4, # equal insertion + deletion of >1 base; a multi-nucleotide polymorphism, e.g. AT>CG
#     delN    = 5, # deletion of >1 base, e.g. AT>-
#     insN    = 6, # insertion of >1 base, e.g. ->AT
#     indelX  = 7, # unequal insertion + deletion of >1 base, e.g. AT>CGT, catch-all for complex events
#     # -----------
#     match   = 8, # no variant, e.g. A>A
#     clipped = 9, # 5'-clipped base
#     lowQual = 10 # masked low quality base
# )
# varTypeColors <- c(
#     snp = "blue",
#     del1 = "red",
#     ins1 = "green",
#     indel0 = "purple",
#     delN = "orange",
#     insN = "cyan",
#     indelX = "brown",
#     match = "gray",
#     clipped = "white",
#     lowQual = "white"
# )
# pileupCodes <- list(
#     CS_MATCH        = ":",
#     CLIPPED         = "!",
#     MASKED_LOW_QUAL = "q"
# )
# pileupCodeVarTypeIs <- c(
#     ":" = varTypeIs$match,
#     "!" = varTypeIs$clipped,
#     "q" = varTypeIs$lowQual
# )
# pileupCodeVals <- names(pileupCodeVarTypeIs)
# varTypes <- names(varTypeIs)
# parseVarType <- function(refBases, altBases){
#     nRefBases <- nchar(refBases)
#     nAltBases <- nchar(altBases)
#     delta <- nAltBases - nRefBases
#     ifelse(
#         refBases %in% pileupCodeVals,
#         pileupCodeVarTypeIs[refBases],
#         ifelse(
#             delta == 0,
#             ifelse(nRefBases == 1, varTypeIs$snp, varTypeIs$indel0),
#             ifelse(
#                 delta == 1,
#                 ifelse(nRefBases == 0, varTypeIs$ins1, varTypeIs$indelX),
#                 ifelse(
#                     delta == -1,
#                     ifelse(nRefBases == 1, varTypeIs$del1, varTypeIs$indelX),
#                     ifelse(
#                         delta > 1,
#                         ifelse(nRefBases == 0, varTypeIs$insN, varTypeIs$indelX),
#                         ifelse(nAltBases == 0, varTypeIs$delN, varTypeIs$indelX)
#                     )
#                 )
#             )
#         )
#     )
# }
# parseVarType_long <- function(varType, refBases, altBases){
#     nRefBases <- nchar(refBases)
#     nAltBases <- nchar(altBases)
#     delta <- nAltBases - nRefBases
#     ifelse(
#         varType %in% c(varTypeIs$snp, varTypeIs$del1, varTypeIs$ins1),
#         paste(refBases, altBases, sep = ">"),
#         ifelse(
#             varType %in% c(varTypeIs$indel0, varTypeIs$delN),
#             nRefBases,
#             ifelse(
#                 varType == varTypeIs$insN,
#                 nAltBases,
#                 pmax(nRefBases, nAltBases)
#             )
#         )
#     )
# }
