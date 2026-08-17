# generic support for data formats and indexed data retrieval using tabix
# these functions expect that files utilized chromIndex1, not chrom

# types
hf3_strands <- list(
    top    = 0,
    bottom = 1
)
hf3_strands_signed <- list(
    top    =  1,
    bottom = -1
)
hf3_junctions <- list(
    typeToBits = c(
        proper        = 0L,
        deletion      = 1L,
        duplication   = 2L,
        inversion     = 4L,
        translocation = 8L
    ),
    bitsToIndex = c(
        1,
        2,
        NA,
        3,
        NA,
        NA,
        NA,
        4
    ),
    indexToType = c(
        "deletion",
        "duplication",
        "inversion",
        "translocation",
        "intergenome" # translocation + 1
    ),
    indexToTypeLabel = c(
        "del",
        "dup",
        "inv",
        "trans",
        "int-gen" # translocation + 1
    ),
    bitsToType = c(
        "deletion",
        "duplication",
        NA,
        "inversion",
        NA,
        NA,
        NA,
        "translocation"
    ),
    bitsToTypeLabel = c(
        "del    ",
        "dup    ",
        NA,
        "inv    ",
        NA,
        NA,
        NA,
        "trans  "
    ),
    typeToRowOffset = list(
        proper          =   NA,
        deletion        = -1/3,
        duplication     =    0,
        inversion       =  1/3,
        translocation   =    0,
        intergenome     =  1/3,
        mixed           =   NA
    ),
    typeToColor = list(
        proper        = CONSTANTS$plotlyColors$black,
        deletion      = CONSTANTS$plotlyColors$blue,
        duplication   = CONSTANTS$plotlyColors$green,
        inversion     = CONSTANTS$plotlyColors$red,
        translocation = CONSTANTS$plotlyColors$orange,
        intergenome   = CONSTANTS$plotlyColors$purple,
        mixed         = CONSTANTS$plotlyColors$black
    ),
    barColors = function(){
        hf3_junctions$typeToColor[1:5 + 1]
    },
    getTypeFromBits = function(bits, is_intergenomic){
        x <- hf3_junctions$bitsToType[bits]
        if(length(x) == 0 || is.na(x)) return("mixed")
        if(x == "translocation" && is_intergenomic == 1) return("intergenome")
        x
    },
    getTypeLabelFromBits = function(bits, is_intergenomic){
        x <- hf3_junctions$bitsToTypeLabel[bits]
        if(length(x) == 0 || is.na(x)) return(NA)
        if(x == "trans  " && is_intergenomic == 1) return("int-gen")
        x
    },
    getColorsFromBits = function(jxn_types, is_intergenomics){
        sapply(seq_along(jxn_types), function(i){
            type <- hf3_junctions$getTypeFromBits(jxn_types[i], is_intergenomics[i])
            hf3_junctions$typeToColor[[type]]
        })
    },
    getColorsFromBits_no_intergenomic = function(jxn_types){
        sapply(seq_along(jxn_types), function(i){
            type <- hf3_junctions$getTypeFromBits(jxn_types[i], 0)
            hf3_junctions$typeToColor[[type]]
        })
    },
    getOffsetsFromBits = function(jxn_types, is_intergenomics){
        sapply(seq_along(jxn_types), function(i){
            type <- hf3_junctions$getTypeFromBits(jxn_types[i], is_intergenomics[i])
            hf3_junctions$typeToRowOffset[[type]]
        })
    },
    getBitsFromTypes = function(jxnTypes){
        hf3_junctions$typeToBits[jxnType]
    },
    getTypesFromBits = function(jxn_types, is_intergenomics){
        sapply(seq_along(jxn_types), function(i){
            hf3_junctions$getTypeFromBits(jxn_types[i], is_intergenomics[i])
        })
    },
    getTypeLabelsFromBits = function(jxn_types, is_intergenomics){
        sapply(seq_along(jxn_types), function(i){
            hf3_junctions$getTypeLabelFromBits(jxn_types[i], is_intergenomics[i])
        })
    }
)
hf3_junctionClasses <- list(
    indexToClass = c(
        "intergenome",
        "artifact",
        "validated",
        "unexpected2",
        "unexpected1"
    ),
    indexToClassLabel = c(
        "intergenome\nall n",
        "artifact\nn=1",
        "validated\nn>=3",
        "unexpected\nn=2",
        "unexpected\nn=1"
    ),
    sizedLabels = c(
        "artifact\nn=1",
        "validated\nn>=3",
        "unexpected\nn=2",
        "unexpected\nn=1"
    ),
    classToIndex = list(
        intergenome = 1L,
        artifact    = 2L,
        validated   = 3L,
        unexpected2 = 4L,
        unexpected1 = 5L
    ),
    classToColor = list(
        intergenome = CONSTANTS$plotlyColors$purple,
        artifact    = CONSTANTS$plotlyColors$orange,
        validated   = CONSTANTS$plotlyColors$green,
        unexpected2 = CONSTANTS$plotlyColors$red,
        unexpected1 = CONSTANTS$plotlyColors$red
    ),
    getClassColor = function(indices){
        sapply(indices, function(i){
            if(i == 0) CONSTANTS$plotlyColors$black
            else hf3_junctionClasses$classToColor[[hf3_junctionClasses$indexToClass[i]]]
        })
    },
    getJxnClassLabels = function(is){
        sapply(is, function(i){
            hf3_junctionClasses$indexToClassLabel[i]
        })
    }
)
hf3_junctionStrata <- list(
    indexToStratum = c(
        "unexpected1",
        "unexpected2",
        "validated"
    ),
    indexToStratumLabel = c(
        "unexpected\nn=1",
        "unexpected\nn=2",
        "validated\nn>=3"
    ),
    stratumToIndex = list(
        unexpected1 = 1L,
        unexpected2 = 2L,
        validated   = 3L
    ),
    getJxnStratumLabels = function(is){
        sapply(is, function(i){
            hf3_junctionStrata$indexToStratumLabel[i]
        })
    }
)
hf3_alignments <- list(
    typeToColor = list(
        alignment  = CONSTANTS$plotlyColors$black,
        projection = CONSTANTS$plotlyColors$grey
    ) 
)

# column definitions
hf3_bgzColumns <- list(
    filteringSitesBgz = c(
        chrom     = "character",
        sitePos1  = "integer",
        inSilico  = "integer",
        nObserved = "integer"
    ),
    svAlignmentsBgz = c(
        chrom_index1   = "integer",
        ref_pos5       = "integer",
        ref_pos3       = "integer",
        ref_proj3      = "integer",
        strand_index0  = "integer",
        jxn_types      = "integer", # middle floating alignments have two bit-encoded types
        n_jxns         = "integer",
        aln_i          = "integer",
        n_observed     = "integer"
    ),
    svJunctions1Bgz = c(
        chrom_index1_1   = "integer",
        ref_pos1_1       = "integer",
        strand_index0_1  = "integer",
        chrom_index1_2   = "integer",
        ref_pos1_2       = "integer",
        strand_index0_2  = "integer",
        offset           = "integer", 
        jxn_seq          = "character",
        # ------------------------------
        jxn_type         = "integer",
        strands          = "integer",
        sv_size          = "integer",
        # ------------------------------
        q_names          = "character",
        insert_sizes     = "character",
        outer_node1s     = "character",
        outer_node2s     = "character",
        channels         = "character",
        n_jxns           = "character",
        is_duplicates    = "character",
        is_duplexes      = "character",
        # ------------------------------
        aln5_is          = "character",
        qry_pos1_aln5_end3s = "character",
        # ------------------------------
        jxn_orientations = "character",
        jxn_failure_flags= "character",
        aln_failure_flags= "character",
        # ------------------------------
        min_stem_lengths = "character",
        min_mapqs        = "character",
        max_divergences  = "character",
        # ------------------------------
        n_instances      = "integer",
        n_reads          = "integer",
        n_instances_dedup= "integer",
        n_reads_dedup    = "integer",
        # ------------------------------
        has_multi_jxn_read = "integer",
        has_multi_instance_read = "integer",
        has_duplex_read  = "integer",
        is_bidirectional = "integer",
        jxn_failure_flag = "integer",
        aln_failure_flag = "integer",
        any_was_chimeric = "integer",
        any_was_not_chimeric = "integer",
        min_stem_length  = "integer",
        max_min_mapq     = "integer",
        min_max_divergence = "double",
        # ------------------------------
        is_intergenomic  = "integer",
        target_1         = "character",
        target_dist_1    = "integer",
        target_2         = "character",
        target_dist_2    = "integer",
        genes_1          = "character",
        gene_dists_1     = "character",
        genes_2          = "character",
        gene_dists_2     = "character",
        is_excluded_1    = "integer",
        is_excluded_2    = "integer",
        # ------------------------------
        sample_bits      = "integer",
        n_samples        = "integer",
        # ------------------------------
        bkpt_coverage_1 = "integer",
        bkpt_coverage_2 = "integer"
    ),
    svReadPaths = c(
        qname             = "character",
        read_len          = "integer",
        insert_size       = "integer",
        has_passed_jxn    = "integer",
        # ------------------------------
        chroms            = "character",
        pos1s             = "character",
        strand0s          = "character",
        n_ref_bases       = "character",
        qry_start0s       = "character",
        qry_end1s         = "character",
        block_ns          = "character",
        mapqs             = "character",
        divergences       = "character",
        aln_failure_flags = "character",
        jxn_failure_flags = "character",
        cigars            = "character",
        # ------------------------------
        seq_strand0        = "integer",
        seq               = "character",
        qual              = "character"
    ),
    variantsBgz = c( # one entry per allowed clonal or subclonal variant
        chrom_index1     = "integer", # 1
        ref_pos0         = "integer", # 2
        tgt_pos0         = "integer",
        tgt_bases        = "character", 
        alt_bases        = "character", 
        is_indel         = "integer",
        start0           = "integer",
        end1             = "integer",
        haplotype        = "integer", 
        n_repeat_bases   = "integer",
        context_base_left  = "character",
        context_base_right = "character",
        n_matching_reads = "integer",
        n_haplotype_reads = "integer",
        n_reads           = "integer",
        n_multivariant_reads = "integer",
        sample_bits      = "integer",
        n_samples        = "integer",
        clonal           = "integer",
        matches_clonal   = "integer",
        max_min_qual     = "integer",
        qnames           = "character"
    ),
    variantReadsBgz = c( # one entry per read with a subclonal variant
        qname           = "character", # first since indexed by qname
        chrom_index1    = "integer",
        start0          = "integer", # for the RE fragment
        end1            = "integer",
        haplotype       = "integer",
        n_repeat_bases  = "integer",
        sample_bit      = "integer",
        n_bases         = "integer",
        n_haplotype_reads = "integer",
        n_reads         = "integer",
        n_variants      = "integer",
        n_low_qual      = "integer",
        n_snv           = "integer",
        n_mnv           = "integer",
        n_del           = "integer",
        n_ins           = "integer",
        n_complex       = "integer",
        variants        = "character"
    ),
    fragmentsOnReferenceBgz = c( # one entry per fully analyzed ReFragment+Haplotype read set
        chrom_index1   = "integer", # 1 # two files, one as aligned to ref, one as aligned to haplotypes
        start0         = "integer", # 2
        end1           = "integer",
        haplotype      = "integer",
        n_repeat_bases = "integer",
        n_reads        = "integer",
        n_multivariant_reads = "integer",
        n_match        = "integer",
        n_alt          = "integer",
        n_del          = "integer",
        n_ins          = "integer",
        n_masked       = "integer",  
        n_variants     = "integer",  
        sample_bits    = "integer",
        n_matches      = "character",
        n_alts         = "character",
        n_dels         = "character",
        n_inss         = "character",
        n_maskeds      = "character",
        n_variantss    = "character",
        sample_bitss   = "character",
        encodings      = "character",
        insertions     = "character",  
        qnames         = "character",
        seq            = "character",
        hap_vs_ref     = "character"
    )
)
hf3_bgzColumns$svJunctions2Bgz <- hf3_bgzColumns$svJunctions1Bgz
hf3_bgzColumns$fragmentsOnHaplotypeBgz <- hf3_bgzColumns$fragmentsOnReferenceBgz

# column display definitions
hf3_bgzColumns_display <- list(
    svJunctions1Bgz = c(
        chrom_index1_1  = "chrom1",
        ref_pos1_1      = "bkptPos1",
        strand_index0_1 = "strand1",
        chrom_index1_2  = "chrom2",
        ref_pos1_2      = "bkptPos2",
        strand_index0_2 = "strand2",
        sv_size         = "svSize",
        jxn_type        = "jxnType",
        is_intergenomic = "interGen",
        n_samples       = "nSmp",
        n_instances_dedup  = "nObs",
        # is_expected     = "expected",
        max_min_mapq       = "mapQ",
        min_max_divergence = "deTag",
        has_duplex_read    = "hasDpx",
        # site_dist       = "siteDist",
        offset          = "alnOffset",
        target_1        = "target1",
        target_2        = "target2",
        bkpt_coverage_1 = "bkptCov1",
        bkpt_coverage_2 = "bkptCov2",
        # insert_sizes    = "insSizes", 
        # stem_lengths    = "stemLens",
        NULL
    )
)

# sample bit to sample name conversion
hf3_sampleIndex <- list()
hf3_getSampleNames <- Vectorize(function(sourceId, sampleBits_, as_string = TRUE){
    if(is.null(hf3_sampleIndex[[sourceId]])) {
        hf3_sampleIndex[[sourceId]] <<- fread(getSourceFilePath(sourceId, "samplesFile"))
    }
    sample_names <- hf3_sampleIndex[[sourceId]][bitwAnd(sample_bit, sampleBits_) > 0, sample_name]
    if(as_string) sample_names <- paste0(sample_names, collapse = " ")
    sample_names
})
hf3_sample_bits <- function(sourceId){
    if(is.null(hf3_sampleIndex[[sourceId]])) {
        hf3_sampleIndex[[sourceId]] <<- fread(getSourceFilePath(sourceId, "samplesFile"))
    }  
    hf3_sampleIndex[[sourceId]]$sample_bit
}

# chrom to chromIndex conversion for bgz queries
# simple cache, sessionCache not used
hf3_chromIndex <- list()  
hf3_getChromNames <- Vectorize(function(sourceId, chromIndex1_){
    if(is.null(hf3_chromIndex[[sourceId]])) {
        hf3_chromIndex[[sourceId]] <<- fread(getSourceFilePath(sourceId, "chromsFile"))
    }
    hf3_chromIndex[[sourceId]][chrom_index1 == chromIndex1_, chrom]
})
hf3_getChromIndex <- function(sourceId, chrom_){
    if(is.null(hf3_chromIndex[[sourceId]])) {
        hf3_chromIndex[[sourceId]] <<- fread(getSourceFilePath(sourceId, "chromsFile"))
    }
    hf3_chromIndex[[sourceId]][chrom == chrom_, chrom_index1]
}
hf3_getChromSize <- function(sourceId, chromIndex1_){
    if(is.null(hf3_chromIndex[[sourceId]])) {
        hf3_chromIndex[[sourceId]] <<- fread(getSourceFilePath(sourceId, "chromsFile"))
    }
    hf3_chromIndex[[sourceId]][chrom_index1 == chromIndex1_, chrom_size]
}

# insert size and stem length reference distributions
hf3_insertSizes <- list()  
hf3_getInsertSizes <- function(sourceId){
    if(is.null(hf3_insertSizes[[sourceId]])) {
        hf3_insertSizes[[sourceId]] <<- list()
        file <- getSourceFilePath(sourceId, "filteredInsertSizes")
        hf3_insertSizes[[sourceId]]$insertSizes <<- if (length(file) > 0) {
            sizes <- fread(file, col.names = c("type", "bin", "count"))
            nonSvSizes <- sizes[type == "nonSV.projected"][order(bin)]
            if(nrow(nonSvSizes) == 0) nonSvSizes <- sizes[type == "nonSV.actual"][order(bin)]
            nonSvSizes$freq <- nonSvSizes$count / sum(nonSvSizes$count)
            nonSvSizes
        } else NULL
    }
    hf3_insertSizes[[sourceId]]
}

# genome excluded regions
hf3_excludedRegions <- list()
hf3_getExcludedRegions <- function(sourceId){
    if(is.null(hf3_excludedRegions[[sourceId]])) {
        hf3_excludedRegions[[sourceId]] <<- fread(getSourceFilePath(sourceId, "exclusionsBed"), fill = TRUE)[, 1:3]
        setnames(hf3_excludedRegions[[sourceId]], c("chrom", "start0", "end1"))
    }
    hf3_excludedRegions[[sourceId]]
}

# general track data retrieval, uses tabix random access not sessionCache
hf3_getTrackData_bgz <- function(sourceId, fileType, coord, use_chrom = FALSE, debug = FALSE, end_col = 3){
    if(!use_chrom) coord$chromosome <- hf3_getChromIndex(sourceId, coord$chromosome)
    bgzFile <- getSourceFilePath(sourceId, fileType)
    # debug <- TRUE
    if (debug) {
        dmsg()
        dmsg(fileType)
        dmsg(bgzFile)
        dmsg(file.exists(bgzFile))
        dstr(coord)
        # dstr(fread(bgzFile))
    }
    if(!isTruthy(bgzFile) || !file.exists(bgzFile)) return(data.table()) # no data of this type
    getCachedTabix(bgzFile, create = debug, force = debug, end_col = end_col) %>% 
    getTabixRangeData(
        coord, 
        col.names  =  names(hf3_bgzColumns[[fileType]]), 
        colClasses = unname(hf3_bgzColumns[[fileType]]), 
        skipChromCheck = TRUE
    )
}
# get and parse RE site, with padding for plotting
hf3_getSites_padded <- function(sourceId, coord){
    coord$start <- coord$start - coord$width
    coord$end   <- coord$end   + coord$width
    if(coord$start < 1) coord$start <- 1L
    hf3_getTrackData_bgz(sourceId, "filteringSitesBgz", coord, use_chrom = TRUE)
}
# get and parse read alignments
hf3_getAlignments <- function(sourceId, coord){
    hf3_getTrackData_bgz(sourceId, "svAlignmentsBgz", coord, end_col = 2)
}
# get and parse unique junction nodes in region ...
hf3_getJunctions <- function(sourceId, coord){
    # svJunctions1Bgz and svJunctions2Bgz are the same row data, just sorted differently
    # having both files allows for recover of all unique SVs with either junction within coord
    # importantly, svJunctions1Bgz and svJunctions2Bgz do NOT! have reversed nodes!
    # svJunctions1Bgz node1 == svJunctions2Bgz node1 NOT(svJunctions2Bgz node2)
    rbind(
        hf3_getTrackData_bgz(sourceId, "svJunctions1Bgz", coord),
        hf3_getTrackData_bgz(sourceId, "svJunctions2Bgz", coord)
    ) %>% unique()
}
# ... additionally filtered the same way as is active the in exploreJunctions step
hf3_getFilteredJunctions <- function(sourceId, coord){
    hf3_applyJunctionFilters(
        hf3_getJunctions(sourceId, coord), 
        app$exploreJunctions$settingsObject, 
        app$exploreJunctions$input
    )
}

# loading entire junction and tally tables, cached in sessionCache
hf3_jxnClassNames <- c("all_jxns", "validated", "artifact", "intergenome", "candidate")
hf3_jxnClasses <- 1:length(hf3_jxnClassNames)
names(hf3_jxnClasses) <- hf3_jxnClassNames
hf3_jxnCreate <- "asNeeded"
hf3_jxnForce  <- FALSE
hf3_getJunctions_all_source <- function(sourceId){
    fileType <- "svJunctions1Bgz"
    filePath <- loadPersistentFile(
        sourceId = sourceId, # ... or source
        contentFileType = fileType,
        force = hf3_jxnForce,
        header = TRUE,
        colClasses = unname(hf3_bgzColumns[[fileType]]), 
        col.names = names(hf3_bgzColumns[[fileType]]), # additional argument passed to fread
        quote = "", # additional argument passed to fread
        postProcess = function(jxns){
            startSpinner(session, message = "post-process junctions")
            jxns[, ":="(
                jxnI = 1:.N,
                sourceId = sourceId,
                sample_names = hf3_getSampleNames(sourceId, sample_bits, as_string = TRUE),
                is_validated = n_instances_dedup >= 3
            )][, ":="(
                is_artifact = n_instances_dedup == 1 & (
                    jxn_type == hf3_junctions$typeToBits["translocation"] |
                    sv_size > 1e6
                )
            )][, ":="(
                jxnClass = ifelse(
                    is_intergenomic, hf3_junctionClasses$classToIndex$intergenome,
                    ifelse(
                        is_artifact, hf3_junctionClasses$classToIndex$artifact,
                        ifelse(
                            is_validated, hf3_junctionClasses$classToIndex$validated,
                            ifelse(
                                n_instances_dedup == 2, hf3_junctionClasses$classToIndex$unexpected2,
                                hf3_junctionClasses$classToIndex$unexpected1
                            )
                        )
                    )
                ),
                jxnStratum = ifelse(
                    is_validated, hf3_junctionStrata$stratumToIndex$validated,
                    ifelse(
                        n_instances_dedup == 2, hf3_junctionStrata$stratumToIndex$unexpected2,
                        hf3_junctionStrata$stratumToIndex$unexpected1
                    )
                )
            )]
            jxns
        }
    )
    persistentCache[[filePath]]$data
}

# get a single SV-containing read path by qName
getTaskFile <- function(sourceId, suffix){
    outputDir <- getSourcePackageOption(sourceId, "output", "output-dir")
    dataName  <- getSourcePackageOption(sourceId, "output", "data-name")
    taskDir   <- file.path(outputDir, dataName)
    fileName  <- paste0(dataName, suffix)
    file.path(taskDir, fileName) %>% Sys.glob # wildcard expansion
}
hf3_getReadPath <- function(sourceId, qName){
    req(sourceId, qName)
    bgzFile <- getTaskFile(sourceId, '.*.read_paths.txt.bgz')
    tabix   <- getCachedTabix(bgzFile)
    gRange <- GenomicRanges::GRanges(
        seqnames = qName, # qName masquerades as rName for this purpose
        ranges = IRanges::IRanges(start = 0, end = 1e6) # coordinates are irrelevant
    )
    lines <- Rsamtools::scanTabix(tabix, param = gRange)[[1]]
    if(length(lines) == 1) lines <- paste0(lines, "\n")
    fread(
        text = lines,
        col.names  =  names(hf3_bgzColumns$svReadPaths), 
        colClasses = unname(hf3_bgzColumns$svReadPaths), 
        header = FALSE,
        sep = "\t"
    )
}
