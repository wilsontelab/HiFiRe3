#----------------------------------------------------------------------
# server components for the exploreFragments appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
exploreFragmentsServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'exploreFragments'
appStepDir <- getAppStepDir(module)
options <- setDefaultOptions(options, stepModuleInfo[[module]])
settings <- activateMdiHeaderLinks( # uncomment as needed
    session,
    # url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
    # dir = appStepDir, # for terminal emulator
    envir = environment(), # for R console
    baseDirs = appStepDir, # for code viewer/editor
    settings = file.path(app$sources$suiteGlobalDir, "settings", "jxn_filters.yml"), #id, # for step-level settings
    size = "m"
)

#----------------------------------------------------------------------
# adjustable color palette
#----------------------------------------------------------------------
fullBaseColors <- list( # generally follow IGV base color conventions
    M = rgb(0.75, 0.75, 0.75),   # any base match = light grey
    "=" = rgb(0.75, 0.75, 0.75),   # any base match = light grey

    A = rgb(0,    0.8,    0),    # green  
    C = rgb(0,    0,    1),    # blue
    G = rgb(0.82, 0.43, 0),      # orange
    T = rgb(1,    0,    0),    # red
    N = rgb(0.9, 0.9, 0.9),      # N, treated as M
    
    a = rgb(0.5,    0.8,    0.5),    # green  
    c = rgb(0.6,    0,    0.6),    # blue
    g = rgb(0.82, 0.63, 0.4),      # orange
    t = rgb(1,    0.6,    0.6),    # red
    n = rgb(0.9, 0.9, 0.9),      # N, treated as M

    D = rgb(0.1, 0.1, 0.1),      # deleted/missing = black
    I = rgb(0.75,   0,    0.75),  # insertion = purple
    "-" = rgb(0.1, 0.1, 0.1),      # deleted/missing = black
    "+" = rgb(0.75,   0,    0.75),  # insertion = purple
    d = rgb(0.9, 0.9, 0.9), # masked indel encodings
    i = rgb(0.9, 0.9, 0.9)
)
#---------------------------------------
noMaskColors <- fullBaseColors
noMaskColors$N <- noMaskColors$M
noMaskColors$n <- noMaskColors$M
noMaskColors$d <- noMaskColors$M
noMaskColors$i <- noMaskColors$M
#---------------------------------------
noIndelColors <- fullBaseColors
noIndelColors$D <- noIndelColors$M
noIndelColors$I <- noIndelColors$M
noIndelColors[["-"]] <- noIndelColors$M
noIndelColors[["+"]] <- noIndelColors$M
noIndelColors$d <- noIndelColors$M
noIndelColors$i <- noIndelColors$M
#---------------------------------------
noMaskOrIndelColors <- noIndelColors
noMaskOrIndelColors$N <- noMaskOrIndelColors$M
noMaskOrIndelColors$n <- noMaskOrIndelColors$M
#---------------------------------------
baseColors <- reactive({
    switch(
        input$colorPalette,
        "show_all"    = fullBaseColors,
        "hide_masked" = noMaskColors,
        "hide_indels" = noIndelColors,
        "hide_masked_and_indels" = noMaskOrIndelColors,
    )
})

#----------------------------------------------------------------------
# load data
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer("dataSource", selection = "single") 
variants <- reactive({
    sourceId <- req(sourceId())
    startSpinner(session, message = "loading variants")
    hf3_getVariants_cached(sourceId)
})
variant_reads <- reactive({
    sourceId <- req(sourceId())
    startSpinner(session, message = "loading variants")
    hf3_getVariantReads_cached(sourceId)
})
reads_on_reference <- reactive({
    sourceId <- req(sourceId())
    startSpinner(session, message = "loading reads on reference")
    hf3_getFragments_cached(sourceId, "reference")
})
reads_on_haplotype <- reactive({
    sourceId <- req(sourceId())
    startSpinner(session, message = "loading reads on haplotype")
    hf3_getFragments_cached(sourceId, "haplotype")
})
called_haplotypes <- reactive({
    sourceId <- req(sourceId())

    startSpinner(session, message = "loading valid subclonal SNVs")
    valid_subclonal_snvs <- hf3_getValidSubclonalSnvs(sourceId)

    startSpinner(session, message = "loading valid haplotypes")
    valid_haplotypes <- hf3_getValidHaplotypes(sourceId)

    startSpinner(session, message = "loading trustworthy haplotypes")
    trustworthy_haplotypes <- hf3_getTrustworthyHaplotypes(
        sourceId, valid_subclonal_snvs, valid_haplotypes
    )
    trustworthy_haplotypes <- trustworthy_haplotypes[n_snv_bases > 0]

    startSpinner(session, message = "loading trustworthy subclonal SNVs")
    trustworthy_subclonal_snvs <- hf3_getTrustworthySubclonalSnvs(
        sourceId, valid_subclonal_snvs, trustworthy_haplotypes
    )
    trustworthy_subclonal_snvs[, has_odd_one := context %in% c(
        "A[C>A]C",
        "G[C>G]C",
        "A[T>A]T",
        "T[T>A]A",
        "G[T>G]T",
        "T[T>G]G"
    )]

    called_haplotypes <- merge(
        trustworthy_haplotypes,
        trustworthy_subclonal_snvs[, 
            .(
                n_called_snvs = .N,
                n_odd_calls = sum(has_odd_one)
            ), 
            keyby = hf3_haplotype_cols
        ],
        all.x = TRUE,
        all.y = FALSE,  
    )
    called_haplotypes <- called_haplotypes[, 
        .(
            fragment = paste0(
                hf3_getChromNames(sourceId, chrom_index1), ":",
                start0, "-", end1
            ),
            haplotype,
            n_bases,
            n_unmasked_bases,
            frac_unmasked = round(frac_unmasked, 2),
            n_snv_bases,
            n_called_snvs = fcase(
                is.na(n_called_snvs), 0L,
                default = n_called_snvs
            ),
            n_odd_calls = fcase(
                is.na(n_odd_calls), 0L,
                default = n_odd_calls
            ),
            n_reads,
            n_valid_reads,
            n_multivariant_reads,
            sample_bits
        )
    ]
    stopSpinner(session)
    called_haplotypes
})   

#----------------------------------------------------------------------
# table of haplotypes with called SNVs
#----------------------------------------------------------------------
haplotypesTable <- bufferedTableServer(
    "haplotypesTable",
    id,
    input,
    called_haplotypes,
    selection = 'single',
    selectionFn = function(selectedRows) {
        called_haplotypes <- called_haplotypes()
        updateTextInput(
            session, 
            "jumpToFragment", 
            value = called_haplotypes[selectedRows, fragment]
        )
    },
    options = list()
)

#----------------------------------------------------------------------
# select a fragment-haplotype
#----------------------------------------------------------------------
fragment <- reactiveVal(NULL)
setFragment <- function(reads_on_ref){
    if (nrow(reads_on_ref) != 1) {
        message("!! NO FRAGMENT TO PLOT !!")
        stopSpinner(session)
    }
    req(nrow(reads_on_ref) == 1)
    variants <- variants()
    reads_on_haplotype <- reads_on_haplotype()
    reads_on_hap1 <- reads_on_haplotype[
        chrom_index1 == reads_on_ref$chrom_index1 &
        start0 == reads_on_ref$start0 &
        end1   == reads_on_ref$end1 & 
        bitwAnd(haplotype, 1) == 1 # thus, hap1==1 and homozygous==3
    ]
    reads_on_hap2 <- reads_on_haplotype[
        chrom_index1 == reads_on_ref$chrom_index1 &
        start0 == reads_on_ref$start0 &
        end1   == reads_on_ref$end1 & 
        haplotype == 2
    ]
    variants <- variants[
        chrom_index1 == reads_on_ref$chrom_index1 &
        ref_pos0 >= reads_on_ref$start0 &
        ref_pos0 <= reads_on_ref$end1
    ]
    stopSpinner(session)
    fragment(
        list(
            reads_on_ref  = reads_on_ref,
            reads_on_hap1 = reads_on_hap1,
            reads_on_hap2 = reads_on_hap2,
            variants      = variants
        )
    )
}
observeEvent(input$anyFragment, {
    reads_on_reference <- req(reads_on_reference())
    startSpinner(session, message = "selecting any fragment")
    reads_on_ref <- reads_on_reference[n_reads >= input$minCoverage][sample(.N, 1)]
    setFragment(reads_on_ref)
})
observeEvent(input$singletonSnv, {
    variants <- req(variants())
    reads_on_reference <- req(reads_on_reference())
    startSpinner(session, message = "selecting singleton SNV")
    snvs <- variants[
        n_reads >= input$minCoverage & 
        is_snv == TRUE & 
        n_matching_reads == 1 & 
        clonal == 0 &
        n_samples == 1
    ]
    snv <- snvs[sample(.N, 1)]
    qname <- strsplit(snv$qname, ",")[[1]][1]
    reads_on_ref <- reads_on_reference[grepl(qname, qnames)][1]
    setFragment(reads_on_ref)
})
observeEvent(input$trueSubclonal, {
    variants <- req(variants())
    reads_on_reference <- req(reads_on_reference())
    startSpinner(session, message = "selecting true subclonal SNV")
    snvs <- variants[
        n_reads >= input$minCoverage & 
        is_snv == TRUE & 
        n_matching_reads >= 1 & 
        clonal == 0 &
        matches_clonal == 0 &
        n_samples == 1 & 
        (
            (
                haplotype == 3 &
                n_multivariant_reads == 0
            ) | 
            (
                haplotype != 3 & 
                n_haplotype_reads >= 3
            )
        )
    ]
    snv <- snvs[sample(.N, 1)]
    qname <- strsplit(snv$qname, ",")[[1]][1]
    reads_on_ref <- reads_on_reference[grepl(qname, qnames)][1]
    setFragment(reads_on_ref)
})
observeEvent(input$snv5Read, {
    variant_reads <- req(variant_reads())
    reads_on_reference <- req(reads_on_reference())
    startSpinner(session, message = "selecting read with 5+SNVs")
    read <- variant_reads[n_snv >= 5][sample(.N, 1)]
    reads_on_ref <- reads_on_reference[
        read$chrom_index1 == chrom_index1 & 
        read$start0 == start0 &
        read$end1 == end1
    ]
    setFragment(reads_on_ref)
})
observeEvent(input$jumpToFragment, {
    jumpToFragment <- trimws(input$jumpToFragment)
    sourceId <- req(jumpToFragment)
    sourceId <- req(sourceId())
    reads_on_reference <- req(reads_on_reference())
    startSpinner(session, message = paste("jumping to", jumpToFragment))
    parts <- strsplit(jumpToFragment, ":")[[1]]
    chrom_index1_ <- hf3_getChromIndex(sourceId, parts[1])
    parts <- as.integer(strsplit(parts[2], "-")[[1]])
    reads_on_ref <- reads_on_reference[
        chrom_index1 == chrom_index1_ &
        start0 == parts[1] &
        end1   == parts[2]
    ]
    setFragment(reads_on_ref)
})
output$fragmentSpan = renderText({
    sourceId <- req(sourceId())
    fragment <- req(fragment())
    chrom <- hf3_getChromNames(sourceId, fragment$reads_on_ref$chrom_index1)
    paste0(
        chrom, ":",
        fragment$reads_on_ref$start0, "-",
        fragment$reads_on_ref$end1
    )
})

#----------------------------------------------------------------------
# encoding plot support
#----------------------------------------------------------------------
dpi <- 96
pointsize <- 7
px_per_read <- 5
px_per_base <- 2
initEncodingPlot <- function(plot, d){
    nStackRows <- floor(d$n_bases / input$windowWidthBases) + 1
    nStrackTracks <- d$n_reads + 1
    ymax <- nStackRows * nStrackTracks
    width_pixels  <- input$pixelsPerBase * input$windowWidthBases
    height_pixels <- input$pixelsPerRead * nStrackTracks * nStackRows
    layout <- list(
        width     = width_pixels,
        height    = height_pixels,
        pointsize = pointsize,
        dpi       = dpi
    )
    png(file = plot$pngFile, width = width_pixels, height = height_pixels, units = "px", 
        pointsize = pointsize, res = dpi, type = "cairo")
    par(mar = c(0, 0, 0, 0))
    plot$initializeFrame(
        layout,
        xlim = c(0, input$windowWidthBases),
        ylim = c(0, ymax),
        xaxs = "i",
        yaxs = "i",
        xaxt = "n",
        yaxt = "n"
    )
    rect(
        xleft   = 0, 
        xright  = input$windowWidthBases, 
        ybottom = 0, 
        ytop    = ymax, 
        col     = fullBaseColors$M, 
        border  = NA
    ) 
    list(
        nStackRows = nStackRows, 
        nStrackTracks = nStrackTracks,
        ymax = ymax,
        layout = layout
    )
}
addEncodingBases <- function(tracks, b0, y1, op, height = 1){
    trackRow1 <- floor(b0 / input$windowWidthBases) + 1
    plotTrackRow0 <- tracks$nStackRows - trackRow1
    x0 <- b0 %% input$windowWidthBases
    y0 <- plotTrackRow0 * tracks$nStrackTracks + y1
    baseColors <- baseColors()
    rect(
        xleft   = x0, 
        xright  = x0 + 1, 
        ybottom = y0,
        ytop    = y0 + height, 
        col     = baseColors[[op]], 
        border  = NA
    ) 
}
addTargetBases <- function(tracks, d){
    seq <- strsplit(d$seq, "")[[1]]
    for (b1 in 1:d$n_bases) {
        addEncodingBases(tracks, b1 - 1, 0, seq[b1])
    } 
}
addEncodings <- function(tracks, d, plotR1, dataR1, encodings){
    encoding <- strsplit(encodings[dataR1], "")[[1]]
    b0 <- 0
    i1 <- 1
    while (i1 <= length(encoding)){
        op <- encoding[i1]
        if (op == "="){
            nMatch <- ""
            while (i1 <= length(encoding) && grepl("[0-9]", encoding[i1 + 1])){
                nMatch <- paste0(nMatch, encoding[i1 + 1])
                i1 <- i1 + 1
            }
            b0 <- b0 + as.integer(nMatch)
        } else {
            addEncodingBases(tracks, b0, plotR1, op)
            b0 <- b0 + 1
        }
        i1 <- i1 + 1
    }  
}
addReadVariants <- function(tracks, d, plotOrder){
    encodings  <- strsplit(d$encodings, ",")[[1]]
    insertions <- strsplit(d$insertions, ",")[[1]]
    for (plotR1 in 1:d$n_reads){
        dataR1 <- plotOrder[plotR1]
        addEncodings(tracks, d, plotR1, dataR1, insertions)
        addEncodings(tracks, d, plotR1, dataR1, encodings)
    } 
}
addHaplotypeClonal <- function(tracks, d){
    encoding <- strsplit(d$hap_vs_ref, "")[[1]]
    b0 <- 0
    i1 <- 1
    while (i1 <= length(encoding)){
        op <- encoding[i1]
        if (op == "="){
            nMatch <- ""
            while (i1 <= length(encoding) && grepl("[0-9]", encoding[i1 + 1])){
                nMatch <- paste0(nMatch, encoding[i1 + 1])
                i1 <- i1 + 1
            }
            b0 <- b0 + as.integer(nMatch)
        } else {
            addEncodingBases(tracks, b0, 1, op, d$n_reads)
            b0 <- b0 + 1
        }
        i1 <- i1 + 1
    } 
}

#----------------------------------------------------------------------
# reads_on_ref plot
#----------------------------------------------------------------------
readsOnRefPlot <- mdiInteractivePlotBoxServer(
    "readsOnRefPlot",
    click  = FALSE,
    brush  = FALSE,
    points = FALSE, # set to TRUE to expose relevant plot options
    lines  = FALSE,
    settings = NULL, # an additional settings template as a list()
    defaults = NULL, # list of default settings values use to inialize settings
    create = function(...) {
        fragment <- req(fragment())
        startSpinner(session, message = "readsOnRefPlot")
        d <- fragment$reads_on_ref
        plotOrder <- if (sum(fragment$variants$clonal == 1) >= 1){
            heterozygous <- fragment$variants[clonal == 1][which.min(abs(vaf - 0.5))]
            het_qnames <- strsplit(heterozygous$qnames, ",")[[1]]
            qnames <- strsplit(d$qnames, ",")[[1]]
            is_het <- qnames %in% het_qnames
            c(which(is_het), which(!is_het))
        } else 1:d$n_reads
        tracks <- initEncodingPlot(readsOnRefPlot, d)
        addTargetBases(tracks, d)
        addReadVariants(tracks, d, plotOrder)
        abline(h = 0:tracks$ymax, lwd = 0.5, col = rgb(0.5,0.5,0.5))
        stopSpinner(session)
        readsOnRefPlot$finishPng(tracks$layout)
    }
)

#----------------------------------------------------------------------
# reads_on_hap plots
#----------------------------------------------------------------------
readsOnHapPlot <- function(id, field){
    plot <- mdiInteractivePlotBoxServer(
        id,
        click  = FALSE,
        brush  = FALSE,
        points = FALSE, # set to TRUE to expose relevant plot options
        lines  = FALSE,
        settings = NULL, # an additional settings template as a list()
        defaults = NULL, # list of default settings values use to inialize settings
        create = function(...) {
            fragment <- req(fragment())
            variant_reads <- req(variant_reads())
            startSpinner(session, message = id)
            tryCatch({
                d <- fragment[[field]]
                re_start0 <- fragment$reads_on_ref$start0
                d$start0 <- d$start0 - re_start0
                qnames <- strsplit(d$qnames, ",")[[1]]
                sample_bits <- variant_reads[qnames, sample_bit]
                plotOrder <- order(sample_bits)
                tracks <- initEncodingPlot(plot, d)
                addTargetBases(tracks, d)
                if (input$showClonal) addHaplotypeClonal(tracks, d)
                addReadVariants(tracks, d, plotOrder)
                abline(h = 0:tracks$ymax, lwd = 0.5, col = rgb(0.5,0.5,0.5)) 
            }, error = function(e) NULL)
            stopSpinner(session)
            plot$finishPng(tracks$layout)
        }
    )
    plot
}
readsOnHap1Plot <- readsOnHapPlot("readsOnHap1Plot", "reads_on_hap1")
readsOnHap2Plot <- readsOnHapPlot("readsOnHap2Plot", "reads_on_hap2")

#----------------------------------------------------------------------
# fragment variants tables
#----------------------------------------------------------------------
variantsTableData <- reactive({
    fragment <- req(fragment())
    if (!input$tableShowClonal) {
        fragment$variants <- fragment$variants[clonal == 0]
    }
    if (!input$tableShowIndels) {
        fragment$variants <- fragment$variants[is_indel == 0]
    }
    comp <- c("A" = "T", "C" = "G", "G" = "C", "T" = "A")
    fragment$variants[, context := fcase(
        is_indel == TRUE,               NA_character_,
        is.na(context_base_left),       NA_character_,
        grepl("N", context_base_left),  NA_character_,
        grepl("N", context_base_right), NA_character_,
        tgt_bases %in% c("C", "T"), paste0(
            context_base_left,      
            "[", paste0(     tgt_bases,  ">",      alt_bases), "]",
            context_base_right
        ),
        tgt_bases %in% c("A", "G"), paste0(
            comp[context_base_right],
            "[", paste0(comp[tgt_bases], ">", comp[alt_bases]), "]",
            comp[context_base_left]
        )
    )]
    variants <- fragment$variants[, 
        .(
            pos = ref_pos0 + 1,
            tgt_bases = tgt_bases,
            alt_bases = alt_bases,
            context = context,
            is_indel = is_indel,
            haplotype = haplotype,
            count = n_matching_reads,
            coverage = n_reads,
            sample_bits = sample_bits,
            n_samples = n_samples,
            clonal = clonal,
            matches_clonal = matches_clonal,
            vaf = round(vaf, 3),
            qual = max_min_qual,
            qnames = qnames
        )
    ]
    if (!input$tableShowQnames || input$tableShowClonal) {
        variants$qnames <- NULL
    }
    variants
})
variantsTable <- bufferedTableServer(
    "variantsTable",
    id,
    input,
    variantsTableData,
    selection = 'single',
    selectionFn = function(selectedRows) NULL,
    options = list()
)

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    settings$replace(bm$settings)
    # # updateSelectInput(session, "sampleSet-sampleSet", selected = bm$input[['sampleSet-sampleSet']])
    if(!is.null(bm$outcomes)) {
    #     # outcomes <<- listToReactiveValues(bm$outcomes)
        # sizePlot$settings$replace(bm$outcomes$sizePlotSettings)
        # offsetPlotWide$settings$replace(bm$outcomes$offsetPlotWideSettings)
        # offsetPlotNarrow$settings$replace(bm$outcomes$offsetPlotNarrowSettings)
    }
    bookmarkObserver$destroy()
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    settings = settings$all_,
    outcomes = reactive({ list(
        # sizePlotSettings = sizePlot$settings$all_(),
        # offsetPlotWideSettings = offsetPlotWide$settings$all_(),
        # offsetPlotNarrowSettings = offsetPlotNarrow$settings$all_()
    ) }),
    settingsObject = settings,
    # junctions_filtered = junctions_filtered,
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
