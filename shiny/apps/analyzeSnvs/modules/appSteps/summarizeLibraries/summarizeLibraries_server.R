#----------------------------------------------------------------------
# server components for the summarizeLibraries appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
summarizeLibrariesServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'summarizeLibraries'
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
# load data
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer("dataSource", selection = "single") 
variants <- reactive({
    sourceId <- req(sourceId())
    startSpinner(session, message = "loading variants")
    hf3_getVariants_cached(sourceId)
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

#----------------------------------------------------------------------
# fragment length distribution plot
#----------------------------------------------------------------------
fragLengthPlot <- mdiInteractivePlotBoxServer(
    "fragLength",
    # click = TRUE,
    # brush = TRUE,
    points  = TRUE, # set to TRUE to expose relevant plot options
    lines   = TRUE,
    settings = NULL, # an additional settings template as a list()
    defaults = NULL, # list of default settings values use to inialize settings
    create = function(...) {
        sourceId <- req(sourceId())
        sample_bits <- hf3_sample_bits(sourceId)
        d <- req(reads_on_reference())
        startSpinner(session, message = "rendering length distribution")
        ymax <- 0
        d <- lapply(sample_bits, function(sample_bit){
            dd <- d[bitwAnd(sample_bits, sample_bit) > 0]
            dd <- dd[, 
                .(n_reads = sum(as.integer(unlist(strsplit(sample_bitss, ","))) == sample_bit)), 
                keyby = .(n_bases_bin)
            ]
            dd[, freq := n_reads / sum(n_reads)]
            ymax <<- max(ymax, dd$freq)
            dd
        })
        layout <- fragLengthPlot$initializePng(mar = c(4.1, 4.1, 0.9, 0.9)) %>% 
                  fragLengthPlot$initializeFrame(
            xlim = c(1000, 9000),
            ylim = c(0, ymax * 1.05),
            xlab = "RE Fragment Length",
            ylab = "Frequency (by sample)",
            xaxs = "i",
            yaxs = "i"
        )
        bar_sep <- 250 / length(d)
        offset <- length(d) / 2 * bar_sep
        lapply(1:length(d), function(i){
            dd <- d[[i]]
            fragLengthPlot$addPoints(
                x = dd$n_bases_bin - offset + (i - 1) * bar_sep,
                y = dd$freq,
                typ = "h",
                lwd = 1.5,
                col = i
            )
        })
        stopSpinner(session)
        fragLengthPlot$finishPng(layout)
    }
)

#----------------------------------------------------------------------
# fragment coverage distribution plot
#----------------------------------------------------------------------
fragCoveragePlot <- mdiInteractivePlotBoxServer(
    "fragCoverage",
    # click = TRUE,
    # brush = TRUE,
    points  = TRUE, # set to TRUE to expose relevant plot options
    lines   = TRUE,
    settings = NULL, # an additional settings template as a list()
    defaults = NULL, # list of default settings values use to inialize settings
    create = function(...) {
        d <- req(reads_on_reference())
        startSpinner(session, message = "rendering coverage distribution")
        d <- d[, .(count = .N), keyby = .(n_reads)]
        d[, freq := count / sum(count)]
        layout <- fragCoveragePlot$initializePng(mar = c(4.1, 4.1, 0.9, 0.9)) %>% 
                  fragCoveragePlot$initializeFrame(
            xlim = c(0, d[, max(n_reads)]),
            ylim = c(0, d[, max(freq) * 1.05]),
            xlab = "RE Fragment Coverage",
            ylab = "Frequency",
            xaxs = "i",
            yaxs = "i"
        )
        abline(v = input$minCoverage, col = "blue", lwd = 1.5)
        fragCoveragePlot$addPoints(
            x = d$n_reads,
            y = d$freq,
            typ = "h"
        )
        stopSpinner(session)
        fragCoveragePlot$finishPng(layout)
    }
)

#----------------------------------------------------------------------
# fragment coverage distribution plot
#----------------------------------------------------------------------
lengthVsCoveragePlot <- mdiInteractivePlotBoxServer(
    "lengthVsCoverage",
    # click = TRUE,
    # brush = TRUE,
    points  = TRUE, # set to TRUE to expose relevant plot options
    lines   = TRUE,
    settings = NULL, # an additional settings template as a list()
    defaults = NULL, # list of default settings values use to inialize settings
    create = function(...) {
        d <- req(reads_on_reference())
        startSpinner(session, message = "rendering correlation")
        layout <- lengthVsCoveragePlot$initializePng(mar = c(4.1, 4.1, 0.9, 0.9)) %>% 
                  lengthVsCoveragePlot$initializeFrame(
            xlim = c(1000, 9000),
            ylim = c(0, quantile(d$n_reads, 0.99) * 1.05),
            xlab = "RE Fragment Length",
            ylab = "Read Coverage (all samples)",
            xaxs = "i",
            yaxs = "i"
        )
        d <- d[sample(.N, 5000, replace = FALSE)]
        lengthVsCoveragePlot$addPoints(
            x = jitter(d$n_bases),
            y = jitter(d$n_reads),
            pch = "."
        )
        abline(h = input$minCoverage, col = "blue", lwd = 1.5)
        stopSpinner(session)
        lengthVsCoveragePlot$finishPng(layout)
    }
)

#----------------------------------------------------------------------
# VAF distribution plot
#----------------------------------------------------------------------
vafPlot <- mdiInteractivePlotBoxServer(
    "vafPlot",
    # click = TRUE,
    # brush = TRUE,
    points  = TRUE, # set to TRUE to expose relevant plot options
    lines   = TRUE,
    settings = NULL, # an additional settings template as a list()
    defaults = NULL, # list of default settings values use to inialize settings
    create = function(...) {
        d <- req(variants())
        startSpinner(session, message = "rendering VAF")
        d <- d[
            is_snv == TRUE &
            n_reads >= input$minCoverage &
            clonal == 1
        ]
        d <- d[, .(count = .N), keyby = .(vaf_bin)]
        d[, freq := count / sum(count)]
        layout <- vafPlot$initializePng(mar = c(4.1, 4.1, 0.9, 0.9)) %>% 
                  vafPlot$initializeFrame(
            xlim = c(0, 1),
            ylim = c(0, d[vaf_bin < 1, max(freq) * 1.05]),
            xlab = "Variant Allele Frequency",
            ylab = "Frequency",
            # xaxs = "i",
            yaxs = "i"
        )
        vafPlot$addLines(
            x = d$vaf_bin,
            y = d$freq,
            typ = "h",
            lwd = 1
        )
        stopSpinner(session)
        vafPlot$finishPng(layout)
    }
)

#----------------------------------------------------------------------
# fragment variants tables
#----------------------------------------------------------------------
variantSummaryTableData <- reactive({
    sourceId <- req(sourceId())
    smp_bits <- hf3_sample_bits(sourceId)
    samples  <- hf3_getSampleNames(sourceId, smp_bits, as_string = FALSE)

    startSpinner(session, message = "loading valid subclonal SNVs")
    valid_subclonal_snvs <- hf3_getValidSubclonalSnvs(sourceId)

    startSpinner(session, message = "loading valid haplotypes")
    valid_haplotypes <- hf3_getValidHaplotypes(sourceId)

    startSpinner(session, message = "loading trustworthy haplotypes")
    trustworthy_haplotypes <- hf3_getTrustworthyHaplotypes(
        sourceId, valid_subclonal_snvs, valid_haplotypes
    )

    startSpinner(session, message = "loading trustworthy subclonal SNVs")
    trustworthy_subclonal_snvs <- hf3_getTrustworthySubclonalSnvs(
        sourceId, valid_subclonal_snvs, trustworthy_haplotypes
    )

    startSpinner(session, message = "loading trustworthy haplotype reads")
    trustworthy_haplotype_reads <- hf3_getHaplotypeReads(
        sourceId, trustworthy_haplotypes
    )

    # assemble a table of metric values
    startSpinner(session, message = "building metrics table")
    smp_all <- "all"
    metrics <- c("n_reads", "n_unmasked_bases", "n_subclonal_snvs")
    n_col <- length(metrics)
    d <- data.table(
        metric     = metrics,
        sample_bit = rep(NA,      n_col),
        sample     = rep(smp_all, n_col),
        value = c(
            trustworthy_haplotype_reads[, c(
                .N, 
                sum(n_unmasked_bases) # includes invariant reads in tally
            )],
            nrow(trustworthy_subclonal_snvs)
        )
    )
    for (i in 1:length(smp_bits)){
        d <- rbind(d, data.table(
            metric     = metrics,
            sample_bit = rep(smp_bits[i], n_col),
            sample     = rep(samples[i],  n_col),
            value = c(
                trustworthy_haplotype_reads[sample_bit == smp_bits[i], c(
                    .N, 
                    sum(n_unmasked_bases)
                )],
                nrow(trustworthy_subclonal_snvs[sample_bits == smp_bits[i]]) # since n_samples == 1
            )
        ))
    }

    # cast to one row per sample, metrics in columns
    d <- dcast(d, sample_bit + sample ~ metric)

    # calculate subclonal SNV rates and return the result
    d[, ":="(
        subclonal_snv_rate = sprintf("%.2e", n_subclonal_snvs / n_unmasked_bases)
    )]
    stopSpinner(session)
    d
})
variantSummaryTable <- bufferedTableServer(
    "variantSummaryTable",
    id,
    input,
    variantSummaryTableData,
    selection = 'single',
    selectionFn = function(selectedRows) NULL,
    options = list()
)

#----------------------------------------------------------------------
# mutation signature modeling
#----------------------------------------------------------------------
sigFitShared <- reactive({
    sourceId <- req(sourceId())

    startSpinner(session, message = "loading mut sig data")
    valid_subclonal_snvs <- hf3_getValidSubclonalSnvs(sourceId)
    valid_haplotypes <- hf3_getValidHaplotypes(sourceId)
    trustworthy_haplotypes <- hf3_getTrustworthyHaplotypes(
        sourceId, valid_subclonal_snvs, valid_haplotypes
    )
    trustworthy_subclonal_snvs <- hf3_getTrustworthySubclonalSnvs(
        sourceId, valid_subclonal_snvs, trustworthy_haplotypes
    )

    startSpinner(session, message = "casting mutation contexts")
    trustworthy_subclonal_snvs[, sample := paste0("smp_", sample_bits)]
    mutation_counts <- trustworthy_subclonal_snvs[
        !is.na(context), 
        .N, 
        keyby = .(context, sample)
    ] %>% dcast(context ~ sample, value.var = "N")

    startSpinner(session, message = "filling out mutation contexts")
    all_contexts <- character()
    for (mutation in c("[C>A]","[C>G]","[C>T]","[T>A]","[T>C]","[T>G]")){
        for (left_context_base in c("A","C","G","T")){
            for (right_context_base in c("A","C","G","T")){
                all_contexts <- c(all_contexts, paste0(
                    left_context_base, 
                    mutation, 
                    right_context_base
                ))
            }
        }
    }
    mutation_counts <- merge(
        data.table(context = all_contexts),
        mutation_counts,
        by = "context",
        all.x = TRUE,
        all.y = FALSE,
        sort = FALSE
    )

    startSpinner(session, message = "summing sample counts")
    sample_columns  <- names(mutation_counts)[  names(mutation_counts) !=     "context"]
    treated_columns <- names(mutation_counts)[!(names(mutation_counts) %in% c("context", "smp_1"))]
    mutation_counts[, all := apply(.SD, 1, sum, na.rm = TRUE), .SDcols = sample_columns]
    mutation_counts[, tx  := apply(.SD, 1, sum, na.rm = TRUE), .SDcols = treated_columns]
    mutation_counts[, ctl := smp_1]
    mutation_counts <- mutation_counts[, .(context, ctl, tx, all)]

    startSpinner(session, message = "parsing mutation matrix")
    mutation_count_matrix <- as.matrix(mutation_counts[, -1, with = FALSE])
    rownames(mutation_count_matrix) <- mutation_counts$context

    startSpinner(session, message = "filling missing counts")
    mutation_count_matrix[is.na(mutation_count_matrix)] <- 0

    startSpinner(session, message = "fetching signatures")
    cosmic_signatures <- MutationalPatterns::get_known_signatures(muttype = "snv")
    rownames(cosmic_signatures) <- rownames(mutation_count_matrix)
    cosmic_sigs <- c("SBS1", "SBS2", "SBS5", "SBS13", "SBS18")
    cosmic_signatures <- cosmic_signatures[, cosmic_sigs]

    neighbor_sigs <- fread(file.path(gitStatusData$suite$dir, "resources", "neighbor_sigs.csv"))
    neighbor_signatures <- cbind(cosmic_signatures, as.matrix(neighbor_sigs[, -1, with = FALSE]))
    colnames(neighbor_signatures) <- c(cosmic_sigs, names(neighbor_sigs)[-1])

    list(
        mutation_counts       = mutation_counts,
        mutation_count_matrix = mutation_count_matrix,
        cosmic_signatures     = cosmic_signatures,
        neighbor_signatures   = neighbor_signatures
    )
})
fit_to_signatures <- function(signatures){
    sigFitShared <- req(sigFitShared())

    startSpinner(session, message = paste("fitting", signatures))
    signatures_fit <- MutationalPatterns::fit_to_signatures(
        sigFitShared$mutation_count_matrix, 
        sigFitShared[[signatures]]
    )

    relative_contributions <- prop.table(signatures_fit$contribution, margin = 2)
    order <- order(-relative_contributions[,"all"])
    print(relative_contributions[order,])

    fit <- signatures_fit$reconstructed
    mutation_fracs <- sigFitShared$mutation_counts[, .(
        context = context,
        ctl = round(ctl / sum(ctl, na.rm = TRUE), 5),
        tx  = round(tx  / sum(tx,  na.rm = TRUE),  5),
        all = round(all / sum(all, na.rm = TRUE), 5),
        ctl_fit = round(fit[,"ctl"] / sum(fit[,"ctl"], na.rm = TRUE), 5),
        tx_fit  = round(fit[,"tx"]  / sum(fit[,"tx"],  na.rm = TRUE), 5),
        all_fit = round(fit[,"all"] / sum(fit[,"all"], na.rm = TRUE), 5)
    )]
    mutation_fracs[is.na(mutation_fracs)] <- 0
    print(mutation_fracs)  

    stopSpinner(session)
    list(
        shared = sigFitShared,
        fit    = signatures_fit
    )
}
sigFitWithoutNeighborPlot <- staticPlotBoxServer(
    "sigFitWithoutNeighbor",
    maxHeight = "600px",
    create = function() {
        fit <- fit_to_signatures("cosmic_signatures")
        plot <- MutationalPatterns::plot_compare_profiles(
            fit$shared$mutation_count_matrix[, "all"],
            fit$fit$reconstructed[, "all"],
            profile_names = c("Original", "Fitted"),
            condensed = TRUE,
            profile_ymax = 0.1
        )
        print(plot)
    }
)
sigFitWithNeighborPlot <- staticPlotBoxServer(
    "sigFitWithNeighbor",
    maxHeight = "600px",
    create = function() {
        fit <- fit_to_signatures("neighbor_signatures")
        plot <- MutationalPatterns::plot_compare_profiles(
            fit$shared$mutation_count_matrix[, "all"],
            fit$fit$reconstructed[, "all"],
            profile_names = c("Original", "Fitted"),
            condensed = TRUE,
            profile_ymax = 0.1
        )
        print(plot)
    }
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
        fragLengthPlot$settings$replace(bm$outcomes$fragLengthPlotSettings)
        fragCoveragePlot$settings$replace(bm$outcomes$fragCoveragePlotSettings)
        lengthVsCoveragePlot$settings$replace(bm$outcomes$lengthVsCoveragePlotSettings)
        vafPlot$settings$replace(bm$outcomes$vafPlotSettings)
        sigFitWithoutNeighborPlot$settings$replace(bm$outcomes$sigFitWithoutNeighborPlotSettings)
        sigFitWithNeighborPlot$settings$replace(bm$outcomes$sigFitWithNeighborPlotSettings)
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
        fragLengthPlotSettings = fragLengthPlot$settings$all_(),
        fragCoveragePlotSettings = fragCoveragePlot$settings$all_(),
        lengthVsCoveragePlotSettings = lengthVsCoveragePlot$settings$all_(),
        vafPlotSettings = vafPlot$settings$all_(),
        sigFitWithoutNeighborPlotSettings = sigFitWithoutNeighborPlot$settings$all_(),
        sigFitWithNeighborPlotSettings = sigFitWithNeighborPlot$settings$all_()
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
