#----------------------------------------------------------------------
# UI components for the summarizeLibraries appStep module
#----------------------------------------------------------------------

# module ui function
summarizeLibrariesUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$summarizeLibraries)

    # return the UI contents
    standardSequentialTabItem(

        # page header text
        options$longLabel,
        options$leaderText,

        # page header links, uncomment as needed
        id = id,
        # documentation = TRUE,
        # terminal = TRUE,
        console = serverEnv$IS_DEVELOPER,
        code = serverEnv$IS_DEVELOPER,
        settings = TRUE,

        # appStep UI elements, populate as needed
        dataSourceTableUI(
            ns("dataSource"), 
            "Data Source", 
            width = 12, 
            collapsible = TRUE
        ),
        fluidRow(
            box(
                title = NULL,
                width = 12,
                solidHeader = FALSE,
                # status = "primary",
                collapsible = FALSE,
                column(
                    width = 2,
                    numericInput(
                        ns("minCoverage"), 
                        "Min Coverage", 
                        value = 10, 
                        min = 0, 
                        max = 50,
                        step = 5
                    )
                ),
                NULL
            )
        ),
        fluidRow(
            mdiInteractivePlotBoxUI(
                ns("fragLength"), 
                "Fragment Length Distribution",
                width = 6, 
                collapsible = TRUE, 
                collapsed = TRUE
            ),
            mdiInteractivePlotBoxUI(
                ns("fragCoverage"), 
                "Fragment Coverage Distribution",
                width = 6, 
                collapsible = TRUE, 
                collapsed = TRUE
            )
        ),
        fluidRow(
            mdiInteractivePlotBoxUI(
                ns("lengthVsCoverage"), 
                "Fragment Length vs. Coverage",
                width = 6, 
                collapsible = TRUE, 
                collapsed = TRUE
            ),
            mdiInteractivePlotBoxUI(
                ns("vafPlot"), 
                "Variant VAF Distribution",
                width = 6, 
                collapsible = TRUE, 
                collapsed = TRUE
            )
        ),
        fluidRow(
            bufferedTableUI(
                ns("variantSummaryTable"),
                title = "VariantSummary Table",
                width = 12, 
                collapsible = TRUE, 
                collapsed = TRUE
            )
        ),
        fluidRow(
            staticPlotBoxUI(
                ns("cosmicAllVsFittedPlot"),
                "All vs. Fitted",
                width = 6, 
                collapsible = TRUE, 
                collapsed = TRUE
            ),
            staticPlotBoxUI(
                ns("cosmicCtlVsFittedPlot"),
                "Control vs. Fitted",
                width = 6, 
                collapsible = TRUE, 
                collapsed = TRUE
            )
        ),
        fluidRow(
            staticPlotBoxUI(
                ns("cosmicTxVsFittedPlot"),
                "Treated vs. Fitted",
                width = 6, 
                collapsible = TRUE, 
                collapsed = TRUE
            )
        ),
        fluidRow(
            staticPlotBoxUI(
                ns("snvCtlVsTxPlot"),
                "Control vs. Treated, SNVs",
                width = 6, 
                collapsible = TRUE, 
                collapsed = FALSE
            ),
            staticPlotBoxUI(
                ns("indelCtlVsTxPlot"),
                "Control vs. Treated, 1-base indels",
                width = 6, 
                collapsible = TRUE, 
                collapsed = FALSE
            )
        ),
        NULL
    )
}
