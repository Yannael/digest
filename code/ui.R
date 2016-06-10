library(shiny)
library(DT)
library(queryBuildR)
library(shinyBS)
library(shinyjs)

shinyUI(
  fluidPage(
    useShinyjs(),
    includeCSS('www/style.css'),
    div(
      fluidRow(
        img(src="mgbck.jpg", height = 150, width = 1000)
      ),
      tags$script("$(document).ready(function() {
                    sampleid=window.location.search.split('?sample_id=')
                    $('#sampleid').val(sampleid)
                });"),
      tags$input(id = 'sampleid', type = 'text', style = 'display:none;'),
      tags$div(class="extraspace2"),
      fluidRow(
        shiny::column(12,
                      uiOutput("loginUI"),
                      tabsetPanel(id="tabset",
                                  tabPanel("Phenotype manager", 
                                           tags$div(class="extraspace2"),
                                           fluidRow(
                                             column(12,
                                                    uiOutput("filterPhenotype"),
                                                    div(
                                                      actionButton("getIDButtonPhenotypesGroup", label = "Get sample IDs"),
                                                      downloadButton('downloadPhenotypesSelection', label = "Download selection (CSV)",class = NULL),
                                                      align="right"),
                                                    bsModal("getIDPhenotypesGroup", "List of sample IDs", "getIDButtonPhenotypesGroup", 
                                                            size = "large",textOutput('listSamplesIDs')
                                                    ),
                                                    uiOutput("showVarPhenotypesUI"),
                                                    dataTableOutput('phenotypesTable'),
                                                    hr(),
                                                    h5(strong("Pivot table")),
                                                    rpivotTableOutput("pivotTablePhenotypes"),
                                                    tags$div(class="extraspace1")
                                             )
                                           )
                                  ),
                                  tabPanel("Gene & variant filtering manager", 
                                           tags$div(class="extraspace2"),
                                           fluidRow(
                                             column(12,
                                                    uiOutput("filterVariant"),
                                                    div(downloadButton('downloadVariantsSelection', label = "Download selection (CSV)",class = NULL),
                                                        align="right"),
                                                    h5(htmlOutput("nbRowsExceededWarningMessage")),
                                                    uiOutput("showVarVariantsUI"),
                                                    dataTableOutput('variantsTable'),
                                                    hr(),
                                                    h5(strong("Pivot table")),
                                                    rpivotTableOutput("pivotTableVariants"),
                                                    tags$div(class="extraspace1")
                                             )
                                           )
                                  ),
                                  tabPanel("DiGeST launcher", 
                                           tags$div(class="extraspace2"),
                                           fluidRow(
                                             column(3,
                                                    h3("1) Variants group(s)"),
                                                    uiOutput("selectSampleControlUI"),
                                                    textInput("caseGroupMAF","MAF threshold","0.5"),
                                                    uiOutput("selectSampleCaseUI"),
                                                    textInput("controlGroupMAF","MAF threshold","0.5")
                                             ),
                                             column(4,offset=1,
                                                    h3("2) Scoring parameters"),
                                                    radioButtons("rankingScale", "Ranking scale",
                                                                 c("Gene" = "gene",
                                                                   "Variant" = "variant"
                                                                 )),
                                                    radioButtons("rankingScope", "Scope",
                                                                 c("Monogenic" = "monogenic",
                                                                   "Digenic" = "digenic"
                                                                 ))
                                                    ,
                                                    radioButtons("rankingCriterion", "Scoring function",
                                                                 c("Sum patients" = "sumpatients",
                                                                   "Sum alleles" = "sumalleles"
                                                                 ),
                                                                 selected=c("sumpatients"))
                                             ),
                                             column(3,
                                                    h3("3) Results collection"),
                                                    textInput("analysisName","Analysis name",""),
                                                    bsAlert("alertStartAnalysis")
                                             )
                                           ),
                                           hr(),
                                           fluidRow(
                                             column(2,offset=5,
                                                    div(actionButton("startAnalysisButton","Start analysis",class="btn btn-primary",disabled = TRUE),align="center"),
                                                    tags$div(class="extraspace1")
                                             )
                                           )
                                  ),
                                  tabPanel("Results explorer", 
                                           tags$div(class="extraspace2"),
                                           fluidRow(
                                             column(3,
                                                    uiOutput("selectAnalysisUI"),
                                                    actionButton("refreshResultsButton","Refresh")
                                             )
                                           ),
                                           hr(),
                                           fluidRow(
                                             column(12,
                                                    uiOutput("resultsPanel"),
                                                    tags$div(class="extraspace1")
                                             )
                                             
                                           )
                                  )
                      )
        )
      ),
      textOutput("done")
    )
  )
  
)



