library(shiny)
library(DT)
library(plyr)
library(ggplot2)
library(queryBuildR)
library(shinyBS)
library(rpivotTable)
library(httr)
library(shinyjs)
require(jsonlite)
library(RCurl)

source("filterPhenotypes.R")
source("filterVariants.R")

shinyServer(function(input, output,session) {
  sessionvalues <- reactiveValues()
  data<-loadData("")
  sessionvalues$variants<-data$data
  sessionvalues$nbRowsExceededWarningMessage<-data$nbRowsExceededWarningMessage
  sessionvalues$phenotypes<-loadPhenotypes("")
  
  sessionvalues$analysesNames<-analysesNames
  sessionvalues$variantDataGene<-NULL
  
  sessionvalues$logged_user<-""
  
  #Add UI components for filtering phenotypes/variant tabs
  output<-createFilterPhenotype(input,output,session,sessionvalues)
  output<-createFilterVariant(input,output,session,sessionvalues)
  
  observe({
    query<-input$sampleid
    if (length(query)>0) {
      query<-substr(query,2,nchar(query))
      if (query!="") {
        updateTabsetPanel(session, "tabset", selected = "Gene & variant filtering manager")
      }
    }
  })
  
  
  ####################################################
  #User login and UI controls
  ####################################################
  
  disableButtons<-function() {
    shinyjs::disable("filterPhenotypeLoad")
    shinyjs::disable("filterPhenotypeSave")
    shinyjs::disable("filterPhenotypeDelete")
    shinyjs::disable("filterVariantLoad")
    shinyjs::disable("filterVariantSave")
    shinyjs::disable("filterVariantDelete")
    shinyjs::disable("startAnalysisButton")
  }
  
  enableButtons<-function() {
    shinyjs::enable("filterPhenotypeLoad")
    shinyjs::enable("filterPhenotypeSave")
    shinyjs::enable("filterPhenotypeDelete")
    shinyjs::enable("filterVariantLoad")
    shinyjs::enable("filterVariantSave")
    shinyjs::enable("filterVariantDelete")
    shinyjs::enable("startAnalysisButton")
  }
  
  output$loginUI<-renderUI({
    if (length(sessionvalues$logged_user)>0) {
      if (sessionvalues$logged_user=="") {
        disableButtons()
        div(div(
          strong("Not logged in. "),
          actionButton("loginButton", icon("log-in","fa-2x",lib="glyphicon"),tooltip="Log in"),
          bsAlert("alertLogin"),
          align="right"),
          bsModal("modalLogin", "Login", "loginButton", 
                  size = "small",
                  textInput("userIDLogin", "User ID:", value = ""),
                  passwordInput("passwordLogin", "Password:", value = ""),
                  actionButton("confirmLogin", label = "Log in")
          ))
      }
      else {
        enableButtons()
        div(div(
          strong(paste0("Logged in as ",sessionvalues$logged_user)),
          actionButton("logoutButton", icon("log-out","fa-2x",lib="glyphicon")),  
          align="right"))
      }
    }
  })
  
  observe({
    input$filterPhenotypeSave
    input$filterPhenotypeDelete
    input$filterVariantSave
    input$filterVariantDelete
    if (sessionvalues$logged_user=="") {
      disableButtons()
    }
    else {
      enableButtons()
    }
  })
  
  observe({
    if (length(input$confirmLogin)>0) {
      if (input$confirmLogin!=0) {
        isolate({
          users<-read.table("users.csv",header=T,stringsAsFactors=F,colClasses=c("character","character"))
          i.users<-which(users$login==input$userIDLogin)
          error_msg<-""
          if (length(i.users)>0) {
            if (input$passwordLogin==users$password[i.users]) {
              sessionvalues$logged_user<-input$userIDLogin
              runjs("$('.modal-backdrop').remove();")
              runjs('$("body").removeClass("modal-open")')
            }
            else {
              error_msg<-"Incorrect password"
            }
          }
          else {
            error_msg<-"Unknown login"
          }
          if (error_msg!="") {
            toggleModal(session, "modalLogin", toggle = "close")
            createAlert(session, "alertLogin", title = "Oops",
                        content = error_msg, append = FALSE)
          }
        })
      }
    }
  }) 
  
  observe({
    if (length(input$logoutButton)>0) {
      if (input$logoutButton) 
        sessionvalues$logged_user<-""
    }
  }) 
  
  ####################################################
  #Phenotypes group manager
  ####################################################
  
  
  #Get sample IDs button
  output$listSamplesIDs<-renderText({
    listIDs<-sessionvalues$phenotypes$Sample_ID
    result<-paste(listIDs,sep="",collapse=" , ")
    result
  })
  
  #If apply filter, update phenotype table
  observe({
    if (length(input$filterPhenotypeQueryBuilderSQL)>0)
      sessionvalues$phenotypes<-loadPhenotypes(input$filterPhenotypeQueryBuilderSQL)
  })
  
  output$showVarPhenotypesUI<-renderUI({
    isolate({
      niceNames<-as.vector(sapply(colnames(sessionvalues$phenotypes),idToName))
      selectInput('showVarPhenotypes', 'Select variables to display', niceNames, 
                  selected=niceNames,multiple=TRUE, selectize=TRUE,width='1050px')
    })
  })
  
  output$phenotypesTable<-DT::renderDataTable({
    if (length(input$showVarPhenotypes)>0) {
      data<-sessionvalues$phenotypes[,sapply(input$showVarPhenotypes,nameToId)]
      data[is.na(data)]<-''
      colnames(data)<-input$showVarPhenotypes
      getWidgetTable(data,session)
    }
  },server=T)
  
  output$downloadPhenotypesSelection <- downloadHandler(
    filename = function() {
      paste('phenotypes.zip', sep='')
    },
    content = function(con) {
      write.csv(sessionvalues$phenotypes, file="phenotypes.csv", row.names=F,quote=T)
      zip(con,c('phenotypes.csv'))
    }
  )
  
  output$pivotTablePhenotypes<-renderRpivotTable({
    if (length(input$showVarPhenotypes)>0) {
      data<-sessionvalues$phenotypes[,sapply(input$showVarPhenotypes,nameToId)]
      colnames(data)<-input$showVarPhenotypes
      rpivotTable(data,width='1050px')
    }
  })
  
  ####################################################
  #Variant group manager
  ####################################################
  
  output$nbRowsExceededWarningMessage<-renderText({
    sessionvalues$nbRowsExceededWarningMessage
  })
  
  #Apply filters
  observe({
    if (length(input$filterVariantQueryBuilderSQL)) {
      withProgress(min=1, max=3, expr={
        setProgress(message = 'Retrieving data, please wait...',
                    value=2)
        data<-loadData(input$filterVariantQueryBuilderSQL)
        sessionvalues$variants<-data$data
        sessionvalues$nbRowsExceededWarningMessage<-data$nbRowsExceededWarningMessage
      })
    }
  })
  
  output$showVarVariantsUI<-renderUI({
    isolate({
      niceNames<-as.vector(sapply(colnames(sessionvalues$variants),idToName))
      niceNames<-colnames(sessionvalues$variants)
      selectInput('showVarVariants', 'Select variables to display', niceNames, 
                  selected=niceNames[c(1,2:7)],multiple=TRUE, selectize=TRUE,width='1050px')
    })
  })
  
  output$variantsTable<-DT::renderDataTable({
    if (length(input$showVarVariants)>0) {
      data<-sessionvalues$variants[,input$showVarVariants]
      data[is.na(data)]<-''
      colnames(data)<-input$showVarVariants
      getWidgetTable(data,session)
    }
  },server=T)
  
  output$downloadVariantsSelection <- downloadHandler(
    filename = function() {
      paste('variantSelection.zip', sep='')
    },
    content = function(con) {
      write.csv(sessionvalues$variants, file="variantSelection.csv", row.names=F,quote=T)
      zip(con,c('variantSelection.csv'))
    }
  )
  output$pivotTableVariants<-renderRpivotTable({
    if (length(input$showVarVariants)>0) {
      data<-sessionvalues$variants[,sapply(input$showVarVariants,nameToId)]
      colnames(data)<-input$showVarVariants
      rpivotTable(data,width='1050px')
    }
  })
  
  ####################################################
  #Scoring tool
  ####################################################
  
  output$selectSampleControlUI<-renderUI({
    input$filterVariantConfirmSave
    variantsGroup<-read.table(sessionvalues$variantGroupFile,header=T,stringsAsFactors=F,colClasses=c("character","character"))
    selectInput('selectSampleControl', 'Control group', 
                choices = list("Groups"=variantsGroup$Name), 
                selected=variantsGroup$Name[1],
                selectize = FALSE)
  })
  
  output$selectSampleCaseUI<-renderUI({
    input$filterVariantConfirmSave
    variantsGroup<-read.table(sessionvalues$variantGroupFile,header=T,stringsAsFactors=F,colClasses=c("character","character"))
    selectInput('selectSampleCase', 'Case group', 
                choices = list("Groups"=variantsGroup$Name), 
                selected=variantsGroup$Name[1],
                selectize = FALSE)
  })
  
  observe({
    if (input$startAnalysisButton>0) {
      isolate({
        analysisName<-input$analysisName
        scope<-input$rankingScope
        scale<-input$rankingScale
      })
      
      if (analysisName!='') {
        analysis<-list()
        variantsGroup<-read.table(sessionvalues$variantGroupFile,header=T,stringsAsFactors=F,colClasses=c("character","character"))
        
        isolate({
          sampleControlName<-input$selectSampleControl
          sampleCaseName<-input$selectSampleCase
          
          selectSampleControlIndex<-which(variantsGroup$Name==sampleControlName)
          sqlControl<-variantsGroup$SQL[selectSampleControlIndex]
          sqlControl<-preprocSQL(sqlControl)
          
          selectSampleCaseIndex<-which(variantsGroup$Name==sampleCaseName)
          sqlCase<-variantsGroup$SQL[selectSampleCaseIndex]
          sqlCase<-preprocSQL(sqlCase)
          
          controlMAF<-input$controlGroupMAF
          caseMAF<-input$caseGroupMAF
          
          scoringFunction<-input$rankingCriterion
          
          jobArguments<-rbind(analysisName,scope,scale,sqlControl,sqlCase,sampleControlName,sampleCaseName,controlMAF,caseMAF,VARIANTS,scoringFunction)
          setwd("spark")
          write.table(file="jobsArguments.conf",jobArguments,quote=F,col.names=F,row.names=F)
          system(paste0("./run_local.sh ",analysisName," ../users/analyses"))
          setwd("..")
        })
      }  
      else {
        createAlert(session, "alertStartAnalysis", title = "Oops",
                    content = "Specify a name for your analysis", append = FALSE)
      }
    }
  })
  
  
  ####################################################
  #Results explorer
  ####################################################
  
  #Refreshing list of analyses
  observe({
    input$refreshResultsButton
    sessionvalues$analysesNames<-getAnalysesNames(sessionvalues$logged_user)
  })
  
  #UI widget for selecting analysis
  output$selectAnalysisUI<-renderUI({
    selectInput('selectAnalysis', 'Select analysis', choices = sessionvalues$analysesNames
                , selected=sessionvalues$analysesNames[1],selectize = FALSE)
  })
  
  #Display selected variables for an analysis
  output$showVarResultsUI<-renderUI({
    nameAnalysis<-input$selectAnalysis
    
    if (length(sessionvalues$results)>0) {
      
      niceNames<-as.vector(sapply(colnames(sessionvalues$results$scoreSummary),idToName))
      initialSelect<-niceNames
      
      selectInput('showVarResults', 'Select variables to display', niceNames, 
                  selected=initialSelect,multiple=TRUE, selectize=TRUE,width='1050px')
    }
  })
  
  #Display ranking (table) for an analysis
  output$resultsTable<-DT::renderDataTable({
    if ((length(input$showVarResults)>0) & (length(sessionvalues$results)>0)) {
      nameAnalysis<-input$selectAnalysis
      if (length(setdiff(sapply(input$showVarResults,nameToId),colnames(sessionvalues$results$scoreSummary)))==0) {
        data<-sessionvalues$results$scoreSummary[,sapply(input$showVarResults,nameToId)]
        colnames(data)<-input$showVarResults
        getWidgetTable(data,session,selection='single')
      }
    }
  },server=TRUE)
  
  #Load data for an analysis
  observe({
    nameAnalysis<-input$selectAnalysis
    if (length(nameAnalysis)>0) {
      sessionvalues$results<-procRes(folderAnalyses,nameAnalysis)
      if (length(input$resultsTable_rows_selected)) {
        withProgress(min=1, max=4, expr={
          setProgress(message = 'Retrieving control data, please wait...',
                      value=2)
          
          if (sessionvalues$results$scale=="variant") {
            if (sessionvalues$results$scope=="monogenic") {
              
              variantsData<-sessionvalues$results$scoreSummary[input$resultsTable_rows_selected,'Variant_ID']
              variantsData<-strsplit(as.character(variantsData),":")[[1]]
              
              sqlControl<-sessionvalues$results$sqlControl
              sqlControl<-paste0(sqlControl," and chr='",variantsData[1],"'"," and pos=",variantsData[2]," and ref='",variantsData[3],"'"," and alt='",variantsData[4],"'")
              
              variantsControl<-loadData(sqlControl,noLimit=T,preproc=F)$data
              nControl<-nrow(variantsControl)
              
              setProgress(message = 'Retrieving case data, please wait...',
                          value=3)
              
              sqlCase<-sessionvalues$results$sqlCase
              sqlCase<-paste0(sqlCase," and chr='",variantsData[1],"'"," and pos=",variantsData[2]," and ref='",variantsData[3],"'"," and alt='",variantsData[4],"'")
              variantsCase<-loadData(sqlCase,noLimit=T,preproc=F)$data
              variants<-rbind(variantsControl,variantsCase)
              variants<-cbind("Group"=c(rep(sessionvalues$results$controlGroupName,nControl),rep(sessionvalues$results$caseGroupName,nrow(variantsCase))),variants)
              sessionvalues$variantDataGene<-variants
            }
            
            if (sessionvalues$results$scope=="digenic") {
              
              variantsData1<-sessionvalues$results$scoreSummary[input$resultsTable_rows_selected,'Variant_ID1']
              variantsData1<-strsplit(as.character(variantsData1),":")[[1]]
              variantsData2<-sessionvalues$results$scoreSummary[input$resultsTable_rows_selected,'Variant_ID2']
              variantsData2<-strsplit(as.character(variantsData2),":")[[1]]
              
              sqlControl<-sessionvalues$results$sqlControl
              sqlControl<-paste0(sqlControl," and ((chr='",variantsData1[1],"'"," and pos=",variantsData1[2]," and ref='",variantsData1[3],"'"," and alt='",variantsData1[4],"') or (chr='",variantsData2[1],"'"," and pos=",variantsData2[2]," and ref='",variantsData2[3],"'"," and alt='",variantsData2[4],"'))")
              
              variantsControl<-loadData(sqlControl,noLimit=T,preproc=F)$data
              nControl<-nrow(variantsControl)
              
              setProgress(message = 'Retrieving case data, please wait...',
                          value=3)
              
              sqlCase<-sessionvalues$results$sqlCase
              sqlCase<-paste0(sqlCase," and ((chr='",variantsData1[1],"'"," and pos=",variantsData1[2]," and ref='",variantsData1[3],"'"," and alt='",variantsData1[4],"') or (chr='",variantsData2[1],"'"," and pos=",variantsData2[2]," and ref='",variantsData2[3],"'"," and alt='",variantsData2[4],"'))")
              variantsCase<-loadData(sqlCase,noLimit=T,preproc=F)$data
              variants<-rbind(variantsControl,variantsCase)
              variants<-cbind("Group"=c(rep(sessionvalues$results$controlGroupName,nControl),rep(sessionvalues$results$caseGroupName,nrow(variantsCase))),variants)
              sessionvalues$variantDataGene<-variants
            }
          }
          
          if (sessionvalues$results$scale=="gene") {
            if (sessionvalues$results$scope=="monogenic") {
              geneID<-sessionvalues$results$scoreSummaryRaw[input$resultsTable_rows_selected,'Gene_Symbol']
              sqlControl<-sessionvalues$results$sqlControl
              sqlControl<-paste0(sqlControl," and gene_symbol='",geneID,"'")
              variantsControl<-loadData(sqlControl,noLimit=T,preproc=F)$data
              nControl<-nrow(variantsControl)
              
              setProgress(message = 'Retrieving case data, please wait...',
                          value=3)
              sqlCase<-sessionvalues$results$sqlCase
              sqlCase<-paste0(sqlCase," and gene_symbol='",geneID,"'")
              variantsCase<-loadData(sqlCase,noLimit=T,preproc=F)$data
              variants<-rbind(variantsControl,variantsCase)
              variants<-cbind("Group"=c(rep(sessionvalues$results$controlGroupName,nControl),rep(sessionvalues$results$caseGroupName,nrow(variantsCase))),variants)
              sessionvalues$variantDataGene<-variants
            }
            if (sessionvalues$results$scope=="digenic") {
              geneID1<-sessionvalues$results$scoreSummaryRaw[input$resultsTable_rows_selected,'Gene_Symbol1']
              geneID2<-sessionvalues$results$scoreSummaryRaw[input$resultsTable_rows_selected,'Gene_Symbol2']
              sqlControl<-sessionvalues$results$sqlControl
              sqlControl<-paste0(sqlControl," and (gene_symbol='",geneID1,"' or gene_symbol='",geneID2,"')")
              variantsControl<-loadData(sqlControl,noLimit=T,preproc=F)$data
              nControl<-nrow(variantsControl)
              
              setProgress(message = 'Retrieving case data, please wait...',
                          value=3)
              sqlCase<-sessionvalues$results$sqlCase
              sqlCase<-paste0(sqlCase," and (gene_symbol='",geneID1,"' or gene_symbol='",geneID2,"')")
              variantsCase<-loadData(sqlCase,noLimit=T,preproc=F)$data
              variants<-rbind(variantsControl,variantsCase)
              variants<-cbind("Group"=c(rep(sessionvalues$results$controlGroupName,nControl),rep(sessionvalues$results$caseGroupName,nrow(variantsCase))),variants)
              sessionvalues$variantDataGene<-variants
            }
          }
        }
        )
      }
    }
  })
  
  #Display UI widget for selecting variant metadata variables
  output$showVarMetadataUI<-renderUI({
    nameAnalysis<-input$selectAnalysis
    if (length(sessionvalues$variantDataGene)>0) {
      niceNames<-as.vector(sapply(colnames(sessionvalues$variantDataGene),idToName))
      selectInput('showVarMetadata', 'Select variables to display', niceNames, 
                  selected=niceNames[c(1:8,24:26,17)],multiple=TRUE, selectize=TRUE,width='1050px')
    }
  })
  
  #Display variant metadata (as a table)
  output$variantsMetadataTable<-DT::renderDataTable({
    nameAnalysis<-input$selectAnalysis
    if (length(sessionvalues$variantDataGene)>0 & length(input$showVarMetadata)>0) {
      data<-sessionvalues$variantDataGene[,sapply(input$showVarMetadata,nameToId)]
      getWidgetTable(data,session)
    }
  },server=TRUE)
  
  #Display variant metadat (as a pivot table)
  output$variantsMetadataPivotTable<-renderRpivotTable({
    nameAnalysis<-input$selectAnalysis
    if (length(sessionvalues$variantDataGene)>0 & length(input$showVarMetadata)>0) {
      data<-sessionvalues$variantDataGene[,sapply(input$showVarMetadata,nameToId)]
      rpivotTable(data)
    }
  })
  
  #Display analysis metadata
  output$resultsMetadata<-renderUI({
    fluidRow(
      column(3,
             strong("Control group: "),br(),
             strong("Pathological group: "),br(),
             strong("Start time: "),br(),
             strong("End time: "),br(),
             strong("Total run time: ")
      ),
      column(4,
             paste0(sessionvalues$results$group1name," (n=",length(sessionvalues$results$controlSampleID),")"),br(),
             paste0(sessionvalues$results$group2name," (n=",length(sessionvalues$results$caseSampleID),")"),br(),
             sessionvalues$results$start_time,br(),
             sessionvalues$results$end_time,br(),
             paste(sessionvalues$results$run_time, "seconds")
      )
      
    )
  })
  
  #Download score list as CSV
  output$downloadScoreTable <- downloadHandler(
    filename = function() {
      paste0(input$selectAnalysis,'.zip')
    },
    content = function(con) {
      filename=paste0(input$selectAnalysis,'.csv')
      write.csv(sessionvalues$results$scoreSummaryRaw, file=filename, row.names=F,quote=F)
      zip(con,c(filename))
    }
  )
  
  #UI for results
  output$resultsPanel<-renderUI({
    fluidRow(
      column(12,
             uiOutput("resultsMetadata"),
             hr(),
             h3("Scoring results"),
             div(
               downloadButton('downloadScoreTable', label = "Download score table (CSV)",class = NULL),
               align="right"),
             uiOutput("showVarResultsUI"),
             DT::dataTableOutput('resultsTable'),
             hr(),
             fluidRow(
               column(12,
                      fluidRow(
                        uiOutput("showVarMetadataUI"),
                        dataTableOutput('variantsMetadataTable')
                      ),
                      fluidRow(
                        h4("Pivot Table"),
                        rpivotTableOutput('variantsMetadataPivotTable')
                      )
               )
             )
      )
    )
    
  })
  
})

