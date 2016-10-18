createFilterVariant<-function(input,output,session,sessionvalues) {
  
  observe({
    sessionvalues$variantGroupFile<-paste0("users/filterVariant.csv")
    variantGroupFolder<-paste0("users/")
    if (!file.exists(variantGroupFolder)) {
      dir.create(variantGroupFolder)
    }
    if (!file.exists(sessionvalues$variantGroupFile)) {
      data<-rbind(c("All",""))
      colnames(data)<-c("Name","SQL")
      write.table(file=sessionvalues$variantGroupFile,data,row.names=F,sep="\t")
    }
  })
  
  output$filterVariantQueryBuilderWidget<-renderQueryBuildR({
    load("filterVariantSpec.Rdata")
    rules<-NULL
    queryBuildR(rules,filters)
  })
  
  output$filterVariantSelectLoadUI<-renderUI({
    input$filterVariantConfirmSave
    input$filterVariantConfirmDelete
    data<-read.table(sessionvalues$variantGroupFile,header=T,stringsAsFactors=F,colClasses=c("character","character"))
    selectInput('filterVariantSelectLoad', 'Select filter', 
                choices = data$Name, 
                selected=data$Name[1],
                selectize = FALSE)
  })
  
  output$filterVariantSelectDeleteUI<-renderUI({
    input$filterVariantConfirmSave
    input$filterVariantConfirmDelete
    data<-read.table(sessionvalues$variantGroupFile,header=T,stringsAsFactors=F,colClasses=c("character","character"))
    selectInput('filterVariantSelectDelete', 'Select filter', 
                choices = data$Name[-1], 
                selected=data$Name[2],
                multiple=T,
                selectize = FALSE)
  })
  
  #Confirm load
  observe({
    input$filterVariantConfirmLoad
    isolate({
      if (length(input$filterVariantSelectLoad)>0) {
        data<-read.table(sessionvalues$variantGroupFile,header=T,stringsAsFactors=F,colClasses=c("character","character"))
        selected<-input$filterVariantSelectLoad
        sqlQuery<-data$SQL[which(data$Name==selected)]
        if (sqlQuery=="") sqlQuery<-"reset"
        session$sendCustomMessage(type='filterVariantCallbackHandlerLoad', sqlQuery)
        toggleModal(session, "filterVariantModalLoad", toggle = "close")
      }
    })
  })
  
  #Confirm Save
  observe({
    input$filterVariantConfirmSave
    isolate({
      if ((length(input$filterVariantNameSave)>0) & length(input$filterVariantQueryBuilderSQL)>0) {
        data<-data.frame(input$filterVariantNameSave,input$filterVariantQueryBuilderSQL)
        colnames(data)<-c("group","sql")
        write.table(file=(sessionvalues$variantGroupFile),data,append=T,row.names=F,col.names=F,sep="\t")
        toggleModal(session, "filterVariantModalSave", toggle = "close")
      }
    })
  })
  
  #Confirm Delete
  observe({
    input$filterVariantConfirmDelete
    isolate({
      
      if (length(input$filterVariantSelectDelete)>0) {
        data<-read.table(sessionvalues$variantGroupFile,header=T,stringsAsFactors=F,colClasses=c("character","character"))
        selected<-input$filterVariantSelectDelete
        data<-data[-match(selected,data$Name),]
        colnames(data)<-c("Name","SQL")
        write.table(file=sessionvalues$variantGroupFile,data,row.names=F,col.names=T,sep="\t")
        toggleModal(session, "filterVariantModalDelete", toggle = "close")
      }
    })
  })
  
  output$filterVariant<-renderUI({
    fluidRow(shiny::column(12,
                    actionButton("filterVariantLoad", label = "Load filter"),  
                    actionButton("filterVariantSave", label = "Save filter"),
                    actionButton("filterVariantDelete", label = "Delete filter(s)"),
                    tags$div(class="extraspace2"),
                    queryBuildROutput("filterVariantQueryBuilderWidget",width="970px",height="100%"),
                    tags$div(class="extraspace2"),
                    actionButton("filterVariantApply", label = "Apply filter"),
                    tags$script('
                                function filterVariantGetSQLStatement() {
                                var sql = $("#filterVariantQueryBuilderWidget").queryBuilder("getSQL", false);
                                Shiny.onInputChange("filterVariantQueryBuilderSQL", sql);
                                };
                                document.getElementById("filterVariantApply").onclick = function() {filterVariantGetSQLStatement()}
                                '),
                    tags$script('            
                                Shiny.addCustomMessageHandler("filterVariantCallbackHandlerLoad",  function(sqlQuery) {
                                if (sqlQuery=="reset") $("#filterVariantQueryBuilderWidget").queryBuilder("reset")
                                else $("#filterVariantQueryBuilderWidget").queryBuilder("setRulesFromSQL",sqlQuery);
                                });
                                '),
                    tags$script('
                  document.getElementById("filterVariantSave").onclick = function() {filterVariantGetSQLStatement()}
                  '),
                    bsModal("filterVariantModalLoad", "Load filter", "filterVariantLoad", 
                            size = "small",
                            uiOutput("filterVariantSelectLoadUI"),
                            actionButton("filterVariantConfirmLoad", label = "Load")
                    ),
                    bsModal("filterVariantModalSave", "Save filter", "filterVariantSave", 
                            size = "small",
                            textInput("filterVariantNameSave", "Save filter as:", value = ""),
                            actionButton("filterVariantConfirmSave", label = "Save")
                    ),
                    bsModal("filterVariantModalDelete", "Delete filter", "filterVariantDelete", 
                            size = "small",
                            uiOutput("filterVariantSelectDeleteUI"),
                            actionButton("filterVariantConfirmDelete", label = "Delete")
                    )
    )
    )
  })
  
  output
}
