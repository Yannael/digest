library('shiny')
library('shinyBS')
library('shinyjs')

library('RMySQL')
library('RJDBC')

library("jsonlite")
library('httr')
library('RCurl')
library('ggplot2')
library('DT')
library('queryBuildR')
library('rpivotTable')
library('plyr')

VARIANTS<-"/home/shiny/variants"
#VARIANTS<-"/Users/yalb/Projects/BridgeIRIS/digest/exomes_1000g_p_subset"

SPARK_HOME<-"/home/shiny/spark"
#SPARK_HOME<-"/Users/yalb/spark"
Sys.setenv(SPARK_HOME=SPARK_HOME)
Sys.setenv(PATH=paste0(SPARK_HOME,"/bin:",SPARK_HOME,"/sbin:",Sys.getenv("PATH")))

USE_CLUSTER<-FALSE

folderAnalyses<-"users/analyses/"

if (USE_CLUSTER) {
  
  IMPALA_CLASSPATH<-"impala-jdbc-0.5-2"
  IMPALA_SERVER<-"jdbc:hive2://127.0.0.1:21050/;auth=noSasl"
  
  drv <- JDBC(driverClass = "org.apache.hive.jdbc.HiveDriver",
              classPath = list.files(IMPALA_CLASSPATH,pattern="jar$",full.names=T),
              identifier.quote="`")
  
  VARIANTS_TABLE<-'digest.exomes_1000g_p'
  
} else {
  .libPaths(c(file.path(Sys.getenv("SPARK_HOME"), "R", "lib"), .libPaths()))
  library(SparkR)
  
  sparkEnvir <- list('spark.sql.parquet.binaryAsString'='true') #Needed to read strings from Parquet
  sparkR.session(master="local[4]",sparkEnvir=sparkEnvir)
  
  VARIANTS_TABLE<-"variants"
  
  df <- read.df(VARIANTS, "parquet")
  createOrReplaceTempView(df, VARIANTS_TABLE);
  
}

if (!file.exists("/tmp/spark-events")) {
  dir.create("/tmp/spark-events")
}

#updatePhenotypeTable()

#Retrieve phenotype table from local DB
loadPhenotypes<-function(sql) {
  
  if (sql!="") {
    sql<-paste0(" where ",sql)
    sql<-gsub(',',"','",sql)
  }
  condb<-dbConnect(RSQLite::SQLite(), paste0("phenotypesSubset.db"))
  data<-dbGetQuery(condb,paste0("select * from phenotypes ",sql))
  dbDisconnect(condb)
  data
}

######################
#Generic functions
#####################

getWidgetTable<-function(data,session,selection='none') {
  action <- dataTableAjax(session, data,rownames=F)
  widget<-datatable(data, 
                    rownames=F,
                    escape=F,
                    selection = selection,
                    options = list(
                      ajax = list(url = action),
                      dom= 'lipt',
                      lengthMenu = list(c(20, 100, 1000), c('10', '100','1000')),pageLength = 20,
                      columnDefs = list(
                        # list(
                        #   targets = c(1),
                        #   render = JS(
                        #     "function(data, type, row, meta) {",
                        #     "return type === 'display' && data.length > 15 ?",
                        #     "'<span title=\"' + data + '\">' + data.substr(0, 11) + '...</span>' : data;",
                        #     "}")
                        # ),
                        list(className="dt-right",targets="_all")
                      )
                    )
  )
  widget
}

preprocSQL<-function(sql) {
  if (sql!="") {
    sql<-paste0(" where ",sql)
    sql<-gsub(' *, *',"','",sql)
  }
  sql
}

loadData<-function(sql,noLimit=F,maxRows=1000,preproc=T) {
  if (preproc) {
    sql<-preprocSQL(sql)
  }
  if (USE_CLUSTER) {
    condb <- dbConnect(drv,IMPALA_SERVER )
    nbrows<-dbGetQuery(condb,paste0("select count(*) from ",VARIANTS_TABLE," ",sql))
  }
  else {
    nbrows<-collect(sql(paste0("select count(*) from ",VARIANTS_TABLE," ",sql)))
  }
  
  if (noLimit) limit<-""
  else limit<-paste0(" limit ",maxRows)
  nbRowsExceededWarningMessage<-""
  if (nbrows>maxRows) {
    nbRowsExceededWarningMessage<-paste0("Warning: Query returns <b>",nbrows," records</b>. First ",maxRows," retrieved.")
  }
  
  if (USE_CLUSTER) {
    data<-dbGetQuery(condb,paste0("select * from ",VARIANTS_TABLE," ",sql,limit))
    dbDisconnect(condb)
  }
  else {
    data<-collect(sql(paste0("select * from ",VARIANTS_TABLE," ",sql,limit)))
  }
  
  results<-list(data=data,nbRowsExceededWarningMessage=nbRowsExceededWarningMessage)
  results
}

######################
#Process results
#####################

procRes<-function(folderAnalyses,nameAnalysis) {
  json_file<-paste0(folderAnalyses,nameAnalysis,"_metadata.json")
  csv_file<-paste0(folderAnalyses,nameAnalysis,"_ranking.csv")
  
  metadata<-fromJSON(json_file)
  
  res<-list()
  res$name<-metadata[[1]]
  res$scale<-metadata[[2]]
  res$scope<-metadata[[3]]
  res$sqlControl<-metadata[[4]]
  res$sqlCase<-metadata[[5]]
  res$controlSampleID=metadata[[6]]
  res$caseSampleID=metadata[[7]]
  res$controlGroupName=metadata[[8]]
  res$caseGroupName=metadata[[9]]
  res$ntests=metadata[[10]]
  res$start_time<-as.POSIXct(metadata[[11]],origin = "1970-01-01",tz="Europe/Brussels")
  res$end_time<-as.POSIXct(metadata[[12]],origin = "1970-01-01",tz="Europe/Brussels")
  res$run_time<-round(sum(metadata[[13]]),digits=2)
  
  res$scores<-read.table(csv_file,stringsAsFactors=F)
  
  if (res$scale=="variant") {
    if (res$scope=="monogenic") {
      variantID<-res$scores[,1]
      geneID<-res$scores[,2]
      geneID_Link<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",geneID,"' target='_blank'>",geneID,"</a>")
      res$scores[,6]<-format(-log(as.numeric(res$scores[,4]))/log(10),scientific=FALSE,digits=3)
      res$scoreSummaryRaw<-cbind(variantID,geneID,res$scores[,c(-1,-2)])
      colnames(res$scoreSummaryRaw)<-c("Variant_ID","Gene_Symbol","Ratio_Difference","-log10_P_Value","Total_Variants","Score_Case","Score_Control")
      res$scoreSummary<-cbind(variantID,geneID_Link,res$scores[,c(-1,-2)])
      colnames(res$scoreSummary)<-c("Variant_ID","Gene_Symbol","Ratio_Difference","-log10_P_Value","Total_Variants","Score_Case","Score_Control")
      
    }
    
    if (res$scope=="digenic") {
      variantID1<-res$scores[,1]
      geneID1<-res$scores[,2]
      variantID2<-res$scores[,3]
      geneID1_Link<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",geneID1,"' target='_blank'>",geneID1,"</a>")
      geneID2<-res$scores[,4]
      geneID2_Link<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",geneID2,"' target='_blank'>",geneID2,"</a>")
      res$scores[,5]<-format(as.numeric(res$scores[,5]),scientific=FALSE,digits=3)
      res$scores[,6]<-format(-log(as.numeric(res$scores[,6]))/log(10),scientific=FALSE,digits=3)
      res$scoreSummaryRaw<-cbind(variantID1,geneID1,variantID2,geneID2,res$scores[,-c(1:4)])
      colnames(res$scoreSummaryRaw)<-c("Variant_ID1","Gene_Symbol1","Variant_ID2","Gene_Symbol2","Ratio_Difference","log10_P_Value","Total_Variants","Score_Case","Score_Control")
      res$scoreSummary<-cbind(variantID1,geneID1_Link,variantID2,geneID2_Link,res$scores[,-c(1:4)])
      colnames(res$scoreSummary)<-c("Variant_ID1","Gene_Symbol1","Variant_ID2","Gene_Symbol2","Ratio_Difference","log10_P_Value","Total_Variants","Score_Case","Score_Control")
    }
  }
  if (res$scale=="gene") {
    if (res$scope=="monogenic") {
      geneID<-res$scores[,1]
      geneID_Link<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",geneID,"' target='_blank'>",geneID,"</a>")
      res$scores[,2]<-format(as.numeric(res$scores[,2]),scientific=FALSE,digits=3)
      res$scores[,3]<-format(-log(as.numeric(res$scores[,3]))/log(10),scientific=FALSE,digits=3)
      
      res$scoreSummaryRaw<-cbind(geneID,res$scores[,-1])
      colnames(res$scoreSummaryRaw)<-c("Gene_Symbol","Ratio_Difference","log10_P_Value","Total_Variants","Score_Case","Score_Control")
      res$scoreSummary<-cbind(geneID_Link,res$scores[,-1])
      colnames(res$scoreSummary)<-c("Gene_Symbol","Ratio_Difference","log10_P_Value","Total_Variants","Score_Case","Score_Control")
      
    }
    
    if (res$scope=="digenic") {
      geneID1<-res$scores[,1]
      geneID1_Link<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",geneID1,"' target='_blank'>",geneID1,"</a>")
      geneID2<-res$scores[,2]
      geneID2_Link<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",geneID2,"' target='_blank'>",geneID2,"</a>")
      
      res$scores[,3]<-format(as.numeric(res$scores[,3]),scientific=FALSE,digits=3)
      res$scores[,4]<-format(-log(as.numeric(res$scores[,4]))/log(10),scientific=FALSE,digits=3)
      
      res$scoreSummaryRaw<-cbind(geneID1,geneID2,res$scores[,-c(1:2)])
      colnames(res$scoreSummaryRaw)<-c("Gene_Symbol1","Gene_Symbol2","Ratio_Difference","log10_P_Value","Total_Variants","Score_Case","Score_Control")
      res$scoreSummary<-cbind(geneID1_Link,geneID2_Link,res$scores[,-c(1:2)])
      colnames(res$scoreSummary)<-c("Gene_Symbol1","Gene_Symbol2","Ratio_Difference","log10_P_Value","Total_Variants","Score_Case","Score_Control")
    }
  }
  res
}

getAnalysesNames<-function(username="") {
  analysesFiles<-basename(Sys.glob("users/analyses/*.json"))
  analysesNames<-as.vector(unlist(sapply(analysesFiles,strsplit,'_metadata.json')))
  analysesNames
}

analysesNames<-getAnalysesNames()

updateVariantFilter<-function() {
  
  condb <- dbConnect(drv,IMPALA_SERVER )
  rs<-dbSendQuery(condb,paste0("select ",fields_select," from ",VARIANTS_TABLE," limit 1"))
  column_info<-dbColumnInfo(rs)
  
  data<-collect(sql(sqlContext,paste0("select * from ",VARIANTS_TABLE," limit 1")))
  data.type<-sapply(data,class)
  column_info<-cbind(names(data.type),as.vector(data.type))
  
  column_info[,2]<-as.character(column_info[,2])
  column_info[,1]<-as.character(column_info[,1])
  
  column_info<-column_info[c(nrow(column_info),1:(nrow(column_info)-1)),]
  column_info[2,2]<-"factor"
  
  column_info<-rbind(column_info[26,],column_info[1:25,])
  
  for (i in 2:nrow(column_info)) {
    if (column_info[i,2]=="character") {
      print(i)
      nbFields<-collect(sql(sqlContext,paste0("select count(distinct ",column_info[i,1],") from ",VARIANTS_TABLE)))
      if (nbFields<30) column_info[i,2]<-"factor"
    }
  }
  
  filters<-list()
  for (i in 1:nrow(column_info)) {
    filterCol<-
      switch(column_info[i,2],
             character=list(
               id= column_info[i,1],
               label= column_info[i,1],
               type= 'string',
               default_value="",
               operators=list('equal','not_equal','contains', 'in', 'not_in','begins_with', 'ends_with','is_null', 'is_not_null')),
             factor={
               values<-sort(collect(sql(sqlContext,paste0("select distinct ",column_info[i,1]," from ",VARIANTS_TABLE)))[,1])
               if (length(values)==1) values<-c("",values)
               list(
                 id= column_info[i,1],
                 label= column_info[i,1],
                 type= 'string',
                 input='select',
                 values=values,
                 default_value=values[1],
                 operators=list('equal','not_equal','contains', 'in', 'not_in','is_null', 'is_not_null'))
             },
             numeric=list(
               id= column_info[i,1],
               label= column_info[i,1],
               type= 'integer',
               default_value="",
               operators=list('equal','not_equal','less', 'less_or_equal', 'greater','greater_or_equal','between','in', 'not_in','is_null', 'is_not_null')),
             integer=list(
               id= column_info[i,1],
               label= column_info[i,1],
               type= 'integer',
               default_value="",
               operators=list('equal','not_equal','less', 'less_or_equal', 'greater','greater_or_equal','between','in', 'not_in','is_null', 'is_not_null'))
     )
    filters<-c(filters,list(filterCol))
  }
  save(file="filterVariantSpec.Rdata",filters)
}



