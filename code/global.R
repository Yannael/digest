
require('RMySQL')
require("jsonlite")
require('rpivotTable')
require('plyr')


require('shiny')
require('RJDBC')
require('RCurl')
require('queryBuildR')

#VARIANTS<-"/home/shiny/variantsulb"
VARIANTS<-"/Users/yalb/Projects/Github/variants/variantsulb"

#SPARK_HOME<-"/home/shiny/spark"
SPARK_HOME<-"/Users/yalb/spark"
Sys.setenv(SPARK_HOME=SPARK_HOME)
Sys.setenv(PATH=paste0(SPARK_HOME,"/bin:",SPARK_HOME,"/sbin:",Sys.getenv("PATH")))

USE_CLUSTER<-TRUE

folderAnalyses<-"users/analyses/"

if (USE_CLUSTER) {
  
  IMPALA_CLASSPATH<-"impala-jdbc-0.5-2"
  IMPALA_SERVER<-"jdbc:hive2://127.0.0.1:21050/;auth=noSasl"
  
  drv <- JDBC(driverClass = "org.apache.hive.jdbc.HiveDriver",
              classPath = list.files(IMPALA_CLASSPATH,pattern="jar$",full.names=T),
              identifier.quote="`")
  
  VARIANTS_TABLE<-'digest.exomes_hc_ulb'
  #VARIANTS_TABLE<-'digest.exomes_hc_1000g'
  
} else {
  .libPaths(c(file.path(Sys.getenv("SPARK_HOME"), "R", "lib"), .libPaths()))
  library(SparkR)
  require('shiny') #Necessary for column function
  
  sparkEnvir <- list('spark.sql.parquet.binaryAsString'='true') #Needed to read strings from Parquet
  sc<-sparkR.init(master="local[4]",sparkEnvir=sparkEnvir)
  sqlContext <- sparkRSQL.init(sc) #sparkRHive.init(sc)#
  
  VARIANTS_TABLE<-"variants"
  
  df <- read.df(sqlContext, VARIANTS, "parquet")
  registerTempTable(df, VARIANTS_TABLE);
  
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
  condb<-dbConnect(RSQLite::SQLite(), paste0("phenotypes.db"))
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
                        list(
                          targets = c(1),
                          render = JS(
                            "function(data, type, row, meta) {",
                            "return type === 'display' && data.length > 15 ?",
                            "'<span title=\"' + data + '\">' + data.substr(0, 11) + '...</span>' : data;",
                            "}")
                        ),
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
    nbrows<-collect(sql(sqlContext,paste0("select count(*) from ",VARIANTS_TABLE," ",sql)))
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
    data<-collect(sql(sqlContext,paste0("select * from ",VARIANTS_TABLE," ",sql,limit)))
  }
  
  results<-list(data=data,nbRowsExceededWarningMessage=nbRowsExceededWarningMessage)
  results
}

######################
#Process results
#####################


get4<-function(l) {l[[4]]}
get1<-function(l) {l[[1]]}
get2<-function(l) {l[[2]]}

procRes<-function(results) {
  res<-list()
  res$name<-results[[1]]
  res$scale<-results[[2]]
  res$scope<-results[[3]]
  res$sqlCase<-results[[4]]
  res$sqlControl<-results[[5]]
  res$start_time<-as.POSIXct(results[[4]],origin = "1970-01-01",tz="Europe/Brussels")
  res$end_time<-as.POSIXct(results[[5]],origin = "1970-01-01",tz="Europe/Brussels")
  res$run_time<-round(results[[6]],digits=2)
  res$scores<-sapply(results[[7]],get2)
  res$locus<-sapply(results[[7]],get1)
  res$caseSampleID=results[[8]]
  res$controlSampleID=results[[9]]
  res$group1name=results[[10]]
  res$group2name=results[[11]]
  
  if (res$scale=="variant") {
    if (res$scope=="monogenic") {
      variantsID<-res$locus[1,]
      scoreSummary<-cbind(t(data.frame(sapply(variantsID,strsplit,':'))),res$locus[2,])
      colnames(scoreSummary)<-c("Chr","Position","Reference","Alternative",'Gene_Symbol')
      rownames(scoreSummary)<-NULL
      scoreSummary[,'Gene_Symbol']<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",scoreSummary[,'Gene_Symbol'],"' target='_blank'>",scoreSummary[,'Gene_Symbol'],"</a>")
      #if (results[[2]]=="pairVariantsMonogenic") {
      #  uniqueid2<-apply(res$locus[5:7,],2,paste,collapse=":")
      #  to.keep<-match(uniqueid2,dbvariants$uniqueid)
      #  res$infovariants2<-dbvariants[to.keep,]
      #  
      #  scoreSummary2<-cbind(res$locus[5,],res$infovariants2[,'reference'],res$infovariants2[,'alternative'],res$infovariants1[,'gene_symbol'])
      #  colnames(scoreSummary2)<-c("Locus2","Reference2", "Alternative2","Gene symbol")
      #}
    }
    
    if (results[[2]]=="pairVariantsDigenic") {
      uniqueid2<-apply(res$locus[6:8,],2,paste,collapse=":")
      to.keep<-match(uniqueid2,dbvariants$uniqueid)
      res$infovariants2<-dbvariants[to.keep,]
      scoreSummary2<-cbind(res$infovariants1[,'gene_symbol'],res$locus[6,],res$infovariants2[,'reference'],res$infovariants2[,'alternative'],res$infovariants2[,'gene_symbol'])
      colnames(scoreSummary2)<-c("Gene symbol1","Locus2","Reference2", "Alternative2","Gene symbol2")
      
    }
    
    res$scoreSummary<-cbind(Score=t(res$scores),scoreSummary)
    colnames(res$scoreSummary)[1:3]<-c("Score","Score_Case","Score_Control")
    
  }
  
  if (res$scale=="gene") {
    if (res$scope=="monogenic") {
      #browser()
      geneID<-res$locus
      geneID_Link<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",geneID,"' target='_blank'>",geneID,"</a>")
      res$scores[1,]<-format(as.numeric(res$scores[1,]),scientific=FALSE,digits=3)
      res$scores[2,]<-format(as.numeric(res$scores[2,]),scientific=TRUE,digits=3)
      res$scores[3,]<-format(as.numeric(res$scores[3,]),scientific=FALSE,digits=3)
      res$scores[4,]<-format(as.numeric(res$scores[4,]),scientific=FALSE,digits=3)
      res$scoreSummaryRaw<-cbind(t(res$scores),geneID)
      colnames(res$scoreSummaryRaw)<-c("Ratio_Difference","P_Value","Ratio_Case","Ratio_Control","Score_Case","Score_Control","Gene_Symbol")
      res$scoreSummary<-cbind(t(res$scores),geneID_Link)
      colnames(res$scoreSummary)<-c("Ratio_Difference","P_Value","Ratio_Case","Ratio_Control","Score_Case","Score_Control","Gene_Symbol")
    }
    
    if (res$scope=="digenic") {
      genes<-ldply(res$scores[1,])
      res$scores<-data.matrix(data.frame(res$scores[2:7,]))
      res$scores[1,]<-format(as.numeric(res$scores[1,]),scientific=FALSE,digits=3)
      res$scores[2,]<-format(as.numeric(res$scores[2,]),scientific=TRUE,digits=3)
      res$scores[3,]<-format(as.numeric(res$scores[3,]),scientific=FALSE,digits=3)
      res$scores[4,]<-format(as.numeric(res$scores[4,]),scientific=FALSE,digits=3)
      
      geneID1<-genes[,1]
      geneID1_Link<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",geneID1,"' target='_blank'>",geneID1,"</a>")
      geneID2<-genes[,2]
      geneID2_Link<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",geneID2,"' target='_blank'>",geneID2,"</a>")
      res$scoreSummaryRaw<-cbind(t(res$scores),geneID1,geneID2)
      colnames(res$scoreSummaryRaw)<-c("Ratio_Difference","P_Value","Ratio_Case","Ratio_Control","Score_Case","Score_Control","Gene_Symbol1","Gene_Symbol2")
      res$scoreSummary<-cbind(t(res$scores),geneID1_Link,geneID2_Link)
      colnames(res$scoreSummary)<-c("Ratio_Difference","P_Value","Ratio_Case","Ratio_Control","Score_Case","Score_Control","Gene_Symbol1","Gene_Symbol2")
    }
    #browser()
  }
  res
}

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
  
  if (res$scale=="gene") {
    if (res$scope=="monogenic") {
      #browser()
      geneID<-res$scores[,1]
      geneID_Link<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",geneID,"' target='_blank'>",geneID,"</a>")
      res$scores[,2]<-format(as.numeric(res$scores[,2]),scientific=FALSE,digits=3)
      res$scores[,3]<-format(as.numeric(res$scores[,3]),scientific=TRUE,digits=3)
      res$scores[,4]<-format(as.numeric(res$scores[,4]),scientific=FALSE,digits=3)
      res$scores[,5]<-format(as.numeric(res$scores[,5]),scientific=FALSE,digits=3)
      
      res$scoreSummaryRaw<-cbind(geneID,res$scores[,-1])
      colnames(res$scoreSummaryRaw)<-c("Gene_Symbol","Ratio_Difference","P_Value","Ratio_Case","Ratio_Control","Score_Case","Score_Control")
      res$scoreSummary<-cbind(geneID_Link,res$scores[,-1])
      colnames(res$scoreSummary)<-c("Gene_Symbol","Ratio_Difference","P_Value","Ratio_Case","Ratio_Control","Score_Case","Score_Control")
      
    }
    
    if (res$scope=="digenic") {
      geneID1<-res$scores[,1]
      geneID1_Link<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",geneID1,"' target='_blank'>",geneID1,"</a>")
      geneID2<-res$scores[,2]
      geneID2_Link<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",geneID2,"' target='_blank'>",geneID2,"</a>")
      
      res$scores[,3]<-format(as.numeric(res$scores[,3]),scientific=FALSE,digits=3)
      res$scores[,4]<-format(as.numeric(res$scores[,4]),scientific=TRUE,digits=3)
      res$scores[,5]<-format(as.numeric(res$scores[,5]),scientific=FALSE,digits=3)
      res$scores[,6]<-format(as.numeric(res$scores[,6]),scientific=FALSE,digits=3)
      
      res$scoreSummaryRaw<-cbind(geneID1,geneID2,res$scores)
      colnames(res$scoreSummaryRaw)<-c("Gene_Symbol1","Gene_Symbol2","Ratio_Difference","P_Value","Ratio_Case","Ratio_Control","Score_Case","Score_Control")
      res$scoreSummary<-cbind(geneID1_Link,geneID2_Link,res$scores)
      colnames(res$scoreSummary)<-c("Gene_Symbol1","Gene_Symbol2","Ratio_Difference","P_Value","Ratio_Case","Ratio_Control","Score_Case","Score_Control")
    }
  }
  res
}

getAnalysesNames<-function() {
analysesFiles<-basename(Sys.glob("users/analyses/*.json"))
analysesNames<-as.vector(unlist(sapply(analysesFiles,strsplit,'_metadata.json')))
analysesNames
}

analysesNames<-getAnalysesNames()

dummy<-function() {
  
  condb<-dbConnect(RSQLite::SQLite(), paste0("digest/phenotypes.db"))
  data<-dbGetQuery(condb,paste0("select * from phenotypes "))
  dbDisconnect(condb)
  
  data<-data[which(data[,1]=="1000 Genomes"),]
  dataEUR<-data[(which(data[,5]=="EUR"))[1:25],]
  dataEAS<-data[(which(data[,5]=="EAS"))[1:25],]
  
  phenotypes<-rbind(dataEUR,dataEAS)
  
  condb<-dbConnect(RSQLite::SQLite(), paste0("digest/phenotypesdemo.db"))
  dbWriteTable(condb,"phenotypes",phenotypes,overwrite=T,row.names=F)
  dbDisconnect(condb)
  
  data[which(data[,1]=="ULB"),5]<-""
  write.table(file="phenotypes.csv",data,col.names=T,row.names=F,sep=',')
  
  gene_list<-read.table("Guillaume_genes_new.txt",stringsAsFactor=F)[,1]
  gene_list_str<-paste0(gene_list,collapse=',')
  
  gene_list<-read.table("DIDA_genes.txt",stringsAsFactor=F)[,1]
  gene_list_str<-paste0(gene_list,collapse="','")
  
  #ARHGAP11B,ASPM,ATR,ATRIP,BLM,BRAT1,C7orf27,BUB1B,CASC5,CASK,CCDC7,CDC6,CDK5RAP2,CDT1,CENPF,CENPJ,CEP135,CEP152,CEP250,CEP63,CIT,COX7B,DYRK1A,EFTUD2,EIF2AK3,ERCC3,ERCC4,ERCC5,ERCC6,ERCC8,IER3IP1,KIF11,KMT2B,MLL4,MLL2,LIG4,MCPH1,MYCN,NBN,NDE1,NIN,NIPBL,ORC1,ORC1L,ORC4,ORC4L,ORC6,ORC6L,PCNT,PHC1,PLK4,PNKP,PPM1D,RAD50,RBBP8,RNU4ATAC,SASS6,SLC25A19,SLC9A6,SMC1A,SMC3,STAMBP,STIL,TRMT10A,RG9MTD2,TUBA1A,TUBB,TUBB2B,TUBB3,TUBG1,TUBGCP4,TUBGCP6,UBE3A,WDR62,ZEB2
  
  
}



