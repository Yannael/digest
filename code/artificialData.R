
#s: vector of samples (e.g, 1:100)
#v: number of variants
#g: number of genes
genData<-function(s,v,g) {
  
  v.s<-rep(s,each=(v*g))
  v.chr<-rep(1,v*length(s)*g)
  v.pos<-rep(1:(v*g),length(s))
  v.ref<-rep("A",v*length(s)*g)
  v.alt<-rep("C",v*length(s)*g)
  v.g<-rep(1:g,length(s),each=v)
  v.geno<-base::sample(c(1,2),v*length(s)*g,rep=T)
  
  data<-data.frame(sample_id=v.s,chr=v.chr,pos=v.pos,ref=v.ref,alt=v.alt,gene_symbol=v.g,zygosity=v.geno)
  
  data
  
}

score<-function(v) {
  1
}

dummy<-function() {
  
  i<-"100_100_10000"
  
  system.time(data<-genData(1:100,100,10000))
  system.time(write.table(file=paste0("data",i,".csv"),data,col.names=F,row.names=F,quote=F,sep=","))
  
  object.size(data)
  
  Sys.setenv('SPARKR_SUBMIT_ARGS'='"--packages" "com.databricks:spark-csv_2.10:1.0.3" "sparkr-shell"')
  #Sys.setenv('SPARKR_SUBMIT_ARGS'='"sparkr-shell"')
  
  Sys.setenv(SPARK_HOME="/Users/yalb/spark")
  .libPaths(c(file.path(Sys.getenv("SPARK_HOME"), "R", "lib"), .libPaths()))
  library(SparkR)
  
  sparkEnvir <- list('spark.sql.parquet.binaryAsString'='true') #Needed to read strings from Parquet
  sc<-sparkR.init(master="local[4]",sparkEnvir=sparkEnvir)
  sqlContext <- sparkRSQL.init(sc)
  
  df1 <- read.df(sqlContext, "output", "com.databricks.spark.csv",delimiter=",")
  
  df <- read.df(sqlContext, "variantsulb", "parquet")
  registerTempTable(df, "df")
  
  
  
  library(magrittr)
  df2 <- group_by(df,"chr","position") 
  
  df22<-lapply(df,score)
  
  #df3 <- summarize(df2,count = sum(df$position))
  #%>% summarize(count = n(flights$dest))
  
  
  RDDcase = sql(sqlContext,"SELECT sample_id,chr,position,reference,alternative,gene_symbol,zygosity FROM df ")
  
  
  
  
  
  df <- sql(sqlContext, "SELECT * FROM variants limit 1")
  sql(sqlContext, "describe sqlContext")
  
  
  
}