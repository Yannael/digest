dummy<-function() {

data<-read.table("genoMatTrios.txt")
#40273 variants
nVarPerGene<-tapply(data[,1],data[,1],length)
length(nVarPerGene)
#15329 genes

dd<-data.frame(nVarPerGene=nVarPerGene)
ggplot(data=dd,aes(nVarPerGene)) + geom_histogram()

#De novo
time<-c(0.7,1.4,2.2,2.5,4.6,8.8)
nsamples<-as.factor(c(30,60,120,30,60,120))
cond<-as.factor(c("200000","200000","200000","880000","880000","880000"))
DF<-data.frame(time=time,nsamples=nsamples,cond=cond)

ggplot(data=DF, aes(x=nsamples, y=time, fill=cond)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black")+
xlab("Number of samples") +
  ylab("Computation time") +
  scale_fill_hue(name="Number of variants")

#Univariate gene
#1 core

times<-c(11,76,116,73,262,420,125,403,824)
nvar<-as.factor(rep(c(100000,500000,1000000),3))
nsamples<-as.factor(rep(c(100,500,1000),each=3))
DF<-data.frame(times=times,nsamples=nsamples,nvar=nvar)

ggplot(data=DF, aes(x=nvar, y=times, fill=nsamples)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black")+
  xlab("Number of variants") +
  ylab("Computation time (s)") +
  scale_fill_hue(name="Number of samples")+
  scale_y_continuous(limits = c(0, 1000))
 # theme(legend.position="bottom")

times<-c(2,8,13,10,34,61,16,57,113)
nvar<-as.factor(rep(c(100000,500000,1000000),3))
nsamples<-as.factor(rep(c(100,500,1000),each=3))
DF<-data.frame(times=times,nsamples=nsamples,nvar=nvar)

ggplot(data=DF, aes(x=nvar, y=times, fill=nsamples)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black")+
  xlab("Number of variants") +
  ylab("Computation time (s)") +
  scale_fill_hue(name="Number of samples")+
  scale_y_continuous(limits = c(0, 1000))

times<-c(11.5,57,102,67,230,472,102,515,1020)
times<-c(6,16,23,13,36,54,24,59,114)
nvar<-as.factor(rep(c(1000,5000,10000),3))
nsamples<-as.factor(rep(c(100,500,1000),each=3))
DF<-data.frame(times=times,nsamples=nsamples,nvar=nvar)

ggplot(data=DF, aes(x=nvar, y=times, fill=nsamples)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black")+
  xlab("Number of genes") +
  ylab("Computation time (s)") +
  scale_fill_hue(name="Number of samples")+
  scale_y_continuous(limits = c(0, 1050))

times<-c(9,16,28,46,24,49,96,168,31,77,142,296)
times<-c(14,35,72,132,66,168,300,600,108,324,660,1140)
nvar<-as.factor(rep(c(100,200,300,400),3))
nsamples<-as.factor(rep(c(100,500,1000),each=4))
DF<-data.frame(times=times,nsamples=nsamples,nvar=nvar)

ggplot(data=DF, aes(x=nvar, y=times, fill=nsamples)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black")+
  xlab("Number of genes") +
  ylab("Computation time (s)") +
  scale_fill_hue(name="Number of samples")+
  scale_y_continuous(limits = c(0, 1150))


#Gene list
genes<-read.table("ID_genelist.txt",stringsAsFactors=F)[,1]
genelist<-paste0("(",paste0(genes,collapse=","),")")

}

dummy<-function() {
  
  #Speedup for monogenic ranking gene, variant
  
  cores<-c(1,5,10,25,50,1,5,10,25,50)
  timesG<-13.3/c(13.3,3.2,1.8,0.75,0.45)
  timesV<-13.6/c(13.6,3.5,1.9,0.9,0.5)
  times<-c(timesG,timesV)
  scale<-c("Gene","Gene","Gene","Gene","Gene","Variant","Variant","Variant","Variant","Variant")
  
  DF<-data.frame(times=times,cores=cores,scale=scale)
  
  ggplot(data=DF, aes(x=cores, y=times,group=scale,colour=scale)) +
    geom_line()+
    xlab("Number of cores") +
    ylab("Speedup") 
  #scale_fill_hue(name="Number of samples")+
    #scale_y_continuous(limits = c(0, 1150))
  
  
}

geneBenchmark<-function() {
  
  #Speedup for monogenic ranking gene, variant
  
  genes<-c(1000,5000,10000,1000,5000,10000,1000,5000,10000)
  times1<-c(13.3,3.2,1.8,0.75,0.45)
  times5<-c(13.3,3.2,1.8,0.75,0.45)
  times50<-c(13.3,3.2,1.8,0.75,0.45)
  times<-c(times1,times5,times50)
  scale<-c("1","1","1","5","5","5","50","50","50")
  
  DF<-data.frame(times=times,cores=cores,scale=scale)
  
  ggplot(data=DF, aes(x=cores, y=times,group=scale,colour=scale)) +
    geom_line()+
    xlab("Number of cores") +
    ylab("Speedup") 
  #scale_fill_hue(name="Number of samples")+
  #scale_y_continuous(limits = c(0, 1150))
  
  
}

plotNbVariantsPerSample_Gene<-function() {
  
  
  #Query from Impala:
  #select count(pos),sample_id from digest.exomes_1000g group by sample_id 
  #Then download as CSV. Result in nbVariantsperSample in BMC folder.
  
  data<-read.csv("~/Projects/Publications/Journals/2016bmc_template/nbVariantsperSample.csv")
  colnames(data)<-c("nb","id")
  
  plot1<-ggplot(data=data, aes(x=nb)) +
    geom_histogram()+
    xlab("Number of variants") +
    ylab("# Samples") 

  #Query from Impala:
  #select count(distinct chr,pos,ref,alt),gene_symbol from digest.exomes_1000g group by gene_symbol
  #Then download as CSV. Result in nbVariantsperGene in BMC folder.
  
  data<-read.csv("~/Projects/Publications/Journals/2016bmc_template/nbVariantsperGene.csv")
  data<-data[which(data[,1]<1000),]
  #data[,1]<-log(data[,1])
  
  colnames(data)<-c("nb","id")
  
  plot2<-ggplot(data=data, aes(x=nb)) +
    geom_histogram()+
    xlab("Number of variants") +
    ylab("# Genes") 
  
}

createSpec<-function(numCol) {
  
  colnamesStr<-paste("Gene VARCHAR, Variant VARCHAR,",paste(paste('V',1:numCol," INT",sep=""),collapse=","),collapse="")
  
  gene count
  
  sqlContext = SQLContext(sc)
  sqlContext.sql("SET spark.sql.parquet.binaryAsString=true")
  
  #parquetFile = sqlContext.read.parquet("hdfs://node001:8020/user/hive/warehouse/highlander.db/exomes_hc")
  parquetFile = sqlContext.read.parquet("hdfs://node001:8020/user/hive/warehouse/1000g.db/exomes_1000g")
  parquetFile.registerTempTable("variantData");
  
  variants = sqlContext.sql("SELECT chr,pos,reference,alternative,gene_symbol,count(patient) FROM variantData group by chr,pos,reference,alternative,gene_symbol")
  variants = variants.groupBy("gene_symbol")
  #vtable=variants.count().collect()
  
  variants.count().write.save("file:///home/yleborgn/1000g", format='com.databricks.spark.csv')
  
  data2<-read.table("geneCount1000g.csv",header=F,sep=",")
  
  
}

