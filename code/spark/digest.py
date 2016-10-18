
# coding: utf-8

# In[2]:

from pyspark import SparkContext, SparkConf
from pyspark.sql import SQLContext
from pyspark.sql import HiveContext
import json
import time
import sys
from  scipy.stats import fisher_exact, ttest_ind
import csv

content = [line.rstrip() for line in open('jobsArguments.conf')]

analysisName=content[0]
scope=content[1]
scale=content[2]
sqlControl=content[3]
sqlCase=content[4]
controlGroupName=content[5]
caseGroupName=content[6]
pathVariants=content[7]

nPartitions=8
conf = (SparkConf()
         .setMaster("local["+str(nPartitions)+"]")
       )
#sc.stop()
sc = SparkContext(conf=conf)



# In[3]:

sqlContext = SQLContext(sc)
sqlContext.sql("SET spark.sql.parquet.binaryAsString=true")

parquetFile = sqlContext.read.parquet(pathVariants)
parquetFile.registerTempTable("variantData");


# In[4]:

parquetFile.take(1)


# In[5]:

#Input is vector patient, chr, pos, ref, alt, gene_symbol, zygosity
def createKey_VariantGene(variantData):
    #ID is chr:pos:ref:alt
    ID=str(variantData[1])+":"+str(variantData[2])+":"+variantData[3]+":"+variantData[4]
    
    #return ID, gene_symbol, patient, zygosity
    zygosity=1
    if variantData[6]=="Homozygous":
    #if variantData[6]==2:
        zygosity=2
    patientsID_dictionnary=patientsID_dictionnary_b.value
    patientIndex=patientsID_dictionnary[variantData[0]]
    
    gene_symbol=variantData[5]
    if gene_symbol is None:
        gene_symbol='NotDefined'
    
    return ((ID,gene_symbol),(patientIndex,zygosity))

#variantGeneEntry: key is (variantID,gene), value is (patientIndex,zygosity)
def geneAsKey(variantGeneEntry):    
    return (variantGeneEntry[0][1],(variantGeneEntry[0][0],variantGeneEntry[1]))

def getVariantID(key_VariantGene):
    return key_VariantGene[0]

def f(splitIndex ,v): 
    return [(splitIndex,list(v))]

def toSQLString(strList):
    strList="','".join(strList)
    return "('"+strList+"')"


# In[6]:

def scoreVariant(ID_genotypeDataList):
    ((variantID,geneID),genotypeDataList)=ID_genotypeDataList
    
    n_cases=global_n_cases.value
    n_controls=global_n_controls.value
    
    genotypeDataList=list(genotypeDataList)
    
    sumCase=sum([y for (x,y) in genotypeDataList if x<n_cases])
    sumControl=sum([y for (x,y) in genotypeDataList if x>=n_cases])
    
    score=sumCase-sumControl
    
    pvalue=fisher_exact([[sumCase,sumControl],[n_cases*2-sumCase,n_controls*2-sumControl]])[1]
    
    score=((variantID,geneID),score,pvalue,1,sumCase,sumControl)

    return score


# In[7]:

def scoreGene(geneID_variantList):
    (geneID,variantList)=geneID_variantList
    
    n_cases=global_n_cases.value
    n_controls=global_n_controls.value
    
    sumCase=0
    sumControl=0
    
    variantList=list(variantList)
    n_variants=len(variantList)
    
    for i in range(n_variants):
        (variantID,genotypeDataList)=variantList[i]
        genotypeDataList=list(genotypeDataList)
    
        sumCase=sumCase+sum([y for (x,y) in genotypeDataList if x<n_cases])
        sumControl=sumControl+sum([y for (x,y) in genotypeDataList if x>=n_cases])

        
    score=sumCase-sumControl
    
    pvalue=fisher_exact([[sumCase,sumControl],[n_cases*n_variants*2-sumCase,n_controls*n_variants*2-sumControl]])[1]
    
    score=((geneID),score,pvalue,n_variants,sumCase,sumControl)

    return score


# In[8]:

def scorePartition(partition_i):
    partition_i=list(partition_i)
    len_i=len(partition_i)
    
    scores=[]
    
    for it in range(len_i):
        if global_scale.value=="variant":
            scores.append((scoreVariant(partition_i[it])))
        if global_scale.value=="gene":
            scores.append((scoreGene(partition_i[it])))
    
    return scores


# In[9]:

def scoreVariantPair(ID_genotypeDataList1,ID_genotypeDataList2):
    ((variantID1,geneID1),genotypeDataList1)=ID_genotypeDataList1
    ((variantID2,geneID2),genotypeDataList2)=ID_genotypeDataList2
    
    n_cases=global_n_cases.value
    n_controls=global_n_controls.value
    
    genotypeDataList1=list(genotypeDataList1)
    genotypeDataList2=list(genotypeDataList2)
    
    sumCase=sum([y for (x,y) in genotypeDataList1 if x<n_cases])
    sumCase=sumCase+sum([y for (x,y) in genotypeDataList2 if x<n_cases])
    sumControl=sum([y for (x,y) in genotypeDataList1 if x>=n_cases])
    sumControl=sumControl+sum([y for (x,y) in genotypeDataList2 if x>=n_cases])
    
    score=sumCase-sumControl
    
    pvalue=fisher_exact([[sumCase,sumControl],[n_cases*2*2-sumCase,n_controls*2*2-sumControl]])[1]
   
    
    score=((variantID1,geneID1,variantID2,geneID2),score,pvalue,2,sumCase,sumControl)
    
    return score


# In[10]:

def scoreGenePair(geneID_scoreList1,geneID_scoreList2):
    (geneID1,score1,pvalue1,n_variants1,sumCase1,sumControl1)=geneID_scoreList1
    (geneID2,score2,pvalue2,n_variants2,sumCase2,sumControl2)=geneID_scoreList2
    
    n_cases=global_n_cases.value
    n_controls=global_n_cases.value
    
    sumCase=0
    sumControl=0
    
    sumCase=sumCase1+sumCase2
    sumControl=sumControl1+sumControl2

    score=score1+score2

    pvalue=fisher_exact([[sumCase,sumControl],[n_cases*(n_variants1+n_variants2)*2-sumCase,n_controls*(n_variants1+n_variants2)*2-sumControl]])[1]

    
    score=((geneID1,geneID2),score,pvalue,n_variants1+n_variants2,sumCase,sumControl)
    
    return score


# In[11]:

def scorePartitionPair(keyPair,partitionPair):
    (i,k)=keyPair
    partitionPair=list(partitionPair)
    partition_i=partitionPair[0]
    if i==k:
        partition_k=partition_i
    else:
        partition_k=partitionPair[1]
    
    len_i=len(partition_i)
    len_k=len(partition_k)
    scores=[]

    start_k=0
    skip_last=0
    
    #Special case: If intra-block computation,
    if i==k:
        skip_last=1
        #If block has size one, no interaction to compute
        if len_i==1:
            len_i=0
    
    if len_i>0 and len_k>0:
        for it_i in range(0,len_i-skip_last):
            if i==k:
                start_k=it_i+1
            for it_k in range(start_k,len_k):
                if global_scale.value=="variant":
                    score=scoreVariantPair(partition_i[it_i],partition_k[it_k])
                if global_scale.value=="gene":
                    score=scoreGenePair(partition_i[it_i],partition_k[it_k])
                scores.append((score))

    return scores


# In[12]:

def pairBlocks(k,v,list_blocksPairs):
    listPartitionPairs=[]
    for i in range(0,len(list_blocksPairs)):
            if k in list_blocksPairs[i]:
                listPartitionPairs.append((list_blocksPairs[i],v))
                
    return listPartitionPairs


# In[13]:

def ranking(scale,scope,p,genotypeMatrixRDD_indexed):

    if scope=='monogenic':
        scores=genotypeMatrixRDD_indexed.flatMap(lambda (k,v):scorePartition(v)).filter(lambda x:x is not None).takeOrdered(1000, key=lambda (k,v1,v2,v3,v4,v5): -v1)

    if scope=='digenic':
        
        if scale=='gene':
            genotypeMatrixRDD_indexed=genotypeMatrixRDD_indexed.map(lambda (k,v):(k,scorePartition(v)))
        
        list_blocksPairs=[]
        for i in range(0,p):
            for j in range(i,p):
                list_blocksPairs.append((i,j))
            
        pairedGenotypeMatrixRDD=genotypeMatrixRDD_indexed.flatMap(lambda (k,v):pairBlocks(k,v,list_blocksPairs)).groupByKey()
        scores=pairedGenotypeMatrixRDD.flatMap(lambda (k,v):scorePartitionPair(k,v)).takeOrdered(1000, key=lambda (k,v1,v2,v3,v4,v5): -v1)
            
    return scores


# In[14]:

def getGenotypeMatrixRDD(sqlCase,sqlControl,p,scale):
    
    variants_case = sqlContext.sql("SELECT sample_id,chr,pos,ref,alt,gene_symbol,zygosity FROM variantData "+sqlCase)
    variants_control= sqlContext.sql("SELECT sample_id,chr,pos,ref,alt,gene_symbol,zygosity FROM variantData "+sqlControl)

    variants=variants_control.unionAll(variants_case)
    
    genotypeMatrixRDD=variants.rdd.map(createKey_VariantGene).groupByKey(p)
    
    if scale=='gene':
        genotypeMatrixRDD=genotypeMatrixRDD.map(geneAsKey).groupByKey(p)

    return genotypeMatrixRDD


# In[15]:

p=10

runtimes=[]
start_time=time.time()

patientsID_case = sqlContext.sql("SELECT distinct sample_id FROM variantData "+sqlCase).collect()
patientsID_case = [patients[0] for patients in patientsID_case]

patientsID_control = sqlContext.sql("SELECT distinct sample_id FROM variantData "+sqlControl).collect()
patientsID_control = [patients[0] for patients in patientsID_control]

patientsID=patientsID_case+patientsID_control
patientsID_dictionnary=dict(zip(patientsID,range(len(patientsID))))

patientsID_dictionnary_b = sc.broadcast(patientsID_dictionnary)

global_n_cases = sc.broadcast(len(patientsID_case))
global_n_controls = sc.broadcast(len(patientsID_control))

global_scale = sc.broadcast(scale)

genotypeMatrixRDD=getGenotypeMatrixRDD(sqlCase,sqlControl,p,scale)
genotypeMatrixRDD_indexed=genotypeMatrixRDD.mapPartitionsWithIndex(lambda splitIndex,v: [(splitIndex,list(v))])

scores=ranking(scale,scope,p,genotypeMatrixRDD_indexed)

end_time=time.time()
all_times=end_time-start_time

ntests=0



# In[16]:

scores


# In[17]:

metadata=[analysisName,scale,scope,sqlControl,sqlCase,patientsID_control,patientsID_case,controlGroupName,caseGroupName,ntests,start_time,end_time,all_times]

with open(analysisName+'_metadata.json', 'w') as outfile:
    json.dump(metadata, outfile)


# In[18]:

if scale=="gene" and scope=="monogenic":
    scores2=scores
else:
    scores2=[list(k1)+[k2,k3,k4,k5,k6] for (k1,k2,k3,k4,k5,k6) in scores]
with open(analysisName+'_ranking.csv', 'w') as csvfile:
    csvwriter = csv.writer(csvfile, delimiter=' ',
                            quotechar='"', quoting=csv.QUOTE_MINIMAL)
    csvwriter.writerows(scores2)


# In[35]:

sc.stop()

