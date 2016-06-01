
# coding: utf-8

# In[15]:

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
controlMAF=content[7]
caseMAF=content[8]
pathVariants=content[9]

nPartitions=8
conf = (SparkConf()
         .setMaster("local["+str(nPartitions)+"]")
#         .setAppName(analysisName)
       )
#sc.stop()
sc = SparkContext(conf=conf)



# In[16]:

#sqlContext = HiveContext(sc) #sqlContext._get_hive_ctx() #HiveContext(sc) 
sqlContext = SQLContext(sc)
sqlContext.sql("SET spark.sql.parquet.binaryAsString=true")

#pathVariants='/user/hive/warehouse/digest.db/exomes_hc_ulb'
parquetFile = sqlContext.read.parquet(pathVariants)
parquetFile.registerTempTable("variantData");


# In[17]:

#Input is vector patient, chr, pos, ref, alt, gene_symbol, zygosity
def createKey_VariantGene(variantData):
    #ID is chr:pos:ref:alt
    ID=variantData[1]+":"+str(variantData[2])+":"+variantData[3]+":"+variantData[4]
    
    #return ID, gene_symbol, patient, zygosity
    zygosity=1
    if variantData[6]=="Homozygous":
    #if variantData[6]==2:
        zygosity=2
    patientsID_dictionnary=patientsID_dictionnary_b.value
    patientIndex=patientsID_dictionnary[variantData[0]]
    return ((ID,variantData[5]),(patientIndex,zygosity))

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


# In[18]:

#Transform sparse data (list of (sample_id,zygozity)) into vector z_i
def vectorize(genotypeDataList):
    genotypeDataList=list(genotypeDataList)
    genotypeVector=[0]*len(patientsID_dictionnary_b.value)
    if len(genotypeDataList)>0:
        for j in range(0,len(genotypeDataList)):
            genotypeVector[genotypeDataList[j][0]]=genotypeDataList[j][1]
        
        sumCase=float(sum([int(x>0) for x in genotypeVector[0:patientsID_split_index_b.value]]))
        sumControl=float(sum([int(x>0) for x in genotypeVector[patientsID_split_index_b.value:len(patientsID_dictionnary_b.value)]]))
    
        ratioCase=sumCase/patientsID_split_index_b.value
        ratioControl=sumControl/(len(patientsID_dictionnary_b.value)-patientsID_split_index_b.value)
        
        if (ratioCase>float(caseMAF_b.value)) or (ratioControl>float(controlMAF_b.value)):
            genotypeVector=[0]*len(patientsID_dictionnary_b.value)
        
    return genotypeVector        


# In[19]:

#Compute burden for variantList
def burden(geneID_variantList):
    (geneID,variantList)=geneID_variantList
    variantList=list(variantList)
    burden=[0]*len(patientsID_dictionnary_b.value)
    
    if len(variantList)>0:
        #Go through list of variants
        for i in range(0,len(variantList)):
            #Get variant ID, and list of sample_index,genotype
            (variantID,genotypeDataList)=variantList[i]
            #if genotypeDataList.__class__==tuple:
            #    genotypeDataList=[genotypeDataList]
            #else:
            #    genotypeDataList=list(genotypeDataList)
            
            #Get genotype vector for current variantID
            genotypeDataVector=vectorize(genotypeDataList)
            #And sum with previous genotype vectors
            burden=[x+y for x,y in zip(burden,genotypeDataVector)]
    
    return (geneID,burden)


# In[20]:

#variantList is [(locusID,[genotype])]
def scoreVariant(ID_genotypeDataList):
    ((variantID,geneID),genotypeDataList)=ID_genotypeDataList
    #genotypeList=list(value_GenotypeList)
    
    patientsID_dictionnary=patientsID_dictionnary_b.value
    patientsID_split_index=patientsID_split_index_b.value
    
    genotypeDataVector=vectorize(genotypeDataList)
    
    #sumCase=float(sum([int(x>0) for x in genotypeDataVector[0:patientsID_split_index]]))
    #sumControl=float(sum([int(x>0) for x in genotypeDataVector[patientsID_split_index:len(patientsID_dictionnary)]]))
    sumCase=float(sum([x for x in genotypeDataVector[0:patientsID_split_index]]))
    sumControl=float(sum([x for x in genotypeDataVector[patientsID_split_index:len(patientsID_dictionnary)]]))
    
    ratioCase=sumCase/patientsID_split_index
    ratioControl=sumControl/(len(patientsID_dictionnary)-patientsID_split_index)
        
    score=ratioCase-ratioControl
    #pvalue=fisher_exact([[sumCase,patientsID_split_index-sumCase],[sumControl,len(patientsID_dictionnary)-patientsID_split_index]],'greater')[1]
    pvalue=ttest_ind(genotypeDataVector[0:patientsID_split_index],genotypeDataVector[patientsID_split_index:len(patientsID_dictionnary)])[1]/2
    
    
    if score>0:
        return (variantID,(score,pvalue,ratioCase,ratioControl,sumCase,sumControl))


# In[21]:

def scoreVariantPair(variantIDpair,value_GenotypeListPair):
    
    genotypeListPair=list(value_GenotypeListPair)
    
    patientsID_dictionnary=patientsID_dictionnary_b.value
    patientsID_split_index=patientsID_split_index_b.value
    
    score=0
    if len(genotypeListPair)==2:
        (variantID,genotypeList1)=genotypeListPair[0]
        (variantID,genotypeList2)=genotypeListPair[1]
        
        variantID1=variantID[0]
        variantID2=variantID[1]
        
        genotypeList1=list(genotypeList1)
        genotypeList2=list(genotypeList2)
        
        genotypeVector1=getGenotypeVector(genotypeList1)
        genotypeVector2=getGenotypeVector(genotypeList2)
        
        genotypeVector=[int(x>0 and y>0) for x,y in zip(genotypeVector1,genotypeVector2)]
        
        sumCase=float(sum([int(x>0) for x in genotypeVector[0:patientsID_split_index]]))
        ratioCase=sumCase/patientsID_split_index
        sumControl=float(sum([int(x>0) for x in genotypeVector[(patientsID_split_index+1):len(patientsID_dictionnary)]]))
        ratioControl=sumControl/(len(patientsID_dictionnary)-patientsID_split_index)
        
        score=ratioCase-ratioControl
        pvalue=fisher_exact([[sumCase,patientsID_split_index-sumCase],[sumControl,len(patientsID_dictionnary)-patientsID_split_index]],'greater')[1]
        
        #if score>0:
        return (variantIDpair,((variantID1,variantID2),score,pvalue,ratioCase,ratioControl,sumCase,sumControl))



# In[30]:

#variantList is [(locusID,[sample_index,genotype])]
def scoreGene(geneID_burden):
    (geneID,burden)=geneID_burden
    
    patientsID_dictionnary=patientsID_dictionnary_b.value
    patientsID_split_index=patientsID_split_index_b.value
    
    #sumCase=float(sum([int(x>0) for x in burden[0:patientsID_split_index]]))
    #sumControl=float(sum([int(x>0) for x in burden[patientsID_split_index:len(patientsID_dictionnary)]]))
    sumCase=float(sum([int(x) for x in burden[0:patientsID_split_index]]))
    sumControl=float(sum([int(x) for x in burden[patientsID_split_index:len(patientsID_dictionnary)]]))
    
    ratioCase=sumCase/patientsID_split_index
    ratioControl=sumControl/(len(patientsID_dictionnary)-patientsID_split_index)
        
    score=ratioCase-ratioControl
    #pvalue=fisher_exact([[sumCase,patientsID_split_index-sumCase],[sumControl,len(patientsID_dictionnary)-patientsID_split_index]],'greater')[1]
    pvalue=ttest_ind(burden[0:patientsID_split_index],burden[patientsID_split_index:len(patientsID_dictionnary)])[1]/2
        
    if score>=0:
        return (geneID,score,pvalue,ratioCase,ratioControl,sumCase,sumControl)


# In[38]:

def scoreGenePair(block_i,block_k,i,k):
    block_k=list(block_k)
    len_i=len(block_i)
    len_k=len(block_k)
    scores=[]

    patientsID_dictionnary=patientsID_dictionnary_b.value
    patientsID_split_index=patientsID_split_index_b.value
    
    start_k=0
    skip_last=0
    if i==k:
        skip_last=1
        if len_i==1:
            len_i=0
    
    if len_i>0 and len_k>0:
        for it_i in range(0,len_i-skip_last):
            if i==k:
                start_k=it_i+1
            for it_k in range(start_k,len_k):
                listLoadBlock_i=block_i[it_i]
                listLoadBlock_k=block_k[it_k]
                #genoSum=[int(x>0 and y>0) for x,y in zip(listLoadBlock_i[1],listLoadBlock_k[1])]
                #sumCase=float(sum([int(x>0) for x in genoSum[0:patientsID_split_index]]))
                #sumControl=float(sum([int(x>0) for x in genoSum[(patientsID_split_index):len(patientsID_dictionnary)]]))
                genoSum=[int(x+y) for x,y in zip(listLoadBlock_i[1],listLoadBlock_k[1])]
                sumCase=float(sum([int(x) for x in genoSum[0:patientsID_split_index]]))
                sumControl=float(sum([int(x) for x in genoSum[(patientsID_split_index):len(patientsID_dictionnary)]]))
        
                ratioCase=sumCase/patientsID_split_index
                ratioControl=sumControl/(len(patientsID_dictionnary)-patientsID_split_index)
        
                score=ratioCase-ratioControl
                #pvalue=fisher_exact([[sumCase,patientsID_split_index-sumCase],[sumControl,len(patientsID_dictionnary)-patientsID_split_index]],'greater')[1]
                pvalue=ttest_ind(genoSum[0:patientsID_split_index],genoSum[patientsID_split_index:len(patientsID_dictionnary)])[1]/2
    
                if score>=0:
                    scores.append((listLoadBlock_i[0],listLoadBlock_k[0],score,pvalue,ratioCase,ratioControl,sumCase,sumControl))
    return scores


# In[51]:

def ranking(sqlCase,sqlControl,scale,scope,p):
    start_time = time.time()

    variants_case = sqlContext.sql("SELECT sample_id,chr,pos,ref,alt,gene_symbol,zygosity FROM variantData "+sqlCase)
    variants_control= sqlContext.sql("SELECT sample_id,chr,pos,ref,alt,gene_symbol,zygosity FROM variantData "+sqlControl)

    variants=variants_control.unionAll(variants_case)
    variants_grouped=variants.map(createKey_VariantGene).groupByKey()

    if scope=='monogenic':
        if scale=='variant':
            ntests=variants_grouped.count()
            finish_load_time=time.time()
            runtime_load=finish_load_time - start_time
            scores=variants_grouped.map(scoreVariant).filter(lambda x:x is not None).takeOrdered(1000, key=lambda (k,(v1,v2,v3,v4,v5,v6)): -v1)
        
        if scale=='gene':
            variants_grouped_by_gene=variants_grouped.map(geneAsKey).groupByKey(p)
            ntests=variants_grouped_by_gene.count()
            finish_load_time=time.time()
            runtime_load=finish_load_time - start_time
            burden_by_gene=variants_grouped_by_gene.map(burden)
            burden_by_gene.count()
            finish_burden_time=time.time()
            runtime_burden=finish_burden_time - finish_load_time
            scores=burden_by_gene.map(scoreGene).filter(lambda x:x is not None).takeOrdered(1000, key=lambda (k,v1,v2,v3,v4,v5,v6): -v1)

    if scope=='digenic':
        if scale=='variant':
            variantsID=variants_grouped.keys().map(getVariantID).collect()
            variants_grouped_by_pairs=variants_grouped.flatMap(lambda (k,v):createPairs(k[0],v,variantsID)).groupByKey()
            scores=variants_grouped_by_pairs.map(lambda (k,v):scoreVariantPair(k,v)).filter(lambda x:x is not None).takeOrdered(1000, key=lambda (k,(variants,v1,v2,v3,v4,v5,v6)): -v1)
            ntests=len(variantsID)*(len(variantsID)+1)/2
   
        if scale=='gene':
            variants_grouped_by_gene=variants_grouped.map(geneAsKey).groupByKey(p)
            ntests=variants_grouped_by_gene.count()
            finish_load_time=time.time()
            runtime_load=finish_load_time - start_time
            burden_by_gene=variants_grouped_by_gene.map(burden)
            burden_by_gene.count()
            finish_burden_time=time.time()
            runtime_burden=finish_burden_time - finish_load_time
            burden_by_gene_with_partitions=burden_by_gene.mapPartitionsWithIndex(lambda splitIndex,v: [(splitIndex,list(v))])
            #burden_by_gene_with_partitions.cache()
            scores=[]
            for i in range(0,p):
                block_i=burden_by_gene_with_partitions.filter(lambda (k,v):k==i).collect()[0][1]
                score=burden_by_gene_with_partitions.filter(lambda (k,v):k>=i).flatMap(lambda (k,v):scoreGenePair(block_i,v,i,k)).takeOrdered(1000, key=lambda (gene1,gene2,v1,v2,v3,v4,v5,v6): -v1)
                scores=scores+score
            scores=sc.parallelize(scores,p).takeOrdered(1000, key=lambda (gene1,gene2,v1,v2,v3,v4,v5,v6): -v1)
            ntests=ntests*(ntests+1)/2
    
    end_time=time.time()
    if scale=='variant':
        runtime_score=end_time - finish_load_time
        all_times=[runtime_load,0,runtime_score]
    if scale=="gene":
        runtime_score=end_time - finish_burden_time
        all_times=[runtime_load,runtime_burden,runtime_score]
        
    return (all_times,scores,ntests)


# In[52]:

p=100
runtimes=[]
start_time=time.time()

patientsID_case = sqlContext.sql("SELECT distinct sample_id FROM variantData "+sqlCase).collect()
patientsID_case = [patients[0] for patients in patientsID_case]

patientsID_control = sqlContext.sql("SELECT distinct sample_id FROM variantData "+sqlControl).collect()
patientsID_control = [patients[0] for patients in patientsID_control]

patientsID=patientsID_case+patientsID_control
patientsID_dictionnary=dict(zip(patientsID,range(len(patientsID))))

patientsID_split_index_b = sc.broadcast(len(patientsID_case))
patientsID_dictionnary_b = sc.broadcast(patientsID_dictionnary)
                     
controlMAF_b=sc.broadcast(controlMAF)
caseMAF_b=sc.broadcast(caseMAF)

end_time=time.time()
(all_times,scores,ntests)=ranking(sqlCase,sqlControl,scale,scope,p)


# In[62]:

metadata=[analysisName,scale,scope,sqlControl,sqlCase,patientsID_control,patientsID_case,controlGroupName,caseGroupName,ntests,start_time,end_time,all_times]

with open(analysisName+'_metadata.json', 'w') as outfile:
    json.dump(metadata, outfile)


# In[58]:

with open(analysisName+'_ranking.csv', 'w') as csvfile:
    csvwriter = csv.writer(csvfile, delimiter=' ',
                            quotechar='"', quoting=csv.QUOTE_MINIMAL)
    csvwriter.writerows(scores)


# In[35]:

sc.stop()


# In[ ]:



