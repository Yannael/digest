# DiGeST : Configuration

Configuration parameters are in the global.R file

## Connection to CliniPhenome

Set CliniPhenomeAPI, e.g.,
```
CliniPhenomeAPI<-"http://bridgeiris.ulb.ac.be:81/bridgeirisportal/search.php"
```

## Cluster or Standalone mode

Set USE_CLUSTER to TRUE for cluster use, false for standalone (local) use, e.g.,
```
USE_CLUSTER<-TRUE
```

###Cluster mode

Set location of Impala JAR files, and Impala server, with IMPALA_CLASSPATH and IMPALA_SERVER, e.g.,
```
IMPALA_CLASSPATH<-"impala-jdbc-0.5-2"
IMPALA_SERVER<-"jdbc:hive2://127.0.0.1:21050/;auth=noSasl"
```

The variant table is set with VARIANTS_TABLE, e.g.,
```
VARIANTS_TABLE<-'highlander.exomes_hc'
```

###Standalone mode

Set SPARK_HOME to your Spark installation folder, e.g.,
```
SPARK_HOME<-"/Users/yalb/spark"
```

Set VARIANTS to the folder containing the SNPs parquet files, e.g.,
```
VARIANTS<-"/Users/yalb/Projects/Github/digest/variantsulb"
```

