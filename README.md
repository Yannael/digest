# DiGeST :Distributed Gene/Variant Scoring Tool

DiGeST is a gene/variant scoring tool, that takes as input two populations (case/control) of SNPs, and outputs a ranking of gene/variant according to their suspected pathogeneticity. It is composed of a front-end which allows a user to browse, filter and create groups of case and control SNPs, and a back-end that performs scoring and ranking in a distributed way, using the Map/Reduce programming model. Scalable storage and computing rely on Parquet file format over HDFS, and the Spark distributed computing framework, respectively. 


Documentation:

* [Installation](../master/doc/install.md)
* [Configuration](../master/doc/configuration.md)
* [Files](../master/doc/files.md)
* [API](../master/doc/API.md)
* [Spark](../master/doc/spark.md)




 