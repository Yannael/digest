# Spark DiGeST

The core of DiGeST is the Spark python code in /code/spark. Input variables come from the jobsArgument.conf file, which is filled when starting an analysis. It consists of the following variables:

* analysisName: The name fo the analysis
* scope: gene or variant
* scale: monogenic or digenic
* sqlControl: The SQL query for the control group
* sqlCase: The SQL query for the case group
* group1name: The name for group 1 (control)
* group2name: The name for group 2 (case)
* controlMAF: The MAF threshold for the control group
* caseMAF: The MAF threshold for the case group
* pathVariants: The path where variant data are stored (in parquet files)


