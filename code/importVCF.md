# Import VCF - Example using 1000 genome data, chromosome 22

Digest reads variant data from the Parquet format. VCF files can be annotated with Highlander to have annotations from SnpEff for example, and then be converted to parquet using the TSV2Parquet script.

## 1) Data

The file 22.exome.vcf.tgz is a VCF file containing variants from the [1000 genome project data](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/), chromosome 22. 

Uncompress with 

 tar xvzf 22.exome.vcf.tgz

## 1) Highlander annotation

Download Highlander from http://sites.uclouvain.be/highlander/download.html (dbBuilder).

Convert 22.exome.vcf to a Highlander annotated TSV format, for a given sample (e.g., HG00096), with

```
java -jar dbBuilder.jar --tool totsv --vcfpath 22.exome.vcf --project test -b 1 -c settings.xml -S HG00096
```

The script extractTSV.sh can be used to convert variant data from the 20 samples provided as examples in http://bridgeiris.ulb.ac.be/digest.

Highlander annotated TSV are in the folder exome.chr22.samples.

## 2) Parquet conversion

Use the notebook TSV2Parquet.ipynb. The files will be converted to parquet in the  folder.
