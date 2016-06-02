# Docker container for DiGeST

This Dockerfile sets up a container for running DiGeST. It is based on rocker/shiny and further installs

* Java 1.8
* Spark 1.6 for Scala 2.10
* Python numpy and scipy
* The following R libraries:
..* From CRAN: 'jsonlite','plyr', 'ggplot2','devtools','htmlwidgets','shinyBS', 'shinyjs''RMySQL','RSQLite','RJDBC','RCurl','roxygen2'
..* From Github: 'rstudio/DT','yannael/queryBuildR','smartinsightsfromdata/rpivotTable'


##Install and run

Clone the repository and build the container

 git clone https://github.com/Yannael/digest
 cd docker
 docker build -t digest .
 cd ..
 
Download example variant data

 wget http://litpc45.ulb.ac.be/variant.tgz
 tar xvzf variant.tgz
 
Run 

 docker run -v `pwd`/code:/srv/shiny-server/digest -v `pwd`/variantsulb:/home/shiny/variants -p 3838:3838 -it digest bash

##Run from dockerhub

Run
 
 docker run yannael/digest
 

 
 