#FROM rocker/shiny
FROM quantumobject/docker-shiny

#Install Java
RUN apt-get update 
RUN apt-get install -y openjdk-7-jdk

ENV JAVA_HOME /usr/lib/jvm/java-7-openjdk-amd64/jre
ENV PATH $PATH:$JAVA_HOME/bin

#Necessary ton install later rJava
RUN R CMD javareconf

#Install packages for MySQL, SSL, git and nano, SSH. net-tools and whois required for SSH
RUN apt-get -y install libmysqlclient-dev libxml2-dev libssl-dev git nano openssh-server  net-tools whois 

#Install R packages required for DiGeST
#Basic
RUN R -e "install.packages(c('jsonlite','plyr', 'ggplot2','devtools'), repos='http://cran.rstudio.com/')"
#HTML widgets
RUN R -e "install.packages(c('htmlwidgets','shinyBS', 'shinyjs'), repos='http://cran.rstudio.com/')"
RUN R -e "devtools::install_github('rstudio/DT')"
#Database connectors
RUN R -e "install.packages(c('RMySQL','RSQLite','RJDBC','RCurl'), repos='http://cran.rstudio.com/')"
#Doc, for Spark
RUN R -e "install.packages(c('roxygen2'), repos='http://cran.rstudio.com/')"

# queryBuildR and rpivotTable need to be installed from GitHub
WORKDIR /home/docker
RUN git clone https://github.com/Yannael/queryBuildR.git
RUN R CMD INSTALL queryBuildR
RUN git clone https://github.com/smartinsightsfromdata/rpivotTable.git
RUN R CMD INSTALL rpivotTable

#Set password for shiny user
RUN echo shiny:shiny1239 | /usr/sbin/chpasswd

#Install Spark
USER shiny
WORKDIR /home/shiny
RUN wget http://apache.belnet.be/maven/maven-3/3.3.9/binaries/apache-maven-3.3.9-bin.tar.gz
RUN tar xvzf apache-maven-3.3.9-bin.tar.gz
RUN mv apache-maven-3.3.9 apache-maven
ENV M2_HOME=/home/shiny/apache-maven
ENV PATH=/home/shiny/apache-maven/bin:$PATH
RUN wget http://d3kbcqa49mib13.cloudfront.net/spark-1.6.1-bin-hadoop2.6.tgz
RUN tar xvzf spark-1.6.1-bin-hadoop2.6.tgz
RUN mv spark-1.6.1-bin-hadoop2.6 spark
WORKDIR /home/shiny/spark
#ENV MAVEN_OPTS="-Xmx2g -XX:MaxPermSize=512M -XX:ReservedCodeCacheSize=512m"
#RUN mvn -Pyarn -Phadoop-2.6 -Dhadoop.version=2.6.0 -Dyarn.version=2.6.0 -DskipTests -Psparkr -Phive -Phive-thriftserver clean package
ENV SPARK_HOME /home/shiny/spark
ENV PATH $SPARK_HOME/bin:$SPARK_HOME/sbin:$PATH

#Needed for statistical tests in ranking
USER root
RUN apt-get install -y python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose


#Cluster
RUN mkdir /home/shiny/.ssh
ADD id_rsa /home/shiny/.ssh/id_rsa
ADD known_hosts /home/shiny/.ssh/known_hosts
RUN chown -R shiny:shiny /home/shiny/.ssh
RUN chmod 700 /home/shiny/.ssh
RUN chmod 600 /home/shiny/.ssh/id_rsa
RUN chmod 600 /home/shiny/.ssh/known_hosts
##################################


#Make silent SSH connection to cluster, and start Shiny server
ADD startupscript.sh /usr/bin/startupscript.sh

WORKDIR /home/shiny

#Expose ports for OpenVPN, Impala, SSH
EXPOSE 1194
EXPOSE 21050
EXPOSE 22

