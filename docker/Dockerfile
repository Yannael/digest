FROM rocker/shiny

#Install Java
RUN apt-get update 
RUN apt-get install -y openjdk-8-jdk

#Necessary to install rJava
RUN R CMD javareconf

#Install packages for MySQL, SSL, git and nano, SSH. net-tools and whois required for SSH
RUN apt-get -y install libmysqlclient-dev libxml2-dev libssl-dev git nano openssh-server  net-tools whois man 

#Create guest user. IMPORTANT: Change here UID 1000 to your host UID if you plan to share folders.
RUN echo shiny:shiny | chpasswd
USER shiny

ENV HOME /home/shiny
WORKDIR $HOME

#Install Spark
RUN wget http://litpc45.ulb.ac.be/spark1.6_2.10.tgz
RUN tar xvzf spark1.6_2.10.tgz
ENV SPARK_HOME $HOME/spark

USER root

#Needed for statistical tests in ranking
RUN apt-get install -y python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose

#Install R packages required for DiGeST
RUN R -e "install.packages(c('jsonlite','plyr', 'ggplot2','devtools'), repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('htmlwidgets','shinyBS', 'shinyjs'), repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('RMySQL','RSQLite','RJDBC','RCurl'), repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('roxygen2'), repos='http://cran.rstudio.com/')"
RUN R -e "devtools::install_github('rstudio/DT')"
RUN R -e "devtools::install_github('yannael/queryBuildR')"
RUN R -e "devtools::install_github('smartinsightsfromdata/rpivotTable')"

##################################
#For in-house cluster connection
USER root
RUN mkdir /home/shiny/.ssh
ADD id_rsa /home/shiny/.ssh/id_rsa
ADD known_hosts /home/shiny/.ssh/known_hosts
RUN chown -R shiny:shiny /home/shiny/.ssh
RUN chmod 700 /home/shiny/.ssh
RUN chmod 600 /home/shiny/.ssh/id_rsa
RUN chmod 600 /home/shiny/.ssh/known_hosts

#Make silent SSH connection to cluster, and start Shiny server
ADD startupscript.sh /usr/bin/startupscript.sh
##################################

WORKDIR $HOME

#RUN echo "export PATH=$SPARK_HOME/bin:$SPARK_HOME/sbin:$HOME/anaconda2/bin:$PATH" >>.bashrc
RUN usermod -s /bin/bash shiny

#Expose ports for Impala and SSH
EXPOSE 21050
EXPOSE 22

