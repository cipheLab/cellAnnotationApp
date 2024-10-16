FROM rocker/rstudio:4.4.1
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libxt-dev \
    python3-pip

# install R packages 

RUN R -e 'install.packages("BiocManager")'
RUN R -e 'library(BiocManager)'
RUN R -e 'BiocManager::install("flowCore")'


RUN R -e "install.packages(c('shiny', 'shinydashboard', 'shinyjs', 'shinybusy', 'bslib', 'shinycssloaders', 'xgboost', 'dplyr', 'tidyr', 'plyr', 'gtools', 'igraph', 'Rcpp', 'DT', 'reticulate', 'openxlsx', 'parallel','stringr'), repos='http://cran.rstudio.com/')"

# python packages 

RUN pip3 install scyan
RUN pip3 install pandas
RUN pip3 install numpy


# copy the app to the image
RUN mkdir /root/cellAnnotationApp
COPY * /root/cellAnnotationApp

# redirect o  to cel
RUN ls -la /root/cellAnnotationApp/*


# configure port 

EXPOSE 3838

# clean install  
RUN apt-get clean\
  && apt-get remove --yes --purge build-essential

RUN mkdir -p /mnt


CMD ["R", "-e", "shiny::runApp(appDir = '/root/cellAnnotationApp', host='0.0.0.0', port=3838)"]
