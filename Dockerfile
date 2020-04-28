FROM rocker/r-ver:3.6.3

# Add dependencies
RUN apt-get update && apt-get install -y \
  libxml2 \
  libxml2-dev \
  libcurl4-gnutls-dev \
  r-cran-rgl \
  libssl-dev

ADD Rprofile /.Rprofile

RUN R -e "install.packages(c('mcclust', 'KernSmooth', 'ks', 'ggplot2', 'gridExtra', 'reshape2'))"

# Add dpclust(3p) dependencies
RUN R -e "install.packages('BiocManager'); BiocManager::install('VariantAnnotation')"

# Install the dpclust3p package
COPY dpclust3p_v1.0.6.tar.gz /opt/galaxy/tools/dpclust3p_v1.0.6.tar.gz 
RUN R -q -e 'install.packages("/opt/galaxy/tools/dpclust3p_v1.0.6.tar.gz", type="source", repos=NULL)'

# Intall the dpclust package
COPY dpclust_v2.2.5.tar.gz /opt/galaxy/tools/dpclust_v2.2.5.tar.gz
RUN R -q -e 'install.packages("/opt/galaxy/tools/dpclust_v2.2.5.tar.gz", type="source", repos=NULL)'

# DPClust pipeline
ADD dpc.R /opt/galaxy/tools/dpclust/dpc.R
ADD dpclust_cellularity.R /opt/galaxy/tools/dpclust/dpclust_cellularity.R
