FROM r-base

# Add dependencies
RUN apt-get update && apt-get install -y libxml2 libxml2-dev libcurl4-gnutls-dev

# Add dpclust(3p) dependencies
RUN R -q -e 'source("http://bioconductor.org/biocLite.R"); biocLite(c("VariantAnnotation","mcclust","KernSmooth","ks","lattice","ggplot2","gridExtra","reshape2"))'

# Preprocessing
ADD dpclust3p_dev6b2f713c9d.tar.gz /opt/galaxy/tools/dpclust3p_dev6b2f713c9d.tar.gz 
RUN R -q -e 'install.packages("/opt/galaxy/tools/dpclust3p_dev6b2f713c9d.tar.gz", type="source", repos=NULL)'

# DPClust
ADD dpclust_dev8468de6df2.tar.gz /opt/galaxy/tools/dpclust_dev8468de6df2.tar.gz
RUN R -q -e 'install.packages("/opt/galaxy/tools/dpclust_dev8468de6df2.tar.gz", type="source", repos=NULL)'

# DPClust pipeline
ADD dpc.R /opt/galaxy/tools/dpclust/dpc.R
ADD dpclust_cellularity.R /opt/galaxy/tools/dpclust/dpclust_cellularity.R
