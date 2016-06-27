FROM r-base

# Add dependencies
RUN apt-get update && apt-get install -y libxml2 libxml2-dev libcurl4-gnutls-dev

# Add dpclust3p dependencies
RUN R -q -e 'source("http://bioconductor.org/biocLite.R"); biocLite("VariantAnnotation")'

# Preprocessing
ADD dpclust3p_v1.0.1.tar.gz /opt/galaxy/tools/dpclust3p_v1.0.1.tar.gz
RUN R -q -e 'install.packages("/opt/galaxy/tools/dpclust3p_v1.0.1.tar.gz", type="source", repos=NULL)'

# DPClust pipeline
ADD dpc.R /opt/galaxy/tools/dpclust/dpc.R
ADD dpclust_cellularity.R /opt/galaxy/tools/dpclust/dpclust_cellularity.R
