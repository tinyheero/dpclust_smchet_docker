# DPClust Docker

This repository is a fork of https://github.com/Wedge-Oxford/dpclust_smchet_docker.

It provides an updated Docker image that is is more stable and uses an updated
version of the dpclust3p and dpclust packages. 

## How to Install

First we grab the most recent versions of dpclust3p and dpclust:

```bash
wget https://github.com/Wedge-Oxford/dpclust3p/archive/v1.0.8.tar.gz \
  -O dpclust_v1.0.8.tar.gz

wget https://github.com/Wedge-Oxford/dpclust/archive/v2.2.8.tar.gz \
  -O dpclust_v2.2.8.tar.gz
```

The `dpc.`R pipeline script takes 4 parameters and assumes the sex of the donor 
is male:

 * Mutect VCF file
 * There can be multiple samples in the VCF file, supply the number of this sample
 * Battenberg copy number output file
 * Battenberg purity output file
