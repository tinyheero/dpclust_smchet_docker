# DPClust Docker

This repository is a fork of https://github.com/Wedge-Oxford/dpclust_smchet_docker.

It provides an updated Docker image that builds properly.

## How to Install

```bash
docker build -t dpclust_smchet .
```

The `dpc.R` pipeline script takes 4 parameters and assumes the sex of the donor 
is male:

 * Mutect VCF file
 * There can be multiple samples in the VCF file, supply the number of this sample
 * Battenberg copy number output file
 * Battenberg purity output file

## Note

* The `dpc.R` script has been tested to work with only v1.0.6 of the dpclust3p 
    package. It does not work with v1.0.8.
