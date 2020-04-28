# DPClust Docker

This repository is a fork of https://github.com/Wedge-Oxford/dpclust_smchet_docker.

It provides an updated Docker image that builds properly.

## How to Install

```bash
docker build -t dpclust_smchet .
```

The `dpc.R` pipeline script takes 5 parameters and assumes the sex of the donor 
is male:

1. Mutect VCF file
1. There can be multiple samples in the VCF file, supply the number of this sample
1. Battenberg copy number output file
1. Battenberg purity output file
1. Output directory

## How to Run

To run:

```bash
docker run --rm -v /:/host/ -w "/host/$(pwd)" -it dpclust_smchet \
    Rscript /opt/galaxy/tools/dpclust/dpc.R \
        /host/path/to/T2-noXY.truth.scoring_vcf.vcf \
        1 \
        /host/path/to/T2-noXY/T2-noXY.battenberg.txt \
        /host/path/to/T2-noXY/T2-noXY.cellularity_ploidy.txt \
        T2_output
```

## Note

* The `dpc.R` script has been tested to work with only v1.0.6 of the dpclust3p 
    package. It does not work with v1.0.8.
