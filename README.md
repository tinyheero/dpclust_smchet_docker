# DPClust Docker

This repository is a fork of https://github.com/Wedge-Oxford/dpclust_smchet_docker.

It provides an updated Docker image that is is more stable and uses an updated
version of the dpclust3p and dpclust packages.

The `dpc.`R pipeline script takes 4 parameters and assumes the sex of the donor 
is male:

 * Mutect VCF file
 * There can be multiple samples in the VCF file, supply the number of this sample
 * Battenberg copy number output file
 * Battenberg purity output file
