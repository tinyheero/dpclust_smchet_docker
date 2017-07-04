# DPClust Docker

This repo contains both a Dockerfile and a Galaxy XML file for DPClust that plug into [SMC-het](https://www.synapse.org/#!Synapse:syn2813581).

The dpc.R pipeline script takes 4 parameters and assumes the sex of the donor is male:
 * Mutect VCF file
 * There can be multiple samples in the VCF file, supply the number of this sample
 * Battenberg copy number output file
 * Battenberg purity output file
