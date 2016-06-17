# DPClust Docker
Legacy version of DPClust from before May 2014. This repo contains both a Dockerfile and a Galaxy XML file.

The dpc.R script takes 4 parameters and assumes the sex of the donor is male:
 * Mutect VCF file
 * There can be multiple samples in the VCF file, supply the number of this sample
 * Battenberg copy number output file
 * Cellularity as a float
