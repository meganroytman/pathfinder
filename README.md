# Pathfinder

# Description
Pathfinder is a statistical fine-mapping framework that integrates genome, chromatin, and gene expression data in order to predict the causal SNP and mark within a gene region that are influencing gene expression. The code runs on a single locus at a time and takes as input the following data for each region:
1. Summary association statistics (Z-scores) between SNPs and marks, as well as marks and expression.
2. SNP LD and mark correlations matrices

Pathfinder outputs the posterior probability of causality for each possible SNP-mark-expression path throughout the region.

# Usage
Download pathfinder.R into your desired target directory. To run the software using R, type:

$ Rscript SampleData/assoc_1.txt SampleData/ld_snps_1.txt SampleData/ld_peaks_1.txt outfile.txt 5

The first line of outfile.txt (after the header) should read: 1 1 4.47714061045043e-12

For detailed information on input files and command line flags, see pathfinder.R.
