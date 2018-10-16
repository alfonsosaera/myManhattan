# myManhattan
I present here an R function to generate Manhattan plots using `ggplot`. It also returns a list of significant SNPs, according to different thresholds, if desired.

From [Wikipedia](https://en.wikipedia.org/wiki/Manhattan_plot), the free encyclopedia

A Manhattan plot is a type of scatter plot, usually used to display data with a large number of data-points - many of non-zero amplitude, and with a distribution of higher-magnitude values, for instance in genome-wide association studies (GWAS). In GWAS Manhattan plots, genomic coordinates are displayed along the X-axis, with the negative logarithm of the association P-value for each single nucleotide polymorphism (SNP) displayed on the Y-axis, meaning that each dot on the Manhattan plot signifies a SNP. Because the strongest associations have the smallest P-values (e.g., 10−15), their negative logarithms will be the greatest (e.g., 15).

The `myManhattanFunction.R` file contains the function and the `myManhattan.R` illustrates how to use it.
