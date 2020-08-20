# Organization of R code

This folder contains the scripts for initializing datasets, generating results
and figures in the manuscript. The code is separated into folders by the respective
type of analysis conducted:

- [`init_reference_data`](https://github.com/ryurko/MAGMA-correspondence/blob/master/R/init_reference_data) - contains code necessary for preparing the 1000 Genomes reference data, randomly selecting genes
for comparing the different methods considered for computing gene-level p-values, and the follow-up resampling step,
- [`init_simulations`](https://github.com/ryurko/MAGMA-correspondence/blob/master/R/init_simulations) - contains all code for generating the null gene simulations as well as approximating the two-sided test covariance under the implication of Brown's original approximation,
- [`hmagma`](https://github.com/ryurko/MAGMA-correspondence/blob/master/R/hmagma) - contains all code for initializing the SNP-gene datasets with five phenotypes, processing through MAGMA, and then generating the corrected Monte Carlo-based p-values.
- [`create_figures`](https://github.com/ryurko/MAGMA-correspondence/blob/master/R/create_figures) - contains all code relevant to reproducing all figures in the manuscript, both the two main figures and each of the supplementary figures.

All code is presented under the assumption of working in a R Project version of this repository and that the necessary initial datasets are available in the `data` folder. See the README in the `data` folder for more information regarding the required initial datasets.

Additionally, the simulation scripts in [`init_simulations`](https://github.com/ryurko/MAGMA-correspondence/blob/master/R/init_simulations) rely on the usage of the [`snpcombineR`](https://github.com/ryurko/snpcombineR) package. `snpcombineR` implements [`Rcpp`](http://dirk.eddelbuettel.com/code/rcpp.html) code for the various methods to compute gene-level
p-values and the functions necessary for simulating GWAS summary statistics given
the reference data example genotypes.