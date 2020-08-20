# Overview of manuscript repository

This repository contains the code for reproducing the results and figures in
our manuscript _H-MAGMA, inheriting a shaky statistical foundation, yields excess false positives_. We explain the flaws of computing gene-level p-values based on two-sided summary statistics with the [MAGMA](https://ctg.cncr.nl/software/magma) software  and its implications on a recent publication [H-MAGMA](https://www.nature.com/articles/s41593-020-0603-0).

The folders are organized in the following manner:

- [`R`](https://github.com/ryurko/HMAGMA-comment/blob/master/R) - all scripts for initializing datasets, generating results and figures in manuscript,
- [`figures`](https://github.com/ryurko/HMAGMA-comment/blob/master/figures) - files for final figures displayed in the manuscript (including supplementary materials),
- [`data`](https://github.com/ryurko/HMAGMA-comment/blob/master/data) - folder
containing files necessary for the analysis.

Note that the [`data`](https://github.com/ryurko/HMAGMA-comment/blob/master/data)
folder does not contain any files in this repository beside the README files
explaining the layout necessary for reproducing the results using the code in the
[`R`](https://github.com/ryurko/HMAGMA-comment/blob/master/R) folder. This
is because the accessed sources of data are from the 1000 Genomes project 
[available to download from the MAGMA software site](https://ctg.cncr.nl/software/magma), publicly available GWAS summary statistics, and the H-MAGMA SNP-gene annotation files [which are available from their public GitHub repository](https://github.com/thewonlab/H-MAGMA).

Additionally, the simulation scripts in [`R`](https://github.com/ryurko/HMAGMA-comment/blob/master/R) rely on the usage of the [`snpcombineR`](https://github.com/ryurko/snpcombineR) package.

## Contact

Ron Yurko: [ryurko@andrew.cmu.edu](mailto:ryurko@andrew.cmu.edu)

## References

- Sey, Nancy YA, et al. A computational tool (H-MAGMA) for improved prediction of brain-disorder risk genes by incorporating brain chromatin interaction profiles. _Nature Neuroscience_. DOI: [10.1038/s41593-020-0603-0](https://www.nature.com/articles/s41593-020-0603-0) (2020).
- de Leeuw, C. A., Mooij, J. M., Heskes, T. & Posthuma, D.  MAGMA: Generalized Gene-Set Analysis of GWAS Data. _PLOS Comput. Biol_. __11__, 1–19, DOI: [10.1371/journal.pcbi.1004219](10.1371/journal.pcbi.1004219) (2015).
- 1000 Genomes Project Consortium and others. An integrated map of genetic variation from 1,092 human genomes. _Nature_ __491__, 56 (2012).
