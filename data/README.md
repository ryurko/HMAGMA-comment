# Overview of data

This folder is meant to contain the files necessary for analysis, as well as 
store the results generated by the code. The following folders are assumed to exist:

- `reference_data` - containing folders `g1000_eur`, `pruned_genotype`, and `resampled_genotypes`,
- `updated_null_sims` - containing folders for the example gene simulations,
- `gwas` - contains the publicly available GWAS summary statistics for the considered phenotypes ADHD, ASD, BD, MDD, and SCZ,
- `hmagma` - contains folders `magma_input`, `output`, `reference_genotypes`, and `snp_annotation` for generating the H-MAGMA results with the corrected p-values.

The following datasets are __assumed to be in `reference_data` folder__ to start 
the analysis: (1) the 1000 Genomes reference dataset files `g1000_eur.fam`, `g1000_eur.bim`, `g1000_eur.bed` (available to download
[here](https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip)) and (2) the NCBI locations
of protein coding genes with build 37 (saved as `NCBI37.3/NCBI37.3.gene.loc` and available to download [here](https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI37.3.zip)).

Additionally, each of the necessary GWAS summary statistics are available from their respective sources:

- ADHD: Demontis, D, et al. Discovery of the first genome-wide significant risk loci for attention deficit/hyperactivity disorder. _Nature genetics_. __51__, 63–75 (2019).
- ASD: Grove, J, et al. Identification of common genetic risk variants for autism spectrum disorder. _Nature genetics_. __51__, 431–444 (2019).
- BD: Stahl, EA, .et al. Genome-wide association study identifies 30 loci associated with bipolar disorder. _Nature genetics_. __51__, 793–803 (2019).
- MDD: Howard, DM, .et al. Genome-wide association study of depression phenotypes in UK biobank identifies variants in excitatory synaptic pathways. _Nature communications_. __9__, 1–10 (2018).
- SCZ: Pardiñas, AF, .et al. Common schizophrenia alleles are enriched in mutation-intolerant genes and in regions under strong255background selection. _Nature genetics_. __50__, 381–389 (2018).