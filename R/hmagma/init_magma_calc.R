# PURPOSE: Generate the results for each phenotype with MAGMA for adult & fetal,
#          will then use this as a guide for ensuring the same number of SNPs
#          based on the synonyms file the MC version

library(tidyverse)
library(data.table)

# Load and prep the ASD file to be MAGMA input ----------------------------
asd_gwas_data <- fread("gunzip -c data/gwas/iPSYCH-PGC_ASD_Nov2017.gz")

# Save a version of this file that has just the SNP and p-value columns
asd_gwas_data[, .(SNP, P)] %>%
  .[, SIZE := rep(1000, nrow(asd_gwas_data))] %>%
  readr::write_tsv(path = "data/hmagma/magma_input/asd_gwas_snps_pvals.txt",
                   col_names = TRUE)

# Run MAGMA using the adult:
# magma --bfile Documents/asd-genome-wide-region-analysis/data/g1000_eur/g1000_eur --gene-annot Documents/asd-genome-wide-region-analysis/data/hmagma/snp_annotation/Adult_brain.genes.annot --pval Documents/asd-genome-wide-region-analysis/data/hmagma/magma_input/asd_gwas_snps_pvals.txt ncol = SIZE --genes-only --out Documents/asd-genome-wide-region-analysis/data/hmagma/output/test_magma/asd/adult_results
# For fetal:
# magma --bfile Documents/asd-genome-wide-region-analysis/data/g1000_eur/g1000_eur --gene-annot Documents/asd-genome-wide-region-analysis/data/hmagma/snp_annotation/Fetal_brain.genes.annot --pval Documents/asd-genome-wide-region-analysis/data/hmagma/magma_input/asd_gwas_snps_pvals.txt ncol = SIZE --genes-only --out Documents/asd-genome-wide-region-analysis/data/hmagma/output/test_magma/asd/fetal_results


# Repeat for BD -----------------------------------------------------------
rm(list = ls())

bd_gwas_data <- fread("gunzip -c data/gwas/bd/bd_2018.gz")

# Save a version of this file that has just the SNP and p-value columns
bd_gwas_data[, .(SNP, P)] %>%
  .[, SIZE := rep(1000, nrow(bd_gwas_data))] %>%
  readr::write_tsv(path = "data/hmagma/magma_input/bd_gwas_snps_pvals.txt",
                   col_names = TRUE)

# Run MAGMA using the adult:
# magma --bfile Documents/asd-genome-wide-region-analysis/data/g1000_eur/g1000_eur --gene-annot Documents/asd-genome-wide-region-analysis/data/hmagma/snp_annotation/Adult_brain.genes.annot --pval Documents/asd-genome-wide-region-analysis/data/hmagma/magma_input/bd_gwas_snps_pvals.txt ncol = SIZE --genes-only --out Documents/asd-genome-wide-region-analysis/data/hmagma/output/test_magma/bd/adult_results
# For fetal:
# magma --bfile Documents/asd-genome-wide-region-analysis/data/g1000_eur/g1000_eur --gene-annot Documents/asd-genome-wide-region-analysis/data/hmagma/snp_annotation/Fetal_brain.genes.annot --pval Documents/asd-genome-wide-region-analysis/data/hmagma/magma_input/bd_gwas_snps_pvals.txt ncol = SIZE --genes-only --out Documents/asd-genome-wide-region-analysis/data/hmagma/output/test_magma/bd/fetal_results

rm(list = ls())

# Repeat for ADHD ---------------------------------------------------------

adhd_gwas_data <- fread("data/gwas/adhd/adhd_jul2017")

# Save a version of this file that has just the SNP and p-value columns
adhd_gwas_data[, .(SNP, P)] %>%
  .[, SIZE := rep(1000, nrow(adhd_gwas_data))] %>%
  readr::write_tsv(path = "data/hmagma/magma_input/adhd_gwas_snps_pvals.txt",
                   col_names = TRUE)

# Run MAGMA using the adult:
# magma --bfile Documents/asd-genome-wide-region-analysis/data/g1000_eur/g1000_eur --gene-annot Documents/asd-genome-wide-region-analysis/data/hmagma/snp_annotation/Adult_brain.genes.annot --pval Documents/asd-genome-wide-region-analysis/data/hmagma/magma_input/adhd_gwas_snps_pvals.txt ncol = SIZE --genes-only --out Documents/asd-genome-wide-region-analysis/data/hmagma/output/test_magma/adhd/adult_results
# For fetal:
# magma --bfile Documents/asd-genome-wide-region-analysis/data/g1000_eur/g1000_eur --gene-annot Documents/asd-genome-wide-region-analysis/data/hmagma/snp_annotation/Fetal_brain.genes.annot --pval Documents/asd-genome-wide-region-analysis/data/hmagma/magma_input/adhd_gwas_snps_pvals.txt ncol = SIZE --genes-only --out Documents/asd-genome-wide-region-analysis/data/hmagma/output/test_magma/adhd/fetal_results

rm(list = ls())

# Repeat for MDD ----------------------------------------------------------

mdd_gwas_data <- fread("unzip -c data/gwas/mdd/mdd_gwas.zip", skip = 1)

# Save a version of this file that has just the MarkerName and p-value columns
mdd_gwas_data %>%
  .[, .(SNP = MarkerName, P)] %>%
  #.[, .(SNP, P)] %>%
  .[, SIZE := rep(1000, nrow(mdd_gwas_data))] %>%
  readr::write_tsv(path = "data/hmagma/magma_input/mdd_gwas_snps_pvals.txt",
                   col_names = TRUE)

# Run MAGMA using the adult:
# magma --bfile Documents/asd-genome-wide-region-analysis/data/g1000_eur/g1000_eur --gene-annot Documents/asd-genome-wide-region-analysis/data/hmagma/snp_annotation/Adult_brain.genes.annot --pval Documents/asd-genome-wide-region-analysis/data/hmagma/magma_input/mdd_gwas_snps_pvals.txt ncol = SIZE --genes-only --out Documents/asd-genome-wide-region-analysis/data/hmagma/output/test_magma/mdd/adult_results
# For fetal:
# magma --bfile Documents/asd-genome-wide-region-analysis/data/g1000_eur/g1000_eur --gene-annot Documents/asd-genome-wide-region-analysis/data/hmagma/snp_annotation/Fetal_brain.genes.annot --pval Documents/asd-genome-wide-region-analysis/data/hmagma/magma_input/mdd_gwas_snps_pvals.txt ncol = SIZE --genes-only --out Documents/asd-genome-wide-region-analysis/data/hmagma/output/test_magma/mdd/fetal_results

rm(list = ls())


# Repeat for SCZ & create its dataset with rs IDs -------------------------

scz_gwas_data <- fread("gunzip -c data/gwas/scz/hmagma_scz_gwas.gz")
# According to the manuscript the paper uses hg19 (GRCh37)
scz_gwas_data[, hg19_id := paste0("chr", CHR, ":", BP, "-", BP)]

# Now need to convert the SNP IMPUTE2 IDs to rs IDs - the simplest way to do
# this is by re-accessing the reference data, using those positions for constructing
# the IDs

library(snpStats)
g1000_eur_fam_path <- "data/g1000_eur/g1000_eur.fam"
g1000_eur_bim_path <- "data/g1000_eur/g1000_eur.bim"
g1000_eur_bed_path <- "data/g1000_eur/g1000_eur.bed"

# Read in the PLINK data using snpStats:
ref_snps_plink <- read.plink(g1000_eur_bed_path,
                             g1000_eur_bim_path,
                             g1000_eur_fam_path)
# Obtain the SNP information from ref_snps_plink list
ref_genobim <- ref_snps_plink$map
rm(ref_snps_plink)
colnames(ref_genobim) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")
ref_genobim <- as.data.table(ref_genobim)
ref_genobim[, hg19_id := paste0("chr", chr, ":", position, "-", position)]

# Only keep the reference data locations for those that are in the GWAS data:
gwas_ref_genobim <- ref_genobim[hg19_id %in% scz_gwas_data$hg19_id,]
nrow(gwas_ref_genobim)
# [1] 7618409
rm(ref_genobim)

# Next need to check that these have the same alleles. First rename the allele
# columns in the reference data:
setnames(gwas_ref_genobim, old = "A1", new = "ref_a1")
setnames(gwas_ref_genobim, old = "A2", new = "ref_a2")
setnames(gwas_ref_genobim, old = "SNP", new = "rsid")

# Only consider the genomic locations in the reference data and check to make
# sure that the alleles match (in either order):
scz_gwas_data <- scz_gwas_data %>%
  .[, A1 := toupper(A1)] %>%
  .[, A2 := toupper(A2)]
scz_gwas_data <- scz_gwas_data %>%
  .[hg19_id %in% gwas_ref_genobim$hg19_id,] %>%
  merge(gwas_ref_genobim[, .(hg19_id, ref_a1, ref_a2, rsid)],
        all.x = TRUE, by = "hg19_id") %>%
  # Create a new column denoting if the alleles match (in either direction):
  .[, ref_match := as.numeric((A1 == ref_a1 & A2 == ref_a2) |
                                (A1 == ref_a2 & A2 == ref_a1))]
table(scz_gwas_data$ref_match)
#    0       1
# 2617 7615891
scz_match_gwas_data <- scz_gwas_data %>%
  .[ref_match == 1,]
rm(scz_gwas_data)

# Save this dataset for use later:
fwrite(scz_match_gwas_data,
       "data/gwas/scz/scz_match_gwas_data.csv")

# Save a version of this file that has just the rsid and p-value columns
scz_match_gwas_data %>%
  .[, .(SNP = rsid, P)] %>%
  .[, SIZE := rep(1000, nrow(scz_match_gwas_data))] %>%
  readr::write_tsv(path = "data/hmagma/magma_input/scz_gwas_snps_pvals.txt",
                   col_names = TRUE)

# Run MAGMA using the adult:
# magma --bfile Documents/asd-genome-wide-region-analysis/data/g1000_eur/g1000_eur --gene-annot Documents/asd-genome-wide-region-analysis/data/hmagma/snp_annotation/Adult_brain.genes.annot --pval Documents/asd-genome-wide-region-analysis/data/hmagma/magma_input/scz_gwas_snps_pvals.txt ncol = SIZE --genes-only --out Documents/asd-genome-wide-region-analysis/data/hmagma/output/test_magma/scz/adult_results
# For fetal:
# magma --bfile Documents/asd-genome-wide-region-analysis/data/g1000_eur/g1000_eur --gene-annot Documents/asd-genome-wide-region-analysis/data/hmagma/snp_annotation/Fetal_brain.genes.annot --pval Documents/asd-genome-wide-region-analysis/data/hmagma/magma_input/scz_gwas_snps_pvals.txt ncol = SIZE --genes-only --out Documents/asd-genome-wide-region-analysis/data/hmagma/output/test_magma/scz/fetal_results



