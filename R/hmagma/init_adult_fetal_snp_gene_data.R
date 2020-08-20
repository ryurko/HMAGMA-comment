# PURPOSE: Initialize the SNP-gene datasets with the various GWAS results to
#          ensure the same number of SNPs that MAGMA assigns from the synonyms
#          file are also used in mcMAGMA.

library(tidyverse)
library(data.table)

# Load the H-MAGMA annotations --------------------------------------------

# First for the adult brain annotations:
tidy_adult_brain_gene_snp_data <- read_lines("data/hmagma/snp_annotation/Adult_brain.genes.annot") %>%
  lapply(
    function(gene) {
      gene_info <- unlist(str_split(gene, "\\t"))
      # Extract the SNPs:
      gene_snps <- gene_info[3:length(gene_info)]
      # Extract the gene location:
      gene_loc <- unlist(str_split(gene_info[2], ":"))
      tibble(snp_id = gene_snps) %>%
        mutate(gene_id = gene_info[1],
               gene_chr = gene_loc[1],
               gene_start = gene_loc[2],
               gene_end = gene_loc[3])
    }
  ) %>%
  bind_rows()

# Repeat for the fetal brain annotations:
tidy_fetal_brain_gene_snp_data <- read_lines("data/hmagma/snp_annotation/Fetal_brain.genes.annot") %>%
  lapply(
    function(gene) {
      gene_info <- unlist(str_split(gene, "\\t"))
      # Extract the SNPs:
      gene_snps <- gene_info[3:length(gene_info)]
      # Extract the gene location:
      gene_loc <- unlist(str_split(gene_info[2], ":"))
      tibble(snp_id = gene_snps) %>%
        mutate(gene_id = gene_info[1],
               gene_chr = gene_loc[1],
               gene_start = gene_loc[2],
               gene_end = gene_loc[3])
    }
  ) %>%
  bind_rows()


# Load the synonym file ---------------------------------------------------

# Read in the synonyms file - creating a list of the synonym SNP sets
synon_snp_lines <- read_lines("data/g1000_eur/g1000_eur.synonyms",
                              skip = 2) %>%
  lapply(., function(x) str_split(x, " "))
#names(synon_snp_lines) <- 1:length(synon_snp_lines)
synon_snp_dt <- rbindlist(synon_snp_lines, idcol = TRUE)
setnames(synon_snp_dt, old = c('.id', 'V1'), new = c('synon_id', 'snp_id'))
synon_snp_dt[, snp_id := paste0("rs", snp_id)]
rm(synon_snp_lines)
# Save this dataset so far:
fwrite(synon_snp_dt,
       "data/hmagma/snp_annotation/ref_synon_data.csv")

#synon_snp_df <- map_df(synon_snp_lines, ~as_tibble(.x), .id = "id")
# Bam it worked
#synon_snp_df <- synon_snp_df %>%
#  mutate(snp_id = paste0("rs", value))
#synon_snp_df <- dplyr::select(-value)

# Load reference data -----------------------------------------------------

library(snpStats)
g1000_eur_fam_path <- "data/g1000_eur/g1000_eur.fam"
g1000_eur_bim_path <- "data/g1000_eur/g1000_eur.bim"
g1000_eur_bed_path <- "data/g1000_eur/g1000_eur.bed"

# Read in the PLINK data using snpStats:
ref_snps_plink <- read.plink(g1000_eur_bed_path,
                             g1000_eur_bim_path,
                             g1000_eur_fam_path)

# Obtain the SnpMatrix object (genotypes) table from ref_snps_plink list
ref_genotypes <- ref_snps_plink$genotypes

# Obtain the SNP information from ref_snps_plink list
ref_genobim <- ref_snps_plink$map
rm(ref_snps_plink)
colnames(ref_genobim) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")
ref_genobim <- as.data.table(ref_genobim)


# Mark which SNPs are in reference data -----------------------------------

# For all three datasets (adult, fetal, synons) mark which SNPs are in the
# reference dataset:
tidy_adult_brain_gene_snp_data <- tidy_adult_brain_gene_snp_data %>%
  mutate(in_ref_data = as.numeric(snp_id %in% ref_genobim$SNP))
tidy_fetal_brain_gene_snp_data <- tidy_fetal_brain_gene_snp_data %>%
  mutate(in_ref_data = as.numeric(snp_id %in% ref_genobim$SNP))
synon_snp_dt[, in_ref_data := as.numeric(snp_id %in% ref_genobim$SNP)]


# Next need to mark which SNPs are in each GWAS ---------------------------


# First ASD:
asd_gwas_data <- fread("gunzip -c data/gwas/iPSYCH-PGC_ASD_Nov2017.gz")

# Left join for each of these the ASD p-values:
tidy_adult_brain_gene_snp_data <- tidy_adult_brain_gene_snp_data %>%
  left_join({
    asd_gwas_data[, list(SNP, asd_pval = P)] %>%
      as_tibble()
  }, by = c("snp_id" = "SNP")) %>%
  mutate(in_asd_data = as.numeric(!is.na(asd_pval)))
tidy_fetal_brain_gene_snp_data <- tidy_fetal_brain_gene_snp_data %>%
  left_join({
    asd_gwas_data[, list(SNP, asd_pval = P)] %>%
      as_tibble()
  }, by = c("snp_id" = "SNP")) %>%
  mutate(in_asd_data = as.numeric(!is.na(asd_pval)))
synon_snp_dt <- merge(synon_snp_dt,
                      asd_gwas_data[, list(snp_id = SNP, asd_pval = P)],
                      all.x = TRUE,
                      by = "snp_id")
synon_snp_dt[, in_asd_data := as.numeric(!is.na(asd_pval))]

rm(asd_gwas_data)

# Next BD:
bd_gwas_data <- fread("gunzip -c data/gwas/bd/bd_2018.gz")

# Left join for each of these the BD p-values:
tidy_adult_brain_gene_snp_data <- tidy_adult_brain_gene_snp_data %>%
  left_join({
    bd_gwas_data[, list(SNP, bd_pval = P)] %>%
      as_tibble()
  }, by = c("snp_id" = "SNP")) %>%
  mutate(in_bd_data = as.numeric(!is.na(bd_pval)))
tidy_fetal_brain_gene_snp_data <- tidy_fetal_brain_gene_snp_data %>%
  left_join({
    bd_gwas_data[, list(SNP, bd_pval = P)] %>%
      as_tibble()
  }, by = c("snp_id" = "SNP")) %>%
  mutate(in_bd_data = as.numeric(!is.na(bd_pval)))
synon_snp_dt <- merge(synon_snp_dt,
                      bd_gwas_data[, list(snp_id = SNP, bd_pval = P)],
                      all.x = TRUE,
                      by = "snp_id")
synon_snp_dt[, in_bd_data := as.numeric(!is.na(bd_pval))]

rm(bd_gwas_data)

# Next MDD:
mdd_gwas_data <- fread("unzip -c data/gwas/mdd/mdd_gwas.zip", skip = 1)

# Left join for each of these the MDD p-values:
tidy_adult_brain_gene_snp_data <- tidy_adult_brain_gene_snp_data %>%
  left_join({
    mdd_gwas_data[, list(MarkerName, mdd_pval = P)] %>%
      as_tibble()
  }, by = c("snp_id" = "MarkerName")) %>%
  mutate(in_mdd_data = as.numeric(!is.na(mdd_pval)))
tidy_fetal_brain_gene_snp_data <- tidy_fetal_brain_gene_snp_data %>%
  left_join({
    mdd_gwas_data[, list(MarkerName, mdd_pval = P)] %>%
      as_tibble()
  }, by = c("snp_id" = "MarkerName")) %>%
  mutate(in_mdd_data = as.numeric(!is.na(mdd_pval)))
synon_snp_dt <- merge(synon_snp_dt,
                      mdd_gwas_data[, list(snp_id = MarkerName, mdd_pval = P)],
                      all.x = TRUE,
                      by = "snp_id")
synon_snp_dt[, in_mdd_data := as.numeric(!is.na(mdd_pval))]

rm(mdd_gwas_data)

# Next ADHD:
adhd_gwas_data <- fread("data/gwas/adhd/adhd_jul2017")

# Left join for each of these the adhd p-values:
tidy_adult_brain_gene_snp_data <- tidy_adult_brain_gene_snp_data %>%
  left_join({
    adhd_gwas_data[, list(SNP, adhd_pval = P)] %>%
      as_tibble()
  }, by = c("snp_id" = "SNP")) %>%
  mutate(in_adhd_data = as.numeric(!is.na(adhd_pval)))
tidy_fetal_brain_gene_snp_data <- tidy_fetal_brain_gene_snp_data %>%
  left_join({
    adhd_gwas_data[, list(SNP, adhd_pval = P)] %>%
      as_tibble()
  }, by = c("snp_id" = "SNP")) %>%
  mutate(in_adhd_data = as.numeric(!is.na(adhd_pval)))
synon_snp_dt <- merge(synon_snp_dt,
                      adhd_gwas_data[, list(snp_id = SNP, adhd_pval = P)],
                      all.x = TRUE,
                      by = "snp_id")
synon_snp_dt[, in_adhd_data := as.numeric(!is.na(adhd_pval))]

rm(adhd_gwas_data)

# Finally SCZ:

scz_gwas_data <- fread("data/gwas/scz/scz_match_gwas_data.csv")

# Left join for each of these the scz p-values:
tidy_adult_brain_gene_snp_data <- tidy_adult_brain_gene_snp_data %>%
  left_join({
    scz_gwas_data[, list(rsid, scz_pval = P)] %>%
      as_tibble()
  }, by = c("snp_id" = "rsid")) %>%
  mutate(in_scz_data = as.numeric(!is.na(scz_pval)))
tidy_fetal_brain_gene_snp_data <- tidy_fetal_brain_gene_snp_data %>%
  left_join({
    scz_gwas_data[, list(rsid, scz_pval = P)] %>%
      as_tibble()
  }, by = c("snp_id" = "rsid")) %>%
  mutate(in_scz_data = as.numeric(!is.na(scz_pval)))
synon_snp_dt <- merge(synon_snp_dt,
                      scz_gwas_data[, list(snp_id = rsid, scz_pval = P)],
                      all.x = TRUE,
                      by = "snp_id")
synon_snp_dt[, in_scz_data := as.numeric(!is.na(scz_pval))]

rm(scz_gwas_data)



# Determine which SNPs to use for each phenotype --------------------------

# First load the MAGMA generated output to get the number of SNPs for each

asd_adult_magma_results <-
  fread("data/hmagma/output/test_magma/asd/adult_results.genes.out")
asd_fetal_magma_results <-
  fread("data/hmagma/output/test_magma/asd/fetal_results.genes.out")

bd_adult_magma_results <-
  fread("data/hmagma/output/test_magma/bd/adult_results.genes.out")
bd_fetal_magma_results <-
  fread("data/hmagma/output/test_magma/bd/fetal_results.genes.out")

mdd_adult_magma_results <-
  fread("data/hmagma/output/test_magma/mdd/adult_results.genes.out")
mdd_fetal_magma_results <-
  fread("data/hmagma/output/test_magma/mdd/fetal_results.genes.out")

adhd_adult_magma_results <-
  fread("data/hmagma/output/test_magma/adhd/adult_results.genes.out")
adhd_fetal_magma_results <-
  fread("data/hmagma/output/test_magma/adhd/fetal_results.genes.out")

scz_adult_magma_results <-
  fread("data/hmagma/output/test_magma/scz/adult_results.genes.out")
scz_fetal_magma_results <-
  fread("data/hmagma/output/test_magma/scz/fetal_results.genes.out")

# Next determine the number of SNPs for each IGNORING the synonyms file at first
# - then identifying the genes where the number of SNPs differ will guide which
# ones to fix for each phenotype

# Start with ASD:
asd_adult_no_syn_summary <- tidy_adult_brain_gene_snp_data %>%
  #filter(in_ref_data == 1, in_asd_data == 1) %>%
  group_by(gene_id) %>%
  summarize(n_snps = sum(in_ref_data * in_asd_data)) %>%
  dplyr::left_join(dplyr::select(asd_adult_magma_results,
                                 GENE, NSNPS), by = c("gene_id" = "GENE")) %>%
  # Remove the missing MAGMA ones since that is not resolved with synonyms
  filter(!is.na(NSNPS))
# Which ones don't match
asd_adult_mismatch <- asd_adult_no_syn_summary %>%
  filter(n_snps < NSNPS)

# Repeat for fetal
asd_fetal_no_syn_summary <- tidy_fetal_brain_gene_snp_data %>%
  group_by(gene_id) %>%
  summarize(n_snps = sum(in_ref_data * in_asd_data)) %>%
  dplyr::left_join(dplyr::select(asd_fetal_magma_results,
                                 GENE, NSNPS), by = c("gene_id" = "GENE")) %>%
  filter(!is.na(NSNPS))
asd_fetal_mismatch <- asd_fetal_no_syn_summary %>%
  filter(n_snps < NSNPS)

# Next BD:
bd_adult_no_syn_summary <- tidy_adult_brain_gene_snp_data %>%
  group_by(gene_id) %>%
  summarize(n_snps = sum(in_ref_data * in_bd_data)) %>%
  dplyr::left_join(dplyr::select(bd_adult_magma_results,
                                 GENE, NSNPS), by = c("gene_id" = "GENE")) %>%
  filter(!is.na(NSNPS))
bd_adult_mismatch <- bd_adult_no_syn_summary %>%
  filter(n_snps < NSNPS)
bd_fetal_no_syn_summary <- tidy_fetal_brain_gene_snp_data %>%
  group_by(gene_id) %>%
  summarize(n_snps = sum(in_ref_data * in_bd_data)) %>%
  dplyr::left_join(dplyr::select(bd_fetal_magma_results,
                                 GENE, NSNPS), by = c("gene_id" = "GENE")) %>%
  filter(!is.na(NSNPS))
bd_fetal_mismatch <- bd_fetal_no_syn_summary %>%
  filter(n_snps < NSNPS)

# MDD:
mdd_adult_no_syn_summary <- tidy_adult_brain_gene_snp_data %>%
  group_by(gene_id) %>%
  summarize(n_snps = sum(in_ref_data * in_mdd_data)) %>%
  dplyr::left_join(dplyr::select(mdd_adult_magma_results,
                                 GENE, NSNPS), by = c("gene_id" = "GENE")) %>%
  filter(!is.na(NSNPS))
mdd_adult_mismatch <- mdd_adult_no_syn_summary %>%
  filter(n_snps < NSNPS)
mdd_fetal_no_syn_summary <- tidy_fetal_brain_gene_snp_data %>%
  group_by(gene_id) %>%
  summarize(n_snps = sum(in_ref_data * in_mdd_data)) %>%
  dplyr::left_join(dplyr::select(mdd_fetal_magma_results,
                                 GENE, NSNPS), by = c("gene_id" = "GENE")) %>%
  filter(!is.na(NSNPS))
mdd_fetal_mismatch <- mdd_fetal_no_syn_summary %>%
  filter(n_snps < NSNPS)

# ADHD:
adhd_adult_no_syn_summary <- tidy_adult_brain_gene_snp_data %>%
  group_by(gene_id) %>%
  summarize(n_snps = sum(in_ref_data * in_adhd_data)) %>%
  dplyr::left_join(dplyr::select(adhd_adult_magma_results,
                                 GENE, NSNPS), by = c("gene_id" = "GENE")) %>%
  filter(!is.na(NSNPS))
adhd_adult_mismatch <- adhd_adult_no_syn_summary %>%
  filter(n_snps < NSNPS)
adhd_fetal_no_syn_summary <- tidy_fetal_brain_gene_snp_data %>%
  group_by(gene_id) %>%
  summarize(n_snps = sum(in_ref_data * in_adhd_data)) %>%
  dplyr::left_join(dplyr::select(adhd_fetal_magma_results,
                                 GENE, NSNPS), by = c("gene_id" = "GENE")) %>%
  filter(!is.na(NSNPS))
adhd_fetal_mismatch <- adhd_fetal_no_syn_summary %>%
  filter(n_snps < NSNPS)

# SCZ:
scz_adult_no_syn_summary <- tidy_adult_brain_gene_snp_data %>%
  group_by(gene_id) %>%
  summarize(n_snps = sum(in_ref_data * in_scz_data)) %>%
  dplyr::left_join(dplyr::select(scz_adult_magma_results,
                                 GENE, NSNPS), by = c("gene_id" = "GENE")) %>%
  filter(!is.na(NSNPS))
scz_adult_mismatch <- scz_adult_no_syn_summary %>%
  filter(n_snps < NSNPS)
scz_fetal_no_syn_summary <- tidy_fetal_brain_gene_snp_data %>%
  group_by(gene_id) %>%
  summarize(n_snps = sum(in_ref_data * in_scz_data)) %>%
  dplyr::left_join(dplyr::select(scz_fetal_magma_results,
                                 GENE, NSNPS), by = c("gene_id" = "GENE")) %>%
  filter(!is.na(NSNPS))
scz_fetal_mismatch <- scz_fetal_no_syn_summary %>%
  filter(n_snps < NSNPS)

# Construct a dataset for each synonym group that has the first SNP with the
# different pieces of information.

# Load the original synonymous data that has the order to use:
synon_order_table <- fread("data/hmagma/snp_annotation/ref_synon_data.csv") %>%
  as_tibble() %>%
  group_by(synon_id) %>%
  mutate(synon_i = 1:n()) %>%
  ungroup()

# Join the order
synon_group_snps <- synon_snp_dt %>%
  as_tibble() %>%
  # Join the order in the group:
  dplyr::left_join(dplyr::select(synon_order_table, -in_ref_data),
                   by = c("synon_id", "snp_id")) %>%
  group_by(synon_id) %>%
  # Sort by the group index (in increasing order)
  arrange(synon_i) %>%
  ungroup() %>%
  arrange(synon_id)

# Save this dataset:
write_csv(synon_group_snps,
          "data/hmagma/snp_annotation/synon_ref_gwas_group_data.csv")
# Remove the other synon datasets for memory:
rm(synon_order_table, synon_snp_dt)

# Make datasets that each piece of data (ref and gwas) for each group that has
# coverage for it (preserving the order from synon_i), starting with reference data:
synon_ref_snps <- synon_group_snps %>%
  filter(in_ref_data == 1) %>%
  group_by(synon_id) %>%
  summarize(syn_ref_snp_id = first(snp_id))

# ASD coverage:
synon_asd_snps <- synon_group_snps %>%
  filter(in_asd_data == 1) %>%
  group_by(synon_id) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(synon_id, snp_id, asd_pval) %>%
  dplyr::rename(syn_asd_snp_id = snp_id, syn_asd_pval = asd_pval)

# BD coverage:
synon_bd_snps <- synon_group_snps %>%
  filter(in_bd_data == 1) %>%
  group_by(synon_id) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(synon_id, snp_id, bd_pval) %>%
  dplyr::rename(syn_bd_snp_id = snp_id, syn_bd_pval = bd_pval)

# MDD coverage:
synon_mdd_snps <- synon_group_snps %>%
  filter(in_mdd_data == 1) %>%
  group_by(synon_id) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(synon_id, snp_id, mdd_pval) %>%
  dplyr::rename(syn_mdd_snp_id = snp_id, syn_mdd_pval = mdd_pval)

# ADHD coverage:
synon_adhd_snps <- synon_group_snps %>%
  filter(in_adhd_data == 1) %>%
  group_by(synon_id) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(synon_id, snp_id, adhd_pval) %>%
  dplyr::rename(syn_adhd_snp_id = snp_id, syn_adhd_pval = adhd_pval)

# SCZ coverage:
synon_scz_snps <- synon_group_snps %>%
  filter(in_scz_data == 1) %>%
  group_by(synon_id) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(synon_id, snp_id, scz_pval) %>%
  dplyr::rename(syn_scz_snp_id = snp_id, syn_scz_pval = scz_pval)


# Next proceed to construct separate datasets for each adult/fetal phenotype
# for the SNPs in the genes mismatching in each in order to identify which
# synonyms to use. First make indicator variables denoting if the SNP is in the
# synonym dataset:

tidy_adult_brain_gene_snp_data <- tidy_adult_brain_gene_snp_data %>%
  mutate(in_syn_data = as.numeric(snp_id %in% synon_snp_dt$snp_id))
tidy_fetal_brain_gene_snp_data <- tidy_fetal_brain_gene_snp_data %>%
  mutate(in_syn_data = as.numeric(snp_id %in% synon_snp_dt$snp_id))

# Create tables using these SNPs to join their respective synon_id:
synon_snp_table <- synon_group_snps %>%
  filter(snp_id %in% c(pull(filter(tidy_adult_brain_gene_snp_data, in_syn_data == 1),
                            snp_id),
                       pull(filter(tidy_fetal_brain_gene_snp_data, in_syn_data == 1),
                            snp_id))) %>%
  dplyr::select(snp_id, synon_id)

# Left join these datasets to both sets of results:
tidy_adult_brain_gene_snp_data <- tidy_adult_brain_gene_snp_data %>%
  dplyr::left_join(synon_snp_table, by = "snp_id")
tidy_fetal_brain_gene_snp_data <- tidy_fetal_brain_gene_snp_data %>%
  dplyr::left_join(synon_snp_table, by = "snp_id")


# Now join the columns with the syn group summaries:
tidy_adult_brain_gene_snp_data <- tidy_adult_brain_gene_snp_data %>%
  dplyr::left_join(synon_ref_snps, by = "synon_id") %>%
  dplyr::left_join(synon_asd_snps, by = "synon_id") %>%
  dplyr::left_join(synon_bd_snps, by = "synon_id") %>%
  dplyr::left_join(synon_mdd_snps, by = "synon_id") %>%
  dplyr::left_join(synon_scz_snps, by = "synon_id") %>%
  dplyr::left_join(synon_adhd_snps, by = "synon_id")
tidy_fetal_brain_gene_snp_data <- tidy_fetal_brain_gene_snp_data %>%
  dplyr::left_join(synon_ref_snps, by = "synon_id") %>%
  dplyr::left_join(synon_asd_snps, by = "synon_id") %>%
  dplyr::left_join(synon_bd_snps, by = "synon_id") %>%
  dplyr::left_join(synon_mdd_snps, by = "synon_id") %>%
  dplyr::left_join(synon_scz_snps, by = "synon_id") %>%
  dplyr::left_join(synon_adhd_snps, by = "synon_id")

# Make the columns with the necessary SNP and p-values to use for each dataset:
tidy_adult_brain_gene_snp_data <- tidy_adult_brain_gene_snp_data %>%
  mutate(use_ref_snp_id = ifelse(in_ref_data == 1, snp_id, syn_ref_snp_id),
         # ASD
         use_asd_snp_id = ifelse(in_asd_data == 1, snp_id, syn_asd_snp_id),
         use_asd_pval = ifelse(in_asd_data == 1, asd_pval, syn_asd_pval),
         # BD
         use_bd_snp_id = ifelse(in_bd_data == 1, snp_id, syn_bd_snp_id),
         use_bd_pval = ifelse(in_bd_data == 1, bd_pval, syn_bd_pval),
         # MDD
         use_mdd_snp_id = ifelse(in_mdd_data == 1, snp_id, syn_mdd_snp_id),
         use_mdd_pval = ifelse(in_mdd_data == 1, mdd_pval, syn_mdd_pval),
         # SCZ
         use_scz_snp_id = ifelse(in_scz_data == 1, snp_id, syn_scz_snp_id),
         use_scz_pval = ifelse(in_scz_data == 1, scz_pval, syn_scz_pval),
         # ADHD
         use_adhd_snp_id = ifelse(in_adhd_data == 1, snp_id, syn_adhd_snp_id),
         use_adhd_pval = ifelse(in_adhd_data == 1, adhd_pval, syn_adhd_pval))

tidy_fetal_brain_gene_snp_data <- tidy_fetal_brain_gene_snp_data %>%
  mutate(use_ref_snp_id = ifelse(in_ref_data == 1, snp_id, syn_ref_snp_id),
         # ASD
         use_asd_snp_id = ifelse(in_asd_data == 1, snp_id, syn_asd_snp_id),
         use_asd_pval = ifelse(in_asd_data == 1, asd_pval, syn_asd_pval),
         # BD
         use_bd_snp_id = ifelse(in_bd_data == 1, snp_id, syn_bd_snp_id),
         use_bd_pval = ifelse(in_bd_data == 1, bd_pval, syn_bd_pval),
         # MDD
         use_mdd_snp_id = ifelse(in_mdd_data == 1, snp_id, syn_mdd_snp_id),
         use_mdd_pval = ifelse(in_mdd_data == 1, mdd_pval, syn_mdd_pval),
         # SCZ
         use_scz_snp_id = ifelse(in_scz_data == 1, snp_id, syn_scz_snp_id),
         use_scz_pval = ifelse(in_scz_data == 1, scz_pval, syn_scz_pval),
         # ADHD
         use_adhd_snp_id = ifelse(in_adhd_data == 1, snp_id, syn_adhd_snp_id),
         use_adhd_pval = ifelse(in_adhd_data == 1, adhd_pval, syn_adhd_pval))


# Now check for each of these updated gene sets the coverage:

# Start with ASD:
asd_adult_w_syn_summary <- tidy_adult_brain_gene_snp_data %>%
  group_by(gene_id) %>%
  summarize(n_snps = sum(as.numeric(!is.na(use_ref_snp_id)) *
                           as.numeric(!is.na(use_asd_pval)))) %>%
  dplyr::left_join(dplyr::select(asd_adult_magma_results,
                                 GENE, NSNPS), by = c("gene_id" = "GENE")) %>%
  # Remove the missing MAGMA ones since that is not resolved with synonyms
  filter(!is.na(NSNPS))
# Any with less?
asd_adult_w_syn_summary %>% filter(n_snps < NSNPS)
# NO

# Repeat for fetal
asd_fetal_w_syn_summary <- tidy_fetal_brain_gene_snp_data %>%
  group_by(gene_id) %>%
  summarize(n_snps = sum(as.numeric(!is.na(use_ref_snp_id)) *
                           as.numeric(!is.na(use_asd_pval)))) %>%
  dplyr::left_join(dplyr::select(asd_fetal_magma_results,
                                 GENE, NSNPS), by = c("gene_id" = "GENE")) %>%
  filter(!is.na(NSNPS))
asd_fetal_w_syn_summary %>% filter(n_snps < NSNPS)
# NO

# Next BD:
bd_adult_w_syn_summary <- tidy_adult_brain_gene_snp_data %>%
  group_by(gene_id) %>%
  summarize(n_snps = sum(as.numeric(!is.na(use_ref_snp_id)) *
                           as.numeric(!is.na(use_bd_pval)))) %>%
  dplyr::left_join(dplyr::select(bd_adult_magma_results,
                                 GENE, NSNPS), by = c("gene_id" = "GENE")) %>%
  filter(!is.na(NSNPS))
bd_adult_w_syn_summary %>% filter(n_snps < NSNPS)
# NO
bd_fetal_w_syn_summary <- tidy_fetal_brain_gene_snp_data %>%
  group_by(gene_id) %>%
  summarize(n_snps = sum(as.numeric(!is.na(use_ref_snp_id)) *
                           as.numeric(!is.na(use_bd_pval)))) %>%
  dplyr::left_join(dplyr::select(bd_fetal_magma_results,
                                 GENE, NSNPS), by = c("gene_id" = "GENE")) %>%
  filter(!is.na(NSNPS))
bd_fetal_w_syn_summary %>% filter(n_snps < NSNPS)
# NO

# MDD:
mdd_adult_w_syn_summary <- tidy_adult_brain_gene_snp_data %>%
  group_by(gene_id) %>%
  summarize(n_snps = sum(as.numeric(!is.na(use_ref_snp_id)) *
                           as.numeric(!is.na(use_mdd_pval)))) %>%
  dplyr::left_join(dplyr::select(mdd_adult_magma_results,
                                 GENE, NSNPS), by = c("gene_id" = "GENE")) %>%
  filter(!is.na(NSNPS))
mdd_adult_w_syn_summary %>% filter(n_snps < NSNPS)
# NO
mdd_fetal_w_syn_summary <- tidy_fetal_brain_gene_snp_data %>%
  group_by(gene_id) %>%
  summarize(n_snps = sum(as.numeric(!is.na(use_ref_snp_id)) *
                           as.numeric(!is.na(use_mdd_pval)))) %>%
  dplyr::left_join(dplyr::select(mdd_fetal_magma_results,
                                 GENE, NSNPS), by = c("gene_id" = "GENE")) %>%
  filter(!is.na(NSNPS))
mdd_fetal_w_syn_summary %>% filter(n_snps < NSNPS)
# NO

# ADHD:
adhd_adult_w_syn_summary <- tidy_adult_brain_gene_snp_data %>%
  group_by(gene_id) %>%
  summarize(n_snps = sum(as.numeric(!is.na(use_ref_snp_id)) *
                           as.numeric(!is.na(use_adhd_pval)))) %>%
  dplyr::left_join(dplyr::select(adhd_adult_magma_results,
                                 GENE, NSNPS), by = c("gene_id" = "GENE")) %>%
  filter(!is.na(NSNPS))
adhd_adult_w_syn_summary %>% filter(n_snps < NSNPS)
# NO
adhd_fetal_w_syn_summary <- tidy_fetal_brain_gene_snp_data %>%
  group_by(gene_id) %>%
  summarize(n_snps = sum(as.numeric(!is.na(use_ref_snp_id)) *
                           as.numeric(!is.na(use_adhd_pval)))) %>%
  dplyr::left_join(dplyr::select(adhd_fetal_magma_results,
                                 GENE, NSNPS), by = c("gene_id" = "GENE")) %>%
  filter(!is.na(NSNPS))
adhd_fetal_w_syn_summary %>% filter(n_snps < NSNPS)
# NO

# SCZ:
scz_adult_w_syn_summary <- tidy_adult_brain_gene_snp_data %>%
  group_by(gene_id) %>%
  summarize(n_snps = sum(as.numeric(!is.na(use_ref_snp_id)) *
                           as.numeric(!is.na(use_scz_pval)))) %>%
  dplyr::left_join(dplyr::select(scz_adult_magma_results,
                                 GENE, NSNPS), by = c("gene_id" = "GENE")) %>%
  filter(!is.na(NSNPS))
scz_adult_w_syn_summary %>% filter(n_snps < NSNPS)
# NO
scz_fetal_w_syn_summary <- tidy_fetal_brain_gene_snp_data %>%
  group_by(gene_id) %>%
  summarize(n_snps = sum(as.numeric(!is.na(use_ref_snp_id)) *
                           as.numeric(!is.na(use_scz_pval)))) %>%
  dplyr::left_join(dplyr::select(scz_fetal_magma_results,
                                 GENE, NSNPS), by = c("gene_id" = "GENE")) %>%
  filter(!is.na(NSNPS))
scz_fetal_w_syn_summary %>% filter(n_snps < NSNPS)
# NO


# Finally construct the dataset with no duplicates ------------------------

tidy_no_dup_adult_data <- tidy_adult_brain_gene_snp_data %>%
  dplyr::select(gene_id, contains("use")) %>%
  filter(!is.na(use_ref_snp_id)) %>%
  distinct()
# Now check for each ref SNP gene combo there is only one p-value for each phenotype:
no_dup_adult_summary <- tidy_no_dup_adult_data %>%
  group_by(gene_id, use_ref_snp_id) %>%
  summarize(n_asd = length(unique(use_asd_snp_id)),
            n_bd = length(unique(use_bd_snp_id)),
            n_mdd = length(unique(use_mdd_snp_id)),
            n_scz = length(unique(use_scz_snp_id)),
            n_adhd = length(unique(use_adhd_snp_id)))
# Based on the count there appear to be a small number with duplicates, inspect these...
# There are a small number of SNPs with multiple p-values for ASD, BD, and ADHD
# (which was noticed earlier when the count increased via the left-join)
# What happens when I drop duplicate p-values for the SNPs that may be the problem here:
tidy_no_dup_pval_adult_data <- tidy_no_dup_adult_data %>%
  dplyr::select(gene_id, use_ref_snp_id, contains("pval")) %>%
  distinct()
# Hmm thats not enough... there are BD ones with different p-values so the simplest
# solution is to just take the first since that preserves the synonymous ordering
tidy_final_pval_adult_data <- tidy_no_dup_pval_adult_data %>%
  group_by(gene_id, use_ref_snp_id) %>%
  slice(1)

# Just repeat the final steps for the fetal data:
tidy_final_pval_fetal_data <- tidy_fetal_brain_gene_snp_data %>%
  dplyr::select(gene_id, contains("use")) %>%
  filter(!is.na(use_ref_snp_id)) %>%
  distinct() %>%
  dplyr::select(gene_id, use_ref_snp_id, contains("pval")) %>%
  group_by(gene_id, use_ref_snp_id) %>%
  slice(1)

# And check to see if all of the genes from each are in the MAGMA results:
all(unique(pull(filter(tidy_final_pval_adult_data, !is.na(use_asd_pval)),
                gene_id)) %in% asd_adult_magma_results$GENE)
# [1] TRUE
all(unique(pull(filter(tidy_final_pval_adult_data, !is.na(use_bd_pval)),
                gene_id)) %in% bd_adult_magma_results$GENE)
# [1] TRUE
all(unique(pull(filter(tidy_final_pval_adult_data, !is.na(use_mdd_pval)),
                gene_id)) %in% mdd_adult_magma_results$GENE)
# [1] TRUE
all(unique(pull(filter(tidy_final_pval_adult_data, !is.na(use_scz_pval)),
                gene_id)) %in% scz_adult_magma_results$GENE)
# [1] TRUE
all(unique(pull(filter(tidy_final_pval_adult_data, !is.na(use_adhd_pval)),
                gene_id)) %in% adhd_adult_magma_results$GENE)
# [1] TRUE


all(unique(pull(filter(tidy_final_pval_fetal_data, !is.na(use_asd_pval)),
                gene_id)) %in% asd_fetal_magma_results$GENE)
# [1] TRUE
all(unique(pull(filter(tidy_final_pval_fetal_data, !is.na(use_bd_pval)),
                gene_id)) %in% bd_fetal_magma_results$GENE)
# [1] TRUE
all(unique(pull(filter(tidy_final_pval_fetal_data, !is.na(use_mdd_pval)),
                gene_id)) %in% mdd_fetal_magma_results$GENE)
# [1] TRUE
all(unique(pull(filter(tidy_final_pval_fetal_data, !is.na(use_scz_pval)),
                gene_id)) %in% scz_fetal_magma_results$GENE)
# [1] TRUE
all(unique(pull(filter(tidy_final_pval_fetal_data, !is.na(use_adhd_pval)),
                gene_id)) %in% adhd_fetal_magma_results$GENE)
# [1] TRUE

# And the other way?
all(asd_adult_magma_results$GENE %in%
      unique(pull(filter(tidy_final_pval_adult_data, !is.na(use_asd_pval)),
                gene_id)))
# [1] TRUE
all(bd_adult_magma_results$GENE %in%
      unique(pull(filter(tidy_final_pval_adult_data, !is.na(use_bd_pval)),
                  gene_id)))
# [1] TRUE
all(mdd_adult_magma_results$GENE %in%
      unique(pull(filter(tidy_final_pval_adult_data, !is.na(use_mdd_pval)),
                  gene_id)))
# [1] TRUE
all(adhd_adult_magma_results$GENE %in%
      unique(pull(filter(tidy_final_pval_adult_data, !is.na(use_adhd_pval)),
                  gene_id)))
# [1] TRUE
all(scz_adult_magma_results$GENE %in%
      unique(pull(filter(tidy_final_pval_adult_data, !is.na(use_scz_pval)),
                  gene_id)))
# [1] TRUE

all(asd_fetal_magma_results$GENE %in%
      unique(pull(filter(tidy_final_pval_fetal_data, !is.na(use_asd_pval)),
                  gene_id)))
# [1] TRUE
all(bd_fetal_magma_results$GENE %in%
      unique(pull(filter(tidy_final_pval_fetal_data, !is.na(use_bd_pval)),
                  gene_id)))
# [1] TRUE
all(mdd_fetal_magma_results$GENE %in%
      unique(pull(filter(tidy_final_pval_fetal_data, !is.na(use_mdd_pval)),
                  gene_id)))
# [1] TRUE
all(adhd_fetal_magma_results$GENE %in%
      unique(pull(filter(tidy_final_pval_fetal_data, !is.na(use_adhd_pval)),
                  gene_id)))
# [1] TRUE
all(scz_fetal_magma_results$GENE %in%
      unique(pull(filter(tidy_final_pval_fetal_data, !is.na(use_scz_pval)),
                  gene_id)))
# [1] TRUE


# Now update the column names for these final datasets and save them prior to
# cleaning environment and creating their respective genotype matrices for fast
# server generation of their p-values

tidy_final_pval_adult_data <- tidy_final_pval_adult_data %>%
  dplyr::rename(snp_id = use_ref_snp_id,
                asd_pval = use_asd_pval, bd_pval = use_bd_pval,
                mdd_pval = use_mdd_pval, scz_pval = use_scz_pval,
                adhd_pval = use_adhd_pval)

tidy_final_pval_fetal_data <- tidy_final_pval_fetal_data %>%
  dplyr::rename(snp_id = use_ref_snp_id,
                asd_pval = use_asd_pval, bd_pval = use_bd_pval,
                mdd_pval = use_mdd_pval, scz_pval = use_scz_pval,
                adhd_pval = use_adhd_pval)

# Save these both:
write_csv(tidy_final_pval_adult_data,
          "data/hmagma/snp_annotation/tidy_final_adult_gene_snp_pval_data.csv")
write_csv(tidy_final_pval_fetal_data,
          "data/hmagma/snp_annotation/tidy_final_fetal_gene_snp_pval_data.csv")

# Clear the environment except for these datasets and the reference data
rm(list = setdiff(ls(), c("tidy_final_pval_adult_data",
                          "tidy_final_pval_fetal_data",
                          "ref_genobim", "ref_genotypes")))

# Initialize the reference genotype matrices -------------------------------

# Just a check are all the SNPs in the reference data:
all(tidy_final_pval_adult_data$snp_id %in% ref_genobim$SNP)
# [1] TRUE
all(tidy_final_pval_fetal_data$snp_id %in% ref_genobim$SNP)
# [1] TRUE
# GOOOOOOOOOOD


# Filter ref_genobim to the SNP positions for the genes of interest for each
# dataset separately:
adult_snp_ref_genobim <-
  ref_genobim[SNP %in% unique(c(tidy_final_pval_adult_data$snp_id)),]
fetal_snp_ref_genobim <-
  ref_genobim[SNP %in% unique(c(tidy_final_pval_fetal_data$snp_id)),]
rm(ref_genobim)

# Create two versions of the genotype data:
adult_snp_ref_genotypes <- ref_genotypes[, adult_snp_ref_genobim$SNP]
fetal_snp_ref_genotypes <- ref_genotypes[, fetal_snp_ref_genobim$SNP]
# Remove the large version:
rm(ref_genotypes)

# Create the numeric versions of these matrices with alleles encoded as 0, 1, or 2:
numeric_adult_snp_ref_genotypes <- as(adult_snp_ref_genotypes, "numeric")
rm(adult_snp_ref_genotypes)

numeric_fetal_snp_ref_genotypes <- as(fetal_snp_ref_genotypes, "numeric")
rm(fetal_snp_ref_genotypes)

# Because these are so big - will NOT flip genotype data yet - but instead
# will do so after filtering to a group's SNPs and working with a much smaller
# matrix inside the test stat calculation

# Remove the row names:
rownames(numeric_adult_snp_ref_genotypes) <- NULL
# Change the colnames to be SNPs ids instead for consistency:
colnames(numeric_adult_snp_ref_genotypes) <- adult_snp_ref_genobim$SNP
rm(adult_snp_ref_genobim)

# Repeat for fetal:
rownames(numeric_fetal_snp_ref_genotypes) <- NULL
colnames(numeric_fetal_snp_ref_genotypes) <- fetal_snp_ref_genobim$SNP
rm(fetal_snp_ref_genobim)


# Create the server versions of these datasets to use ---------------------

# For both the fetal and adult datasets, assign the genes to 15 different folds
# First for adult
set.seed(1971)
adult_gene_server_fold <- tidy_final_pval_adult_data %>%
  ungroup() %>%
  dplyr::select(gene_id) %>%
  distinct() %>%
  mutate(server_fold = sample(rep(1:15, length.out = n())))
# Left-join this to the tidy-dataset then save the tidy data in the server folder:
tidy_final_pval_adult_data <- tidy_final_pval_adult_data %>%
  ungroup() %>%
  dplyr::left_join(adult_gene_server_fold,
                   by = "gene_id")
write_csv(tidy_final_pval_adult_data,
          "data/hmagma/snp_annotation/server_input/tidy_final_adult_gene_snp_pval_data.csv")

# Repeat for fetal:
set.seed(1979)
fetal_gene_server_fold <- tidy_final_pval_fetal_data %>%
  ungroup() %>%
  dplyr::select(gene_id) %>%
  distinct() %>%
  mutate(server_fold = sample(rep(1:15, length.out = n())))
tidy_final_pval_fetal_data <- tidy_final_pval_fetal_data %>%
  ungroup() %>%
  dplyr::left_join(fetal_gene_server_fold,
                   by = "gene_id")
write_csv(tidy_final_pval_fetal_data,
          "data/hmagma/snp_annotation/server_input/tidy_final_fetal_gene_snp_pval_data.csv")



# Now for each fold - grab the SNPs, and save that subset's genotype matrix:
for (fold_i in 1:max(adult_gene_server_fold$server_fold)) {

  adult_fold_snps <- tidy_final_pval_adult_data %>%
    filter(server_fold == fold_i) %>%
    pull(snp_id)
  fetal_fold_snps <- tidy_final_pval_fetal_data %>%
    filter(server_fold == fold_i) %>%
    pull(snp_id)

  # Save the genotype matrices:
  write_csv(as_tibble(numeric_adult_snp_ref_genotypes[, adult_fold_snps]),
            paste0("data/hmagma/reference_genotypes/adult/fold",
                   fold_i, "_ref_genotypes.csv"))

  write_csv(as_tibble(numeric_fetal_snp_ref_genotypes[, fetal_fold_snps]),
            paste0("data/hmagma/reference_genotypes/fetal/fold",
                   fold_i, "_ref_genotypes.csv"))

}



