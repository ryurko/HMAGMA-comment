# PURPOSE: Compare the updated counts in the number of genes I'm using for
#          H-MAGMA as compared to what they report

library(tidyverse)
library(data.table)

tidy_fetal_brain_gene_snp_data <-
  read_csv("data/hmagma/snp_annotation/server_input/tidy_final_fetal_gene_snp_pval_data.csv")
tidy_adult_brain_gene_snp_data <-
  read_csv("data/hmagma/snp_annotation/server_input/tidy_final_adult_gene_snp_pval_data.csv")

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(sheet_i) readxl::read_excel(filename,
                                                           sheet = sheet_i))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  return(x)
}

# Load the adult and fetal brain H-MAGMA results:
hmagma_adult_brain_results <-
  read_excel_allsheets("data/hmagma/output/reported/H-MAGMA_Adult_brain_output.xlsx")
hmagma_fetal_brain_results <-
  read_excel_allsheets("data/hmagma/output/reported/H-MAGMA_Fetal_brain_output.xlsx")

# Compare the count for each phenotype:

# ASD - adult
length(setdiff(
  hmagma_adult_brain_results$`H-MAGMA_ASD`$GENE,
  pull(filter(tidy_adult_brain_gene_snp_data, !is.na(asd_pval)), gene_id)))
# [1] 2
# ASD - fetal
length(setdiff(
  hmagma_fetal_brain_results$`H-MAGMA_ASD`$GENE,
  pull(filter(tidy_fetal_brain_gene_snp_data, !is.na(asd_pval)), gene_id)))
# [1] 1

# ADHD - adult
length(setdiff(
  hmagma_adult_brain_results$`H-MAGMA_ADHD`$GENE,
  pull(filter(tidy_adult_brain_gene_snp_data, !is.na(adhd_pval)), gene_id)))
# [1] 41
# ADHD - fetal
length(setdiff(
  hmagma_fetal_brain_results$`H-MAGMA_ADHD`$GENE,
  pull(filter(tidy_fetal_brain_gene_snp_data, !is.na(adhd_pval)), gene_id)))
# [1] 45

# SCZ - adult
length(setdiff(
  hmagma_adult_brain_results$`H-MAGMA_SCZ`$GENE,
  pull(filter(tidy_adult_brain_gene_snp_data, !is.na(scz_pval)), gene_id)))
# [1] 0
# SCZ - fetal
length(setdiff(
  hmagma_fetal_brain_results$`H-MAGMA_SCZ`$GENE,
  pull(filter(tidy_fetal_brain_gene_snp_data, !is.na(scz_pval)), gene_id)))
# [1] 1

# MDD - adult
length(setdiff(
  hmagma_adult_brain_results$`H-MAGMA_MDD`$GENE,
  pull(filter(tidy_adult_brain_gene_snp_data, !is.na(mdd_pval)), gene_id)))
# [1] 0
# MDD - fetal
length(setdiff(
  hmagma_fetal_brain_results$`H-MAGMA_MDD`$GENE,
  pull(filter(tidy_fetal_brain_gene_snp_data, !is.na(mdd_pval)), gene_id)))
# [1] 0

# BD - adult
length(setdiff(
  hmagma_adult_brain_results$`H-MAGMA_BD`$GENE,
  pull(filter(tidy_adult_brain_gene_snp_data, !is.na(bd_pval)), gene_id)))
# [1] 0
# BD - fetal
length(setdiff(
  hmagma_fetal_brain_results$`H-MAGMA_BD`$GENE,
  pull(filter(tidy_fetal_brain_gene_snp_data, !is.na(bd_pval)), gene_id)))
# [1] 0

# And now check with the regular use of MAGMA given the provided datasets
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

# ASD - adult
length(setdiff(
  asd_adult_magma_results$GENE,
  pull(filter(tidy_adult_brain_gene_snp_data, !is.na(asd_pval)), gene_id)))
# [1] 0
# ASD - fetal
length(setdiff(
  asd_fetal_magma_results$GENE,
  pull(filter(tidy_fetal_brain_gene_snp_data, !is.na(asd_pval)), gene_id)))
# [1] 0

# ADHD - adult
length(setdiff(
  adhd_adult_magma_results$GENE,
  pull(filter(tidy_adult_brain_gene_snp_data, !is.na(adhd_pval)), gene_id)))
# [1] 0
# ADHD - fetal
length(setdiff(
  adhd_fetal_magma_results$GENE,
  pull(filter(tidy_fetal_brain_gene_snp_data, !is.na(adhd_pval)), gene_id)))
# [1] 0

# SCZ - adult
length(setdiff(
  scz_adult_magma_results$GENE,
  pull(filter(tidy_adult_brain_gene_snp_data, !is.na(scz_pval)), gene_id)))
# [1] 0
# SCZ - fetal
length(setdiff(
  scz_fetal_magma_results$GENE,
  pull(filter(tidy_fetal_brain_gene_snp_data, !is.na(scz_pval)), gene_id)))
# [1] 0

# MDD - adult
length(setdiff(
  mdd_adult_magma_results$GENE,
  pull(filter(tidy_adult_brain_gene_snp_data, !is.na(mdd_pval)), gene_id)))
# [1] 0
# MDD - fetal
length(setdiff(
  mdd_fetal_magma_results$GENE,
  pull(filter(tidy_fetal_brain_gene_snp_data, !is.na(mdd_pval)), gene_id)))
# [1] 0

# BD - adult
length(setdiff(
  bd_adult_magma_results$GENE,
  pull(filter(tidy_adult_brain_gene_snp_data, !is.na(bd_pval)), gene_id)))
# [1] 0
# BD - fetal
length(setdiff(
  bd_fetal_magma_results$GENE,
  pull(filter(tidy_fetal_brain_gene_snp_data, !is.na(bd_pval)), gene_id)))
# [1] 0





