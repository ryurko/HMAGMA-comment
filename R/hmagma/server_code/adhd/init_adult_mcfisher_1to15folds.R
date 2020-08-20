# PURPOSE: Initialize the H-MAGMA adhd mcMAGMA p-values

# SERVER SCRIPT: folds 1 - 15

# AUTHOR: Ron Yurko

# ------------------------------------------------------------------------------

# Access necessary packages

library(backports, lib.loc = "/home/ryurko/Rpackages")
library(vctrs, lib.loc = "/home/ryurko/Rpackages")
library(dplyr, lib.loc = "/home/ryurko/Rpackages")
library(magrittr, lib.loc = "/home/ryurko/Rpackages")
library(data.table, lib.loc = "/home/ryurko/Rpackages")
library(snpcombineR, lib.loc = "/home/ryurko/Rpackages")
library(readr, lib.loc = "/home/ryurko/Rpackages")
library(parallel, lib.loc = "/home/ryurko/Rpackages")


# ------------------------------------------------------------------------------




# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


flip_matrix <- function(x) {
  zero2 <- which(x == 0)
  two0 <- which(x == 2)
  x[zero2] <- 2
  x[two0] <- 0
  return(x)
}

# Repeat the following for the three different folds:
tidy_adult_brain_gene_snp_data <-
  read_csv("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/tidy_snp_gene_data/tidy_final_adult_gene_snp_pval_data.csv") %>%
  filter(!is.na(adhd_pval),
         server_fold %in% 1)
numeric_snp_ref_genotypes <- as.matrix(
  data.table::fread("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/reference_data/adult/fold1_ref_genotypes.csv"))
gene_ids <- unique(tidy_adult_brain_gene_snp_data$gene_id)


tidy_group_pval_summary <-
  mclapply(1:length(gene_ids), mc.cores = 30,
           FUN = function(i) {

             group_x <- gene_ids[i]
             # Filter to the SNPs belonging to group X in the dataset:
             group_x_snp_rows <- tidy_adult_brain_gene_snp_data %>%
               filter(gene_id == group_x)

             # If only 1 SNP then just make a dataset with
             # that SNP's single p-value and test statistics
             # just for reference:
             if (nrow(group_x_snp_rows) == 1) {

               group_results <- group_x_snp_rows %>%
                 dplyr::select(gene_id,
                               adhd_pval) %>%
                 mutate(fisher_stat = -2 * log(adhd_pval),
                        fisher_pval = adhd_pval,
                        n_snps = 1) %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             } else {

               # Pull the SNPs:
               group_x_snps <- pull(group_x_snp_rows, snp_id)

               # Only look at these SNPs' genotype data - with the genotype
               # columns in the correct order:
               group_geno_ref_matrix <- numeric_snp_ref_genotypes[, group_x_snps]
               # Flip the matrix:
               group_geno_ref_matrix <- flip_matrix(group_geno_ref_matrix)

               # Compute the group's correlation matrix:
               group_cor_matrix <- compute_cor_matrix(group_geno_ref_matrix)

               # Compute the simulated truth statistics with this correlation matrix
               # using 300000 simulations to start:
               group_sim_truth <-
                 generate_sim_truth_data(group_cor_matrix,
                                         120000, 10, .Machine$double.eps)

               # The first column is the VEGAS data, second column is
               # the Fisher data - use these along with the observed
               # test statistics to compute empirical p-values and
               # return:
               group_results <- group_x_snp_rows %>%
                 group_by(gene_id) %>%
                 summarize(fisher_stat = -2 * sum(log(adhd_pval)),
                           n_snps = n()) %>%
                 ungroup() %>%
                 mutate(fisher_pval = mean(as.numeric(group_sim_truth[,2] > fisher_stat)))  %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             }

             group_results
           }) %>%
  bind_rows()
write_csv(tidy_group_pval_summary,
          "/home/ryurko/asd_genome_wide_analysis/hmagma/data/output/adhd/adult/tidy_gene_pval_summary_1.csv")

# ------------------------------------------------------------------------------

# Remove the datasets and repeat for the second fold:
rm(list = ls())


flip_matrix <- function(x) {
  zero2 <- which(x == 0)
  two0 <- which(x == 2)
  x[zero2] <- 2
  x[two0] <- 0
  return(x)
}

# Repeat the following for the three different folds:
tidy_adult_brain_gene_snp_data <-
  read_csv("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/tidy_snp_gene_data/tidy_final_adult_gene_snp_pval_data.csv") %>%
  filter(!is.na(adhd_pval),
         server_fold %in% 2)
numeric_snp_ref_genotypes <- as.matrix(
  data.table::fread("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/reference_data/adult/fold2_ref_genotypes.csv"))
gene_ids <- unique(tidy_adult_brain_gene_snp_data$gene_id)


tidy_group_pval_summary <-
  mclapply(1:length(gene_ids), mc.cores = 30,
           FUN = function(i) {

             group_x <- gene_ids[i]
             # Filter to the SNPs belonging to group X in the dataset:
             group_x_snp_rows <- tidy_adult_brain_gene_snp_data %>%
               filter(gene_id == group_x)

             # If only 1 SNP then just make a dataset with
             # that SNP's single p-value and test statistics
             # just for reference:
             if (nrow(group_x_snp_rows) == 1) {

               group_results <- group_x_snp_rows %>%
                 dplyr::select(gene_id,
                               adhd_pval) %>%
                 mutate(fisher_stat = -2 * log(adhd_pval),
                        fisher_pval = adhd_pval,
                        n_snps = 1) %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             } else {

               # Pull the SNPs:
               group_x_snps <- pull(group_x_snp_rows, snp_id)

               # Only look at these SNPs' genotype data - with the genotype
               # columns in the correct order:
               group_geno_ref_matrix <- numeric_snp_ref_genotypes[, group_x_snps]
               # Flip the matrix:
               group_geno_ref_matrix <- flip_matrix(group_geno_ref_matrix)

               # Compute the group's correlation matrix:
               group_cor_matrix <- compute_cor_matrix(group_geno_ref_matrix)

               # Compute the simulated truth statistics with this correlation matrix
               # using 300000 simulations to start:
               group_sim_truth <-
                 generate_sim_truth_data(group_cor_matrix,
                                         120000, 10, .Machine$double.eps)

               # The first column is the VEGAS data, second column is
               # the Fisher data - use these along with the observed
               # test statistics to compute empirical p-values and
               # return:
               group_results <- group_x_snp_rows %>%
                 group_by(gene_id) %>%
                 summarize(fisher_stat = -2 * sum(log(adhd_pval)),
                           n_snps = n()) %>%
                 ungroup() %>%
                 mutate(fisher_pval = mean(as.numeric(group_sim_truth[,2] > fisher_stat)))  %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             }

             group_results
           }) %>%
  bind_rows()
write_csv(tidy_group_pval_summary,
          "/home/ryurko/asd_genome_wide_analysis/hmagma/data/output/adhd/adult/tidy_gene_pval_summary_2.csv")

# ------------------------------------------------------------------------------

# Remove the datasets and repeat for the second fold:
rm(list = ls())


flip_matrix <- function(x) {
  zero2 <- which(x == 0)
  two0 <- which(x == 2)
  x[zero2] <- 2
  x[two0] <- 0
  return(x)
}

# Repeat the following for the three different folds:
tidy_adult_brain_gene_snp_data <-
  read_csv("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/tidy_snp_gene_data/tidy_final_adult_gene_snp_pval_data.csv") %>%
  filter(!is.na(adhd_pval),
         server_fold %in% 3)
numeric_snp_ref_genotypes <- as.matrix(
  data.table::fread("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/reference_data/adult/fold3_ref_genotypes.csv"))
gene_ids <- unique(tidy_adult_brain_gene_snp_data$gene_id)


tidy_group_pval_summary <-
  mclapply(1:length(gene_ids), mc.cores = 30,
           FUN = function(i) {

             group_x <- gene_ids[i]
             # Filter to the SNPs belonging to group X in the dataset:
             group_x_snp_rows <- tidy_adult_brain_gene_snp_data %>%
               filter(gene_id == group_x)

             # If only 1 SNP then just make a dataset with
             # that SNP's single p-value and test statistics
             # just for reference:
             if (nrow(group_x_snp_rows) == 1) {

               group_results <- group_x_snp_rows %>%
                 dplyr::select(gene_id,
                               adhd_pval) %>%
                 mutate(fisher_stat = -2 * log(adhd_pval),
                        fisher_pval = adhd_pval,
                        n_snps = 1) %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             } else {

               # Pull the SNPs:
               group_x_snps <- pull(group_x_snp_rows, snp_id)

               # Only look at these SNPs' genotype data - with the genotype
               # columns in the correct order:
               group_geno_ref_matrix <- numeric_snp_ref_genotypes[, group_x_snps]
               # Flip the matrix:
               group_geno_ref_matrix <- flip_matrix(group_geno_ref_matrix)

               # Compute the group's correlation matrix:
               group_cor_matrix <- compute_cor_matrix(group_geno_ref_matrix)

               # Compute the simulated truth statistics with this correlation matrix
               # using 300000 simulations to start:
               group_sim_truth <-
                 generate_sim_truth_data(group_cor_matrix,
                                         120000, 10, .Machine$double.eps)

               # The first column is the VEGAS data, second column is
               # the Fisher data - use these along with the observed
               # test statistics to compute empirical p-values and
               # return:
               group_results <- group_x_snp_rows %>%
                 group_by(gene_id) %>%
                 summarize(fisher_stat = -2 * sum(log(adhd_pval)),
                           n_snps = n()) %>%
                 ungroup() %>%
                 mutate(fisher_pval = mean(as.numeric(group_sim_truth[,2] > fisher_stat)))  %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             }

             group_results
           }) %>%
  bind_rows()
write_csv(tidy_group_pval_summary,
          "/home/ryurko/asd_genome_wide_analysis/hmagma/data/output/adhd/adult/tidy_gene_pval_summary_3.csv")

# ------------------------------------------------------------------------------

# Remove the datasets and repeat for the second fold:
rm(list = ls())


flip_matrix <- function(x) {
  zero2 <- which(x == 0)
  two0 <- which(x == 2)
  x[zero2] <- 2
  x[two0] <- 0
  return(x)
}

# Repeat the following for the three different folds:
tidy_adult_brain_gene_snp_data <-
  read_csv("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/tidy_snp_gene_data/tidy_final_adult_gene_snp_pval_data.csv") %>%
  filter(!is.na(adhd_pval),
         server_fold %in% 4)
numeric_snp_ref_genotypes <- as.matrix(
  data.table::fread("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/reference_data/adult/fold4_ref_genotypes.csv"))
gene_ids <- unique(tidy_adult_brain_gene_snp_data$gene_id)


tidy_group_pval_summary <-
  mclapply(1:length(gene_ids), mc.cores = 30,
           FUN = function(i) {

             group_x <- gene_ids[i]
             # Filter to the SNPs belonging to group X in the dataset:
             group_x_snp_rows <- tidy_adult_brain_gene_snp_data %>%
               filter(gene_id == group_x)

             # If only 1 SNP then just make a dataset with
             # that SNP's single p-value and test statistics
             # just for reference:
             if (nrow(group_x_snp_rows) == 1) {

               group_results <- group_x_snp_rows %>%
                 dplyr::select(gene_id,
                               adhd_pval) %>%
                 mutate(fisher_stat = -2 * log(adhd_pval),
                        fisher_pval = adhd_pval,
                        n_snps = 1) %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             } else {

               # Pull the SNPs:
               group_x_snps <- pull(group_x_snp_rows, snp_id)

               # Only look at these SNPs' genotype data - with the genotype
               # columns in the correct order:
               group_geno_ref_matrix <- numeric_snp_ref_genotypes[, group_x_snps]
               # Flip the matrix:
               group_geno_ref_matrix <- flip_matrix(group_geno_ref_matrix)

               # Compute the group's correlation matrix:
               group_cor_matrix <- compute_cor_matrix(group_geno_ref_matrix)

               # Compute the simulated truth statistics with this correlation matrix
               # using 300000 simulations to start:
               group_sim_truth <-
                 generate_sim_truth_data(group_cor_matrix,
                                         120000, 10, .Machine$double.eps)

               # The first column is the VEGAS data, second column is
               # the Fisher data - use these along with the observed
               # test statistics to compute empirical p-values and
               # return:
               group_results <- group_x_snp_rows %>%
                 group_by(gene_id) %>%
                 summarize(fisher_stat = -2 * sum(log(adhd_pval)),
                           n_snps = n()) %>%
                 ungroup() %>%
                 mutate(fisher_pval = mean(as.numeric(group_sim_truth[,2] > fisher_stat)))  %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             }

             group_results
           }) %>%
  bind_rows()
write_csv(tidy_group_pval_summary,
          "/home/ryurko/asd_genome_wide_analysis/hmagma/data/output/adhd/adult/tidy_gene_pval_summary_4.csv")

# ------------------------------------------------------------------------------

# Remove the datasets and repeat for the second fold:
rm(list = ls())


flip_matrix <- function(x) {
  zero2 <- which(x == 0)
  two0 <- which(x == 2)
  x[zero2] <- 2
  x[two0] <- 0
  return(x)
}

# Repeat the following for the three different folds:
tidy_adult_brain_gene_snp_data <-
  read_csv("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/tidy_snp_gene_data/tidy_final_adult_gene_snp_pval_data.csv") %>%
  filter(!is.na(adhd_pval),
         server_fold %in% 5)
numeric_snp_ref_genotypes <- as.matrix(
  data.table::fread("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/reference_data/adult/fold5_ref_genotypes.csv"))
gene_ids <- unique(tidy_adult_brain_gene_snp_data$gene_id)

# Iterate through the groups to compute the group level test statistics:
tidy_group_pval_summary <-
  mclapply(1:length(gene_ids), mc.cores = 30,
           FUN = function(i) {

             group_x <- gene_ids[i]
             # Filter to the SNPs belonging to group X in the dataset:
             group_x_snp_rows <- tidy_adult_brain_gene_snp_data %>%
               filter(gene_id == group_x)

             # If only 1 SNP then just make a dataset with
             # that SNP's single p-value and test statistics
             # just for reference:
             if (nrow(group_x_snp_rows) == 1) {

               group_results <- group_x_snp_rows %>%
                 dplyr::select(gene_id,
                               adhd_pval) %>%
                 mutate(fisher_stat = -2 * log(adhd_pval),
                        fisher_pval = adhd_pval,
                        n_snps = 1) %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             } else {

               # Pull the SNPs:
               group_x_snps <- pull(group_x_snp_rows, snp_id)

               # Only look at these SNPs' genotype data - with the genotype
               # columns in the correct order:
               group_geno_ref_matrix <- numeric_snp_ref_genotypes[, group_x_snps]
               # Flip the matrix:
               group_geno_ref_matrix <- flip_matrix(group_geno_ref_matrix)

               # Compute the group's correlation matrix:
               group_cor_matrix <- compute_cor_matrix(group_geno_ref_matrix)

               # Compute the simulated truth statistics with this correlation matrix
               # using 300000 simulations to start:
               group_sim_truth <-
                 generate_sim_truth_data(group_cor_matrix,
                                         120000, 10, .Machine$double.eps)

               # The first column is the VEGAS data, second column is
               # the Fisher data - use these along with the observed
               # test statistics to compute empirical p-values and
               # return:
               group_results <- group_x_snp_rows %>%
                 group_by(gene_id) %>%
                 summarize(fisher_stat = -2 * sum(log(adhd_pval)),
                           n_snps = n()) %>%
                 ungroup() %>%
                 mutate(fisher_pval = mean(as.numeric(group_sim_truth[,2] > fisher_stat)))  %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             }

             group_results
           }) %>%
  bind_rows()
write_csv(tidy_group_pval_summary,
          "/home/ryurko/asd_genome_wide_analysis/hmagma/data/output/adhd/adult/tidy_gene_pval_summary_5.csv")

# ------------------------------------------------------------------------------

# Remove the datasets and repeat for the third fold:
rm(list = ls())


flip_matrix <- function(x) {
  zero2 <- which(x == 0)
  two0 <- which(x == 2)
  x[zero2] <- 2
  x[two0] <- 0
  return(x)
}

# Repeat the following for the three different folds:
tidy_adult_brain_gene_snp_data <-
  read_csv("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/tidy_snp_gene_data/tidy_final_adult_gene_snp_pval_data.csv") %>%
  filter(!is.na(adhd_pval),
         server_fold %in% 6)
numeric_snp_ref_genotypes <- as.matrix(
  data.table::fread("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/reference_data/adult/fold6_ref_genotypes.csv"))
gene_ids <- unique(tidy_adult_brain_gene_snp_data$gene_id)

# Iterate through the groups to compute the group level test statistics:
tidy_group_pval_summary <-
  mclapply(1:length(gene_ids), mc.cores = 30,
           FUN = function(i) {

             group_x <- gene_ids[i]
             # Filter to the SNPs belonging to group X in the dataset:
             group_x_snp_rows <- tidy_adult_brain_gene_snp_data %>%
               filter(gene_id == group_x)

             # If only 1 SNP then just make a dataset with
             # that SNP's single p-value and test statistics
             # just for reference:
             if (nrow(group_x_snp_rows) == 1) {

               group_results <- group_x_snp_rows %>%
                 dplyr::select(gene_id,
                               adhd_pval) %>%
                 mutate(fisher_stat = -2 * log(adhd_pval),
                        fisher_pval = adhd_pval,
                        n_snps = 1) %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             } else {

               # Pull the SNPs:
               group_x_snps <- pull(group_x_snp_rows, snp_id)

               # Only look at these SNPs' genotype data - with the genotype
               # columns in the correct order:
               group_geno_ref_matrix <- numeric_snp_ref_genotypes[, group_x_snps]
               # Flip the matrix:
               group_geno_ref_matrix <- flip_matrix(group_geno_ref_matrix)

               # Compute the group's correlation matrix:
               group_cor_matrix <- compute_cor_matrix(group_geno_ref_matrix)

               # Compute the simulated truth statistics with this correlation matrix
               # using 300000 simulations to start:
               group_sim_truth <-
                 generate_sim_truth_data(group_cor_matrix,
                                         120000, 10, .Machine$double.eps)

               # The first column is the VEGAS data, second column is
               # the Fisher data - use these along with the observed
               # test statistics to compute empirical p-values and
               # return:
               group_results <- group_x_snp_rows %>%
                 group_by(gene_id) %>%
                 summarize(fisher_stat = -2 * sum(log(adhd_pval)),
                           n_snps = n()) %>%
                 ungroup() %>%
                 mutate(fisher_pval = mean(as.numeric(group_sim_truth[,2] > fisher_stat)))  %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             }

             group_results
           }) %>%
  bind_rows()
write_csv(tidy_group_pval_summary,
          "/home/ryurko/asd_genome_wide_analysis/hmagma/data/output/adhd/adult/tidy_gene_pval_summary_6.csv")

# ------------------------------------------------------------------------------

rm(list = ls())


flip_matrix <- function(x) {
  zero2 <- which(x == 0)
  two0 <- which(x == 2)
  x[zero2] <- 2
  x[two0] <- 0
  return(x)
}

# Repeat the following for the three different folds:
tidy_adult_brain_gene_snp_data <-
  read_csv("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/tidy_snp_gene_data/tidy_final_adult_gene_snp_pval_data.csv") %>%
  filter(!is.na(adhd_pval),
         server_fold %in% 7)
numeric_snp_ref_genotypes <- as.matrix(
  data.table::fread("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/reference_data/adult/fold7_ref_genotypes.csv"))
gene_ids <- unique(tidy_adult_brain_gene_snp_data$gene_id)

# Iterate through the groups to compute the group level test statistics:
tidy_group_pval_summary <-
  mclapply(1:length(gene_ids), mc.cores = 30,
           FUN = function(i) {

             group_x <- gene_ids[i]
             # Filter to the SNPs belonging to group X in the dataset:
             group_x_snp_rows <- tidy_adult_brain_gene_snp_data %>%
               filter(gene_id == group_x)

             # If only 1 SNP then just make a dataset with
             # that SNP's single p-value and test statistics
             # just for reference:
             if (nrow(group_x_snp_rows) == 1) {

               group_results <- group_x_snp_rows %>%
                 dplyr::select(gene_id,
                               adhd_pval) %>%
                 mutate(fisher_stat = -2 * log(adhd_pval),
                        fisher_pval = adhd_pval,
                        n_snps = 1) %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             } else {

               # Pull the SNPs:
               group_x_snps <- pull(group_x_snp_rows, snp_id)

               # Only look at these SNPs' genotype data - with the genotype
               # columns in the correct order:
               group_geno_ref_matrix <- numeric_snp_ref_genotypes[, group_x_snps]
               # Flip the matrix:
               group_geno_ref_matrix <- flip_matrix(group_geno_ref_matrix)

               # Compute the group's correlation matrix:
               group_cor_matrix <- compute_cor_matrix(group_geno_ref_matrix)

               # Compute the simulated truth statistics with this correlation matrix
               # using 300000 simulations to start:
               group_sim_truth <-
                 generate_sim_truth_data(group_cor_matrix,
                                         120000, 10, .Machine$double.eps)

               # The first column is the VEGAS data, second column is
               # the Fisher data - use these along with the observed
               # test statistics to compute empirical p-values and
               # return:
               group_results <- group_x_snp_rows %>%
                 group_by(gene_id) %>%
                 summarize(fisher_stat = -2 * sum(log(adhd_pval)),
                           n_snps = n()) %>%
                 ungroup() %>%
                 mutate(fisher_pval = mean(as.numeric(group_sim_truth[,2] > fisher_stat)))  %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             }

             group_results
           }) %>%
  bind_rows()
write_csv(tidy_group_pval_summary,
          "/home/ryurko/asd_genome_wide_analysis/hmagma/data/output/adhd/adult/tidy_gene_pval_summary_7.csv")


# ------------------------------------------------------------------------------

# Remove the datasets and repeat for the third fold:
rm(list = ls())


flip_matrix <- function(x) {
  zero2 <- which(x == 0)
  two0 <- which(x == 2)
  x[zero2] <- 2
  x[two0] <- 0
  return(x)
}

# Repeat the following for the three different folds:
tidy_adult_brain_gene_snp_data <-
  read_csv("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/tidy_snp_gene_data/tidy_final_adult_gene_snp_pval_data.csv") %>%
  filter(!is.na(adhd_pval),
         server_fold %in% 8)
numeric_snp_ref_genotypes <- as.matrix(
  data.table::fread("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/reference_data/adult/fold8_ref_genotypes.csv"))
gene_ids <- unique(tidy_adult_brain_gene_snp_data$gene_id)

# Iterate through the groups to compute the group level test statistics:
tidy_group_pval_summary <-
  mclapply(1:length(gene_ids), mc.cores = 30,
           FUN = function(i) {

             group_x <- gene_ids[i]
             # Filter to the SNPs belonging to group X in the dataset:
             group_x_snp_rows <- tidy_adult_brain_gene_snp_data %>%
               filter(gene_id == group_x)

             # If only 1 SNP then just make a dataset with
             # that SNP's single p-value and test statistics
             # just for reference:
             if (nrow(group_x_snp_rows) == 1) {

               group_results <- group_x_snp_rows %>%
                 dplyr::select(gene_id,
                               adhd_pval) %>%
                 mutate(fisher_stat = -2 * log(adhd_pval),
                        fisher_pval = adhd_pval,
                        n_snps = 1) %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             } else {

               # Pull the SNPs:
               group_x_snps <- pull(group_x_snp_rows, snp_id)

               # Only look at these SNPs' genotype data - with the genotype
               # columns in the correct order:
               group_geno_ref_matrix <- numeric_snp_ref_genotypes[, group_x_snps]
               # Flip the matrix:
               group_geno_ref_matrix <- flip_matrix(group_geno_ref_matrix)

               # Compute the group's correlation matrix:
               group_cor_matrix <- compute_cor_matrix(group_geno_ref_matrix)

               # Compute the simulated truth statistics with this correlation matrix
               # using 300000 simulations to start:
               group_sim_truth <-
                 generate_sim_truth_data(group_cor_matrix,
                                         120000, 10, .Machine$double.eps)

               # The first column is the VEGAS data, second column is
               # the Fisher data - use these along with the observed
               # test statistics to compute empirical p-values and
               # return:
               group_results <- group_x_snp_rows %>%
                 group_by(gene_id) %>%
                 summarize(fisher_stat = -2 * sum(log(adhd_pval)),
                           n_snps = n()) %>%
                 ungroup() %>%
                 mutate(fisher_pval = mean(as.numeric(group_sim_truth[,2] > fisher_stat)))  %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             }

             group_results
           }) %>%
  bind_rows()
write_csv(tidy_group_pval_summary,
          "/home/ryurko/asd_genome_wide_analysis/hmagma/data/output/adhd/adult/tidy_gene_pval_summary_8.csv")


# ------------------------------------------------------------------------------
rm(list = ls())


flip_matrix <- function(x) {
  zero2 <- which(x == 0)
  two0 <- which(x == 2)
  x[zero2] <- 2
  x[two0] <- 0
  return(x)
}

# Repeat the following for the three different folds:
tidy_adult_brain_gene_snp_data <-
  read_csv("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/tidy_snp_gene_data/tidy_final_adult_gene_snp_pval_data.csv") %>%
  filter(!is.na(adhd_pval),
         server_fold %in% 9)
numeric_snp_ref_genotypes <- as.matrix(
  data.table::fread("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/reference_data/adult/fold9_ref_genotypes.csv"))
gene_ids <- unique(tidy_adult_brain_gene_snp_data$gene_id)


tidy_group_pval_summary <-
  mclapply(1:length(gene_ids), mc.cores = 30,
           FUN = function(i) {

             group_x <- gene_ids[i]
             # Filter to the SNPs belonging to group X in the dataset:
             group_x_snp_rows <- tidy_adult_brain_gene_snp_data %>%
               filter(gene_id == group_x)

             # If only 1 SNP then just make a dataset with
             # that SNP's single p-value and test statistics
             # just for reference:
             if (nrow(group_x_snp_rows) == 1) {

               group_results <- group_x_snp_rows %>%
                 dplyr::select(gene_id,
                               adhd_pval) %>%
                 mutate(fisher_stat = -2 * log(adhd_pval),
                        fisher_pval = adhd_pval,
                        n_snps = 1) %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             } else {

               # Pull the SNPs:
               group_x_snps <- pull(group_x_snp_rows, snp_id)

               # Only look at these SNPs' genotype data - with the genotype
               # columns in the correct order:
               group_geno_ref_matrix <- numeric_snp_ref_genotypes[, group_x_snps]
               # Flip the matrix:
               group_geno_ref_matrix <- flip_matrix(group_geno_ref_matrix)

               # Compute the group's correlation matrix:
               group_cor_matrix <- compute_cor_matrix(group_geno_ref_matrix)

               # Compute the simulated truth statistics with this correlation matrix
               # using 300000 simulations to start:
               group_sim_truth <-
                 generate_sim_truth_data(group_cor_matrix,
                                         120000, 10, .Machine$double.eps)

               # The first column is the VEGAS data, second column is
               # the Fisher data - use these along with the observed
               # test statistics to compute empirical p-values and
               # return:
               group_results <- group_x_snp_rows %>%
                 group_by(gene_id) %>%
                 summarize(fisher_stat = -2 * sum(log(adhd_pval)),
                           n_snps = n()) %>%
                 ungroup() %>%
                 mutate(fisher_pval = mean(as.numeric(group_sim_truth[,2] > fisher_stat)))  %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             }

             group_results
           }) %>%
  bind_rows()
write_csv(tidy_group_pval_summary,
          "/home/ryurko/asd_genome_wide_analysis/hmagma/data/output/adhd/adult/tidy_gene_pval_summary_9.csv")

# ------------------------------------------------------------------------------

# Remove the datasets and repeat for the second fold:
rm(list = ls())


flip_matrix <- function(x) {
  zero2 <- which(x == 0)
  two0 <- which(x == 2)
  x[zero2] <- 2
  x[two0] <- 0
  return(x)
}

# Repeat the following for the three different folds:
tidy_adult_brain_gene_snp_data <-
  read_csv("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/tidy_snp_gene_data/tidy_final_adult_gene_snp_pval_data.csv") %>%
  filter(!is.na(adhd_pval),
         server_fold %in% 10)
numeric_snp_ref_genotypes <- as.matrix(
  data.table::fread("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/reference_data/adult/fold10_ref_genotypes.csv"))
gene_ids <- unique(tidy_adult_brain_gene_snp_data$gene_id)


tidy_group_pval_summary <-
  mclapply(1:length(gene_ids), mc.cores = 30,
           FUN = function(i) {

             group_x <- gene_ids[i]
             # Filter to the SNPs belonging to group X in the dataset:
             group_x_snp_rows <- tidy_adult_brain_gene_snp_data %>%
               filter(gene_id == group_x)

             # If only 1 SNP then just make a dataset with
             # that SNP's single p-value and test statistics
             # just for reference:
             if (nrow(group_x_snp_rows) == 1) {

               group_results <- group_x_snp_rows %>%
                 dplyr::select(gene_id,
                               adhd_pval) %>%
                 mutate(fisher_stat = -2 * log(adhd_pval),
                        fisher_pval = adhd_pval,
                        n_snps = 1) %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             } else {

               # Pull the SNPs:
               group_x_snps <- pull(group_x_snp_rows, snp_id)

               # Only look at these SNPs' genotype data - with the genotype
               # columns in the correct order:
               group_geno_ref_matrix <- numeric_snp_ref_genotypes[, group_x_snps]
               # Flip the matrix:
               group_geno_ref_matrix <- flip_matrix(group_geno_ref_matrix)

               # Compute the group's correlation matrix:
               group_cor_matrix <- compute_cor_matrix(group_geno_ref_matrix)

               # Compute the simulated truth statistics with this correlation matrix
               # using 300000 simulations to start:
               group_sim_truth <-
                 generate_sim_truth_data(group_cor_matrix,
                                         120000, 10, .Machine$double.eps)

               # The first column is the VEGAS data, second column is
               # the Fisher data - use these along with the observed
               # test statistics to compute empirical p-values and
               # return:
               group_results <- group_x_snp_rows %>%
                 group_by(gene_id) %>%
                 summarize(fisher_stat = -2 * sum(log(adhd_pval)),
                           n_snps = n()) %>%
                 ungroup() %>%
                 mutate(fisher_pval = mean(as.numeric(group_sim_truth[,2] > fisher_stat)))  %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             }

             group_results
           }) %>%
  bind_rows()
write_csv(tidy_group_pval_summary,
          "/home/ryurko/asd_genome_wide_analysis/hmagma/data/output/adhd/adult/tidy_gene_pval_summary_10.csv")

# ------------------------------------------------------------------------------

# Remove the datasets and repeat for the second fold:
rm(list = ls())


flip_matrix <- function(x) {
  zero2 <- which(x == 0)
  two0 <- which(x == 2)
  x[zero2] <- 2
  x[two0] <- 0
  return(x)
}

# Repeat the following for the three different folds:
tidy_adult_brain_gene_snp_data <-
  read_csv("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/tidy_snp_gene_data/tidy_final_adult_gene_snp_pval_data.csv") %>%
  filter(!is.na(adhd_pval),
         server_fold %in% 11)
numeric_snp_ref_genotypes <- as.matrix(
  data.table::fread("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/reference_data/adult/fold11_ref_genotypes.csv"))
gene_ids <- unique(tidy_adult_brain_gene_snp_data$gene_id)


tidy_group_pval_summary <-
  mclapply(1:length(gene_ids), mc.cores = 30,
           FUN = function(i) {

             group_x <- gene_ids[i]
             # Filter to the SNPs belonging to group X in the dataset:
             group_x_snp_rows <- tidy_adult_brain_gene_snp_data %>%
               filter(gene_id == group_x)

             # If only 1 SNP then just make a dataset with
             # that SNP's single p-value and test statistics
             # just for reference:
             if (nrow(group_x_snp_rows) == 1) {

               group_results <- group_x_snp_rows %>%
                 dplyr::select(gene_id,
                               adhd_pval) %>%
                 mutate(fisher_stat = -2 * log(adhd_pval),
                        fisher_pval = adhd_pval,
                        n_snps = 1) %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             } else {

               # Pull the SNPs:
               group_x_snps <- pull(group_x_snp_rows, snp_id)

               # Only look at these SNPs' genotype data - with the genotype
               # columns in the correct order:
               group_geno_ref_matrix <- numeric_snp_ref_genotypes[, group_x_snps]
               # Flip the matrix:
               group_geno_ref_matrix <- flip_matrix(group_geno_ref_matrix)

               # Compute the group's correlation matrix:
               group_cor_matrix <- compute_cor_matrix(group_geno_ref_matrix)

               # Compute the simulated truth statistics with this correlation matrix
               # using 300000 simulations to start:
               group_sim_truth <-
                 generate_sim_truth_data(group_cor_matrix,
                                         120000, 10, .Machine$double.eps)

               # The first column is the VEGAS data, second column is
               # the Fisher data - use these along with the observed
               # test statistics to compute empirical p-values and
               # return:
               group_results <- group_x_snp_rows %>%
                 group_by(gene_id) %>%
                 summarize(fisher_stat = -2 * sum(log(adhd_pval)),
                           n_snps = n()) %>%
                 ungroup() %>%
                 mutate(fisher_pval = mean(as.numeric(group_sim_truth[,2] > fisher_stat)))  %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             }

             group_results
           }) %>%
  bind_rows()
write_csv(tidy_group_pval_summary,
          "/home/ryurko/asd_genome_wide_analysis/hmagma/data/output/adhd/adult/tidy_gene_pval_summary_11.csv")

# ------------------------------------------------------------------------------

# Remove the datasets and repeat for the second fold:
rm(list = ls())



flip_matrix <- function(x) {
  zero2 <- which(x == 0)
  two0 <- which(x == 2)
  x[zero2] <- 2
  x[two0] <- 0
  return(x)
}

# Repeat the following for the three different folds:
tidy_adult_brain_gene_snp_data <-
  read_csv("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/tidy_snp_gene_data/tidy_final_adult_gene_snp_pval_data.csv") %>%
  filter(!is.na(adhd_pval),
         server_fold %in% 12)
numeric_snp_ref_genotypes <- as.matrix(
  data.table::fread("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/reference_data/adult/fold12_ref_genotypes.csv"))
gene_ids <- unique(tidy_adult_brain_gene_snp_data$gene_id)

# Iterate through the groups to compute the group level test statistics:
tidy_group_pval_summary <-
  mclapply(1:length(gene_ids), mc.cores = 30,
           FUN = function(i) {

             group_x <- gene_ids[i]
             # Filter to the SNPs belonging to group X in the dataset:
             group_x_snp_rows <- tidy_adult_brain_gene_snp_data %>%
               filter(gene_id == group_x)

             # If only 1 SNP then just make a dataset with
             # that SNP's single p-value and test statistics
             # just for reference:
             if (nrow(group_x_snp_rows) == 1) {

               group_results <- group_x_snp_rows %>%
                 dplyr::select(gene_id,
                               adhd_pval) %>%
                 mutate(fisher_stat = -2 * log(adhd_pval),
                        fisher_pval = adhd_pval,
                        n_snps = 1) %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             } else {

               # Pull the SNPs:
               group_x_snps <- pull(group_x_snp_rows, snp_id)

               # Only look at these SNPs' genotype data - with the genotype
               # columns in the correct order:
               group_geno_ref_matrix <- numeric_snp_ref_genotypes[, group_x_snps]
               # Flip the matrix:
               group_geno_ref_matrix <- flip_matrix(group_geno_ref_matrix)

               # Compute the group's correlation matrix:
               group_cor_matrix <- compute_cor_matrix(group_geno_ref_matrix)

               # Compute the simulated truth statistics with this correlation matrix
               # using 300000 simulations to start:
               group_sim_truth <-
                 generate_sim_truth_data(group_cor_matrix,
                                         120000, 10, .Machine$double.eps)

               # The first column is the VEGAS data, second column is
               # the Fisher data - use these along with the observed
               # test statistics to compute empirical p-values and
               # return:
               group_results <- group_x_snp_rows %>%
                 group_by(gene_id) %>%
                 summarize(fisher_stat = -2 * sum(log(adhd_pval)),
                           n_snps = n()) %>%
                 ungroup() %>%
                 mutate(fisher_pval = mean(as.numeric(group_sim_truth[,2] > fisher_stat)))  %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             }

             group_results
           }) %>%
  bind_rows()
write_csv(tidy_group_pval_summary,
          "/home/ryurko/asd_genome_wide_analysis/hmagma/data/output/adhd/adult/tidy_gene_pval_summary_12.csv")

# ------------------------------------------------------------------------------

# Remove the datasets and repeat for the third fold:
rm(list = ls())


flip_matrix <- function(x) {
  zero2 <- which(x == 0)
  two0 <- which(x == 2)
  x[zero2] <- 2
  x[two0] <- 0
  return(x)
}

# Repeat the following for the three different folds:
tidy_adult_brain_gene_snp_data <-
  read_csv("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/tidy_snp_gene_data/tidy_final_adult_gene_snp_pval_data.csv") %>%
  filter(!is.na(adhd_pval),
         server_fold %in% 13)
numeric_snp_ref_genotypes <- as.matrix(
  data.table::fread("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/reference_data/adult/fold13_ref_genotypes.csv"))
gene_ids <- unique(tidy_adult_brain_gene_snp_data$gene_id)

# Iterate through the groups to compute the group level test statistics:
tidy_group_pval_summary <-
  mclapply(1:length(gene_ids), mc.cores = 30,
           FUN = function(i) {

             group_x <- gene_ids[i]
             # Filter to the SNPs belonging to group X in the dataset:
             group_x_snp_rows <- tidy_adult_brain_gene_snp_data %>%
               filter(gene_id == group_x)

             # If only 1 SNP then just make a dataset with
             # that SNP's single p-value and test statistics
             # just for reference:
             if (nrow(group_x_snp_rows) == 1) {

               group_results <- group_x_snp_rows %>%
                 dplyr::select(gene_id,
                               adhd_pval) %>%
                 mutate(fisher_stat = -2 * log(adhd_pval),
                        fisher_pval = adhd_pval,
                        n_snps = 1) %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             } else {

               # Pull the SNPs:
               group_x_snps <- pull(group_x_snp_rows, snp_id)

               # Only look at these SNPs' genotype data - with the genotype
               # columns in the correct order:
               group_geno_ref_matrix <- numeric_snp_ref_genotypes[, group_x_snps]
               # Flip the matrix:
               group_geno_ref_matrix <- flip_matrix(group_geno_ref_matrix)

               # Compute the group's correlation matrix:
               group_cor_matrix <- compute_cor_matrix(group_geno_ref_matrix)

               # Compute the simulated truth statistics with this correlation matrix
               # using 300000 simulations to start:
               group_sim_truth <-
                 generate_sim_truth_data(group_cor_matrix,
                                         120000, 10, .Machine$double.eps)

               # The first column is the VEGAS data, second column is
               # the Fisher data - use these along with the observed
               # test statistics to compute empirical p-values and
               # return:
               group_results <- group_x_snp_rows %>%
                 group_by(gene_id) %>%
                 summarize(fisher_stat = -2 * sum(log(adhd_pval)),
                           n_snps = n()) %>%
                 ungroup() %>%
                 mutate(fisher_pval = mean(as.numeric(group_sim_truth[,2] > fisher_stat)))  %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             }

             group_results
           }) %>%
  bind_rows()
write_csv(tidy_group_pval_summary,
          "/home/ryurko/asd_genome_wide_analysis/hmagma/data/output/adhd/adult/tidy_gene_pval_summary_13.csv")

# ------------------------------------------------------------------------------

# Remove the datasets and repeat for the third fold:
rm(list = ls())


flip_matrix <- function(x) {
  zero2 <- which(x == 0)
  two0 <- which(x == 2)
  x[zero2] <- 2
  x[two0] <- 0
  return(x)
}

# Repeat the following for the three different folds:
tidy_adult_brain_gene_snp_data <-
  read_csv("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/tidy_snp_gene_data/tidy_final_adult_gene_snp_pval_data.csv") %>%
  filter(!is.na(adhd_pval),
         server_fold %in% 14)
numeric_snp_ref_genotypes <- as.matrix(
  data.table::fread("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/reference_data/adult/fold14_ref_genotypes.csv"))
gene_ids <- unique(tidy_adult_brain_gene_snp_data$gene_id)

# Iterate through the groups to compute the group level test statistics:
tidy_group_pval_summary <-
  mclapply(1:length(gene_ids), mc.cores = 30,
           FUN = function(i) {

             group_x <- gene_ids[i]
             # Filter to the SNPs belonging to group X in the dataset:
             group_x_snp_rows <- tidy_adult_brain_gene_snp_data %>%
               filter(gene_id == group_x)

             # If only 1 SNP then just make a dataset with
             # that SNP's single p-value and test statistics
             # just for reference:
             if (nrow(group_x_snp_rows) == 1) {

               group_results <- group_x_snp_rows %>%
                 dplyr::select(gene_id,
                               adhd_pval) %>%
                 mutate(fisher_stat = -2 * log(adhd_pval),
                        fisher_pval = adhd_pval,
                        n_snps = 1) %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             } else {

               # Pull the SNPs:
               group_x_snps <- pull(group_x_snp_rows, snp_id)

               # Only look at these SNPs' genotype data - with the genotype
               # columns in the correct order:
               group_geno_ref_matrix <- numeric_snp_ref_genotypes[, group_x_snps]
               # Flip the matrix:
               group_geno_ref_matrix <- flip_matrix(group_geno_ref_matrix)

               # Compute the group's correlation matrix:
               group_cor_matrix <- compute_cor_matrix(group_geno_ref_matrix)

               # Compute the simulated truth statistics with this correlation matrix
               # using 300000 simulations to start:
               group_sim_truth <-
                 generate_sim_truth_data(group_cor_matrix,
                                         120000, 10, .Machine$double.eps)

               # The first column is the VEGAS data, second column is
               # the Fisher data - use these along with the observed
               # test statistics to compute empirical p-values and
               # return:
               group_results <- group_x_snp_rows %>%
                 group_by(gene_id) %>%
                 summarize(fisher_stat = -2 * sum(log(adhd_pval)),
                           n_snps = n()) %>%
                 ungroup() %>%
                 mutate(fisher_pval = mean(as.numeric(group_sim_truth[,2] > fisher_stat)))  %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             }

             group_results
           }) %>%
  bind_rows()
write_csv(tidy_group_pval_summary,
          "/home/ryurko/asd_genome_wide_analysis/hmagma/data/output/adhd/adult/tidy_gene_pval_summary_14.csv")




# Remove the datasets and repeat for the third fold:
rm(list = ls())


flip_matrix <- function(x) {
  zero2 <- which(x == 0)
  two0 <- which(x == 2)
  x[zero2] <- 2
  x[two0] <- 0
  return(x)
}

# Repeat the following for the three different folds:
tidy_adult_brain_gene_snp_data <-
  read_csv("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/tidy_snp_gene_data/tidy_final_adult_gene_snp_pval_data.csv") %>%
  filter(!is.na(adhd_pval),
         server_fold %in% 15)
numeric_snp_ref_genotypes <- as.matrix(
  data.table::fread("/home/ryurko/asd_genome_wide_analysis/hmagma/data/input/reference_data/adult/fold15_ref_genotypes.csv"))
gene_ids <- unique(tidy_adult_brain_gene_snp_data$gene_id)

# Iterate through the groups to compute the group level test statistics:
tidy_group_pval_summary <-
  mclapply(1:length(gene_ids), mc.cores = 30,
           FUN = function(i) {

             group_x <- gene_ids[i]
             # Filter to the SNPs belonging to group X in the dataset:
             group_x_snp_rows <- tidy_adult_brain_gene_snp_data %>%
               filter(gene_id == group_x)

             # If only 1 SNP then just make a dataset with
             # that SNP's single p-value and test statistics
             # just for reference:
             if (nrow(group_x_snp_rows) == 1) {

               group_results <- group_x_snp_rows %>%
                 dplyr::select(gene_id,
                               adhd_pval) %>%
                 mutate(fisher_stat = -2 * log(adhd_pval),
                        fisher_pval = adhd_pval,
                        n_snps = 1) %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             } else {

               # Pull the SNPs:
               group_x_snps <- pull(group_x_snp_rows, snp_id)

               # Only look at these SNPs' genotype data - with the genotype
               # columns in the correct order:
               group_geno_ref_matrix <- numeric_snp_ref_genotypes[, group_x_snps]
               # Flip the matrix:
               group_geno_ref_matrix <- flip_matrix(group_geno_ref_matrix)

               # Compute the group's correlation matrix:
               group_cor_matrix <- compute_cor_matrix(group_geno_ref_matrix)

               # Compute the simulated truth statistics with this correlation matrix
               # using 300000 simulations to start:
               group_sim_truth <-
                 generate_sim_truth_data(group_cor_matrix,
                                         120000, 10, .Machine$double.eps)

               # The first column is the VEGAS data, second column is
               # the Fisher data - use these along with the observed
               # test statistics to compute empirical p-values and
               # return:
               group_results <- group_x_snp_rows %>%
                 group_by(gene_id) %>%
                 summarize(fisher_stat = -2 * sum(log(adhd_pval)),
                           n_snps = n()) %>%
                 ungroup() %>%
                 mutate(fisher_pval = mean(as.numeric(group_sim_truth[,2] > fisher_stat)))  %>%
                 dplyr::select(gene_id,
                               fisher_stat,
                               fisher_pval,
                               n_snps)

             }

             group_results
           }) %>%
  bind_rows()
write_csv(tidy_group_pval_summary,
          "/home/ryurko/asd_genome_wide_analysis/hmagma/data/output/adhd/adult/tidy_gene_pval_summary_15.csv")




