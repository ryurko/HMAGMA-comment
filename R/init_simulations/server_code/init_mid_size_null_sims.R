# SERVER VERSION

library(backports, lib.loc = "/home/ryurko/Rpackages")
library(vctrs, lib.loc = "/home/ryurko/Rpackages")
library(dplyr, lib.loc = "/home/ryurko/Rpackages")
library(magrittr, lib.loc = "/home/ryurko/Rpackages")
library(data.table, lib.loc = "/home/ryurko/Rpackages")
library(snpcombineR, lib.loc = "/home/ryurko/Rpackages")
library(readr, lib.loc = "/home/ryurko/Rpackages")
library(stringr, lib.loc = "/home/ryurko/Rpackages")
library(parallel, lib.loc = "/home/ryurko/Rpackages")

# Gene ID: 169355

# Load the resampled data:
resampled_data <- as.matrix(
  data.table::fread(
    "/home/ryurko/asd_gene_analysis/combine_snps_sims/resampled_genotypes/ex_gene169355_chr8_nsnps101_genotype.csv"))

# Compute the gene correlation matrix:
gene_cor_matrix <- compute_cor_matrix(resampled_data)

# Compute the simulated truth statistics 
gene_sim_truth_matrix <- generate_sim_truth_data(gene_cor_matrix,
                                                 120000, 10, .Machine$double.eps)


sim_list <-
  mclapply(1:25, mc.cores = 25,
           FUN = function(i) {
             
             null_sims_results <- sim_gene_cor_gwas_gene(resampled_data,
                                                         gene_cor_matrix,
                                                         gene_sim_truth_matrix,
                                                         120000,
                                                         FALSE,
                                                         # These parameters are ignored
                                                         c(1:5), 1.5,
                                                         # .5 case rate
                                                         0.5,
                                                         # Use  core
                                                         1)
             write_csv(as.data.frame(null_sims_results),
                       paste0(
                         "/home/ryurko/asd_gene_analysis/combine_snps_sims/update_resampled_sims/more_null/gene169355/null_sim_",
                         i, ".csv"))
             
             return(i)
             
           })

rm(list = ls())



# Gene ID: 80243

# Load the resampled data:
resampled_data <- as.matrix(
  data.table::fread(
    "/home/ryurko/asd_gene_analysis/combine_snps_sims/resampled_genotypes/ex_gene80243_chr8_nsnps218_genotype.csv"))

# Compute the gene correlation matrix:
gene_cor_matrix <- compute_cor_matrix(resampled_data)

# Compute the simulated truth statistics 
gene_sim_truth_matrix <- generate_sim_truth_data(gene_cor_matrix,
                                                 120000, 10, .Machine$double.eps)


sim_list <-
  mclapply(1:25, mc.cores = 25,
           FUN = function(i) {
             
             null_sims_results <- sim_gene_cor_gwas_gene(resampled_data,
                                                         gene_cor_matrix,
                                                         gene_sim_truth_matrix,
                                                         120000,
                                                         FALSE,
                                                         # These parameters are ignored
                                                         c(1:5), 1.5,
                                                         # .5 case rate
                                                         0.5,
                                                         # Use  core
                                                         1)
             write_csv(as.data.frame(null_sims_results),
                       paste0(
                         "/home/ryurko/asd_gene_analysis/combine_snps_sims/update_resampled_sims/more_null/gene80243/null_sim_",
                         i, ".csv"))
             
             return(i)
             
           })



rm(list = ls())



# Gene ID: 26064

# Load the resampled data:
resampled_data <- as.matrix(
  data.table::fread(
    "/home/ryurko/asd_gene_analysis/combine_snps_sims/resampled_genotypes/ex_gene80243_chr8_nsnps218_genotype.csv"))

# Compute the gene correlation matrix:
gene_cor_matrix <- compute_cor_matrix(resampled_data)

# Compute the simulated truth statistics 
gene_sim_truth_matrix <- generate_sim_truth_data(gene_cor_matrix,
                                                 120000, 10, .Machine$double.eps)


sim_list <-
  mclapply(1:25, mc.cores = 25,
           FUN = function(i) {
             
             null_sims_results <- sim_gene_cor_gwas_gene(resampled_data,
                                                         gene_cor_matrix,
                                                         gene_sim_truth_matrix,
                                                         120000,
                                                         FALSE,
                                                         # These parameters are ignored
                                                         c(1:5), 1.5,
                                                         # .5 case rate
                                                         0.5,
                                                         # Use  core
                                                         1)
             write_csv(as.data.frame(null_sims_results),
                       paste0(
                         "/home/ryurko/asd_gene_analysis/combine_snps_sims/update_resampled_sims/more_null/gene26064/null_sim_",
                         i, ".csv"))
             
             return(i)
             
           })



