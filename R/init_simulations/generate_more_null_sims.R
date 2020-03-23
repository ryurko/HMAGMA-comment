# PURPOSE: Generate more null simulations (using the already sampled genotype
#          matrix) for the multiple testing simulations. This script can be called
#          twice, just changing the files for saving.
# NOTE: This script was distributed across a server - you may need to adjust
#       the number of cores/threads depending on the machine/server considered

# Access necessary packages:
library(tidyverse)
library(snpcombineR)
library(future)
library(furrr)

# ------------------------------------------------------------------------------

# Create the vector of gene files:
ex_gene_pruned_file_paths <- list.files(
  "data/reference_data/pruned_genotype",
  full.names = TRUE) 

# Just select the randomly selected example genes:
ex_gene_pruned_file_paths <- 
  str_subset(
    ex_gene_pruned_file_paths,
    "(gene51306_)|(gene3848_)|(gene659_)|(gene80243_)|(gene169355_)|(gene26064_)|(gene6263_)|(gene123624_)|(gene11122_)")

# Create a data frame of the gene info:
ex_gene_info_df <- map_dfr(ex_gene_pruned_file_paths,
                           function(gene_file) {
                             # Get the file info summary:
                             file_info <- 
                               str_remove(gene_file,
                                          "data/reference_data/pruned_genotype") %>%
                               str_remove("_genotype\\.csv")
                             
                             # Split this:
                             file_split_info <- unlist(str_split(file_info, 
                                                                 "_"))
                             
                             # Extract the gene, chromosome and number of SNPs:
                             data.frame(
                               "gene_name" = str_remove(file_split_info[2],
                                                        "gene"),
                               "gene_chr" = as.numeric(str_remove(file_split_info[3],
                                                                  "chr")),
                               "gene_n_snps" = as.numeric(str_remove(file_split_info[4],
                                                                     "nsnps")),
                               "file" = gene_file
                             )
                           }) %>%
  arrange(gene_n_snps)

plan(tweak(multiprocess, workers = 9))
results <- future_map_dfr(1:nrow(ex_gene_info_df),
                          function(gene_i) {
                            ex_gene_file_path <- ex_gene_info_df$file[gene_i]
                            ex_gene_file_info <- str_remove(ex_gene_file_path,
                                                            "data/reference_data/pruned_genotype") %>%
                              str_remove("_genotype\\.csv")
                            # Split this:
                            file_split_info <- unlist(str_split(ex_gene_file_info,
                                                                "_"))
                            
                            # Extract the gene, chromosome and number of SNPs:
                            gene_name <- ex_gene_info_df$gene_name[gene_i]
                            gene_chr <- ex_gene_info_df$gene_chr[gene_i]
                            gene_n_snps <- ex_gene_info_df$gene_n_snps[gene_i]
                            
                            # Load the resampled data:
                            resampled_data <- as.matrix(
                              read_csv(
                                paste0(
                                  "data/reference_data/resampled_genotypes/",
                                  ex_gene_file_info, "_genotype.csv")))
                            
                            # Compute the gene correlation matrix:
                            gene_cor_matrix <- compute_cor_matrix(resampled_data)
                            
                            # Compute the simulated truth statistics
                            gene_sim_truth_matrix <- generate_sim_truth_data(gene_cor_matrix,
                                                                             10, 10000, .Machine$double.eps)
                            
                            null_sims_results <- sim_gene_cor_gwas_gene(resampled_data,
                                                                        gene_cor_matrix,
                                                                        gene_sim_truth_matrix,
                                                                        450000,
                                                                        FALSE,
                                                                        # These parameters are ignored
                                                                        c(1:5), 1.5,
                                                                        # .5 case rate
                                                                        0.5,
                                                                        # Use 30 cores
                                                                        30)
                            write_csv(as.data.frame(null_sims_results),
                                      paste0("data/simulations/null/",
                                             ex_gene_file_info, 
                                             "_sims_2.csv"))
                                             #"_sims_3.csv"))
 
                            # Just return the return:
                            ex_gene_info_df[gene_i,]
                            
                          })
