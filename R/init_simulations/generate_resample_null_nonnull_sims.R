# PURPOSE: Generate null and non-null simulations for the randomly picked genes,
#          after initially re-sampling the datasets.
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

# Set up a dataframe with alternative settings:
odds_ratios <- seq(1.1, 2, by = .1)
causal_perc <- c(0, .1)  # will use pmax so that 0 corresponds to 1 causal SNP
alt_settings <- expand.grid(odds_ratios, causal_perc)
colnames(alt_settings) <- c("odds_ratios", "perc_causal_snps")

# Now in parallel resample the genotype matrices, and generate each gene's 
# null and non-null simulation results:
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

                            # Make a copy of the alt_settings:
                            gene_alt_settings <- alt_settings
                            # Now compute the n_snps taking the maximum of either
                            # 1 or the perc_causal_snps * gene_n_snps:
                            gene_alt_settings <- gene_alt_settings %>%
                               mutate(n_snps = pmax(1, ceiling(gene_n_snps * perc_causal_snps))) %>%
                               dplyr::select(-perc_causal_snps) %>%
                               distinct()

                            # First load the data:
                            example_gene_ref_genotypes <- suppressMessages(as.matrix(
                              readr::read_csv(ex_gene_file_path)))

                            # Next resample it 5000 times
                            resampled_data <- example_gene_ref_genotypes[sample(1:nrow(example_gene_ref_genotypes),
                                                                                size = 5000, replace = TRUE),]

                            # Save this resampled data:
                            write_csv(as.data.frame(resampled_data),
                                      paste0("data/reference_data/resampled_genotypes/",
                                             ex_gene_file_info, "_genotype.csv"))

                            # Compute the gene correlation matrix:
                            gene_cor_matrix <- compute_cor_matrix(resampled_data)

                            # Compute the simulated truth statistics with this correlation matrix:
                            gene_sim_truth_matrix <- generate_sim_truth_data(gene_cor_matrix,
                                                                             10, 10000, .Machine$double.eps)
                            
                            # Next generate the null sim results:
                            null_sims_results <- sim_gene_cor_gwas_gene(resampled_data,
                                                                        gene_cor_matrix,
                                                                        gene_sim_truth_matrix,
                                                                        100000, FALSE,
                                                                        # These parameters are ignored
                                                                        c(1:5), 1.5,
                                                                        # .5 case rate
                                                                        0.5,
                                                                        # Use 30 cores
                                                                        30)
                            write_csv(as.data.frame(null_sims_results),
                                      paste0("data/simulations/null/",
                                             ex_gene_file_info, 
                                             "_sims_1.csv"))

                            # What are the distinct choices for number of causal SNPs:
                            distinct_n_snps <- unique(gene_alt_settings$n_snps)

                            for (i in 1:length(distinct_n_snps)) {
                              n_causal_snps <- distinct_n_snps[i]

                              # Randomly pick the causal SNPs (subtract 1 for C++)
                              causal_snp_i <- sample(1:ncol(gene_cor_matrix), n_causal_snps) - 1

                              # Next loop through the signal strengths:
                              filter_gene_settings <- gene_alt_settings %>%
                                dplyr::filter(n_snps == n_causal_snps)

                              for (j in 1:nrow(filter_gene_settings)) {

                                # Odds ratio:
                                or_setting <- filter_gene_settings$odds_ratios[j]

                                # Now generate the non-null results:
                                alt_sims_results <- sim_gene_cor_gwas_gene(resampled_data,
                                                                           gene_cor_matrix,
                                                                           gene_sim_truth_matrix,
                                                                           100000,
                                                                           # non-null settings:
                                                                           TRUE,
                                                                           causal_snp_i,
                                                                           or_setting,
                                                                           # .5 case rate
                                                                           0.5,
                                                                           # Use 30 cores
                                                                           30)
                                
                                write_csv(as.data.frame(alt_sims_results),
                                          paste0("data/simulations/nonnull/",
                                                 ex_gene_file_info, "_nsnps",
                                                 n_causal_snps, "_or", or_setting * 10,
                                                 "_sims.csv"))


                              }


                            }
                            
                            # Just return the return:
                            ex_gene_info_df[gene_i,]

                          })

