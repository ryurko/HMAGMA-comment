# PURPOSE: Create figure comparing power from non-null simulations

library(tidyverse)
library(latex2exp)
library(cowplot)

# ------------------------------------------------------------------------------

null_sim_files <- list.files("data/simulations/null/",
                             full.names = TRUE) %>%
  str_remove("_sims_(1|2|3)\\.csv") %>%
  unique()

# For each gene - need to go through and find the alpha that leads to control
# at the original target of 0.05
null_sim_adj_alpha_summary <- map_dfr(null_sim_files,
                                      function(sim_file) {
                                        
                                        # Get the file info summary:
                                        file_info <- str_remove(sim_file,
                                                                "data/simulations/null/")
                                        
                                        # Split this:
                                        file_split_info <- unlist(str_split(file_info, 
                                                                            "_"))
                                        
                                        # Extract the gene, chromosome and number of SNPs:
                                        gene_name <- str_remove(file_split_info[2],
                                                                "gene")
                                        
                                        gene_chr <- as.numeric(str_remove(file_split_info[3],
                                                                          "chr"))
                                        gene_n_snps <- as.numeric(str_remove(file_split_info[4],
                                                                             "nsnps"))
                                        
                                        print(paste0("gene ", gene_name,
                                                     ", CHR ", gene_chr,
                                                     ", ", gene_n_snps, " SNPs"))
                                        # Next extract just the columns with the p-values
                                        # from the data (load all null sim results)
                                        gene_pvals <- suppressMessages(
                                          suppressWarnings(
                                            dplyr::bind_rows(read_csv(paste0(sim_file, "_sims_1.csv")),
                                                             read_csv(paste0(sim_file, "_sims_2.csv")),
                                                             read_csv(paste0(sim_file, "_sims_3.csv")))))
                                        
                                        gene_pvals <- gene_pvals[, c(4, 8, 12, 16, 18, 20)]
                                        colnames(gene_pvals) <- 
                                          c("magma_fudge", "magma_no_fudge",
                                            "browns_two_sided", "browns_one_sided",
                                            "fisher", "vegas")
                                        
                                        # Next go through each column of the dataset,
                                        # corresponding to each of the different tests,
                                        # and find when .05 is within the t1 error rate interval:
                                        correct_alpha_results <- map_dfr(1:ncol(gene_pvals),
                                                                         function(gene_test_i) {
                                                                           test_type <- colnames(gene_pvals)[gene_test_i]
                                                                           print(paste0("Using method ", test_type))
                                                                           # Initialize the adjustment term:
                                                                           alpha_adj <- 0.0001
                                                                           
                                                                           # Contains target alpha?
                                                                           reach_target <- FALSE
                                                                           loop_i <- 0
                                                                           
                                                                           if (test_type %in% c("fisher", "vegas")) {
                                                                             target_alpha <- 0.05
                                                                           } else {
                                                                             # Loop until target alpha 0.05 reached:
                                                                             while(!reach_target) {
                                                                               
                                                                               target_alpha <- 0.05 - (loop_i * alpha_adj)
                                                                               #print(target_alpha)
                                                                               gene_target_results <- data.frame("pval" = gene_pvals[[gene_test_i]]) %>%
                                                                                 mutate(t1_error = ifelse(pval <= target_alpha, 1, 0)) %>%
                                                                                 summarize(t1_error_rate = mean(t1_error, na.rm = TRUE),
                                                                                           n_sims = n()) %>%
                                                                                 mutate(t1_error_se = sqrt((t1_error_rate * (1- t1_error_rate)) /
                                                                                                             n_sims),
                                                                                        t1_error_lower = ifelse(t1_error_rate - 2 * t1_error_se < 0,
                                                                                                                0,
                                                                                                                t1_error_rate - 2 * t1_error_se),
                                                                                        t1_error_upper = ifelse(t1_error_rate + 2 * t1_error_se > 1,
                                                                                                                1,
                                                                                                                t1_error_rate + 2 * t1_error_se))
                                                                               
                                                                               # Is alpha 0.05 within the interval?
                                                                               reach_target <- (gene_target_results$t1_error_lower <= 0.05) &
                                                                                 (gene_target_results$t1_error_upper >= 0.05)
                                                                               loop_i <- loop_i + 1
                                                                               
                                                                             }
                                                                           }
                                                                           
                                                                           # Return a dataset with the target alpha used and the method:
                                                                           data.frame("target_alpha" = target_alpha,
                                                                                      "method" = test_type)
                                                                         }) %>%
                                          mutate(gene = gene_name,
                                                 chr = gene_chr,
                                                 n_snps = gene_n_snps)
                                        return(correct_alpha_results)
                                      })


# ------------------------------------------------------------------------------

# Make a vector of the non-null sim files - with the suffix removed 
nonnull_sim_files <- list.files("data/simulations/nonnull",
                                full.names = TRUE) %>%
  str_remove("_sims\\.csv") %>%
  unique()

nonnull_sim_result_summary <- map_dfr(nonnull_sim_files,
                                      function(sim_file) {
                                        
                                        # Get the file info summary:
                                        file_info <- str_remove(sim_file,
                                                                "data/simulations/nonnull/")
                                        
                                        # Split this:
                                        file_split_info <- unlist(str_split(file_info, 
                                                                            "_"))
                                        
                                        # Extract the gene, chromosome and number of SNPs:
                                        gene_name <- str_remove(file_split_info[2],
                                                                "gene")
                                        gene_chr <- as.numeric(str_remove(file_split_info[3],
                                                                          "chr"))
                                        gene_n_snps <- as.numeric(str_remove(file_split_info[4],
                                                                             "nsnps"))
                                        n_cause <- as.numeric(str_remove(file_split_info[5],
                                                                         "nsnps"))
                                        odds_ratio <- as.numeric(str_remove(file_split_info[6],
                                                                            "or"))
                                        
                                        # Next extract just the columns with the p-values
                                        # from the data:
                                        gene_pvals <-  suppressMessages(
                                          suppressWarnings(
                                            dplyr::bind_rows(read_csv(paste0(sim_file, "_sims.csv")))))
                                        gene_pvals <- gene_pvals[, c(4, 8, 12, 16, 18, 20)]
                                        colnames(gene_pvals) <- 
                                          c("magma_fudge", "magma_no_fudge",
                                            "browns_two_sided", "browns_one_sided",
                                            "fisher", "vegas")
                                        
                                        gene_results <- gene_pvals %>%
                                          mutate(sim_index = 1:n()) %>%
                                          gather(method, pval, -sim_index) %>%
                                          mutate(gene = gene_name,
                                                 chr = gene_chr,
                                                 n_snps = gene_n_snps) %>%
                                          left_join(null_sim_adj_alpha_summary,
                                                    by = c("method", "gene",
                                                           "chr", "n_snps")) %>%
                                          mutate(tp = ifelse(pval <= target_alpha, 1, 0)) %>%
                                          group_by(method) %>%
                                          summarize(power = mean(tp, na.rm = TRUE),
                                                    n_sims = n()) %>%
                                          mutate(power_se = sqrt((power * (1- power)) /
                                                                   n_sims),
                                                 power_lower = ifelse(power - 2 * power_se < 0,
                                                                      0,
                                                                      power - 2 * power_se),
                                                 power_upper = ifelse(power + 2 * power_se > 1,
                                                                      1,
                                                                      power + 2 * power_se)) %>%
                                          mutate(gene = gene_name,
                                                 chr = gene_chr,
                                                 n_snps = gene_n_snps,
                                                 n_causal_snps = n_cause,
                                                 or = odds_ratio)
                                        return(gene_results)
                                      })

# Create the power plot result using only 1 causal SNP:
power_sim_plot <- nonnull_sim_result_summary %>%
  filter(n_causal_snps == 1) %>%
  mutate(method = fct_relevel(method,
                              "browns_one_sided", "magma_no_fudge",
                              "magma_fudge", "browns_two_sided",
                              "fisher", 
                              "vegas"),
         method = fct_recode(method,
                             `MAGMA: paper` = "browns_one_sided",
                             `MAGMA: rho^2` = "magma_no_fudge",
                             `MAGMA: code` = "magma_fudge",
                             `Two-sided approx.` = "browns_two_sided",
                             `MC Fisher` = "fisher",
                             `VEGAS` = "vegas"),
         gene_label = paste0(n_snps, " SNPs"),
         gene_label = fct_reorder(gene_label, n_snps)) %>%
  ggplot(aes(x = as.factor(or / 10), color = method)) + 
  geom_point(aes(y = power),
             position = position_dodge(width = .85)) +
  geom_errorbar(aes(ymin = power_lower, ymax = power_upper),
                position = position_dodge(width = .85)) +
  theme_bw() + 
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1),
                              nrow = 1))+
  labs(x = "Odds ratio", color = "Method", y = "Power") +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~gene_label, ncol = 3) +
  ggthemes::scale_color_colorblind(labels = c(`MAGMA: paper` = "MAGMA: paper",
                                              `MAGMA: rho^2` = parse(text = TeX('MAGMA: $\\rho^2$')),
                                              `MAGMA: code` = "MAGMA: code",
                                              `Two-sided approx.` = "Two-sided approx.",
                                              `MC Fisher` = "MC Fisher",
                                              `VEGAS` = "VEGAS")) +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 24),
        strip.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18))


# save_plot("figures/si/sf_power_results.pdf",
#           power_sim_plot, ncol = 3, nrow = 3)
# save_plot("nonpdf_figures/si/sf_power_results.jpeg",
#           power_sim_plot, ncol = 3, nrow = 3)

