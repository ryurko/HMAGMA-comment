# PURPOSE: Create supplementary figures for multiple testing results

library(tidyverse)
library(latex2exp)
library(cowplot)

# ------------------------------------------------------------------------------

null_sim_files <- list.files("data/simulations/null/",
                             full.names = TRUE) %>%
  str_remove("_sims_(1|2|3)\\.csv") %>%
  unique()

null_mt_pvalues <- map_dfr(null_sim_files,
                           function(sim_file) {
                             
                             # Get the file info summary:
                             file_info <- str_remove(sim_file,
                                                     "data/gene_test_sims/update_resampled_sims/null/")
                             
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
                             # Next extract just the columns with the p-values
                             # from the data:
                             gene_pvals <- suppressMessages(
                               suppressWarnings(
                                 dplyr::bind_rows(read_csv(paste0(sim_file, "_test_1.csv")),
                                                  read_csv(paste0(sim_file, "_test_2.csv")),
                                                  read_csv(paste0(sim_file, "_test_3.csv")),
                                                  read_csv(paste0(sim_file, "_test_4.csv")))))
                             gene_pvals <- gene_pvals[, c(4, 8, 12, 16, 18, 20)]
                             colnames(gene_pvals) <- 
                               c("magma_fudge", "magma_no_fudge",
                                 "browns_two_sided", "browns_one_sided",
                                 "fisher", "vegas")
                             
                             # Randomly assing 1000 to each set so in total there
                             # are 9000 genes in each set:
                             gene_pvals <- gene_pvals %>%
                               mutate(sim_index = 1:n(),
                                      mt_set = rep(1:(n() / 1000), 1000),
                                      mt_set = sample(mt_set, n(), replace = FALSE),
                                      gene = gene_name,
                                      chr = gene_chr,
                                      n_snps = gene_n_snps)
                             return(gene_pvals)
                           }) %>%
  dplyr::select(-gene, -chr, -n_snps, -sim_index) %>%
  gather(test_type, pval, -mt_set) %>%
  group_by(mt_set, test_type) %>%
  mutate(bf_pval = p.adjust(pval, method = "bonferroni"),
         bh_pval = p.adjust(pval, method = "BH"),
         bf_error = as.numeric(bf_pval <= 0.05),
         bh_error = as.numeric(bh_pval <= 0.05)) %>%
  summarize(n_genes = n(),
            obs_fwer = ifelse(any(bf_error == 1), 1, 0),
            obs_fdr = ifelse(any(bh_error == 1), 1, 0),
            obs_bf_v = sum(bf_error),
            obs_bh_v = sum(bh_error)) %>%
  ungroup()

# Create plot displaying the expected number of false positives by method
bf_fp_sim_plot <- null_mt_pvalues %>%
  group_by(test_type) %>%
  summarize(mean_v_bf = mean(obs_bf_v),
            n_sims = n(),
            se_v_bf = sd(obs_bf_v) / sqrt(n_sims)) %>%
  mutate(obs_v_lower = pmax(0, mean_v_bf - 2 * se_v_bf),
         obs_v_upper = mean_v_bf + 2 * se_v_bf,
         test_type = fct_relevel(test_type,
                                 "browns_one_sided", 
                                 "magma_no_fudge",
                                 "magma_fudge", 
                                 "browns_two_sided", 
                                 "fisher",
                                 "vegas"),
         test_type = fct_rev(fct_recode(test_type,
                                        `MAGMA: paper` = "browns_one_sided",
                                        `MAGMA: rho^2` = "magma_no_fudge",
                                        `MAGMA: code` = "magma_fudge",
                                        `Two-sided approx.` = "browns_two_sided",
                                        `MC Fisher` = "fisher",
                                        `VEGAS` = "vegas"))) %>%
  ggplot(aes(x = test_type, y = mean_v_bf)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkred") +
  geom_errorbar(aes(ymin = obs_v_lower, ymax = obs_v_upper)) +
  theme_bw() + 
  labs(x = "Method", y = "Number of false positives") +
  scale_x_discrete(labels = c(`MAGMA: paper` = "MAGMA: paper",
                              `MAGMA: rho^2` = parse(text = TeX('MAGMA: $\\rho^2$')),
                              `MAGMA: code` = "MAGMA: code",
                              `Two-sided approx.` = "Two-sided approx.",
                              `MC Fisher` = "MC Fisher",
                              `VEGAS` = "VEGAS")) +
  coord_flip() + 
  theme(axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12))

save_plot("figures/si/sf_bf_mt_fp_sims.pdf",
          bf_fp_sim_plot, ncol = 1, nrow = 1)
save_plot("nonpdf_figures/si/sf_bf_mt_fp_sims.jpeg",
          bf_fp_sim_plot, ncol = 1, nrow = 1)
