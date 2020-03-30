# PURPOSE: Create the supplementary gene-set analysis figures

library(tidyverse)
library(latex2exp)
library(cowplot)

# ------------------------------------------------------------------------------

null_sim_files <- list.files("data/simulations/null/",
                             full.names = TRUE) %>%
  str_remove("_sims_(1|2|3)\\.csv") %>%
  unique()

null_gene_set_pvalues <- map_dfr(null_sim_files,
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
                                   # Next extract just the columns with the p-values
                                   # from the data:
                                   gene_pvals <- suppressMessages(
                                     suppressWarnings(
                                       dplyr::bind_rows(read_csv(paste0(sim_file, "_sims_1.csv")))))
                                   gene_pvals <- gene_pvals[, c(4, 8, 12, 16, 18, 20)]
                                   colnames(gene_pvals) <- 
                                     c("magma_fudge", "magma_no_fudge",
                                       "browns_two_sided", "browns_one_sided",
                                       "fisher", "vegas")
                                   
                                   gene_pvals <- gene_pvals %>%
                                     mutate(sim_index = 1:n(),
                                            gene_set = rep(1:(n() / 5), 5),
                                            gene_set = sample(gene_set, n(), replace = FALSE),
                                            gene = gene_name,
                                            chr = gene_chr,
                                            n_snps = gene_n_snps)
                                   return(gene_pvals)
                                 })
null_gene_set_sim_results <- null_gene_set_pvalues %>%
  dplyr::select(-sim_index, -gene, -chr, -n_snps) %>%
  gather(test_type, pval, -gene_set) %>%
  mutate(probit_stat = qnorm(1 - pval),
         probit_stat = ifelse(is.infinite(probit_stat) & (pval < 0.5),
                              qnorm(1 - .Machine$double.eps),
                              probit_stat),
         probit_stat = ifelse(is.infinite(probit_stat) & (pval > 0.5),
                              qnorm(1 - (1 - .Machine$double.eps)),
                              probit_stat),
         neg_log_pval = -log(pval),
         neg_log_pval = ifelse(is.infinite(neg_log_pval) & (pval < 0.5),
                               -log(.Machine$double.eps),
                               neg_log_pval),
         neg_log_pval = ifelse(is.infinite(neg_log_pval) & (pval > 0.5),
                               -log(1 - .Machine$double.eps),
                               neg_log_pval)) %>%
  group_by(gene_set, test_type) %>%
  summarize(fisher_test_stat = 2 * sum(neg_log_pval),
            fisher_df = 2 * n(),
            stouffer_test_stat = sum(probit_stat) / sqrt(n()),
            magma_gs_test_stat = mean((probit_stat)),
            magma_gs_se = sd((probit_stat)) / sqrt(n())) %>%
  mutate(fisher_gs_pval = pchisq(q = fisher_test_stat,
                                 df = fisher_df, lower.tail = FALSE),
         stouffer_gs_pval = pnorm(stouffer_test_stat, lower.tail = FALSE),
         magma_gs_pval = pnorm(magma_gs_test_stat / magma_gs_se, lower.tail = FALSE)) %>%
  ungroup() 

# Compute the gene-set type 1 error rates:
null_gene_set_t1_results <- null_gene_set_sim_results %>%
  dplyr::select(-fisher_test_stat, -fisher_df, -stouffer_test_stat,
                -magma_gs_test_stat, -magma_gs_se) %>%
  gather(gene_set_test, gene_set_pval, -gene_set, -test_type) %>%
  mutate(t1_error = ifelse(gene_set_pval <= .05, 1, 0)) %>%
  group_by(gene_set_test, test_type) %>%
  summarize(t1_error_rate = mean(t1_error, na.rm = TRUE),
            n_sims = n()) %>%
  mutate(t1_se = sqrt((t1_error_rate * (1- t1_error_rate)) /
                        n_sims),
         t1_lower = ifelse(t1_error_rate - 2 * t1_se < 0,
                           0,
                           t1_error_rate - 2 * t1_se),
         t1_upper = ifelse(t1_error_rate + 2 * t1_se > 1,
                           1,
                           t1_error_rate + 2 * t1_se))

# Create a figure displaying the type 1 error results for Fisher and Stouffer's method
fisher_gs_t1_error_sim_plot <- null_gene_set_t1_results %>%
  filter(gene_set_test == "fisher_gs_pval") %>%
  mutate(test_type = fct_relevel(test_type,
                                 "browns_one_sided", 
                                 "magma_no_fudge",
                                 "magma_fudge", 
                                 "browns_two_sided", 
                                 "fisher",
                                 "vegas"),
         test_type = fct_recode(test_type,
                                `MAGMA: paper` = "browns_one_sided",
                                `MAGMA: rho^2` = "magma_no_fudge",
                                `MAGMA: code` = "magma_fudge",
                                `Two-sided approx.` = "browns_two_sided",
                                `MC Fisher` = "fisher",
                                `VEGAS` = "vegas"),
         test_type = fct_rev(test_type)) %>%
  ggplot(aes(x = test_type)) + 
  geom_point(aes(y = t1_error_rate),
             size = 2) +
  geom_errorbar(aes(ymin = t1_lower, ymax = t1_upper)) +
  scale_x_discrete(labels = c(`MAGMA: paper` = "MAGMA: paper",
                              `MAGMA: rho^2` = parse(text = TeX('MAGMA: $\\rho^2$')),
                              `MAGMA: code` = "MAGMA: code",
                              `Two-sided approx.` = "Two-sided approx.",
                              `MC Fisher` = "MC Fisher",
                              `VEGAS` = "VEGAS")) +
  scale_y_continuous(limits = c(0, 0.625)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "darkred") +
  theme_bw() + 
  labs(x = "Method", y = "Gene-set type 1 error rate (Fisher's)") +
  coord_flip() + 
  theme(axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12))
stouffer_gs_t1_error_sim_plot <- null_gene_set_t1_results %>%
  filter(gene_set_test == "stouffer_gs_pval") %>%
  mutate(test_type = fct_relevel(test_type,
                                 "browns_one_sided", 
                                 "magma_no_fudge",
                                 "magma_fudge", 
                                 "browns_two_sided", 
                                 "fisher",
                                 "vegas"),
         test_type = fct_recode(test_type,
                                `MAGMA: paper` = "browns_one_sided",
                                `MAGMA: rho^2` = "magma_no_fudge",
                                `MAGMA: code` = "magma_fudge",
                                `Two-sided approx.` = "browns_two_sided",
                                `MC Fisher` = "fisher",
                                `VEGAS` = "vegas"),
         test_type = fct_rev(test_type)) %>%
  ggplot(aes(x = test_type)) + 
  geom_point(aes(y = t1_error_rate),
             size = 2) +
  geom_errorbar(aes(ymin = t1_lower, ymax = t1_upper)) +
  scale_x_discrete(labels = c(`MAGMA: paper` = "MAGMA: paper",
                              `MAGMA: rho^2` = parse(text = TeX('MAGMA: $\\rho^2$')),
                              `MAGMA: code` = "MAGMA: code",
                              `Two-sided approx.` = "Two-sided approx.",
                              `MC Fisher` = "MC Fisher",
                              `VEGAS` = "VEGAS")) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "darkred") +
  scale_y_continuous(limits = c(0, .1)) +
  theme_bw() + 
  labs(x = "Method", y = "Gene-set type 1 error rate (Stouffer's)") +
  coord_flip() + 
  theme(axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12))


# Save an arrangement of these figures together:
other_gs_plots <- plot_grid(fisher_gs_t1_error_sim_plot,
                            stouffer_gs_t1_error_sim_plot,
                            ncol = 1, 
                            labels = c("a", "b"),
                            label_fontface = "bold",
                            align = "hv")
save_plot("figures/si/sf_other_gene_set_t1_error_sims.pdf",
          other_gs_plots, ncol = 1, nrow = 2)
save_plot("nonpdf_figures/si/sf_other_gene_set_t1_error_sims.jpeg",
          other_gs_plots, ncol = 1, nrow = 2)

# ------------------------------------------------------------------------------

# Next create a figure with the p-value distributions for each of the approaches
gene_set_pval_hist <- null_gene_set_sim_results %>%
  dplyr::select(-fisher_test_stat, -fisher_df, -stouffer_test_stat,
                -magma_gs_test_stat, -magma_gs_se) %>%
  gather(gene_set_test, gene_set_pval, -gene_set, -test_type) %>%
  mutate(gene_set_test = fct_relevel(gene_set_test,
                                     "magma_gs_pval", "fisher_gs_pval",
                                     "stouffer_gs_pval"),
         gene_set_test = fct_recode(gene_set_test,
                                    "MAGMA~gene-set" = "magma_gs_pval",
                                    "Fisher~gene-set" = "fisher_gs_pval",
                                    "Stouffer~gene-set" = "stouffer_gs_pval"),
         test_type = fct_relevel(test_type,
                                 "browns_one_sided", 
                                 "magma_no_fudge",
                                 "magma_fudge", 
                                 "browns_two_sided", 
                                 "fisher",
                                 "vegas"),
         test_type = factor(test_type,
                            labels = c("MAGMA:~paper",
                                       "MAGMA:~rho^2",
                                       "MAGMA:~code",
                                       "Two-sided~approx.",
                                       "MC~Fisher",
                                       "VEGAS"))) %>%
  ggplot(aes(x = gene_set_pval)) +
  geom_histogram(breaks = seq(0, 1, by = 0.05),
                 aes(fill = test_type), size = .1,
                 color = "white", alpha = 1,
                 position = "identity") +
  # add expected count lines based on assumption of uniform distribution with 
  # equal probability in each bin:
  geom_hline(yintercept = qbinom(0.025, max(null_gene_set_sim_results$gene_set), 1 / 20),
             color = "darkred", linetype = "dashed", alpha = 1,
             size = 1) +
  geom_hline(yintercept = qbinom(0.975, max(null_gene_set_sim_results$gene_set), 1 / 20),
             color = "darkred", linetype = "dashed", alpha = 1,
             size = 1) +
  ggthemes::scale_fill_colorblind() +
  facet_grid(test_type ~ gene_set_test, 
             scales = "free_y", labeller = label_parsed) +
  theme_bw() +
  labs(x = "Gene-set analysis p-value", y = "Count") + 
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text.y = element_text(size = 24),
        strip.text.x = element_text(size = 24),
        axis.title = element_text(size = 34),
        axis.text = element_text(size = 18))

save_plot("figures/si/sf_gene_set_pval_hist.pdf",
          gene_set_pval_hist, ncol = 3, nrow = 6)
save_plot("nonpdf_figures/si/sf_gene_set_pval_hist.jpeg",
          gene_set_pval_hist, ncol = 3, nrow = 6)

# ------------------------------------------------------------------------------

# Next create a grid of the MAGMA test stat histograms for each gene-level p-value

magma_gs_stat_hist_plots <- null_gene_set_sim_results %>%
  mutate(test_type = fct_relevel(test_type,
                                 "browns_one_sided", 
                                 "magma_no_fudge",
                                 "magma_fudge", 
                                 "browns_two_sided", 
                                 "fisher",
                                 "vegas"),
         test_type = factor(test_type,
                            labels = c("MAGMA:~paper",
                                       "MAGMA:~rho^2",
                                       "MAGMA:~code",
                                       "Two-sided~approx.",
                                       "MC~Fisher",
                                       "VEGAS"))) %>%
  mutate(magma_gs_std_stat = magma_gs_test_stat / magma_gs_se) %>%
  ggplot(aes(x = magma_gs_std_stat)) +
  geom_histogram(aes(y = ..density..,
                     fill = test_type), 
                 breaks = seq(-5, 5, by = 0.25),
                 size = .1,
                 color = "white", alpha = 1,
                 position = "identity") +
  stat_function(fun = dnorm, n = 100000, 
                args = list(mean = 0, 
                            sd = 1),
                color = "darkred", alpha = 0.8, size = .75) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "darkred", alpha = 0.8, size = .75) +
  labs(y = "Density", x = TeX('MAGMA gene-set $z_s$')) +
  facet_wrap(~test_type, labeller = label_parsed,
             ncol = 3) +
  ggthemes::scale_fill_colorblind(guide = FALSE) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 24),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))

# Save these histograms
save_plot("figures/si/sf_magma_gs_test_stat.pdf",
          magma_gs_stat_hist_plots, ncol = 3, nrow = 2)
save_plot("nonpdf_figures/si/sf_magma_gs_test_stat.jpeg",
          magma_gs_stat_hist_plots, ncol = 3, nrow = 2)





