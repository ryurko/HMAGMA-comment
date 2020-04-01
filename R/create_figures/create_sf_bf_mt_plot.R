# PURPOSE: Create supplementary figures for multiple testing results

library(tidyverse)
library(latex2exp)
library(cowplot)

# ------------------------------------------------------------------------------

null_sim_files <- list.files("data/simulations/null/",
                             full.names = TRUE) %>%
  str_remove("_sims_(1|2|3)\\.csv") %>%
  unique()

set.seed(1389)
null_mt_pvalues <- map_dfr(null_sim_files,
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
                                 dplyr::bind_rows(read_csv(paste0(sim_file, "_sims_1.csv")),
                                                  read_csv(paste0(sim_file, "_sims_2.csv")),
                                                  read_csv(paste0(sim_file, "_sims_3.csv")))))
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


# Create figure comparing the distributions of the number of false positives
# with each method:
all_fp_distr <- null_mt_pvalues %>%
  mutate(test_type = fct_relevel(test_type,
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
  ggplot() +
  geom_quasirandom(aes(x = test_type, y = obs_bf_v, color = test_type),
                   alpha = .2) +
  geom_boxplot(aes(x = test_type, y = obs_bf_v), fill = "white",
               color = "black", alpha = 0.8, width = .2,
               outlier.alpha = 0) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkred") +
  theme_bw() + 
  labs(x = "Method", 
       y = TeX('Number of false positives')) +
  scale_x_discrete(labels = c(`MAGMA: paper` = "MAGMA: paper",
                              `MAGMA: rho^2` = parse(text = TeX('MAGMA: $\\rho^2$')),
                              `MAGMA: code` = "MAGMA: code",
                              `Two-sided approx.` = "Two-sided approx.",
                              `MC Fisher` = "MC Fisher",
                              `VEGAS` = "VEGAS")) +
  scale_color_manual(guide = FALSE, 
                     values = rev(ggthemes::colorblind_pal()(6))[1:6]) +
  coord_flip() + 
  theme(axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12),
        panel.grid.major.y = element_blank())


# Next one that displays the proportion of simulations with a number of 
# false positives, but just for the non 'MAGMA: paper' approaches:
fp_facet_prop_points <- null_mt_pvalues %>%
  filter(test_type != "browns_one_sided") %>%
  group_by(test_type, obs_bf_v) %>%
  count() %>%
  ungroup() %>%
  bind_rows(data.frame(test_type = c(rep("magma_fudge", 2), 
                                     rep("browns_two_sided", 2),
                                     rep("fisher", 4),
                                     rep("vegas", 3)),
                       obs_bf_v = c(4, 5, 4, 5, 2, 3, 4, 5, 3, 4, 5),
                       n = rep(0, 11))) %>%
  filter(obs_bf_v < 6) %>%
  mutate(test_type = fct_relevel(test_type,
                                 "magma_no_fudge",
                                 "magma_fudge", 
                                 "browns_two_sided", 
                                 "fisher",
                                 "vegas"),
         test_type = fct_recode(test_type,
                                `MAGMA:~rho^2` = "magma_no_fudge",
                                `MAGMA:~code` = "magma_fudge",
                                `Two-sided~approx.` = "browns_two_sided",
                                `MC~Fisher` = "fisher",
                                `VEGAS` = "vegas")) %>%
  mutate(prop_sims = n / 1000,
         se_prop = sqrt((prop_sims * (1 - prop_sims)) / 1000),
         prop_lower = ifelse(prop_sims - 2 * se_prop < 0,
                             0,
                             prop_sims - 2 * se_prop),
         prop_upper = ifelse(prop_sims + 2 * se_prop > 1,
                             1,
                             prop_sims + 2 * se_prop),
         obs_bf_v = obs_bf_v) %>%
  ggplot(aes(x = obs_bf_v)) +
  geom_line(aes(y = prop_sims, group = test_type),
            color = "gray80", linetype = "dashed", size = .75) +
  geom_point(aes(y = prop_sims, color = test_type),
             position = position_dodge(width = .5),
             size = 2) +
  geom_errorbar(aes(ymin = prop_lower, ymax = prop_upper,
                    color = test_type),
                position = position_dodge(width = .5),
                width = .3) +
  labs(color = "Method", y = "Proportion of simulations",
       x = TeX('Number of false positives')) +
  facet_wrap(~ test_type, ncol = 3, labeller = label_parsed) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = 0:5) +
  scale_color_manual(values = rev(rev(ggthemes::colorblind_pal()(6))[1:5]),
                     labels = rev(c(`MAGMA: rho^2` = parse(text = TeX('MAGMA: $\\rho^2$')),
                                    `MAGMA: code` = "MAGMA: code",
                                    `Two-sided approx.` = "Two-sided approx.",
                                    `MC Fisher` = "MC Fisher",
                                    `VEGAS` = "VEGAS")),
                     guide = FALSE) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "bottom",
        panel.border = element_blank())

# Make the distribution chart the missing panel:
fp_facet_prop_points_combined_chart <-
  ggdraw() +
  draw_plot(plot_grid(fp_facet_prop_points + 
                        theme(axis.title = element_text(size = 20),
                              axis.text.x = element_text(size = 16),
                              axis.text.y = element_text(size = 16),
                              panel.grid.minor.x = element_blank(),
                              legend.position = "none",
                              strip.background = element_blank(),
                              strip.text = element_text(size = 22)),
                      labels = "a", label_fontface = "bold",
                      label_size = 24, vjust = 1.2)) +
  draw_plot(plot_grid(all_fp_distr + 
                        theme(axis.title.x = element_text(size = 10),
                              axis.text.x = element_text(size = 8),
                              axis.title.y = element_blank(),
                              axis.text.y = element_text(size = 8),
                              panel.grid.major.y = element_blank()), 
                      labels = "b", label_size = 24,
                      label_fontface = "bold", vjust = 0.5, hjust = -0.25), 
            x = 0.675, y = 0.05, width = .325, height = .4)

# Save the figure:
save_plot("nonpdf_figures/si/sf_fp_facet_prop_points_distr.jpeg",
          fp_facet_prop_points_combined_chart, ncol = 3, nrow = 2)
