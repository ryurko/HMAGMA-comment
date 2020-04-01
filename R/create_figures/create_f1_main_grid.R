# PURPOSE: Generate the grid of results for Figure 1

library(tidyverse)
library(latex2exp)
library(cowplot)

# ------------------------------------------------------------------------------

# FIGURE A

# First create the figure displaying a comparison in the covariance approximations

# Brown's approximation is the following:
# cov_{jk} = \begin{cases} \rho_{jk} \cdot (3.25 + 0.75 \rho_{jk}) \text{ if } \rho_{jk} \geq 0\\
#                          \rho_{jk} \cdot (3.27 + 0.71 \rho_{jk}) \end{cases}
# Meanwhile the two-sided approximation is:
# cov_{jk} =  3.9081 \rho_{jk}^2 + 0.0313 \rho_{jk}^4 + 0.1022 \rho_{jk}^6 - 0.1378 \rho_{jk}^8 + 0.0941 \rho_{jk}^10 

# Create a display of these different two approximations as a function of rho_{ij}
cov_approx_line_plot <- tibble(x = c(-1, 1)) %>%
  ggplot(aes(x)) +
  stat_function(fun = function(.x) 3.9081*(.x)^2 + 0.0313*(.x)^4 + 0.1022*(.x)^6 -
                  0.1378*(.x)^8 + 0.0941*(.x)^(10),
                color = "#009E73", size = 1.5) +
  stat_function(fun = function(.x) ifelse(.x >= 0, 
                                          .x * (3.25 + 0.75 * .x),
                                          .x * (3.27 + 0.71 * .x)),
                linetype = "dashed", color = "#000000",
                size = 1.5) +
  stat_function(fun = function(.x) (.x)^2 * (3.25 + 0.75 * (.x)^2),
                color = "#E69F00", linetype = "dashed",
                size = 1.5) +
  annotate("text", x = 0, y = 2.5, size = 10,
           label = "Two-sided\n approximation",
           color = "#009E73") +
  annotate("text", x = 0, y = -1.5, size = 10,
           label = "MAGMA: paper\n(One-sided approximation)",
           color = "#000000") +
  annotate("text", x = -0.75, y = .5, size = 10,
           label = TeX('MAGMA: $\\rho^2$'), 
           color = "#E69F00") +
  geom_hline(yintercept = 0, color = "gray", linetype = "dotted") +
  labs(x = TeX('$\\rho$'),
       y = TeX('$\\approx$ covariance')) +
  theme_cowplot() +
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 34))

# ------------------------------------------------------------------------------

# FIGURE B

# Next a comparison of the Type 1 error rate results

# Make a vector of the null sim files - with the test suffix removed - this
# will make it easier to then load both datasets
null_sim_files <- list.files("data/simulations/null/",
                             full.names = TRUE) %>%
  str_remove("_sims_(1|2|3)\\.csv") %>%
  unique()

# Create a data frame with the type 1 error rate simulations to visualize:
null_sim_result_summary <- map_dfr(null_sim_files,
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
                                     
                                     # Compute the type 1 error rate for each gene-level test:
                                     gene_results <- gene_pvals %>%
                                       mutate(sim_index = 1:n()) %>%
                                       gather(test_type, pval, -sim_index) %>%
                                       mutate(t1_error = ifelse(pval <= .05, 1, 0)) %>%
                                       group_by(test_type) %>%
                                       summarize(t1_error_rate = mean(t1_error, na.rm = TRUE),
                                                 n_sims = n()) %>%
                                       mutate(t1_error_se = sqrt((t1_error_rate * (1- t1_error_rate)) /
                                                                   n_sims),
                                              t1_error_lower = ifelse(t1_error_rate - 2 * t1_error_se < 0,
                                                                      0,
                                                                      t1_error_rate - 2 * t1_error_se),
                                              t1_error_upper = ifelse(t1_error_rate + 2 * t1_error_se > 1,
                                                                      1,
                                                                      t1_error_rate + 2 * t1_error_se)) %>%
                                       mutate(gene = gene_name,
                                              chr = gene_chr,
                                              n_snps = gene_n_snps)
                                     return(gene_results)
                                   })

# Visualize the type 1 error rate sims with plus/minus two standard errors
t1_error_sim_plot <- null_sim_result_summary %>%
  mutate(test_type = fct_relevel(test_type,
                                 "browns_one_sided", "magma_no_fudge",
                                 "magma_fudge", "browns_two_sided",
                                 "fisher", 
                                 "vegas"),
         test_type = fct_recode(test_type,
                                `MAGMA: paper` = "browns_one_sided",
                                `MAGMA: rho^2` = "magma_no_fudge",
                                `MAGMA: code` = "magma_fudge",
                                `Two-sided approx.` = "browns_two_sided",
                                `MC Fisher` = "fisher",
                                `VEGAS` = "vegas"),
         test_type = fct_rev(test_type),
         gene_label = paste0(n_snps, " SNPs"),
         gene_label = fct_reorder(gene_label, n_snps)) %>%
  ggplot(aes(x = test_type)) + 
  scale_x_discrete(labels = c(`MAGMA: paper` = "MAGMA: paper",
                              `MAGMA: rho^2` = parse(text = TeX('MAGMA: $\\rho^2$')),
                              `MAGMA: code` = "MAGMA: code",
                              `Two-sided approx.` = "Two-sided approx.",
                              `MC Fisher` = "MC Fisher",
                              `VEGAS` = "VEGAS")) +
  geom_point(aes(y = t1_error_rate), size = 2) +
  geom_errorbar(aes(ymin = t1_error_lower, ymax = t1_error_upper)) +
  theme_bw() + 
  labs(x = "Method", y = "Type 1 error rate") +
  scale_y_continuous(limits = c(0, .175)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "darkred") +
  facet_wrap(~gene_label, ncol = 3) +
  coord_flip() +
  theme(strip.background = element_blank(),
        plot.title = element_text(size = 32),
        plot.subtitle = element_text(size = 24),
        strip.text = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 18))

# ------------------------------------------------------------------------------

# FIGURES C AND D

# C: qq-plot for example gene with 772 SNPs
ex_gene_null_file <- null_sim_files %>%
  str_subset("nsnps772")

# Load the data for this example gene:
inflated_gene_pvals <- suppressMessages(
  suppressWarnings(
    dplyr::bind_rows(read_csv(paste0(ex_gene_null_file, "_sims_1.csv")),
                     read_csv(paste0(ex_gene_null_file, "_sims_2.csv")),
                     read_csv(paste0(ex_gene_null_file, "_sims_3.csv")))))

# Name the columns accordingly:
inflated_gene_pvals <- inflated_gene_pvals[, c(4, 8, 12, 16, 18, 20)]
colnames(inflated_gene_pvals) <- 
  c("magma_fudge", "magma_no_fudge",
    "browns_two_sided", "browns_one_sided",
    "fisher", "vegas")

# Generate the qq-plot comparing each gene-level p-value calculation:
ex_gene_qqplots_colored <- inflated_gene_pvals %>%
  gather(test_type, pval) %>%
  mutate(test_type = fct_relevel(test_type,
                                 "browns_one_sided", 
                                 "magma_no_fudge",
                                 "magma_fudge", "browns_two_sided", 
                                 "fisher",
                                 "vegas"),
         test_type = fct_recode(test_type,
                                `MAGMA: paper` = "browns_one_sided",
                                `MAGMA: rho^2` = "magma_no_fudge",
                                `MAGMA: code` = "magma_fudge",
                                `Two-sided approx.` = "browns_two_sided",
                                `MC Fisher` = "fisher",
                                `VEGAS` = "vegas")) %>%
  group_by(test_type) %>%
  arrange(pval) %>%
  mutate(pval_index = 1:n(),
         neglog_p = -log10(pval),
         neglog_p = ifelse(is.infinite(neglog_p) & neglog_p > 0,
                           -log10(.Machine$double.eps), neglog_p),
         neglog_p = ifelse(is.infinite(neglog_p) & neglog_p < 0,
                           -log10(1 - .Machine$double.eps), neglog_p),
         exp_neglog_p = -log10(ppoints(n())),
         upper_beta = -log10(qbeta(.05 / 2, 
                                   pval_index, n() - pval_index)),
         lower_beta = -log10(qbeta(1 - (.05 / 2), 
                                   pval_index, n() - pval_index))) %>%
  ungroup() %>%
  ggplot(aes(x = exp_neglog_p)) +
  geom_line(aes(y = upper_beta), alpha = 0.5,
            linetype = "dashed", color = "darkred") +
  geom_line(aes(y = lower_beta), alpha = 0.5,
            linetype = "dashed", color = "darkred") +
  geom_point(aes(y = neglog_p, color = test_type),
             alpha = 0.75, size = 0.75) +
  ggthemes::scale_color_colorblind(labels = c(`MAGMA: paper` = "MAGMA: paper",
                                              `MAGMA: rho^2` = parse(text = TeX('MAGMA: $\\rho^2$')),
                                              `MAGMA: code` = "MAGMA: code",
                                              `Two-sided approx.` = "Two-sided approx.",
                                              `MC Fisher` = "MC Fisher",
                                              `VEGAS` = "VEGAS")) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.3,
              linetype = "dashed", color = "gray") +
  scale_y_continuous(limits = c(0, 7)) +
  theme_cowplot() +
  guides(color = guide_legend(override.aes = list(size = 10, alpha = 1)))+
  labs(x = TeX('Expected -log$_{10}$(gene-level p-value)'),
       y = TeX('Observed -log$_{10}$(gene-level p-value)'),
       color = "Method") + 
  theme(legend.position = c(0.125, 0.725),
        legend.direction = "vertical",
        plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        legend.background = element_rect(fill = alpha("white", alpha = 0.5)))


# D: histograms for example gene

# Create a list of the different histograms for each method - this is just to
# make it easier for setting the y-axes compared to using facets
method_types <- c("browns_one_sided", "magma_no_fudge", "magma_fudge", 
                  "browns_two_sided", "fisher", "vegas")
method_labels <- c("MAGMA: paper", "MAGMA: rho^2", "MAGMA: code",
                   "Two-sided approx.", "MC Fisher", "VEGAS")

# Loop through this to make a list of the plots:
null_hist_list <- lapply(1:length(method_types),
                         function(method_i) {
                           hist_plot <- inflated_gene_pvals %>%
                             gather(test_type, pval) %>%
                             filter(test_type == method_types[method_i]) %>%
                             ggplot(aes(x = pval)) +
                             geom_histogram(breaks = seq(0, 1, by = 0.05),
                                            fill = ggthemes::colorblind_pal()(6)[method_i], 
                                            size = .1,
                                            color = "white", alpha = 1,
                                            position = "identity") +
                             geom_hline(yintercept = qbinom(0.025, nrow(inflated_gene_pvals), 1 / 20),
                                        color = "darkred", linetype = "dashed", alpha = .5,
                                        size = .5) +
                             geom_hline(yintercept = qbinom(0.975, nrow(inflated_gene_pvals), 1 / 20),
                                        color = "darkred", linetype = "dashed", alpha = .5,
                                        size = .5) +
                             theme_bw() +
                             theme(plot.margin = unit(c(0, 0.1, 0, 0.1), "cm"),
                                   strip.background = element_blank(),
                                   axis.title = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.text.y = element_blank(),
                                   axis.text.x = element_text(size = 6))
                           
                           # If not Brown's one-sided then use the same y-axis limits:
                           if (method_i > 1) {
                             hist_plot <- hist_plot +
                               scale_y_continuous(limits = c(0, 62500))
                           }
                           
                           # Next add the titles depending on if it is the rho^2:
                           if (method_types[method_i] == "magma_no_fudge") {
                             hist_plot <- hist_plot +
                               labs(title = TeX('MAGMA: $\\rho^2$'))
                           } else {
                             hist_plot <- hist_plot +
                               labs(title = method_labels[method_i])
                           }
                           
                           hist_plot +
                             theme(plot.title = element_text(hjust = 0.5, vjust = 0,
                                                             size = 10))
                           
                         })

# Generate the grid of these histograms:
null_hist_plots <- plot_grid(plot_grid(plotlist = null_hist_list, ncol = 3, align = "hv"),
                             ggdraw() + 
                               draw_label(
                                 "Gene-level p-value",
                                 fontface = 'plain',
                                 x = 0,
                                 hjust = 0.5,
                                 size = 10
                               ) + theme(
                                 # add margin so label is centered
                                 plot.margin = margin(0, 0, 0, 11.75, unit = "cm")
                               ), ncol = 1, rel_heights = c(1, .08))


# Place this image in the bottom right corner using cowplot
ex_gene_null_qqplot_hist <-
  ggdraw() +
  draw_plot(plot_grid(ex_gene_qqplots_colored +
                        theme(plot.title = element_text(size = 24),
                              plot.subtitle = element_text(size = 18),
                              axis.title = element_text(size = 20),
                              axis.text = element_text(size = 14),
                              legend.text = element_text(size = 24),
                              legend.title = element_blank()),
                      labels = "c", label_fontface = "bold",
                      label_size = 32, vjust = .5)) +
  draw_plot(plot_grid(null_hist_plots, 
                      labels = "d", label_size = 32,
                      label_fontface = "bold", vjust = 0.5, hjust = -0.25), 
            x = 0.52, y = .09, width = .46, height = .3)

# ------------------------------------------------------------------------------

# FIGURES E AND F

# E: FWER error rate control:

# Set-up the dataset for plotting with multiple testing sets of 9000 genes,
# where there are 1000 from each gene:
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
         bf_error = as.numeric(bf_pval <= 0.05)) %>%
  summarize(n_genes = n(),
            obs_fwer = ifelse(any(bf_error == 1), 1, 0),
            obs_bf_v = sum(bf_error)) %>%
  ungroup()

# Create the plot displaying comparison of FWER control:
mt_fwer_sim_plot <- null_mt_pvalues %>%
  group_by(test_type) %>%
  summarize(obs_fwer_bf = mean(obs_fwer),
            n_sims = n()) %>%
  mutate(obs_fwer_se = sqrt((obs_fwer_bf * (1- obs_fwer_bf)) /
                              n_sims)) %>%
  mutate(obs_fwer_lower = pmax(0, obs_fwer_bf - 2 * obs_fwer_se),
         obs_fwer_upper = obs_fwer_bf + 2 * obs_fwer_se,
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
  ggplot(aes(x = test_type, y = obs_fwer_bf)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = obs_fwer_lower, ymax = obs_fwer_upper)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "darkred") +
  theme_bw() + 
  labs(x = "Method", y = "Observed FWER") +
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

# F: Gene set type 1 error rate
set.seed(1389)
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


# Just zoom in on the MAGMA gene set test results:
magma_gs_t1_error_sim_plot <- null_gene_set_t1_results %>%
  filter(gene_set_test == "magma_gs_pval") %>%
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
  scale_y_continuous(limits = c(0, 0.06)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "darkred") +
  theme_bw() + 
  labs(x = "Method", y = "Gene-set type 1 error rate") +
  coord_flip() + 
  theme(axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12))


# ------------------------------------------------------------------------------

# Now arrange in a grid and save:

magma_plot_grid <- plot_grid(plot_grid(cov_approx_line_plot, t1_error_sim_plot,
                                       ncol = 2, rel_widths = c(1, 2),
                                       labels = c("a", "b"), label_fontface = "bold",
                                       label_size = 32),
                             ggdraw() + 
                               draw_label(
                                 "",
                                 x = 0,
                                 hjust = 0
                               ),
                             plot_grid(ex_gene_null_qqplot_hist,
                                       plot_grid(mt_fwer_sim_plot + theme(axis.title.x = element_text(size = 24),
                                                                          axis.text.x = element_text(size = 18),
                                                                          axis.title.y = element_blank(),
                                                                          axis.text.y = element_text(size = 18)), 
                                                 magma_gs_t1_error_sim_plot + theme(axis.title.x = element_text(size = 24),
                                                                                    axis.text.x = element_text(size = 18),
                                                                                    axis.title.y = element_blank(),
                                                                                    axis.text.y = element_text(size = 18)),
                                                 ncol = 1, align = "hv", labels = c("e", "f"), vjust = .5, hjust = -0.3,
                                                 label_fontface = "bold", label_size = 32),
                                       ncol = 2, rel_widths = c(2, 1), rel_heights = c(1, .8)),
                             ncol = 1, rel_heights = c(1, .05, 1))

# save_plot("figures/main/f1_grid.pdf",
#           magma_plot_grid, ncol = 5, nrow = 5)
# save_plot("nonpdf_figures/main/f1_grid.jpeg",
#           magma_plot_grid, ncol = 5, nrow = 5)







