# PURPOSE: Create null test stat distribution figure for the example gene 
#          highlighted in the qq-plot and histograms

library(tidyverse)
library(latex2exp)
library(cowplot)

# ------------------------------------------------------------------------------

# Get the example gene file name:
ex_gene_null_file <- null_sim_files %>%
  str_subset("nsnps772")

# Load the data for this example gene:
example_null_gene_stats <- suppressMessages(
  suppressWarnings(
    dplyr::bind_rows(read_csv(paste0(ex_gene_null_file, "_sims_1.csv")),
                     read_csv(paste0(ex_gene_null_file, "_sims_2.csv")),
                     read_csv(paste0(ex_gene_null_file, "_sims_3.csv")))))


# Get the parameters for MAGMA and the correct two-sided approximation:
ex_gene_twosided_params <- 
  data.frame("scale" = as.numeric(example_null_gene_stats[1, 10]),
             "dof" = as.numeric(example_null_gene_stats[1, 11]))
ex_gene_magma_params <- 
  data.frame("scale" = as.numeric(example_null_gene_stats[1, 2]),
             "dof" = as.numeric(example_null_gene_stats[1, 3]))
ex_gene_brown_params <-
  data.frame("scale" = as.numeric(example_null_gene_stats[1, 14]),
             "dof" = as.numeric(example_null_gene_stats[1, 15]))

null_test_stat_distr_plot <- data.frame("fisher_sim" = example_null_gene_stats[[1]]) %>%
  ggplot(aes(x = fisher_sim)) +
  geom_histogram(color = "gray", fill = "gray", 
                 alpha = 0.5, bins = 50,
                 aes(y = ..density..)) +
  stat_function(fun = dgamma, n = 1000000, 
                args = list(shape = ex_gene_magma_params$dof / 2, 
                            scale = 2 * ex_gene_magma_params$scale),
                color = "#E69F00") +
  stat_function(fun = dgamma, n = 1000000, 
                args = list(shape = ex_gene_twosided_params$dof / 2, 
                            scale = 2 * ex_gene_twosided_params$scale),
                color = "#009E73") +
  stat_function(fun = dgamma, n = 1000000, 
                args = list(shape = ex_gene_brown_params$dof / 2, 
                            scale = 2 * ex_gene_brown_params$scale),
                color = "black") +
  annotate("text", x = 2100, y = 0.001,
           label = "Two-sided approx.",
           color = "#009E73") +
  annotate("text", x = 1900, y = 0.00275,
           label = "MAGMA: paper",
           color = "black") +
  annotate("text", x = 2000, y = 0.0015,
           label = TeX("MAGMA: $\\rho^2$"),
           color = "#E69F00") +
  theme_bw() +
  labs(y = "Density",
       x = TeX('T = -2 $\\sum_i^m log p_j$'))

# save_plot("figures/si/sf1_null_test_stat_example.pdf",
#           null_test_stat_distr_plot)
# save_plot("nonpdf_figures/si/sf1_null_test_stat_example.jpeg",
#           null_test_stat_distr_plot)
