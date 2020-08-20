# PURPOSE: Create figure 1 displaying the simulation based results 

library(tidyverse)
library(latex2exp)
library(cowplot)

# -------------------------------------------------------------------------

# Fig 1a - multivariate normal example ------------------------------------

# Create a function to generate blocks of p-values that are underlying correlated
# Gaussians, compute the test statistics - then compute the chi-squared approximation
# based on the actual correlation matrix
create_single_rho_block <- function(rho, block_dim) {
  block <- matrix(rho, nrow = block_dim, ncol = block_dim)
  diag(block) <- 1
  return(block)
}

# Generate the distribution of test statistics for a single correlation value,
# will use correlation of 0.5 with block size of 50:
set.seed(1865)
sim_pval_test_stats <- map_dfr(1:1000000,
                               function(i) {
                                 test_cor_block <- create_single_rho_block(.5, 50)
                                 MASS::mvrnorm(1, rep(0, 50), 
                                               Sigma = test_cor_block) %>%
                                   sapply(., function(z) -2 * log(2 * pnorm(-abs(z)))) %>%
                                   sum() %>%
                                   tibble("fisher_stat" = .)
                               })

# Save this dataset to access later:
write_csv(sim_pval_test_stats, "data/mvnorm_block_sims.csv")
#sim_pval_test_stats <- read_csv("data/mvnorm_block_sims.csv")

# Function to return just the upper diagonal values of the correlation matrix as a vector:
create_single_rho_upper_diag_vector <- function(rho, block_dim) {
  block <- matrix(rho, nrow = block_dim, ncol = block_dim)
  diag(block) <- 0 # since we dont want the diagonal elements
  return(block[upper.tri(block, diag = FALSE)])
}


# Write a function that returns both the scale and degrees of freedom for a given
# block size and correlation value
compute_browns_approx <- function(rho, block_dim, type) {
  
  # Next what's the expected value of the test stat:
  test_expected_value <- 2 * block_dim
  
  test_cor_vector <- create_single_rho_upper_diag_vector(rho, block_dim)
  
  if (type == "browns") {
    cov_vector <- ifelse(test_cor_vector >= 0,
                         test_cor_vector * (3.25 + 0.75 * test_cor_vector),
                         test_cor_vector * (3.27 + 0.71 * test_cor_vector))
  } else if (type == "squared") {
    cov_vector <- (test_cor_vector^2) * (3.25 + 0.75 * (test_cor_vector^2))
  } else {
    cov_vector <- 
      3.9081 * (test_cor_vector^2) + 
      0.0313 * (test_cor_vector^4) + 
      0.1022 * (test_cor_vector^6) - 
      0.1378 * (test_cor_vector^8) + 
      0.0941 * (test_cor_vector^10)
  }
  cov_sum <- 2 * sum(cov_vector, na.rm = TRUE)
  
  test_variance <- 4 * block_dim + cov_sum
  
  test_c <- test_variance / (2 * test_expected_value)
  test_dof <- (2 * test_expected_value^2) / test_variance
  
  return(list("c" = test_c,
              "dof" = test_dof))
  
}

# Table of the different approximations with rho = 0.5 and m = 50 SNPs
magma_approx_tbl <- map_dfr(c("browns", "squared", "twosided"),
                            function(cov_type) {
                              as_tibble(compute_browns_approx(0.5,
                                                              50, cov_type)) %>%
                                mutate(type = cov_type)
                            })

# Will be using the Nature color palette:
ggsci::pal_npg()(5)
# [1] "#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF"

null_test_stat_ex_plot <- sim_pval_test_stats %>%
  ggplot(aes(x = fisher_stat)) +
  geom_histogram(color = "white", fill = "#F39B7FFF",
                 bins = 50, closed = "left", size = 0.1,
                 aes(y = ..density..)) +
  #geom_density(fill = "#F39B7FFF", 
  #             alpha = .8) +
  stat_function(fun = dgamma, n = 1000000, 
                args = list(shape = pull(filter(magma_approx_tbl, type == "squared"), dof) / 2, 
                            scale = 2 * pull(filter(magma_approx_tbl, type == "squared"), c)),
                color = "#4DBBD5FF") +
  stat_function(fun = dgamma, n = 1000000, 
                args = list(shape = pull(filter(magma_approx_tbl, type == "twosided"), dof) / 2, 
                            scale = 2 * pull(filter(magma_approx_tbl, type == "twosided"), c)),
                color = "#00A087FF") +
  stat_function(fun = dgamma, n = 1000000, 
                args = list(shape = pull(filter(magma_approx_tbl, type == "browns"), dof) / 2, 
                            scale = 2 * pull(filter(magma_approx_tbl, type == "browns"), c)),
                color = "#E64B35FF") +
  scale_x_continuous(limits = c(0, 450)) +
  annotate("text", x = 200, y = 0.00625,
           size = 6,
           label = "Two-sided approx.",
           color = "#00A087FF") +
  annotate("text", x = 275, y = 0.0025,
           size = 6,
           label = "MAGMA: paper",
           color = "#E64B35FF") +
  annotate("text", x = 150, y = 0.01,
           size = 6,
           label = TeX("MAGMA: $\\rho^2$"),
           color = "#4DBBD5FF") +
  annotate("text", x = 175, y = 0.01375,
           size = 6,
           label = TeX("Correct null distribution"),
           color = "#F39B7FFF") +
  theme_cowplot() +
  labs(y = "Density",
       x = TeX('T = -2 $\\sum_j^m log p_j$'))


# Save the individual plot
save_plot("figures/pdf/main/f1a_null_test_stat_example.pdf",
          null_test_stat_ex_plot, base_asp = 2)
save_plot("figures/nonpdf/main/f1a_null_test_stat_example.jpg",
          null_test_stat_ex_plot, base_asp = 2)


# Just density curve based version:
null_test_stat_ex_plot_curve <- sim_pval_test_stats %>%
  ggplot(aes(x = fisher_stat)) +
  geom_density(color = "#F39B7FFF") +
  stat_function(fun = dgamma, n = 1000000, 
                args = list(shape = pull(filter(magma_approx_tbl, type == "squared"), dof) / 2, 
                            scale = 2 * pull(filter(magma_approx_tbl, type == "squared"), c)),
                color = "#4DBBD5FF") +
  stat_function(fun = dgamma, n = 1000000, 
                args = list(shape = pull(filter(magma_approx_tbl, type == "twosided"), dof) / 2, 
                            scale = 2 * pull(filter(magma_approx_tbl, type == "twosided"), c)),
                color = "#00A087FF") +
  stat_function(fun = dgamma, n = 1000000, 
                args = list(shape = pull(filter(magma_approx_tbl, type == "browns"), dof) / 2, 
                            scale = 2 * pull(filter(magma_approx_tbl, type == "browns"), c)),
                color = "#E64B35FF") +
  scale_x_continuous(limits = c(0, 450)) +
  annotate("text", x = 200, y = 0.00625,
           size = 6,
           label = "Two-sided approx.",
           color = "#00A087FF") +
  annotate("text", x = 275, y = 0.0025,
           size = 6,
           label = "MAGMA: paper",
           color = "#E64B35FF") +
  annotate("text", x = 150, y = 0.01,
           size = 6,
           label = TeX("MAGMA: $\\rho^2$"),
           color = "#4DBBD5FF") +
  annotate("text", x = 175, y = 0.01375,
           size = 6,
           label = TeX("Correct null distribution"),
           color = "#F39B7FFF") +
  theme_cowplot() +
  labs(y = "Density",
       x = TeX('T = -2 $\\sum_j^m log p_j$'))
save_plot("figures/pdf/main/f1a_null_test_stat_example_curves.pdf",
          null_test_stat_ex_plot_curve, base_asp = 2)
save_plot("figures/nonpdf/main/f1a_null_test_stat_example_curves.jpg",
          null_test_stat_ex_plot_curve, base_asp = 2)

# Will use the curves

# Fig 1b - type 1 error rate results --------------------------------------

# Load the smallest gene (1 million sims)
# Gene: 3848
set.seed(3848)
gene3848_type1_null_sims <-
  map_dfr(list.files("data/updated_null_sims/gene3848_1/", full.names = TRUE),
          read_csv) %>%
  sample_n(1000000, replace = FALSE) %>%
  # Randomly assign them
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(gene = "Small gene",
         fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))

# Load the largest gene (again 1 million sims)
# Gene: 123624
set.seed(123624)
gene123624_type1_null_sims <-
  map_dfr(list.files("data/updated_null_sims/gene123624/", full.names = TRUE),
          read_csv) %>%
  sample_n(1000000, replace = FALSE) %>%
  # Randomly assign them
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(gene = "Large gene",
         fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))


# Stack these together - selecting only the necessary p-value columns after
# renaming them:
gene_t1_sim_results <- 
  bind_rows(mutate(gene3848_type1_null_sims, 
                   sim_index = 1:n()),
            mutate(gene123624_type1_null_sims, 
                   sim_index = 1:n())) %>%
  dplyr::select(-vegas) %>%
  pivot_longer(cols = magma_fudge:fisher,
               names_to = "method", values_to = "pval") %>%
  mutate(t1_error = ifelse(pval <= .05, 1, 0)) %>%
  group_by(gene, method) %>%
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

t1_error_sim_plot <- gene_t1_sim_results %>%
  ungroup() %>%
  mutate(method = fct_relevel(method,
                              "browns_one_sided", 
                              "magma_no_fudge",
                              "browns_two_sided", 
                              "magma_fudge", 
                              "fisher"),
         method = fct_recode(method,
                             `MAGMA: paper` = "browns_one_sided",
                             `MAGMA: rho^2` = "magma_no_fudge",
                             `MAGMA: code` = "magma_fudge",
                             `Two-sided approx.` = "browns_two_sided",
                             `Corrected` = "fisher"),
         method = fct_rev(method),
         gene = fct_relevel(gene, "Small gene")) %>%
  ggplot(aes(x = method)) + 
  scale_x_discrete(labels = c(`MAGMA: paper` = "MAGMA: paper",
                              `MAGMA: rho^2` = parse(text = TeX('MAGMA: $\\rho^2$')),
                              `MAGMA: code` = "MAGMA: code",
                              `Two-sided approx.` = "Two-sided approx.",
                              `Corrected` = "Corrected")) +
  geom_point(aes(y = t1_error_rate), size = 2) +
  geom_errorbar(aes(ymin = t1_error_lower, ymax = t1_error_upper)) +
  theme_bw() + 
  labs(x = "Method", y = "Type 1 error rate") +
  scale_y_continuous(limits = c(.025, .175)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "darkred") +
  facet_wrap(~gene, ncol = 3) +
  coord_flip() +
  theme(strip.background = element_blank(),
        plot.title = element_text(size = 32),
        plot.subtitle = element_text(size = 24),
        strip.text = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 18))

# Save the individual plot
save_plot("figures/pdf/main/f1b_t1_error.pdf",
          t1_error_sim_plot, ncol = 2, nrow = 1)
save_plot("figures/nonpdf/main/f1b_t1_error.jpg",
          t1_error_sim_plot, ncol = 2, nrow = 1)


# Fig 1c - example testing gene simulations -------------------------------

# Remove the large simulation data files that are no longer needed:
rm(list = c("sim_pval_test_stats", 
            "gene3848_type1_null_sims", "gene123624_type1_null_sims"))

# Gene: 3848
set.seed(3848)
gene3848_null_sims <-
  map_dfr(c(list.files("data/updated_null_sims/gene3848_1/", full.names = TRUE),
            list.files("data/updated_null_sims/gene3848_2/", full.names = TRUE),
            list.files("data/updated_null_sims/gene3848_3/", full.names = TRUE)),
          read_csv) %>%
  sample_n(14663000, replace = FALSE) %>%
  # Randomly assign them
  mutate(sim_group = sample(rep(1:1000, 14663))) %>%
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))

# Gene: 51306
set.seed(51306)
gene51306_null_sims <-
  map_dfr(c(list.files("data/updated_null_sims/gene51306_1/", full.names = TRUE),
            list.files("data/updated_null_sims/gene51306_2/", full.names = TRUE),
            list.files("data/updated_null_sims/gene51306_3/", full.names = TRUE)),
          read_csv) %>% 
  sample_n(14663000, replace = FALSE) %>%
  # Randomly assign them
  mutate(sim_group = sample(rep(1:1000, 14663))) %>%
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))

# Gene: 659
set.seed(659)
gene659_null_sims <-
  map_dfr(c(list.files("data/updated_null_sims/gene659_1/", full.names = TRUE),
            list.files("data/updated_null_sims/gene659_2/", full.names = TRUE),
            list.files("data/updated_null_sims/gene659_3/", full.names = TRUE)),
          read_csv) %>%  
  sample_n(14663000, replace = FALSE) %>%
  # Randomly assign them
  mutate(sim_group = sample(rep(1:1000, 14663))) %>%
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))

# Gene: 169355
set.seed(169355)
gene169355_null_sims <-
  map_dfr(list.files("data/updated_null_sims/gene169355/", full.names = TRUE),
          read_csv) %>%
  sample_n(2827000, replace = FALSE) %>%
  # Randomly assign them
  mutate(sim_group = sample(rep(1:1000, 2827))) %>%
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))

# Gene: 80243
set.seed(80243)
gene80243_null_sims <-
  map_dfr(list.files("data/updated_null_sims/gene80243/", full.names = TRUE),
          read_csv) %>%
  sample_n(2827000, replace = FALSE) %>%
  # Randomly assign them
  mutate(sim_group = sample(rep(1:1000, 2827))) %>%
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))

# Gene: 26064
set.seed(26064)
gene26064_null_sims <-
  map_dfr(list.files("data/updated_null_sims/gene26064/", full.names = TRUE),
          read_csv) %>%
  sample_n(2827000, replace = FALSE) %>%
  # Randomly assign them
  mutate(sim_group = sample(rep(1:1000, 2827))) %>%
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))

# Gene: 6263
set.seed(6263)
gene6263_null_sims <-
  map_dfr(list.files("data/updated_null_sims/gene6263/", full.names = TRUE),
          read_csv) %>%
  sample_n(177000, replace = FALSE) %>%
  # Randomly assign them
  mutate(sim_group = sample(rep(1:1000, 177))) %>%
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE))  %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))

# Gene: 11122
set.seed(11122)
gene11122_null_sims <-
  map_dfr(list.files("data/updated_null_sims/gene11122/", full.names = TRUE),
          read_csv) %>%
  sample_n(177000, replace = FALSE) %>%
  # Randomly assign them
  mutate(sim_group = sample(rep(1:1000, 177))) %>%
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))

# Gene: 123624
set.seed(123624)
gene123624_null_sims <-
  map_dfr(list.files("data/updated_null_sims/gene123624/", full.names = TRUE),
          read_csv) %>%
  sample_n(177000, replace = FALSE) %>%
  # Randomly assign them
  mutate(sim_group = sample(rep(1:1000, 177))) %>%
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))


# Stack these together - selecting only the necessary p-value columns after
# renaming them:
gene_mt_sims <- bind_rows(mget(str_subset(ls(), "_null_sims"))) 

# Remove the additional large gene sets:
rm(list = str_subset(ls(), "_null_sims"))

# Just use the first sim group to display its p-value distributions:
# Randomly pick a set:
set.seed(1389)
pick_sim_group <- sample(1:max(gene_mt_sims$sim_group), 1)

ex_mt_gene_set <- gene_mt_sims %>% 
  filter(sim_group == pick_sim_group) %>%
  dplyr::select(-vegas)

# Create a list of the different histograms for each method - this is just to
# make it easier for setting the y-axes compared to using facets
method_types <- c("browns_one_sided", "magma_no_fudge",  
                  "browns_two_sided", "magma_fudge",
                  "fisher")
method_labels <- c("MAGMA: paper", "MAGMA: rho^2", 
                   "Two-sided approx.", "MAGMA: code", "Corrected")
method_colors <-  ggsci::pal_npg()(5)


# Loop through this to make a list of the plots:
null_hist_list <- lapply(1:length(method_types),
                         function(method_i) {
                           hist_plot <- ex_mt_gene_set %>%
                             pivot_longer(cols = magma_fudge:fisher,
                                          names_to = "method", values_to = "pval")  %>%
                             filter(method == method_types[method_i]) %>%
                             ggplot(aes(x = pval)) +
                             geom_histogram(breaks = seq(0, 1, by = 0.05),
                                            closed = "left",
                                            fill = method_colors[method_i], 
                                            size = .1,
                                            color = "white", alpha = 1,
                                            position = "identity") +
                             geom_hline(yintercept = qbinom(0.025, nrow(ex_mt_gene_set), 1 / 20),
                                        color = "darkred", linetype = "dashed", alpha = .5,
                                        size = .5) +
                             geom_hline(yintercept = qbinom(0.975, nrow(ex_mt_gene_set), 1 / 20),
                                        color = "darkred", linetype = "dashed", alpha = .5,
                                        size = .5) +
                             theme_bw() +
                             theme(plot.margin = unit(c(0, 0.1, 0, 0.1), "cm"),
                                   strip.background = element_blank(),
                                   strip.text = element_text(size = 12),
                                   axis.title = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.text.y = element_blank(),
                                   axis.text.x = element_text(size = 10))
                           
                           # If not Brown's one-sided then use the same y-axis limits:
                           if (method_i > 1) {
                             hist_plot <- hist_plot +
                               scale_y_continuous(limits = c(0, 3500))
                           }
                           
                           # Next add the titles depending on if it is the rho^2:
                           if (method_types[method_i] == "magma_no_fudge") {
                             hist_plot <- hist_plot +
                               labs(title = TeX('MAGMA: $\\rho^2$'))
                           } else {
                             hist_plot <- hist_plot +
                               labs(title = method_labels[method_i])
                           }
                           
                           # Add x-axis label for bottom:
                           if (method_types[method_i] == "fisher") {
                             hist_plot <- hist_plot +
                               labs(x = "Gene-level p-value") +
                               theme(axis.title.x = element_text(size = 18))
                           }
                           
                           hist_plot +
                             theme(plot.title = element_text(hjust = 0.5, vjust = 0,
                                                             size = 18))
                           
                         })

null_hist_plot_all <-
  plot_grid(plotlist = null_hist_list, ncol = 1, align = "hv")

# Save the individual plot
save_plot("figures/pdf/main/f1c_null_hist_all.pdf",
          null_hist_plot_all, ncol = 1, nrow = 5)
save_plot("figures/nonpdf/main/f1c_null_hist_all.jpg",
          null_hist_plot_all, ncol = 1, nrow = 5)


# Create another option which has the corrected version displayed behind each one

# Loop through this to make a list of the plots:
null_hist_list_upd <- lapply(1:(length(method_types) - 1),
                         function(method_i) {
                           hist_data <- ex_mt_gene_set %>%
                             pivot_longer(cols = magma_fudge:fisher,
                                          names_to = "method", values_to = "pval")
                           hist_plot <- hist_data %>%
                             filter(method == "fisher") %>%
                             ggplot(aes(x = pval)) +
                             geom_histogram(breaks = seq(0, 1, by = 0.05),
                                            closed = "left",
                                            fill = method_colors[5], 
                                            size = .1,
                                            color = "white", alpha = 1,
                                            position = "identity") +
                             geom_histogram(data = 
                                              filter(hist_data,
                                                     method %in% method_types[method_i]),
                                            breaks = seq(0, 1, by = 0.05),
                                            closed = "left",
                                            alpha = 0.75,
                                            fill = method_colors[method_i], 
                                            size = .1,
                                            color = "white", 
                                            position = "identity") +
                             geom_hline(yintercept = qbinom(0.025, nrow(ex_mt_gene_set), 1 / 20),
                                        color = "darkred", linetype = "dashed", alpha = .5,
                                        size = .5) +
                             geom_hline(yintercept = qbinom(0.975, nrow(ex_mt_gene_set), 1 / 20),
                                        color = "darkred", linetype = "dashed", alpha = .5,
                                        size = .5) +
                             theme_bw() +
                             theme(plot.margin = unit(c(0, 0.1, 0, 0.1), "cm"),
                                   strip.background = element_blank(),
                                   strip.text = element_text(size = 12),
                                   axis.title = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.text.y = element_blank(),
                                   axis.text.x = element_text(size = 12))
                           
                           # If not Brown's one-sided then use the same y-axis limits:
                           if (method_i > 1) {
                             hist_plot <- hist_plot +
                               scale_y_continuous(limits = c(0, 3500))
                           }
                           
                           # Next add the titles depending on if it is the rho^2:
                           if (method_types[method_i] == "magma_no_fudge") {
                             hist_plot <- hist_plot +
                               labs(title = TeX('MAGMA: $\\rho^2$'))
                           } else {
                             hist_plot <- hist_plot +
                               labs(title = method_labels[method_i])
                           }
                           
                           # Add x-axis label for bottom:
                           if (method_i == 4) {
                             hist_plot <- hist_plot +
                               labs(x = "Gene-level p-value") +
                               theme(axis.title.x = element_text(size = 20))
                           }
                           
                           # Annotate the first one:
                           if (method_i == 1) {
                             hist_plot <- hist_plot +
                               annotate("text", x = .35, y = 3350,
                                        size = 6,
                                        label = "Corrected null distribution",
                                        color = "#F39B7FFF")
                           }
                           
                           hist_plot +
                             theme(plot.title = element_text(hjust = 0.5, vjust = 0,
                                                             size = 20))
                           
                         })
null_hist_plot_background <-
  plot_grid(plotlist = null_hist_list_upd, ncol = 1, align = "hv")
# Save the individual plot
save_plot("figures/pdf/main/f1c_null_hist_background.pdf",
          null_hist_plot_background, ncol = 1, nrow = 4)
save_plot("figures/nonpdf/main/f1c_null_hist_background.jpg",
          null_hist_plot_background, ncol = 1, nrow = 4)


# Another version - make a summary of the bins across all 1000 simulations:

pval_bin_sim_summary <- gene_mt_sims %>%
  dplyr::select(-vegas) %>%
  mutate(fisher = pmin(1 - 1e-15, fisher)) %>%
  pivot_longer(cols = -sim_group,
               names_to = "method",
               values_to = "pval") %>%
  mutate(pval_group = cut(pval, breaks = seq(0, 1, by = 0.05), right = FALSE)) %>%
  group_by(method, pval_group, sim_group) %>%
  summarize(n_genes = n()) %>%
  ungroup() %>%
  group_by(method, pval_group) %>%
  summarize(mean_genes = mean(n_genes),
            se_genes = sd(n_genes) / sqrt(n()),
            sd_genes = sd(n_genes)) %>%
  mutate(upper_genes = mean_genes + 2 * se_genes,
         lower_genes = mean_genes - 2 * se_genes,
         upper_genes_sd = mean_genes + sd_genes,
         lower_genes_sd = mean_genes - sd_genes) 




# Loop through this to make a list of the plots:
null_hist_list_summary <- lapply(1:length(method_types),
                                 function(method_i) {
                                   bar_plot <- pval_bin_sim_summary %>%
                                     filter(method == method_types[method_i]) %>%
                                     ggplot(aes(x = pval_group)) +
                                     geom_bar(aes(y = mean_genes), stat = "identity",
                                              fill = method_colors[method_i],
                                              alpha = 0.8) +
                                     geom_line(aes(y = upper_genes_sd, group = 1),
                                               color = method_colors[method_i]) +
                                     geom_line(aes(y = lower_genes_sd, group = 1),
                                               color = method_colors[method_i]) +
                                     geom_hline(yintercept = qbinom(0.5, nrow(ex_mt_gene_set), 1 / 20),
                                                color = "black", linetype = "dashed", alpha = .5,
                                                size = 1) +
                                     # geom_hline(yintercept = qbinom(0.975, nrow(ex_mt_gene_set), 1 / 20),
                                     #            color = "darkred", linetype = "dashed", alpha = .5,
                                     #            size = .5) +
                                     theme_bw() +
                                     scale_x_discrete(labels = seq(0, 1, by = 0.05)) +
                                     theme(plot.margin = unit(c(0, 0.1, 0, 0.1), "cm"),
                                           strip.background = element_blank(),
                                           strip.text = element_text(size = 12),
                                           axis.title = element_blank(),
                                           axis.ticks.y = element_blank(),
                                           axis.text.y = element_blank(),
                                           axis.text.x = element_text(size = 8, angle = 45))
                                   
                                   # If not Brown's one-sided then use the same y-axis limits:
                                   if (method_i > 1) {
                                     bar_plot <- bar_plot +
                                       scale_y_continuous(limits = c(0, 3500))
                                   }
                                   
                                   # Next add the titles depending on if it is the rho^2:
                                   if (method_types[method_i] == "magma_no_fudge") {
                                     bar_plot <- bar_plot +
                                       labs(title = TeX('MAGMA: $\\rho^2$'))
                                   } else {
                                     bar_plot <- bar_plot +
                                       labs(title = method_labels[method_i])
                                   }
                                   
                                   # Add x-axis label for bottom:
                                   if (method_types[method_i] == "fisher") {
                                     bar_plot <- bar_plot +
                                       labs(x = "Gene-level p-value") +
                                       theme(axis.title.x = element_text(size = 18))
                                   }
                                   
                                   bar_plot +
                                     theme(plot.title = element_text(hjust = 0.5, vjust = 0,
                                                                     size = 18))
                                   
                                 })

null_hist_plot_summary <-
  plot_grid(plotlist = null_hist_list_summary, ncol = 1, align = "hv")

# Save the individual plot
save_plot("figures/pdf/main/f1c_null_hist_summary.pdf",
          null_hist_plot_summary, ncol = 1, nrow = 5)
save_plot("figures/nonpdf/main/f1c_null_hist_summary.jpg",
          null_hist_plot_summary, ncol = 1, nrow = 5)


# Add the monte carlo to the background:
# Loop through this to make a list of the plots:
null_hist_list_summary_upd <- lapply(1:(length(method_types) - 1),
                                     function(method_i) {
                                       bar_plot <- pval_bin_sim_summary %>%
                                         filter(method == "fisher") %>%
                                         ggplot(aes(x = pval_group)) +
                                         geom_bar(aes(y = mean_genes), stat = "identity",
                                                  fill = method_colors[5],
                                                  alpha = 1) +
                                         geom_line(aes(y = upper_genes_sd, group = 1),
                                                   color = method_colors[5]) +
                                         geom_line(aes(y = lower_genes_sd, group = 1),
                                                   color = method_colors[5]) +
                                         geom_bar(data = filter(pval_bin_sim_summary,
                                                                method == method_types[method_i]),
                                                  aes(y = mean_genes), stat = "identity",
                                                  fill = method_colors[method_i],
                                                  alpha = 0.8) +
                                         geom_line(data = filter(pval_bin_sim_summary,
                                                                 method == method_types[method_i]),
                                                   aes(y = upper_genes_sd, group = 1),
                                                   color = method_colors[method_i]) +
                                         geom_line(data = filter(pval_bin_sim_summary,
                                                                 method == method_types[method_i]),
                                                   aes(y = lower_genes_sd, group = 1),
                                                   color = method_colors[method_i]) +
                                         geom_hline(yintercept = qbinom(0.5, nrow(ex_mt_gene_set), 1 / 20),
                                                    color = "black", linetype = "dashed", alpha = .5,
                                                    size = 1) +
                                         # geom_hline(yintercept = qbinom(0.975, nrow(ex_mt_gene_set), 1 / 20),
                                         #            color = "darkred", linetype = "dashed", alpha = .5,
                                         #            size = .5) +
                                         theme_bw() +
                                         scale_x_discrete(labels = seq(0, 1, by = 0.05)) +
                                         theme(plot.margin = unit(c(0, 0.1, 0, 0.1), "cm"),
                                               strip.background = element_blank(),
                                               strip.text = element_text(size = 12),
                                               axis.title = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               axis.text.x = element_text(size = 8, angle = 45))
                                       
                                       # If not Brown's one-sided then use the same y-axis limits:
                                       if (method_i > 1) {
                                         bar_plot <- bar_plot +
                                           scale_y_continuous(limits = c(0, 3500))
                                       }
                                       
                                       # Next add the titles depending on if it is the rho^2:
                                       if (method_types[method_i] == "magma_no_fudge") {
                                         bar_plot <- bar_plot +
                                           labs(title = TeX('MAGMA: $\\rho^2$'))
                                       } else {
                                         bar_plot <- bar_plot +
                                           labs(title = method_labels[method_i])
                                       }
                                       
                                       # Add x-axis label for bottom:
                                       if (method_i == 4) {
                                         bar_plot <- bar_plot +
                                           labs(x = "Gene-level p-value") +
                                           theme(axis.title.x = element_text(size = 18))
                                       }
                                       
                                       # Annotate the first one:
                                       if (method_i == 1) {
                                         bar_plot <- bar_plot +
                                           annotate("text", x = "[0.3,0.35)", y = 3350,
                                                    size = 6,
                                                    label = "Corrected null distribution",
                                                    color = "#F39B7FFF")
                                       }
                                       
                                       bar_plot +
                                         theme(plot.title = element_text(hjust = 0.5, vjust = 0,
                                                                         size = 18))
                                       
                                     })

null_hist_plot_summary_upd <-
  plot_grid(plotlist = null_hist_list_summary_upd, ncol = 1, align = "hv")

# Save the individual plot
save_plot("figures/pdf/main/f1c_null_hist_summary_background.pdf",
          null_hist_plot_summary_upd, ncol = 1, nrow = 4)
save_plot("figures/nonpdf/main/f1c_null_hist_summary_background.jpg",
          null_hist_plot_summary_upd, ncol = 1, nrow = 4)


# Another option - display each bar from the simulations stacked on top of each
# other with small alpha to visualize the uncertainty


pval_bin_sims <- gene_mt_sims %>%
  dplyr::select(-vegas) %>%
  mutate(fisher = pmin(1 - 1e-15, fisher)) %>%
  pivot_longer(cols = -sim_group,
               names_to = "method",
               values_to = "pval") %>%
  mutate(pval_group = cut(pval, breaks = seq(0, 1, by = 0.05), right = FALSE)) %>%
  group_by(method, pval_group, sim_group) %>%
  summarize(n_genes = n()) %>%
  ungroup()

pval_bin_sims <- pval_bin_sims %>%
  separate(pval_group, c("left_end", "right_end"), sep = ",", remove = FALSE) %>%
  mutate(left_end = as.numeric(str_remove(left_end, "\\[")),
         right_end = as.numeric(str_remove(right_end, "\\)")))

pval_bin_sim_summary <- pval_bin_sim_summary %>%
  separate(pval_group, c("left_end", "right_end"), sep = ",", remove = FALSE) %>%
  mutate(left_end = as.numeric(str_remove(left_end, "\\[")),
         right_end = as.numeric(str_remove(right_end, "\\)")))
  

# Add the monte carlo to the background:
# Loop through this to make a list of the plots:
null_hist_list_all_sims <- lapply(1:(length(method_types) - 1),
                                     function(method_i) {
                                       bar_plot <- pval_bin_sims %>%
                                         filter(method == "fisher") %>%
                                         ggplot() +
                                         geom_segment(#stat = "identity",
                                           aes(x = left_end, xend = right_end,
                                               y = n_genes, yend = n_genes,
                                               group = sim_group),
                                           color = method_colors[5],
                                           #position = "identity",
                                           alpha = 0.05) +
                                         geom_rect(data = filter(pval_bin_sim_summary, method == "fisher"),
                                                   aes(xmin = left_end, xmax = right_end,
                                                       ymin = 0, ymax = mean_genes),
                                                   fill = method_colors[5], alpha = 0.5, 
                                                   color = "white") +
                                         # geom_segment(data = filter(pval_bin_sim_summary, method == "fisher"),
                                         #   aes(x = left_end, xend = right_end,
                                         #       y = mean_genes, yend = mean_genes),
                                         #   color = "white",
                                         #   #position = "identity",
                                         #   alpha = 1) +
                                         geom_segment(data = filter(pval_bin_sims,
                                                                    method == method_types[method_i]),
                                                      aes(x = left_end, xend = right_end,
                                                          y = n_genes, yend = n_genes,
                                                          group = sim_group),
                                                      #stat = "identity",
                                                      #position = "identity",
                                                      color = method_colors[method_i],
                                                      alpha = 0.05) +
                                         geom_rect(data = filter(pval_bin_sim_summary, 
                                                                 method == method_types[method_i]),
                                                   aes(xmin = left_end, xmax = right_end,
                                                       ymin = 0, ymax = mean_genes),
                                                   fill = method_colors[method_i], alpha = 0.5, 
                                                   color = "white") +
                                         geom_hline(yintercept = qbinom(0.025, nrow(ex_mt_gene_set), 1 / 20),
                                                    color = "black", linetype = "dashed", alpha = .5,
                                                    size = 1) +
                                         geom_hline(yintercept = qbinom(0.975, nrow(ex_mt_gene_set), 1 / 20),
                                                    color = "black", linetype = "dashed", alpha = .5,
                                                    size = 1) +
                                         theme_bw() +
                                         #scale_x_discrete(labels = seq(0, 1, by = 0.05)) +
                                         theme(plot.margin = unit(c(0, 0.1, 0, 0.1), "cm"),
                                               strip.background = element_blank(),
                                               strip.text = element_text(size = 12),
                                               axis.title = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               axis.text.x = element_text(size = 8))
                                       
                                       # If not Brown's one-sided then use the same y-axis limits:
                                       if (method_i > 1) {
                                         bar_plot <- bar_plot +
                                           scale_y_continuous(limits = c(0, 3500))
                                       }
                                       
                                       # Next add the titles depending on if it is the rho^2:
                                       if (method_types[method_i] == "magma_no_fudge") {
                                         bar_plot <- bar_plot +
                                           labs(title = TeX('MAGMA: $\\rho^2$'))
                                       } else {
                                         bar_plot <- bar_plot +
                                           labs(title = method_labels[method_i])
                                       }
                                       
                                       # Add x-axis label for bottom:
                                       if (method_i == 4) {
                                         bar_plot <- bar_plot +
                                           labs(x = "Gene-level p-value") +
                                           theme(axis.title.x = element_text(size = 18))
                                       }
                                       
                                       # Annotate the first one:
                                       if (method_i == 1) {
                                         bar_plot <- bar_plot +
                                           annotate("text", x = .3, y = 3350,
                                                    size = 6,
                                                    label = "Corrected null distribution",
                                                    color = "#F39B7FFF")
                                       }
                                       
                                       bar_plot +
                                         theme(plot.title = element_text(hjust = 0.5, vjust = 0,
                                                                         size = 18))
                                       
                                     })

null_hist_plot_all_sims <-
  plot_grid(plotlist = null_hist_list_all_sims, ncol = 1, align = "hv")

# Save the individual plot
save_plot("figures/pdf/main/f1c_null_hist_all_sims.pdf",
          null_hist_plot_all_sims, ncol = 1, nrow = 4)
save_plot("figures/nonpdf/main/f1c_null_hist_all_sims.jpg",
          null_hist_plot_all_sims, ncol = 1, nrow = 4)



# Add the monte carlo to the background:
# Loop through this to make a list of the plots:
null_hist_list_all_ave <- lapply(1:length(method_types),
                                  function(method_i) {
                                    bar_plot <- pval_bin_sim_summary %>%
                                      filter(method == method_types[method_i]) %>%
                                      ggplot(aes(xmin = left_end, xmax = right_end,
                                                 ymin = 0, ymax = mean_genes)) +
                                      geom_rect(fill = method_colors[method_i], 
                                                alpha = 0.8, 
                                                color = method_colors[method_i]) +
                                      geom_segment(y = qbinom(0.5, nrow(ex_mt_gene_set), 1 / 20),
                                                   yend = qbinom(0.5, nrow(ex_mt_gene_set), 1 / 20),
                                                   x = 0, xend = 1,
                                                 color = "black", alpha = 1, size = 1) +
                                      theme_bw() +
                                      #scale_x_discrete(labels = seq(0, 1, by = 0.05)) +
                                      theme(plot.margin = unit(c(0, 0.1, 0, 0.1), "cm"),
                                            strip.background = element_blank(),
                                            strip.text = element_text(size = 12),
                                            axis.title = element_blank(),
                                            axis.ticks.y = element_blank(),
                                            axis.text.y = element_blank(),
                                            axis.text.x = element_text(size = 8))
                                    
                                    # If not Brown's one-sided then use the same y-axis limits:
                                    if (method_i > 1) {
                                      bar_plot <- bar_plot +
                                        scale_y_continuous(limits = c(0, 3500))
                                    }
                                    
                                    # Next add the titles depending on if it is the rho^2:
                                    if (method_types[method_i] == "magma_no_fudge") {
                                      bar_plot <- bar_plot +
                                        labs(title = TeX('MAGMA: $\\rho^2$'))
                                    } else {
                                      bar_plot <- bar_plot +
                                        labs(title = method_labels[method_i])
                                    }
                                    
                                    # Add x-axis label for bottom:
                                    if (method_i == 5) {
                                      bar_plot <- bar_plot +
                                        labs(x = "Gene-level p-value") +
                                        theme(axis.title.x = element_text(size = 18))
                                    }
                                    

                                    bar_plot +
                                      theme(plot.title = element_text(hjust = 0.5, vjust = 0,
                                                                      size = 18))
                                    
                                  })

null_hist_plot_all_ave <-
  plot_grid(plotlist = null_hist_list_all_ave, ncol = 1, align = "hv")

# Save the individual plot
save_plot("figures/pdf/main/f1c_null_hist_all_ave.pdf",
          null_hist_plot_all_ave, ncol = 1, nrow = 5)
save_plot("figures/nonpdf/main/f1c_null_hist_all_ave.jpg",
          null_hist_plot_all_ave, ncol = 1, nrow = 5)


# Fig 1d - FDR control simulation results ---------------------------------


# Now summarize each group's FWER, FDR, # FP by each approach:
gene_mt_results <- gene_mt_sims %>%
  pivot_longer(cols = -sim_group,
               names_to = "method",
               values_to = "pval") %>%
  group_by(sim_group, method) %>%
  mutate(bf_pval = p.adjust(pval, method = "bonferroni"),
         bf_error = as.numeric(bf_pval <= 0.05),
         bh_pval = p.adjust(pval, method = "BH"),
         bh_error = as.numeric(bh_pval <= 0.05)) %>%
  summarize(n_genes = n(),
            obs_fwer = as.numeric(any(bf_error == 1)),
            obs_bf_v = sum(bf_error),
            obs_fdp = as.numeric(any(bh_error == 1)),
            obs_bh_v = sum(bh_error)) %>%
  ungroup()

# Display the FDR control with standard errors:
mt_fdr_sim_plot <- gene_mt_results %>%
  filter(method != "vegas") %>%
  group_by(method) %>%
  summarize(fdr = mean(obs_fdp),
            n_sims = n()) %>%
  mutate(fdp_se = sqrt((fdr * (1- fdr)) /
                         n_sims)) %>%
  mutate(fdp_lower = pmax(0, fdr - 2 * fdp_se),
         fdp_upper = pmin(1, fdr + 2 * fdp_se),
         method = fct_relevel(method,
                              "browns_one_sided", 
                              "magma_no_fudge",
                              "browns_two_sided", 
                              "magma_fudge", 
                              "fisher"),
         method = fct_rev(fct_recode(method,
                                     `MAGMA: paper` = "browns_one_sided",
                                     `MAGMA: rho^2` = "magma_no_fudge",
                                     `MAGMA: code` = "magma_fudge",
                                     `Two-sided approx.` = "browns_two_sided",
                                     `Corrected` = "fisher"))) %>%
  ggplot(aes(x = method, y = fdr)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = fdp_lower, ymax = fdp_upper)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "darkred") +
  theme_bw() + 
  labs(x = "Method", y = "Observed FDP") +
  scale_x_discrete(labels = c(`MAGMA: paper` = "MAGMA: paper",
                              `MAGMA: rho^2` = parse(text = TeX('MAGMA: $\\rho^2$')),
                              `MAGMA: code` = "MAGMA: code",
                              `Two-sided approx.` = "Two-sided approx.",
                              `Corrected` = "Corrected")) +
  coord_flip() + 
  theme(strip.background = element_blank(),
        plot.title = element_text(size = 32),
        plot.subtitle = element_text(size = 24),
        strip.text = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 18))

# Save the individual plot
save_plot("figures/pdf/main/f1d_fdp_sims.pdf",
          mt_fdr_sim_plot, ncol = 1, nrow = 1)
save_plot("figures/nonpdf/main/f1d_fdp_sims.jpg",
          mt_fdr_sim_plot, ncol = 1, nrow = 1)


# Arrange the figures together for a single plot --------------------------

# First using just the single example histogram:
f1_display <- plot_grid(
  plot_grid(null_test_stat_ex_plot_curve, t1_error_sim_plot, mt_fdr_sim_plot,
            ncol = 1, labels = c("a", "b", "d"),
            label_fontface = "bold", 
            label_size = 20),
  null_hist_plot_all_ave, ncol = 2, labels = c("", "c"),
  label_fontface = "bold", label_size = 20)
save_plot("figures/pdf/main/f1_display_1.pdf",
          f1_display, ncol = 3, nrow = 5)
save_plot("figures/nonpdf/main/f1_display_1.jpg",
          f1_display, ncol = 3, nrow = 5)


f1_display_2 <- plot_grid(
  plot_grid(null_test_stat_ex_plot_curve,
            t1_error_sim_plot + facet_wrap(~gene, ncol = 1) +
              scale_y_continuous(limits = c(0.04, .16),
                                 breaks = seq(0.05, .15, by = 0.05)), 
            mt_fdr_sim_plot,
            ncol = 1, labels = c("a", "b", "d"),
            label_fontface = "bold", 
            rel_heights = c(2, 2, 1),
            label_size = 24),
  null_hist_plot_all_ave, ncol = 2, labels = c("", "c"),
  label_fontface = "bold", label_size = 24)
save_plot("figures/pdf/main/f1_display_2.pdf",
          f1_display_2, ncol = 3, nrow = 5)
save_plot("figures/nonpdf/main/f1_display_2.jpg",
          f1_display_2, ncol = 3, nrow = 5)


# Now with the histogram that displays the simulations:
f1_display_3 <- plot_grid(
  plot_grid(null_test_stat_ex_plot_curve,
            t1_error_sim_plot + facet_wrap(~gene, ncol = 1) +
              scale_y_continuous(limits = c(0.04, .16),
                                 breaks = seq(0.05, .15, by = 0.05)), 
            mt_fdr_sim_plot,
            ncol = 1, labels = c("a", "b", "d"),
            label_fontface = "bold", 
            rel_heights = c(2, 2, 1),
            label_size = 20),
  null_hist_plot_all_sims, ncol = 2, labels = c("", "c"),
  label_fontface = "bold", label_size = 20)
save_plot("figures/pdf/main/f1_display_3.pdf",
          f1_display_3, ncol = 3, nrow = 4)
save_plot("figures/nonpdf/main/f1_display_3.jpg",
          f1_display_3, ncol = 3, nrow = 4)

# With the corrected separately
f1_display_4 <- plot_grid(
  plot_grid(null_test_stat_ex_plot_curve,
            t1_error_sim_plot + facet_wrap(~gene, ncol = 1) +
              scale_y_continuous(limits = c(0.04, .16),
                                 breaks = seq(0.05, .15, by = 0.05)), 
            mt_fdr_sim_plot,
            ncol = 1, labels = c("a", "b", "d"),
            label_fontface = "bold", 
            rel_heights = c(2, 2, 1),
            label_size = 20),
  null_hist_plot_all, ncol = 2, labels = c("", "c"),
  label_fontface = "bold", label_size = 20)
save_plot("figures/pdf/main/f1_display_4.pdf",
          f1_display_4, ncol = 3, nrow = 5)
save_plot("figures/nonpdf/main/f1_display_4.jpg",
          f1_display_4, ncol = 3, nrow = 5)




