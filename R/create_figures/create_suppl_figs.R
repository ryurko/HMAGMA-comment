# PURPOSE: Create the supplementary figures

library(tidyverse)
library(latex2exp)
library(cowplot)

# -------------------------------------------------------------------------

# Will be using the Nature color palette:
ggsci::pal_npg()(5)
# [1] "#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF"


# Suppl Fig 1 - covariance approximations ---------------------------------

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
                color = "#00A087FF", size = .5) +
  stat_function(fun = function(.x) ifelse(.x >= 0, 
                                          .x * (3.25 + 0.75 * .x),
                                          .x * (3.27 + 0.71 * .x)),
                linetype = "dashed", color = "#E64B35FF",
                size = .5) +
  stat_function(fun = function(.x) (.x)^2 * (3.25 + 0.75 * (.x)^2),
                color = "#4DBBD5FF", linetype = "dashed",
                size = .5) +
  annotate("text", x = 0, y = 2.5, size = 3,
           label = "Two-sided\n approximation",
           color = "#00A087FF") +
  annotate("text", x = 0, y = -1.5, size = 3,
           label = "MAGMA: paper\n(One-sided approximation)",
           color = "#E64B35FF") +
  annotate("text", x = -0.75, y = .5, size = 3,
           label = TeX('MAGMA: $\\rho^2$'), 
           color = "#4DBBD5FF") +
  geom_hline(yintercept = 0, color = "gray", linetype = "dotted") +
  labs(x = TeX('$\\rho$'),
       y = TeX('$\\approx$ covariance')) +
  theme_cowplot() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))


# Save the individual plot
save_plot("figures/pdf/suppl/sf1_cov_approx.pdf",
          cov_approx_line_plot, base_asp = 1)
save_plot("figures/nonpdf/suppl/sf1_cov_approx.jpg",
          cov_approx_line_plot, base_asp = 1)



# Suppl Fig 2 - remaining gene type 1 error plots -------------------------


# Gene: 51306
set.seed(51306)
gene51306_null_sims <-
  map_dfr(c(list.files("data/updated_null_sims/gene51306_1/", full.names = TRUE)),
          read_csv) %>% 
  sample_n(1000000, replace = FALSE) %>%
  # Randomly assign them
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(gene = "24 SNPs",
         fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))

# Gene: 659
set.seed(659)
gene659_null_sims <-
  map_dfr(c(list.files("data/updated_null_sims/gene659_1/", full.names = TRUE)),
          read_csv) %>%  
  sample_n(1000000, replace = FALSE) %>%
  # Randomly assign them
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(gene = "38 SNPs",
         fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))

# Gene: 169355
set.seed(169355)
gene169355_null_sims <-
  map_dfr(list.files("data/updated_null_sims/gene169355/", full.names = TRUE),
          read_csv) %>%
  sample_n(1000000, replace = FALSE) %>%
  # Randomly assign them
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(gene = "101 SNPs",
         fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))

# Gene: 80243
set.seed(80243)
gene80243_null_sims <-
  map_dfr(list.files("data/updated_null_sims/gene80243/", full.names = TRUE),
          read_csv) %>%
  sample_n(1000000, replace = FALSE) %>%
  # Randomly assign them
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(gene = "218 SNPs",
         fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))

# Gene: 26064
set.seed(26064)
gene26064_null_sims <-
  map_dfr(list.files("data/updated_null_sims/gene26064/", full.names = TRUE),
          read_csv) %>%
  sample_n(1000000, replace = FALSE) %>%
  # Randomly assign them
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(gene = "220 SNPs",
         fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))

# Gene: 6263
set.seed(6263)
gene6263_null_sims <-
  map_dfr(list.files("data/updated_null_sims/gene6263/", full.names = TRUE),
          read_csv) %>%
  sample_n(1000000, replace = FALSE) %>%
  # Randomly assign them
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(gene = "772 SNPs",
         fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))

# Gene: 11122
set.seed(11122)
gene11122_null_sims <-
  map_dfr(list.files("data/updated_null_sims/gene11122/", full.names = TRUE),
          read_csv) %>%
  sample_n(1000000, replace = FALSE) %>%
  # Randomly assign them
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(gene = "1,026 SNPs",
         fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))



# Stack these together - selecting only the necessary p-value columns after
# renaming them:
gene_t1_sim_results <- 
  bind_rows(mget(str_subset(ls(), "_null_sims"))) %>%
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

sf_t1_error_sim_plot <- gene_t1_sim_results %>%
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
         gene = fct_relevel(gene, "24 SNPs", "38 SNPs", "101 SNPs", "218 SNPs",
                            "220 SNPs", "772 SNPs", "1,026 SNPs")) %>%
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
  facet_wrap(~gene, ncol = 2) +
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
save_plot("figures/pdf/suppl/sf2_t1_error.pdf",
          sf_t1_error_sim_plot, ncol = 2, nrow = 4)
save_plot("figures/nonpdf/suppl/sf2_t1_error.jpg",
          sf_t1_error_sim_plot, ncol = 2, nrow = 4)


# Suppl Fig 3 - multiple testing FWER plots -------------------------------
rm(list = str_subset(ls(), "null_sims"))


# Following the code in the create_f1_sim_results.R script

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
mt_fwer_sim_plot <- gene_mt_results %>%
  filter(method != "vegas") %>%
  group_by(method) %>%
  summarize(fwer = mean(obs_fwer),
            n_sims = n()) %>%
  mutate(fwer_se = sqrt((fwer * (1 - fwer)) /
                         n_sims)) %>%
  mutate(fwer_lower = pmax(0, fwer - 2 * fwer_se),
         fwer_upper = pmin(1, fwer + 2 * fwer_se),
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
  ggplot(aes(x = method, y = fwer)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = fwer_lower, ymax = fwer_upper)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "darkred") +
  theme_bw() + 
  labs(x = "Method", y = "Observed FWER") +
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
save_plot("figures/pdf/suppl/sf3_fwer_sims.pdf",
          mt_fwer_sim_plot, ncol = 1, nrow = 1)
save_plot("figures/nonpdf/suppl/sf3_fwer_sims.jpg",
          mt_fwer_sim_plot, ncol = 1, nrow = 1)

# Suppl Fig 4 - FWER number of FPs ----------------------------------------


# Create figure comparing the distributions of the number of false positives
# with each method:
library(ggbeeswarm)
all_fp_distr <- gene_mt_results %>%
  filter(method != "vegas") %>%
  mutate(method = fct_relevel(method,
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
  ggplot() +
  geom_quasirandom(aes(x = method, y = obs_bf_v, color = method),
                   alpha = .2) +
  geom_boxplot(aes(x = method, y = obs_bf_v), fill = "white",
               color = "black", alpha = 0.8, width = .2,
               outlier.alpha = 0) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkred") +
  theme_bw() + 
  labs(x = "Method", 
       y = TeX('Number of false positives')) +
  scale_x_discrete(labels = c(`MAGMA: paper` = "MAGMA: paper",
                              `MAGMA: rho^2` = parse(text = TeX('MAGMA: $\\rho^2$')),
                              `Two-sided approx.` = "Two-sided approx.",
                              `MAGMA: code` = "MAGMA: code",
                              `Corrected` = "Corrected")) +
  scale_color_manual(guide = FALSE, 
                     values = rev(ggsci::pal_npg()(5))[1:5]) +
  coord_flip() + 
  theme(axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12),
        panel.grid.major.y = element_blank())


# Next one that displays the proportion of simulations with a number of 
# false positives, but just for the non 'MAGMA: paper' approaches:
fp_facet_prop_points <- gene_mt_results %>%
  filter(!(method %in% c("vegas", "browns_one_sided"))) %>%
  group_by(method, obs_bf_v) %>%
  count() %>%
  ungroup() %>%
  bind_rows(data.frame(method = c(rep("magma_fudge", 3), 
                                     rep("browns_two_sided", 1),
                                     rep("fisher", 5)),
                       obs_bf_v = c(4, 5, 6, 6, 2:6),
                       n = rep(0, 9))) %>%
  mutate(method = fct_relevel(method,
                              "magma_no_fudge",
                              "browns_two_sided", 
                              "magma_fudge", 
                              "fisher"),
         method = fct_recode(method,
                             `MAGMA:~rho^2` = "magma_no_fudge",
                             `Two-sided~approx.` = "browns_two_sided",
                             `MAGMA:~code` = "magma_fudge",
                             `Corrected` = "fisher")) %>%
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
  geom_line(aes(y = prop_sims, group = method),
            color = "gray80", linetype = "dashed", size = .75) +
  geom_point(aes(y = prop_sims, color = method),
             position = position_dodge(width = .5),
             size = 2) +
  geom_errorbar(aes(ymin = prop_lower, ymax = prop_upper,
                    color = method),
                position = position_dodge(width = .5),
                width = .3) +
  labs(color = "Method", y = "Proportion of simulations",
       x = TeX('Number of false positives')) +
  facet_wrap(~ method, ncol = 2, labeller = label_parsed) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = 0:6) +
  scale_color_manual(values = rev(rev(ggsci::pal_npg()(5))[1:4]),
                     guide = FALSE) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.border = element_blank())


fwer_n_fp_plots <-
  plot_grid(fp_facet_prop_points, all_fp_distr,
          ncol = 2, labels = c("a", "b"),
          rel_widths = c(2, 1.25),
          label_fontface = "bold", label_size = 24)


# Save the individual plot
save_plot("figures/pdf/suppl/sf4_mt_fp.pdf",
          fwer_n_fp_plots, ncol = 3, nrow = 2)
save_plot("figures/nonpdf/suppl/sf4_mt_fp.jpg",
          fwer_n_fp_plots, ncol = 3, nrow = 2)


# Suppl Fig 5 to 7 - gene-set analysis ------------------------------------


# Construct the gene-sets of size 45 with 5 from each and 200,000 simulations

# Gene: 3848
set.seed(3848 + 1)
gene3848_null_gs_sims <-
  map_dfr(c(list.files("data/updated_null_sims/gene3848_1/", full.names = TRUE),
            list.files("data/updated_null_sims/gene3848_2/", full.names = TRUE),
            list.files("data/updated_null_sims/gene3848_3/", full.names = TRUE)),
          read_csv) %>%
  sample_n(1000000, replace = FALSE) %>%
  # Randomly assign them
  mutate(sim_group = sample(rep(1:200000, 5))) %>%
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))

# Gene: 51306
set.seed(51306 + 1)
gene51306_null_gs_sims <-
  map_dfr(c(list.files("data/updated_null_sims/gene51306_1/", full.names = TRUE),
            list.files("data/updated_null_sims/gene51306_2/", full.names = TRUE),
            list.files("data/updated_null_sims/gene51306_3/", full.names = TRUE)),
          read_csv) %>% 
  sample_n(1000000, replace = FALSE) %>%
  # Randomly assign them
  mutate(sim_group = sample(rep(1:200000, 5))) %>%
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))

# Gene: 659
set.seed(659 + 1)
gene659_null_gs_sims <-
  map_dfr(c(list.files("data/updated_null_sims/gene659_1/", full.names = TRUE),
            list.files("data/updated_null_sims/gene659_2/", full.names = TRUE),
            list.files("data/updated_null_sims/gene659_3/", full.names = TRUE)),
          read_csv) %>%  
  sample_n(1000000, replace = FALSE) %>%
  # Randomly assign them
  mutate(sim_group = sample(rep(1:200000, 5))) %>%
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))

# Gene: 169355
set.seed(169355 + 1)
gene169355_null_gs_sims <-
  map_dfr(list.files("data/updated_null_sims/gene169355/", full.names = TRUE),
          read_csv) %>%
  sample_n(1000000, replace = FALSE) %>%
  # Randomly assign them
  mutate(sim_group = sample(rep(1:200000, 5))) %>%
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))

# Gene: 80243
set.seed(80243 + 1)
gene80243_null_gs_sims <-
  map_dfr(list.files("data/updated_null_sims/gene80243/", full.names = TRUE),
          read_csv) %>%
  sample_n(1000000, replace = FALSE) %>%
  # Randomly assign them
  mutate(sim_group = sample(rep(1:200000, 5))) %>%
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))

# Gene: 26064
set.seed(26064 + 1)
gene26064_null_gs_sims <-
  map_dfr(list.files("data/updated_null_sims/gene26064/", full.names = TRUE),
          read_csv) %>%
  sample_n(1000000, replace = FALSE) %>%
  # Randomly assign them
  mutate(sim_group = sample(rep(1:200000, 5))) %>%
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))

# Gene: 6263
set.seed(6263 + 1)
gene6263_null_gs_sims <-
  map_dfr(list.files("data/updated_null_sims/gene6263/", full.names = TRUE),
          read_csv) %>%
  sample_n(1000000, replace = FALSE) %>%
  # Randomly assign them
  mutate(sim_group = sample(rep(1:200000, 5))) %>%
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))

# Gene: 11122
set.seed(11122 + 1)
gene11122_null_gs_sims <-
  map_dfr(list.files("data/updated_null_sims/gene11122/", full.names = TRUE),
          read_csv) %>%
  sample_n(1000000, replace = FALSE) %>%
  # Randomly assign them
  mutate(sim_group = sample(rep(1:200000, 5))) %>%
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))

# Gene: 123624
set.seed(123624 + 1)
gene123624_null_gs_sims <-
  map_dfr(list.files("data/updated_null_sims/gene123624/", full.names = TRUE),
          read_csv) %>%
  sample_n(1000000, replace = FALSE) %>%
  # Randomly assign them
  mutate(sim_group = sample(rep(1:200000, 5))) %>%
  rename(magma_fudge = V4, magma_no_fudge = V8,
         browns_two_sided = V12, browns_one_sided = V16,
         fisher = V18, vegas = V20) %>%
  dplyr::select(-contains("V", ignore.case = FALSE)) %>%
  # Make conservative adjustment to fisher and vegas:
  mutate(fisher = (1200000 * fisher + 1) / (1200001),
         vegas = (1200000 * vegas + 1) / (1200001))


# Stack these together - selecting only the necessary p-value columns after
# renaming them:
gene_gs_sims <- bind_rows(mget(str_subset(ls(), "_null_gs_sims"))) 

# Remove the additional large gene sets:
rm(list = str_subset(ls(), "_null_gs_sims"))

# Now compute the type 1 error rates for each of the three different gene-set
# analysis approaches considered:

gs_sim_pvals <- gene_gs_sims %>%
  dplyr::select(-vegas) %>%
  pivot_longer(cols = -sim_group,
               names_to = "method",
               values_to = "pval") %>%
  # Cap the p-values just in case before computing the various stats:
  mutate(pval = pmax(1e-15, pmin(1 - 1e-15, pval)),
         probit_stat = qnorm(1 - pval),
         neg_log_pval = -log(pval)) %>%
  group_by(sim_group, method) %>%
  summarize(fisher_test_stat = 2 * sum(neg_log_pval),
            fisher_df = 2 * n(),
            stouffer_test_stat = sum(probit_stat) / sqrt(n()),
            magma_gs_se = sd(probit_stat) / sqrt(n()),
            magma_gs_test_stat = mean(probit_stat) / magma_gs_se) %>%
  mutate(fisher_gs_pval = pchisq(q = fisher_test_stat,
                                 df = fisher_df, lower.tail = FALSE),
         stouffer_gs_pval = pnorm(stouffer_test_stat, lower.tail = FALSE),
         magma_gs_pval = pnorm(magma_gs_test_stat, lower.tail = FALSE),
         magma_gs_t_pval = pt(magma_gs_test_stat, df = 44, lower.tail = FALSE)) %>%
  ungroup() 

gs_sim_pvals %>%
  filter(method == "fisher") %>%
  ggplot(aes(x = magma_gs_t_pval)) +
  geom_histogram(breaks = seq(0, 1, by = 0.05),
                 closed = "left")

gs_sim_pvals %>%
  filter(method == "fisher") %>%
  ggplot(aes(x = magma_gs_test_stat)) +
  geom_histogram()


gs_sim_results <- gs_sim_pvals %>%
  dplyr::select(sim_group, method, fisher_gs_pval,
                stouffer_gs_pval, magma_gs_t_pval) %>%
  pivot_longer(fisher_gs_pval:magma_gs_t_pval,
               names_to = "gs_method", values_to = "gs_pval") %>%
  mutate(t1_error = ifelse(gs_pval <= .05, 1, 0)) %>%
  group_by(method, gs_method) %>%
  summarize(t1_error_rate = mean(t1_error, na.rm = TRUE),
            n_sims = n()) %>%
  mutate(t1_se = sqrt((t1_error_rate * (1- t1_error_rate)) /
                        n_sims),
         t1_lower = pmax(0, t1_error_rate - 2 * t1_se),
         t1_upper = pmin(1, t1_error_rate + 2 * t1_se))

# Display plot with comparison of each method's type 1 error rate

gs_t1_error_rate_plot <- gs_sim_results %>%
  ungroup() %>%
  mutate(method = fct_relevel(method,
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
                                     `Corrected` = "fisher")),
         gs_method = fct_relevel(gs_method,
                                 "magma_gs_t_pval",
                                 "fisher_gs_pval",
                                 "stouffer_gs_pval"),
         gs_method = fct_recode(gs_method,
                                `MAGMA: gene-set` = "magma_gs_t_pval",
                                `Fisher's` = "fisher_gs_pval",
                                `Stouffer's` = "stouffer_gs_pval")) %>%
  ggplot(aes(x = method, y = t1_error_rate)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = t1_lower, ymax = t1_upper)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "darkred") +
  theme_bw() + 
  labs(x = "Method", y = "Gene-set type 1 error rate") +
  scale_x_discrete(labels = c(`MAGMA: paper` = "MAGMA: paper",
                              `MAGMA: rho^2` = parse(text = TeX('MAGMA: $\\rho^2$')),
                              `MAGMA: code` = "MAGMA: code",
                              `Two-sided approx.` = "Two-sided approx.",
                              `Corrected` = "Corrected")) +
  facet_wrap(~ gs_method, ncol = 3, scales = "free_x") +
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
save_plot("figures/pdf/suppl/sf5_gs_t1_error.pdf",
          gs_t1_error_rate_plot, ncol = 3, nrow = 1)
save_plot("figures/nonpdf/suppl/sf5_gs_t1_error.jpg",
          gs_t1_error_rate_plot, ncol = 3, nrow = 1)

# Next create a figure with the p-value distributions for each of the approaches
gene_set_pval_hist <- gs_sim_pvals %>%
  dplyr::select(sim_group, method, fisher_gs_pval,
                stouffer_gs_pval, magma_gs_t_pval) %>%
  pivot_longer(fisher_gs_pval:magma_gs_t_pval,
               names_to = "gs_method", values_to = "gs_pval") %>%
  mutate(method = fct_relevel(method,
                              "browns_one_sided", 
                              "magma_no_fudge",
                              "browns_two_sided", 
                              "magma_fudge", 
                              "fisher"),
         method = fct_recode(method,
                                     `MAGMA:~paper` = "browns_one_sided",
                                     `MAGMA:~rho^2` = "magma_no_fudge",
                                     `MAGMA:~code` = "magma_fudge",
                                     `Two-sided~approx.` = "browns_two_sided",
                                     `Corrected` = "fisher"),
         gs_method = fct_relevel(gs_method,
                                 "magma_gs_t_pval",
                                 "fisher_gs_pval",
                                 "stouffer_gs_pval"),
         gs_method = fct_recode(gs_method,
                                `MAGMA:~gene-set` = "magma_gs_t_pval",
                                `Fisher~gene-set` = "fisher_gs_pval",
                                `Stouffer~gene-set` = "stouffer_gs_pval")) %>%
  ggplot(aes(x = gs_pval)) +
  geom_histogram(breaks = seq(0, 1, by = 0.05),
                 closed = "left", aes(fill = method), 
                 size = .1, color = "white", alpha = 1,
                 position = "identity") +
  # add expected count lines based on assumption of uniform distribution with 
  # equal probability in each bin:
  geom_hline(yintercept = qbinom(0.025, max(gs_sim_pvals$sim_group), 1 / 20),
             color = "darkred", linetype = "dashed", alpha = 1,
             size = 1) +
  geom_hline(yintercept = qbinom(0.975, max(gs_sim_pvals$sim_group), 1 / 20),
             color = "darkred", linetype = "dashed", alpha = 1,
             size = 1) +
  ggsci::scale_fill_npg() +
  #ggthemes::scale_fill_colorblind() +
  facet_grid(method ~ gs_method, 
             scales = "free_y", labeller = label_parsed) +
  theme_bw() +
  labs(x = "Gene-set analysis p-value", y = "Count") + 
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text.y = element_text(size = 24),
        strip.text.x = element_text(size = 24),
        axis.title = element_text(size = 34),
        axis.text = element_text(size = 18))


# Save the individual plot
save_plot("figures/pdf/suppl/sf6_gs_pval_hist.pdf",
          gene_set_pval_hist, ncol = 3, nrow = 5)
save_plot("figures/nonpdf/suppl/sf6_gs_pval_hist.jpg",
          gene_set_pval_hist, ncol = 3, nrow = 5)

# And finally the MAGMA test stat distribution comparison:

magma_gs_stat_hist_plots <- gs_sim_pvals %>%
  mutate(method = fct_relevel(method,
                              "browns_one_sided", 
                              "magma_no_fudge",
                              "browns_two_sided", 
                              "magma_fudge", 
                              "fisher"),
         method = fct_recode(method,
                             `MAGMA:~paper` = "browns_one_sided",
                             `MAGMA:~rho^2` = "magma_no_fudge",
                             `MAGMA:~code` = "magma_fudge",
                             `Two-sided~approx.` = "browns_two_sided",
                             `Corrected` = "fisher")) %>%
  ggplot(aes(x = magma_gs_test_stat)) +
  geom_histogram(aes(y = ..density..,
                     fill = method), 
                 breaks = seq(-5, 5, by = 0.25),
                 size = .1,
                 color = "white", alpha = 1,
                 position = "identity") +
  stat_function(fun = dt, n = 1000000, 
                args = list(df = 44),
                color = "black", alpha = 1, size = .75) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "black", alpha = 1, size = .75) +
  labs(y = "Density", x = "MAGMA: gene-set test statistic") +
  facet_wrap(~method, labeller = label_parsed,
             ncol = 3) +
  ggsci::scale_fill_npg(guide = FALSE) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 24),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))

# Save these histograms
save_plot("figures/pdf/suppl/sf7_magma_gs_test_stat.pdf",
          magma_gs_stat_hist_plots, ncol = 3, nrow = 2)
save_plot("figures/nonpdf/suppl/sf7_magma_gs_test_stat.jpg",
          magma_gs_stat_hist_plots, ncol = 3, nrow = 2)

