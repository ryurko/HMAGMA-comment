# PURPOSE: Create figure 2 displaying the H-MAGMA results

library(tidyverse)
library(latex2exp)
library(cowplot)


# Load the H-MAGMA data ---------------------------------------------------

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(sheet_i) readxl::read_excel(filename,
                                                           sheet = sheet_i))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  return(x)
}

# Load the adult and fetal brain H-MAGMA results:
hmagma_adult_brain_results <-
  read_excel_allsheets("data/hmagma/output/reported/H-MAGMA_Adult_brain_output.xlsx")
hmagma_fetal_brain_results <-
  read_excel_allsheets("data/hmagma/output/reported/H-MAGMA_Fetal_brain_output.xlsx")


# Get the H-MAGMA fetal discoveries ---------------------------------------

asd_fetal_hmagma_bh_results <-
  hmagma_fetal_brain_results$`H-MAGMA_ASD`$GENE[
    which(
      p.adjust(hmagma_fetal_brain_results$`H-MAGMA_ASD`$P,
               method = "BH") <= 0.05)
  ]
scz_fetal_hmagma_bh_results <-
  hmagma_fetal_brain_results$`H-MAGMA_SCZ`$GENE[
    which(
      p.adjust(hmagma_fetal_brain_results$`H-MAGMA_SCZ`$P,
               method = "BH") <= 0.05)
  ]
adhd_fetal_hmagma_bh_results <-
  hmagma_fetal_brain_results$`H-MAGMA_ADHD`$GENE[
    which(
      p.adjust(hmagma_fetal_brain_results$`H-MAGMA_ADHD`$P,
               method = "BH") <= 0.05)
  ]
bd_fetal_hmagma_bh_results <-
  hmagma_fetal_brain_results$`H-MAGMA_BD`$GENE[
    which(
      p.adjust(hmagma_fetal_brain_results$`H-MAGMA_BD`$P,
               method = "BH") <= 0.05)
  ]
mdd_fetal_hmagma_bh_results <-
  hmagma_fetal_brain_results$`H-MAGMA_MDD`$GENE[
    which(
      p.adjust(hmagma_fetal_brain_results$`H-MAGMA_MDD`$P,
               method = "BH") <= 0.05)
  ]

# Get the H-MAGMA adult discoveries ---------------------------------------

asd_adult_hmagma_bh_results <-
  hmagma_adult_brain_results$`H-MAGMA_ASD`$GENE[
    which(
      p.adjust(hmagma_adult_brain_results$`H-MAGMA_ASD`$P,
               method = "BH") <= 0.05)
  ]
scz_adult_hmagma_bh_results <-
  hmagma_adult_brain_results$`H-MAGMA_SCZ`$GENE[
    which(
      p.adjust(hmagma_adult_brain_results$`H-MAGMA_SCZ`$P,
               method = "BH") <= 0.05)
  ]
adhd_adult_hmagma_bh_results <-
  hmagma_adult_brain_results$`H-MAGMA_ADHD`$GENE[
    which(
      p.adjust(hmagma_adult_brain_results$`H-MAGMA_ADHD`$P,
               method = "BH") <= 0.05)
  ]
bd_adult_hmagma_bh_results <-
  hmagma_adult_brain_results$`H-MAGMA_BD`$GENE[
    which(
      p.adjust(hmagma_adult_brain_results$`H-MAGMA_BD`$P,
               method = "BH") <= 0.05)
  ]
mdd_adult_hmagma_bh_results <-
  hmagma_adult_brain_results$`H-MAGMA_MDD`$GENE[
    which(
      p.adjust(hmagma_adult_brain_results$`H-MAGMA_MDD`$P,
               method = "BH") <= 0.05)
  ]


# Get the H-MAGMA union ---------------------------------------------------

asd_hmagma_bh_union <-
  union(asd_fetal_hmagma_bh_results, asd_adult_hmagma_bh_results)

scz_hmagma_bh_union <-
  union(scz_fetal_hmagma_bh_results, scz_adult_hmagma_bh_results)

adhd_hmagma_bh_union <-
  union(adhd_fetal_hmagma_bh_results, adhd_adult_hmagma_bh_results)

bd_hmagma_bh_union <-
  union(bd_fetal_hmagma_bh_results, bd_adult_hmagma_bh_results)

mdd_hmagma_bh_union <-
  union(mdd_fetal_hmagma_bh_results, mdd_adult_hmagma_bh_results)



# Load the corrected fetal results ----------------------------------------

asd_fetal_mcfisher <-
  map_dfr(list.files("data/hmagma/output/asd/fetal_new/",
                     full.names = TRUE), read_csv) %>%
  # Conservative adjustment to p-values:
  mutate(fisher_pval = (1200000 * fisher_pval + 1) / (1200000))

scz_fetal_mcfisher <-
  map_dfr(list.files("data/hmagma/output/scz/fetal_new/", 
                     full.names = TRUE), read_csv) %>%
  # Conservative adjustment to p-values:
  mutate(fisher_pval = (1200000 * fisher_pval + 1) / (1200000))
adhd_fetal_mcfisher <-
  map_dfr(list.files("data/hmagma/output/adhd/fetal_new/", 
                     full.names = TRUE), read_csv) %>%
  # Conservative adjustment to p-values:
  mutate(fisher_pval = (1200000 * fisher_pval + 1) / (1200000))

mdd_fetal_mcfisher <-
  map_dfr(list.files("data/hmagma/output/mdd/fetal_new/", 
                     full.names = TRUE), read_csv) %>%
  # Conservative adjustment to p-values:
  mutate(fisher_pval = (1200000 * fisher_pval + 1) / (1200000))

bd_fetal_mcfisher <-
  map_dfr(list.files("data/hmagma/output/bd/fetal_new/", 
                     full.names = TRUE), read_csv) %>%
  # Conservative adjustment to p-values:
  mutate(fisher_pval = (1200000 * fisher_pval + 1) / (1200000))



# Load the adult H-MAGMA MC Fisher results --------------------------------

asd_adult_mcfisher <-
  map_dfr(list.files("data/hmagma/output/asd/adult_new/",
                     full.names = TRUE), read_csv) %>%
  # Conservative adjustment to p-values:
  mutate(fisher_pval = (1200000 * fisher_pval + 1) / (1200000))

scz_adult_mcfisher <-
  map_dfr(list.files("data/hmagma/output/scz/adult_new/", 
                     full.names = TRUE), read_csv) %>%
  # Conservative adjustment to p-values:
  mutate(fisher_pval = (1200000 * fisher_pval + 1) / (1200000))
adhd_adult_mcfisher <-
  map_dfr(list.files("data/hmagma/output/adhd/adult_new/", 
                     full.names = TRUE), read_csv) %>%
  # Conservative adjustment to p-values:
  mutate(fisher_pval = (1200000 * fisher_pval + 1) / (1200000))

mdd_adult_mcfisher <-
  map_dfr(list.files("data/hmagma/output/mdd/adult_new/", 
                     full.names = TRUE), read_csv) %>%
  # Conservative adjustment to p-values:
  mutate(fisher_pval = (1200000 * fisher_pval + 1) / (1200000))

bd_adult_mcfisher <-
  map_dfr(list.files("data/hmagma/output/bd/adult_new/", 
                     full.names = TRUE), read_csv) %>%
  # Conservative adjustment to p-values:
  mutate(fisher_pval = (1200000 * fisher_pval + 1) / (1200000))


# Get the fetal H-MAGMA MC Fisher BH results ------------------------------

asd_fetal_mcfisher_bh_results <-
  asd_fetal_mcfisher$gene_id[
    which(p.adjust(asd_fetal_mcfisher$fisher_pval, method = "BH") <= 0.05)
  ]
length(asd_fetal_mcfisher_bh_results)
# [1] 85

scz_fetal_mcfisher_bh_results <-
  scz_fetal_mcfisher$gene_id[
    which(p.adjust(scz_fetal_mcfisher$fisher_pval, method = "BH") <= 0.05)
  ]
length(scz_fetal_mcfisher_bh_results)
# [1] 6062

adhd_fetal_mcfisher_bh_results <-
  adhd_fetal_mcfisher$gene_id[
    which(p.adjust(adhd_fetal_mcfisher$fisher_pval, method = "BH") <= 0.05)
  ]
length(adhd_fetal_mcfisher_bh_results)
# [1] 146

bd_fetal_mcfisher_bh_results <-
  bd_fetal_mcfisher$gene_id[
    which(p.adjust(bd_fetal_mcfisher$fisher_pval, method = "BH") <= 0.05)
  ]
length(bd_fetal_mcfisher_bh_results)
# [1] 792

mdd_fetal_mcfisher_bh_results <-
  mdd_fetal_mcfisher$gene_id[
    which(p.adjust(mdd_fetal_mcfisher$fisher_pval, method = "BH") <= 0.05)
  ]
length(mdd_fetal_mcfisher_bh_results)
# [1] 1580

# Get the adult H-MAGMA MC Fisher BH results ------------------------------

asd_adult_mcfisher_bh_results <-
  asd_adult_mcfisher$gene_id[
    which(p.adjust(asd_adult_mcfisher$fisher_pval, method = "BH") <= 0.05)
  ]
length(asd_adult_mcfisher_bh_results)
# [1] 102

scz_adult_mcfisher_bh_results <-
  scz_adult_mcfisher$gene_id[
    which(p.adjust(scz_adult_mcfisher$fisher_pval, method = "BH") <= 0.05)
  ]
length(scz_adult_mcfisher_bh_results)
# [1] 6469

adhd_adult_mcfisher_bh_results <-
  adhd_adult_mcfisher$gene_id[
    which(p.adjust(adhd_adult_mcfisher$fisher_pval, method = "BH") <= 0.05)
  ]
length(adhd_adult_mcfisher_bh_results)
# [1] 168

bd_adult_mcfisher_bh_results <-
  bd_adult_mcfisher$gene_id[
    which(p.adjust(bd_adult_mcfisher$fisher_pval, method = "BH") <= 0.05)
  ]
length(bd_adult_mcfisher_bh_results)
# [1] 845

mdd_adult_mcfisher_bh_results <-
  mdd_adult_mcfisher$gene_id[
    which(p.adjust(mdd_adult_mcfisher$fisher_pval, method = "BH") <= 0.05)
  ]
length(mdd_adult_mcfisher_bh_results)
# [1] 1664


# Get the union of these sets ---------------------------------------------

asd_mcfisher_bh_union <-
  union(asd_fetal_mcfisher_bh_results, asd_adult_mcfisher_bh_results)

scz_mcfisher_bh_union <-
  union(scz_fetal_mcfisher_bh_results, scz_adult_mcfisher_bh_results)

adhd_mcfisher_bh_union <-
  union(adhd_fetal_mcfisher_bh_results, adhd_adult_mcfisher_bh_results)

bd_mcfisher_bh_union <-
  union(bd_fetal_mcfisher_bh_results, bd_adult_mcfisher_bh_results)

mdd_mcfisher_bh_union <-
  union(mdd_fetal_mcfisher_bh_results, mdd_adult_mcfisher_bh_results)

# Display a plot showing the reduction in discoveries ---------------------

# Make a pie chart version of this

# Make a table where the the counts of which ones that are in HAGMA
# versus in the mcMAGMA are displayed:

hmagma_lost_table <- 
  tibble(disc_set = rep(c("H-MAGMA only", "Corrected"), 5),
         phenotype = c(rep("ASD", 2), rep("SCZ", 2),
                       rep("MDD", 2), rep("ADHD", 2),
                       rep("BD", 2)),
         n_disc = c(length(setdiff(asd_hmagma_bh_union, asd_mcfisher_bh_union)),
                    length(intersect(asd_hmagma_bh_union, asd_mcfisher_bh_union)),
                    length(setdiff(scz_hmagma_bh_union, scz_mcfisher_bh_union)),
                    length(intersect(scz_hmagma_bh_union, scz_mcfisher_bh_union)),
                    length(setdiff(mdd_hmagma_bh_union, mdd_mcfisher_bh_union)),
                    length(intersect(mdd_hmagma_bh_union, mdd_mcfisher_bh_union)),
                    length(setdiff(adhd_hmagma_bh_union, adhd_mcfisher_bh_union)),
                    length(intersect(adhd_hmagma_bh_union, adhd_mcfisher_bh_union)),
                    length(setdiff(bd_hmagma_bh_union, bd_mcfisher_bh_union)),
                    length(intersect(bd_hmagma_bh_union, bd_mcfisher_bh_union))))

hmagma_pie_chart <- hmagma_lost_table %>%
  mutate(phenotype = fct_relevel(phenotype, "ASD"),
         disc_set = fct_relevel(disc_set, "H-MAGMA only")) %>%
  #filter(phenotype %in% c("ASD", "SCZ")) %>%
  group_by(phenotype) %>%
  mutate(disc_prop = n_disc / sum(n_disc),
         pos = cumsum(disc_prop) - disc_prop / 2) %>%
  ungroup() %>%
  ggplot(aes(x = "", y = disc_prop, fill = disc_set)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(disc_prop, 2) * 100, "%"),
                color = disc_set), 
            position = position_stack(vjust = 0.5), size = 5) +
  coord_polar(theta = "y") + 
  facet_wrap(~ phenotype, ncol = 3) +
  scale_fill_manual(values = ggsci::pal_npg()(5)[c(4, 5)]) +
  scale_color_manual(values = c("white", "black"), guide = FALSE) +
  guides(fill = guide_legend(override.aes = list(size = 15))) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        legend.position = c(.85, .3),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_blank())

# Save the individual plot
save_plot("figures/pdf/main/f2b_hmagma_pies.pdf",
          hmagma_pie_chart, ncol = 3, nrow = 2)
save_plot("figures/nonpdf/main/f2b_hmagma_pies.jpg",
          hmagma_pie_chart, ncol = 3, nrow = 2)


# Try proceeding through each phenotype and make the list separately due
# to the strange spacing that is taking place with the facets

pheno_vec <- c("ASD", "ADHD", "BD", "MDD", "SCZ")

pie_list <- lapply(pheno_vec,
                   function(pheno_i) {
                     hmagma_lost_table %>%
                       filter(phenotype == pheno_i) %>%
                       mutate(disc_set = fct_relevel(disc_set, "H-MAGMA only")) %>%
                       mutate(disc_prop = n_disc / sum(n_disc),
                              pos = cumsum(disc_prop) - disc_prop / 2) %>%
                       ggplot(aes(x = "", y = disc_prop, fill = disc_set)) +
                       geom_bar(stat = "identity",
                                aes(color = disc_set, alpha = disc_set)) +
                       geom_text(aes(label = paste0(round(disc_prop, 2) * 100, "%")),
                                 color = "black",
                                 position = position_stack(vjust = 0.5), size = 5) +
                       coord_polar(theta = "y") + 
                       labs(title = pheno_i) +
                       scale_fill_manual(values = ggsci::pal_npg()(5)[c(4, 5)]) +
                       scale_alpha_manual(values = c(.4, .8)) +
                       scale_color_manual(values = ggsci::pal_npg()(5)[c(4, 5)]) +
                       guides(fill = guide_legend(override.aes = list(size = 15))) +
                       labs(fill = "Method", alpha = "Method", color = "Method") +
                       theme_minimal() +
                       theme(axis.text = element_blank(),
                             axis.ticks = element_blank(),
                             panel.grid = element_blank(),
                             axis.title = element_blank(),
                             #legend.position = c(.85, .3),
                             # strip.text = element_text(size = 18),
                             plot.title = element_text(size = 18, hjust = 0.5, vjust = 0),
                             legend.text = element_text(size = 18),
                             legend.title = element_blank())
                   })

display_pie_list <- lapply(pie_list, 
                           function(plot) plot + theme(legend.position = "none",
                                                       plot.margin = unit(c(0,0,0,0), "cm")))

top_row_pie <- plot_grid(display_pie_list[[1]], display_pie_list[[2]], display_pie_list[[3]],
                         ncol = 3)
bottom_row_pie <-
  plot_grid(display_pie_list[[4]], display_pie_list[[5]], get_legend(pie_list[[1]]),
          ncol = 3)

pie_chart_grid <- plot_grid(top_row_pie, bottom_row_pie, ncol = 1)


# Save the individual plot
save_plot("figures/pdf/main/f2b_hmagma_pie_grid.pdf",
          pie_chart_grid, ncol = 3, nrow = 2)
save_plot("figures/nonpdf/main/f2b_hmagma_pie_grid.jpg",
          pie_chart_grid, ncol = 3, nrow = 2)


# Create H-MAGMA versus MC ------------------------------------------------


# Histogram for H-MAGMA ASD:
asd_comp_hist <- 
  bind_rows(hmagma_adult_brain_results$`H-MAGMA_ASD` %>%
              dplyr::select(P) %>%
              rename(asd_pval = P) %>%
              mutate(type = "H-MAGMA"),
            asd_adult_mcfisher %>%
              dplyr::select(fisher_pval) %>%
              rename(asd_pval = fisher_pval) %>%
              mutate(type = "Corrected")) %>%
  #mutate(type = fct_relevel(type, "H-MAGMA")) %>%
  ggplot(aes(x = asd_pval, fill = type, color = type,
# ,
#          color = "white",
#              
             alpha = type)) +
#              color = type, size = type,)) +
  geom_histogram(breaks = seq(0, 1, by = 0.05),
                 position = "identity",
                 closed = "left") +
  #facet_zoom(xlim = c(.5, 1), horizontal = FALSE) +
  scale_alpha_manual(values = c(.8, .4)) +
  #scale_size_manual(values = c(.5, 1)) +
  #ggthemes::scale_fill_colorblind() +
  scale_color_manual(values = ggsci::pal_npg()(5)[c(5, 4)]) +
  scale_fill_manual(values = ggsci::pal_npg()(5)[c(5, 4)]) +
  labs(x = "ASD gene-level p-value",
       y = "Number of genes",
       color = "Method", fill = "Method", alpha = "Method", size = "Method") +
  theme_bw() +
  guides(fill = guide_legend(override.aes = list(size = 15))) +
  theme(#strip.background = element_blank(),
    plot.title = element_text(size = 32),
    plot.subtitle = element_text(size = 24),
    strip.text = element_text(size = 24),
    axis.title.x = element_text(size = 24),
    axis.text.x = element_text(size = 16),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 18),
    legend.position = c(.6, .7),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.direction = "horizontal")

# Save the individual plot
save_plot("figures/pdf/main/f2a_asd_hist.pdf",
          asd_comp_hist, ncol = 1, nrow = 1)
save_plot("figures/nonpdf/main/f2a_asd_hist.jpg",
          asd_comp_hist, ncol = 1, nrow = 1)



# Now the layout:
fig2_layout <-
  plot_grid(asd_comp_hist, pie_chart_grid,
          labels = c("a", "b"), label_fontface = "bold",
          ncol = 2, rel_widths = c(2, 2),
          label_size = 24)
# Save the individual plot
save_plot("figures/pdf/main/f2_hmagma_grid.pdf",
          fig2_layout, ncol = 4, nrow = 2)
save_plot("figures/nonpdf/main/f2_hmagma_grid.jpg",
          fig2_layout, ncol = 4, nrow = 2)


