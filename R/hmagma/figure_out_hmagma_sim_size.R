library(tidyverse)

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


# Check out how many genes of the different sizes are in each of the groups:
hmagma_asd_gene_size_distr <- hmagma_adult_brain_results$`H-MAGMA_ASD` %>%
  mutate(gene_size_group = cut(NSNPS, breaks = c(-Inf, 100, 500, Inf))) %>%
  group_by(gene_size_group) %>%
  count() %>%
  ungroup() %>%
  mutate(prop = n / sum(n))
#     gene_size_group     n    prop
#     <fct>           <int>   <dbl>
#   1 (-Inf,100]      43741 0.826
#   2 (100,500]        8681 0.164
#   3 (500, Inf]        505 0.00954

# ASD results are roughly 53K genes - round these to 83, 16, and 1
hmagma_asd_gene_size_distr %>%
  mutate(sim_n = round(prop, digits = 2) * 53000,
         sim_n_per_gene = sim_n / 3,
         sim_n_gene_total = sim_n_per_gene * 1000)
# A tibble: 3 x 6
#     gene_size_group     n    prop sim_n sim_n_per_gene sim_n_gene_total
#     <fct>           <int>   <dbl> <dbl>          <dbl>            <dbl>
#   1 (-Inf,100]      43741 0.826   43990         14663.        14663333.
#   2 (100,500]        8681 0.164    8480          2827.         2826667.
#   3 (500, Inf]        505 0.00954   530           177.          176667.

# Okay - so could do 15 million sims for each of the smaller size genes,
# 3 million for the mid-size, and then only 180K for the largest (done)


