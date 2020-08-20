# PURPOSE: Initialize resampled genotype matrices to use for the null gene
#          simulations

library(tidyverse)

# Gene 3848 ---------------------------------------------------------------

# Load the pruned genotype data
gene3848_pruned_genotype <-  as.matrix(
  data.table::fread("data/reference_data/pruned_genotypes/ex_gene3848_chr12_nsnps17_genotype.csv"))

# Resample 5000 times:
set.seed(3848)
gene3848_resampled_genotype <- 
  gene3848_pruned_genotype[sample(1:nrow(gene3848_pruned_genotype),
                                  size = 5000, replace = TRUE),]

# Save it:
write_csv(as_tibble(gene3848_resampled_genotype),
          "data/reference_data/resampled_genotypes/ex_gene3848_chr12_nsnps17_genotype.csv")

rm(list = ls())

# Gene 51306 --------------------------------------------------------------

# Load the pruned genotype data
gene51306_pruned_genotype <-  as.matrix(
  data.table::fread("data/reference_data/pruned_genotypes/ex_gene51306_chr5_nsnps24_genotype.csv"))

# Resample 5000 times:
set.seed(51306)
gene51306_resampled_genotype <- 
  gene51306_pruned_genotype[sample(1:nrow(gene51306_pruned_genotype),
                                  size = 5000, replace = TRUE),]

# Save it:
write_csv(as_tibble(gene51306_resampled_genotype),
          "data/reference_data/resampled_genotypes/ex_gene51306_chr5_nsnps24_genotype.csv")

rm(list = ls())


# Gene 659 ----------------------------------------------------------------

# Load the pruned genotype data
gene659_pruned_genotype <-  as.matrix(
  data.table::fread("data/reference_data/pruned_genotypes/ex_gene659_chr2_nsnps38_genotype.csv"))

# Resample 5000 times:
set.seed(659)
gene659_resampled_genotype <- 
  gene659_pruned_genotype[sample(1:nrow(gene659_pruned_genotype),
                                   size = 5000, replace = TRUE),]

# Save it:
write_csv(as_tibble(gene659_resampled_genotype),
          "data/reference_data/resampled_genotypes/ex_gene659_chr2_nsnps38_genotype.csv")

rm(list = ls())


# Gene 169355 -------------------------------------------------------------

# Load the pruned genotype data
gene169355_pruned_genotype <-  as.matrix(
  data.table::fread("data/reference_data/pruned_genotypes/ex_gene169355_chr8_nsnps101_genotype.csv"))

# Resample 5000 times:
set.seed(169355)
gene169355_resampled_genotype <- 
  gene169355_pruned_genotype[sample(1:nrow(gene169355_pruned_genotype),
                                 size = 5000, replace = TRUE),]

# Save it:
write_csv(as_tibble(gene169355_resampled_genotype),
          "data/reference_data/resampled_genotypes/ex_gene169355_chr8_nsnps101_genotype.csv")

rm(list = ls())


# Gene 80243 --------------------------------------------------------------

# Load the pruned genotype data
gene80243_pruned_genotype <-  as.matrix(
  data.table::fread("data/reference_data/pruned_genotypes/ex_gene80243_chr8_nsnps218_genotype.csv"))

# Resample 5000 times:
set.seed(80243)
gene80243_resampled_genotype <- 
  gene80243_pruned_genotype[sample(1:nrow(gene80243_pruned_genotype),
                                    size = 5000, replace = TRUE),]

# Save it:
write_csv(as_tibble(gene80243_resampled_genotype),
          "data/reference_data/resampled_genotypes/ex_gene80243_chr8_nsnps218_genotype.csv")

rm(list = ls())


# Gene 26064 --------------------------------------------------------------

# Load the pruned genotype data
gene26064_pruned_genotype <-  as.matrix(
  data.table::fread("data/reference_data/pruned_genotypes/ex_gene26064_chr5_nsnps220_genotype.csv"))

# Resample 5000 times:
set.seed(26064)
gene26064_resampled_genotype <- 
  gene26064_pruned_genotype[sample(1:nrow(gene26064_pruned_genotype),
                                   size = 5000, replace = TRUE),]

# Save it:
write_csv(as_tibble(gene26064_resampled_genotype),
          "data/reference_data/resampled_genotypes/ex_gene26064_chr5_nsnps220_genotype.csv")

rm(list = ls())


# Gene 6263 ---------------------------------------------------------------

# Load the pruned genotype data
gene6263_pruned_genotype <-  as.matrix(
  data.table::fread("data/reference_data/pruned_genotypes/ex_gene6263_chr15_nsnps772_genotype.csv"))

# Resample 5000 times:
set.seed(6263)
gene6263_resampled_genotype <- 
  gene6263_pruned_genotype[sample(1:nrow(gene6263_pruned_genotype),
                                   size = 5000, replace = TRUE),]

# Save it:
write_csv(as_tibble(gene6263_resampled_genotype),
          "data/reference_data/resampled_genotypes/ex_gene6263_chr15_nsnps772_genotype.csv")

rm(list = ls())

# Gene 11122 --------------------------------------------------------------

# Load the pruned genotype data
gene11122_pruned_genotype <-  as.matrix(
  data.table::fread("data/reference_data/pruned_genotypes/ex_gene11122_chr20_nsnps1026_genotype.csv"))

# Resample 5000 times:
set.seed(11122)
gene11122_resampled_genotype <- 
  gene11122_pruned_genotype[sample(1:nrow(gene11122_pruned_genotype),
                                  size = 5000, replace = TRUE),]

# Save it:
write_csv(as_tibble(gene11122_resampled_genotype),
          "data/reference_data/resampled_genotypes/ex_gene11122_chr20_nsnps1026_genotype.csv")

rm(list = ls())


# Gene 123624 -------------------------------------------------------------

# Load the pruned genotype data
gene123624_pruned_genotype <-  as.matrix(
  data.table::fread("data/reference_data/pruned_genotypes/ex_gene123624_chr15_nsnps1068_genotype.csv"))

# Resample 5000 times:
set.seed(123624)
gene123624_resampled_genotype <- 
  gene123624_pruned_genotype[sample(1:nrow(gene123624_pruned_genotype),
                                   size = 5000, replace = TRUE),]

# Save it:
write_csv(as_tibble(gene123624_resampled_genotype),
          "data/reference_data/resampled_genotypes/ex_gene123624_chr15_nsnps1068_genotype.csv")

rm(list = ls())





