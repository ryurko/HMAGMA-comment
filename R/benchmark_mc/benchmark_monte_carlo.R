# PURPOSE: Benchmark the run-time for generating Monte Carlo-based null

library(tidyverse)
library(snpcombineR)
library(microbenchmark)

# Load small gene example -------------------------------------------------

# Gene ID: 659

# Load the original reference data:
small_genotype_data <-
  as.matrix(
    data.table::fread(
    "data/reference_data/pruned_genotypes/ex_gene659_chr2_nsnps38_genotype.csv")
  )


# Load medium gene example ------------------------------------------------

# Gene ID: 169355

# Load the original reference data:
medium_genotype_data <- as.matrix(
    data.table::fread(
      "data/reference_data/pruned_genotypes/ex_gene169355_chr8_nsnps101_genotype.csv")
    )


# Load large gene example -------------------------------------------------

# Load the resampled data:
large_genotype_data <- as.matrix(
  data.table::fread(
    "data/reference_data/pruned_genotypes/ex_gene123624_chr15_nsnps1068_genotype.csv")
  )


# Compute correlation matrices --------------------------------------------

compute_cor_benchmark <- microbenchmark(
  "small" = {
    small_cor_matrix <- compute_cor_matrix(small_genotype_data)
  },
  "medium" = {
    medium_cor_matrix <- compute_cor_matrix(medium_genotype_data)
  },
  "large" = {
    large_cor_matrix <- compute_cor_matrix(large_genotype_data)
  },
  times = 1, control = list(order = "inorder")
)

compute_cor_benchmark
# Unit: microseconds
# expr           min         lq       mean     median         uq        max neval
# small      374.595    374.595    374.595    374.595    374.595    374.595     1
# medium    2462.000   2462.000   2462.000   2462.000   2462.000   2462.000     1
# large   251882.981 251882.981 251882.981 251882.981 251882.981 251882.981     1


# Generate null test statistics -------------------------------------------

monte_carlo_null_benchmark <- microbenchmark(
  "small" = {
    small_null_results <- generate_sim_truth_data(small_gene_cor_matrix,
                                                  1000000, 1, 
                                                  .Machine$double.eps)
  },
  "medium" = {
    medium_null_results <- generate_sim_truth_data(medium_cor_matrix,
                                                   1000000, 1, 
                                                   .Machine$double.eps)
  },
  "large" = {
    large_null_results <- generate_sim_truth_data(large_cor_matrix,
                                                  1000000, 1, 
                                                  .Machine$double.eps)
  },
  times = 5, control = list(order = "inorder")
)

# Unit: seconds
# expr       min         lq       mean     median         uq        max neval
# small   3.06932   3.085098   3.094981   3.086144   3.092263   3.142082     5
# medium  10.13034  10.151468  10.218496  10.154909  10.177478  10.478282     5
# large 566.72926 570.248120 574.944150 571.078743 572.218497 594.446128     5

