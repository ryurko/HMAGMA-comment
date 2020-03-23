# PURPOSE: Randomly pick genes among the pruned genes, stratified by the number
#          of SNPs (post-pruning generated via prep_reference_data.R)

# Access tidyverse
library(tidyverse)

# ------------------------------------------------------------------------------

# First get the example gene genotype data file paths:
ex_pruned_gene_files <- list.files("data/reference_data/pruned_genotype", 
                                   full.names = TRUE)

ex_pruned_gene_info <- map_dfr(ex_pruned_gene_files,
                               function(sim_file) {
                                 # Get the file info summary:
                                 file_info <- str_remove(sim_file,
                                                         "data/reference_data/pruned_genotype/") %>%
                                   str_remove("_genotype\\.csv")
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
                                 
                                 data.frame("gene_id" = gene_name,
                                            "chr" = gene_chr,
                                            "n_snps" = gene_n_snps)
                               })
set.seed(1389)
sample_ex_pruned_gene_info <- ex_pruned_gene_info %>%
  # Remove genes with less than 10 SNPs
  filter(n_snps > 9) %>%
  # Need groups to reflect gene sizes - 10-100, 100-500, and then 500 and up
  mutate(gene_size_group = cut(n_snps, breaks = c(9, 100, 500, Inf))) %>%
  # Next sample by gene size:
  group_by(gene_size_group) %>%
  sample_n(3) %>%
  ungroup()

#     gene_id   chr n_snps gene_size_group
#     <chr>   <dbl>  <dbl> <fct>          
#   1 51306       5     24 (9,100]        
#   2 3848       12     17 (9,100]        
#   3 659         2     38 (9,100]        
#   4 80243       8    218 (100,500]      
#   5 169355      8    101 (100,500]      
#   6 26064       5    220 (100,500]      
#   7 6263       15    772 (500,Inf]      
#   8 123624     15   1068 (500,Inf]      
#   9 11122      20   1026 (500,Inf] 

# Will then use the following vector to select the genes for generating the 
# simulated data:
paste0("(gene", paste0(sample_ex_pruned_gene_info$gene_id, collapse = "_)(gene"), "_)")
#[1] "(gene51306_)|(gene3848_)|(gene659_)|(gene80243_)|(gene169355_)|(gene26064_)|(gene6263_)|(gene123624_)|(gene11122_)"
