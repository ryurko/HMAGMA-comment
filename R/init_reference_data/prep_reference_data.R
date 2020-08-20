# PURPOSE: Perform pre-processing of the 1000 Genomes reference data.  
#          First determine which SNPs and individuals meet basic criteria 
#          (such as appropriate MAF values, etc). Then will assign the remaining 
#          SNPs to genes using the MAGMA software and prepare pruned genotype
#          to choose from in R/init_data/random_ex_genes.R

# Access packages
library(snpStats)
library(tidyverse)
library(bigsnpr)

# ------------------------------------------------------------------------------

# Load the 1000 Genomes reference data

# Get the file paths for the reference data: -----
g1000_eur_fam_path <- "data/reference_data/g1000_eur/g1000_eur.fam"
g1000_eur_bim_path <- "data/reference_data/g1000_eur/g1000_eur.bim"
g1000_eur_bed_path <- "data/reference_data/g1000_eur/g1000_eur.bed"

# Read in the PLINK data -----
ref_snps_plink <- read.plink(g1000_eur_bed_path, 
                             g1000_eur_bim_path,
                             g1000_eur_fam_path)


# Obtain the SnpMatrix object (genotypes) table from ref_snps_plink list
ref_genotypes <- ref_snps_plink$genotypes
print(ref_genotypes)
# A SnpMatrix with  503 rows and  22665064 columns
# Row names:  HG00096 ... NA20832 
# Col names:  rs537182016 ... rs781880 

#Obtain the SNP information from ref_snps_plink list
ref_genobim <- ref_snps_plink$map
colnames(ref_genobim) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")
print(head(ref_genobim))

# Remove raw file to free up memory
rm(ref_snps_plink)

# ------------------------------------------------------------------------------

# Filter the data based on traditional pre-processing steps

# Create SNP summary statistics (MAF, call rate, etc.)
snp_sum_col <- col.summary(ref_genotypes)

# Set thresholds for both the MAF and call rate
min_call_rate <- 0.95
min_maf <- 0.05

# Filter SNPs on MAF and call rate
use_snps_i <- with(snp_sum_col, 
                   (!is.na(MAF) & MAF >= min_maf) & Call.rate >= min_call_rate)
use_snps_i[is.na(use_snps_i)] <- FALSE                # Remove NA's as well

cat(ncol(ref_genotypes) - sum(use_snps_i),
    "SNPs will be removed due to low MAF or call rate.\n") 
# 16445832 SNPs will be removed due to low MAF or call rate.
rm(snp_sum_col)

# Subset genotype and SNP summary data for SNPs that pass call rate and MAF criteria
ref_genotypes <- ref_genotypes[, use_snps_i]
rm(use_snps_i)

# Next filter the individuals based on their call rate (analogous to SNP call rate)
min_sample_call_rate <- 0.95 
snp_sum_row <- row.summary(ref_genotypes)
sample_use_i <- with(snp_sum_row, !is.na(Call.rate) & Call.rate >= min_sample_call_rate)
sample_use_i[is.na(sample_use_i)] <- FALSE   
cat(nrow(ref_genotypes) - sum(sample_use_i), 
    "subjects will be removed due to low sample call rate.\n")
# Oh good that is nice to see - and not surprising since this is reference data

# Now filter ref_genobim 
ref_genobim <- ref_genobim[which(ref_genobim$SNP %in% colnames(ref_genotypes)),]

# Save a version of this file that is tab delimited, with no header and contains
# only the SNP ID, CHR, and BP (in that order):
readr::write_tsv(ref_genobim[, c("SNP", "chr", "position")],
                 path = "data/reference_data/g1000_eur/ref_data_snps.txt",
                 col_names = FALSE)

# The following MAGMA command is used to generate annotations with 10kb padding 
# upstream and downstream
# (NOTE: need to replace the file paths with your correct ones to generate the 
#        annotations in their respective paths, since this is not a R command
#        but rather expected to be copied and pasted in a command-line terminal:
# magma --annotate window=10 --snp-loc /data/reference_data/g1000_eur/ref_data_snps.txt --gene-loc /data/reference_data/NCBI37.3/NCBI37.3.gene.loc --out /data/reference_data/g1000_eur/ref_data_genes_10kb_pad

# 6219232 SNP locations read from file
# of those, 2893021 (46.52%) mapped to at least one gene
# Writing annotation to file /data/reference_data/g1000_eur/ref_data_genes_10kb_pad.genes.annot
# for chromosome  1, 7 genes are empty (out of 2016)
# for chromosome  2, 5 genes are empty (out of 1226)
# for chromosome  8, 1 gene is empty (out of 669)
# for chromosome  9, 3 genes are empty (out of 775)
# for chromosome 10, 3 genes are empty (out of 723)
# for chromosome 12, 2 genes are empty (out of 1009)
# for chromosome 15, 4 genes are empty (out of 586)
# for chromosome 16, 5 genes are empty (out of 817)
# for chromosome 17, 2 genes are empty (out of 1147)
# for chromosome 19, 1 gene is empty (out of 1389)
# for chromosome  X, 21 genes are empty (out of 805)
# for chromosome  Y, 47 genes are empty (out of 47)
# at least one SNP mapped to each of a total of 19326 genes (out of 19427)
# Only need to focus on the 2,893,021 SNPs

# ------------------------------------------------------------------------------

# Load the protein coding gene annotations generated using MAGMA with 10kb padding 
ref_gene_annotation_data <- read_lines("data/reference_data/g1000_eur/ref_data_genes_10kb_pad.genes.annot",
                                       skip = 2) # skip the first two lines without genes
# Create a dataset to use -------
ref_gene_info_data <- purrr::map_dfr(ref_gene_annotation_data,
                                     function(gene) {
                                       gene_info <- unlist(str_split(gene, "\\t"))
                                       # Extract the SNPs:
                                       gene_snps <- gene_info[3:length(gene_info)]
                                       gene_snps <- gene_snps[str_detect(gene_snps, "rs")]
                                       # Extract the CHR:
                                       gene_chr <- unlist(str_split(gene_info[2], ":"))[1]
                                       data.frame("gene_id" = gene_info[1],
                                                  "gene_chr" = gene_chr,
                                                  "gene_n_snps" = length(gene_snps),
                                                  "gene_snps" = paste0(gene_snps, collapse = "_"))
                                     }) 

# Remove genes with 0 SNPs: ------
ref_gene_info_data <- ref_gene_info_data %>%
  filter(gene_n_snps > 0)
# Remove the annotation data -----
rm(ref_gene_annotation_data)

# What does this distribution look like?
hist(ref_gene_info_data$gene_n_snps)
# Incredibly skewed
table(ref_gene_info_data$gene_chr)
# 1   10   11   12   13   14   15   16   17   18   19    2   20   21   22    3    4    5    6    7 
# 2009  720 1275 1007  320  595  582  812 1145  271 1388 1221  527  215  442 1050  745  856 1016  906 
# 8    9    X 
# 668  772  784
summary(ref_gene_info_data$gene_n_snps)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.0    47.0    86.0   166.4   167.0 12911.0 


# Consider genes of sizes between 9, 100, 500+, then group by the size group and
# the chromosome, randomly select five from each::
set.seed(1389)
example_genes_info_data <- ref_gene_info_data %>%
  # Remove the sex chromosome and remove genes with less than 10 SNPs
  filter(gene_chr != "X", gene_n_snps > 9) %>%
  # Need groups to reflect gene sizes - don't want to use quantiles exactly
  # since there really isn't a need to vary the group sizes between 1-30 and 
  # then 30-40 and so on. Just use 10-100, 100-500, and then 500 and up
  mutate(gene_size_group = cut(gene_n_snps, breaks = c(9, 100, 500, Inf))) %>%
  # Next sample by both the chromosome and gene size:
  group_by(gene_chr, gene_size_group) %>%
  # For now pick 5 from each group, chromosome and gene size group:
  slice(sample(1:n(), 5)) %>%
  ungroup()
rm(ref_gene_info_data)

# Extract the unique SNPs from this dataset and only use these for the genotype data
unique_ex_snps <- purrr::map({(
  example_genes_info_data %>%
    dplyr::pull(gene_snps))}, 
  function(snps) {
    unlist(stringr::str_split(snps, pattern = "_"))
  }) %>%
  unlist() %>%
  unique()  
ref_genotypes <- ref_genotypes[, which(colnames(ref_genotypes) %in%
                                         unique_ex_snps)]
# Remove the vector:
rm(unique_ex_snps)

# Create the numeric version of this matrix with alleles encoded as 0, 1, or 2:
numeric_ref_genotypes <- as(ref_genotypes, "numeric")
# Flip it:
flip_matrix <- function(x) {
  zero2 <- which(x == 0)
  two0 <- which(x == 2)
  x[zero2] <- 2
  x[two0] <- 0
  return(x)
}
numeric_ref_genotypes <- flip_matrix(numeric_ref_genotypes)

# Now will proceed through each of these genes, saving pruned genotype data to use
# A for loop should suffice here - pruned versions using bigsnpr
for (gene_i in 1:nrow(example_genes_info_data)) {
  # Gene id:
  gene_id <- example_genes_info_data$gene_id[gene_i]
  print(paste0("Starting gene: ", gene_id))
  
  # What's the chromosome?
  gene_chr <- example_genes_info_data$gene_chr[gene_i]
  
  # What are these gene's SNPs?
  gene_snps <- unlist(str_split(example_genes_info_data$gene_snps[gene_i], "_"))
  
  # Only look at these SNPs' genotype data:
  gene_ref_genotypes <- numeric_ref_genotypes[, 
                                              which(colnames(numeric_ref_genotypes) %in%
                                                      gene_snps)]
  # Remove the row names:
  row.names(gene_ref_genotypes) <- NULL
  
  # Turn this is into a fake SNP dataset using bigsnpr to perform LD clumping 
  # with MAF as the variable of importancw:
  fake_gene_ref <- snp_fake(nrow(gene_ref_genotypes), 
                            ncol(gene_ref_genotypes))
  fake_gene_ref$genotypes[] <- gene_ref_genotypes
  
  fake_snp_keep_i <- snp_clumping(fake_gene_ref$genotypes,
                                  infos.chr = fake_gene_ref$map$chromosome, 
                                  thr.r2 = 0.95, # Removes extreme LD values
                                  size = ncol(gene_ref_genotypes))
  
  # Save the pruned matrix with a name including the chromosome it is from:
  write_csv(as.data.frame(gene_ref_genotypes[, fake_snp_keep_i]),
            paste0("data/reference_data/pruned_genotype/ex_gene", 
                   gene_id,
                   "_chr", 
                   gene_chr, "_nsnps", length(fake_snp_keep_i), "_genotype.csv"))
  
}

