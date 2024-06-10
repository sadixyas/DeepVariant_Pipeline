setwd("/bigdata/stajichlab/sadikshs/GEN242/challenge_project/SadikshyaSharma_project/varseq")

library(dplyr) #Load library
annotated_vcf_dir <- "results/vcf_anno" # Define the directory containing the annotated VCF files

# List all annotated .tsv files in the directory
annotated_vcf_files <- list.files(path = annotated_vcf_dir, pattern = "_annotated\\.tsv$", full.names = TRUE, recursive = TRUE)

# Check if there are any annotated VCF files found
if (length(annotated_vcf_files) == 0) {
  stop("No annotated .tsv files found in the specified directory.")
}

# Read all annotated files and combine them into one data frame
combined_annotations <- do.call(rbind, lapply(annotated_vcf_files, read.delim))

# Save the combined annotations to a file
write.table(combined_annotations, "results/combineDF_all.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

cat("Combination of annotated VCF files is complete.\n")
