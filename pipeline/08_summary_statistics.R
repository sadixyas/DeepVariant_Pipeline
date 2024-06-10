setwd("/bigdata/stajichlab/sadikshs/GEN242/challenge_project/SadikshyaSharma_project/varseq")

library(dplyr) # Load library
combined_file <- "results/combineDF_all.tsv" # Define the path to the combined annotation file
combined_df <- read.delim(combined_file) # Read the combined file

# Generate summary statistics
summary_stats <- combined_df %>%
  group_by(coding_change) %>%
  summarise(count = n())

# Save the summary statistics to a file
write.table(summary_stats, "results/variantStats_all.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

cat("Summary statistics of variants are generated.\n")

