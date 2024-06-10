setwd("/bigdata/stajichlab/sadikshs/GEN242/challenge_project/SadikshyaSharma_project/varseq")

# Load necessary libraries
library(VariantAnnotation)
library(GenomicFeatures)
library(GenomicRanges)

filtered_vcf_dir <- "results/vcf_filter" # Define the directory containing the filtered VCF files

# Recursively list all filtered .vcf files in the directory
filtered_vcf_files <- list.files(path = filtered_vcf_dir, pattern = "_filtered\\.vcf$", full.names = TRUE, recursive = TRUE)

# Check if there are any filtered VCF files found
if (length(filtered_vcf_files) == 0) {
  stop("No filtered .vcf files found in the specified directory.")
}

txdb <- loadDb("./data/tair10.sqlite") # Load the annotation database
fa <- FaFile("data/tair10.fasta") # Ensure the fasta file is indexed
indexFa(fa)

# Function to annotate a single VCF file
annotate_vcf <- function(vcf_file, txdb, fa, output_dir) {
  cat("Annotating VCF file:", vcf_file, "\n")
  
  # Read the VCF file
  vcf <- tryCatch({
    readVcf(vcf_file, "A. thaliana")
  }, error = function(e) {
    message("Error reading VCF file: ", vcf_file, "\n", e)
    return(NULL)
  })
  
  if (is.null(vcf)) return(NULL)
  
  # Trim out-of-bound ranges
  vcf <- trim(vcf)
  
  # Annotate variants
  located_variants <- tryCatch({
    locateVariants(vcf, txdb, CodingVariants())
  }, error = function(e) {
    message("Error locating variants: ", vcf_file, "\n", e)
    return(NULL)
  })
  
  if (is.null(located_variants)) return(NULL)
  
  # Print summary of located variants for debugging
  cat("Located variants summary:\n")
  print(summary(located_variants))
  
  # Predict coding changes
  coding_changes <- tryCatch({
    predictCoding(vcf, txdb, seqSource = fa)
  }, error = function(e) {
    message("Error predicting coding changes: ", vcf_file, "\n", e)
    return(NULL)
  })
  
  if (is.null(coding_changes)) return(NULL)
  
  # Print summary of coding changes for debugging
  cat("Coding changes summary:\n")
  print(summary(coding_changes))
  
  # Convert to data frames
  annotated_variants_df <- as.data.frame(located_variants)
  coding_changes_df <- as.data.frame(coding_changes)
  
  # Ensure both data frames have the same length
  if (nrow(annotated_variants_df) != nrow(coding_changes_df)) {
    cat("Error: Mismatch in row counts between annotated_variants_df and coding_changes_df for ", vcf_file, "\n")
    stop("Lengths of annotated_variants_df and coding_changes_df do not match.")
  }
  
  # Combine the annotations into a data frame
  variant_report <- tryCatch({
    data.frame(
      seqnames = annotated_variants_df$seqnames,
      start = annotated_variants_df$start,
      end = annotated_variants_df$end,
      width = annotated_variants_df$width,
      strand = annotated_variants_df$strand,
      variant_type = annotated_variants_df$LOCATION,
      gene = annotated_variants_df$GENEID,
      coding_change = coding_changes_df$CONSEQUENCE
    )
  }, error = function(e) {
    message("Error creating variant report data frame: ", vcf_file, "\n", e)
    return(NULL)
  })
  
  if (is.null(variant_report)) return(NULL)
  
  # Define the output file path
  output_file <- file.path(output_dir, paste0(basename(sub("\\.vcf$", "_annotated", vcf_file)), ".tsv"))
  
  # Save the annotation report
  tryCatch({
    write.table(variant_report, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  }, error = function(e) {
    message("Error writing annotation report: ", vcf_file, "\n", e)
    return(NULL)
  })
  
  return(output_file)
}

# Define the output directory for annotation reports
output_dir <- "results/vcf_anno"

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop through each filtered VCF file and annotate it
annotated_vcf_files <- sapply(filtered_vcf_files, annotate_vcf, txdb = txdb, fa = fa, output_dir = output_dir)

cat("Annotation of all VCF files is complete.\n")

# Example check for the first annotated VCF file
if (length(annotated_vcf_files) > 0 && !is.null(annotated_vcf_files[1])) {
  cat("Example annotation result:\n")
  print(head(read.delim(annotated_vcf_files[1])))
}

