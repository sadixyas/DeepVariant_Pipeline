setwd("/bigdata/stajichlab/sadikshs/GEN242/challenge_project/SadikshyaSharma_project/varseq")

library(VariantAnnotation) #Load library

vcf_dir <- "results" # Define the directory containing the VCF files

# Recursively list all .vcf files in the directory
vcf_files <- list.files(path = vcf_dir, pattern = "\\.vcf$", full.names = TRUE, recursive = TRUE)

# Check if there are any VCF files found
if (length(vcf_files) == 0) {
  stop("No .vcf files found in the specified directory.")
}

# Function to filter a single VCF file
filter_vcf <- function(vcf_file, filter_criteria, output_dir) {
  cat("Filtering VCF file:", vcf_file, "\n")
  
  # Read the VCF file
  vcf <- tryCatch({
    readVcf(vcf_file, "A. thaliana")
  }, error = function(e) {
    message("Error reading VCF file: ", vcf_file, "\n", e)
    return(NULL)
  })
  
  if (is.null(vcf)) return(NULL)
  
  # Print summary information for debugging
  cat("VCF dimensions:", dim(vcf), "\n")
  cat("VCF rowRanges length:", length(rowRanges(vcf)), "\n")
  cat("VCF geno fields:", names(geno(vcf)), "\n")
  
  # Convert the VCF object to a VRanges object
  vr <- tryCatch({
    as(vcf, "VRanges")
  }, error = function(e) {
    message("Error converting VCF to VRanges: ", vcf_file, "\n", e)
    return(NULL)
  })
  
  if (is.null(vr)) return(NULL)
  
  # Print VRanges summary for debugging
  cat("VRanges length:", length(vr), "\n")
  
  # Apply the filter criteria
  filtered_vr <- tryCatch({
    subset(vr, eval(parse(text=filter_criteria)))
  }, error = function(e) {
    message("Error applying filter criteria: ", vcf_file, "\n", e)
    return(NULL)
  })
  
  if (is.null(filtered_vr)) return(NULL)
  
  # Ensure filtered_vr has valid indices for subsetting the VCF object
  filtered_indices <- which(start(rowRanges(vcf)) %in% start(filtered_vr) & end(rowRanges(vcf)) %in% end(filtered_vr) & seqnames(rowRanges(vcf)) %in% seqnames(filtered_vr))
  
  if (length(filtered_indices) == 0) {
    message("No variants passed the filter criteria in file: ", vcf_file)
    return(NULL)
  }
  
  # Print indices length for debugging
  cat("Filtered indices length:", length(filtered_indices), "\n")
  cat("Filtered VRanges length:", length(filtered_vr), "\n")
  
  # Create a new VCF object with the filtered variants
  filtered_vcf <- tryCatch({
    vcf[filtered_indices, , drop=FALSE]
  }, error = function(e) {
    message("Error creating filtered VCF object: ", vcf_file, "\n", e)
    return(NULL)
  })
  
  if (is.null(filtered_vcf)) return(NULL)
  
  # Define the output file path
  output_file <- file.path(output_dir, basename(sub("\\.vcf$", "_filtered.vcf", vcf_file)))
  
  # Write the filtered VCF file
  tryCatch({
    writeVcf(filtered_vcf, output_file)
  }, error = function(e) {
    message("Error writing filtered VCF file: ", vcf_file, "\n", e)
    return(NULL)
  })
  
  return(output_file)
}

# Define the filter criteria
filter_criteria <- "totalDepth(vr) >= 20 & (altDepth(vr) / totalDepth(vr) >= 0.8)"

# Define the output directory for filtered VCF files
output_dir <- "results/vcf_filter"

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop through each VCF file and filter it
filtered_vcf_files <- sapply(vcf_files, filter_vcf, filter_criteria=filter_criteria, output_dir=output_dir)

cat("Filtering of all VCF files is complete.\n")

# Check filtering outcome for one sample
check_filter <- function(vcf_file, filtered_vcf_file) {
  original_vr <- as(readVcf(vcf_file, "A. thaliana"), "VRanges")
  filtered_vr <- as(readVcf(filtered_vcf_file, "A. thaliana"), "VRanges")
  
  cat("Original number of variants:", length(original_vr), "\n")
  cat("Filtered number of variants:", length(filtered_vr), "\n")
}

# Example check for the first VCF file
if (length(filtered_vcf_files) > 0 && !is.null(filtered_vcf_files[1])) {
  check_filter(vcf_files[1], filtered_vcf_files[1])
}

library(VariantAnnotation)
vcf_file <- "results/PlantA1/PlantA1..vcf"
vcf <- readVcf(vcf_file, "A. thaliana")
print(vcf)
