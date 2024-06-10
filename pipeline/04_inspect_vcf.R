setwd("/bigdata/stajichlab/sadikshs/GEN242/challenge_project/SadikshyaSharma_project/varseq")

library(VariantAnnotation) #Load library

vcf_dir <- "results" # Define the directory containing the VCF files

# Recursively list all .vcf files in the directory
vcf_files <- list.files(path = vcf_dir, pattern = "\\.vcf$", full.names = TRUE, recursive = TRUE)

# Check if there are any VCF files found
if (length(vcf_files) == 0) {
  stop("No .vcf files found in the specified directory.")
}

# Function to inspect a single VCF file
inspect_vcf <- function(vcf_file) {
  cat("Inspecting VCF file:", vcf_file, "\n")
  
  # Read the VCF file
  vcf <- readVcf(vcf_file, "A. thaliana")
  
  # Capture the summary information
  summary_info <- capture.output({
    print(vcf)
    print(as(vcf, "VRanges"))
    summary(vcf)
  })
  
  # Define the output file path
  summary_file <- sub("\\.vcf$", "_summary.txt", vcf_file)
  
  # Write the summary information to the output file
  writeLines(summary_info, con = summary_file)
}

# Loop through each VCF file and inspect it
for (vcf_file in vcf_files) {
  inspect_vcf(vcf_file)
}

cat("Inspection of all VCF files is complete.\n")
