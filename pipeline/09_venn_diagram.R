setwd("/bigdata/stajichlab/sadikshs/GEN242/challenge_project/SadikshyaSharma_project/varseq")

# Load libraries
library(systemPipeR)
library(VennDiagram)


annotated_vcf_dir <- "results/vcf_anno" # Define the directory containing the annotated VCF files

annotated_vcf_files <- list.files(path = annotated_vcf_dir, pattern = "_annotated\\.tsv$", full.names = TRUE, recursive = TRUE) # List all annotated .tsv files in the directory

# Check if there are any annotated VCF files found
if (length(annotated_vcf_files) == 0) {
  stop("No annotated .tsv files found in the specified directory.")
}

# Make a list of variants from the first three samples
varlist <- lapply(annotated_vcf_files[1:3], function(file) {
  as.character(read.delim(file)$seqnames)
})

names(varlist) <- basename(annotated_vcf_files[1:3])

# Create Venn diagram using VennDiagram package
venn.plot <- venn.diagram(
  x = varlist,
  category.names = names(varlist),
  filename = "results/vennplot_var.png",
  output = TRUE,
  imagetype = "png",
  height = 1000,  # Increased height
  width = 1000,   # Increased width
  resolution = 300,
  compression = "lzw",
  col = c("red", "blue", "green"),
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cex = 0.5, # Text size
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.5, # Category text size
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cat.pos = c(-30, 30, 180), # Label positions
  cat.dist = c(0.1, 0.1, 0.1), # Increased distance from circles
  cat.just = list(c(0.5, 0.5), c(0.5, 0.5), c(0.5, 0.5)),
  main = "Venn Plot of First 3 Samples",
  main.cex = 0.8, # Main title text size
  margin = 0.2 # Increased margin
)

cat("Venn diagram created at results/vennplot_var.png.\n")
