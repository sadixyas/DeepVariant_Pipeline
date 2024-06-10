setwd("/bigdata/stajichlab/sadikshs/GEN242/challenge_project/SadikshyaSharma_project/varseq")
list.files("data")

#To load necessary libraries and modules
library(systemPipeR)
moduleload("sratoolkit/3.0.0")
system("vdb-config --prefetch-to-cwd")

#To define SRA IDs for my project 
sraidv <- paste("SRR1051", 389:415, sep="")

#Defining download function
getSRAfastq <- function(sraid, threads=1) {                                                                                                                                         
  system(paste("prefetch", sraid)) # makes download faster                                                                                                                        
  system(paste("vdb-validate", sraid)) # checks integrity of the downloaded SRA file                                                                   
  system(paste("fastq-dump --split-files --gzip", sraid)) # gzip option makes it slower but saves storage space                                                                   
  unlink(x=sraid, recursive = TRUE, force = TRUE) # deletes sra download directory                                                                                                
}  

#Run the download
mydir <- getwd(); setwd("data")
for(i in sraidv) getSRAfastq(sraid=i)
setwd(mydir)

#download reference genome and annotation
downloadRefs <- function(rerun=FALSE) {
  if(rerun==TRUE) {
    library(Biostrings)
    download.file("https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz", "./data/tair10.fasta.gz")
    R.utils::gunzip("./data/tair10.fasta.gz")
    download.file("https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.59.gff3.gz", "./data/tair10.gff.gz")
    R.utils::gunzip("./data/tair10.gff.gz")
    txdb <- GenomicFeatures::makeTxDbFromGFF(file = "data/tair10.gff", format = "gff", dataSource = "TAIR", organism = "Arabidopsis thaliana")
    AnnotationDbi::saveDb(txdb, file="./data/tair10.sqlite")
    download.file("https://cluster.hpcc.ucr.edu/~tgirke/Teaching/GEN242/data/tair10_functional_descriptions", "./data/tair10_functional_descriptions")
  }
}


#Execute reference genome download
downloadRefs(rerun=TRUE) 
