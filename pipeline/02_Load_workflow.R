setwd("/bigdata/stajichlab/sadikshs/GEN242/challenge_project/SadikshyaSharma_project/varseq")
library(systemPipeR)
sal <- SPRproject()
sal <- importWF(sal, file_path = "systemPipeVARseq.Rmd")



stepsWF(sal)

resources <- list(
  conffile = ".batchtools.conf.R",
  template = "batchtools.slurm.tmpl",
  Njobs = 18,
  walltime = 180,
  ntasks = 1,
  ncpus = 4,
  memory = 4096,
  partition = "gen242"
)

sal <- addResources(sal, step = c("preprocessingg", "trimming", "hisat2_mapping"), resources = resources)

sal <- runWF(sal)

sal <- renderReport(sal)
