## setting 1Gb memory for each core
options('future.globals.maxSize' = 1014*1024^2)

## loading required packages:
library(tidyverse) #for data wrangling
library(xscss) #for ChIP-Seq read count simulations
# Instructions for downloading xscss: https://bioinf.wehi.edu.au/csaw/
# as part of Lun and Smyth's Paper on CSAW: https://academic.oup.com/nar/article/44/5/e45/2464481?login=false
library(csaw) #for xscss dependencies, BiocManager::install("csaw")
library(edgeR) #for xscss dependencies, BiocManager::install("edgeR")
library(psych) #for geometric mean calculations
library(truncnorm) #for truncated normal distribution calculation

library(MAnorm2) #for MAnorm2 normalization and db analysis, install.packages("MAnorm2")
library(DiffBind) #for bulk of norm functions + db analysis, BiocManager::install("DiffBind")
library(fuzzyjoin) #for data wrangling to join results to ground truth via start, end, and chromosome information

## setting up structure:

source("big_function.R")

equal.lib.size <- TRUE
prior.df <- 10
replicates <- 3
prop.consensus <- 0.75
fraglen <- 100
npeaks <- 10000
base.mu.jitter <- 5
back.const.down <- 0.33
back.const.up <- 0.67
back.vary.down <- 0
back.vary.up <- 0.20
FDR.threshold = 0.05

symmetry <- TRUE
equal.DNA.bind <- TRUE
constant.background <- TRUE


library(future)
library(furrr)

slurm_job_user <- Sys.getenv("SLURM_JOB_USER")
slurm_job_id <- Sys.getenv("SLURM_JOB_ID")

setwd(paste0("/tmpfs/",slurm_job_user,"/",slurm_job_id))

plan(multisession, workers = 60) #workers argument sets number of cores used and affects total memory.
# try to lower the workers if you run into memory problems

start_time <- Sys.time()
param_vals <- tidyr::expand_grid(iters = 1:100, 
                                 prop.do = c(seq(0.05, 0.95, by = 0.05))) 

big_possibly = purrr::possibly(.f = big_function, otherwise = NULL)

results_parallel <- param_vals |>
  future_pmap(big_possibly, symmetry = symmetry,
              equal.DNA.bind = equal.DNA.bind,
              constant.background = constant.background,
              equal.lib.size = equal.lib.size, 
              replicates = replicates, prop.consensus = prop.consensus,
              fraglen = fraglen,
              prior.df = prior.df,
              npeaks = npeaks, back.const.down = back.const.down,
              base.mu.jitter = base.mu.jitter,
              back.const.up = back.const.up, back.vary.down = back.vary.down,
              back.vary.up = back.vary.up,
              FDR.threshold = FDR.threshold) |>
  list_rbind()

write_csv(results_parallel, file = paste0("res_symm", symmetry, "_back", 
                                 constant.background,"_DNA",equal.DNA.bind,".csv"))
end_time <- Sys.time()

sink(paste0("time_symm", symmetry, "_back", 
            constant.background, "_DNA",equal.DNA.bind,".csv"))
print(end_time - start_time)
sink()

