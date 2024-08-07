# require(tidyverse) #for data wrangling
# require(xscss) #for simulation functions
# # Instructions for downloading xscss: https://bioinf.wehi.edu.au/csaw/
# # as part of Lun and Smyth's Paper on CSAW: https://academic.oup.com/nar/article/44/5/e45/2464481?login=false
# require(csaw) #for xscss dependencies, BiocManager::install("csaw")
# require(edgeR) #for xscss dependencies, BiocManager::install("edgeR")
# require(psych) #for geometric mean
# require(truncnorm) #for truncated normal distribution
# require(MAnorm2) #for MAnorm2 normalization functions
# require(DiffBind) #for bulk of norm functions + DB analysis, BiocManager::install("DiffBind")
# require(fuzzyjoin) #for data wrangling of truth to sim data


## these R scripts all live within the repo

source("simulate_reads_proportions.R") #to ChIP-Seq simulate read counts (w/ background)
source("read_analysis.R") #to conduct db analysis on simulated ChIP-Seq read counts
source("dba.normalize.sim.R") #to make DiffBind norm specific to technique
source("MAnorm2.sim.R") #for MAnorm2 normalization and differential binding analysis


big_function <- function(iters, prop.do = 0.05,
                         symmetry = TRUE, equal.DNA.bind = TRUE,
                         constant.background = TRUE, equal.lib.size = TRUE,
                         replicates = 2, fraglen = 100,
                         prior.df = 10, npeaks = 10000,
                         base.mu.jitter = 5,
                         back.const.down, back.const.up,
                         back.vary.down, back.vary.up,
                         prop.consensus = 0.75,
                         FDR.threshold = 0.05) 
{
  
### simulate ChIP-Seq read counts   
  simulate_reads_proportions(iters = iters, prop.do = prop.do, 
                             symmetry = symmetry, 
                             equal.DNA.bind = equal.DNA.bind,
                             constant.background = constant.background,
                             equal.lib.size = equal.lib.size,
                             replicates = replicates, 
                             fraglen = fraglen,
                             prior.df = prior.df, npeaks = npeaks, 
                             base.mu.jitter = base.mu.jitter, 
                             back.const.down = back.const.down,
                             back.const.up = back.const.up, 
                             back.vary.down = back.vary.down,
                             back.vary.up = back.vary.up)
  
### analyze the simulated read counts using DiffBind and MAnorm2
  df_all_final <- read_analysis(iters = iters, prop.do = prop.do,
                                symmetry = symmetry,
                                equal.DNA.bind = equal.DNA.bind,
                                constant.background = constant.background,
                                equal.lib.size = equal.lib.size,
                                replicates = replicates,
                                prop.consensus = prop.consensus,
                                FDR.threshold = FDR.threshold) 
  
## For the sake of memory, we unlink (i.e., delete) each folder after creating the final data frame:
 
     unlink(paste0("./sim_folder_",iters,"_do_",prop.do,"/"), recursive = TRUE)

  return(df_all_final)
}

  
  