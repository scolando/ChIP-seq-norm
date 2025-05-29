# The code used to simulate ChIP-Seq reads borrows substantially from
# the code created by Lun and Smyth (2016) and Evans et al. (2016)

# require(tidyverse)
# require(xscss)
# require(csaw) 
# require(edgeR) 
# require(psych)
# require(truncnorm) 

source("addBackground_new.R")

# iters: integer; the iteration number used to index folders and files
# prop.do: numeric; proportion of peaks that have differential DNA occupancy across A and B
# symmetry: logical; indicates whether the same proportion of peaks have more DNA occupancy across A and B
# constant.background: logical; indicates whether the background binding **between** A and B for the same genomic region is constant
# equal.DNA.bind: logical; indicates whether there is the same amount of DNA Occupancy in the simulated peaks across A and B
# equal.lib.size: library size is equal across the two experimental states in the called regions
#      TRUE = 600,000, FALSE = 1,032,000 for A and 602,000 for B 
#      logic: keep above 600,000 so that the peak reads are substantially above the background reads
# replicates: integer; the number of .bam files generated within each experimental state (i.e., for A and B)
# fraglen: integer; half of the width of the peak
# npeaks: integer; number of simulated peaks in each replicate (some of which have differential DNA occupancy, some of which don't)
# prior.df: integer; connected to dispersion, used in calculating reads from means using negative binomial distribution
# true.width: artifact of original simulation code created by Lun and Smyth (2016)
#      parameter is only relevant in their histone marker simulation -- however, we are simulating transcription factors
# base.mu.jitter: integer; when sampling the base mu values, we `sample(c(base.mu_A-base.mu.jitter,base.mu_A+base.mu.jitter),npeaks, replace = TRUE)`
# back.const.down and back.const.up: numeric; multiplicative factors (on base.mu) to set the range
#      of background noise across the two experimental states
# back.vary.down and back.vary.up: numeric; multiplicative factors (on base.mu) to set the range
#      of the second experimental state (B) to have a different background amount

simulate_reads_proportions <- function(iters, prop.do = 0.05, symmetry = TRUE, 
                                       equal.DNA.bind = TRUE,
                                       constant.background = TRUE,
                                       equal.lib.size = TRUE,
                                       replicates = 2,
                                       fraglen = 100, true.width = 500,
                                       npeaks = 10000, prior.df = 50, #maybe make less peaks maybe 500-1000 or something
                                       base.mu.jitter = 5, back.const.down = 0.33,
                                       back.const.up = 0.67, back.vary.down = 0,
                                       back.vary.up = 0.20)

{
  
  ####################################SETUP#####################################
  
  it <- iters
  dir.create(paste0("./sim_folder_",it,"_do_",prop.do,"/"))
  
  dispersion <- 1/prior.df # for the sake of possible flexibility later; right now, it just cancels out.
  disp <- prior.df*dispersion/rchisq(npeaks, df = prior.df)
  
  grouping <- c(rep("A", replicates), rep("B", replicates))
  
  radius <- fraglen
  
  ## setting sequencing depth
  sizeRatio <- truncnorm::rtruncnorm(length(grouping), a = 0.6, b = 1.4, mean = 1, sd = 0.2)
  
  ## setting the expected library size for each replicate
  
  if(equal.lib.size){
    expected.lib.size <- rep(600000, length(grouping)) 
  } else{
    expected.lib.size <- 860000 * c(rep(1.2, sum(grouping == "A")), 
                                    rep(0.7, sum(grouping == "B")))
  }
  
  ## setting how many peaks will have differential DNA occupancy
  ndo <- prop.do*npeaks
  
  if(symmetry & equal.DNA.bind){
    up.A <- 1:(ndo/2) 
    up.B <- (ndo/2 + 1):ndo
    fc.up.A <- 2 
    fc.up.B <- 2 
  } else if(symmetry & !(equal.DNA.bind)){
    up.A <- 1:(ndo/2) 
    up.B <- (ndo/2 + 1):ndo
    fc.up.A <- 4 
    fc.up.B <- 2 
  } else if(!(symmetry) & equal.DNA.bind){
    up.A <- 1:(ndo/8)
    up.B <- floor((ndo/8)+1):ndo
    fc.up.A <- 8 
    fc.up.B <- 2 
  }else{ #asymmetric and different total DNA occupancy
    up.A <- 1:(ndo/8) 
    up.B <- floor((ndo/8)+1):ndo
    fc.up.A <- 2 
    fc.up.B <- 2 
  }
  
  ## creating log2 fold change vector for convenience later on
  fc_vec <- rep(0, times = npeaks)
  fc_vec[up.A] <- log(fc.up.A, base = 2) - 0
  fc_vec[up.B] <- 0 - log(fc.up.B, base = 2)
  
  
  ####################################SIMULATING READS#######################################
  
  ## Single chromosome, for convenience.
  distances<-round(runif(npeaks, 10000, 20000))
  pos.1 <- cumsum(distances)	
  sizes <- c(chrA = max(pos.1) + 10000) # a single number that represents length of the chromosome
  chrs <- rep("chrA", npeaks) # simulating all peaks on one chromosome for simplicity
  
  fnames <- list()
  input_controls <- list()
  basesums <- c()
  
  ## given we want the baseline props to retain the correct fc,
  # we use the minimum library size to define base.mu and cur.mu 
  ## using the minimum library size allows for greater variability in cur.mu
  #and thus baseline props, holding base.mu.jitter constant
  
  base.mu <- min(expected.lib.size)/npeaks
  
  if (base.mu.jitter == 0) {
    cur.mu_input <- c(base.mu, base.mu)
  } else {
    cur.mu_input <- (base.mu-base.mu.jitter):(base.mu+base.mu.jitter)
  }
  
  cur.mu <- sample(cur.mu_input, npeaks, replace = TRUE) # samples npeaks times
  
  ## iterating through the .sam files to generate read counts in the peak regions  
  
  for (lib in 1:length(grouping)) {
    baseline_props <- cur.mu/sum(cur.mu) #resetting so we don't write over the already edited proportions 
    
    fname <- paste0("./sim_folder_", it, "_do_", prop.do,"/tf_out_", lib,"_sim_",it,"_do_",prop.do,"_cond", grouping[lib],"_exp.sam")
    input_control <- paste0("./sim_folder_", it, "_do_", prop.do,"/tf_out_", lib,"_sim_",it,"_do_",prop.do,"_cond", grouping[lib],"_input.sam")
    
    if (grouping[lib] == "A") { 
      baseline_props[up.A] <- fc.up.A*baseline_props[up.A]

      # need basesum_A to calculate the oracle size factor
      basesums[lib] <- sum(baseline_props)
      prop.mu_A <- baseline_props / sum(baseline_props)
      mu.used.A <- prop.mu_A * expected.lib.size[lib] * sizeRatio[lib]
      
      xscss::peakFile(fname, chrs=chrs, pos=pos.1, mu=mu.used.A, disp=disp,
                      sizes=sizes, fraglen=fraglen, width=true.width, tf=TRUE)
      
      xscss::peakFile(input_control, chrs=chrs, pos=pos.1, mu=0, disp=0,
                      sizes=sizes, fraglen=fraglen, width=true.width, tf=TRUE)
      
    } else { #this is for grouping B
      baseline_props[up.B] <- fc.up.B*baseline_props[up.B]

      # need basesum_B to calculate the oracle size factor
      # to allow for correct normalization even with different library sizes,
      #we multiply by the ratio of the expected library sizes 
      # when equal.library.size = TRUE, this ratio is just 1, so basesums[lib] = sum(baseline_props)
      
      basesums[lib] <- sum(baseline_props)*(expected.lib.size[grouping == "A"][1]/expected.lib.size[lib]) # lib corresponds to grouping == "B"
      prop.mu_B <- baseline_props / sum(baseline_props)
      mu.used.B <- prop.mu_B * expected.lib.size[lib] * sizeRatio[lib]
      
      xscss::peakFile(fname, chrs=chrs, pos=pos.1, mu=mu.used.B, disp=disp,
                      sizes=sizes, fraglen=fraglen, width=true.width, tf=TRUE)
      
      
      xscss::peakFile(input_control, chrs=chrs, pos=pos.1, mu=0, disp=0,
                      sizes=sizes, fraglen=fraglen, width=true.width, tf=TRUE)
      
    }
    
    fnames[[lib]] <- fname
    input_controls[[lib]] <- input_control
    
  }
  
  ####################################CREATING ORACLE SIZE FACTORS#####################################
  
  basesums <- t(basesums)
  
  
  normFactors <- sizeRatio/basesums[1,]
  SizeFactors <- as_tibble(normFactors/psych::geometric.mean(normFactors))
  
  write_csv(SizeFactors, file = paste0("./sim_folder_",it,"_do_",prop.do,"/OracleFactors_sim_",it,"_do_",prop.do,".txt"))
  
  ####################################ADDING BACKGROUND BINDING########################################
  
  fnames <- unlist(fnames)
  input_controls <- unlist(input_controls)
  
  # first, we add the background binding for the input controls
  more.background.cond <- FALSE 
  
  background.mu <- c(base.mu*back.const.down, base.mu*back.const.up)
  bdispersion <- (1/mean((background.mu)/base.mu))*(dispersion)
  
  addBackground_new(input_controls, sizes = sizes, width = 1000,
                    back.mu.min = background.mu[1],
                    back.mu.max = background.mu[2],
                    back.added.min = 0,
                    back.added.max = 0,
                    rlen = 10, 
                    dispersion.base = bdispersion,
                    dispersion.added = bdispersion,
                    more.background.cond = more.background.cond,
                    prior.df = prior.df, append = TRUE)
  
  # now we add background to experimental conditions
  if (constant.background) {
    
    more.background.cond <- FALSE 
    
    background.mu <- c(base.mu*back.const.down, base.mu*back.const.up)
    bdispersion <- (1/mean((background.mu)/base.mu))*(dispersion)
    
    addBackground_new(fnames, sizes = sizes, width = 1000,
                      back.mu.min = background.mu[1],
                      back.mu.max = background.mu[2],
                      back.added.min = 0,
                      back.added.max = 0,
                      rlen = 10, 
                      dispersion.base = bdispersion,
                      dispersion.added = bdispersion,
                      more.background.cond = more.background.cond,
                      prior.df = prior.df, append = TRUE)
    
  } else {
    
    ## one possibility is randomly choosing which sample has more background binding
    ## do so, we use the following code:  
    
    more.background.cond <- sample(c("A", "B"), 1)
    
    ## to hold constant, like we had in earlier versions of code, we can do:
    
    #more.background.cond <- "B"
    
    background.mu <- c(base.mu*back.const.down, base.mu*back.const.up)
    bdispersion.base <- (1/mean((background.mu)/base.mu))*(dispersion)
    background.added <- c(base.mu*back.vary.down, base.mu*back.vary.up)
    bdispersion.added <-(1/mean((background.mu + background.added)/base.mu))*(dispersion)
    
    addBackground_new(fnames, sizes = sizes, width = 1000,
                      back.mu.min = background.mu[1], 
                      back.mu.max = background.mu[2],
                      back.added.min = background.added[1], 
                      back.added.max = background.added[2],
                      rlen = 10, 
                      dispersion.base = bdispersion.base, 
                      dispersion.added = bdispersion.added,
                      more.background.cond = more.background.cond,
                      prior.df = prior.df, append = TRUE)
    
  }
  
  ####################################CONVERTING TO .BAM FILES#####################################
  sam.files <- c(fnames, input_controls)
  bam.files <- xscss::crunch2BAM(sam.files, dir = paste0("./sim_folder_", it,"_do_",prop.do,"/"))
  unlink(sam.files)
  
  lfile <- paste0("./sim_folder_",it,"_do_",prop.do,"/tf_log_sim_",it,"_do_",prop.do,".txt")
  
    write.table(file=lfile, data.frame(Chr= chrs[up.A],
                                       Start = pos.1[up.A]-radius,
                                       End = pos.1[up.A]+radius,
                                       log2FC = fc_vec[up.A], 
                                       more.background.cond = more.background.cond),
                row.names=FALSE, sep="\t", quote=FALSE)
    
    write.table(file=lfile, data.frame(Chr = chrs[up.B],
                                       Start = pos.1[up.B]-radius,
                                       End = pos.1[up.B]+radius,
                                       log2FC=fc_vec[up.B],
                                       more.background.cond = more.background.cond),
                row.names=FALSE, sep="\t",  quote=FALSE, append=TRUE, col.names=FALSE)
    
    lfile2 <- paste0("./sim_folder_",it,"_do_",prop.do,"/tf_log_sim_",it,"_do_",prop.do,"_true_peaks.txt")
    
    write.table(file=lfile2, data.frame(Chr= chrs,
                                       Start = pos.1-radius,
                                       End = pos.1+radius,
                                       log2FC = fc_vec, 
                                       more.background.cond = more.background.cond),
                row.names=FALSE, sep="\t", quote=FALSE)
    
  
  ####################################PEAK CALLING WITH MACS2#####################################
  
  macs2_script <- paste('sh macs2.sh ', it, prop.do)
  system(macs2_script)
  
  ####################################CLEAN UP###################################################
  
  narrowpeaks_files <- vector(length = length(grouping)) 
  summit_files <- vector(length = length(grouping)) 
  
  for (i in 1:length(grouping)){
    narrowpeaks_files[i] <- list.files(paste0("./sim_folder_",it,"_do_",prop.do,"/tf_out_",i,"_sim_",it,"_do_",prop.do,"_cond",grouping[i],"_macs_peak/"), pattern = "*.narrowPeak")
    summit_files[i] <- list.files(paste0("./sim_folder_",it,"_do_",prop.do,"/tf_out_",i,"_sim_",it,"_do_",prop.do,"_cond",grouping[i],"_macs_peak/"), pattern = "*bed")
    
    file.copy(from = paste0("./sim_folder_",it,"_do_",prop.do,"/tf_out_",i,"_sim_",it,"_do_",prop.do,"_cond",grouping[i],"_macs_peak/", narrowpeaks_files[i]),
              to = paste0("./sim_folder_",it,"_do_",prop.do, "/"))
    
    file.copy(from = paste0("./sim_folder_",it,"_do_",prop.do,"/tf_out_",i,"_sim_",it,"_do_",prop.do,"_cond",grouping[i],"_macs_peak/", summit_files[i]),
              to = paste0("./sim_folder_",it,"_do_",prop.do, "/"))
    
    # to conserve memory, we delete the macs2 peak folders after extracting the .narrowPeak files we need for downstream analysis
    unlink(paste0("./sim_folder_",it,"_do_",prop.do,"/tf_out_",i,"_sim_",it,"_do_",prop.do,"_cond",grouping[i],"_macs_peak/"), recursive = TRUE) 
    
  }
  
}

