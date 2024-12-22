# The code used to simulate ChIP-Seq reads borrows substantially from
# the code created by Lun and Smyth (2016) and Evans et al. (2016)

# require(tidyverse)
# require(xscss)
# require(csaw) 
# require(edgeR) 
# require(psych)
# require(truncnorm) 

# prop.do: numeric; proportion of peaks that have differential DNA occupancy across A and B
# symmetry: logical; indicates whether the same proportion of peaks have more DNA occupancy across A and B
# constant.background: logical; indicates whether the background binding **between** A and B for the same genomic region is constant
# equal.DNA.bind: logical; indicates whether there is the same amount of DNA Occupancy in the simulated peaks across A and B
# equal.lib.size: library size is equal across the two experimental states in the called regions
#      TRUE = 600,000, FALSE = 1,032,000 for A and 602,000 for B 
#      logic: keep above 600,000 so that the peak reads are substantially above the background reads
# npeaks: integer; number of simulated peaks in each replicate (some of which have differential DNA occupancy, some of which don't)
# prior.df: integer; connected to dispersion, used in calculating reads from means using negative binomial distribution
# base.mu.jitter: integer; when sampling the base mu values, we `sample(c(base.mu_A-base.mu.jitter,base.mu_A+base.mu.jitter),npeaks, replace = TRUE)`
# back.const.down and back.const.up: numeric; multiplicative factors (on base.mu) to set the range
#      of background noise across the two experimental states
# back.vary.down and back.vary.up: numeric; multiplicative factors (on base.mu) to set the range
#      of the second experimental state (B) to have a different background amount

simulate_means <- function(prop.do = 0.05, symmetry = TRUE,
                           equal.DNA.bind = TRUE,
                           constant.background = TRUE,
                           equal.lib.size = TRUE,
                           npeaks = 20000, prior.df = 50, 
                           base.mu.jitter = 0, back.const.down = 0.33,
                           back.const.up = 0.67, back.vary.down = 0,
                           back.vary.up = 0.20)
  
{
  
  ####################################SETUP#####################################
  
  dispersion <- 1/prior.df
  disp <- prior.df*dispersion/rchisq(npeaks, df = prior.df)
  
  grouping <- c("A", "B")
  
  ## setting sequencing depth
  sizeRatio <- truncnorm::rtruncnorm(length(grouping), a = 0.6, b = 1.4, mean = 1, sd = 0.2)
  
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
  
  ####################################SIMULATING READS#######################################
  
  ## Single chromosome, for convenience.
  distances <- round(runif(npeaks, 10000, 20000))
  pos.1 <- cumsum(distances)	
  sizes <- c(chrA = max(pos.1) + 10000) # a single number that represents length of the chromosome
  chrs <- rep("chrA", npeaks) # simulating all peaks on one chromosome for simplicity
  
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
  
  cur.mu <- sample(cur.mu_input, npeaks, replace = TRUE)
  
  ####################################################################################################
  mu_A_df <- data.frame(peak_index = c(), position = c(), mu_A = c())
  mu_B_df <- data.frame(peak_index = c(), position = c(), mu_B = c())
  
  
  for (lib in 1:length(grouping)) {
    baseline_props <- cur.mu/sum(cur.mu)
    
    if (grouping[lib] == "A") { 
      baseline_props[up.A] <- fc.up.A*baseline_props[up.A]
      
      prop.mu_A <- baseline_props / sum(baseline_props)
      mu.used.A <- prop.mu_A * expected.lib.size[lib] * sizeRatio[lib]
      
      mu_A_df <- mu_A_df |> rbind(data.frame(peak_index = 1:npeaks, position = pos.1, mu_A = mu.used.A))
      
    } else { #this is for grouping B
      baseline_props[up.B] <- fc.up.B*baseline_props[up.B]
      
      prop.mu_B <- baseline_props / sum(baseline_props)
      mu.used.B <- prop.mu_B * expected.lib.size[lib] * sizeRatio[lib]
      
      mu_B_df <- mu_B_df |> rbind(data.frame(peak_index = 1:npeaks, position = pos.1, mu_B = mu.used.B))
      
    }
    
  }
  
  simulation_mu_df <- inner_join(mu_A_df, mu_B_df, by = c('peak_index', 'position'))
  
  width <- 1000
  nbins <- as.integer(sizes[['chrA']]/width)
  ## adding background based on parameter value...
  if (constant.background) {
    
    more.background.cond <- FALSE 
    
    background.mu <- c(base.mu*back.const.down, base.mu*back.const.up)
    background.added <- c(0,0)
    
  } else {

    more.background.cond <- "B"
    
    background.mu <- c(base.mu*back.const.down, base.mu*back.const.up)
    background.added <- c(base.mu*back.vary.down, base.mu*back.vary.up)
    
  }
  
  bmu.A <- runif(nbins, background.mu[1], background.mu[2])
  bmu.B <- bmu.A + runif(nbins, background.added[1], background.added[2])

  pos.back.bin <- seq(width/2, sizes - width/2, by = width)
  
  simulation_back_df <- data.frame(back_A = bmu.A, back_B = bmu.B, back_pos = pos.back.bin)
  
  return(list(back = simulation_back_df, mu = simulation_mu_df))
  
  }

