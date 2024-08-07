# amending addBackground function that is in xscss package
# to allow for variable background binding across experimental states

addBackground_new<- function(fnames, sizes, width = 1000, 
                             back.mu.min, back.mu.max,
                             back.added.min, back.added.max, 
                             rlen = 10, prior.df = 50,  
                             more.background.cond,
                             dispersion.base = 0.05, 
                             dispersion.added = 0.05,
                             append = FALSE) 
{
  if (!append) {
    for (fname in fnames) {
      .addHeader(fname, sizes)
    }
  }
  
  width <- as.integer(width + 0.5)
  rlen <- as.integer(rlen + 0.5)
  
  for (x in names(sizes)) {
    nbins <- as.integer(sizes[[x]]/width)
    
    if (more.background.cond == "A"){
      bmu.B <-  runif(nbins, back.mu.min, back.mu.max)
      bdisp.B <- dispersion.base
      bmu.A <-  bmu.B + runif(nbins, back.added.min, back.added.max)
      bdisp.A <- dispersion.added
    } else { # covers if more.background.cond == "B" or is FALSE (as it is when constant.background = TRUE)
      bmu.A <-  runif(nbins, back.mu.min, back.mu.max)
      bdisp.A <- dispersion.base
      bmu.B <-  bmu.A + runif(nbins, back.added.min, back.added.max)
      bdisp.B <- dispersion.added
    } 
    
   #bdisp <- dispersion * prior.df/rchisq(nbins, df = prior.df)
    
    for(fname in fnames) {
      if((grepl("condA", fname) == TRUE)){
        back.y.A <- as.integer(rnbinom(nbins, mu = bmu.A, size = 1/bdisp.A) + 0.5)
        out <- .Call("R_dump_background", fname, x, width,
                     back.y.A, rlen, PACKAGE = "xscss")
        if (is.character(out)) {
          stop(out)
        }
      }
      if((grepl("condB", fname) == TRUE)){
        back.y.B <- as.integer(rnbinom(nbins, mu = bmu.B, size = 1/bdisp.B) + 0.5)
        out <- .Call("R_dump_background", fname, x, width,
                     back.y.B, rlen, PACKAGE = "xscss")
        if (is.character(out)) {
          stop(out)
        }
      }
    }
  }
  
  return(invisible()) 
}


###############################################################################################