# require(DiffBind)
# this is adapted from DiffBind's dba.normalize() in order to allow the norm_method to be input as a string
# captures the 7 normalization methods outlined on pg. 46 of DiffBind's vignette (Section 7.5)

dba.normalize.sim <- function(DBA, method = DBA$config$AnalysisMethod,
                              norm_method = NULL, library = NULL, 
                              background = FALSE, spikein = FALSE, offsets = FALSE,
                              libFun=mean, bRetrieve=FALSE, ...){
  
  # Jo from this help file: https://support.bioconductor.org/p/p133921/
  # set the number of cores being used to normalize
  
  ## 3 below normalize via reads in peaks
  
  if (norm_method == "LIB_PEAKS") {
    DBA$config$cores <- 1
    res_norm <- dba.normalize(DBA, normalize = DBA_NORM_LIB, library = DBA_LIBSIZE_PEAKREADS,
                              background=FALSE, spikein=FALSE)
  }
  
  else if (norm_method == "TMM_PEAKS") {
    DBA$config$cores <- 1
    res_norm <- dba.normalize(DBA, normalize = DBA_NORM_TMM, library = DBA_LIBSIZE_PEAKREADS,
                              background=FALSE, spikein=FALSE)
  }
  
  else if (norm_method == "RLE_PEAKS") {
    DBA$config$cores <- 1
    res_norm <- dba.normalize(DBA, normalize = DBA_NORM_RLE, library = DBA_LIBSIZE_PEAKREADS,
                              background=FALSE, spikein=FALSE)
  }
  
  ## 3 methods below normalize via background bins
  
  else if (norm_method == "LIB_BINS") {
    DBA$config$cores <- 1
    res_norm <- dba.normalize(DBA, normalize = DBA_NORM_LIB, library = DBA_LIBSIZE_BACKGROUND,
                              background=TRUE, spikein=FALSE)
  }
  
  else if (norm_method == "TMM_BINS") {
    DBA$config$cores <- 1
    res_norm <- dba.normalize(DBA, normalize = DBA_NORM_TMM, library = DBA_LIBSIZE_BACKGROUND,
                              background=TRUE, spikein=FALSE)
  }
  
  else if (norm_method == "RLE_BINS") {
    DBA$config$cores <- 1
    res_norm <- dba.normalize(DBA, normalize = DBA_NORM_RLE, library = DBA_LIBSIZE_BACKGROUND,
                              background=TRUE, spikein=FALSE)
  }
  
  # offset normalization (using loess on reads in peaks). Offsets will be adjusted for library size and mean centering pre-dba  
  
  else if (norm_method == "DBA_NORM_OFFSETS_ADJUST") {
    DBA$config$cores <- 1
    res_norm <- dba.normalize(DBA, normalize = DBA_NORM_OFFSETS_ADJUST, library = DBA_LIBSIZE_PEAKREADS,
                              background=FALSE, spikein=FALSE, offsets = TRUE)
  }
  
  ## if normalization method string is not recognized/not provided, then will use DiffBind's default normalization method.
  
  else{
    
    warning("normalization method string wasn't recognized, using DiffBind's default...")
    DBA$config$cores <- 1
    res_norm <- dba.normalize(DBA, normalize = DBA_NORM_DEFAULT) #uses default (more info in dba.normalize's vingette)
    
  }
  
  
  return(res_norm)
  #return()
  print("finished!")
  
}