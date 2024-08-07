# require(MAnorm2)
# require(DiffBind)
# require(truncnorm)
# require(fuzzyjoin)
# require(tidyverse)

source("MAnorm2.sim.R")
source("dba.normalize.sim.R")

read_analysis <- function(iters, prop.do = 0.05, symmetry = TRUE,
                          equal.DNA.bind = TRUE,
                          constant.background = TRUE,
                          equal.lib.size = TRUE,
                          replicates = 2, prop.consensus = 0.75,
                          FDR.threshold = 0.05){
  
  it <- iters
  
  
  ####################################DIFFBIND: CREATING INITIAL DBA OBJECT##################################### 
  
  ## First, we create a diffbind object with the bams, reads, and peak files
  
  bamReads <- as_tibble(paste(paste0("./sim_folder_",it,"_do_",prop.do,"/"), list.files(paste0("./sim_folder_",it,"_do_",prop.do,"/"), pattern= "*.bam$")[!str_detect(list.files(paste0("./sim_folder_",it,"_do_",prop.do,"/"), pattern="*.bam$"), "bam.bai")], sep="")) |>
    mutate(Condition = str_extract(value, "cond(.)", group = TRUE),
           number = as.numeric(str_extract(value, "out_(.*)_sim", group = TRUE))) |>
    arrange(Condition, number) |>
    mutate(Replicate= rep(1:replicates, 2),
           SampleID = paste0("cond",Condition,"_rep",Replicate))
  
  peaks <- as_tibble(paste(paste0("./sim_folder_",it,"_do_",prop.do,"/"), list.files(paste0("./sim_folder_",it,"_do_",prop.do,"/"), pattern="*.narrowPeak"), sep="")) |>
    mutate(Condition = str_extract(value, "cond(.)", group = TRUE),
           sample_id = as.numeric(str_extract(value, "out_(.*)_sim", group = TRUE))) |>
    arrange(Condition, sample_id) |>
    mutate(Replicate= rep(1:replicates, 2))
  
  
  sample_sheet<-data.frame(
    SampleID= bamReads$SampleID,
    Condition= bamReads$Condition,
    Replicate= bamReads$Replicate,
    bamReads= bamReads$value,
    Peaks= peaks$value, 
    PeakCaller= rep("narrow", 2*replicates)
  )
  
  res_s <- dba(sampleSheet=sample_sheet,config=data.frame(AnalysisMethod=DBA_DESEQ2, th=FDR.threshold,
                                                          DataType=DBA_DATA_GRANGES, RunParallel=FALSE, 
                                                          bUsePval=FALSE, design=TRUE,
                                                          doBlacklist=FALSE, doGreylist=FALSE))
  
  ####################################DIFFBIND: CREATING CONSENSUS PEAK SET##################################### 
  
  ## adds consensus peaks: those that overlap with at least (prop.consensus)*(minimum number of reps)
  
  grouping <- c(rep("A", replicates), rep("B", replicates))
  
  #first we calculate our minOverlap:
  minReps <- ceiling(prop.consensus*min(table(grouping)))
  
  res_consensus <- dba.peakset(res_s, consensus = DBA_CONDITION, 
                               minOverlap = minReps)
  
  ## removes other peak sets, focuses only on the consensus peaks
  res_consensus <- dba(res_consensus, mask=res_consensus$masks$Consensus,
                       minOverlap=1)
  
  ## creates an object containing the consensus peaks, and writes it out as a bed file
  consensus_peaks <- dba.peakset(res_consensus, bRetrieve=TRUE)
  
  ## finding the proportion of peaks with differential DNA occupancy peaks that are are actually called
  ## we start by reading in the .txt file created during data generation (i.e., in simulate_reads_proportions.R) that identifies the regions with differential DNA occupancy 
  
  true_do <- readr::read_table(paste0("./sim_folder_",it,"_do_",prop.do,"/tf_log_sim_",it, "_do_",prop.do,".txt"))
  consensus_peaks_df <- data.frame(Start = consensus_peaks@ranges@start, Width = consensus_peaks@ranges@width) |>
    mutate(End = Start + Width) |> mutate(Chr = "chrA")
  
  called_true_do_peaks <- genome_inner_join(consensus_peaks_df, true_do, by = c("Chr", "Start", "End"))
  
  n_true_do <- nrow(true_do)
  n_consensus_peaks <- nrow(consensus_peaks_df)
  n_true_do_called <- nrow(called_true_do_peaks) # intersection of true_do and consensus_peaks 
  
  ####################################DIFFBIND: COUNTING READS##################################### 
  
  ## Finding summit number
  
  summary <- summary(res_s$binding[,3]-res_s$binding[,2])
  min <- unname(0.5*(summary[1]))
  quart_1 <- unname(0.5*(summary[2]))
  summits <- ceiling((min+quart_1)/2)
  
  ## Now we are creating the dba counts/contrasts to be used in DiffBind normalization
  # Jo added bParallel=FALSE (and below) to help with running on the HPC
  # from this help file: https://www.biostars.org/p/9512204/#9548140
  
  res_c <- dba.count(res_s, peaks = consensus_peaks, summits = summits, bParallel=FALSE)
  
  res_c <- dba.contrast(res_c, contrast=c("Condition", 'A', 'B'))
  
  # Jo from this help file: https://support.bioconductor.org/p/p133921/
  # set the number of cores being used to normalize
  res_c$config$cores <- 1
  
  ## creating a directory for created csvs, which are useful when running code manually
  
  dir.create(paste0("./sim_folder_",it,"_do_",prop.do,"/norm_csvs_",it,"_do_",prop.do,"/"))
  
  more.background.cond <- true_do$more.background.cond[1] #retrieving condition with more background. If background is constant, this is FALSE
  
  ####################################ORACLE DBA ANALYSIS#####################################
  
  technique <- "Oracle"
  
  ## we start by reading in the .txt file that was created durng data generation (in simulate_reads_proportions.R) 
  
  OracleFactors <- read_csv(paste0("./sim_folder_",it,"_do_",prop.do, "/OracleFactors_sim_",it,"_do_",prop.do,".txt"))
  OracleFactors <- OracleFactors$value
  
  res_Oracle <- dba.normalize(res_c, normalize = OracleFactors)
  res_Oracle <- dba.analyze(res_Oracle)
  dba.peakset(res_Oracle, bRetrieve=TRUE, 
              writeFile=paste0("./sim_folder_",it,"_do_",prop.do,"/norm_csvs_",it,"_do_",prop.do,"/","Oracle_",it,"_do_",prop.do,".csv"))
  
  norm_factor_oracle <- c()
  sample_id_oracle <- c()
  
  for (rep in 1:(2*replicates)){
    norm_factor_oracle[rep] <- OracleFactors[rep]
    sample_id_oracle[rep] <- res_Oracle[["samples"]][["SampleID"]][rep]
  }
  
  ## calculating FP, FN, TP, TNs, and Power for the Oracle
  
  info_df <- res_Oracle$peaks[[1]] |> mutate(id = 1:nrow(res_Oracle$peaks[[1]]))
  called_regions <- inner_join(info_df, res_Oracle$contrasts[[1]]$DESeq2$de, by = "id")
  
  true_positives <- genome_inner_join(true_do, called_regions, 
                                      by = c("Chr", "Start", "End")) |> filter(padj < FDR.threshold) |> nrow()
  
  false_positives <- genome_anti_join(called_regions, true_do, 
                                      by = c("Chr", "Start", "End")) |> filter(padj < FDR.threshold) |> nrow()
  
  false_negatives <- genome_left_join(true_do, called_regions, by = c("Chr", "Start", "End")) |>
    mutate(padj = replace_na(padj, 1)) |> filter(padj >= FDR.threshold) |> nrow()
  
  true_negatives <- genome_anti_join(called_regions, true_do, by = c("Chr", "Start", "End")) |>
    filter(padj >= FDR.threshold) |> nrow()
  
  FDR <- false_positives/(false_positives+true_positives)
  
  ## creating data frames
  
  df_all <- data.frame(technique, true_positives, false_positives, true_negatives, 
                       false_negatives, FDR, n_true_do_called, 
                       n_consensus_peaks, n_true_do, it, prop.do,
                       symmetry, constant.background,
                       equal.DNA.bind, equal.lib.size,
                       prop.consensus, more.background.cond, FDR.threshold)
  
  
  df_norm_factors <- data.frame(technique, t(setNames(norm_factor_oracle, sample_id_oracle)))
  
  ####################################DIFFBIND DBA ANALYSIS#####################################
  
  norm_techniques <- c("LIB_PEAKS","TMM_PEAKS","RLE_PEAKS","LIB_BINS","TMM_BINS","RLE_BINS", "DBA_NORM_OFFSETS_ADJUST")
  norm_factor_tech <- c()
  sample_id_tech <- c()
  
  for(tech in norm_techniques) {
    
    technique <- paste0(tech)
    
    res_norm <- dba.normalize.sim(res_c, norm_method = tech)
    
      norm_factors <- as.numeric(dba.normalize(res_norm, bRetrieve = TRUE)$norm.factors)
      for (rep in 1:(2*replicates)){
        norm_factor_tech[rep] <- norm_factors[rep]
        sample_id_tech[rep] <- res_norm[["samples"]][["SampleID"]][rep]
      }
    
    print(paste0("Performing DiffBind analysis for ",tech))
    res_norm <- dba.analyze(res_norm)
    dba.peakset(res_norm, bRetrieve=TRUE, 
                writeFile=paste0("./sim_folder_",it,"_do_",prop.do,"/norm_csvs_",it,"_do_",prop.do,"/",tech,"_",it,"_do_",prop.do,".csv"))
    
    ## calculating FP, FN, TP, TNs, and Power for the normalization technique
    
    info_df <- res_norm$peaks[[1]] |> mutate(id = 1:nrow(res_norm$peaks[[1]]))
    called_regions <- inner_join(info_df, res_norm$contrasts[[1]]$DESeq2$de, by = "id") 
    
    true_positives <- genome_inner_join(true_do, called_regions, by = c("Chr", "Start", "End")) |>
      filter(padj < FDR.threshold) |> nrow()
    
    false_positives<- genome_anti_join(called_regions, true_do, by = c("Chr", "Start", "End")) |>
      filter(padj < FDR.threshold) |> nrow()
    
    false_negatives <- genome_left_join(true_do, called_regions, by = c("Chr", "Start", "End")) |>
      mutate(padj = replace_na(padj, 1)) |> filter(padj >= FDR.threshold) |> nrow()
    
    true_negatives <- genome_anti_join(called_regions, true_do, by = c("Chr", "Start", "End")) |>
      filter(padj >= FDR.threshold) |> nrow()
    
    FDR <- (false_positives)/((false_positives)+(true_positives))
    
    ## creating data frames   
    
    df <- data.frame(technique, true_positives, false_positives, true_negatives, 
                     false_negatives, FDR, n_true_do_called, 
                     n_consensus_peaks, n_true_do, it, prop.do,
                     symmetry, constant.background,
                     equal.DNA.bind, equal.lib.size,
                     prop.consensus, more.background.cond, FDR.threshold)
    
    df_all <- rbind(df_all, df)
    
    df_norm <- data.frame(technique, t(setNames(norm_factor_tech, sample_id_tech)))
    
    
    df_norm_factors <- rbind(df_norm_factors, df_norm)
    
  }
  
  ####################################MANORM2 DBA ANALYSIS#####################################
  
  ## MAnorm2 is a non-DiffBind supported norm technique so we create a new folder and move relevant files into it and rename them for convenience 
  
  dir.create(paste0("./sim_folder_",it,"_do_",prop.do,"/MAnorm2/"))
  
  bam_files <- list.files(path = paste0("./sim_folder_",it,"_do_",prop.do,"/"),
                          pattern = ".bam$")
  
  for(i in 1:length(bam_files)){
    file.copy(from = paste0("./sim_folder_",it,"_do_",prop.do,"/",bam_files[i]),
              to = paste0("./sim_folder_",it,"_do_",prop.do,"/MAnorm2/",bam_files[i]))
    new_filename <- gsub(paste0("_sim_",it,"_do_",prop.do,"_cond."), "", paste0("./sim_folder_",it,"_do_",prop.do,"/MAnorm2/",bam_files[i]))
    file.rename(from = paste0("./sim_folder_",it,"_do_",prop.do,"/MAnorm2/",bam_files[i]), 
                to = new_filename)
  }
  
  narrowPeak_files <- list.files(path = paste0("./sim_folder_",it,"_do_",prop.do,"/"),
                                 pattern = ".narrowPeak$")
  
  for(i in 1:length(narrowPeak_files)){
    file.copy(from = paste0("./sim_folder_",it,"_do_",prop.do,"/",narrowPeak_files[i]),
              to = paste0("./sim_folder_",it,"_do_",prop.do,"/MAnorm2/",narrowPeak_files[i]))
    new_filename <- gsub(paste0("_sim_",it,"_do_",prop.do,"_cond."), "", paste0("./sim_folder_",it,"_do_",prop.do,"/MAnorm2/",narrowPeak_files[i]))
    new_filename <- gsub(".bam_macs_peak_peaks", "", new_filename)
    file.rename(from = paste0("./sim_folder_",it,"_do_",prop.do,"/MAnorm2/",narrowPeak_files[i]), 
                to = new_filename)
  }
  
  #### now that our files are conveniently placed and named we will run our shell script which creates the .xls file MAnorm2 uses as its input
  
  MAnorm2_script <- paste('MAnorm2.sh ', it, prop.do, 2*replicates)
  system(MAnorm2_script)
  MAnorm2 <- read.table(paste0("./sim_folder_",it,"_do_",prop.do,"/MAnorm2/","cond_A_and_B_profile_bins.xls"), sep = "\t", header = TRUE)
  
  ### MAnorm2.sim performs the MAnorm2 analysis steps that are outlined in the MAnorm2 package vignette
  
  df_all_final <- MAnorm2.sim(MAnorm2, iters = it, df_all = df_all, prop.do = prop.do,
                              equal.DNA.bind = equal.DNA.bind,
                              symmetry = symmetry, constant.background = constant.background,
                              equal.lib.size = equal.lib.size,
                              consensus_peaks = consensus_peaks, prop.consensus = prop.consensus,
                              FDR.threshold = FDR.threshold)
  
  ####################################RETURN STATEMENT#####################################
  
  return(full_join(df_all_final, df_norm_factors, by = "technique"))
  
}