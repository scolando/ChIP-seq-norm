# require(MAnorm2)
# require(DiffBind)


MAnorm2.sim <- function(df, iters, df_all = NULL, prop.do = 0.05,
                        equal.DNA.bind = TRUE,
                        symmetry = TRUE, constant.background = TRUE, 
                        equal.lib.size = TRUE,
                        consensus_peaks, replicates, prop.consensus,
                        FDR.threshold = 0.05) {
  
  it <- iters
  
  technique <- "MAnorm2"
  
## reading in table with regions that have (true) differential DNA occupancy
  true_do <- readr::read_table(paste0("./sim_folder_",it,"_do_",prop.do,"/tf_log_sim_",it, "_do_",prop.do,".txt")) 
  more.background.cond <- true_do$more.background.cond[1] # adding in the condition with more background. FALSE indicates the conditions have the same amount
  
## making sure the norm is only wrt to consensus peaks used for rest of methods for a more direct comparison
  
  consensus_starts <- consensus_peaks@ranges@start
  
  df <- df |>
    dplyr::filter(start %in% consensus_starts) |>
    dplyr::rename(Chr = chrom,
                  Start = start,
                  End = end)
  
  called_true_do_peaks <- genome_inner_join(df, true_do, by = c("Chr", "Start", "End"))
  
  n_true_do <- nrow(true_do)
  n_consensus_peaks <- nrow(df) # this is number of consensus peaks that live in MAnorm2 file
  n_true_do_called <- nrow(called_true_do_peaks) # intersection of true_do and consensus_peaks from MAnorm2 

## MAnorm2's within-group normalization
  
  norm <- MAnorm2::normalize(df, count = grep("s1.read", names(df)), 
                             occupancy = grep("s1.occupancy", names(df)))
  norm <- MAnorm2::normalize(norm, count = grep("s2.read", names(df)), 
                             occupancy = grep("s2.occupancy", names(df)))
  
## constructing a bioCond object for both A and B
  conds <- list(cond_A = bioCond(norm[grep("s1.read", names(df))[1]:tail(grep("s1.read", names(df)), n = 1)], 
                                 norm[grep("s1.occupancy", names(df))[1]:tail(grep("s1.occupancy", names(df)), n = 1)],
                                 name = "cond_A"),
                cond_B = bioCond(norm[grep("s2.read", names(df))[1]:tail(grep("s2.read", names(df)), n = 1)],
                                 norm[grep("s2.occupancy", names(df))[1]:tail(grep("s2.occupancy", names(df)), n = 1)],
                                 name = "cond_B"))
  
## MAnorm2's between-group normalization
  # common.peak.regions is left unspecified (i.e., NULL) so that they are defined using MAnorm2's default
  conds <- normBioCond(conds, common.peak.regions = NULL) 
  
  # fittng mean variance curve, set init.coef to c(0.1, 10) based on suggestion in vignette for fitMeanVarCurve function
  conds <- fitMeanVarCurve(conds, method = "parametric", occupy.only = TRUE, init.coef = c(0.1, 10))
  
## creating data frame of regions and their associated adjusted p-values (which are calculated using diffTest function)
  called_regions <- MAnorm2::diffTest(conds[[1]], conds[[2]]) |>
    cbind(df[1:3]) |>
    select(Chr, Start, End, everything()) 
  
## calculatng true pos, false pos, etc. for MAnorm2
  
  true_positives <- genome_inner_join(true_do, called_regions, by = c("Chr", "Start", "End")) |>
    dplyr::filter(padj < FDR.threshold) |> nrow()
  
  false_positives<- genome_anti_join(called_regions, true_do, by = c("Chr", "Start", "End")) |>
    dplyr::filter(padj < FDR.threshold) |> nrow()
  
  false_negatives <- genome_left_join(true_do, called_regions, by = c("Chr", "Start", "End")) |>
    mutate(padj = replace_na(padj, 1)) |> dplyr::filter(padj >= FDR.threshold) |> nrow()
  
  true_negatives <- genome_anti_join(called_regions, true_do, by = c("Chr", "Start", "End")) |>
    mutate(padj = replace_na(padj, 1)) |> dplyr::filter(padj >= FDR.threshold) |> nrow()
  
  FDR <- (false_positives)/((false_positives)+(true_positives))
  
## creating a df with false pos., true pos. etc for MAnorm2
  
  df <- data.frame(technique, true_positives, false_positives, true_negatives, 
                   false_negatives, FDR, n_true_do_called, 
                   n_consensus_peaks, n_true_do, it, prop.do,
                   symmetry, constant.background,
                   equal.DNA.bind, equal.lib.size,
                   prop.consensus, more.background.cond, FDR.threshold)
  
## binding to df_all which has rest of normalization methods' results  
  
  rbind(df_all, df)
  
}