#' @title DSBInt
#' @description simulates DNA DSBs across virtual chromosomes. Requires input hotspot map
#' @param chrom specifies which chromosome to simulate
#' @param ActiveChromatids the number of chromatids to simulate. This has no effect on the global number of DSBs, which is always scaled to match a situation where all 4 chromatids exist. If trans interference is on, trans interference will be applied across all chromatids.
#' @param RepN number of DSBs to simulate. Represents the number of global DSBs, across all chromosomes and chromatids (as if there are 4 chromatids per chromosome)
#' @param WTmap the DSB map with interference (i.e wild-type). Used to calculate the real ratio
#' @param tel1Dmap the DSB map without interference (i.e tel1D). Used for DSB site selection and calculation of both real and sim ratios
#' @param script_name name of the script. Affixed to the end of the output folder name
#' @param output_dir name of the directory to output simulated data. Will be added to subdirectory entitled "sim_parts".
#' @param expname if given, puts all of the sim parts in a subfolder with this name.
#' @param Winmethod the window function to use for simulation. Defaults to "Hann", also supports "Exponential" and "Tukey".
#' @param winr the r value used for the Tukey window. Only has an effect if Winmethod is set to "Tukey".
#' @param Windows the window widths for the interference function (in kb). Wider windows will spread further from the DSB. Supports lists, simulator will loop for each value given.
#' @param trans when set to "ON", DSBs will cause interference across all simulated chromatids.
#' @param cisStrength the strength of cis interference.
#' @param transStrength the strength of trans interference. Supports lists, simulator will loop for each value given.
#' @param fail_on if set to "TRUE", DSBs will fail. Failed DSBs will generate interference but not contribute to the SimRatio or the DSB number. The rate of failure can be set with failrate.
#' @param failrate the chance of a DSB to fail to mature. Requires fail_on to be TRUE to have any effect.
#' @param failStrength The strength of DSB failure interference, as a proportion of successful DSB interference strength. Only functions if fail_on is set to "TRUE" and failrate is greater than 0.
#' @param res the width of bins in bp.
#' @param G the depth of the sim as total number of DSBs to simulate (as a function of global). This ensures that all sims are the same depth, regardless of the number of DSBs simulated.
#' @param min_prob_threshold the minimum possible probability of a DSB forming in any given bin. Ensures the simulator will always run to completion, no matter how strong interference is.
#' @param Smooth the strength of the smoother used to produce intermediate plots. Has no effect on the simulator output.
#' @param ChrSizes the sizes of the chromosomes. Defaults to cerevisiae chromosomes. The chromosome to simulate can be set with chr.
#' @param CEN the positions of centromeres. Only affects the intermediate plots, does not affect simulator output.
#' @param ncores the number of cores to run the simulator across in parallel.

DSBInt <- function(
    chrom = 15,
    ActiveChromatids = 2,
    RepN = 200,
    tel1Dmap,
    WTmap,
    script_name = "DSBInt_v23T",
    output_dir = "./Output",
    expname = NULL,
    Winmethod = "Hann",
    winr = 0.5,
    Windows = c(10, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900, 1000),
    trans = "ON",
    cisStrength = 1,
    transStrength = seq(0, 1, by = 0.1),
    fail_on = FALSE,
    failrate = 0.5,
    failStrength = 1,
    res = 100,
    G = 1e7,
    min_prob_threshold = 1e-10,
    Smooth = 0.1,
    ChrSizes = c(230218,813184,(316620+1173+3077),1531933,576874,270161,
                 1090940,562643,439888,745751,666816,1078177,
                 924431,784333,1091291,948066),
    CEN = c(151523.5, 238265, 114443, 449766,152045.5, 148568.5,
            496979,105644.5, 355687, 436366,440187.5, 150887.5,
            268090, 628816.5, 326643, 556015),
    ncores = max(1, parallel::detectCores() - 1)
){
  # Try to set working directory based on RStudio, otherwise fall back
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  } else {
    args_full <- commandArgs(trailingOnly = FALSE)
    script_path <- sub("--file=", "", args_full[grep("--file=", args_full)])
    if (length(script_path) > 0) {
      setwd(dirname(normalizePath(script_path)))
    }
  }
  
  #list of packages
  packages = c('e1071', 'tictoc', 'data.table', 'plyr', 'dplyr', 'bio3d', 'ggplot2', 'viridis',
               'gridExtra', 'doParallel', 'foreach', 'doFuture', 'ggnewscale', 'bspec')
  
  #load/install
  package.check <- lapply(
    packages,
    FUN = function(x) {
      if (!require(x, character.only = TRUE)) {
        install.packages(x, dependencies = TRUE)
        library(x, character.only = TRUE)
      }
    }
  )
  
  
  #libraries
  library(e1071)
  library(tictoc)
  library(data.table)
  library(plyr)
  library(dplyr)
  library(bio3d)
  library(ggplot2)
  library(viridis)
  library(gridExtra)
  library(doParallel)
  library(foreach)
  library(doFuture)
  library(ggnewscale)
  library(bspec)
  
  packagelist<- c('data.table', 'plyr', 'ggplot2', 'e1071', 'tictoc', 'bio3d', 
                  'viridis', 'gridExtra', 'doParallel', 'foreach', 'doFuture',
                  'dplyr', 'bspec')
  # Helper function to handle file path OR dataframe
  load_hotspot <- function(x) {
    if (is.data.frame(x)) {
      # already a data.frame in memory
      return(x)
    } else if (is.character(x) && length(x) == 1 && file.exists(x)) {
      # file path given, load it
      return(read.table(x, header = TRUE, stringsAsFactors = FALSE))
    } else {
      stop("Invalid hotspot input: must be a data.frame or an existing file path")
    }
  }
  # Load A/B hotspot data (file path or dataframe)
  A.hotspot <- load_hotspot(tel1Dmap)
  B.hotspot <- load_hotspot(WTmap)
  
  #reformatting for newer HS tables
  if ("NormHpMChr" %in% names(A.hotspot) && !"NormHpChr" %in% names(A.hotspot)) {
    A.hotspot$NormHpChr <- A.hotspot$NormHpMChr
    #message("Renamed 'NormHpMChr' -> 'NormHpChr' in A.hotspot")
  }
  A.hotspot <- A.hotspot[, c("Chr", "Midpoint", "NormHpChr"), drop = FALSE]
  
  tic() # total timer
  
  #calculate parameters
  E1= max(RepN)
  RepTN<- round(G/RepN)
  maxsim<- max(RepTN)
  Gratio = sum(subset(A.hotspot, A.hotspot$Chr == chrom)$NormHpM)/sum(A.hotspot$NormHpM)
  
  #Exp window function
  exponential_window <- function(N, target_width) {
    n <- 0:(N-1)
    center <- (N-1)/2
    mu <- log(2) / (target_width / 2)
    w <- exp(-mu * abs(n - center))
    return(w)
  }
  
  #exp windows precompute
  if(Winmethod == 'Exponential'){
    exp_windows <- list()
    for (w in Windows) {
      maxlength = max(Windows)*1000/res
      exp_windows[[as.character(w)]] <- exponential_window(w*1000/res, w*1000/res)
    }
  }
  
  fail_label <- sprintf("%.2f", failrate)
  
  # directories
  tpart = '_varying_trans'
  t <- Sys.time()
  t <- strsplit(as.character(Sys.time()), ' ')
  t <- paste(t[[1]][1], gsub("\\.", "", t[[1]][2]), sep = '.')
  t <- gsub('.{9}$', '', t)
  t <- gsub(':', '\\.', t)
  
  # format failrate
  if (exists("fail_on") && fail_on) {
    fail_label <- sprintf("failrate_%.2f_", failrate)
  } else {
    fail_label <- ""
  }
  
  # build folder name
  # directories
  tpart = '_varying_trans'
  t <- Sys.time()
  t <- strsplit(as.character(Sys.time()), ' ')
  t <- paste(t[[1]][1], gsub("\\.", "", t[[1]][2]), sep = '.')
  t <- gsub('.{9}$', '', t)
  t <- gsub(':', '\\.', t)
  
  # format failrate nicely
  if (exists("fail_on") && fail_on) {
    fail_label <- sprintf("failrate_%.2f_", failrate)
  } else {
    fail_label <- ""
  }
  
  # build folder name (with failrate inserted before script)
  # directories
  tpart = '_varying_trans'
  t <- Sys.time()
  t <- strsplit(as.character(Sys.time()), ' ')
  t <- paste(t[[1]][1], gsub("\\.", "", t[[1]][2]), sep = '.')
  t <- gsub('.{9}$', '', t)
  t <- gsub(':', '\\.', t)
  
  # format failrate nicely
  if (exists("fail_on") && fail_on) {
    fail_label <- sprintf("failrate_%.2f_", failrate)
  } else {
    fail_label <- ""
  }
  
  # build folder name (with failrate inserted before script)
  fname <- paste(
    ActiveChromatids, 'chroms_',
    as.numeric(G) / 1000000, 'M_',
    min(RepN), '-', max(RepN), 'DSBs_',
    Winmethod, min(Windows), '-', max(Windows),
    '_Chrom', chrom, tpart, '_',
    fail_label,
    t, '_',
    script_name,
    sep = ''
  )
  # Base folder (for simulation parts)
  sim_parts_base <- file.path(output_dir, "sim_parts")
  
  # If expname is specified, nest one level deeper
  if (!is.null(expname) && nzchar(expname)) {
    sim_parts_base <- file.path(sim_parts_base, expname)
  }
  
  # Create full simulation subfolder
  sdir <- file.path(sim_parts_base, fname)
  
  # Make directories if needed
  if (!dir.exists(sdir)) dir.create(sdir, recursive = TRUE)
  
  # Similarly adjust other directories
  pdir  <- file.path(output_dir, "plots", fname)
  psdir <- file.path(output_dir, "sim_plots", fname)
  if (!file.exists(pdir)) dir.create(file.path(pdir))
  if (!file.exists(sdir)) dir.create(file.path(sdir))
  if (!file.exists(psdir)) dir.create(file.path(psdir))
  
  #Subset hotspots
  needed_cols <- c("Chr", "Midpoint", "NormHpChr")
  A1 <- subset(A.hotspot, Chr == chrom)[, intersect(needed_cols, names(A.hotspot)), drop = FALSE]
  B1 <- subset(B.hotspot, Chr == chrom)[, intersect(needed_cols, names(B.hotspot)), drop = FALSE]
  A1$Midpoint=round(A1$Midpoint/res)
  B1$Midpoint=round(B1$Midpoint/res)
  
  A1_dt <- as.data.table(A1)
  B1_dt <- as.data.table(B1)
  #Keep only numeric columns before summing
  numcols_A1 <- names(A1_dt)[sapply(A1_dt, is.numeric)]
  A1_agg <- A1_dt[, lapply(.SD, sum, na.rm=TRUE), by = Midpoint, .SDcols = numcols_A1]
  
  numcols_B1 <- names(B1_dt)[sapply(B1_dt, is.numeric)]
  B1_agg <- B1_dt[, lapply(.SD, sum, na.rm=TRUE), by = Midpoint, .SDcols = numcols_B1]
  A1 <- A1_agg
  B1 <- B1_agg
  
  L=round(ChrSizes[chrom]/res)
  A=rep(0,L)
  A[A1$Midpoint]=A1$NormHpChr
  B=rep(0,L)
  B[B1$Midpoint]=B1$NormHpChr
  
  ### FIX: ensure integer, at least 1
  E_raw <- L/(12500000/res)*E1/4*ActiveChromatids
  E <- as.integer(round(E_raw))
  if (is.na(E) || E < 1) E <- 1
  
  Repseq = floor((E/E1)*RepN)
  G1=seq(1000,1000000, by=1000)
  AllResult<- NULL
  
  plan(multisession, workers = ncores)
  foreach(k = Windows, .options.future = list(packages = packagelist, seed = TRUE), .combine = rbind) %dofuture%{
    tic()
    for(tr in transStrength){
      W1=k
      W = W1*1000/res
      W <- as.integer(round(W))
      W <- max(1, min(W, 2*L))   ### FIX cap window
      
      if(Winmethod == 'Hann'){
        C_full = hanning.window(W)
        C = 1 - C_full
        winmethod = Winmethod
      }else if(Winmethod == 'Exponential'){
        C <- 1 - exp_windows[[as.character(k)]]
        winmethod = Winmethod
      }else{
        C = 1 - tukeywindow(W, winr)
        winmethod = paste(Winmethod, winr, sep = '_')
      }
      
      H=vector("list", E)
      for (idx in seq_len(E)) H[[idx]] <- integer(L)
      
      F=integer(E)
      # updated cellcount fields to store success fraction and counts
      cellcount <- data.table(DSBs = integer(0), GDSBs = integer(0), Cells = integer(0),
                              success_frac = numeric(0), success_count = integer(0), failed_count = integer(0))
      
      for (i in 1:maxsim){
        DSBN = 0; DSBS = 0; DSBF = 0; RepIndex = 1
        D <- vector("list", ActiveChromatids)
        for (d_idx in seq_len(ActiveChromatids)) D[[d_idx]] <- rep(1, L)
        
        if (i %in% G1){
          cat("\r", "GlobalDSBs:", E1, "/ ChromDSBs:", E, "/ BinRes:", res, "/ Chr:", chrom, "/ WindowWidth:", k,
              "/ Cell:", i, " of ", G, "/ ")
          flush.console()
        }
        
        DSB <- integer(E)
        for (j in 1:E){
          if(i <= RepTN[RepIndex]){
            Dact <- sample.int(ActiveChromatids, 1)
            weighted_probs <- A * D[[Dact]]
            if (sum(weighted_probs) <= 0) next
            DSB[j] <- sample.int(L, 1, replace = TRUE, prob = weighted_probs)
            
            # Determine failure (if enabled)
            is_failed <- FALSE
            if (fail_on) {
              is_failed <- (runif(1) < failrate)
            }
            
            # Update counters
            DSBN <- DSBN + 1
            if (is_failed) DSBF <- DSBF + 1 else DSBS <- DSBS + 1
            
            F[j] <- F[j] + 1
            if(!is_failed){
              H[[j]][DSB[j]] <- H[[j]][DSB[j]] + 1
            }
            # interference window
            C1 <- rep(1, L)
            left <- max(1, DSB[j] - floor(W/2))
            right <- min(L, DSB[j] + floor(W/2))
            
            offsetC <- 1 + (left - (DSB[j] - floor(W/2)))
            offsetC <- as.integer(offsetC)
            offsetC_start <- max(1, offsetC)
            offsetC_end   <- min(W, offsetC + (right - left))
            
            c_sub <- C[offsetC_start:offsetC_end]
            if (length(c_sub) != (right - left + 1)){
              target_len <- (right - left + 1)
              if (length(c_sub) < target_len){
                c_sub <- c(c_sub, rep(tail(c_sub, 1), target_len - length(c_sub)))
              } else {
                c_sub <- c_sub[1:target_len]
              }
            }
            C1[left:right] <- c_sub
            
            # Apply interference — scale the C1 by failStrength for failed attempts so interference strength is adjustable
            scale_factor <- ifelse(is_failed, failStrength, 1)
            
            if (tr > 0) {
              for (ii in setdiff(seq_len(ActiveChromatids), Dact)) {
                # apply scaled interference for trans chromatids
                D[[ii]] <- D[[ii]] * (((C1 * scale_factor) + ((1 / tr) - 1)) * tr)
                D[[ii]] <- pmax(D[[ii]], min_prob_threshold)
              }
            }
            # cis effect on the active chromatid (scaled if failed)
            D[[Dact]] <- D[[Dact]] * (((C1 * scale_factor) + ((1 / cisStrength) - 1)) * cisStrength)
            D[[Dact]] <- pmax(D[[Dact]], min_prob_threshold)
            
            # Snapshot cellcount at replication sequence points
            if (j %in% Repseq & i <= RepTN[RepIndex]) {
              # compute fraction of successes in this cell up to this point
              success_frac_current <- if (DSBN > 0) DSBS / DSBN else 0
              success_count_current <- round(success_frac_current * RepN[RepIndex])
              failed_count_current  <- RepN[RepIndex] - success_count_current
              cellcount <- rbindlist(list(cellcount, data.table(
                DSBs = DSBN,
                GDSBs = RepN[RepIndex],
                Cells = i,
                success_frac = success_frac_current,
                success_count = success_count_current,
                failed_count = failed_count_current
              )), use.names=TRUE)
              RepIndex <- RepIndex + 1
            }
          }
        }
      }
      
      F1 <- cumsum(F) / maxsim
      timer <- toc()
      
      H1 <- vector("list", E)
      H1[[1]] <- H[[1]]
      for (h in 2:E) H1[[h]] <- H1[[h-1]] + H[[h]]
      for (h in seq_along(H1)){
        total <- sum(H1[[h]])
        if (total > 0){
          H1[[h]] <- H1[[h]] / total * 1e6
        } else {
          H1[[h]] <- rep(0, L)
        }
      }
      
      H1=list()
      for (h in 1:E) H1[[h]]=Reduce("+",H[1:h])
      for (h in 1:E) H1[[h]]=H1[[h]]/sum(H1[[h]])*1000000
      
      for(h in Repseq){
        Result=data.frame(NULL)
        Result[1:length(A),"Pos"]=(1:length(A))*res/1000
        Result[1:length(A),"tel1D"]=A
        Result[1:length(B),"WT"]=B
        Result[1:length(H1[[h]]),"sim"]=H1[[h]]
        Result["SimRatio"]=log2(Result["sim"]/Result["tel1D"])
        Result["RealRatio"]=log2(Result["WT"]/Result["tel1D"])
        Result=na.omit(Result)
        Result=Result[!is.infinite(rowSums(Result)),]
        
        SimR=lowess(Result$Pos, Result$SimRatio, f=Smooth)
        RealR=lowess(Result$Pos, Result$RealRatio, f=Smooth)
        DevR=lowess(Result$Pos, Result$SimRatio-Result$RealRatio, f=Smooth)
        rmsd=mean((Result$SimRatio-Result$RealRatio)^2)^0.5
        rmsdA=mean((DevR$y)^2)^0.5
        
        cellsub<- subset(cellcount, cellcount$GDSBs == round_any(h/(E/E1), 50))
        Result$method = winmethod
        Result$window = k
        Result$N.attempted = round_any(h/(E/E1), 20)
        # now include both successful and failed DSB counts (averaged over samples if present)
        Result$successful = if(nrow(cellsub)>0) mean(cellsub$success_count) else 0
        Result$failed     = if(nrow(cellsub)>0) mean(cellsub$failed_count) else 0
        Result$success_frac = if(nrow(cellsub)>0) mean(cellsub$success_frac) else 0
        Result$transstrength = tr
        
        fwrite(Result, paste(sdir ,'/SimTable_', k, '_window_', tr, '_transStrength_', round_any(h/(E/E1), 20), '_attempted_DSBs_part.csv', sep = ''))
      }
    }
  }
  
  files<- list.files(sdir)
  files<- paste(sdir, files, sep = '/')
  AllResult <- NULL
  for (f in files){
    tmp <- fread(f)
    if (is.null(AllResult)) AllResult <- tmp else AllResult <- rbindlist(list(AllResult, tmp), use.names=TRUE, fill=TRUE)
  }
  unlink(paste(sdir, '*', sep = '/'))
  fwrite(AllResult, paste(sdir, '/simtable_', ActiveChromatids, '_chromatids_', G, '_N_', res, '_resolution.csv', sep = ''))
  
  simsum<- NULL
  for(i in levels(as.factor(AllResult$method))){
    tsub1<- subset(AllResult, AllResult$method == i)
    for(j in levels(as.factor(tsub1$window))){
      tsub2<- subset(tsub1, tsub1$window == j)
      for(k in levels(as.factor(tsub2$transstrength))){
        tsub3<- subset(tsub2, tsub2$transstrength == k)
        simsum<- rbind(simsum, data.frame(method = i, window = j, transstrength = k,
                                          successful.DSBs = mean(tsub3$successful, na.rm=TRUE),
                                          failed.DSBs = mean(tsub3$failed, na.rm=TRUE)))
      }
    }
  }
  
  print('Sim generation done')
}
