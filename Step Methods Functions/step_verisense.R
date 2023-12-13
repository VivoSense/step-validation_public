# Adapted from Maylor et al 2022
# DOI: 10.3390/s22249984
step_verisense <- function(input_data = signal, export_event = T,
                           sf = 128, 
                           k = 4, # Number of samples considered per peak
                           period_min = 4, # 
                           period_max = 20, # 
                           sim_thres= -1,# similarity threshold
                           cont_win_size = 4,  # continuity window size
                           cont_thres = 4,     # continuity threshold
                           var_thres = .01,# variance threshold
                           mag_thres = 1.25, # threshold on VM magnitude to consider a step
                           mag_thres1_sdacc = 0.025) {
  # adapted from https://github.com/Maylor8/Verisense-Toolbox/blob/master/Verisense_step_algorithm/verisense_count_steps.R
  acc = input_data$VM
  
  if (sd(acc) < mag_thres1_sdacc) {
    # acceleration too low, no steps
    num_seconds = round(length(acc) / sf)
    steps_per_sec = rep(0,num_seconds)
  } else {
    # Search for steps
    
    # find the peak rms value is every range of k
    half_k <- round(k/2)
    segments <- floor(length(acc) / k)
    peak_info <- matrix(NA,nrow=segments,ncol = 5)
    # peak_info[,1] - peak location
    # peak_info[,2] - acc magnitude
    # peak_info[,3] - periodicity (samples)
    # peak_info[,4] - similarity
    # peak_info[,5] - continuity
    
    # for each segment find the peak location
    for (i in 1:segments) {
      start_idx <- (i-1) * k + 1
      end_idx <- start_idx + (k-1)
      tmp_loc_a <- which.max(acc[start_idx:end_idx])
      tmp_loc_b <- (i-1) * k + tmp_loc_a
      # only save if this is a peak value in range of -k/2:+K/2
      start_idx_ctr <- tmp_loc_b - half_k
      if (start_idx_ctr < 1) {
        start_idx_ctr <- 1
      }
      end_idx_ctr <- tmp_loc_b + half_k
      if (end_idx_ctr > length(acc)) {
        end_idx_ctr <- length(acc)
      }
      check_loc <- which.max(acc[start_idx_ctr:end_idx_ctr])
      if (check_loc == (half_k + 1)) {
        peak_info[i,1] <- tmp_loc_b
        peak_info[i,2] <- max(acc[start_idx:end_idx])
      }
    }
    peak_info <- peak_info[is.na(peak_info[,1])!=TRUE,] # get rid of na rows
    peak_info_screen = data.frame(peak_info[NULL,c(1:3)]) # will always use column 3 for ID reason for exclusion

    peak_info <-peak_info[peak_info[,2] > mag_thres,] 
    
    if (length(peak_info) > 10) {  # there must be at least two steps
      num_peaks <- length(peak_info[,1])
      
      no_steps = FALSE
      if (num_peaks > 2) {
        # Calculate Features (periodicity, similarity, continuity)
        peak_info[1:(num_peaks-1),3] <- diff(peak_info[,1]) # calculate periodicity - RM for VivoSense modification to include
        
        peak_info <- peak_info[peak_info[,3] > period_min,] # filter peaks based on period_min

        peak_info <- peak_info[peak_info[,3] < period_max,]   # filter peaks based on period_max 
      } else {
        no_steps = TRUE
      }
    } else {
      no_steps = TRUE
    }
    
    if ( length(peak_info)==0 || length(peak_info) == sum(is.na(peak_info)) || no_steps == TRUE) {
      # no steps found
      num_seconds = round(length(acc) / sf)
      steps_per_sec = rep(0,num_seconds)
    } else {
      # calculate similarity
      num_peaks <- length(peak_info[,1])
      peak_info[1:(num_peaks-2),4] <- -abs(diff(peak_info[,2],2)) # calculate similarity
      
      peak_info <- peak_info[peak_info[,4] > sim_thres,]  # filter based on sim_thres
      peak_info <- peak_info[is.na(peak_info[,1])!=TRUE,] # previous statement can result in an NA in col-1
      
      # calculate continuity
      if (length(peak_info[,3]) > 5) {
        end_for <- length(peak_info[,3])-1
        for (i in cont_thres:end_for) {
          # for each bw peak period calculate acc var
          v_count <- 0 # count how many windows were over the variance threshold
          for (x in 1:cont_thres) {
            if (var(acc[peak_info[i-x+1,1]:peak_info[i-x+2,1]]) > var_thres) {
              v_count = v_count + 1
            }
          }
          if (v_count >= cont_win_size) {
            peak_info[i,5] <- 1 # set continuity to 1, otherwise, 0
          } else {
            peak_info[i,5] <- 0
          }
        }
      } 
      
      peak_info <- peak_info[peak_info[,5]==1,1] # continuity test - only keep locations after this
      peak_info <- peak_info[is.na(peak_info)!=TRUE] # previous statement can result in an NA in col-1
      
      if (length(peak_info)==0) {
        # no steps found
        num_seconds = round(length(acc) / sf)
        steps_per_sec = rep(0,num_seconds)
      } else {
        
        # debug plot
        # is_plot = F
        # if (is_plot) {
        #   library(ggplot2)
        #   library(plotly)
        #   acc.df <- data.frame(acc=acc, det_step=integer(length(acc)))
        #   acc.df$det_step[peak_info] <- 1  # to plot annotations, prepare a 0/1 column on dataframe
        #   acc.df$idx <- as.numeric(row.names(acc.df))
        #   pl <- ggplot(data=acc.df,aes(x=idx,y=acc)) 
        #   pl2 <- pl + geom_line()
        #   pl3 <- pl2 + geom_point(data=subset(acc.df,det_step==1),aes(x=idx,y=acc),color='red',size=1,alpha=0.7)
        #   pl4 <- ggplotly(pl3)
        #   print(pl4)  
        # }
        
        # for GGIR, output the number of steps in 1 second chunks
        start_idx_vec <- seq(from=1,to=length(acc),by=sf)
        steps_per_sec <- table(factor(findInterval(peak_info, start_idx_vec), levels = seq_along(start_idx_vec)))
        steps_per_sec <- as.numeric(steps_per_sec)
      }
    }
  }
  
  peak_info_screen = peak_info_screen %>% dplyr::rename(index = X1,
                                                        magnitude = X2,
                                                        reason = X3) %>%
    dplyr::mutate(steps = F)
  
  peak_info_screen$Timestamp = input_data$Timestamp[peak_info_screen$index]
  # Added by RM for VivoSense----
  if(export_event == T){
    export_data = data.frame(index = peak_info,
                             Timestamp = input_data$Timestamp[peak_info],
                             steps = TRUE)
    
    export_data = bind_rows(export_data, peak_info_screen %>% dplyr::select(index, steps,reason))
  } else {
    # Append Timestamp to export data
    export_data = data.frame(index = 1:length(steps_per_sec),
                             Timestamp = input_data$Timestamp[seq(1, nrow(input_data), by = sf)[1:length(steps_per_sec)]],
                             steps = steps_per_sec)
  }
  return(export_data)
}