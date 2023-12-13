# Adapted from Kang et al 2018
# DOI: 10.3390/s18010297

step_fft = function(signal, sf = 128, 
                    window_secs = 3.25, # 
                    slide_secs = 1.25, # selected based on average step duration
                    min_sdvm_thresh = 0.025,  # min value of avg spectral amplitude. >=10 Derived from gyro data in seminal paper, needs to be tuned to ACC data. Using pseudo value for now
                    spec_smooth_weight = 0.8,
                    export_suppl = F,
                    export_event = F){
  
  # Stage 0: Identify window indices----
  indices = seq(1, nrow(signal), by = slide_secs*sf)
  
  export_data = data.frame(obs_n = 1:length(indices),
                           index = indices,
                           Timestamp = signal$Timestamp[indices],
                           activity_state = NA,
                           steps = NA)
  
  export_data_suppl = data.frame(obs_n = 1:length(indices),
                                 index = indices,
                                 Timestamp = signal$Timestamp[indices],
                                 window_start =NA,
                                 window_end = NA,
                                 w_0 = NA,
                                 w_c = NA,
                                 sensitive_axis = NA,
                                 sens_ax_criteria = NA,
                                 walk_freq = NA)
  
  # Implement steps 1 (walk detection) and 2 (step count alg)----
  for(i in 1:(length(indices)-1)){
    window_start = indices[i]
    
    window_end = ifelse(i != length(indices) & (indices[i] + window_secs*sf) <= nrow(signal), (indices[i] + window_secs*sf), nrow(signal))
    
    sd_vm = sd(signal$VM[window_start:window_end], na.rm = T)
    
    temp_data = signal[window_start:window_end,]
    
    # Stage 1a: Select sensitive axis based on magnitude. ----
    # Since it was developed from Gyro data, all axes have expected value of 0 during non-movement
    # Kang et al identified sensitive axis as the one with the largest angular velocity value
    # Unfortunately, ACC detect gravity so depending on orientation always a non-zero value present
    # Thus, to detect axis most sensitive to movement, will examine the axis that contributes most to VM summary
    # Current Issue: ID'ing axis contributing most to VM is still preferential to axis parallel to gravity vector unless movement is substantially large
    # OR to simplify this process, can use the VM but for now will implement sensitive axis detection
    # OR if we want to maintain sensitive axis, can use sd to identify which axis demonstrates the greatest variability
    sens_axis = c(mean(abs(temp_data$AxisX)/temp_data$VM, na.rm = T),
                  mean(abs(temp_data$AxisY)/temp_data$VM, na.rm = T),
                  mean(abs(temp_data$AxisZ)/temp_data$VM, na.rm = T))
    
    axis_index = which.max(sens_axis)
    
    sens_axis_data = switch(axis_index, 
                            '1' = temp_data$AxisX,
                            '2' = temp_data$AxisY,
                            '3' = temp_data$AxisZ)
    
    # Detrend sens_axis_data
    sens_axis_data = sens_axis_data-mean(sens_axis_data)
    
    # Variance in signal may be 0 if device was docked resulting in repeated values. If so, skip this window since no relevant data
    if(sd(sens_axis_data) == 0){
      # Update export data
      export_data$activity_state[i] = 'non-locomotion'
      export_data$steps[i] = 0
      
      # Update export supplementary data
      export_data_suppl$window_start[i] = window_start
      export_data_suppl$window_end[i] = window_end
      export_data_suppl$activity_state[i] = 'non-locomotion'
      export_data_suppl$w_0[i] = 0
      export_data_suppl$w_c[i] = 0
      export_data_suppl$sensitive_axis[i] = axis_index
      export_data_suppl$sens_ax_criteria[i] = max(sens_axis, na.rm = T)
      export_data_suppl$walk_freq[i] = 0
      
      next
    }
    
    # Stage 1b: Determine if walking is occurring in window----
    # Identified as the average spectral amplitude between 0.6 to 2 Hz (w_c) is > the average spectral amplitude between 0 to 0.6 Hz (w_0)
    spec = seewave::spec(sens_axis_data, f= sf, fftw = T, flim = c(0, .01),plot = F)
    spec[,1] = spec[,1]*1000
    
    w_0 = mean(spec[which(spec[,1] <0.6),2], na.rm = T)
    w_c = mean(spec[which(spec[,1] >=0.6 & spec[,1] <=2),2], na.rm = T)
    
    activity_state = ifelse(w_0 > w_c, 'non-locomotion', # locomotion-related spectral density is not high enough 
                            ifelse(sd_vm > min_sdvm_thresh, 'locomotion','non-locomotion'))
    
    # Update export supplementary data
    export_data_suppl$window_start[i] = window_start
    export_data_suppl$window_end[i] = window_end
    export_data_suppl$activity_state[i] = activity_state
    export_data_suppl$w_0[i] = w_0
    export_data_suppl$w_c[i] = w_c
    export_data_suppl$sensitive_axis[i] = axis_index
    export_data_suppl$sens_ax_criteria[i] = max(sens_axis, na.rm = T)
    
    if(activity_state == 'non-locomotion'){
      # Update export data
      export_data$activity_state[i] = activity_state
      export_data$steps[i] = 0
      
      next
    }
    
    export_data$activity_state[i] = activity_state
    
    # Step 2: Determine Walking Frequency and determine step count in window----
    curve_data = data.frame(spec[which(spec[,1] >= 0.6 & spec[,1] <=2),])
    colnames(curve_data)= c('freq','amp')
    
    walk_freq = curve_data$freq[which.max(curve_data$amp)]
    
    if(i != 1){
      # Besides the first sample in the data, determine if current sample is the beginning of a walking bout
      # If so, store the walk_freq as is
      if(export_data$activity_state[i-1] != 'locomotion'){ # Beginning of walking bout, can't smooth walking freq from prior datapoint
        steps = round(walk_freq*slide_secs)
        export_data$steps[i] = steps
        
        # Update export supplementary data
        export_data_suppl$walk_freq[i] = walk_freq
      } else {
        # Otherwise if current window is a continuation of walking, smooth the walking freq to minimize spurious spikes in walking frequency
        walk_freq = spec_smooth_weight*export_data_suppl$walk_freq[i-1]+(1-spec_smooth_weight)*walk_freq
        
        steps = round(walk_freq*slide_secs)
        export_data$steps[i] = steps
        
        # Update export supplementary data
        export_data_suppl$walk_freq[i] = walk_freq
      }
      
      
    } else {
      # First obs in dataset, store walk freq and steps estimated
      steps = round(walk_freq*slide_secs)
      export_data$steps[i] = steps
      
      # Update export supplementary data
      export_data_suppl$walk_freq[i] = walk_freq
    }
  }
  
  # If export_event is false, then function will produce second-by-second step counts
  if(export_event == F & ('Timestamp' %in% colnames(signal))){
    export_data = export_data %>% dplyr::mutate(Timestamp_rounded = lubridate::floor_date(Timestamp, unit = 'seconds')) %>%
      dplyr::group_by(Timestamp_rounded) %>%
      dplyr::summarize(steps = sum(steps)) %>% dplyr::rename(Timestamp = Timestamp_rounded)
    
    export_data = left_join(data.frame(Timestamp = ymd_hms(floor_date(seq(min(signal$Timestamp, na.rm = T), max(signal$Timestamp, na.rm = T), by = 'secs')))),
                            export_data)
    
    # Fill NA periods in steps with 0
    export_data$steps = nafill(export_data$steps, type = 'const', fill = 0)
  }
  # Wrap up with supplementary export data----
  if(export_suppl == T){
    export_data = left_join(export_data, export_data_suppl)
  }
  
  return(export_data)

}