# Method adapted from Magistro et al 2018
# DOI: 10.1371/journal.pone.0190753
step_adamo = function(signal, apply_bpfilt = T, threshold = 'median', 
                      chunks_samples = 80, bp_low = 0.25, bp_high = 2.5,
                      cross_slope = -1, sf = 128, 
                      locomotion_classifier = F, export_event = F,
                      loco_data = NULL){
  # signal                Tri-axial signal with VM included as a column name in data frame
  # apply_bpfilt          Logical indicating whether bandpass filter should be applied or not. Default is TRUE
  # threshold             The threshold used to determine crossings. Default is median crossing for ADAMO, but can accept absolute numerical inputs such as 0
  # chunks_samples        Number of samples processed at a time
  # bp_low                High-pass filter frequency for constructing bandpass filter. Used in step_bpfilt function
  # bp_high               Low-pass filter frequency for constructing bandpass filter. Used in step_bpfilt function
  # cross_slope           Identify crossings that have a specific slope. Inputs of -1, 0, and 1 correspond to negative, all, and positive slopes
  # sf                    Sampling Frequency
  # locomotion_classifier Logical indicating if steps should only be considered during classified locomotion periods
  # loco_data             External data of locomotion classifications occuring every 1.25 seconds
  
  # TO DO: Paper mentions a 2-tap FIR filter but no info on coefficients. using a BP filter for now
  
  temp_data = data.frame(index = 1:nrow(signal),
                         signal = signal$VM,
                         step = F)
  
  if(locomotion_classifier== T){
    if(is.null(loco_data))
      loco_data = classify_locomotion_kang(signal)
    temp_loco_data = rep(loco_data$activity_state, each = (1.25*sf))
    signal$locomotion = temp_loco_data[1:nrow(signal)]
    
    temp_data$locomotion = signal$locomotion
    
  }
  
  # Applies the custom bandpass filter from ots_methods.R
  if(apply_bpfilt){
    signal = step_bpfilt(signal, bp_low = bp_low, bp_high = bp_high)
    temp_data$signal = signal$VM
  }
  
  
  chunks = seq(1, nrow(temp_data), by = chunks_samples)
  
  for(j in 1:(length(chunks)-1)){
    temp2 = temp_data[chunks[j]:pmin((chunks[j+1]-1), nrow(temp_data)),]
    
    if(locomotion_classifier== T){
      # Check if majority window is locomotion. If not, do not detect steps
      if(sum(temp2$locomotion == 'locomotion', na.rm = T) < (nrow(temp2)/2))
        next
    }
    
    if(threshold == 'median'){
      cross_threshold = (max(temp2$signal)+min(temp2$signal))/2
      temp2$signal = temp2$signal-cross_threshold # offset signal by median threshold so now 0 is equivalent to median
    }
    
    temp2 = temp2 %>% 
      dplyr::mutate(lag = dplyr::lag(signal),
                    lead = dplyr::lead(signal))
    
    temp_cross_index = temp2$index[which((sign(temp2$signal) < sign(temp2$lag)) & # Implement slope requirement
                                           (sign(temp2$signal) <= 0 & sign(temp2$lag) >= 0))] 
    if(length(temp_cross_index)>0){
      temp_data$step[temp_cross_index] = T
    }
  }
  
  temp_data = temp_data %>% dplyr::filter(step == T)
  cross_indices = temp_data$index
  
  # export either the step event indices in seconds or provide the actual timestamps
  if('Timestamp' %in% colnames(signal)){
    export_data = data.frame(index = cross_indices,
                             Timestamp = signal$Timestamp[cross_indices],
                             time = cross_indices/sf,
                             signal_value = temp_data$signal[cross_indices])  
  } else {
    export_data = data.frame(index = cross_indices,
                             time = cross_indices/sf,
                             signal_value = temp_data$signal[cross_indices])  
  }
  
  # If export_event is false, then function will produce second-by-second step counts
  if(export_event == F & ('Timestamp' %in% colnames(signal))){
    export_data = export_data %>% dplyr::mutate(Timestamp_rounded = lubridate::floor_date(Timestamp, unit = 'seconds')) %>%
      dplyr::group_by(Timestamp_rounded) %>%
      dplyr::summarize(steps = n()) %>% dplyr::rename(Timestamp = Timestamp_rounded)
    
    export_data = left_join(data.frame(Timestamp = ymd_hms(floor_date(seq(min(signal$Timestamp, na.rm = T), max(signal$Timestamp, na.rm = T), by = 'secs')))),
                            export_data)
    
    # Fill NA periods in steps with 0
    export_data$steps = nafill(export_data$steps, type = 'const', fill = 0)
  }
  
  return(export_data)
}

# Function to apply bandpass filter to signal
step_bpfilt = function(signal, filt_order = 4, bp_low = 0.25, bp_high = 2.5, sf = 128){
  # signal    Tri-axial accelerometer data with AxisX, AxisY, AxisZ and VM as column names
  # filt_order  Order of the bandpass filter. Default is 4
  # bp_low      High-pass filter frequency of bandpass filter. Default is 2.5 Hz
  # bp_high     Low-pass filter frequency of bandpass filter. Default is 2.5 Hz
  # sf          Sampling frequency of the data. default is 128
  
  bf <- butter(filt_order, c(bp_low, bp_high)/(sf/2), type="pass")
  
  signal$AxisX = filtfilt(bf, signal$AxisX)
  signal$AxisY = filtfilt(bf, signal$AxisY)
  signal$AxisZ = filtfilt(bf, signal$AxisZ)
  signal$VM = filtfilt(bf, signal$VM)
  
  return(signal)
}

# Function to apply Kang 2018 locomotion classification step to data and export locomotion estimates
classify_locomotion_kang = function(signal, sf = 128, 
                                    window_secs = 3.25, # 
                                    slide_secs = 1.25, # selected based on average step duration
                                    min_sdvm_thresh = 0.025,
                                    locofreq_LL = 0.6,
                                    locofreq_UL = 2.0){
  
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
  
  # Locomotion classifier implementation
  for(i in 1:(length(indices)-1)){
    window_start = indices[i]
    
    window_end = ifelse(i != length(indices) & (indices[i] + window_secs*sf) <= nrow(signal), (indices[i] + window_secs*sf), nrow(signal))
    
    sd_vm = sd(signal$VM[window_start:window_end], na.rm = T)
    
    temp_data = signal[window_start:window_end,]
    
    # Stage 1a: Select sensitive axis based on magnitude. ----
    # Since it was developed from Gyro data, all axes have expected value of 0 during non-movement
    # Kang et al identified sensitive axis as the one with the largest angular velocity value
    # Unfortunately, ACC detect gravity so depending on orientation always a non-zero value present
    # Thus, to detect axis most sensitive to movement, we examine the axis that contributes most to VM summary
    # Current Issue: ID'ing axis contributing most to VM is still preferential to axis parallel to gravity vector unless movement is substantially large
    # OR to simplify this process, can use the VM but for now will implement sensitive axis detection
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
    
    # Custom fft function
    # spec = freq_spec(sens_axis_data, scale= F)
    
    # seewave spec function exports frequency in kHz, change to Hz
    
    
    w_0 = mean(spec[which(spec[,1] <locofreq_LL),2], na.rm = T)
    w_c = mean(spec[which(spec[,1] >=locofreq_LL & spec[,1] <=locofreq_UL),2], na.rm = T)
    
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
  }
  return(export_data)
}