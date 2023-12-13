# Adapted from Ducharme et al 2021
# DOI: 10.1123/jmpb.2021-0011
step_ducharme = function(signal, apply_bpfilt =T, 
                          threshold = .0359, 
                          bp_LL = 0.25,
                          bp_UL = 2.5, 
                          sf = 128, export_reject = F, 
                          locomotion_classifier = F, loco_data = NULL,
                          export_event = F, ...){
  # signal                Single channel acceleration signal provided as a vector
  # threshold             A numeric value indicating the threshold (in gravitational units) used to determine peaks detected
  # sf                    Sampling Frequency
  # apply_bpfilt          Logical indicating whether bandpass filter should be applied or not. Default is TRUE
  # bp_low                High-pass filter frequency for constructing bandpass filter. Used in step_bpfilt function
  # bp_high               Low-pass filter frequency for constructing bandpass filter. Used in step_bpfilt function
  # locomotion_classifier Logical indicating if steps should only be considered during classified locomotion periods
  # loco_data             External data of locomotion classifications occuring every 1.25 seconds
  
  
  temp_data = data.frame(index = seq(1, nrow(signal), by = 1),
                         signal = signal$VM,
                         step = F)
  
  if(locomotion_classifier== T){
    # loco_data = classify_locomotion_kang(signal)
    temp_loco_data = rep(loco_data$activity_state, each = 160)
    signal$locomotion = temp_loco_data[1:nrow(signal)]
    
    temp_data$locomotion = signal$locomotion
    
  } else {
    temp_data$locomotion = NULL
  }
  
  # Applies the custom bandpass filter from ots_methods.R
  if(apply_bpfilt){
    signal = step_bpfilt(signal, bp_low = bp_LL, bp_high = bp_UL)
    temp_data$signal = signal$VM
  }
  
  peak_indices = which(diff(sign(diff(temp_data$signal)))==-2)+1
  
  # export either the step event indices in seconds or provide the actual timestamps
  if('Timestamp' %in% colnames(signal)){
    export_data = data.frame(index = peak_indices,
                             Timestamp = signal$Timestamp[peak_indices],
                             time = peak_indices/sf,
                             signal_value = temp_data$signal[peak_indices])  
  } else {
    export_data = data.frame(index = peak_indices,
                             time = peak_indices/sf,
                             signal_value = temp_data$signal[peak_indices])  
  }
  
  export_data$step = export_data$signal_value >= threshold
  
  # If loco classifier is enabled, keep only loco-related steps and all other peaks stored in reject_step dataframe
  if(locomotion_classifier==T){
    export_data = left_join(export_data, temp_data %>% dplyr::select(index, locomotion))
    reject_steps = export_data %>% dplyr::filter(step == F | locomotion != 'locomotion')
    
    if(export_reject ==F){
      export_data = export_data %>% dplyr::filter(step == T, locomotion == 'locomotion')
    }
    
  } else {
    reject_steps = export_data %>% dplyr::filter(step == F)
    if(export_reject ==F){
      export_data = export_data %>% dplyr::filter(step == T)
    }
  }
  
  # If export_event is false, then function will produce second-by-second step counts
  if(export_event == F & ('Timestamp' %in% colnames(signal))){
    export_data = export_data %>% dplyr::mutate(Timestamp_rounded = lubridate::floor_date(Timestamp, unit = 'seconds')) %>%
      dplyr::group_by(Timestamp_rounded) %>%
      dplyr::summarize(steps = sum(step, na.rm = T)) %>% dplyr::rename(Timestamp = Timestamp_rounded)
    
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