measure_metrics_with_antennas <- function (ind_movement){
  require (dplyr)
  require (data.table)
  
  ant_data <- list ()
  last_comp <- nrow(ind_movement)-1
  for (i in 1:last_comp){
    range_movement <- range (ind_movement[i,'cum_movement'],
                             ind_movement[i+1,'cum_movement'])
    range_mov_ant <- data.frame (anten_location = sort(anten_location),
                                 initial = range_movement[1],
                                 final = range_movement[2])
    value <- range_mov_ant$anten_location
    range_mov_ant <- range_mov_ant %>% 
      mutate(detection = value >= initial  & value <= final)
    range_mov_ant$comp <- i
    ant_data [[i]] <- range_mov_ant
  }
  
  ant_data <- do.call(rbind.data.frame, ant_data)
  
  #Total distance
  detected <- subset (ant_data, detection == TRUE)
  for (i in 1:nrow(detected)){
    detected$dist_to_initial[i] <- dist(c(detected$anten_location[i],
                                          detected$initial[i]))}
  detected_simplified <- setDT (detected)[, .SD[which.max(dist_to_initial)],by=comp]
  distances <- vector ()
  for (i in 1:nrow(detected_simplified)){
    distances[i] <- dist (c(detected_simplified$anten_location[i],
                            detected_simplified$anten_location[i+1]))
  }
  total_distance <- sum (distances, na.rm = TRUE)
  
  # Amount of time without move
  DidItMove <- aggregate (detection~comp,data=ant_data, FUN=any)
  Amount_time_Wmove <- table(DidItMove$detection)[1]
  
  # CV of distance tracked
  Zero_vector <- rep (0, Amount_time_Wmove)
  CV_movement <- sd (c(distances, Zero_vector),na.rm = TRUE)/
    mean (c(distances, Zero_vector),na.rm = TRUE) * 100
  
  # Average of of distance tracked 
  Mean_movement <- mean (c(distances, Zero_vector),na.rm = TRUE)
  
  # Last time recorded
  last_time_recorded <- subset (DidItMove, detection == 'TRUE')
  last_time_recorded <- tail (last_time_recorded[order (last_time_recorded$comp),],1)[,1]
  last_time_recorded <- n_days - last_time_recorded  
  
  metrics <- c(total_distance, Amount_time_Wmove, Mean_movement, 
               CV_movement, last_time_recorded)
  names(metrics) <- c('Total_Distance', 'Amount_Time_Wmove', 'Mean_Movement', 
                      'CV_Movement', 'Last_Time_Recorded')
  return (metrics)
  }

measure_metrics_with_antennas(ind_movement)


