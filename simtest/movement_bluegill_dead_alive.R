# Packages
require (fishmove)
require (ggplot2)
require (truncnorm)
require (MASS)
require (GGally)
require (randomForest)

# Libraries

# functions 
source ('measure_metrics_with_antennas.R')

# Parameters
n_days = 30
n_individuals = 500
anten_location = c(0,
                   -15.24, 15.24,
                   -30.48, 30.48,
                   -76.2, 76.2,
                   -152.4,152.4) # distance from central point (meters)
avg_size <- 180  # mm
sd_size <- 15 
body_sizes <- rnorm(n=n_individuals, mean = avg_size, sd = sd_size)
plot (density(
  body_sizes)) # visualize body size distribution
p_stat <- 0.67 # 0.67 default 
stream_order <- 6

avg_movement_dead_fish <- 0
sd_movement_dead_fish <- 10
plot (density(
rnorm(n=n_individuals, 
      mean = avg_movement_dead_fish, 
      sd = sd_movement_dead_fish))) #visualize movement of dead fish

#general scale movement
bluegill_30days <- fishmove(species="Lepomis macrochirus",T=30,
                     L=avg_size, SO=stream_order)
tiff(file=paste ('move30d_size',avg_size,'.tiff',sep=''),
     width=2000,height=2000,units = "px",
     compression="lzw",res=300)
pdk (bluegill_30days)
dev.off()


mov_upper95<- bluegill_30days$pred.fishmove[6] #sigma_mob upper

#individual simulation 
results <- list ()
for (ind_i in 1:length(body_sizes)){
  bluegill <- fishmove(species="Lepomis macrochirus",T=1,
                       L=body_sizes[ind_i], SO=stream_order)
  distance <- seq(from=0, to=1000, by=1)
  cum_prob <- rep(NA,length(distance))
  
  for (dist_i in 1:length(distance)){
    cum_prob[dist_i] <- fishmove.query (bluegill,dist=distance[dist_i],p=p_stat)
  } #calculate cumulative probability
  
  probs <- rep (NA, length(cum_prob))
  
  for (prob_i in 1:length(cum_prob)){
    cum_prob2 <- c(cum_prob,0)
    probs[prob_i] <- cum_prob [prob_i] - cum_prob2 [prob_i+1]
  } #calculate probability of each case
  
  distance <- c(distance*-1,distance)
  probs <- rep (probs,2)
  
  movement_fish_alive <- sample (x = distance, size = n_days,
                                 rep = TRUE, prob = probs)
  movement_fish_dead <- rnorm(n=n_days, 
                              mean = avg_movement_dead_fish, 
                              sd = sd_movement_dead_fish)
  chance_to_disapear <-sample (c(1,0),size=n_days,
                               prob=c(0.9,0.1),replace=T)
  disappear <- FALSE
  for (p in 1:length(chance_to_disapear)){
    if (disappear == FALSE){
      chance_to_disapear[p] <- chance_to_disapear[p]
      if (chance_to_disapear[p]==0){
        disappear <- TRUE
        movement_fish_dead[p] <- 0}
      else {
        movement_fish_dead[p] <- movement_fish_dead[p]
      }
    } else {
      movement_fish_dead[p] <- 0
    }
    }

  cum_movement_alive <- c(0,
    cumsum(movement_fish_alive))
  
  cum_movement_dead <- c(0,
                         cumsum(movement_fish_dead))
  
  movement_fish_alive <- c(0,movement_fish_alive)
  movement_fish_dead <- c(0,movement_fish_dead)
  
  simulations <- data.frame (ind = rep (NA, n_days+1),
                             time = rep (NA, n_days+1),
                             size = rep (NA, n_days+1),
                             cum_movement_alive = rep (NA, n_days+1),
                             movement_alive = rep (NA, n_days+1),
                             cum_movement_dead = rep (NA, n_days+1),
                             movement_dead = rep (NA, n_days+1)
                             )
  simulations$ind <- ind_i
  simulations$size <- body_sizes[ind_i]
  simulations$time <- 0:n_days
  simulations$movement_alive <- movement_fish_alive
  simulations$cum_movement_alive <- cum_movement_alive
  simulations$movement_dead <- movement_fish_dead
  simulations$cum_movement_dead <- cum_movement_dead
  
  results[[ind_i]] <- simulations
  
}

tiff(file=paste ('sim_alive_size',avg_size,'.tiff',sep=''),
     width=2000,height=2000,units = "px",
     compression="lzw",res=300)
par(mai=c(1,1,1,1))
plot(0,type='n',
     ylim=c(-1.05*mov_upper95,1.05*mov_upper95),xlim=c(1,n_days),
     ylab="Distance (m)",xlab="Time (days)",main=paste('p=',p_stat,
                                                       ' n=',n_individuals,
                                                       ' size=', avg_size,
                                                       ' river order=', stream_order),
     cex.lab=1.2,cex.axis=1.1, cex.main=1.2)
rect(xleft=0, ybottom=anten_location[1], xright=30, ytop=anten_location[3],
     border = "aquamarine3", col = "aquamarine3")
rect(xleft=0, ybottom=anten_location[2], xright=30, ytop=anten_location[1],
     border = "azure3", col = "azure3")  
rect(xleft=0, ybottom=anten_location[3], xright=30, ytop=anten_location[5],
     border = "azure3", col = "azure3") 
rect(xleft=0, ybottom=anten_location[4], xright=30, ytop=anten_location[2],
     border = "aquamarine3", col = "aquamarine3")
rect(xleft=0, ybottom=anten_location[5], xright=30, ytop=anten_location[7],
     border = "aquamarine3", col = "aquamarine3")
rect(xleft=0, ybottom=anten_location[6], xright=30, ytop=anten_location[4],
     border = "azure3", col = "azure3") 
rect(xleft=0, ybottom=anten_location[7], xright=30, ytop=anten_location[9],
     border = "azure3", col = "azure3") 
rect(xleft=0, ybottom=anten_location[8], xright=30, ytop=anten_location[6],
     border = "aquamarine3", col = "aquamarine3")

metrics_alive <- data.frame (Total_Distance= numeric(), Amount_Time_Wmove= numeric(),
                             Mean_Movement= numeric(), CV_Movement= numeric(),
                             Last_Time_Recorded= numeric())
for (ind in 1:length(results)){
  ind_movement <- results [[ind]]
  ind_movement <- ind_movement [,c(1:5)]
  colnames (ind_movement)[4:5] <- c('cum_movement','movement')
  lines (ind_movement$cum_movement~ind_movement$time)
  metrics_alive [ind,] <- measure_metrics_with_antennas (ind_movement)
}
dev.off()

tiff(file=paste ('sim_dead',avg_size,'.tiff',sep=''),
     width=2000,height=2000,units = "px",
     compression="lzw",res=300)
par(mai=c(1,1,1,1))
plot(0,type='n',
     ylim=c(-1.05*mov_upper95,1.05*mov_upper95),xlim=c(1,n_days),
     ylab="Distance (m)",xlab="Time (days)",main=paste('p=',p_stat,
                                                       ' n=',n_individuals,
                                                       ' size=', avg_size,
                                                       ' river order=', stream_order),
     cex.lab=1.2,cex.axis=1.1, cex.main=1.2)
rect(xleft=0, ybottom=anten_location[1], xright=30, ytop=anten_location[3],
     border = "aquamarine3", col = "aquamarine3")
rect(xleft=0, ybottom=anten_location[2], xright=30, ytop=anten_location[1],
     border = "azure3", col = "azure3")  
rect(xleft=0, ybottom=anten_location[3], xright=30, ytop=anten_location[5],
     border = "azure3", col = "azure3") 
rect(xleft=0, ybottom=anten_location[4], xright=30, ytop=anten_location[2],
     border = "aquamarine3", col = "aquamarine3")
rect(xleft=0, ybottom=anten_location[5], xright=30, ytop=anten_location[7],
     border = "aquamarine3", col = "aquamarine3")
rect(xleft=0, ybottom=anten_location[6], xright=30, ytop=anten_location[4],
     border = "azure3", col = "azure3") 
rect(xleft=0, ybottom=anten_location[7], xright=30, ytop=anten_location[9],
     border = "azure3", col = "azure3") 
rect(xleft=0, ybottom=anten_location[8], xright=30, ytop=anten_location[6],
     border = "aquamarine3", col = "aquamarine3")

metrics_dead <- data.frame (Total_Distance= numeric(), Amount_Time_Wmove= numeric(),
                             Mean_Movement= numeric(), CV_Movement= numeric(),
                             Last_Time_Recorded= numeric())
for (ind in 1:length(results)){
  ind_movement <- results [[ind]]
  ind_movement <- ind_movement [,c(1:3,6:7)]
  colnames (ind_movement)[4:5] <- c('cum_movement','movement')
  lines (ind_movement$cum_movement~ind_movement$time)
  metrics_dead [ind,] <- measure_metrics_with_antennas (ind_movement)
}
dev.off()

df <- data.frame (Fate = rep(c('alive','dead'),each=n_individuals),
                  Total_Distance = c(metrics_alive$Total_Distance,
                                     metrics_dead$Total_Distance),
                  Amount_Time_Wmove = c(metrics_alive$Amount_Time_Wmove,
                                        metrics_dead$Amount_Time_Wmove),
                  Mean_Movement = c(metrics_alive$Mean_Movement,
                                    metrics_dead$Mean_Movement),
                  CV_Movement = c(metrics_alive$CV_Movement,
                                        metrics_dead$CV_Movement),
                  Last_Time_Recorded = c(metrics_alive$Last_Time_Recorded,
                                         metrics_dead$Last_Time_Recorded))

tiff(file=paste ('corr',avg_size,'.tiff',sep=''),
     width=3000,height=3000,units = "px",
     compression="lzw",res=300)
ggpairs(df, upper = list(continuous = "points", combo = "box"), 
        lower = list(continuous = "points", combo = "box"),
        ggplot2::aes(colour=Fate))
dev.off()

# Checking model accuracy

# Check sample effect
df$Fate <- ifelse (df$Fate=='alive',1,0)
df$Fate <- factor(df$Fate)

df <- df[complete.cases(df), ] #remove NAs

results_n <- data.frame (n= numeric(), Accuracy_lda = numeric(),
                         Accuracy_rf= numeric())
count <- 1
n_begin <- 3
n_end <- 100
for (n in n_begin:n_end){
  # lda
  df_temp <- do.call( rbind, lapply( split(df, df$Fate) ,
                                     function(df) df[sample(nrow(df), n) , ] )
  )
  df_temp <- df_temp[,-4]
  fit_cv_temp <- lda(Fate ~ Total_Distance+Amount_Time_Wmove+
                       CV_Movement+Last_Time_Recorded, 
                     data=df_temp,CV=TRUE)
  ct <- table(df_temp$Fate, fit_cv_temp$class)
  results_n [count,1]<- n
  results_n [count,2]<- sum(diag(prop.table(ct)))
  
  # Random forest
  rf_pred <- c()
  n_total <- n*2
  id <- 1:n_total
  for(k in id){
    #loo
    rf <- randomForest(Fate ~ ., data=df_temp[id!=k,])
    rf_pred[k] <- predict(rf, newdata=df_temp[id==k,])
    
  }
  rf_ct <- table(df_temp$Fate, rf_pred)
  results_n [count,3]<- sum(diag(prop.table(rf_ct)))
  count <- count + 1
}

# just doing some rearranging for the plot
results_n <- data.frame (n= rep (results_n$n, 2),
                         Accuracy=  c(results_n$Accuracy_lda,
                                      results_n$Accuracy_rf),
                         model= c(rep ('lda',nrow(results_n)),
                                  rep ('rf',nrow(results_n))))

#plot 
tiff(file=paste ('accuracy_n_relationship',avg_size,'.tiff',sep=''),
     width=2000,height=2000,units = "px",
     compression="lzw",res=300)
ggplot (results_n ,aes(x=n, y=Accuracy, color=model))+
  geom_point(size=2)+geom_smooth(method="loess", alpha=0.1, 
                                 aes(fill = model))+
  theme_bw()
dev.off()

