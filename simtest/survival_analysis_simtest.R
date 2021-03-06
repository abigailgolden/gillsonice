## Survival analysis simulation test
## Abigail Golden


# Setup -------------------------------------------------------------------

# Clear workspace
rm(list = ls())

options(scipen = 999)

library(tidyverse)
library(survival)
library(survminer)


source("simtest/survival_functions.R")
  
# Set up parameters for data simulation ---------------------------------------------

# set up 50 time steps
t <- 1:50
# Weibull shape and scale parameters; when p <1, baseline hazard decreases over time
p <- 0.5
lambda <- 0.1
# baseline hazard h0(t) is the probability of the event of interest occurring when all covariates are at zero
# here, I'm modeling it as a decreasing Weibull distribution
h0 <- p*lambda*((lambda*t)^(p-1))

plot(t, h0)

# range of winter temperatures in celsius
temp <- seq(-11,0, by = 0.5)
# range of handling times in seconds, from 30 seconds to 10 minutes
hand_time <- seq(30,600, by = 30)

# center these variables so that I can calculate a baseline hazard at the mean for all variables
ctemp <- scale(temp, center = TRUE, scale = TRUE)
chand_time <- scale(hand_time, center = TRUE, scale = TRUE)

# define beta coefficients for each variable and for the interaction between them
btemp <- -0.8
bhand_time <- 0.8
b_int <- 0.5
b_cage <- -0.5  # fish released in predator-control cage are less likely to die


# Simulate data -----------------------------------------------------------

n <- 50
pcage <- 0.4

simulate <- function(n, t, pcage){
  # assign capture-related covariates to fish
  cap <- data.frame(id = 1:n, 
                    ctemp = sample(ctemp, size = n, replace = T),
                    chand_time = sample(chand_time, size = n, replace = T),
                    cage = rbinom(n = n, size = 1, prob = pcage))
  # make output dataframe for survival histories
  dat <- array(data = NA, dim = c(length(t), n))
  # simulate survival histories
  for (i in t){
    # at first time step
    if(i == 1){
      dat[i,] <- rbinom(n,1, 
                       prob = h0[i]*exp(btemp*cap[,2] + bhand_time*cap[,3] + 
                                          b_int*cap[,2]*cap[,3] + b_cage*cap[,4]))
    } else {
      # at subsequent time steps, death depends on whether the fish was dead in the prev time step
      dat[i,] <- ifelse(dat[i-1,] == 1, 1, rbinom(n,1, 
                        prob = h0[i]*exp(btemp*cap[,2] + bhand_time*cap[,3] + 
                                           b_int*cap[,2]*cap[,3] + b_cage*cap[,4])))
    }
  }
  # determine time of death. NA values indicate censoring (aka, death not observed before end of study)
  survtime <- as.data.frame(dat) %>% 
    summarize_all(first.1)
  
  # combine capture and time of death data
  # include indicator column showing whether the data is censored
  data <- cbind(cap, t(survtime))
  colnames(data)[5] <- "survtime"
  data$cens <- ifelse(is.na(data$survtime),0, 1)
  data$survtime <- as.numeric(ifelse(is.na(data$survtime), max(t), survtime))
  return(data)
}


# write function to iterate this process many times for a given sample size and number of timesteps

iterate <- function(nsims, N, timesteps, pcage){
  n_iters <- nsims
  out <- data.frame(ts = rep(NA, n_iters),
                    samplesize = rep(NA, n_iters),
                    b1_sig = rep(NA, n_iters),
                    b2_sig = rep(NA, n_iters),
                    b3_sig = rep(NA, n_iters),
                    b4_sig = rep(NA, n_iters))
  for (j in 1:nsims){
    data <- simulate(n = N, t = timesteps, pcage = pcage)
    mod <- fit.mod(data)
    temp <- test.sig(summary(mod))
    out[j,1] <- max(timesteps)
    out[j,2] <- N
    out[j,3:6] <- temp
  }
  prop_sig <- out %>% 
    group_by(ts, samplesize) %>% 
  summarize(sig0 = (length(which(rowSums(.[,3:6]) == 0)))/nsims,
            sig1 = (length(which(rowSums(.[,3:6]) == 1)))/nsims,
            sig2 = (length(which(rowSums(.[,3:6]) == 2)))/nsims,
            sig3 = (length(which(rowSums(.[,3:6]) == 3)))/nsims,
            sig4 = (length(which(rowSums(.[,3:6]) == 4)))/nsims)
  return(prop_sig)
  
  
    
}

#test <- iterate(nsims = 100, N = )

# and now to run a power analysis across sample sizes

n_power_analysis <- function(samplesize){
  n_iters <- length(samplesize)
  out <- data.frame(timesteps = rep(NA, n_iters),
                    samplesize = rep(NA, n_iters),
                    sig0 = rep(NA, n_iters),
                    sig1 = rep(NA, n_iters),
                    sig2 = rep(NA, n_iters),
                    sig3 = rep(NA, n_iters),
                    sig4 = rep(NA, n_iters))
  for (k in seq_along(samplesize)){
    dat <- iterate(nsims = 100, N = samplesize[k], pcage = pcage)
    out[k,] <- dat
  }
  return(out)
}

# write function to run power analysis across a varying number of timesteps at a given sample size

Ts <- seq(7, 49, by = 7)
samplesize <- 100
timesteps <- Ts

t_power_analysis <- function(timesteps, samplesize){
  n <- samplesize
  n_iters <- length(timesteps)
  out <- data.frame(ts = rep(NA, n_iters),
                    samplesize = rep(NA, n_iters),
                    sig0 = rep(NA, n_iters),
                    sig1 = rep(NA, n_iters),
                    sig2 = rep(NA, n_iters),
                    sig3 = rep(NA, n_iters),
                    sig4 = rep(NA, n_iters))
  for (k in 1:length(timesteps)){
    t <- 1:timesteps[k]
    dat <- iterate(nsims = 100, N = n, timesteps = t, pcage = pcage)
    out[k,] <- dat
  }
  return(out)
}


# Sample size power analysis ----------------------------------------------------

ns<- seq(30, 500, by = 25)

test2 <- n_power_analysis(samplesize = ns, timesteps = max(t))

reshape_test <- pivot_longer(test2, cols = 3:7, names_to = "number_significant",
                             values_to = "prop_sims")

p <- ggplot(reshape_test, aes(x = samplesize, y = prop_sims))+
  geom_line(aes(color = number_significant))+
  labs(x = "Sample size", y = "Proportion of simulations")+
  scale_color_discrete(name = "Number of significant\ncoefficient estimates", labels = c("0", "1", "2", "3", "4"))+
  annotate(geom = "text", x = 300, y = 0.7, label = "h(t) = h(0)exp(b1*temp +\nb2*handling time + b3*cage + b4*temp*handling time)")+
  annotate(geom = "text", x = 100, y = 0.8, label = paste0("Proportion caged:\n", pcage))+
  scale_x_continuous(breaks = seq(0,max(ns), by = 25))+
  theme_classic()

print(p)
 
ggsave("survival_power_analysis_4vars.png", p, width = 10, height = 5, units = "in")



# Time steps power analysis -----------------------------------------------

Ts <- seq(7, 49, by = 7)
n <- 100

tpa <- t_power_analysis(samplesize = n, timesteps = Ts)
tpa_long <- pivot_longer(tpa, cols = 3:7, names_to = "number_significant",
                         values_to = "prop_sims")

p <- ggplot(tpa_long, aes(x = ts, y = prop_sims))+
  geom_line(aes(color = number_significant))+
  labs(x = "Number of timesteps", y = "Proportion of simulations")+
  scale_color_discrete(name = "Number of significant\ncoefficient estimates", labels = c("0", "1", "2", "3", "4"))+
  annotate(geom = "text", x =28, y = 0.3, label = "h(t) = h(0)exp(b1*temp +\nb2*handling time + b3*cage + b4*temp*handling time)")+
  annotate(geom = "text", x = 20, y = 0.8, label = paste0("Proportion caged: ", pcage))+
  annotate(geom = "text", x = 45, y = 0.8, label = paste0("Sample size = ", n))+
  scale_x_continuous(breaks = seq(0, max(Ts), by = 7))+
  theme_classic()

print(p)
