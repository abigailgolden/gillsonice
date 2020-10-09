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
btemp <- -0.75
bhand_time <- 0.8
b_int <- 0.5


# Simulate data -----------------------------------------------------------

#n <- seq(10,100, by = 5)
n <- 50

simulate <- function(n, t){
  # assign capture-related covariates to fish
  cap <- data.frame(id = 1:n, 
                    ctemp = sample(ctemp, size = n, replace = T),
                    chand_time = sample(chand_time, size = n, replace = T))
  # make output dataframe for survival histories
  dat <- array(data = NA, dim = c(length(t), n))
  # simulate survival histories
  for (i in t){
    if(i == 1){
      dat[i,] <- rbinom(n,1, 
                       prob = h0[i]*exp(btemp*cap[,2] + bhand_time*cap[,3] + 
                                          b_int*cap[,2]*cap[,3]))
    } else {
      dat[i,] <- ifelse(dat[i-1,] == 1, 1, rbinom(n,1, 
                        prob = h0[i]*exp(btemp*cap[,2] + bhand_time*cap[,3] + 
                                           b_int*cap[,2]*cap[,3])))
    }
  }
  # determine time of death. NA values indicate censoring (aka, death not observed before end of study)
  survtime <- as.data.frame(dat) %>% 
    summarize_all(first.1)
  
  # combine capture and time of death data
  # include indicator column showing whether the data is censored
  data <- cbind(cap, t(survtime))
  colnames(data)[4] <- "survtime"
  data$cens <- ifelse(is.na(data$survtime),0, 1)
  data$survtime <- as.numeric(ifelse(is.na(data$survtime), max(t), survtime))
  return(data)
}

simdat <- simulate(n = 100, t = t)

test.sig(fit.mod(simdat))

# write function to iterate this process across many sample sizes
Ns <- seq(10, 50, by = 5)
Ts <- seq(10,50, by = 5)


iterate <- function(Ns, Ts){
 
  n_iters <- length(Ns)*length(Ts) 
  out <- data.frame(timesteps = rep(NA, n_iters),
                    samplesize = rep(NA, n_iters),
                    b1_sig = rep(NA, n_iters),
                    b2_sig = rep(NA, n_iters),
                    b3_sig = rep(NA, n_iters))
  r <- 1
  
  for (k in 1:length(Ts)){
    t <- seq(1, Ts[k], by = 1)
      for (l in 1:length(Ns)){
        data <- simulate(n = Ns[l], t = t)
        mod <- fit.mod(data)
        temp <- test.sig(summary(mod))
        out[r,1] <- Ts[k]
        out[r,2] <- Ns[l]
        out[r,3:5] <- temp
        r <- r + 1
    }
  }
 
  return(out)
}

power_analysis <- iterate(Ns = Ns, Ts = Ts)
