## Survival analysis simulation test
## Abigail Golden


# Setup -------------------------------------------------------------------

# Clear workspace
rm(list = ls())

options(scipen = 999)

library(tidyverse)
library(survival)
library(survminer)


# Define functions --------------------------------------------------------

# count the first positive occurrence of a binary variable (in this case, alive vs dead) in a time series
first.1 <- function(x){
  match(1, x)

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
btemp <- -0.5
bhand_time <- 0.8
b_int <- 0.1


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
      dat[i,] <- ifelse(t > 1 & dat[i-1,] == 1, 1, rbinom(n,1, 
                        prob = h0[i]*exp(btemp*cap[,2] + bhand_time*cap[,3] + 
                                           b_int*cap[,2]*cap[,3])))
    }
  }
  # determine time of death. NA values indicate censoring (aka, death not observed before end of study)
  survtime <- as.data.frame(dat) %>% 
    summarize_all(first.1)
  
  # combine capture and time of death data
  data <- cbind(cap, t(survtime))
  colnames(data)[4] <- "survtime"
  return(data)
}

simdat <- simulate(n = n, t = t)

## Fit survival model and see if it detects significant coefficient values for temp and handling time


# write function to iterate this across many sample sizes