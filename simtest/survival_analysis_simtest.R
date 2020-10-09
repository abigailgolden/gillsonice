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

# write function to iterate this process many times for a given sample size

#ts <- seq(10,50, by = 5)




iterate <- function(nsims, N){
  n_iters <- nsims
  out <- data.frame(timesteps = rep(NA, n_iters),
                    samplesize = rep(NA, n_iters),
                    b1_sig = rep(NA, n_iters),
                    b2_sig = rep(NA, n_iters),
                    b3_sig = rep(NA, n_iters))
  for (j in 1:nsims){
    data <- simulate(n = N, t = t)
    mod <- fit.mod(data)
    temp <- test.sig(summary(mod))
    out[j,1] <- max(t)
    out[j,2] <- N
    out[j,3:5] <- temp
  }
  prop_sig <- out %>% 
    group_by(timesteps, samplesize) %>% 
  summarize(sig0 = (length(which(rowSums(.[,3:5]) == 0)))/nsims,
            sig1 = (length(which(rowSums(.[,3:5]) == 1)))/nsims,
            sig2 = (length(which(rowSums(.[,3:5]) == 2)))/nsims,
            sig3 = (length(which(rowSums(.[,3:5]) == 3)))/nsims)
  return(prop_sig)
  
  
    
}

power_analysis <- function(samplesize){
  n_iters <- length(samplesize)
  out <- data.frame(timesteps = rep(NA, n_iters),
                    samplesize = rep(NA, n_iters),
                    sig0 = rep(NA, n_iters),
                    sig1 = rep(NA, n_iters),
                    sig2 = rep(NA, n_iters),
                    sig3 = rep(NA, n_iters))
  for (k in seq_along(samplesize)){
    dat <- iterate(nsims = 100, N = samplesize[k])
    out[k,] <- dat
  }
  return(out)
}

ns<- seq(15, 100, by = 5)

test2 <- power_analysis(samplesize = ns)

reshape_test <- pivot_longer(test2, cols = 3:6, names_to = "number_significant",
                             values_to = "prop_sims")

p <- ggplot(reshape_test, aes(x = samplesize, y = prop_sims))+
  geom_line(aes(color = number_significant))+
  labs(x = "Sample size", y = "Proportion of simulations")+
  scale_color_discrete(name = "Number of significant\ncoefficient estimates", labels = c("0", "1", "2", "3"))+
  annotate(geom = "text", x = 50, y = 0.9, label = "h(t) = h(0)exp(b1*temp +\nb2*handling time + b3*temp*handling time)")+
  theme_classic()
 
ggsave("survival_power_analysis.png", p)
