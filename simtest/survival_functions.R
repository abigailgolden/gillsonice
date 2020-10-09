# Define functions --------------------------------------------------------

# count the first positive occurrence of a binary variable (in this case, alive vs dead) in a time series
first.1 <- function(x){
  match(1, x)
}

## Fit a cox proportional hazards model

fit.mod <- function(d){
  fit <- coxph(Surv(survtime, cens) ~ ctemp + chand_time + ctemp*chand_time, data = d)
  return(fit)
  
}

#see if the model detects significant coefficient values for temp and handling time
test.sig <- function(model_obj){
  fit <- model_obj
  stemp <- ifelse(fit$coefficients[1,5] < 0.05, 1,0)
  shand_time <- ifelse(fit$coefficients[2,5] < 0.05, 1,0)
  s_int <- ifelse(fit$coefficients[3,5] < 0.05, 1,0)
  sig <- c(stemp, shand_time, s_int)
  return(sig)
}

# calculate proportion of total simulations

prop <- function(x){
  x/nsims
}
