# Discriminant function analysis

# Libraries
require (MASS)
require (truncnorm)
require (GGally)
require (lattice)

# functions
confusion <- function(actual, predicted, names = NULL, printit = TRUE,
                      prior = NULL) {
  if (is.null(names))
    names <- levels(actual)
  tab <- table(actual, predicted)
  acctab <- t(apply(tab, 1, function(x) x/sum(x)))
  dimnames(acctab) <- list(Actual = names, "Predicted (cv)" = names)
  if (is.null(prior)) {
    relnum <- table(actual)
    prior <- relnum/sum(relnum)
    acc <- sum(tab[row(tab) == col(tab)])/sum(tab)
  }
  else {
    acc <- sum(prior * diag(acctab))
    names(prior) <- names
  }
  if (printit)
    print(round(c("Overall accuracy" = acc, "Prior frequency" = prior),
                4))
  if (printit) {
    cat("\nConfusion matrix", "\n")
    print(round(acctab, 4))
  }
  invisible(acctab)
} # Function to create confusion matrices

# Simulate data
df <- data.frame (Fate = factor(rep(c('Dead','Alive'),each=15)),
                  Home_range = NA,
                  Time_still = NA,
                  Lst_tm_rec = NA)

for (i in 1:nrow(df)){
  if (df [i,"Fate"]=="Alive"){
    df [i,"Home_range"] <- rtruncnorm(n=1, a=0, b=50, mean = 20, sd = 20) 
    df [i,"Time_still"] <- rtruncnorm(n=1, a=0, b=30, mean = 3, sd = 1) 
    df [i,"Lst_tm_rec"] <- rtruncnorm(n=1, a=0, b=30, mean = 3, sd = 3) 
    } 
  else {
    df [i,"Home_range"] <- rtruncnorm(n=1, a=0, b=50, mean = 7, sd = 7)  
    df [i,"Time_still"] <- rtruncnorm(n=1, a=0, b=30, mean = 6, sd = 2) 
    df [i,"Lst_tm_rec"] <- rtruncnorm(n=1, a=0, b=30, mean = 5, sd = 5)
    }
}

#Association between variables
ggpairs(df, upper = list(continuous = "points", combo = "box"), 
        lower = list(continuous = "points", combo = "box"),
        ggplot2::aes(colour=Fate))

#Model
fit <- lda(Fate ~ ., data=df,CV=FALSE) #All data
fit
fit.hat <- predict(fit)
confusion(df$Fate,fit.hat$class)

densityplot(~fit.hat$x, groups = df$Fate) # Plot discriminant function
plot(fit) # Plot discriminant function

#leave-one-out cross-validation - Similar to jack-knife estimation (Efron, 1982)
fit_cv <- lda(Fate ~ ., data=df,CV=TRUE) 
confusion(df$Fate,fit_cv$class)
