# Discriminant function analysis

# Libraries
require (MASS)
require (truncnorm)
require (GGally)
require (lattice)
require (randomForest)

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

df <- data.frame (Fate = factor(rep(c('Dead','Alive'),each=1000)),
                  Home_range = NA,
                  Time_still = NA,
                  Lst_tm_rec = NA)

# Simulate data
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

ggpairs(df, upper = list(continuous = "points", combo = "box"), 
        lower = list(continuous = "points", combo = "box"),
        ggplot2::aes(colour=Fate))

#Model
fit <- lda(Fate ~ ., data=df,CV=FALSE) #All data
fit.hat <- predict(fit)
conf_matrix<-confusion(df$Fate,fit.hat$class)

densityplot(~fit.hat$x, groups = df$Fate) # Plot discriminant function
plot(fit) # Plot discriminant function

#leave-one-out cross-validation - Similar to jack-knife estimation (Efron, 1982)
fit_cv <- lda(Fate ~ ., data=df,CV=TRUE) 
confusion(df$Fate,fit_cv$class)

# Check sample effect
results_n <- data.frame (n= numeric(), Accuracy_lda = numeric(),
                         Accuracy_rf= numeric())
count <- 1
n_begin <- 5
n_end <- 100
for (n in n_begin:n_end){
  # lda
  df_temp <- df[sample(nrow(df), n), ]
  fit_cv_temp <- lda(Fate ~ ., data=df_temp,CV=TRUE)
  ct <- table(df_temp$Fate, fit_cv_temp$class)
  results_n [count,1]<- n
  results_n [count,2]<- sum(diag(prop.table(ct)))
  
  # Random forest
  rf_pred <- c()
  id <- 1:n
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
ggplot (results_n ,aes(x=n, y=Accuracy, color=model))+
  geom_point(size=2)+geom_smooth(method="loess", alpha=0.1, 
                           aes(fill = model))+
  theme_bw()

