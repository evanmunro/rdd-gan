source("estimators.R")

library(future.apply)
library(kableExtra)
plan(multiprocess)

dgp <- function() { 
  n = nrow(data)
  X = data$x[sample(1:n,n,replace=F)]
  y = sapply(X,FUN = function(x) {
    mu = 0.4905 + 0.5648*x+ 0.0562*x^2 -0.2450*x^3   
    sigma = 0.1375
    return(rnorm(1,mean=mu,sd=sigma))
  } )
  return(data.frame(x=X,y=y))
}

ate_sample <- function(M) {
  df <- dgp() 
  print(sum(is.na(df)))
  return(c(rdd_IK(df$y,df$x)$ate,rdd_AK(df$y,df$x,M=M)$ate,rdd_IW(df$y,df$x,M=M)$ate,rdd_QD(df$y,df$x)$ate))
}

make_table <- function(samples,gt=0.0) {
  bias <- apply(samples,MARGIN=1,FUN=function(x){mean(x-gt)})
  sd <- apply(samples,MARGIN=1,FUN=sd)
  rmse <- apply(samples,MARGIN=1,FUN=function(x){sqrt(mean((x-gt)^2))})
  table <- data.frame(estimate = rowMeans(samples),rmse=rmse,bias=bias,sd=sd)
  rownames(table) <- c("rdd_IK","rdd_AK","rdd_IW","rdd_QD")
  return(table)
}

data <- read.csv("data/cleaned/lee.csv")

#tuning parameter based on original data
M1 = RDHonest::NPR_MROT.fit(RDHonest::RDData(data[,c("y","x")], cutoff=0))
#tuning parameter based on sample of generated data
M2 = RDHonest::NPR_MROT.fit(RDHonest::RDData(dgp(), cutoff=0))
#tuning parameter based on real max deriv
M3 = 1.58
#tuning parameter based on deriv at threshold 
M4 = 0.1124

samples_M1 <- future_replicate(1000,ate_sample(M1)) 
samples_M2 <- future_replicate(1000,ate_sample(M2))
samples_M3 <- future_replicate(1000,ate_sample(M3)) 
samples_M4 <- future_replicate(1000,ate_sample(M4))

result <- make_table(samples_M1)
kable(result,"latex",booktabs=T)


result <- make_table(samples_M2)
kable(result,"latex",booktabs=T)


result <- make_table(samples_M3)
kable(result,"latex",booktabs=T)

result <- make_table(samples_M4)
kable(result,"latex",booktabs=T)
  
#run a different of sample size equal to lee each time 