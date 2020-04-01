source("estimators.R")
library(ggplot2)
library(future.apply)
library(kableExtra)
plan(multiprocess(workers=20))

dgp <- function() { 
  n = nrow(data)*10
  X = data$x[sample(1:nrow(data),n,replace=T)]
  y = sapply(X,FUN = function(x) {
    mu = 0.4905 + 0.5648*x+ 0.0562*x^2 -0.2450*x^3   
    sigma = 0.1375
    return(rnorm(1,mean=mu,sd=sigma))
  } )
  return(data.frame(x=X,y=y))
}

ate_sample <- function(M) {
  df <- dgp() 
  M= RDHonest::NPR_MROT.fit(RDHonest::RDData(df[,c("y","x")], cutoff=0))
  return(c(rdd_IK(df$y,df$x)$ate,rdd_AK(df$y,df$x,M=M)$ate,rdd_IW(df$y,df$x,M=M)$ate,rdd_QD(df$y,df$x)$ate,M))
}

make_table <- function(samples,gt=0.0) {
  #gt = rdd_IK(gen$y,gen$x)$ate
  bias <- apply(samples[1:4,],MARGIN=1,FUN=function(x){mean(x-gt)})
  sd <- apply(samples[1:4,],MARGIN=1,FUN=sd)
  rmse <- apply(samples[1:4,],MARGIN=1,FUN=function(x){sqrt(mean((x-gt)^2))})
  table <- data.frame(estimate = rowMeans(samples[1:4,]),rmse=rmse,bias=bias,sd=sd)
  rownames(table) <- c("rdd_IK","rdd_AK","rdd_IW","rdd_QD")
  return(table)
}

data <- read.csv("data/cleaned/lee.csv")

#tuning parameter based on original data
M1 = RDHonest::NPR_MROT.fit(RDHonest::RDData(data[,c("y","x")], cutoff=0))
#tuning parameter based on sample of generated data
M2 = RDHonest::NPR_MROT.fit(RDHonest::RDData(dgp()[,c("y","x")], cutoff=0))
#tuning parameter based on real max deriv
M3 = 1.58
#tuning parameter based on deriv at threshold 
M4 = 0.1124

samples_M1 <- future_replicate(1000,ate_sample(M1)) 
#samples_M2 <- future_replicate(1000,ate_sample(M2))
#samples_M3 <- future_replicate(1000,ate_sample(M3)) 
#samples_M4 <- future_replicate(1000,ate_sample(M4))


Ms = samples_M1[5,] 
print("M statistics:")
print(mean(Ms))
print(sd(Ms))

ggplot(data.frame(M=Ms), aes(x=M)) + 
  geom_histogram(binwidth=0.5)
ggsave("m_poly_hist.pdf")

result <- make_table(samples_M1)
print(kable(result,"latex",booktabs=T)) 


#result <- make_table(samples_M2)
#kable(result,"latex",booktabs=T)

#result <- make_table(samples_M3)
#kable(result,"latex",booktabs=T)
#
#result <- make_table(samples_M4)
#kable(result,"latex",booktabs=T)
  
#run a different of sample size equal to lee each time 