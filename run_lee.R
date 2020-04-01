source("estimators.R")
library(future.apply)
library(ggplot2)
library(kableExtra)
plan(multiprocess(workers=4))

data = read.csv("data/cleaned/lee.csv")
gen <- feather::read_feather("data/generated/lee_gen.feather")

real.da <- data[data['x']>0,]
real.db <- data[data['x']<=0,]

ate_sample <- function(dfa,dfb,na,nb) {
  dfa <- dfa[sample(1:nrow(dfa),na),]
  dfb <- dfb[sample(1:nrow(dfb),nb),] 
  df <- data.frame(rbind(dfa,dfb))
  M= RDHonest::NPR_MROT.fit(RDHonest::RDData(df[,c("y","x")], cutoff=0))
  return(c(rdd_IK(df$y,df$x)$ate,rdd_AK(df$y,df$x,M)$ate,rdd_IW(df$y,df$x,M)$ate,rdd_QD(df$y,df$x)$ate,M))
}


make_table <- function(samples) {
  gt = rdd_IK(gen$y,gen$x)$ate 
  bias <- apply(samples[1:4,],MARGIN=1,FUN=function(x){mean(x-gt)})
  sd <- apply(samples[1:4,],MARGIN=1,FUN=sd)
  rmse <- apply(samples[1:4,],MARGIN=1,FUN=function(x){mean((x-gt)^2)})
  table <- data.frame(estimate = rowMeans(samples[1:4,]),rmse=rmse,bias=bias,sd=sd)
  rownames(table) <- c("rdd_IK","rdd_AK","rdd_IW","rdd_QD")
  return(table) 
}

samples <- future_replicate(20,ate_sample(gen[gen$x>0,],gen[gen$x<=0,],nrow(real.da),nrow(real.db)))

result <- make_table(samples) 
print(kable(result,"latex",booktabs=T))

Ms = samples[5,] 

ggplot(data.frame(M=Ms), aes(x=M)) + 
  geom_histogram(binwidth=0.5)
ggsave("m_hist.pdf")
