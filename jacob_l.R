source("estimators.R")

data <- read.csv("data/cleaned/jl_math.csv")

real.da <- data[data['x']>0,]
real.db <- data[data['x']<=0,]


Y= data$y
X = data$x

fit = rdd::RDestimate(y~x,gen[sample(1:nrow(gen),70000),],cutpoint=0)
M=RDHonest::NPR_MROT.fit(RDHonest::RDData(data[,c("y","x")], cutoff=0))



#W = as.numeric(data$x<=threshold)
#model <-  optrdd::optrdd(data$x,data$y,W,
 #                        max.second.derivative=1.74,verbose=F,estimation.point=0,
  #                       try.elnet.for.sigma.sq=T,optimizer='mosek')
#model <- RDHonest::RDHonest(y~x,data=data,cutoff=0,M=1.74,opt.criterion="MSE",kern="optimal",bw.equal=TRUE)

results = rdd_IK(Y,X)
results = rbind(results,rdd_IW(Y,X))
results = rbind(results,rdd_AK(Y,X))
results = rbind(results,rdd_QD(Y,X))


gen <- feather::read_feather("data/generated/jl_math_generated.feather")
print(head(gen))
#X = gen$x
#Y=gen$y
#obs = sample(1:nrow(gen),6000)
#rdd_QD(Y[obs],X[obs])

ate_sample <- function(dfa,dfb,na,nb) {
  dfsa <- dfa[sample(1:nrow(dfa),na),]
  dfsb <- dfb[sample(1:nrow(dfb),nb),]
  df <- data.frame(rbind(dfsa,dfsb))
  print(dim(df))
  #print(sum(is.na(df)))
  #return(c(rdd_IK(df$y,df$x)$ate,rdd_AK(df$y,df$x)$ate,rdd_IW(df$y,df$x)$ate,rdd_QD(df$y,df$x)$ate))
  return(c(rdd_IK(df$y,df$x)$ate,rdd_AK(df$y,df$x)$ate,rdd_QD(df$y,df$x)$ate))
}
bw_sample <- function(dfa,dfb,na,nb) {
  dfa <- dfa[sample(1:nrow(dfa),na),]
  dfb <- dfb[sample(1:nrow(dfb),nb),]
  df <- data.frame(rbind(dfa,dfb))
  print(sum(is.na(df)))
  return(c(rdd_IK(df$y,df$x)$bw,rdd_AK(df$y,df$x)$bw,rdd_QD(df$y,df$x)$bw))
}

library(future.apply)
library(kableExtra)
plan(multiprocess)


gen$x = round(gen$x,1)
samples <- future_replicate(100,ate_sample(gen[gen$x>0,],gen[gen$x<=0,],nrow(real.da),nrow(real.db)))

make_table <- function(samples,gt=-0.0892) {
  bias <- apply(samples,MARGIN=1,FUN=function(x){mean(x-gt)})
  sd <- apply(samples,MARGIN=1,FUN=sd)
  mse <- apply(samples,MARGIN=1,FUN=function(x){mean((x-gt)^2)})
  table <- data.frame(estimate = rowMeans(samples),mse=mse,bias=bias,sd=sd)
  #rownames(table) <- c("rdd_IK","rdd_AK","rdd_IW","rdd_QD")
  rownames(table) <- c("rdd_IK","rdd_AK","rdd_QD")
  return(table)
}

result <- make_table(samples)
kable(result,"latex",booktabs=T)


model.gen <- rdd::RDestimate(y~x,gen)

model.gen  <- RDHonest::RDHonest(y~x,data=gen,cutoff=0,M=1.74,opt.criterion="MSE")
#model.base <- RDHonest::RDHonest(y~x,data,cutoff=0,M=10,opt.criterion="MSE")

model.gen <- rdrobust::rdrobust(gen$y,gen$x,c=0,p=2)
model.base <-rdrobust::rdrobust(data$y,data$x,c=0,p=2)
test <- gen[gen$x>-0.01&gen$x<0.01,]

#-0.0898 is GT for the RDHonest 
