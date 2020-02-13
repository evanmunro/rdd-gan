source("estimators.R")

data <- read.csv("data/cleaned/lee.csv")
X = data$x
Y = data$y 
M=RDHonest::NPR_MROT.fit(RDHonest::RDData(data[,c("y","x")], cutoff=0))

fit = rdd::RDestimate(y~x,data,model=TRUE,frame=TRUE)
bw= fit$bw[[1]]
wts.ik = rdd::kernelwts(data$x,0,bw)
#plot(data$x[wts>0],wts[wts>0],type="p")

W = as.numeric(data$x>0)
model <-  optrdd::optrdd(data$x,data$y,W,
                         max.second.derivative=14.28,verbose=F,estimation.point=0,
                         try.elnet.for.sigma.sq=T,optimizer='mosek')
wts.iw <- abs(model$gamma)/sum(abs(model$gamma))
#plot(data$x[wts>0.000001],(wts/sum(wts))[wts>0.000001])

dat <- data.frame(x = data$x, wts.ik=wts.ik,wts.iw=wts.iw)
dat.m <- reshape2::melt(dat, id.vars = "x")
ggplot(dat.m, aes(x, value, colour = variable)) +
  geom_point(alpha=0.2)# +
  #scale_colour_manual(values = c("red", "blue"))

model <- RDHonest::RDHonest(y~x,data=data,cutoff=0,M=14.28,opt.criterion="MSE",kern="optimal",bw.equal=TRUE)

#results = rdd_IK(Y,X)
##results = rbind(results,rdd_IW(Y,X))
#results = rbind(results,rdd_AK(Y,X))
#results = rbind(results,rdd_QD(Y,X)) 


gen <- feather::read_feather("data/generated/lee_cond_generated.feather")
print(head(gen))
#X = gen$x
#Y=gen$y
#obs = sample(1:nrow(gen),6000)
#rdd_QD(Y[obs],X[obs])

ate_sample <- function(dfa,dfb,na=3818,nb=2740) {
  dfa <- dfa[sample(1:nrow(dfa),na),]
  dfb <- dfb[sample(1:nrow(dfb),nb),] 
  df <- data.frame(rbind(dfa,dfb))
  print(sum(is.na(df)))
  return(c(rdd_IK(df$y,df$x)$ate,rdd_AK(df$y,df$x)$ate,rdd_IW(df$y,df$x)$ate,rdd_QD(df$y,df$x)$ate))
  #return(c(rdd_IK(df$y,df$x)$ate,rdd_AK(df$y,df$x)$ate,rdd_QD(df$y,df$x)$ate))
}
bw_sample <- function(dfa,dfb,na=3818,nb=2740) {
  dfa <- dfa[sample(1:nrow(dfa),na),]
  dfb <- dfb[sample(1:nrow(dfb),nb),] 
  df <- data.frame(rbind(dfa,dfb))
  print(sum(is.na(df)))
  return(c(rdd_IK(df$y,df$x)$bw,rdd_AK(df$y,df$x)$bw,rdd_QD(df$y,df$x)$bw))
}

library(future.apply)
library(kableExtra)
plan(multicore,workers=1)
dfa<- gen[gen$x>0,]
dfb<- gen[gen$x<=0,]

samples <- future_replicate(100,ate_sample(dfa,dfb)) 

make_table <- function(samples,gt=0.0721) {
  bias <- apply(samples,MARGIN=1,FUN=function(x){mean(x-gt)})
  sd <- apply(samples,MARGIN=1,FUN=sd)
  mse <- apply(samples,MARGIN=1,FUN=function(x){mean((x-gt)^2)})
  table <- data.frame(estimate = rowMeans(samples),mse=mse,bias=bias,sd=sd)
  rownames(table) <- c("rdd_IK","rdd_AK","rdd_IW","rdd_QD")
  return(table)
}

result <- make_table(samples)
kable(result,"latex",booktabs=T)

#0.07217 
#0.0711




rowMeans(future_replicate(100,bw_sample(dfa,dfb)))

model.gen <- rdd::RDestimate(y~x,gen)
model.base <- rdd::RDestimate(y~x,data)

dl <- RDHonest::RDData(data,cutoff=0)

model.gen  <- RDHonest::RDHonest(y~x,data=gen,cutoff=0,M=14.28,opt.criterion="MSE")
model.base <- RDHonest::RDHonest(y~x,data,cutoff=0,M=10,opt.criterion="MSE")

model.gen <- rdrobust::rdrobust(gen$y,gen$x,c=0,p=2)
model.base <-rdrobust::rdrobust(data$y,data$x,c=0,p=2)
test <- gen[gen$x>-0.01&gen$x<0.01,]