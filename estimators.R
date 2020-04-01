

rdd_IK <- function(Y,X, c=0) {
  rd<- rddtools::rdd_data(x=X, y=Y, cutpoint=0)
  bw <- rddtools::rdd_bw_ik(rd)
  #print(bw)
  model <- rdd::RDestimate(y~x,data.frame(y=Y,x=X),cutpoint=c,bw=bw)
  ate <- as.numeric(model$est[1])
  se <- model$se[1]
  ci <- model$ci[1,]
  return(list(ate=ate,se=se,ci.lower=ci[1],ci.upper=ci[2],bw=model$bw[1]))
}

rdd_IW <- function(Y,X,M,c=0) {
  W = as.numeric(X>c)
  model <-  optrdd::optrdd(X,Y,W,
              max.second.derivative=M,verbose=F,estimation.point=c,
              try.elnet.for.sigma.sq=T,optimizer='mosek')
  ate <- model$tau.hat
  se <- model$sampling.se
  ci <- model$tau.hat +c(-1,1)*model$tau.plusminus
  return(list(ate=ate,se=se,ci.lower=ci[1],ci.upper=ci[2],bw=0))
}

rdd_AK <- function(Y,X,M,c=0) {
  model <- RDHonest::RDHonest(y~x,data=data.frame(y=Y,x=X),cutoff=c,M=M,kern="triangular",opt.criterion="MSE",sclass="T")
  ate <- model$estimate
  se <- as.numeric(model$sd)
  ci <- c(as.numeric(model$lower),as.numeric(model$upper))
  return(list(ate=ate,se=se,ci.lower=ci[1],ci.upper=ci[2],bw=model$hp))
}

rdd_QD<- function(Y,X,c=0) {
  model <- rdrobust::rdrobust(Y,X,c=c,p=2)
  ate <- model$coef[3]
  se <- model$se[3]
  ci <- ate+c(-1,1)*1.96*se
  return(list(ate=ate,se=se,ci.lower=ci[1],ci.upper=ci[2],bw=model$bws['h','left']))
}


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
  rmse <- apply(samples[1:4,],MARGIN=1,FUN=function(x){sqrt(mean((x-gt)^2))})
  table <- data.frame(estimate = rowMeans(samples[1:4,]),rmse=rmse,bias=bias,sd=sd)
  rownames(table) <- c("rdd_IK","rdd_AK","rdd_IW","rdd_QD")
  return(table)
}
