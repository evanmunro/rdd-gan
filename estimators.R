

rdd_IK <- function(Y,X, c=0) {
  model <- rdd::RDestimate(y~x,data.frame(y=Y,x=X),cutpoint=c)
  ate <- as.numeric(model$est[1])
  se <- model$se[1]
  ci <- model$ci[1,]
  return(list(ate=ate,se=se,ci.lower=ci[1],ci.upper=ci[2]))
}

rdd_IW <- function(Y,X,c=0) {
  model <-  optrdd::optrdd(X,Y,X>c,
              max.second.derivative=10,verbose=F,
              try.elnet.for.sigma.sq=T)
  ate <- model$tau.hat
  se <- model$sampling.se
  ci <- model$tau.hat +c(-1,1)*model$tau.plusminus
  return(list(ate=ate,se=se,ci.lower=ci[1],ci.upper=ci[2]))
}

rdd_AK <- function(Y,X,c=0) {
  model <- RDHonest::RDHonest(y~x,data=data.frame(y=Y,x=X),cutoff=c,M=10,opt.criterion="MSE")
  ate <- model$estimate
  se <- as.numeric(model$sd)
  ci <- c(as.numeric(model$lower),as.numeric(model$upper))
  return(list(ate=ate,se=se,ci.lower=ci[1],ci.upper=ci[2]))
}

rdd_QD<- function(Y,X,c=0) {
  model <- rdrobust::rdrobust(Y,X,c=c,p=2)
  ate <- model$coef[3]
  se <- model$se[3]
  ci <- ate+c(-1,1)*1.96*se
  return(list(ate=ate,se=se,ci.lower=ci[1],ci.upper=ci[2]))
}
