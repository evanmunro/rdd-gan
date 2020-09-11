
#library(JuliaCall)

#julia_eval("include(\"../BayesRDD/src/code/estimators/bayesrdd.jl\")")
#julia_eval("include(\"../BayesRDD/src/code/estimators/gpsimple.jl\")")

global_M <- function(df) { 
  model <- stats::lm(df$y ~ 0 + outer(df$x, 0:4, "^"))
  r1 <- unname(model$coefficients)
  f2 <- function(x) abs(2*r1[3]+6*x*r1[4]+12*x^2*r1[5])
  return(f2(df$x))
}

select_M<- function(da, db, bw) {
  a2 = global_M(da) 
  b2 = global_M(db) 
  Ma = max(a2[abs(da$x)<bw])
  Mb = max(b2[abs(db$x)<bw])
  return(max(Ma, Mb))
}

rddBayesFast <- function(Y, X, c=0) {
  nfolds = ceiling(length(X)/1000)
  folds = sample(1:nfolds, length(X), replace=T)
  ates = sapply(1:nfolds, FUN =function(x) return(julia_eval("gpRDDSimple")(Y[folds==x], X[folds==x])))
  return(list(ate=mean(ates), se=1, ci.lower=1, ci.upper=1, bw=1))
}

rddKRRFast <- function(Y, X, c=0) {
  nfolds = ceiling(length(X)/250)
  folds = sample(1:nfolds, length(X), replace=T)
  print(nfolds)
  ates = sapply(1:nfolds, FUN = function(x) return(rddKRR(Y[folds==x], X[folds==x])$ate))
  print(ates)
  return(list(ate=mean(ates), se=0, ci.lower=0, ci.upper=0, bw=0))
}

rddBayes <- function(Y, X, c=0, discrete=TRUE) {
  #result <- julia_eval("gpRDDSimple")(Y, X)
  #start_time = Sys.time()
  result <- list(0)
  i=1
  while(length(result)<4) {
    if(i>1) print("julia error, trying again")
    result <- julia_eval("bayesRDD")(Y, X)
    i=i+1
  }
 # end_time = Sys.time()
  print("Bayes")
  #print(end_time - start_time)
  return(list(ate=result[[1]], se=result[[2]], ci.lower=result[[3]], ci.upper=result[[4]], bw=1))
}

rddKRR <- function(Y, X, c=0) {
  mu1  = listdtr::krr(X[X>c], Y[X>c])
  mu0 = listdtr::krr(X[X<c], Y[X<c])
  return(list(ate=predict(mu1, c(0))- predict(mu0, c(0)), se=0, ci.lower=0, ci.upper=0, bw=0))
}

rddIK <- function(Y, X, c=0) {
  rd <- rddtools::rdd_data(x=X, y=Y, cutpoint=0)
  bw <- rddtools::rdd_bw_ik(rd)
  #print(bw)
  model <- rdd::RDestimate(y~x, data.frame(y=Y, x=X), cutpoint=c, bw=bw)
  ate <- as.numeric(model$est[1])
  se <- model$se[1]
  ci <- model$ci[1,]
  return(list(ate=ate, se=se, ci.lower=ci[1], ci.upper=ci[2], bw=model$bw[1]))
}

#robust MSE optimal local linear regression
rddLLRM <- function(Y, X, c=0) {
  sink("/dev/null")
  model <- rdrobust::rdrobust(Y, X, c=c, bwselect="mserd")
  sink()
  ate <- model$coef[3]
  se <- model$se[3]
  ci.lower<- model$ci[3, 1]
  ci.upper <- model$ci[3, 2]
  return(list(ate=ate, se=se, ci.lower=ci.lower, ci.upper=ci.upper, bw=model$bws['h','left']))
}

#coverage optimal local linear regression
rddLLRC <- function(Y, X, c=0) {
  sink("/dev/null")
  model <- rdrobust::rdrobust(Y, X, c=c, bwselect="cerrd")
  sink()
  ate <- model$coef[3]
  se <- model$se[3]
  ci.lower<- model$ci[3, 1]
  ci.upper <- model$ci[3, 2]
  return(list(ate=ate, se=se, ci.lower=ci.lower, ci.upper=ci.upper, bw=model$bws['h','left']))
}

rddIW <- function(Y, X, M, c=0) {
  #start_time = Sys.time()
  W = as.numeric(X>c)
  model <-  optrdd::optrdd(X,Y,W,
              max.second.derivative=M,verbose=F,estimation.point=c,
              try.elnet.for.sigma.sq=T,optimizer='mosek')
  ate <- model$tau.hat
  se <- model$sampling.se
  ci <- model$tau.hat +c(-1,1)*model$tau.plusminus
  #end_time = Sys.time()
  #print("IW")
  #print(end_time - start_time)
  return(list(ate=ate,se=se,ci.lower=ci[1],ci.upper=ci[2],bw=0))
}

rddAK <- function(Y, X, M, c=0) {
  #start_time = Sys.time()
  model <- RDHonest::RDHonest(y~x,data=data.frame(y=Y,x=X),cutoff=c,M=M,kern="triangular",opt.criterion="MSE",sclass="T")
  ate <- model$estimate
  se <- as.numeric(model$sd)
  ci <- c(as.numeric(model$lower),as.numeric(model$upper))
  #print("AK")
 # end_time = Sys.time()
  #print(end_time - start_time)
  return(list(ate=ate, se=se, ci.lower=ci[1], ci.upper=ci[2], bw=model$hp))
}

#local quadratic with bias correction, MSE optimal bandwidth
rddQD<- function(Y, X, c=0) {
  sink("/dev/null")
  model <- rdrobust::rdrobust(Y, X, c=c, p=2)
  sink()
  ate <- model$coef[3]
  se <- model$se[3]
  ci.lower<- model$ci[3, 1]
  ci.upper <- model$ci[3, 2]
  return(list(ate=ate, se=se, ci.lower=ci.lower, ci.upper=ci.upper, bw=model$bws['h','left']))
}

estimate_sample <- function(estimators, dfa, dfb, na, nb, gt) {
  dfa <- dfa[sample(1:nrow(dfa), na), ]
  dfb <- dfb[sample(1:nrow(dfb), nb), ]
  df <- data.frame(rbind(dfa, dfb))
  #M <- RDHonest::NPR_MROT.fit(RDHonest::RDData(df[, c("y","x")], cutoff=0))/20
  rd <- rddtools::rdd_data(x=df$x, y=df$y, cutpoint=0)
  bw <- rddtools::rdd_bw_ik(rd)
  M <- select_M(dfa, dfb, bw)
  print(M)
  estimates <- sapply(estimators, FUN = function(x) return(run_estimate(x, df, M)))
  return(estimates)
}

save_generated_sample <- function(gen_path, real_path, name) {
    gen <- read_feather(gen_path)
    real <- read.csv(real_path)
    real.da <- real[real$x>0, ]
    real.db <- real[real$x<=0, ]
    ga <- as.data.frame(gen[gen$x>0,])
    gb <- as.data.frame(gen[gen$x<0,])
    na = nrow(real.da)
    nb = nrow(real.db)
    print(na)
    print(nb)
    dfa <- ga[sample(1:nrow(ga), na), ]
    dfb <- gb[sample(1:nrow(gb), nb), ]
    df <- data.frame(rbind(dfa, dfb))
    print(dim(df))
    write.csv(df, file=name, row.names=F)
}

run_estimate <- function(estimator, df, M) {
  if (estimator %in% c("rddAK", "rddIW")) {
    est <- do.call(estimator, list(df$y, df$x, M=M))
  } else {
    est <- do.call(estimator, list(df$y, df$x))
  }
  return(c(ate=est$ate, ci.l=est$ci.lower, ci.u=est$ci.upper))
}


make_table <- function(samples, gt) {
  d <- dim(samples)[2]
  estimates <- rowMeans(samples[1, 1:d, ])
  bias <- apply(samples[1, 1:d, ], MARGIN=1, FUN=function(x){mean(x-gt)})
  sd <- apply(samples[1, 1:d, ], MARGIN=1, FUN=sd)
  rmse <- apply(samples[1, 1:d, ], MARGIN=1,FUN=function(x){sqrt(mean((x-gt)^2))})
  covg <- apply(samples[2:3, 1:d,], MARGIN=2, FUN = function(x){mean(gt <= x[2,] & gt >= x[1,])})
  width <- apply(samples[2:3, 1:d, ], MARGIN=2, FUN= function(x){mean(x[2,]- x[1,])})
  table <- data.frame(estimate = estimates, rmse=rmse, bias=bias, sd=sd, coverage=covg, width=width)
  rownames(table) <- names(covg)
  return(table)
}
