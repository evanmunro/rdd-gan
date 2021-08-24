
#library(JuliaCall)
#julia_eval("include(\"simulations/ikhonest.jl\")")

#julia_eval("include(\"../BayesRDD/src/code/estimators/bayesrdd.jl\")")
#julia_eval("include(\"../BayesRDD/src/code/estimators/gpsimple.jl\")")
library(mgcv)
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

rddIKHonest <- function(Y, X, c=0) {
  result <- julia_eval("IKHonest")(as.numeric(Y), as.numeric(X))
  print(result)
  return(list(ate = result[[1]], se = result[[2]], ci.lower = result[[3]], ci.upper = result[[4]], bw=0))
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
  #Y = Y[abs(X) < median(abs(X))]
  #X = X[abs(X) < median(abs(X))]
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
rddCLL <- function(Y, X, c=0) {
  sink("/dev/null")
  model <- rdrobust::rdrobust(Y, X, c=c, bwselect="mserd")
  sink()
  ate <- model$coef[3]
  se <- model$se[3]
  ci.lower<- model$ci[3, 1]
  ci.upper <- model$ci[3, 2]
  return(list(ate=ate, se=se, ci.lower=ci.lower, ci.upper=ci.upper, bw=model$bws['h','left']))
}


#sp of 30 beats IK for lee data 
gamfit <- function(y, x, xpred) {
  #model <- gam(y ~ s(x, k=30, bs="re"), data = data.frame(y=y, x=x), method="ML")
  #weights = (max(abs(x)) - abs(x))/max(abs(x))^2
  #sp = 0.005473115
  model <- gam(y ~ s(x, k=30, bs="ts", sp = 30), data = data.frame(y=y, x=x), method="REML")
  print(model$sp)
  #model <- gam(y ~ s(x, k=30, bs="ts"), data = data.frame(y=y, x=x),weights=weights, select=TRUE, method="P-REML")
  #print(gam.check(model))
  return(predict(model, newdata=data.frame(x=xpred), se.fit=TRUE))
}

scamfit <- function(y, x, xpred) {
  model <- scam::scam(y ~ s(x, bs="mpi"),  data = data.frame(y=y, x=x))
  #model <- gam(y ~ s(x, k=30, bs="ts"),)
  #print(gam.check(model))
  return(predict(model, newdata=data.frame(x=xpred), se.fit=TRUE))
}

boostfit <- function(y, x, xpred){
  gam2 <- gamboost(y ~ bbs(x), data  = data.frame(y=y, x=x),
                   ctrl <- boost_control(mstop = 2))
  #cvm <- cvrisk(gam2)
  #model <- gam2[ mstop(cvm) ]
  model = gam2
  return(predict(model, newdata=data.frame(x=xpred)))
}

#alpha=0.8
polyfit <- function(y, x, xpred) {
  xAll = append(x, xpred)
  X = model.matrix(~poly(xAll, 5, raw=TRUE)-1)
  xpred = tail(X, 1)
  x = head(X, -1)
  cvfit <- glmnet::cv.glmnet(x, y, alpha=0.8)
  return(predict(cvfit, newx = xpred, s = "lambda.min"))
}

npr_fit <- function(y, x, xfit) {
  bw = thumbBw(x, y, 1, EpaK, weig = rep(1, length(y)))
  print(bw)
  model = locpol(y~x,bw=2*bw, xeval = c(0), data=data.frame(y=y,x=x))
  #model = npreg(tydat=y, txdat=data.frame(x=x), exdat = data.frame(x=xfit), regtype="ll", bwmethod="cv.aic", gradients=TRUE)
 # print(model$bw)
  return(model$lpFit$y)
}
rddNPR <- function (Y, X, c=0) {
  Y = Y[abs(X) < median(abs(X))]
  X = X[abs(X) < median(abs(X))]
  mu1 = npr_fit(Y[X>0], X[X>0], c(0))
  mu0 = npr_fit(Y[X<0], X[X<0], c(0))
  return(list(ate=mu1-mu0, se=0, ci.lower = 0, ci.upper = 0, bw =0 ))
}
rddBOOST <- function(Y, X, c=0){
  Y = Y[abs(X) < median(abs(X))]
  X = X[abs(X) < median(abs(X))]
  mu1 = boostfit(Y[X>0], X[X>0], c(0))
  mu0 = boostfit(Y[X<0], X[X<0], c(0))
  return(list(ate=mu1-mu0, se=0, ci.lower = 0, ci.upper = 0, bw =0 ))
}

rffit <- function(y, x, xpred) {
  model <- ranger::ranger(y~x, data=data.frame(y=y, x=x))
  return(predict(model, data.frame(x=xpred))$predictions)
}

#
rddRF <- function(Y, X, c=0) {
  Y = Y[abs(X) < median(abs(X))]
  X = X[abs(X) < median(abs(X))]
  mu1 = rffit(Y[X>0], X[X>0], c(0))
  mu0 = rffit(Y[X<0], X[X<0], c(0))
  return(list(ate=mu1-mu0, se=0, ci.lower = 0, ci.upper = 0, bw =0 ))
}

#does not work particularly well (larger sd ) about 0.0096 RMSE
rddPoly <- function(Y, X, c=0) {
  Y = Y[abs(X) < median(abs(X))]
  X = X[abs(X) < median(abs(X))]
  mu1 = polyfit(Y[X>0], X[X>0], c(0))
  mu0 = polyfit(Y[X<0], X[X<0], c(0))
  return(list(ate=mu1-mu0, se=0, ci.lower = 0, ci.upper = 0, bw =0 ))
}

rddSCAM <- function(Y, X, c=0) {
  #Y = Y[abs(X) < median(abs(X))]
  #X = X[abs(X) < median(abs(X))]
  mu1 = scamfit(Y[X>0], X[X>0], c(0))
  mu0 = scamfit(Y[X<0], X[X<0], c(0))
  se = sqrt(mu1$se.fit^2 + mu0$se.fit^2)
  tau = mu1$fit - mu0$fit
  return(list(ate=tau, se=se, ci.lower = tau-se*1.96, ci.upper = tau+se*1.96, bw =0 ))
}

rddGAM <- function(Y, X, c=0) {
  #Y = Y[abs(X) < median(abs(X))]
  #X = X[abs(X) < median(abs(X))]
  mu1 = gamfit(Y[X>0], X[X>0], c(0))
  mu0 = gamfit(Y[X<0], X[X<0], c(0))
  se = sqrt(mu1$se.fit^2 + mu0$se.fit^2)
  tau = mu1$fit - mu0$fit
  return(list(ate=tau, se=se, ci.lower = tau-se*1.96, ci.upper = tau+se*1.96, bw =0 ))
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
  W = as.integer(X>c)
  print(head(W))
  model <-  optrdd::optrdd(as.numeric(X),as.numeric(Y),W,
              max.second.derivative=M,verbose=F,estimation.point=c,
              try.elnet.for.sigma.sq=F)#, optimizer='mosek')
  ate <- model$tau.hat
  se <- model$sampling.se
  ci <- model$tau.hat + c(-1,1)*model$tau.plusminus
  #end_time = Sys.time()
  #print("IW")
  #print(end_time - start_time)
  return(list(ate=ate,se=se,ci.lower=ci[1],ci.upper=ci[2],bw=0))
}

rddAK <- function(Y, X, M, c=0) {
  #start_time = Sys.time()
  model <- RDHonest::RDHonest(y~x,data=data.frame(y=Y,x=X),cutoff=c,M=M,kern="optimal",opt.criterion="MSE",sclass="T")
  ate <- model$estimate
  se <- as.numeric(model$sd)
  ci <- c(as.numeric(model$lower),as.numeric(model$upper))
  #print("AK")
 # end_time = Sys.time()
  #print(end_time - start_time)
  return(list(ate=ate, se=se, ci.lower=ci[1], ci.upper=ci[2], bw=model$hp))
}


nn = function(Y, X, X_eval, probs=F, hidden=c(32, 8),
              l2=0.01, n_batch=128, n_epochs=50, lr=1e-3, ...) {
  if (probs) Y = as.character(Y)
  net = ANN2::neuralnetwork(X, Y, hidden.layers = hidden,
                           regression = !probs,
                           standardize = TRUE,
                           loss.type = ifelse(probs, "log", "squared"),
                           activ.functions = "relu",
                           optim.type = "adam",
                           L2 = l2,
                           batch.size = n_batch,
                           n.epochs = n_epochs,
                           verbose = F,
                           learn.rates = lr, ...)
  if (probs) return(predict(net, X_eval)$probabilities[, "class_1"])
  if (!probs) return(predict(net, X_eval)$predictions)
}

rddNN <- function(Y, X, M) {
  mu1 = nn(Y[X>0], X[X>0], c(0))
  mu0 = nn(Y[X<0], X[X<0], c(0))
  return(list(ate = mu1-mu0, se=0, ci.lower = 0, ci.upper = 0, bw =0 ))
}



#local quadratic with bias correction, MSE optimal bandwidth
rddCLQ<- function(Y, X, c=0) {
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
