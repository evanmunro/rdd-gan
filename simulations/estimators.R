library(JuliaCall)
julia_command('include("../MinMaxRDD/adaptiveminmax.jl")')

calc_max_d1 <- function(Y, X) {
  df  = data.frame(y=Y, x=X)
  model <- stats::lm(df$y ~ 0 + outer(df$x, 0:4, "^"))
  r1 <- unname(model$coefficients)
  f1 <- function(x) abs(r1[2] + 2*r1[3]*x + 3*r1[4]*x^2 + 4*r1[5]*x^3)
  return(max(f1(X)))
}

first_deriv_bound <- function(Y, X) {
  return(max(calc_max_d1(Y[X>0], X[X>0]), calc_max_d1(Y[X<0], X[X<0])))
}

second_deriv_bound <- function(Y, X) {
  df = data.frame(y = Y, x=X)
  M = RDHonest::NPR_MROT.fit(RDHonest::RDData(df[, c("y","x")], cutoff=0))
  return(M)
}

discreteBWSelect <- function(Y, X) {
  hchoice = sort(unique(abs(X)))[1:25]
  amse = rep(0, length(hchoice))
  rd <- rddtools::rdd_data(x=X, y=Y, cutpoint=0)
  amse <- sapply(hchoice, function(h) rddtools::ik_amse(rd, kernel="Triangular", bw=h) )
  hstar = hchoice[which.min(amse)]
  return(hstar)
}

MinMaxAdapt <- function(Y, X, c=0) {
  julia_assign("y", Y)
  julia_assign("x", X)
  julia_command("using Revise")
  julia_command("result = adaptrdd(x, y)")
  ate = julia_eval("result.τ")
  se = julia_eval("result.σ")
  ci.lower = julia_eval("result.ci[1]")
  ci.upper = julia_eval("result.ci[2]")
  bw = 0 
  return(list(ate=ate, se=se, ci.lower = ci.lower, ci.upper = ci.upper, bw= bw))
}
DiscreteIK <- function(Y, X, c=0) {
  bw <- discreteBWSelect(Y, X)
  model <- rdd::RDestimate(y~x, data.frame(y=Y, x=X), cutpoint=c, bw=c(bw, bw))
  ate <- as.numeric(model$est[1])
  se <- model$se[1]
  ci <- model$ci[1,]
  return(list(ate=ate, se=se, ci.lower=ci[1], ci.upper=ci[2], bw=model$bw[1]))
}

LLinearIK <- function(Y, X, c=0) {
  #Y = Y[abs(X) < median(abs(X))]
  #X = X[abs(X) < median(abs(X))]
  rd <- rddtools::rdd_data(x=X, y=Y, cutpoint=0)
  bw <- rddtools::rdd_bw_ik(rd, kernel="Triangular")
  #print(bw)
  model <- rdd::RDestimate(y~x, data.frame(y=Y, x=X), cutpoint=c, bw=c(bw, bw))
  ate <- as.numeric(model$est[1])
  se <- model$se[1]
  ci <- model$ci[1,]
  return(list(ate=ate, se=se, ci.lower=ci[1], ci.upper=ci[2], bw=model$bw[1]))
}

#robust MSE optimal local linear regression
LLinearCatt <- function(Y, X, c=0) {
  sink("/dev/null")
  model <- rdrobust::rdrobust(Y, X, c=c, bwselect="mserd")
  sink()
  ate <- model$coef[3]
  se <- model$se[3]
  ci.lower<- model$ci[3, 1]
  ci.upper <- model$ci[3, 2]
  return(list(ate=ate, se=se, ci.lower=ci.lower, ci.upper=ci.upper, bw=model$bws['h','left']))
}

#local quadratic with bias correction, MSE optimal bandwidth
LQuadCatt <- function(Y, X, c=0) {
  sink("/dev/null")
  model <- rdrobust::rdrobust(Y, X, c=c, p=2, bwselect="mserd")
  sink()
  ate <- model$coef[3]
  se <- model$se[3]
  ci.lower<- model$ci[3, 1]
  ci.upper <- model$ci[3, 2]
  return(list(ate=ate, se=se, ci.lower=ci.lower, ci.upper=ci.upper, bw=model$bws['h','left']))
}


MinMaxIW <- function(Y, X, c=0) {
  M = second_deriv_bound(Y, X)
  W = as.integer(X>c)
  model <-  optrdd::optrdd(as.numeric(X),as.numeric(Y),W,
              max.second.derivative=M,verbose=F,estimation.point=c,
              try.elnet.for.sigma.sq=F)#, optimizer='mosek')
  ate <- model$tau.hat
  se <- model$sampling.se
  ci <- model$tau.hat + c(-1,1)*model$tau.plusminus
  return(list(ate=ate,se=se,ci.lower=ci[1],ci.upper=ci[2],bw=0))
}

MinMaxAK <- function(Y, X, c=0) {
  M = second_deriv_bound(Y, X)
  model <- RDHonest::RDHonest(y~x,data=data.frame(y=Y,x=X),cutoff=c,M=M,kern="optimal",opt.criterion="MSE",sclass="T")
  ate <- model$estimate
  se <- as.numeric(model$sd)
  ci <- c(as.numeric(model$lower),as.numeric(model$upper))
  return(list(ate=ate, se=se, ci.lower=ci[1], ci.upper=ci[2], bw=model$hp))
}

MinMaxIW20 <- function(Y, X, c=0) {
  W = as.integer(X>c)
  M = second_deriv_bound(Y, X)/20
  model <-  optrdd::optrdd(as.numeric(X), as.numeric(Y), W,
              max.second.derivative=M, verbose=F, estimation.point=c,
              try.elnet.for.sigma.sq=F)#, optimizer='mosek')
  ate <- model$tau.hat
  se <- model$sampling.se
  ci <- model$tau.hat + c(-1,1)*model$tau.plusminus
  return(list(ate=ate, se=se, ci.lower=ci[1], ci.upper=ci[2], bw=0))
}

MinMaxAK20 <- function(Y, X, c=0) {
  M = second_deriv_bound(Y, X)/20
  model <- RDHonest::RDHonest(y~x,data=data.frame(y=Y,x=X), cutoff=c, M=M, opt.criterion="MSE")
  ate <- model$estimate
  se <- as.numeric(model$sd)
  ci <- c(as.numeric(model$lower),as.numeric(model$upper))
  return(list(ate=ate, se=se, ci.lower=ci[1], ci.upper=ci[2], bw=model$hp))
}

MinMaxMON <- function(Y, X) {
  M = first_deriv_bound(Y, X)
  W <- X > 0
  Yt <- Y[X>0]
  Yc <- Y[X<0]
  Xt <- X[X>0]
  Xc <- X[X<0]
  ci <- rdadapt::CI_minimax_RD(Yt, Yc, Xt, Xc, M, 1, t.dir = "right", alpha = 0.05)$ci
  return(list(ate=(ci[1] + ci[2])/2, se=0.0, ci.lower=ci[1], ci.upper=ci[2], bw=0.0))
}

estimate_sample <- function(estimators, dfa, dfb, na, nb, gt) {
  dfa <- dfa[sample(1:nrow(dfa), na), ]
  dfb <- dfb[sample(1:nrow(dfb), nb), ]
  df <- data.frame(rbind(dfa, dfb))
  #sub <- abs(df$x) <0.25
 # df <- df[sub, ]
  #rd <- rddtools::rdd_data(x=df$x, y=df$y, cutpoint=0)
  #bw <- rddtools::rdd_bw_ik(rd)
  M <- second_deriv_bound(df$y, df$x)
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
  est <- do.call(estimator, list(df$y, df$x))
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
