library(feather)
source("estimators.R")

args = commandArgs(TRUE)
datapath = args[1]
outpath  = args[2]

print(datapath)
print(outpath)

run_estimate <- function(estimator, df, M) {
  if (estimator %in% c("rddAK", "rddIW")) {
    est <- do.call(estimator, list(df$y, df$x, M=M))
  } else {
    est <- do.call(estimator, list(df$y, df$x))
  }
  return(c(ate=est$ate, ci.l=est$ci.lower, ci.u=est$ci.upper))
}

estimators = c("rddIK", "rddLLRM", "rddLLRC", "rddIW", "rddAK", "rddQD", "rddBayes")
df <- read.csv(datapath)
M <- RDHonest::NPR_MROT.fit(RDHonest::RDData(df[, c("y","x")], cutoff=0))
estimates <- t(sapply(estimators, FUN = function(x) return(run_estimate(x, df, M))))

write.table(estimates, file=outpath, row.names=FALSE, col.names=FALSE, sep=",")
