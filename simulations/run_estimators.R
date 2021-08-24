library(arrow)
source("simulations/estimators.R")

args = commandArgs(TRUE)
datapath = args[1]
outpath  = args[2]

print(datapath)
print(outpath)

estimators = c("rddIK", "rddCLL", "rddCLQ", "rddIW", "rddAK",  "rddGAM")
df <- read_feather(datapath)
print(length(unique(df$x)))
M <- RDHonest::NPR_MROT.fit(RDHonest::RDData(df[, c("y","x")], cutoff=0))
estimates <- t(sapply(estimators, FUN = function(x) return(run_estimate(x, df, M))))

write.table(estimates, file=outpath, row.names=FALSE, col.names=FALSE, sep=",")
