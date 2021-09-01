library(arrow)
source("simulations/estimators.R")

args = commandArgs(TRUE)
name = args[1]
datapath = args[2]
outpath  = args[3]

print(datapath)
print(outpath)
if (name %in% c("jl_math", "mats_math")) {
    estimators = c("LLinearIK", "DiscreteIK", "LLinearCatt", "LQuadCatt", "MinMaxIW", "MinMaxAK", "MinMaxIW20", "MinMaxAK20", "MinMaxMON")
} else{
    estimators = c("LLinearIK", "LLinearCatt", "LQuadCatt", "MinMaxIW", "MinMaxAK", "MinMaxIW20", "MinMaxAK20", "MinMaxMON")
}
df <- read_feather(datapath)
print(length(unique(df$x)))
M <- RDHonest::NPR_MROT.fit(RDHonest::RDData(df[, c("y","x")], cutoff=0))
estimates <- t(sapply(estimators, FUN = function(x) return(run_estimate(x, df, M))))

write.table(estimates, file=outpath, row.names=FALSE, col.names=FALSE, sep=",")
