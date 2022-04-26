library(RDHonest)
setwd("~/Documents/GitHub/rdd-gan")
data = read.csv("data/cleaned/lee.csv")
M = NPR_MROT.fit(RDData(data[, c("y","x")], cutoff=0))
result = RDHonest(y~x, data=data, cutoff=0.0, M=M, kern="optimal", opt.criterion="MSE", sclass="T")

#just look at below cutoff 
result$lff$m

result$lff$p