setwd("/scratch/users/munro/rdd-gan")
source("estimators.R")
library(future.apply)
library(ggplot2)
library(feather) 
library(kableExtra)
plan(multiprocess(workers=2))

generate_tables <- function(real_path,gen_path,n.sims=3) {
  print(real_path)
  data <- read.csv(real_path)
  real.da <- data[data$x>0,]
  real.db <- data[data$x<=0,]
  #estimate on real data 
  M= RDHonest::NPR_MROT.fit(RDHonest::RDData(data[,c("y","x")], cutoff=0))
  estimators = c("RDD_IK","RDD_AK","RDD_IW","RDD_QD")
  estimates = rbind(rdd_IK(data$y,data$x), rdd_AK(data$y,data$x,M), 
                    rdd_IW(data$y,data$x,M), rdd_QD(data$y,data$x)) 
  colnames(estimates) = c("estimate","se","ci.lower","ci.upper","bw")
  rownames(estimates) = c("RDD_IK","RDD_AK","RDD_IW","RDD_QD")
  estimates = apply(data.frame(estimates),MARGIN=2,FUN= function(x) as.numeric(x))
  
  print(kable(data.frame(estimates),"latex",digits=3,booktabs=T))
  
  #then run simulation 
  gen <- read_feather(gen_path)
  gt = rdd_IK(gen$y,gen$x)$ate
  samples <- future_replicate(n.sims,ate_sample(gen[gen$x>0,],gen[gen$x<=0,],nrow(real.da),nrow(real.db),gt))
  result <- make_table(samples,gt) 
  print(kable(result,"latex",digits=4,booktabs=T))
}

generate_tables("data/cleaned/lee.csv", "data/generated/lee_generated.feather")

#generate_tables("data/cleaned/m_math.csv", "data/generated/mats_math_generated.feather")

#generate_tables("data/cleaned/jl_math.csv","data/generated/jl_math_generated.feather") 