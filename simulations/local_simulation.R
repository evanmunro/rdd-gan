startbig = Sys.time() 
setwd("~/Documents/rdd-gan/")
source("simulations/estimators.R")
library(future.apply)
library(ggplot2)
library(feather) 
library(kableExtra)
#plan(multiprocess(workers=2))


table_real_estimates <- function(estimators, data)  {
  M= RDHonest::NPR_MROT.fit(RDHonest::RDData(data[,c("y","x")], cutoff=0))
  estimates = rbind(rddIK(data$y,data$x), rddLLRM(data$y, data$x), 
                    rddLLRC(data$y, data$x), rddAK(data$y,data$x, M), 
                    rddIW(data$y,data$x, M), rddQD(data$y,data$x), rddBayes(data$y, data$x)) 
  colnames(estimates) = c("estimate", "se", "ci.lower", "ci.upper", "bw")
  rownames(estimates) = estimators 
  estimates = apply(data.frame(estimates), MARGIN=2, FUN= function(x) as.numeric(x))
  print(kable(data.frame(estimates),"latex", row.names=T, digits=3, booktabs=T))
}

generate_tables <- function(real_path, gen_path, n.sims=5, digits=NULL, small=NULL) {
  #estimators = c("rddIK", "rddLLRM", "rddLLRC", "rddIW", "rddAK", "rddQD")
  estimators = c("rddIK",  "rddAK")
  data <- read.csv(real_path)
  
  if(!is.null(small)) {
    data$dist = abs(data$x) 
    data = data[order(data$dist), ]
    data = data[1:small, ]
    # estimators = c(estimators, "rddKRR")
  }
  #data <- data[sample(1:nrow(data), 6000, replace=FALSE), ]
  real.da <- data[data$x>0, ]
  real.db <- data[data$x<=0, ]
  
  #data$x = data$x/mean(data$x)
  #estimate on real data 
  #table_real_estimates(estimators, data) 
  
  #then run simulation 
  gen <- read_feather(gen_path)
  gt = rddIK(gen$y,gen$x)$ate
  if(!is.null(digits)){ gen$x <- round(gen$x, digits) } 
  #gen$x = gen$x/mean(gen$x)
  samples <- replicate(n.sims, estimate_sample(estimators, as.data.frame(gen[gen$x>0,]), as.data.frame(gen[gen$x<0,]), floor(nrow(real.da)), floor(nrow(real.db))))
  save(samples,file='samples.RData')
  result <- make_table(samples, gt) 
  print(kable(result, "latex", digits=4, booktabs=T))
}
generate_tables("data/cleaned/lee.csv", "data/generated/lee_generated.feather")
#generate_tables("data/cleaned/m_math.csv", "data/generated/mats_math_generated.feather", digits=0)
#generate_tables("data/cleaned/jl_math.csv","data/generated/jl_math_generated.feather", digits=2) 
 
#