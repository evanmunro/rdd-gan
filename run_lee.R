setwd("/scratch/users/munro/rdd-gan")
source("estimators.R")
library(future.apply)
library(ggplot2)
library(kableExtra)
plan(multiprocess(workers=4))

data = read.csv("data/cleaned/lee.csv")
gen <- feather::read_feather("data/generated/lee_gen.feather")

real.da <- data[data['x']>0,]
real.db <- data[data['x']<=0,]

samples <- future_replicate(20,ate_sample(gen[gen$x>0,],gen[gen$x<=0,],nrow(real.da),nrow(real.db)))

result <- make_table(samples) 
print(kable(result,"latex",booktabs=T))

Ms = samples[5,] 

print("M statistics:")
print(mean(Ms))
print(sd(Ms))

ggplot(data.frame(M=Ms), aes(x=M)) + 
  geom_histogram(binwidth=0.5)
ggsave("m_hist.pdf")

#M ocnvergence
data = read.csv("data/cleaned/lee.csv")
dgp <- function() { 
  n = 10000000
  X = data$x[sample(1:nrow(data),n,replace=T)]
  y = sapply(X,FUN = function(x) {
    mu = 0.4905 + 0.5648*x+ 0.0562*x^2 -0.2450*x^3   
    sigma = 0.1375
    return(rnorm(1,mean=mu,sd=sigma))
  } )
  return(data.frame(x=X,y=y))
}


#
big_poly <- dgp()

m_seq = c(200,400,800,1600,3200,6400,12800,25600,51200,102400,204800,409600,819200,1638400,3276800,6553600,10000000)
convg = sapply(m_seq, FUN= function (x) {
                      rows = sample(1:nrow(big_poly),x,replace=F)
                      M= RDHonest::NPR_MROT.fit(RDHonest::RDData(big_poly[rows,c("y","x")], cutoff=0))
                      return(M) 
              })
