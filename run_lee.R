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
