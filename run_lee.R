source("estimators.R")

library(future.apply)
library(kableExtra)
plan(multiprocess(workers=4))

gen <- feather::read_feather("data/generated/lee_cond_generated.feather")