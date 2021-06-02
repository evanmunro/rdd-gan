library(kableExtra)
library(arrow)
include("simulations/estimators.R")

sample_generate <- function(data.real, data.gen, name) {
  da.r = data.real[data.real$x >0, ]
  db.r = data.real[data.real$x <=0, ]
  na = nrow(da.r)
  nb = nrow(db.r)
  dfa = data.gen[data.gen$x >0, ]
  dfb = data.gen[data.gen$x <=0, ]
  dfa <- dfa[sample(1:nrow(dfa), na), ]
  dfb <- dfb[sample(1:nrow(dfb), nb), ]
  df <- data.frame(rbind(dfa, dfb))
  if (name == "jl_math") {
    df$x = round(df$x, 2)
  }
  if (name == "mats_math") {
    df$x = round(df$x, 0)
  }
  return(df)
}

summarize <- function(variable) {
  avg = round(mean(variable), 2)
  std = round(sd(variable), 2)
  return(paste(avg, " (", std, ")", sep=""))
}

te <- function(data) {
  result = rddIK(data$y, data$x)
  return(paste(round(result$ate, 5), " (", round(result$se,5), ")", sep=""))
}

table_stats = data.frame()
table_te = data.frame()
for (name in c("lee", "jl_math", "mats_math")) {
  data.real  = read.csv(paste("data/cleaned/", name, ".csv", sep=""))
  print(nrow(data.real))
  data.gen = read_feather(paste("data/generated/", name, "_generated.feather", sep=""))
  data.gen = sample_generate(data.real, data.gen, name)
  labelx = paste(name, "running variable")
  labely = paste(name, "outcome")
  labelte = paste(name, "treatment effect")
  table_stats = rbind(table_stats, c(labelx, summarize(data.real$x), summarize(data.gen$x)))
  table_stats = rbind(table_stats, c(labely, summarize(data.real$y), summarize(data.gen$y)))

  table_te = rbind(table_te, c(labelte, te(data.real), te(data.gen)))
}

names(table_te) = c(" ", "Real", "Generated")
names(table_stats) = c(" ", "Real", "Generated")

print(kable(table_te, "latex", booktabs=T))
print(kable(table_stats, "latex", booktabs=T))
