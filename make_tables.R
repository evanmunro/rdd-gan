library(kableExtra)
library(feather)
args = commandArgs(TRUE)
name =  args[1]

rddIK <- function(Y, X, c=0) {
  rd <- rddtools::rdd_data(x=X, y=Y, cutpoint=0)
  bw <- rddtools::rdd_bw_ik(rd)
  #print(bw)
  model <- rdd::RDestimate(y~x, data.frame(y=Y, x=X), cutpoint=c, bw=bw)
  ate <- as.numeric(model$est[1])
  se <- model$se[1]
  ci <- model$ci[1,]
  return(list(ate=ate, se=se, ci.lower=ci[1], ci.upper=ci[2], bw=model$bw[1]))
}


estimators = c("rddIK", "rddLLRM", "rddLLRC", "rddIW", "rddAK", "rddQD", "rddBayes")

outpath  = paste0("output/", name, "/")
tablepath = paste0("tables/", name, "_sims", ".txt")

files = paste0(outpath, list.files(outpath))
samples <- lapply(files,
            FUN = function(x) {
                df <- read.csv(x, header=F)
                df$type <- estimators
                return(df)
              })
samples <- do.call(rbind, samples)
colnames(samples) = c("ate", "ci.l", "ci.u", "type")
print("number of samples: ")
print(nrow(samples)/length(estimators))

df.g <- read_feather(paste0("data/generated/", name, "_generated.feather"))
print(gt)
gt <- rddIK(df.g$y,df.g$x)$ate

results <- t(sapply(estimators, FUN = function (x) {
                    df <- samples[samples$type==x, ]
                    est <- mean(df$ate)
                    rmse <- sqrt(mean((df$ate - gt)^2))
                    bias <- mean(df$ate-gt)
                    sd <- sd(df$ate)
                    covg <- mean(gt <= df$ci.u & gt >= df$ci.l)
                    width <- mean(df$ci.u - df$ci.l)
                    return(c(est, rmse, bias, sd, covg, width))
                }))
colnames(results) <- c("est", "rmse", "bias", "sd", "covg", "width")
result.ltx <- kable(results, "latex", digits=4, booktabs=T)
print(result.ltx)

sink(tablepath)
  print(result.ltx)
sink()
