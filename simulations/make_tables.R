library(kableExtra)
library(arrow)

make_tables <- function(name, gt) {
    print(name)
    estimators = c("rddIK", "rddCLL", "rddCLQ", "rddIW", "rddAK",  "rddGAM")
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
}

make_tables("meyersson", 3.809218406677246)
make_tables("senate", 7.691082000732422)
make_tables("brazil", )
#make_tables("lee", 0.08543315529823303)
#make_tables("jl_math", -0.208665132522583)
#make_tables("mats_math", -0.03234344720840454)
