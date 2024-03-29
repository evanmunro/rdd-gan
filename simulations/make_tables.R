library(kableExtra)
library(arrow)

make_tables <- function(name, gt) {
    print(name)
    if (name %in% c("jl_math", "mats_math")) {
        estimators = c("LLinearIK", "DiscreteIK", "LLinearCatt", "LQuadCatt", "MinMaxIW", "MinMaxAK", "MinMaxIW20", "MinMaxAK20", "MinMaxMON")
    } else{
        estimators = c("LLinearIK", "LLinearCatt", "LQuadCatt", "MinMaxIW", "MinMaxAK", "MinMaxIW20", "MinMaxAK20", "MinMaxMON")
    }
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
    results <- data.frame(results)
    results$rmsep <- results$rmse/min(results$rmse)
    results <- results[, c("est", "rmsep", "rmse", "bias", "sd", "covg", "width")]
    result.ltx <- kable(results, "latex", digits=4, booktabs=T)
    print(result.ltx)

    sink(tablepath)
      print(result.ltx)
    sink()
}

make_tables("curved", 1.02711510658264)
make_tables("mc_curved", 0.5)
make_tables("lee", 0.08543315529823303)
make_tables("jl_math", -0.208665132522583)
make_tables("mats_math", -0.03234344720840454)
make_tables("meyersson", 3.809218406677246)
make_tables("senate", 7.691082000732422)
make_tables("brazil", -4.142387390136719)
