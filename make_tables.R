library(kableExtra)
args = commandArgs(TRUE)
name =  args[1]
estimators = c("rddIK", "rddLLRM", "rddLLRC", "rddIW", "rddAK", "rddQD", "rddBayes")

outpath  = paste0("output/", outname)
tablepath = paste0("tables/", outname, ".txt")

files = list.files(outpath)
print(files)
#samples = rbind(lapply(files, FUN = function(x) return read.csv(x, row.names=estimators))
#colnames(samples) = c("ate", "ci.l", "ci.u")

#collate_samples  <- function(est.name, samples) {
  #

#}
#sink(tablepath)

#sink()

#file_list =
rbind()
