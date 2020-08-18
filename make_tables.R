args = commandArgs(TRUE)
outname =  args[1]
estimators = c("rddIK", "rddLLRM", "rddLLRC", "rddIW", "rddAK", "rddQD", "rddBayes")

outpath  = paste0("output/", outname)
tablepath = paste0("tables/", outname, ".txt")

sink(tablepath)

sink()

file_list =
rbind()
