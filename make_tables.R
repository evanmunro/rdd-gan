args = commandArgs(TRUE)
outname =  args[1]
tablename = args[2]
estimators = c("rddIK", "rddLLRM", "rddLLRC", "rddIW", "rddAK", "rddQD", "rddBayes")

outpath  = paste0("output/", outname)
tablepath = paste0("tables/", tablename)

rbind()
