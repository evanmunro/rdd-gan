library(feather)

args = commandArgs(TRUE)
name = "jl_math"
fraction = 1
runs = 2

function submit_slurm(name, chunkpath, runpath) {

   sink(paste(".job/", name,'.job',sep=""))

   # the basic job submission script is a bash script
   cat("#!/bin/bash\n")
   cat("#SBATCH --job-name=", dataname, ".job\n", sep="")
   cat("#SBATCH --output=.out/", dataname, ".out\n", sep="")
   cat("#SBATCH --error=.out/", dataname, ".err\n", sep="")
   cat("#SBATCH --time=00:05:00\n")
   cat("#SBATCH --mem=8GB\n")
   cat("#SBATCH -p owners\n")
   cat("Rscript run_estimators.R", chunkpath, runpath, sep=" ")

   # Close the sink!
   sink()

   # Submit to run on cluster
   system(paste("sbatch",paste(".job/", name, ".job",sep="")))

}

cleanpath = paste0("data/cleaned/", name, ".csv")
genpath = paste0("data/generated/", name, ".feather")

dfr= read.csv(cleanpath)
na = floor(length(dfr$x[dfr$x>0])*fraction)
nb = floor(length(dfr$x[dfr$x<=0])*fraction)
dfg= read_feather(genpath)
df.ga = dfg[dfg$x>0, ]
df.gb = dfg[dfg$x<=0, ]

for (i in 1:runs) {
  dataname = paste0(name, "_", sample(1:100000))
  chunkpath = paste0("data/chunks", dataname, ".feather")
  outpath = paste0("output/", name, "/", dataname, "run.csv\n")

  dfa <- df.ga[sample(1:nrow(df.ga), na), ]
  dfb <- df.gb[sample(1:nrow(df.gb), nb), ]
  df <- data.frame(rbind(dfa, dfb))

  write_feather(df, chunkpath)
  submit_slurm(name, chunkpath, runpath)
}
