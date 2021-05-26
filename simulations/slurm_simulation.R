library(arrow)

args = commandArgs(TRUE)
name = "lee"
fraction = 1
runs = 2

submit_slurm <- function(name, chunkpath, runpath) {

   sink(paste(".job/", name,'.job',sep=""))

   # the basic job submission script is a bash script
   cat("#!/bin/bash\n")
   cat("#SBATCH --job-name=", dataname, ".job\n", sep="")
   cat("#SBATCH --output=.out/", dataname, ".out\n", sep="")
   cat("#SBATCH --error=.out/", dataname, ".err\n", sep="")
   cat("#SBATCH --time=00:5:00\n")
   cat("#SBATCH --mem=6GB\n")
   cat("#SBATCH -p owners\n")
   cat("module load gcc/10.1.0\n")
   cat("module load R/4.0.2\n")
   cat("Rscript simulations/run_estimators.R", chunkpath, runpath, sep=" ")

   # Close the sink!
   sink()
   # Submit to run on cluster
   system(paste("sbatch",paste(".job/", name, ".job",sep="")))

}

cleanpath = paste0("data/cleaned/", name, ".csv")
genpath = paste0("data/generated/", name, "_generated", ".feather")

dfr= read.csv(cleanpath)
na = floor(length(dfr$x[dfr$x>0])*fraction)
nb = floor(length(dfr$x[dfr$x<=0])*fraction)
dfg= read_feather(genpath)

if (name=="jl_math"){
   dfg$x = round(dfg$x, 2)
}
if (name=="mats_math"){
   dfg$x = round(dfg$x, 0)
}
df.ga = dfg[dfg$x>0, ]
df.gb = dfg[dfg$x<=0, ]

for (i in 1:runs) {
  dataname = paste0(name, "_", sample(1:100000, 1))
  chunkpath = paste0("data/chunks/", dataname, ".feather")
  runpath = paste0("output/", name, "/", dataname, "run.csv\n")

  dfa <- df.ga[sample(1:nrow(df.ga), na), ]
  dfb <- df.gb[sample(1:nrow(df.gb), nb), ]
  df <- data.frame(rbind(dfa, dfb))
  print(chunkpath)
  write_feather(df, chunkpath)
  submit_slurm(name, chunkpath, runpath)
}
