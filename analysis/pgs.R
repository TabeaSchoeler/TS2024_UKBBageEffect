#!/usr/bin/Rscript
 
args = commandArgs(trailingOnly=TRUE)
HOME=args[1]
task=args[2]
pheno=args[3]

source(paste0(HOME, "/analysis/input.R"))
source(paste0(HOME, "/analysis/functions.R"))


start_time <- Sys.time()

con <- file(paste0(HOME, "/output/log/",task, "_", pheno,"_R.log"))
sink(con, append=TRUE, type="output")

print("Remove old file")
file.remove(paste0(OUT, "/output/pgs/scores/", pheno, ".rds"))

print(paste0("Get PGS for ", pheno))
genPRS(pheno=pheno)


print("======= RUNNING TIME =======")
end_time <- Sys.time()
print(end_time - start_time)
print("============================")
print("All data saved on cluster")
