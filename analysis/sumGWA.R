#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
HOME=args[1]
task=args[2]

source(paste0(HOME, "/analysis/input.R"))
source(paste0(HOME, "/analysis/functions.R"))
library(GenomicSEM)

print("Remove old file")
file.remove(paste0(HOME,"/output/rds/clump.rds"))

con <- file(paste0(HOME, "/output/log/",task,"_R.log"))
sink(con, append=TRUE, type="output")

start_time <- Sys.time()

print("Summarize GWA results")
readIn=read.table(paste0(OUT, "/output/regenie/phenofiles/batchNameG"), head = F)$V1
print("Remove phenotypes with unestimated interaction effects")
phenoAvailL=lapply(readIn, function(x) data.table::fread(paste0(OUT, "/output/regenie/genofiles/processed/", x), nrows=20))
phenoAvail=data.frame(pheno=readIn, avail=do.call(rbind, lapply(phenoAvailL, function(x) NROW(x))))
readIn=subset(phenoAvail, avail==20)$pheno
paste0(NROW(subset(phenoAvail, avail==0)$pheno), " not available")

print("Generated GWA statistics")
print(table(readIn))

print("Extract clumped estimates")
pheno=stringr::str_remove(grep("0_interaction", readIn, value = TRUE), "_0_interaction")
clumpL=lapply(pheno, function(x) compareSNP(var=x))
clumpDF=do.call(rbind, Filter(function(x) !is.null(x) && nrow(x) > 0, clumpL))

print("Get estimates from pre-COVID sample")
clumpC=lapply(pheno, function(x) compareSNP(var=x, s="_c"))
clumpCDF=do.call(rbind, Filter(function(x) !is.null(x) && nrow(x) > 0, clumpC))

print("Upload all output")
clumpS=rbind(clumpDF, subset(clumpCDF, type=="interaction_within_c"))
attr(clumpS, "timestamp") <- Sys.time()
saveOutput(object=clumpS, label="clump", upload="yes")


end_time <- Sys.time()

print("======= RUNNING TIME =======")
print(end_time - start_time)
print("============================")

print("All data saved on cluster")

