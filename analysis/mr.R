#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
HOME=args[1]
task=args[2]
source(paste0(HOME, "/analysis/input.R"))
source(paste0(HOME, "/analysis/functions.R"))


library(TwoSampleMR)

con <- file(paste0(HOME, "/output/log/",task,"_R.log"))
sink(con, append=TRUE, type="output")

start_time <- Sys.time()


print(" ==============  START MR (causal effects in TwoSampleMR) ========================= ")
mr1_long=funMR(exposure="bmi", outcome="SBP", model="change", nTop=20)
mr1_cs=funMR(exposure="bmi", outcome="SBP", model="0_interaction", nTop=20)
mr1_mean=funMR(exposure="bmi", outcome="SBP", model="0")
mrOut1=rbind(mr1_long[["mr"]], mr1_cs[["mr"]], mr1_mean[["mr"]])
mrOut1$pval=2*pnorm(-abs(mrOut1$b/mrOut1$se))
mrOut1$uCI=mrOut1$b + 1.96 * mrOut1$se
mrOut1$lCI=mrOut1$b - 1.96 * mrOut1$se

print("Extract SNP effects")
snpsOut=rbind(mr1_long[["snps"]], mr1_cs[["snps"]],mr1_mean[["snps"]])


print("Save results")
attr(mrOut1, "timestamp") <- Sys.time()
saveOutput(object=mrOut, label="mr", upload="yes")

attr(snpsOut, "timestamp") <- Sys.time()
saveOutput(object=snpsOut, label="mrDat", upload="yes")


end_time <- Sys.time()


print("======= RUNNING TIME =======")
print(end_time - start_time)
print("============================")

print("All data saved on cluster")

