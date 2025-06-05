#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
HOME=args[1]
task=args[2]
source(paste0(HOME, "/analysis/input.R"))
source(paste0(HOME, "/analysis/functions.R"))


library(lme4)
library(survey)
library(svylme)

con <- file(paste0(HOME, "/output/log/",task,"_R.log"))
sink(con, append=TRUE, type="output")

start_time <- Sys.time()

print("Reading variable names")
longitudinalVar=readRDS(paste0(OUT,"/output/rds/phenoInc.rds"))

print("Remove old files")
file.remove(paste0(HOME,"/output/rds/ageEff.rds")) 
file.remove(paste0(HOME,"/output/rds/snpEff.rds")) 

print("=========== Prepare data =============")
UKBB_C=readRDS(paste0(OUT,"/output/rds/UKBB_C.rds")) # Use pre-COVID samle
UKBB_CL=lapply(longitudinalVar, function(x) prepareAge(df=UKBB_C, var=x, ID="eid"))
names(UKBB_CL)=longitudinalVar

UKBBFU=readRDS(paste0(OUT,"/output/rds/UKBB.rds")) # use the whole UKB sample
UKBB_L=lapply(longitudinalVar, function(x) prepareAge(df=UKBBFU, var=x, ID="eid"))
names(UKBB_L)=longitudinalVar

print("=========== Get age effects (pre-COVID)=============")
ageEffL_preC=lapply(longitudinalVar, function(x) ageModel(var=x, data=UKBB_CL[[x]], datAge=UKBB_C)) #  pre Covid
ageEff_pc=do.call(rbind, ageEffL_preC)
ageEff_pc$var=paste0(ageEff_pc$var, "_c")

print("=========== Get age effects (all)=============")
ageEffL=lapply(longitudinalVar, function(x) ageModel(var=x, data=UKBB_L[[x]], datAge=UKBBFU)) #  all
ageEff=do.call(rbind, ageEffL)

# age estimates across the two samples
ageCOVID=rbind(subset(ageEff_pc, type=="within_age"), subset(ageEff, type=="within_age"))

print("=========== Get age-varying genetic effects =============")
gene=readRDS(paste0(OUT,"/output/rds/gene.rds"))
clumpA=readRDS(paste0(HOME,"/output/rds/clump.rds"))
clump=subset(clumpA, !type %in% c("marginal", "interaction_within_c") & p_fdr <0.05)
snpEffL=lapply(longitudinalVar, function(x) ageModel(var=x, data=UKBB_L[[x]], datAge=UKBBFU, snp=subset(clump, pheno==x), gene=gene)) #  datH=hse
snpEff=do.call(rbind, snpEffL[!sapply(snpEffL, function(x) is.null(x) || all(is.na(x)))])


print("Save results")
attr(ageEff, "timestamp") <- Sys.time()
saveOutput(object=ageEff, label="ageEff", upload="yes")
attr(snpEff, "timestamp") <- Sys.time()
saveOutput(object=snpEff, label="snpEff", upload="yes")
attr(ageCOVID, "timestamp") <- Sys.time()
saveOutput(object=ageCOVID, label="ageCOVID", upload="yes")


end_time <- Sys.time()

print("======= RUNNING TIME =======")
print(end_time - start_time)
print("============================")

print("All data saved on cluster")
