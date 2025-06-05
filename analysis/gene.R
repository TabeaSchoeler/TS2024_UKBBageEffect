#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
HOME=args[1]
task=args[2]
source(paste0(HOME, "/analysis/input.R"))
source(paste0(HOME, "/analysis/functions.R"))

con <- file(paste0(HOME, "/output/log/",task,"_R.log"))
sink(con, append=TRUE, type="output")

start_time <- Sys.time()

library(bigsnpr)
library(dplyr)
clump=readRDS(paste0(HOME,"/output/rds/clump.rds"))
clumpS=subset(clump, !type %in% c("marginal", "interaction_within_c") & p_fdr <0.05)
alleleL=lapply(levels(as.factor(clumpS$SNP)), function(x) readSNP(df=clumpS, snp=x))
alleleDF=Reduce(function(x,y) merge(x = x, y = y, by = "eid", all=TRUE ),  alleleL)

print("Remove individuals for genetic analyses")
dfIND=alleleL[[1]]
dfIND$pheno_pheno=1
remIND=prepGWAdat(df=dfIND, var="pheno", saveLabel="pheno", save="no", relatives="no")
remIND$eid=remIND$IID
remINDs=subset(remIND, eid>0, select=c("eid", "SEX", "batch", paste0("PC", 1:20)))
alleleSel=merge(alleleDF, remINDs, by="eid", all.y=TRUE)
saveRDS(alleleSel , paste0(OUT,"/output/rds/gene.rds"))

print("Get direction of age-varying genetic effects")
UKB=readRDS(paste0(OUT,"/output/rds/UKBB.rds"))
directionBetweenL=lapply(levels(as.factor(clumpS$pheno)), function(x) directionI(df=UKB, var=x, snpList=subset(clumpS, pheno==x)$SNP, gene=alleleSel, type="interaction_between"))
directionWithinL=lapply(levels(as.factor(clumpS$pheno)), function(x) directionI(df=UKB, var=x, snpList=subset(clumpS, pheno==x)$SNP, gene=alleleSel, type="interaction_within"))
directionSum=rbind(do.call(rbind, directionBetweenL), do.call(rbind, directionWithinL))

print("Save results")
attr(directionSum, "timestamp") <- Sys.time()
saveOutput(object=directionSum, label="direction", upload="yes")


end_time <- Sys.time()
print("======= RUNNING TIME =======")
print(end_time - start_time)
print("============================")

print("All data saved on cluster")