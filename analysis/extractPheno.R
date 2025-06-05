#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
HOME=args[1]
task=args[2]
source(paste0(HOME, "/analysis/input.R"))
source(paste0(HOME, "/analysis/functions.R"))


library("psych") 
library("data.table")
library(stringr)
library(plyr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(ggpubr)

start_time <- Sys.time()
con <- file(paste0(HOME, "/output/log/",task,"_R.log"))
sink(con, append=TRUE, type="output")

print("Remove old file")
file.remove(paste0(OUT,"/output/rds/UKBB.rds"))
file.remove(paste0(OUT,"/output/rds/UKBB_C.rds")) 
file.remove(list.files(paste0(OUT,"/output/regenie/phenofiles"), full.names = TRUE))

print("Select variables")
dfAll=readVarNames()
dfExtract=readVarNames(do="extract")
longDF=readVarNames(do="inLong")
long=longDF$label
baseline=subset(dfAll, baseline=="yes")$label


print("Read in individual-level data")
redIn=c(679895) # most recent (largest file)
UKBBList=lapply(redIn, function(x) fread(paste0(UKBB, "/phenotypes/ukb", x, ".csv"))) # ,  nrows=10000
names(UKBBList)=redIn

print("Process data")
UKBBL=lapply(dfExtract$label, function(x) extractPheno(var=x, varFile=dfAll, dfL=UKBBList, return="df"))
UKBBind=Reduce(function(x,y) merge(x = x, y = y, by = "eid", all=T ),  UKBBL)

print("Check coding of the variables")
checkCoding(varInc=longDF)

print(" ========================================================= ")
print(" ===================== Recode variables ================== ")
print(" ========================================================= ")
recodeL=lapply(levels(as.factor(c(long, baseline))), function(x) recodeVar(df=UKBBind, var=x, varInfo=dfAll))
UKBBr=Reduce(function(x,y) merge(x = x, y = y, by = "eid", suffixes=c("", "_rem"), all=T ),  recodeL)

print("Excluded individuals (out of age exclusion criteria):")
UKBBa=subset(UKBBr, age_0 >= 40 & age_0 <= 69)
print(NROW(UKBBr)-NROW(UKBBa))

print("=========== Merge with sampling weights =============")
weights=readRDS(paste0(OUT, "/data/datHSEUKBB.rds"))
weights$w=weights$propensity.weight.normalized
ws=subset(weights, sampleName=="UKBB", select=c(eid, w))
UKBBw=merge(UKBBa, ws, by="eid", all.x=TRUE)
rm(UKBBr, UKBBa)


print(" ========================================================= ")
print(" ==================== Get follow up data ================= ")
print(" ========================================================= ")
print("All UKB participants")
UKBBfuL=lapply(long, function(x) deriveFU(var=x, df=UKBBw, varDF=dfAll))
UKBBFU=Reduce(function(x,y) merge(x = x, y = y, by = "eid" , suffixes=c("", "_rem"), all=T),  UKBBfuL)
saveRDS(UKBBFU , paste0(OUT,"/output/rds/UKBB.rds"))
 
print("Pre-COVID UKB participants")
UKBB_CL=lapply(long, function(x) deriveFU(var=x, df=UKBBw, varDF=dfAll, sub="COVID"))
UKBB_C=Reduce(function(x,y) merge(x = x, y = y, by = "eid" , suffixes=c("", "_rem"), all=T),  UKBB_CL)
saveRDS(UKBB_C , paste0(OUT,"/output/rds/UKBB_C.rds"))

print(" ========================================================= ")
print(" ======== INFO on follow-up participation ================ ")
print(" ========================================================= ")
apcINFO_a=do.call(rbind, lapply(long, function(x) apcI(df=UKBBFU, var=x)))
print(paste0("Exlcuded phenotypes as <50,000 FU observations: ", paste0(levels(as.factor(subset(apcINFO_a, n_FU<50000)$var)), collapse=", ")))
apcINFO=subset(apcINFO_a, n_FU>50000)
phenoInc=levels(as.factor(subset(apcINFO_a, n_FU>50000)$var))
print(paste0("Number of phenotypes included: ", NROW(phenoInc)))

attr(apcINFO, "timestamp") <- Sys.time()
saveOutput(object=apcINFO, label="apcINFO", upload="yes")
saveOutput(object=phenoInc, label="phenoInc", upload="yes")

print(" ========================================================= ")
print(" ======== Prepare data for change scores ================= ")
print(" ========================================================= ")
print("All UKB participants")
changeDat=lapply(phenoInc, function(x) changeFunc(var=x, df=UKBBFU, return="data"))
names(changeDat)=phenoInc

print("Pre-COVID UKB participants")
changeDat_C=lapply(phenoInc, function(x) changeFunc(var=x, df=UKBB_C, return="data"))
names(changeDat_C)=phenoInc

print(" ========================================================= ")
print(" ==================== Normality plots ==================== ")
print(" ========================================================= ")
normPL=lapply(phenoInc, function(x) normPlot(df=changeDat[[x]], var=x))
normC <- Filter(Negate(is.null), normPL)
normPC=ggpubr::ggarrange(plotlist=normC, ncol=1)
ggsave(file=paste0(HOME, "/output/rds/normP.pdf"), plot = normPC, width = 23, height = 120, units = "cm", bg='transparent', limitsize = FALSE)


print(" ========================================================= ")
print(" ========Prepare phenotype file for regenie ============== ")
print(" ========================================================= ")
print("Save REGENIE files for change")
lapply(phenoInc, function(x) prepGWAdat(df=changeDat[[x]], var=x, saveLabel="change"))
lapply(phenoInc, function(x) prepGWAdat(df=changeDat_C[[x]], var=x, saveLabel="change", sample="_c"))
lapply(phenoInc, function(x) prepGWAdat(df=UKBBFU, var=x, saveLabel="0"))


write.table(data.frame(pheno=c(paste0(phenoInc, "_change"), paste0(phenoInc, "_change_c"), paste0(phenoInc, "_0"))),
            file= paste0(OUT, "/output/regenie/phenofiles/batchName"),
            sep="\t",
            row.names = FALSE,
            col.names=F,
            quote=F)


print("Files in output folders")
print(list.files(paste0(OUT,"/output/regenie/phenofiles/")))


end_time <- Sys.time()

print("======= RUNNING TIME =======")
print(end_time - start_time)
print("============================")

print("All data saved on cluster")