#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
HOME=args[1]
task=args[2]
source(paste0(HOME, "/analysis/input.R"))
source(paste0(HOME, "/analysis/functions.R"))


con <- file(paste0(HOME, "/output/log/bias_",task,"_R.log"))
sink(con, append=TRUE, type="output")


print("Read in UKBB data harmonized with HSE")
UKBBHSE=readRDS( paste0(OUT,"/data/datHSEUKBB.rds"))
UKBBHSE$ID=paste0(UKBBHSE$eid, "_", UKBBHSE$sampleName)

print("Get IDs with available follow up assessment")
phenoInc=fread(paste0(OUT, "/output/pgs/phenofiles/scoreName"), header=FALSE)$V1
UKBBFU=readRDS(paste0(OUT,"/output/rds/UKBB.rds"))
changeDat=lapply(phenoInc, function(x) changeFunc(var=x, df=UKBBFU, return="data"))
names(changeDat)=phenoInc
listID=lapply(phenoInc, function(x) fuIDL(var=x, dfIn=changeDat))
listDF=unique(do.call(rbind, listID))
listDF$miss0=0


print("Merge IDs with dataframe")
UKBB=merge(UKBBFU, listDF, by="eid", all.x=T)
UKBB$missFU = ifelse(is.na(UKBB$miss0)==TRUE, 1, UKBB$miss0)
UKBB$ID=paste0(UKBB$eid, "_UKBB")

print("Add death data")
death=readRDS(paste0(OUT,"/output/rds/death.rds"))
UKBBd=merge(UKBB, subset(death, select=c(eid, death_age)), by="eid", all.x=T)
UKBBd$missDEATH=ifelse(is.na(UKBBd$death_age)==TRUE, 0, 1 ) # 1 = dead
UKBBd$missDEATHc=ifelse(UKBB$missFU==0 & UKBBd$missDEATH==1, 0, UKBBd$missDEATH) # count only death occuring prior to FU assessments


print("Sample with available FU data, compare with HSE")
UKBBHSEFU=rbind(subset(UKBBHSE, sampleName=="HSE"),
                subset(UKBBHSE, ID %in% subset(UKBB, missFU==0)$ID)) # select individuals with available FU

library("readxl")
varList <- readxl::read_excel(paste0(HOME, "/data/variableListWeighting.xlsx"))

print("read in variables for weighting")
vars=readRDS(paste0(HOME, "/data/weightVariables.rds"))

print("================= Start participation (FU versus HSE ) LASSO ======================")
lassoPBFU=runLasso(data=UKBBHSEFU, iteration="participationFU", varInc=vars, variableList=varList)
print("================= Finnished participation LASSO ======================")
     
#saveRDS(lassoPBFU, paste0(OUT,"/output/rds/lassoPBFU.rds"))
#file.remove(paste0(OUT,"/output/rds/lassoPBFU.rds"))
#lassoPBFU=readRDS(paste0(OUT,"/output/rds/lassoPBFU.rds"))

print("Combine weights")
UKBw=subset(UKBBHSE, sample==1)
datW=merge(lassoPBFU[["data"]], data.frame(eid=UKBw$eid, PW_participation_BL=UKBw$propensity.weight.normalized), by="eid", all=TRUE)

print("Save results")
attr(datW, "timestamp") <- Sys.time()
saveOutput(object=datW, label="weights", upload="yes")

print("Get descriptive data")
nAll=NROW(UKBBd)
nFU=NROW(subset(UKBBd, missFU==0))
perFU=(nFU/nAll)*100
nLoss=NROW(subset(UKBBd, missFU==1))
nLossDead=NROW(subset(UKBBd, missFU==1 & missDEATHc==1))
percLossDead=(nLossDead/nLoss)*100

lossINFO=data.frame(label=c("nAll", "perFU","nFU", "nLoss", "nLossDead", "percLossDead"),
           est=c(nAll, perFU, nFU, nLoss, nLossDead, percLossDead))

print("Save results")
attr(lossINFO, "timestamp") <- Sys.time()
saveOutput(object=lossINFO, label="lossINFO", upload="yes")


end_time <- Sys.time()


print("======= RUNNING TIME =======")
print(end_time - start_time)
print("============================")

print("All data saved on cluster")
