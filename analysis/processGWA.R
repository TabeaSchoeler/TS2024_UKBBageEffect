#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
HOME=args[1]
pheno=args[3]
#pheno="bmi_0"
source(paste0(HOME, "/analysis/input.R"))
source(paste0(HOME, "/analysis/functions.R"))
library(GenomicSEM)


print("Process GWA")
con <- file(paste0(HOME, "/output/log/processGWA.",pheno, "_R.log"))
sink(con, append=TRUE, type="output")

print("Remove old file")
file.remove(paste0(OUT, "/output/regenie/genofiles/processed/", pheno))

print('================ Start processing REGENIE output ================')
GWAout=readGWA(model=pheno)
saveGWA(df=GWAout, name=pheno)

print("All GWA processed and saved on cluster")        
nSNP=NROW(GWAout)

print("================ Clump data =========================")
gwaClump=clumpData(clump_input=subset(GWAout, P<5e-8), name=pheno)
gwaClump$nSNP=nSNP

if(ncol(gwaClump) <= 4){
    cols=c("CHR", "POS", "SNP", "A1", "A2", "BETA", "SE", "P", "INFO", "N", "EAF", "MAF", "pheno", "nSNP")
    gwaClump=as.data.frame(matrix(ncol=length(cols), nrow=0))
    colnames(gwaClump)=cols
}
saveRDS(gwaClump, paste0(HOME, "/output/clump/clump_",pheno, ".rds"))

print("================ Munge data =========================")
mungeData(name=pheno, path=paste0(OUT, "/output/regenie/genofiles/processed/"))

print("All output saved")
