
################################################################
# ====================== Create output tables ==================
################################################################

# +++++++++++++++++++++++++ RUN FUNCTIOND AND DEFINE STYLES ++++++++++++++++++++++++++++++++++++++++++++++++++
ColStart=2
RowHeader=2
RowSubheaderStart=3
RowSubheaderEnds=3
RowTable=4

library(openxlsx)

# Create info text
createInfo=function(dataInfoPath){
  datOut=read.csv(dataInfoPath,header=T)
  datOut$X=NULL
  datOut_merged=paste0(datOut[,1],": " ,datOut[,2])
  return(datOut_merged)
}

# Headerstyle
hs1 <- createStyle(halign = "CENTER", textDecoration = "Bold",
                   border = "Bottom", fontColour = "black", fgFill = "white")

h_info <- createStyle(halign = "left", textDecoration = "Bold",
                      border = "Bottom", fontColour = "black", fgFill = "white")

addTable=function(sheet, table){
  writeDataTable(wb, sheet, table, headerStyle=hs1, tableStyle = "TableStyleLight1",
                 startRow = RowTable, startCol = ColStart)
  setColWidths(wb, sheet, cols=2:10, widths = 15)
}

# HEADER
headerFunc=function(TITLE, sheet){
  writeData(wb, sheet = sheet, TITLE,
            colNames = FALSE, rowNames = FALSE,
            startCol = ColStart, startRow = RowHeader)
}
# INFO ROW
InfoFunc=function(TITLE, sheet){
  writeData(wb, sheet = sheet, TITLE,
            colNames = FALSE, rowNames = FALSE,
            startCol = ColStart, startRow = RowSubheaderStart)
}

# Create new workbook
setwd(paste0(HOME,"/results/tables/"))
wb <- openxlsx::createWorkbook()


# ============== Table 1 =================
# gwas significant findings
ageEff_supp=ageEff
ageEff_supp$beta_CI=paste0(round(ageEff_supp$beta, 3), " (", round(ageEff_supp$lCI, 3), "; ", round(ageEff_supp$uCI, 3), ")")
ageEff_supp$P_c=formatC(ageEff_supp$P, 4)
ageEffsupp_inc=subset(ageEff_supp, select=c(l_clean, typeC, beta_CI, P_c, n))
colnames(ageEffsupp_inc)=c("Aging phenotype", "Model / UKBB sample", "beta (95% CI)", "P", "(Effective) sample size")
head(ageEffsupp_inc)
ageEfft="Supplementary Data 1"
addWorksheet(wb, ageEfft)
# Add parameters
title_name="Supplementary Data 1. Cross-sectional and longitudinal convergence of aging effects."
sheet=ageEfft
table=ageEffsupp_inc
Info_text="The source data used in figure 2 (main manuscript). Of note, for estimates obtained from inverse probability weighting, the corresponding n was replaced by the effective sample size that accounts for the unequal contribution per observation."
# Run functions
addTable(sheet, table)
headerFunc(title_name, sheet)
InfoFunc(Info_text, sheet)


# ============== Table 2 =================
# gwas significant findings obtained from REGENIE
clump_supp=clump
clump_supp$beta_CI=paste0(round(clump_supp$BETA, 3), " (", round(clump_supp$lCI, 3), "; ", round(clump_supp$uCI, 3), ")")
clump_supp$P_c=formatC(clump_supp$P, 4)
clump_supp$P_bf=formatC(clump_supp$p_fdr, 4)
clump_supp$typeC=revalCOL(var=clump_supp$type)
clump_supp$N
clump_supp_inc=subset(clump_supp, select=c(typeC, label_clean, SNP, gene, beta_CI, P_c, P_bf, N))
colnames(clump_supp_inc)=c( "Model", "Aging phenotype", "SNP", "Nearest gene","beta (95% CI)", "P", "P (Bonferroni corrected)","Sample size")
ageEfft="Supplementary Data 2"
addWorksheet(wb, ageEfft)
# Add parameters
title_name="Supplementary Data 2. Cross-sectional and longitudinal identification of age-varying genetic effects"
sheet=ageEfft
table=clump_supp_inc
Info_text="The source data used in figure 3 (main manuscript). All estimates were obtained from participants of European genetic ancestry using REGENIE v3.2.6"
# Run functions
addTable(sheet, table)
headerFunc(title_name, sheet)
InfoFunc(Info_text, sheet)

# ============== Table 3 =================
snpEff_supp=snpEff
snpEff_supp$beta_CI=paste0(round(snpEff_supp$beta, 3), " (", round(snpEff_supp$lCI, 3), "; ", round(snpEff_supp$uCI, 3), ")")
snpEff_supp$P_c=formatC(snpEff_supp$P, 4)
snpEff_supp$P_bf=formatC(snpEff_supp$p_fdr, 4)
snpEff_supp$typeC=revalCOL(var=snpEff_supp$type)
snpEff_supp_inc=na.omit(subset(snpEff_supp, select=c(typeC, label_clean, SNP, beta_CI, P_c, n)))
colnames(snpEff_supp_inc)=c( "Model", "Aging phenotype", "SNP", "beta (95% CI)", "P", "Sample size")
snpEffEfft="Supplementary Data 3"
addWorksheet(wb, snpEffEfft)
# Add parameters
title_name="Supplementary Data 3. Sources of bias contributing to discrepancies in age-varying genetic effects"
sheet=snpEffEfft
table=snpEff_supp_inc
Info_text=paste0("The source data used in figure 5 (main manuscript). All estimates were obtained from unrelated participants of European genetic ancestry. ")
                 
# Run functions
addTable(sheet, table)
headerFunc(title_name, sheet)
InfoFunc(Info_text, sheet)


################################################################
# ====================== Export table ==================
################################################################

# Create new styles
s <- createStyle(fgFill = "#FFFFFF")
h_info <- createStyle(halign = "left",
                      border = "BOTTOM", fontColour = "black", fgFill = "white", fontSize=16, textDecoration = "Bold", numFmt="TEXT", borderColour = "black")
info_info <- createStyle(halign = "left",
                         border = NULL, fontColour = "black", fgFill = "white", fontSize=14, textDecoration = NULL, numFmt="TEXT", wrapText=TRUE)
# Run loop
for(curr_sheet in names(wb)){
  addStyle(wb,sheet = curr_sheet, s, cols=1:40, rows=1:2000, gridExpand = TRUE)
  setColWidths(wb, sheet = curr_sheet, cols=1:40, widths = 20)
  addStyle(wb,sheet = curr_sheet, h_info, cols=ColStart:20, rows=RowHeader, gridExpand = TRUE)
  addStyle(wb,sheet = curr_sheet, info_info, cols=ColStart:5, rows=RowSubheaderStart, gridExpand = TRUE)
  mergeCells(wb,sheet = curr_sheet, cols = 2:100, rows = RowSubheaderStart:RowSubheaderEnds)
}

library( openxlsx)
openxlsx::saveWorkbook(wb, paste0(HOME,"/results/tables/geneageUKBB.xlsx"), overwrite = TRUE)
# Open File
openXL(wb)
