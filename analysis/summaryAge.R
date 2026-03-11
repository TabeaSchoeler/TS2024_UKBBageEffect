#######################################################
# =================== FUNCTIONS =======================
#######################################################
rm(list = ls())
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")
HOME=getwd()
source(paste0(HOME, "/analysis/summaryAgeFunc.R"))
library(metafor)
setwd("..");setwd("..")
MD=getwd()
MDP=paste0(MD, "/mountainduck.localized/TS2024_ageEffect")
file.copy(from=list.files(paste0(MDP, "/analysis"), full.names = TRUE), to=paste0(HOME, "/analysis/MD"), overwrite = TRUE)

print("Copy over variable names file")
file.copy(from=paste0(MDP, "/data/variableAge.xlsx"), to=paste0(HOME, "/data/"), overwrite = TRUE)

print("Copy over result files")
copyFiles=c("ageEff", "snpEff","clump", "direction", "apcINFO", "ageCOVID", "rDil", "mr", "mrDat", "lossINFO", "snpPGS")
cDat=lapply(copyFiles, function(x) importFiles(file=x))



#########################################################################################
# ================================= PHENOTYPIC ANALYSES =================================
#########################################################################################
# ================== Age effect ======================
ageEff=readRDS( paste0(HOME,"/output/rds/ageEff.rds"))
ageEff$uCI=ageEff$beta + 1.96 * ageEff$se
ageEff$lCI=ageEff$beta - 1.96 * ageEff$se
ageEffC=recodeLabel(df=ageEff, mergeBy="var")
ageEff$l_clean =ageEffC$label_clean
ageEff$typeC=revalCOL(var=ageEff$type)
ageEff$colour <- factor(ageEff$typeC, labels = 
                          c("darkblue",  "steelblue3", "darkgreen", "cyan4", "darkorange3",  "palevioletred4", "#bc4f5e"))

mOrder=data.frame(var=levels(as.factor(ageEff$var)), ID=seq(1,NROW(levels(as.factor(ageEff$var))), 1))
mOrderL=recodeLabel(df=mOrder, mergeBy="var")
mOrder$label_clean =mOrderL$label_clean

# Plot 1 (longitudinal and cross-sectional age effects)
ageEffP1=merge(mOrder, subset(ageEff, type %in% c("within_age", "between_age")), by="var", all.x=TRUE)

agePlotP1 <- ggplot(ageEffP1, aes(x = beta, y =fct_reorder(label_clean, -ID), colour=typeC)) +
  geom_line(colour="black") +
  geom_pointrange(aes(xmin = lCI, xmax = uCI), size = 1) +
  scale_colour_manual("",values = levels(droplevels(ageEffP1$colour)), label=levels(droplevels(ageEffP1$typeC)) ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Reference line at beta = 0
  labs(x = expression(paste( beta[CS], " / ", beta[L] )), y = "Phenotype", color = "Variable", y = NULL) +
  theme_minimal() +
  themeAll +
  guides(color = guide_legend(nrow = 2)) +
  ggtitle("Model A1/A2")  +
  labs(subtitle = "Age effects") 
agePlotP1

# Plot 2: Cohort effects
ageEffP2=merge(mOrder, subset(ageEff, type %in% c("cohort")), by="var", all.x=TRUE)
agePlotP2 <- ggplot(ageEffP2, aes(x = beta, y =fct_reorder(label_clean, -ID), colour=typeC)) +
  geom_pointrange(aes(xmin = lCI, xmax = uCI), size = 1) +
  scale_colour_manual("",values = levels(droplevels(ageEffP2$colour)), label=levels(droplevels(ageEffP2$typeC)) ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Reference line at beta = 0
  labs(x = expression(paste( beta[C])), y = "", color = "") +
  themeAge +
  themeAll +
  ggtitle("Model A3") +
  labs(subtitle = "Year of birth effects") 
agePlotP2


# Plot 3: Non-linear effects
ageEffP3=merge(mOrder, subset(ageEff, type %in% c("non_linear_age")), by="var", all.x=TRUE)
agePlotP3 <- ggplot(ageEffP3, aes(x = beta, y =fct_reorder(label_clean, -ID), colour=typeC)) +
  geom_pointrange(aes(xmin = lCI, xmax = uCI), size = 1) +
  scale_colour_manual("",values = levels(droplevels(ageEffP3$colour)), label=levels(droplevels(ageEffP3$typeC)) ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Reference line at beta = 0
  labs(x = expression(paste(  beta[Q] )), y = "Phenotype", color = "Variable",  y = NULL) +
  themeAge +
  themeAll +
  guides(color = guide_legend(nrow = 2)) +
  ggtitle("Model A4") +
  labs(subtitle = "Non-linear age effects") 
agePlotP3

# Plot 4: Weighted cross-sectional age effects
ageEffP4=merge(mOrder, subset(ageEff, type %in% c("between_age", "age_weighted", "between_age_FU")), by="var", all.x=TRUE)
agePlotP4 <- ggplot(ageEffP4, aes(x = beta, y =fct_reorder(label_clean, -ID), colour=typeC)) +
  geom_line(colour="black") +
  geom_pointrange(aes(xmin = lCI, xmax = uCI), size = 1) +
  scale_colour_manual("",values = levels(droplevels(ageEffP4$colour)), label=levels(droplevels(ageEffP4$typeC)) ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Reference line at beta = 0
  labs(x = expression(paste( beta[W], " / ", beta[CS] )), y = "Phenotype", color = "Variable",  y = NULL) +
  themeAge +
  themeAll +
  guides(color = guide_legend(nrow = 3)) +
  ggtitle("Model A5") +
  labs(subtitle = "Participation effects") 
agePlotP4


# correlations
ageEff$vi=ageEff$se^2 # get variance
d1=subset(ageEff, type %in% c("age_weighted"), select=c(beta, vi, var)); d1$beta_weighted=d1$beta;d1$beta=NULL;d1$vi_weighted=d1$vi;d1$vi=NULL
d2=subset(ageEff, type %in% c("within_age"), select=c(beta, vi, var)); d2$beta_within=d2$beta;d2$beta=NULL;d2$vi_within=d2$vi;d2$vi=NULL
d3=subset(ageEff, type %in% c("between_age"), select=c(beta, vi, var)); d3$beta_between=d3$beta;d3$beta=NULL;d3$vi_between=d3$vi;d3$vi=NULL
d4=subset(ageEff, type %in% c("cohort"), select=c(beta, vi, var)); d4$beta_cohort=d4$beta;d4$beta=NULL;d4$vi_cohort=d4$vi;d4$vi=NULL
d5=subset(ageEff, type %in% c("between_age_FU"), select=c(beta, vi, var)); d5$beta_between_FU=d5$beta;d5$beta=NULL;d5$vi_between_FU=d5$vi;d5$vi=NULL
d6=subset(ageEff, type %in% c("non_linear_age"), select=c(beta, vi, var)); d6$beta_age2=d6$beta;d6$beta=NULL;d6$vi_age2=d6$vi;d6$vi=NULL

ageComb <- Reduce(function(x, y) merge(x, y, by = "var", all = TRUE), list(d1, d2, d3, d4, d5, d6))
ageComb$beta_diff=ageComb$beta_within - ageComb$beta_between
ageComb$beta_diff_w_fu=ageComb$beta_weighted - ageComb$beta_between_FU
rDil=readRDS( paste0(HOME,"/output/rds/rDil.rds"))
corW=as.numeric(rDil$estimate)
ageComb$vi_diff_w_fu    <-  ageComb$vi_weighted + ageComb$vi_between_FU - 2 * corW * sqrt(ageComb$vi_weighted) * sqrt(ageComb$vi_between_FU)


# R2 per predictor
m3=lm(beta_diff ~ beta_age2 + beta_cohort + beta_diff_w_fu, data=ageComb)
mDiscSTD <- lm.beta::lm.beta(m3)
discDF=as.data.frame(coef(summary(mDiscSTD)))
discDF=data.frame(predictor=rownames(discDF), beta_std=discDF$Standardized, P=discDF$`Pr(>|t|)`)
rel_imp <- relaimpo::calc.relimp(m3, type = "lmg") # Lindeman, Merenda & Gold method
dfR2=data.frame(r2=rel_imp@lmg)
totVar=round(sum(dfR2$r2)*100, 2)
dfR2$predictor=rownames(dfR2)
dfR2$labelT=paste0(round(dfR2$r2 * 100, 1), "%")
r2AddR <- data.frame( r2 = (100-totVar)/100, predictor = "missing_var", labelT="Missing")
dfR2 <- rbind(dfR2, r2AddR)
dfR2$predictor <- factor(dfR2$predictor, levels = c("missing_var", "beta_cohort", "beta_age2","beta_diff_w_fu"))
dfR2$group=1
colR2=c("lightgrey","cyan4","darkorange3", "darkblue")
dfR2m=merge(dfR2, discDF, by="predictor", all.x=TRUE)
saveRDS(dfR2m, paste0(HOME, "/output/rds/expR2pheno.rds"))


# Plot R2
dfR2$labelT=ifelse(dfR2$r2<0.01, "", dfR2$labelT)
R2Plot <- ggplot(dfR2, aes(x = group, y = r2 * 100, fill = predictor)) +
  geom_bar(stat = "identity") +
  labs(title = "Stacked Bar Chart of Explained R²",
       #subtitle = paste0("Total explained \nvariance: ", totVar, "%"),
       x = "",
       y = "Explained variance (%)") +
  theme_minimal() +
  scale_fill_manual("", values = colR2) +
  geom_text(aes(label = labelT),
            position = position_stack(vjust = 0.5),
            size = 5,
            color = "white") +
  scale_y_continuous(position = "right") +  # Moves axis labels, ticks, and title to the right
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major.x = element_blank(), # Remove vertical grid lines
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(hjust = 0.5),  
        axis.text.y = element_text(hjust = 0)) + 
  guides(color = guide_legend(nrow = 3)) +
  ggtitle("Relative \ncontribution")

R2Plot
# ======= FIGURE 2: age effects
ageCombine <- agePlotP1 + agePlotP2  + agePlotP3 + agePlotP4 + R2Plot + plot_layout(ncol = 5, widths = c(1.2, 1, 1,1, 0.4))
saveFigure(fileName=paste0(HOME,"/results/figures/AgePlotSum"), plotName=ageCombine, w=40, h=20 )


###################################################################################
# ============================ REGENIE findings ===================================
###################################################################################
clumpA=readRDS( paste0(HOME,"/output/rds/clump.rds"))
clumpA$uCI=clumpA$BETA + 1.96 * clumpA$SE
clumpA$lCI=clumpA$BETA - 1.96 * clumpA$SE
clumpA$snppheno=paste0(clumpA$SNP, "_", clumpA$pheno)
clumpL=recodeLabel(df=clumpA, mergeBy="pheno")
clumpA$label_clean =clumpL$label_clean
clumpA$label_wide=paste0(clumpA$SNP, " (", clumpA$label_clean, ")")
clumpS=subset(clumpA, !type %in% c("marginal", "interaction_within_c") ) # remove marginal effects and pre-COVID sample estimates
sigSNP=subset(clumpS, p_fdr <0.05)$snppheno
clumpC=subset(clumpS, snppheno %in% sigSNP)

# Map to nearby genes using open Target
clumpOT=openTarget(df=clumpC)
clump=merge(clumpC, subset(clumpOT, select=c(SNP, gene)), by="SNP", all.x=TRUE)
clump$geneS=ifelse(is.na(clump$gene)==TRUE, clump$SNP, clump$gene)
clump$geneC=ifelse(is.na(clump$gene)==TRUE, "", paste0(" / ", clump$gene))
clump$pheno_snp=paste0(clump$SNP, "_", clump$A1, clump$geneC)
clump$desc=0; clump$dimension="g"
clump$type <- factor(clump$type , levels = c("interaction_between", "interaction_within"))
clump$sig_text=ifelse(clump$p_fdr<0.05, "*", "")
cSNP=c( "steelblue3", "palevioletred4", "black")


###################################################################################
# ====================== FREQUENCY OF IDENTIFICATION ==============================
###################################################################################
# age-varying genetic effects
m1R=subset(clump, type=="interaction_between", select=c(pheno,SNP, BETA, label_clean, snppheno, p_fdr, geneS ))
m2R=subset(clump, type=="interaction_within", select=c(BETA, snppheno, p_fdr ))
clumpR=merge(m1R, m2R, by="snppheno",  suffixes = c("_b", "_w") )
clumpR$identification=ifelse(clumpR$p_fdr_b<0.05 & clumpR$p_fdr_w>=0.05 ,"1_between", NA)
clumpR$identification=ifelse(clumpR$p_fdr_b>=0.05 & clumpR$p_fdr_w<0.05 ,"2_within", clumpR$identification)
clumpR$identification=ifelse(clumpR$p_fdr_b<0.05 & clumpR$p_fdr_w<0.05 ,"3_both", clumpR$identification)

# add marginal genetic effects
clumpI=data.frame(snppheno=clumpR$snppheno, label_clean=clumpR$label_clean, gene=clumpR$geneS, identification=clumpR$identification)
clumpM=merge(subset(clumpA, type == c("marginal") & snppheno %in% clumpI$snppheno & p_fdr<0.05), subset(clumpI, select=c(snppheno, gene)), by="snppheno", all.x=TRUE )
clumpM$identification="4_marginal"
clumpMs=data.frame(snppheno=clumpM$snppheno, label_clean=clumpM$label_clean, gene=clumpM$gene, identification=clumpM$identification)
clumpIM=rbind(clumpI, clumpMs)

# ======= FIGURE 3: Identification plot
freqPheno <- clumpIM %>%
  dplyr::group_by(label_clean, identification) %>%
  dplyr::summarise(n = dplyr::n(),  snp = paste0(gene, collapse = "\n")) %>%
  dplyr::arrange(desc(n))

allPheno=data.frame(pheno=levels(as.factor(clumpA$label_clean)))
allPhenoS=subset(allPheno, !pheno %in% freqPheno$label_clean)
freqPhenoR=data.frame(label_clean=allPhenoS$pheno, identification="1_between", n=NA)
freqPhenoC=rbind(freqPheno, freqPhenoR)
freqPhenoC$n=ifelse(freqPhenoC$identification=="4_marginal", freqPhenoC$n*(-1), freqPhenoC$n)
freqPhenoC$snp=ifelse(freqPhenoC$identification=="4_marginal", "", freqPhenoC$snp)

# Plot frequency as a bar chart
freqPhenoC$identification=droplevels(as.factor(freqPhenoC$identification))
levels(as.factor(freqPhenoC$identification))
cSNPM=c(cSNP, "grey")

freqGene <- ggplot(freqPhenoC, aes(x = reorder(label_clean, -abs(n)), y = n, fill = identification)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = snp), 
            position = position_stack(vjust = 0.5),  # Works only if stacked
            size = 1, 
            color = "white", 
            lineheight = 0.8, 
            na.rm = TRUE) +  # Removes missing values
  theme_minimal() +
  scale_y_continuous(labels = abs) +
  labs(title = "",
       x = "",
       y = "Number of identified variants") +
  scale_fill_manual("", values = cSNPM, 
                    labels = c("Cross-sectional", "Longitudinal", "Both (cross-sectional & longitudinal)", "Marginal" ) ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(face = "bold", hjust = 0),
        legend.position = "top",
        plot.margin = margin(t = 0, r = 0, b = 0, l = 2, "cm"))  + ggtitle("B")

saveFigure(fileName=paste0(HOME,"/results/figures/freqGene"), plotName=freqGene, w=22, h=12 )



#########################################################################################
# ============================ SOUCES OF BIAS (SNP LEVEL) ===============================
#########################################################################################
snpEff=readRDS( paste0(HOME,"/output/rds/snpEff.rds"))
snpEff$uCI=snpEff$beta + 1.96 * snpEff$se
snpEff$lCI=snpEff$beta - 1.96 * snpEff$se
snpEff$alpha=ifelse(snpEff$P < 0.05, 1, 0.5)
snpEffL=recodeLabel(df=snpEff, mergeBy="var")
snpEff$label_clean=snpEffL$label_clean
snpEff$label=paste0(snpEff$SNP, "_", snpEff$var, "_", snpEff$type)
snpEff$snp_var=paste0(snpEff$SNP, "_", snpEff$var)
nSNPsFDR=length(levels(as.factor(paste0(snpEff$SNP, "_", snpEff$var))))
snpEff$p_fdr= do.call(rbind, lapply(snpEff$P, function(x) p.adjust(x, "bonferroni",  nSNPsFDR )))
snpEff$label_cleanN =paste0(snpEff$label_clean, " (", snpEff$SNP, ")")
snpEff$type <- factor(snpEff$type, levels = c("snp_interaction_between_lm","snp_age2_lme", "age_snp_weighted", "snp_interaction_between_FU_lm","snp_cohort", "snp_age_lme"))
snpEff$var_snp=paste0(snpEff$var, "_", snpEff$SNP)
snpEff$typeC=revalCOL(var=as.factor(snpEff$type))
snpEff$colour <- factor(snpEff$typeC, labels =  c( "steelblue3", "darkorange3", "darkblue", "darkgreen","cyan4", "palevioletred4"))
mOrderG=data.frame(var=levels(as.factor(snpEff$var)), ID=seq(1,NROW(levels(as.factor(snpEff$var))), 1))
mOrderGL=recodeLabel(df=mOrderG, mergeBy="var")
mOrderG$label_clean =mOrderGL$label_clean

# Plot 1: Gene-by-age effects (cross-sectional versus longitudinal estimates)
snpEffP1=merge(mOrderG, subset(snpEff, type %in% c("snp_age_lme", "snp_interaction_between_lm")), by="var", all.x=TRUE, suffixes = c("", "_rem"))
mergeVar=subset(snpEffP1, select=c(ID, var_snp, label_cleanN))

snpPlotP1 <- ggplot(snpEffP1, aes(x = beta, y =fct_reorder(label_cleanN, -ID), colour=typeC, alpha=alpha, xmin = lCI, xmax = uCI)) +
  geom_pointrange(aes(xmin = lCI, xmax = uCI), size = 1) +
  scale_colour_manual("",values = levels(droplevels(snpEffP1$colour)), label=levels(droplevels(snpEffP1$typeC)) ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(x = expression(paste( delta[CS], " / ",   delta[L])), y = "Phenotype", color = "Variable", y = NULL) +
  theme_minimal() +
  facet_grid(rows = vars(ID), scales = "free_y", space = "free_y") + 
  themeAll +
  guides(color = guide_legend(nrow = 2), alpha = "none" ) +
  ggtitle("Model B1/B2")  +
  labs(subtitle = "Age-varying genetic effects") 


# Plot 2: Cohort-varying genetic effects
snpEffP2=merge(mergeVar, subset(snpEff, type %in% c("snp_cohort")), by="var_snp", all.x=TRUE, suffixes = c("", "_rem"))
snpPlotP2 <- ggplot(snpEffP2, aes(x = beta, y =fct_reorder(label_cleanN, -ID), colour=typeC, alpha=alpha, xmin = lCI, xmax = uCI)) +
  geom_pointrange(size = 1) +
  scale_colour_manual("",values = levels(droplevels(snpEffP2$colour)), label=levels(droplevels(snpEffP2$typeC)) ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") + 
  labs(x = expression(paste( delta[C] )), y = "", color = "") +
  themeAge +
  themeAll+
  facet_grid(rows = vars(ID), scales = "free_y", space = "free_y") +
  guides(color = guide_legend(nrow = 2), alpha = "none" ) +
  ggtitle("Model B3")   +
  labs(subtitle = "Year of birth-varying genetic effects") 
 

# Plot 3: Non-linear age-varying genetic effects
snpEffP3=merge(mergeVar, subset(snpEff, type %in% c("snp_age2_lme")), by="var_snp", all.x=TRUE, suffixes = c("", "_rem"))
snpPlotP3 <- ggplot(snpEffP3, aes(x = beta, y =fct_reorder(label_cleanN, -ID), colour=typeC, alpha=alpha, xmin = lCI, xmax = uCI)) +
  geom_pointrange(size = 1) +
  scale_colour_manual("",values = levels(droplevels(snpEffP3$colour)), label=levels(droplevels(snpEffP3$typeC)) ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") + 
  labs(x = expression(paste( delta[Q] )), y = "Phenotype", color = "Variable", , y = NULL) +
  themeAge +
  themeAll+
  facet_grid(rows = vars(ID), scales = "free_y", space = "free_y") +
  guides(color = guide_legend(nrow = 2), alpha = "none" ) +
  ggtitle("Model B4")  +
  labs(subtitle = "Non-linear age-varying genetic effects") 


# Plot 4: Weighted cross-sectional age effects
snpEffP4=merge(mergeVar, subset(snpEff, type %in% c("snp_interaction_between_lm", "age_snp_weighted", "snp_interaction_between_FU_lm")), by="var_snp", all.x=TRUE, suffixes = c("", "_rem"))
snpPlotP4 <- ggplot(snpEffP4, aes(x = beta, y =fct_reorder(label_cleanN, -ID), colour=typeC, alpha=alpha, xmin = lCI, xmax = uCI)) +
  geom_pointrange(size = 1) +
  scale_colour_manual("",values = levels(droplevels(snpEffP4$colour)), label=levels(droplevels(snpEffP4$typeC)) ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  
  labs(x = expression(paste( delta[W], " / ", delta[CS] )), y = "Phenotype", color = "Variable",  y = NULL) +
  themeAge +
  themeAll+
  facet_grid(rows = vars(ID), scales = "free_y", space = "free_y") +
  guides(color = guide_legend(nrow = 3), alpha = "none") +
  ggtitle("Model B5") +
  labs(subtitle = "Participation effects")


# Explained variance in effect size differences
snpEff$vi=snpEff$se^2 # get variance
d1=subset(snpEff, type %in% c("age_snp_weighted"), select=c(beta, vi, var_snp)); d1$beta_weighted=d1$beta;d1$beta=NULL;d1$vi_weighted=d1$vi;d1$vi=NULL
d2=subset(snpEff, type %in% c("snp_age_lme"), select=c(beta, vi, var_snp)); d2$beta_within=d2$beta;d2$beta=NULL;d2$vi_within=d2$vi;d2$vi=NULL
d3=subset(snpEff, type %in% c("snp_interaction_between_lm"), select=c(beta, vi, var_snp)); d3$beta_between=d3$beta;d3$beta=NULL;d3$vi_between=d3$vi;d3$vi=NULL
d4=subset(snpEff, type %in% c("snp_cohort"), select=c(beta, vi, var_snp)); d4$beta_cohort=d4$beta;d4$beta=NULL;d4$vi_cohort=d4$vi;d4$vi=NULL
d5=subset(snpEff, type %in% c("snp_interaction_between_FU_lm"), select=c(beta, vi, var_snp)); d5$beta_between_FU=d5$beta;d5$beta=NULL;d5$vi_between_FU=d5$vi;d5$vi=NULL
d6=subset(snpEff, type %in% c("snp_age2_lme"), select=c(beta, vi, var_snp)); d6$beta_age2=d6$beta;d6$beta=NULL;d6$vi_age2=d6$vi;d6$vi=NULL

# Combine
snpComb <- Reduce(function(x, y) merge(x, y, by = "var_snp", all = TRUE), list(d1, d2, d3, d4, d5, d6))
snpComb$beta_diff=snpComb$beta_within-snpComb$beta_between
snpComb$beta_diff_fu=snpComb$beta_weighted-snpComb$beta_between_FU
snpComb$vi_diff_fu    <-  snpComb$vi_weighted + snpComb$vi_between_FU - 2 * corW * sqrt(snpComb$vi_weighted) * sqrt(snpComb$vi_between_FU)


# R2 per predictor
m3_snp=lm(beta_diff ~ beta_age2 + beta_cohort + beta_diff_fu, data=snpComb)
as.data.frame(coef(summary(m3_snp)))
rel_imp_snp <- relaimpo::calc.relimp(m3_snp, type = "lmg") # Lindeman, Merenda & Gold method
dfR2snp=data.frame(r2=rel_imp_snp@lmg)
totVar=round(sum(dfR2snp$r2)*100, 2)
dfR2snp$predictor=rownames(dfR2snp)
dfR2snp$labelT=paste0(round(dfR2snp$r2 * 100, 1), "%")
r2AddR <- data.frame( r2 = (100-totVar)/100, predictor = "missing_var", labelT="Missing")
dfR2snp <- rbind(dfR2snp, r2AddR)
dfR2snp$predictor <- factor(dfR2snp$predictor, levels = c("missing_var", "beta_cohort", "beta_age2","beta_diff_w", "beta_diff_fu"))
dfR2snp$predictorC=as.factor(revalCOL(var=dfR2snp$predictor))
dfR2snp$group=1
levels(dfR2snp$predictorC)
colR2=c("lightgrey", "cyan4", "darkorange2", "darkblue", "darkgreen")

mDiscGeneSTD <- lm.beta::lm.beta(m3_snp)
discGeneDF=as.data.frame(coef(summary(mDiscGeneSTD)))
discGDF=data.frame(predictor=rownames(discGeneDF), beta_std=discGeneDF$Standardized, P=discGeneDF$`Pr(>|t|)`)
dfR2Genem=merge(dfR2snp, discGDF, by="predictor", all.x=TRUE)
saveRDS(dfR2Genem, paste0(HOME, "/output/rds/expR2geno.rds"))

R2Plot_snp <- ggplot(dfR2snp, aes(x = group, y = r2 * 100, fill = predictorC)) +
  geom_bar(stat = "identity") +
  labs(title = "",
       x = "",
       y = "Explained variance (%)") +
  theme_minimal() +
  scale_fill_manual("", values = colR2) +
  geom_text(aes(label = labelT),
            position = position_stack(vjust = 0.5),
            size = 5,
            color = "white") +
  scale_y_continuous(position = "right") +  # Moves axis labels, ticks, and title to the right
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major.x = element_blank(), # Remove vertical grid lines
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.y = element_text(hjust = 0)) +
  ggtitle("Relative \ncontribution") 

R2Plot_snp

# ===== FIGURE 5
snpCombine <- snpPlotP1 + snpPlotP2 + snpPlotP3  + snpPlotP4 + R2Plot_snp + plot_layout(ncol = 5, widths = c(1.2, 1, 1,1, 0.3))
saveFigure(fileName=paste0(HOME,"/results/figures/geneAgeBias"), plotName=snpCombine, w=48, h=40 )


###################################################################################
# ======================= DIRECTION OF INTERACTION ================================
###################################################################################
# Direction inferred from visual inspection 
directionDF=readRDS( paste0(HOME,"/output/rds/direction.rds"))
directionDF$label_merge=paste0(directionDF$pheno, "_", directionDF$SNP)
directionL=recodeLabel(df=directionDF, mergeBy="pheno")
directionDF$label_clean =directionL$label_clean
directionDF$typeC=revalCOL(var= directionDF$type)

# Direction inferred from interaction effects
dirSumL=lapply(levels(as.factor(snpEffP1$snp_var)), function(x) getDirection(df=snpEffP1, snp=x))
dirSumA=do.call(rbind, dirSumL)
saveRDS(dirSumA, paste0(HOME,"/output/rds/dirSum.rds"))

# Plot intensification effects

snpIntensification=levels(droplevels(as.factor(subset(dirSumA, direction=="intensification")$label_merge)))
plotInteL=lapply(c("reaction_time_rs181003402", "insomnia_rs113851554"), function(x) getDirection(snp=x, dfD=directionDF, return="plot", topM=1, bottom=0, bg="blue4"))
plotInt=ggarrange(plotlist=plotInteL, nrow=1, ncol=2, common.legend = TRUE,   legend = "top")

# Plot attenuation effects
snpAtten=levels(droplevels(as.factor(subset(dirSumA, direction=="attenuation")$label_merge)))
plotAttenL=lapply(c("bmi_rs6751993", "body_fat_free_mass_rs374892365"), function(x) getDirection(snp=x, dfD=directionDF, return="plot", legend="none", bg="aquamarine4"))
plotAtten=ggarrange(plotlist=plotAttenL, nrow=1, ncol=2)

# Plot crossover effects
snpCrossover=levels(droplevels(as.factor(subset(dirSumA, direction=="crossover")$label_merge)))
plotCOL=lapply(c("falls_last_year_rs73351184", "pulse_rate_rs13175615"), function(x) getDirection(snp=x, dfD=directionDF, return="plot", legend="none", bg="lightblue3"))
plotCOc=ggarrange(plotlist=plotCOL, nrow=1, ncol=2) 

# Plot incosistent effects
snpDiscrepant=levels(droplevels(as.factor(subset(dirSumA, direction=="inconsistent")$label_merge)))
plotDisL=lapply(c("overall_health_rs112054368", "number_medication_rs55730499"), function(x) getDirection(snp=x, dfD=directionDF, return="plot", legend="none", bg="lightcoral"))
plotDisc=ggarrange(plotlist=plotDisL, nrow=1, ncol=2)


# Plot summary of direction of interactions for all phenotype=genotype associations
dirSumA$typeC=revalDIR(var= dirSumA$direction)
dirSumA$colour <- factor(dirSumA$direction, labels = c( "aquamarine4", "lightblue3", "lightcoral", "blue4"))
dirSumA$label_cleanS=paste0(dirSumA$label_clean, " (", dirSumA$SNP, ")")
dirSumS=subset(dirSumA, type_raw=="snp_age_lme")

dirPlot=ggplot(dirSumS, aes(x = as.factor(type), y = label_cleanS, colour = direction)) +
  geom_point(size = 3, shape=15) +
  facet_grid(rows = vars(pheno), scales = "free_y", space = "free_y") +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_colour_manual("",values = levels(droplevels(dirSumS$colour)), label=levels(droplevels(as.factor(dirSumS$typeC))) ) +
  theme_minimal() +
  theme(  legend.position = c(0, 1),         # position at top right
          legend.justification = c(0.5, 0),  
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 8), 
          axis.text.x =element_blank(),
          axis.ticks.x = element_blank(),
          strip.text.x = element_text(size = 10),
          plot.margin = margin(t = 2, r = 0, b = 0, l = 0.3, "cm"),
          strip.text.y =  element_blank(),  # vertical strip text for pheno
          panel.spacing = unit(0.1, "lines"), 
          panel.grid.major.x = element_blank()
  ) + labs( y = "",shape = "", title = "" ) +
  guides(color = guide_legend(nrow = 2), alpha = "none" ) 
dirPlot

# ===== FIGURE 4
plotDirEx=ggarrange(plotInt, plotAtten, plotCOc, plotDisc, nrow=4, common.legend = TRUE)
directionComb=ggarrange(dirPlot, plotDirEx, ncol=2, widths=c(0.6, 0.7), labels=c("A", "B"))
saveFigure(fileName=paste0(HOME,"/results/figures/directionPlot"), plotName=directionComb, w=25, h=30 )


# =================================================== SUPPLEMENT FIGURES ========================================================


#########################################################################################
# =================== SUPPLEMENT: POLYGENIC SCORE INTERACTIONS ==========================
#########################################################################################
snpPGS=readRDS( paste0(HOME,"/output/rds/snpPGS.rds"))
snpPGS$r2_perd=round(snpPGS$r2*100, 2)
snpPGS$uCI=snpPGS$beta + 1.96 * snpPGS$se
snpPGS$lCI=snpPGS$beta - 1.96 * snpPGS$se
snpPGS$alpha=ifelse(snpPGS$P < 0.05, 1, 0.5)
snpPGSL=recodeLabel(df=snpPGS, mergeBy="var")
snpPGS$label_clean=snpPGSL$label_clean
snpPGS$label=paste0(snpPGS$SNP, "_", snpPGS$var, "_", snpPGS$type)
snpPGS$snp_var=paste0(snpPGS$SNP, "_", snpPGS$var)
nSNPsFDR=length(levels(as.factor(paste0(snpPGS$SNP, "_", snpPGS$var))))
snpPGS$p_fdr= do.call(rbind, lapply(snpPGS$P, function(x) p.adjust(x, "bonferroni",  nSNPsFDR )))
snpPGS$type <- factor(snpPGS$type, levels = c("snp_interaction_between_lm","snp_age2_lme", "age_snp_weighted", "snp_interaction_between_FU_lm","snp_cohort", "snp_age_lme"))
snpPGS$typeC=revalCOL(var=as.factor(snpPGS$type))
snpPGS$colour <- factor(snpPGS$typeC, labels =  c( "steelblue3", "darkorange3", "darkblue", "darkgreen","cyan4", "palevioletred4"))
mOrderG=data.frame(var=levels(as.factor(snpPGS$var)), ID=seq(1,NROW(levels(as.factor(snpPGS$var))), 1))
mOrderGL=recodeLabel(df=mOrderG, mergeBy="var")
mOrderG$label_clean =mOrderGL$label_clean

# Plot 1: Gene-by-age effects (cross-sectional versus longitudinal estimates)
snpPGSP1=merge(mOrderG, subset(snpPGS, type %in% c("snp_age_lme", "snp_interaction_between_lm")), by="var", all.x=TRUE, suffixes = c("", "_rem"))
mergeVar=subset(snpPGSP1, select=c(ID, var, label_clean))

pgsPlotP1 <- ggplot(snpPGSP1, aes(x = beta, y =fct_reorder(label_clean, -ID), colour=typeC, alpha=alpha, xmin = lCI, xmax = uCI)) +
  geom_pointrange(aes(xmin = lCI, xmax = uCI), size = 1) +
  scale_colour_manual("",values = levels(droplevels(snpPGSP1$colour)), label=levels(droplevels(snpPGSP1$typeC)) ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(x = expression(paste( delta[CS], " / ",   delta[L])), y = "Phenotype", color = "Variable", y = NULL) +
  theme_minimal() +
  themeAll +
  guides(color = guide_legend(nrow = 2), alpha = "none" ) +
  ggtitle("Model B1/B2")  +
  labs(subtitle = "Age-varying polygenetic effects") 

# Plot 2: Cohort-varying genetic effects
snpPGSP2=merge(mergeVar, subset(snpPGS, type %in% c("snp_cohort")), by="var", all.x=TRUE, suffixes = c("", "_rem"))
pgsPlotP2 <- ggplot(snpPGSP2, aes(x = beta, y =fct_reorder(label_clean, -ID), colour=typeC, alpha=alpha, xmin = lCI, xmax = uCI)) +
  geom_pointrange(size = 1) +
  scale_colour_manual("",values = levels(droplevels(snpPGSP2$colour)), label=levels(droplevels(snpPGSP2$typeC)) ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") + 
  labs(x = expression(paste( delta[C] )), y = "", color = "") +
  themeAge +
  themeAll+
  guides(color = guide_legend(nrow = 2), alpha = "none" ) +
  ggtitle("Model B3")   +
  labs(subtitle = "Year of birth-varying polygenetic effects") 
pgsPlotP2

# Plot 3: Non-linear age-varying genetic effects
snpPGSP3=merge(mergeVar, subset(snpPGS, type %in% c("snp_age2_lme")), by="var", all.x=TRUE, suffixes = c("", "_rem"))
pgsPlotP3 <- ggplot(snpPGSP3, aes(x = beta, y =fct_reorder(label_clean, -ID), colour=typeC, alpha=alpha, xmin = lCI, xmax = uCI)) +
  geom_pointrange(size = 1) +
  scale_colour_manual("",values = levels(droplevels(snpPGSP3$colour)), label=levels(droplevels(snpPGSP3$typeC)) ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") + 
  labs(x = expression(paste( delta[Q] )), y = "Phenotype", color = "Variable", , y = NULL) +
  themeAge +
  themeAll+
  guides(color = guide_legend(nrow = 2), alpha = "none" ) +
  ggtitle("Model B4")  +
  labs(subtitle = "Non-linear age-varying polygenetic effects") 
pgsPlotP3

# Plot 4: Weighted cross-sectional age effects
snpPGSP4=merge(mergeVar, subset(snpPGS, type %in% c("snp_interaction_between_lm", "age_snp_weighted", "snp_interaction_between_FU_lm")), by="var", all.x=TRUE, suffixes = c("", "_rem"))
pgsPlotP4 <- ggplot(snpPGSP4, aes(x = beta, y =fct_reorder(label_clean, -ID), colour=typeC, alpha=alpha, xmin = lCI, xmax = uCI)) +
  geom_pointrange(size = 1) +
  scale_colour_manual("",values = levels(droplevels(snpPGSP4$colour)), label=levels(droplevels(snpPGSP4$typeC)) ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  
  labs(x = expression(paste( delta[W], " / ", delta[CS] )), y = "Phenotype", color = "Variable",  y = NULL) +
  themeAge +
  themeAll+
  guides(color = guide_legend(nrow = 3), alpha = "none") +
  ggtitle("Model B5") +
  labs(subtitle = "Participation effects")
pgsPlotP4

# Explained variance in effect size differences
snpPGS$vi=snpPGS$se^2 # get variance
d1=subset(snpPGS, type %in% c("age_snp_weighted"), select=c(beta, vi, var)); d1$beta_weighted=d1$beta;d1$beta=NULL;d1$vi_weighted=d1$vi;d1$vi=NULL
d2=subset(snpPGS, type %in% c("snp_age_lme"), select=c(beta, vi, var)); d2$beta_within=d2$beta;d2$beta=NULL;d2$vi_within=d2$vi;d2$vi=NULL
d3=subset(snpPGS, type %in% c("snp_interaction_between_lm"), select=c(beta, vi, var)); d3$beta_between=d3$beta;d3$beta=NULL;d3$vi_between=d3$vi;d3$vi=NULL
d4=subset(snpPGS, type %in% c("snp_cohort"), select=c(beta, vi, var)); d4$beta_cohort=d4$beta;d4$beta=NULL;d4$vi_cohort=d4$vi;d4$vi=NULL
d5=subset(snpPGS, type %in% c("snp_interaction_between_FU_lm"), select=c(beta, vi, var)); d5$beta_between_FU=d5$beta;d5$beta=NULL;d5$vi_between_FU=d5$vi;d5$vi=NULL
d6=subset(snpPGS, type %in% c("snp_age2_lme"), select=c(beta, vi, var)); d6$beta_age2=d6$beta;d6$beta=NULL;d6$vi_age2=d6$vi;d6$vi=NULL

# Combine
pgsComb <- Reduce(function(x, y) merge(x, y, by = "var", all = TRUE), list(d1, d2, d3, d4, d5, d6))
pgsComb$beta_diff=pgsComb$beta_within-pgsComb$beta_between
pgsComb$beta_diff_fu=pgsComb$beta_weighted-pgsComb$beta_between_FU
pgsComb$vi_diff_fu    <-  pgsComb$vi_weighted + pgsComb$vi_between_FU - 2 * corW * sqrt(pgsComb$vi_weighted) * sqrt(pgsComb$vi_between_FU)


# R2 per predictor
m3_pgs=lm(beta_diff ~ beta_age2 + beta_cohort + beta_diff_fu, data=pgsComb)
as.data.frame(coef(summary(m3_pgs)))
rel_imp_pgs <- relaimpo::calc.relimp(m3_pgs, type = "lmg") # Lindeman, Merenda & Gold method
dfR2pgs=data.frame(r2=rel_imp_pgs@lmg)
totVarPGS=round(sum(dfR2pgs$r2)*100, 2)
dfR2pgs$predictor=rownames(dfR2pgs)
dfR2pgs$labelT=paste0(round(dfR2pgs$r2 * 100, 1), "%")
r2AddR <- data.frame( r2 = (100-totVarPGS)/100, predictor = "missing_var", labelT="Missing")
dfR2pgs <- rbind(dfR2pgs, r2AddR)
dfR2pgs$predictor <- factor(dfR2pgs$predictor, levels = c("missing_var", "beta_cohort", "beta_age2","beta_diff_w", "beta_diff_fu"))
dfR2pgs$predictorC=as.factor(revalCOL(var=dfR2pgs$predictor))
dfR2pgs$group=1
levels(dfR2pgs$predictorC)
colR2=c("lightgrey", "cyan4", "darkorange2", "darkblue", "darkgreen")

mDiscPGSSTD <- lm.beta::lm.beta(m3_pgs)
discPGSDF=as.data.frame(coef(summary(mDiscPGSSTD)))
discPDF=data.frame(predictor=rownames(discPGSDF), beta_std=discPGSDF$Standardized, P=discPGSDF$`Pr(>|t|)`)
dfR2PGSm=merge(dfR2pgs, discPDF, by="predictor", all.x=TRUE)
saveRDS(dfR2PGSm, paste0(HOME, "/output/rds/expR2pgs.rds"))

R2Plot_pgs <- ggplot(dfR2pgs, aes(x = group, y = r2 * 100, fill = predictorC)) +
  geom_bar(stat = "identity") +
  labs(title = "",
       x = "",
       y = "Explained variance (%)") +
  theme_minimal() +
  scale_fill_manual("", values = colR2) +
  geom_text(aes(label = labelT),
            position = position_stack(vjust = 0.5),
            size = 5,
            color = "white") +
  scale_y_continuous(position = "right") +  # Moves axis labels, ticks, and title to the right
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major.x = element_blank(), # Remove vertical grid lines
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.y = element_text(hjust = 0)) +
  ggtitle("Relative \ncontribution") 

R2Plot_pgs

# ===== SUPPLEMENT FIGURE
pgsCombine <- pgsPlotP1 + pgsPlotP2 + pgsPlotP3  + pgsPlotP4 + R2Plot_pgs + plot_layout(ncol = 5, widths = c(1.2, 1, 1,1, 0.3))
saveFigure(fileName=paste0(HOME,"/results/figures/pgsAgeBias"), plotName=pgsCombine, w=45, h=25 )


#########################################################################################
# ========================= SUPPLEMENT: LONGITUDINAL IPW  ===============================
#########################################################################################

# ======= Cross-sectional models
crosssecIPW=merge(mOrder, subset(ageEff, type %in% c("between_age", "age_weighted")), by="var", all.x=TRUE)
crossIPWplot <- ggplot(crosssecIPW, aes(x = beta, y =fct_reorder(label_clean, -ID), colour=typeC)) +
  geom_line(colour="black") +
  geom_pointrange(aes(xmin = lCI, xmax = uCI), size = 1) +
  scale_colour_manual("",values = levels(droplevels(crosssecIPW$colour)), label=levels(droplevels(crosssecIPW$typeC)) ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") + 
  labs(x = expression(paste( beta[W], " / ", beta )), y = "Phenotype", color = "Variable", y = NULL) +
  theme_minimal() +
  themeAll +
  guides(color = guide_legend(nrow = 2)) +
  ggtitle("A")  +
  labs(subtitle = "Participation effects (cross-sectional)") 

# ======= Longitudinal models
longIPW=merge(mOrder, subset(ageEff, type %in% c("within_age_complete_weights", "within_age_weighted")), by="var", all.x=TRUE) 
# note: within_age => all individuals with longitudinal data; within_age_complete_weights = individuals with both longitudinal data and sampling weigths (used here to ensure comparability)
longIPWplot <- ggplot(longIPW, aes(x = beta, y =fct_reorder(label_clean, -ID), colour=typeC)) +
  geom_line(colour="black") +
  geom_pointrange(aes(xmin = lCI, xmax = uCI), size = 1) +
  scale_colour_manual("",values = levels(droplevels(longIPW$colour)), label=levels(droplevels(longIPW$typeC)) ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Reference line at beta = 0
  labs(x = expression(paste( beta[W], " / ", beta )), y = "Phenotype", color = "Variable",  y = NULL) +
  themeAge +
  themeAll +
  guides(color = guide_legend(nrow = 3)) +
  ggtitle("B") +
  labs(subtitle = "Participation effects (longitudinal)") 
longIPWplot

IPW_diff <- subset(ageEff, type %in% c( "between_age", "age_weighted", "within_age", "within_age_weighted")) %>%
  dplyr::select(var, type, beta) %>%
  tidyr::pivot_wider(names_from = type, values_from = beta) %>%
  dplyr::mutate(diff_long =
                  within_age - within_age_weighted) %>%
  dplyr::mutate(diff_cross =
                  between_age - age_weighted)

IPW_diffL=data.frame(type=c(rep("diff_long", NROW(IPW_diff)), rep("diff_cross", NROW(IPW_diff))),
                     beta=c(IPW_diff$diff_long, IPW_diff$diff_cross),
                     var=c(IPW_diff$var, IPW_diff$var))
IPW_diffo=merge(mOrder, IPW_diffL, by="var", all.x=TRUE)
IPW_diffo$typeC=revalCOL(var=IPW_diffo$type)
saveRDS(IPW_diffo,  paste0(HOME,"/output/rds/ageW.rds"))


diffIPWplot <- ggplot(
  IPW_diffo, aes( x = beta, y = fct_reorder(label_clean, -ID), fill = type, colour = type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_col(position = position_dodge(width = 0.75),  linewidth = 0.9 ) +
  scale_colour_manual(values = c(diff_long  = "#470a0a", diff_cross = "darkblue")) +
  scale_fill_manual("", values = c(diff_long  = "#bc4f5e", diff_cross = "steelblue3"),   labels = c(
    diff_long  = "Longitudinal model",
    diff_cross = "Cross-sectional model"
  )) +
  labs(x = expression(paste(beta, " - ", beta[W])), y = NULL, colour = "Type") +
  themeAge +
  themeAll +
  ggtitle("C") +
  labs(subtitle = "Differences between unweighted and weighted effect estimates") + guides( colour = "none" )
diffIPWplot

# ======= COMBINED FIGURE
IPWplot <- crossIPWplot + longIPWplot + diffIPWplot+ plot_layout(ncol = 3, widths = c(1.2, 1.2, 1.2))
saveFigure(fileName=paste0(HOME,"/results/figures/IPWplot"), plotName=IPWplot, w=40, h=25 )


###################################################################################
# =================== Plot participation differences ==============================
###################################################################################
partBaseline=subset(snpEffP4, select=c(var_snp, beta, P), type=="snp_interaction_between_lm")
names(partBaseline)[names(partBaseline) == 'beta'] <- 'beta_B'
names(partBaseline)[names(partBaseline) == 'P'] <- 'P_B'
partFU=subset(snpEffP4, select=c(var_snp, beta, P), type=="snp_interaction_between_FU_lm")
names(partFU)[names(partFU) == 'beta'] <- 'beta_FU'
names(partFU)[names(partFU) == 'P'] <- 'P_FU'
partBaselineW=subset(snpEffP4, select=c(var_snp, beta, P), type=="age_snp_weighted")
names(partBaselineW)[names(partBaselineW) == 'beta'] <- 'beta_W'
names(partBaselineW)[names(partBaselineW) == 'P'] <- 'P_W'
partDF <- Reduce(function(x, y) merge(x, y, by = "var_snp"), 
                 list(partBaseline, partFU, partBaselineW))
partU <- partDF %>% distinct(var_snp, .keep_all = TRUE)

# Generate plots
scatterPartB=scatterPart(df=partU, x="W", y="B", colLine="aquamarine4")
scatterPartB_FU=scatterPart(df=partU, x="B", y="FU", colLine="gold4")
scatterPartW_FU=scatterPart(df=partU, x="W", y="FU", colLine="forestgreen")
scatterPartSum=ggarrange(scatterPartB, scatterPartB_FU, scatterPartW_FU, ncol=3, labels = c("A", "B", "C"))
saveFigure(fileName=paste0(HOME,"/results/figures/partPlot"), plotName=scatterPartSum, w=30, h=10 )


#########################################################################################
# ============================ Mendelian Randomization =================================
#########################################################################################

# Mendelian Randomization Results
library(TwoSampleMR)
mr=readRDS( paste0(HOME,"/output/rds/mr.rds"))
mrEXP=recodeLabel(df=mr, mergeBy="exposure")
mr$exp_clean=mrEXP$label_clean
mrOUT=recodeLabel(df=mr, mergeBy="outcome")
mr$out_clean=mrOUT$label_clean
mr$path=paste0(mr$exp_clean, " --> ", mr$out_clean, " (", mr$model, ")")
mr$model_c=revalMR(mr$model)
mrCol=c("darkgrey", "steelblue3", "palevioletred4")



plotMR=function(exp, out, m){
  if(m=="change"){
    c="palevioletred4"
  }
  if(m=="0_interaction"){
    c="steelblue3"
  }
  if(m=="0"){
    c="darkgrey"
  }
  
  mr_results=subset(mr, exposure==exp & outcome==out & model==m)
  dat=subset(mrDat, exposure==exp & outcome==out & model==m)
  dat$model_c=revalMR(dat$model)
  l=dat$model_c[1]
  
  mrRes=ggplot(mr_results, aes(x = b, y = model_c, colour=model_c)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = lCI, xmax = uCI), height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_label(
      aes(label = sprintf("%.2f", b)),
      nudge_y = 0.1,
      size = 3,
      label.size = 0.2,
      show.legend = FALSE
    )  +
    theme_minimal() +
    xlim(0, 0.5) +
    scale_colour_manual("", values = c) +
    theme(legend.position = "none",
          plot.margin = margin(t = 0, r = 1, b = 0, l = 0, "cm"),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank())
  #p=mr_scatter_plot(mr_results, dat)
  #print(p)
  mrRes
  combos <- unique(dat[, c("id.exposure", "id.outcome")])
  
  mrres <- lapply(seq_len(nrow(combos)), function(i) {
    d <- dat[dat$id.exposure == combos$id.exposure[i] & dat$id.outcome == combos$id.outcome[i], ]
    d <- subset(d, mr_keep)
    index <- d$beta.exposure < 0
    d$beta.exposure[index] <- d$beta.exposure[index] * -1
    d$beta.outcome[index] <- d$beta.outcome[index] * -1
    mrres <- subset(
      mr_results,
      id.exposure == d$id.exposure[1] & id.outcome == d$id.outcome[1]
    )
    mrres$a <- 0
    
    ggplot2::ggplot( data = d, aes(x = beta.exposure, y = beta.outcome)) +
      geom_errorbar(aes(ymin = beta.outcome - se.outcome,
                        ymax = beta.outcome + se.outcome),
                    colour = "grey",
                    width = 0) +
      geom_errorbar(aes(xmin = beta.exposure - se.exposure,
                        xmax = beta.exposure + se.exposure),
                    colour = "grey",
                    width = 0,
                    orientation = "y") +
      geom_point(colour = c) +
      geom_abline(data = mrres, colour = c, size=1, aes(intercept = a, slope = b), show.legend = TRUE) +
      labs(x = paste("Exposure (BMI)"), y = paste("Outcome (SBP)")) + theme_minimal() +
      theme(legend.direction = "vertical",
            legend.position = "none",
            strip.text.x = element_blank()) + guides(colour = ggplot2::guide_legend(ncol = 2))
  })
  
  mrOut=ggarrange(mrRes, mrres[[1]], ncol=2, widths = c(1, 1.2))
  
  mrOut = annotate_figure(
    mrOut,
    top = text_grob(l, size = 12)
  )
  
  print(mrOut)
  return(mrOut)
}

mrDat=readRDS( paste0(HOME,"/output/rds/mrDat.rds"))
mrP1=plotMR(exp="bmi", out="SBP", m="change")
mrP2=plotMR(exp="bmi", out="SBP", m="0_interaction")
mrP3=plotMR(exp="bmi", out="SBP", m="0")

mrComb <- ggarrange(
  mrP1,
  mrP2,
  mrP3,
  nrow = 3,
  ncol = 1
)

saveFigure(fileName=paste0(HOME,"/results/figures/mrPlot"), plotName=mrComb, w=17, h=20 )



# =============================================================
# ========== Scatter plot (consistency of effects) ============
# =============================================================
m1=subset(snpEff, type=="snp_interaction_between_lm", select=c(var,SNP, beta, label_clean, snp_var ))
m2=subset(snpEff, type=="snp_age_lme", select=c(beta, snp_var ))
clumpW=merge(m1, m2, by="snp_var",  suffixes = c("_b", "_w") )
clumpW$identification=ifelse(sign(clumpW$beta_b)==sign(clumpW$beta_w),"consistent", "inconsistent")
clumpW$beta_gene_diff <- abs(clumpW$beta_b - clumpW$beta_w)
top10 <- clumpW[order(-clumpW$beta_gene_diff), ][1:10, ]

# test if slope is significantly different from 0
corWB=cor.test(clumpW$beta_w, clumpW$beta_b)
mWB=lm(beta_w  ~ beta_b, data=clumpW) 
r2_pheno=summary(mWB)$r.squared
mBW=as.data.frame(coef(summary(mWB)))
intercept <- round(mBW$Estimate[1], 3)
slope <- round(mBW$Estimate[2], 2)
pval <- formatC(mBW$`Pr(>|t|)`[2], format = "e", digits = 2)  # Scientific notation

# test if slope is significantly different from 1
t_value <- (mBW$Estimate[2] - 1) / mBW$`Std. Error`[2]
df <- df.residual(mWB)
p_value_t <- 2 * pt(-abs(t_value), df)
slopeR=data.frame(slope=slope, p_0=pval, p_1=p_value_t, t_value=t_value, r=corWB$estimate, r2=r2_pheno)
saveRDS(slopeR, paste0(HOME, "/output/rds/slopeR.rds"))

minmax=0.007

scatterP <- ggplot(clumpW, aes(x = beta_b, y = beta_w, color = identification)) +
  scale_colour_manual("", values = c("#1b9e77", "lightcoral"), 
                      labels =  c("Consistent direction", "Inconsistent direction")) +
  labs(x = expression(paste(delta[CS])), y=expression(paste(delta[L]))) +
  xlim(-minmax, minmax) + ylim(-minmax, minmax) +
  geom_abline(intercept = intercept, slope = slope, colour="black", linetype=1) +
  geom_abline(intercept = 0, slope = 1, colour="grey") +
  geom_abline(intercept = 0, slope = 0, colour="grey") +
  geom_point(aes(alpha=0.3)) +
  ggrepel::geom_text_repel(data = top10, aes(label = label_clean), size = 3, box.padding = 0.5, max.overlaps = 10) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.text = element_text(size=13),
    plot.title = element_text(face = "bold", hjust = 0),
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5)  # Rotate label, adjust position
  ) +
  guides(alpha = "none", colour = guide_legend(override.aes = list(size = 6))) +
  annotate("text", 
           x = 0.004, y = mBW$Estimate[1] + 0.003 * slope, 
           label = paste0("alpha[1] == ", slope),
           parse = TRUE, hjust = 0, size = 4)

print(scatterP)
saveFigure(fileName=paste0(HOME,"/results/figures/scatterDir"), plotName=scatterP, w=17, h=13 )


###################################################################################
# ================== Compare all FU versus pre-COVID FU ===========================
###################################################################################

# ================== Phdenotypic analyses ======================
ageCOVID=readRDS( paste0(HOME,"/output/rds/ageCOVID.rds"))
ageCOVID$period=ifelse(grepl( "_c", ageCOVID$var)==T, "pre", "all")
ageCOVID$var=stringr::str_remove(ageCOVID$var, "_c")
ageCOVID_w <- ageCOVID %>%
  dplyr::filter(period %in% c("pre", "all")) %>%
  dplyr::select(var, beta, period) %>%
  pivot_wider(names_from = period, values_from = beta)
ageCOVIDC=recodeLabel(df=ageCOVID_w, mergeBy="var")
ageCOVID_w$label_clean =ageCOVIDC$label_clean

scatterCO_p=scatterCOVID(df=ageCOVID_w, pre="pre", all="all", t="A. Phenotypic analyses", minSens=-0.06, maxSens=0.06, colLine="darkred", letter="beta")

# ================== Genotypic analyses ======================
m1Sens=subset(clump, type=="interaction_within", select=c(label_wide, pheno,SNP, BETA, lCI, uCI, label_clean, p_fdr, P, pheno_snp, gene))
m2Sens=subset(clumpA, snppheno %in% sigSNP & type=="interaction_within_c", select=c(label_wide, BETA, lCI, uCI, p_fdr, P))
clumpSens=merge(m1Sens, m2Sens, by="label_wide",  suffixes = c("_all", "_pre") )
scatterCO_g=scatterCOVID(df=clumpSens, pre="BETA_pre", all="BETA_all", t="B. Genotypic analyses", minSens=-0.03, maxSens=0.03, colLine="firebrick2", letter="delta")

# Combine
scatterSens=ggarrange(scatterCO_p, scatterCO_g, ncol=2)
saveFigure(fileName=paste0(HOME,"/results/figures/scatterSens"), plotName=scatterSens, w=23, h=11 )

corSensP=cor.test(ageCOVID_w$all, ageCOVID_w$pre)
corSensG=cor.test(clumpSens$BETA_all, clumpSens$BETA_pre)
corSensPG=rbind(data.frame(type="geno", r=corSensG$estimate, p=corSensG$p.value),
                data.frame(type="pheno", r=corSensP$estimate, p=corSensP$p.value))
saveRDS(corSensPG, paste0(HOME,"/output/rds/corSens.rds"))


#########################################################################################
# ================ ILLUSTRATION NON-LINEAR AGE EFFECTS ==================================
#########################################################################################
age2Pos=plotAge2(dir="pos")
age2Neg=plotAge2(dir="neg")
biasSimP=ggarrange(age2Pos, age2Neg, widths=c(1, 1), labels = c("A", "B"), common.legend = TRUE)
saveFigure(fileName=paste0(HOME,"/results/figures/biasSimP"), plotName=biasSimP, w=22, h=12 )

# =======================================================================
# ================ ILLUSTRATTION  MANHATTAN PLOT ========================
# =======================================================================
# part of figure 2
num_chromosomes <- 22            
snps_per_chromosome <- 1000      
total_snps <- num_chromosomes * snps_per_chromosome

# Generate chromosome, position, and p-values
chr <- rep(1:num_chromosomes, each = snps_per_chromosome)
pos <- unlist(lapply(rep(snps_per_chromosome, num_chromosomes), function(x) sort(runif(x, 1, 1e6))))
pval <- runif(total_snps, 0, 1)

# Create data frame
data <- data.frame(CHR = chr, POS = pos, P = pval)
data$logP <- -log10(data$P)
data$CHRb <- ifelse(data$CHR %% 2 == 0, 1, 0)

# Plot
manH=ggplot(data, aes(x = POS, y = logP, color = as.factor(CHRb))) +
  geom_point(alpha = 0.6, size = 0.7) +
  scale_color_manual(values = c("lightgrey", "darkgrey")) +
  scale_x_continuous(labels = NULL) +
  labs(title = "", x = "", y = expression(-log[10](italic(P)))) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),   
    axis.text.x = element_blank(),    
    axis.text.y = element_blank(),     
    axis.ticks.x = element_blank(),   
    axis.ticks.y = element_blank(),     
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 8),
    panel.spacing = unit(0, "lines")
  ) +
  facet_wrap(~CHR, scales = "free_x", nrow = 1, strip.position = "bottom")
manH
saveFigure(fileName=paste0(HOME,"/results/figures/manH"), plotName=manH, h=10, w=13)


# =======================================================================
# ======================== REGRESSION DILUTION ==========================
# =======================================================================
# ============== Phenotype results ===================
cohortDil=getDil(df=ageComb, dfR2=dfR2, var="cohort")
age2Dil=getDil(df=ageComb, dfR2=dfR2, var="age2")
partDil=getDil(df=ageComb, dfR2=dfR2, var="diff_w_fu")
sumDilpheno=rbind(cohortDil, age2Dil, partDil)
sumDilpheno$r2_scaled <- sumDilpheno$r2_cor / sum(sumDilpheno$r2_cor)
sumDilpheno$r2_cor=sumDilpheno$r2_scaled 
sumDilpheno$predictor=paste0(sumDilpheno$predictor, "_pheno")
sumDilpheno$type="1_pheno"

# ============== Genotype results ===================
cohortDilgeno=getDil(df=snpComb, dfR2=dfR2Genem, var="cohort")
age2Dilgeno=getDil(df=snpComb, dfR2=dfR2Genem, var="age2")
partDilgeno=getDil(df=snpComb, dfR2=dfR2Genem, var="diff_fu")
sumDilgeno=rbind(cohortDilgeno, age2Dilgeno, partDilgeno)
sumDilgeno$predictor=paste0(sumDilgeno$predictor, "_geno")

sumDilgeno$r2_scaled <- sumDilgeno$r2_cor / sum(sumDilgeno$r2_cor)
sumDilgeno$r2_cor=sumDilgeno$r2_scaled 
sumDilgeno$type="2_geno"
sumDil=rbind(subset(sumDilpheno, select=c(predictor, r2, dillution, r2_cor, type)), subset(sumDilgeno, select=c(predictor, r2, dillution, r2_cor, type)))

saveRDS(sumDilgeno, paste0(HOME, "/output/rds/dilution.rds"))


# Combine
sumDil_long <- sumDil %>% pivot_longer(cols = c(r2, r2_cor),  names_to = "Correction",  values_to = "R2") %>%
  mutate( Correction = recode(Correction, r2 = "Uncorrected", r2_cor = "Corrected"), label = paste0("λ = ", round(dillution, 3)))

predictor_order <- c("beta_cohort_pheno", "beta_cohort_geno", "beta_diff_w_fu_pheno", "beta_diff_fu_geno", "beta_age2_pheno", "beta_age2_geno")

sumDil_long <- sumDil_long %>% mutate(type = forcats::fct_recode(type,  "Genotypic analyses" = "2_geno", "Phenotypic analyses" = "1_pheno"),
                                      type = forcats::fct_rev(type),
                                      Correction = forcats::fct_rev(Correction),
                                      predictor = factor(predictor, levels = predictor_order))

math_labels <- c(
  "beta_cohort_pheno"       = expression(paste("Cohort (", beta[C], ")")),
  "beta_cohort_geno"      = expression(paste("Cohort (", delta[C], ")")),
  "beta_age2_pheno"               = expression(paste("Non-linear age (", beta[Q], ")")),
  "beta_diff_w_fu_pheno"         = expression(paste("Participation (", beta[P], ")")),
  "beta_age2_geno"                = expression(paste("Non-linear age (", delta[Q], ")")),
  "beta_diff_fu_geno"             = expression(paste("Participation (", delta[P], ")"))
)


# Plot
dilutionPlot=ggplot(sumDil_long, aes(x = predictor, y = R2, fill = Correction)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.5) +
  facet_grid(cols = vars(type), scales = "free_x", space = "free_y") +
  geom_text(
    aes(label = ifelse(Correction == "Corrected", label, "")),
    position = position_dodge(width = 0.6),
    vjust = -0.5, size = 3.5
  ) +
  labs( x = "", y = expression(R^2), fill = "") +
  theme_minimal() +
  ylim(0,1) +
  scale_fill_manual( values = c("Corrected" = "#1f77b4", "Uncorrected" = "#ff7f0e")) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = math_labels, guide = guide_axis(angle = 45)) 
dilutionPlot

showtext::showtext_auto()  # Automatically use showtext for all text drawing
saveFigure(fileName=paste0(HOME,"/results/figures/dilutionPlot"), plotName=dilutionPlot, w=20, h=15 )


###################################################################################
# ================================== CREATE TABLES ================================
###################################################################################

source(paste0(HOME, "/analysis/createTablesAge.R"))








