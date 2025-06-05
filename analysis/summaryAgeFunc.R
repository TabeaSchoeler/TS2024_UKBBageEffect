#######################################################
# =================== LOAD LIBRARIES ==================
#######################################################
# remove all of the objects that are stored in the global environment
.libPaths(new=c("/Users/tabea/Dropbox/progs/R/library"))
#system('R_LIBS=/Users/tabea/Dropbox/progs/R/library')
# Load and install libraries
load.lib=c('ggplot2', 'lme4','sjPlot', 'sjmisc', 'sjlabelled', 'rsq', 'tidyr', 'MuMIn', 'tidyverse', 'reshape2', 'ggtext',
           'gridExtra', 'scales','grid', 'lattice', 'ggpubr', 'paletteer', 'pcaMethods', 'plyr', 'patchwork', 'latex2exp', 
           'cowplot', "otargen", "ggcorrplot", 'effectsize', 'viridis', 'parameters', 'ggrepel', 'gwasrapidd', 'biomaRt', 'english')

install.lib<-load.lib[!load.lib %in% installed.packages()]


# Install missing packages
for(lib in install.lib) install.packages(pkgs=lib, dependencies=TRUE)

# Load all packages
sapply(load.lib,require,character=TRUE)



print("Read in variable names")
variables <- as.data.frame(readxl::read_excel(paste0(HOME, "/data/variableAge.xlsx")))

recodeLabel=function(df, mergeBy="label"){
  label <- as.data.frame(readxl::read_excel(paste0(HOME, "/data/variableAge.xlsx")))
  lDF=data.frame(label=label$label, label_clean=label$label_clean, dimension=label$dimension)
  lDFC <- na.omit(lDF)
  lDFC=subset(lDFC, as.character(label)!="NA")
  dfm=distinct(lDFC, label, .keep_all = TRUE)
  
  df$ID=seq(1,NROW(df), 1)
  dfsel=merge(df, dfm, by.x=mergeBy, by.y="label", all.x=T, suffixes = c("", "_rem") )
  dfOut <- dfsel[order(dfsel$ID),] 
  dfOut$ID=NULL
  return(dfOut)
}

saveFigure=function(fileName, plotName, w, h){
  print("Save PDF")
  ggsave(file=paste0(fileName, ".pdf"), plot = plotName, width = w, height = h, units = "cm", bg='transparent', limitsize = FALSE)
  ggsave(file=paste0(fileName, ".svg"), plot = plotName, width = w, height = h, units = "cm", bg='transparent', limitsize = FALSE)
}


importFiles=function(file){
  print(file)
  file.copy(from=paste0(MDP, "/output/rds/", file, ".rds"), to=paste0(HOME, "/output/rds/"), overwrite = TRUE)
  df=readRDS(paste0(HOME, "/output/rds/", file, ".rds"))
  aDF=data.frame(date=attributes(df)$timestamp[1])
  df1=readRDS(paste0(MDP, "/output/rds/", file, ".rds"))
  aDF1=data.frame(date=attributes(df1)$timestamp[1])
  print( paste0(file, ", creation date file on mountainduck: ", aDF1$date))
  print( paste0(file, ", creation date copied file: ", aDF$date))
}



openTarget=function(df){
  snps=levels(as.factor(df$SNP))
  listG=list()
  
  for (i in 1:length(snps)) {
    print(snps[i])
    s <- snps[i]
    dfLabel <- paste0(subset(df, SNP == s)$label_clean, collapse = ", ")
    info <- NULL  # Define info outside tryCatch block
    
    tryCatch({
      info <- variantInfo(s)
      pos_GRCh37 = info$position
      pos=geneInfo(info$nearestGene.symbol)
      start_GRCh37=pos$start
      end_GRCh37=pos$end
    }, error = function(e) {
      cat("Error occurred for snp", s, ":", conditionMessage(e), "\n")
    })
    
    # Check if info is NULL before accessing its elements
    if (!is.null(info)) {
      listG[[i]] <- data.frame(SNP = s, gene = info$nearestGene.symbol, pos_GRCh37, start_GRCh37, end_GRCh37)
    } else {
      print("Check GWA catalog")
        info=get_associations(variant_id = s)
        
        if(NROW(info@genes$gene_name)>0){
          gene_common <- names(sort(table(na.omit(info@genes$gene_name)), decreasing = TRUE))[1]
          listG[[i]] <- data.frame(SNP = s, gene = gene_common, pos_GRCh37=NA, start_GRCh37=NA, end_GRCh37=NA)
        } 
        if(NROW(info@genes$gene_name)==0){
          print("Check biomaRt")
          snp_mart <- useEnsembl( biomart = "snp",  dataset = "hsapiens_snp")
          snp_info <- getBM( attributes = c("refsnp_id", "chr_name", "chrom_start", "allele", "associated_gene"),
            filters = "snp_filter",
            values = s, mart = snp_mart)
          
          associated_genes=snp_info$associated_gene
          gene_entries <- associated_genes[associated_genes != "" & !grepl(",", associated_genes)][1]
          
          if(NROW(gene_entries)>0){
            listG[[i]] <- data.frame(SNP = s, gene = gene_entries, pos_GRCh37=NA, start_GRCh37=NA, end_GRCh37=NA)
          }
        }
    }
    
  }
  return(do.call(rbind, listG))
  
}



nonzero_decimal <- function(x) {
  if (x == 0) return("0")
  s <- sprintf("%.10f", x)
  s_trim <- sub("^0\\.", "", s)
  nz_pos <- regexpr("[1-9]", s_trim)
  decimals <- as.integer(nz_pos)
  x_round <- round(x, decimals)
  out <- format(x_round, nsmall=decimals, scientific=FALSE, trim=TRUE)
  out <- sub("0+$", "", sub("\\.$", "", out))
  out
}


revalCOL=function(var){
  variableOutF=as.factor(
    revalue( as.factor(var), 
             c("within"="explained % (within-subject)", 
               "between"="explained % (between-subject)", 
               "within_age" = "Longitudinal age effects (prospective sample)",
               "age_snp_within" = "Longitudinal age-varying genetic effects (prospective sample)",
               "snp_age_lme" = "Longitudinal age-varying genetic effects (prospective sample)",
               "age_snp_between" = "Cross-sectional age-varying genetic effects (baseline sample)",
               "interaction_within" = "Longitudinal",
               "interaction_between" = "Cross-sectional",
               "age_snp_weighted" = "Cross-sectional age-varying genetic effects (weighted baseline sample)",
               "snp_age2_lme" = "Non-linear age-varying genetic effects",
               "snp_cohort" = "Cohort-varying genetic effects",
               "snp_cohort2" = "Non-linear cohort-varying genetic effects",
               "within_age2_effect" = "Non-linear effect (age/cohort)",
               "non_linear_age"  = "Non-linear age effects",
               "beta_age2" = "Non-linear age effects",
               "beta_cohort2" = "Non-linear effects (cohort)",
               "beta_cohort" = "Cohort effect", 
               "cohort" = "Cohort effect", 
               "cohort2" = "Non-linear cohort effect", 
               "beta_diff_w" = "Selective participation",
               "between_age" = "Cross-sectional age effects (baseline sample)",
               "between_age_FU" = "Cross-sectional age effects (prospective sample)",
               "snp_interaction_between_lm" = "Cross-sectional age-varying genetic effects (baseline sample)",
               "age_weighted" ="Cross-sectional age effects (weighted baseline sample)",
               "snp_interaction_between_FU_lm" ="Cross-sectional age-varying genetic effects (prospective sample)",
               "between_age_weighted" = "Cross-sectional age effects (weighted)",
               "between_adjusted_age2" = "Cross-sectional age effects (age^2 adjusted)",
               "between_age2_effect" = "Age^2 effects (cross-sectional data)"))) 
  return(variableOutF)
 
}



plotParticipation=function(df){
  slope1 <- lm(Y ~ X, data = df)
  slope2 <- lm(Y ~ X, data = subset(df, inclusion == 1))
  slope1_intercept <- coef(slope1)[1]
  slope1_slope <- coef(slope1)[2]
  slope2_intercept <- coef(slope2)[1]
  slope2_slope <- coef(slope2)[2]
  
  dfR=sample_n(df, 100)
  dfS=sample_n(subset(df, inclusion == 1), 20)
  
  NROW(subset(df, inclusion == 1))
  
  cR <- "#20A39E"   # Teal for full data regression line
  cS <- "#EF5B5B"   # Coral Red for subset regression line
  
  #labelR="All individuals"
  #labelS="Selected sample"
  # Plot
  p <- ggplot(df, aes(x = X, y = Y)) +
    geom_abline(aes(color = "All individuals"), 
                intercept = slope1_intercept, slope = slope1_slope, linetype = "solid", size = 1) +
    geom_abline(aes(color = "Selected sample"), 
                intercept = slope2_intercept, slope = slope2_slope, linetype = "solid", size = 1) +
    geom_point(data = dfR, aes(color = "All individuals"), alpha = 0.6) +  
    geom_point(data = dfS, aes(color = "Selected sample"), alpha = 0.6, size = 3) +  
    scale_color_manual(name = "", 
                       values = c("All individuals" = cS, "Selected sample" = cR)) +
    labs(ttitle = "Age-differential participation",
         x = "Age (years)",
         y = "") +
    ylim(0, 200) +
    theme_minimal() +
    theme(legend.position = c(0.1, 0.1),  # Adjust x, y position (bottom-left)
          legend.justification = c(0, 0),  # Anchor legend at bottom-left
          legend.background = element_rect(fill = alpha("white", 0.5), color = NA),  # Optional background for readability
          legend.box.margin = margin(5, 5, 5, 5),
          plot.margin = unit(c(top=1, right=0, 0, 0), "cm"))
  
  # Print plot
  print(p)
  
  return(p)
}




plotAge2=function(dir="neg" ){

  # Generate age sequence
  minAgeCS <- 40
  maxAgeCS <- 60
  minAgeL <- 55
  maxAgeL <- 70
  age <- seq(40, 70, by = 0.01)
  age_c <- age - mean(age)  # center around 75

  age2Eff=0.005
  # Adjusted coefficients
  ageEff <- -0.2
  if(dir=="neg"){
    age2Eff <- -1*age2Eff
    label="Negative quadratic effects"
  }
  if(dir=="pos"){
    label="Positive quadratic effects"
  }

  interceptHeight <- 170  # reasonable average height
  
  # Simulated outcomes
  outcome <- interceptHeight + ageEff * age_c + age2Eff * age_c^2
  
  dfN <- data.frame(outcome = outcome, age = age)
  dfNLong <- subset(dfN, age >= minAgeL & age <= maxAgeL)
  dfNCS <- subset(dfN, age >= minAgeCS & age <= maxAgeCS)
  dfComb <- data.frame(age = age, outcome = outcome, type = label)
  
  # Colors for illustration
  colIll <- c("steelblue3", "palevioletred4", "darkgrey")
  
  # Calculate slopes and intercepts
  linesN_L <- data.frame(slope = getSlope(df = dfNLong)$slope, intercept = getSlope(df = dfNLong)$intercept, model = "Longitudinal", min = minAgeL, max = maxAgeL)
  linesN_C <- data.frame(slope = getSlope(df = dfNCS)$slope, intercept = getSlope(df = dfNCS)$intercept, model = "Cross-sectional", min = minAgeCS, max = maxAgeCS)
  linesDF <- rbind(linesN_L, linesN_C)
  
  # Plot with ggplot2
  p <- ggplot(dfComb, aes(x = age, y = outcome)) +
    scale_colour_manual("", values = colIll) +
    geom_line(size = 5, linetype = "solid", alpha = 0.3) +
    geom_segment(data = linesDF, aes(x = min, xend = max, 
                                     y = intercept + slope * min, 
                                     yend = intercept + slope * max,
                                     colour = model), size = 2) +
    labs(x = "Age", y = "Height (cm)", subtitle = label) +
    theme_minimal() +
    theme(legend.position = "top",
          plot.margin = margin(t = 0.9, r = 0, b = 0.1, l = 0, "cm")) +
    guides(color = guide_legend(nrow = 1))
  
  print(p)
  
  return( p)
  
}


themeAll=theme(axis.title.x = element_text(size = 15),
               legend.position="top",
               plot.title=element_text(face="bold", hjust=0.5),
               strip.text.y = element_blank(),
               plot.subtitle=element_text( hjust=0.5),
               axis.title.y = element_blank()) 

plotRES=function(df, c, lim, title, hline){
  
  p=ggplot(df, aes(x = param, y = beta, colour=model, shape=type)) +
    geom_point(size=4, alpha=0.5) +  # Scatter points
    scale_colour_manual("",values = c ) + # guide = "none" labels=c("Age^2", "Cohort effects") 
    scale_shape_manual("", values = c( 16, 17), labels=c("Cross-sectional", "Longitudinal") ) + # c
    labs(title = "",
         x = "",
         y = "") + 
    theme_minimal() +
    theme(legend.position="top") + #          axis.text.x = element_text(angle = 45, hjust = 1)
    geom_hline(yintercept = 0, linetype = "dashed", color = "lightgrey", size = 1) +
    geom_hline(yintercept = hline, color = "black", size = 1, alpha=0.5) +
    #facet_grid(vars(type) , scales="free", space = "free") +
    labs(title="", # title
         x = "",
         y = "Estimated age effects") +
    scale_x_reverse() 
  print(p)
  return(p)
}


themeAge=theme_minimal() +
  theme(legend.position="top",
        plot.title=element_text(face="bold", hjust=0.5),
        plot.subtitle=element_text( hjust=0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),  # Hide y-axis text for the right plot
        axis.ticks.y = element_blank(),
        strip.text.y = element_blank()) 



getDirection=function(df=NULL, snp=NULL, return="data", legend="top", topM=1, bg, bottom=1, dfD=NULL){
  print(snp)
  
  if(return=="data"){
  dfC=readRDS( paste0(HOME,"/output/rds/clump.rds"))
  dfC$var_snp=paste0(dfC$SNP, "_", dfC$pheno)
  df$var_snp=paste0(df$SNP, "_", df$var)
  dfI=subset(df, var_snp==snp)
  dfM=subset(dfC, var_snp==snp & type=="marginal")
  
  if(dfM$p_fdr>0.05){
    dfI$direction="crossover"
  } else {
    dfI$direction=ifelse(sign(dfI$beta)==sign(dfM$BETA), "intensification", "attenuation")
  }
  dfI$direction=ifelse(sign(subset(dfI, type=="snp_interaction_between_lm")$beta)!=sign(subset(dfI, type=="snp_age_lme")$beta), "inconsistent",     dfI$direction)
  dfOut=data.frame(pheno=dfI$var, SNP=dfI$SNP[1], label_merge=snp, label_clean=dfI$label_clean, type=dfI$typeC, type_raw=dfI$type, direction=dfI$direction)
}
  
  if(return=="plot"){
    dfS=subset(dfD, label_merge==snp)
    dfS$label=paste0(dfS$label_clean, "\n(", dfS$SNP, ")")
    var=paste0( dfS$label[1])
    colM=c( "steelblue3", "palevioletred4")
    
    dfOut=ggplot(dfS, aes(x = medianAge, y = BETA, color = typeC, group=typeC)) +
      geom_line(size=1) +
      geom_point(size = 4) +
      theme_minimal() +
      scale_colour_manual("",values = colM, label=levels(droplevels(as.factor(dfD$typeC))) ) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", size=0.3) +  # Reference line at beta = 0
      labs(x = "Age", y = expression(paste( italic(beta) )), title = var) +
      guides(color = guide_legend(nrow = 1)) +
      theme(legend.position=legend,
            plot.title = element_text(size = 8),   
            panel.grid = element_blank(),
            panel.background = element_rect(fill = scales::alpha(bg, 0.3), color = NA),  # more transparent
            plot.background = element_rect(fill = scales::alpha("white", 0), color = NA),   # <- This sets the plot background
            plot.margin = margin(t = topM, r = 0, b = bottom, l = 0.3, "cm"),
            axis.title.y = element_text(angle = 0, vjust = 0.5) )
    print(dfOut)
  }
  return(dfOut)
}

revalDIR=function(var){
  variableOutF=as.factor(
    revalue( as.factor(var), 
             c("intensification"="Intensification", 
               "attenuation"="Attenuation", 
               "inconsistent"="Inconsistent", 
               "crossover" = "Crossover")))
}





scatterCOVID=function(df, pre, all, minSens, maxSens, t, colLine, letter){
  head(df)
  dfI=data.frame(preC=df[[pre]], allC=df[[all]], label=df$label_clean)
  mC=lm(preC  ~ allC, data=dfI)
  mSENS=as.data.frame(coef( summary(mC)))
  
  dfI$betaDiff <- abs(dfI$preC - dfI$allC)
  top10 <- dfI[order(-dfI$betaDiff), ][1:5, ]
  
  interceptSens <- round(mSENS$Estimate[1], 3)
  slopeSens <- round(mSENS$Estimate[2], 2)
  pvalSens <- formatC(mSENS$`Pr(>|t|)`[2], format = "e", digits = 2)  # Scientific notation
  corSens=cor.test(dfI$preC, dfI$allC)
  
  # test if slope is significantly different from 1
  t_value_t <- (mSENS$Estimate[2] - 1) / mSENS$`Std. Error`[2]
  df <- df.residual(mC)
  p_value_t <- 2 * pt(-abs(t_value_t), df)
  
  annot_text <- paste0("Slope = ", slopeSens,
                       "\nP = ", formatC(p_value_t, format = "e", digits = 1), " (t = ", round(t_value_t, 2), ")")
  
  x_subscript = "ALL"
  y_subscript = "pre-COVID"
  
  scatterSens <- ggplot(data=dfI, aes( x = allC, y = preC)) +
    labs(
      x = bquote(italic(.(as.name(letter))[.(x_subscript)])),
      y = bquote(italic(.(as.name(letter))[.(y_subscript)]))
    ) +
    ggrepel::geom_text_repel(data = top10, aes(label = label), size = 3, box.padding = 0.5, max.overlaps = 10) +
    xlim(minSens, maxSens) + ylim(minSens, maxSens) +
    geom_abline(intercept = interceptSens, slope = slopeSens, colour=colLine, linetype=1) +
    geom_abline(intercept = 0, slope = 1, colour="grey") +
    geom_abline(intercept = 0, slope = 0, colour="grey") +
    geom_point(aes(alpha=0.3)) +
    theme_minimal() +
    theme(legend.position = "top",
          plot.title = element_text(face = "bold", hjust = 0)) +
    guides( color = "none", alpha = "none") +
    geom_label(aes(x = minSens+0.005, y = maxSens - 0.005 * slopeSens, label = annot_text),
               hjust = 0, vjust = 1,
               size = 3.5,
               label.size = 0.2,
               fill = "white",
               color = colLine) +
    ggtitle(t)

  scatterSens
  return(scatterSens)
}


scatterPart=function(df, x, y, colLine){
  dfIn=data.frame(x=df[[paste0("beta_", x)]], y=df[[paste0("beta_", y)]], x_p=df[[paste0("P_", x)]], y_p=df[[paste0("P_", y)]], label=df$var_snp)
  dfS=subset(dfIn, x_p<0.05 | y_p <0.05)
  head(dfIn)
  
  mWB=lm(x  ~ y, data=dfS) 
  mBW=as.data.frame(coef(summary(mWB)))
  slope <- round(mBW$Estimate[2], 2)
  intercept <- mBW$Estimate[1]
  
  # test if slope is significantly different from 1
  t_value_t <- (mBW$Estimate[2] - 1) / mBW$`Std. Error`[2]
  df <- df.residual(mWB)
  p_value_t <- 2 * pt(-abs(t_value_t), df)
  
  annot_text <- paste0("Slope = ", slope,
                       "\nP = ", formatC(p_value_t, format = "e", digits = 1), " (t = ", round(t_value_t, 2), ")")
  
  minmax = ceiling(max(abs(dfS$x), abs(dfS$y)) * 1000) / 1000
  x_label=x
  y_label=y
  
  scatterP <- ggplot(dfS, aes(x = x, y = y)) +
    labs(
      x = bquote(delta[.(x_label)]),
      y = bquote(delta[.(y_label)])
    ) +
    xlim(-minmax, minmax) + ylim(-minmax, minmax) +
    geom_abline(intercept = intercept, slope = slope, colour=colLine, linetype=1, size=2) +
    geom_abline(intercept = 0, slope = 1, colour="grey") +
    geom_abline(intercept = 0, slope = 0, colour="grey") +
    geom_point(aes(alpha=0.3)) +
    theme_minimal() +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", hjust = 0),
      axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5)  # Rotate label, adjust position
    ) +
    guides(color = "none", alpha = "none") +
    geom_label(aes(x = -minmax * 0.95, y = minmax * 0.95, label = annot_text),
               hjust = 0, vjust = 1,
               size = 3.5,
               label.size = 0.2,
              fill = "white",
              color = colLine)


  print(scatterP)
  return(scatterP)
}

