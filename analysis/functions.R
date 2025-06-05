
print("Install R packages")
PROGSLOCAL <- file.path(dirname(HOME), "programs")
DATA <- file.path(dirname(OUT), "data")
PROGS <- file.path(dirname(OUT), "programs")
package="lme4_1.1-35.5.tar.gz"
#package="bigsparser_0.7.3.tar.gz"
package="bigstatsr_1.6.1.tar.gz"
package="svylme-master"

#install.packages(paste0(PROGSLOCAL, "/", package), lib=paste0(PROGS, "/R"),repos = NULL)


print("Load packages")
R_libPaths=paste0(PROGS, "/R")
.libPaths(R_libPaths)
library(stringr)

saveOutput=function(object, label, upload="no"){
  saveRDS(object , paste0(OUT,"/output/rds/",label,".rds"))
  if(upload=="yes"){
    saveRDS(object , paste0(HOME,"/output/rds/",label,".rds"))
  }
}



# ======================= extractPheno.R ============================


extractPheno=function(var, varFile, dfL, return){
   # Get data with fewest missing values
   print(var)

  print("Get phenotype data")
   varInfo=subset(varFile, label==var)
   nameDat=selectDat(listIn=dfL, ID=varInfo$ID)

   df=dfL[[nameDat]]
   col=data.frame(str_split_fixed(colnames(df), "-", 2))
   colAvailall=subset(col, X1==varInfo$ID)

   # Chek how many time points available
   split_values <- sapply(colAvailall$X2, function(x) {
      digits <- unlist(strsplit(as.character(x), '\\.'))
      digits <- as.numeric(digits)
      return(digits)
      })
    colAvailall$t= split_values[1, ]
    colAvailall$instance <- split_values[2, ]
    colAvailall$time_instance=paste0(colAvailall$X1, "-", colAvailall$X2)
    print(colAvailall)
  
    minInstance=min(  colAvailall$instance)
    maxInstance=max(  colAvailall$instance)
    minTime=min(  colAvailall$t)
    maxTime=max(  colAvailall$t)
    timePoint=length(levels(as.factor(  colAvailall$t)))-1
  
    print("Get corresponding date")
    dateInfo=subset(varFile, label==varInfo$date1)
    dateDat=selectDat(dfL, ID=dateInfo$ID)
    dfDate=dfL[[dateDat]]
    colDat=data.frame(str_split_fixed(colnames(dfDate), "-", 2))
    dateCols=subset(colDat, X1==dateInfo$ID)

  if(is.na(varInfo$instanceStart)==F){
          print("get instances")
          instances=seq(varInfo$instanceStart, varInfo$instanceEnd, 1) 
          for ( j in 0:(timePoint) ) {
            print(j)
            selDFmean=subset(df, select=c(paste0(varInfo$ID, "-",j,".",instances ))) # get mean
            df[[paste0(varInfo$ID, "-",j,".0")]]=rowMeans(selDFmean, na.rm=T)
            }
          }
                  

       print(paste0("Number of follow ups for ", var, ": ", timePoint))
       # Return phenotype data
       dfSel=subset(df, select=c("eid",paste0(varInfo$ID, "-",seq(minTime,maxTime, 1),".", minInstance)))
       fuPoints=seq(0, timePoint, 1)
       colnames(dfSel)=c("eid",paste0(varInfo$label, "_", fuPoints))

       # add info on first time point and assessment method
       dfAssess=firstAssess(df=dfSel, time=seq(minTime,maxTime, 1), var, fuPoints, date=varInfo$date1)
       dfSelA=cbind(dfSel, dfAssess)

       # Return dates data
       dfDat=subset(dfDate, select=c("eid", paste0(dateInfo$ID, "-",seq(minTime,maxTime, 1),".0")))
       colnames(dfDat)=c("eid",paste0("date_",varInfo$label, "_",fuPoints))
       dfCombined=merge(dfSelA, dfDat, by="eid")


       if(return=="df"){
        print(paste0("N rows extraced: ", NROW(dfCombined)))
        return(dfCombined)
        }
   }




firstAssess=function(df, time, var, fuPoints, date){
         dfOut=as.data.frame(matrix(nrow=NROW(df), ncol=length(time)))
         colnames(dfOut)=paste0("first_", var,  "_", fuPoints)

         addF=seq(1, length(fuPoints), 1)

            for ( i in 1:length(addF) ) {
              dfOut[,i]= time[i]
            }
            head(dfOut)
            
        dfOutA=as.data.frame(matrix(nrow=NROW(df), ncol=length(time)))
        colnames(dfOutA)=paste0(var,  "_", fuPoints)

        dfOutA$remove=0
        dfOutD <- as.data.frame(dfOutA) %>% mutate_all(list(assessment = ~ date))
        dfOutE=subset(dfOutD, select= paste0(var,  "_", fuPoints, "_assessment"))
        dfOutR=cbind(dfOut, dfOutE)
        return(dfOutR)
     }


readVarNames=function(do="none"){
    print("Get variable ID")
    variableList <- as.data.frame(readxl::read_excel(paste0(HOME, "/data/variableAge.xlsx")))

    if(do=="none"){
      return(variableList)
    } else {
          varInc=subset(variableList, variableList[do]=="yes")
         return(varInc)
    }

}


selectDat=function(listIn, ID, dfOutL=list()) {

  for ( i in 1:length(listIn) ) {
    #print(head(listIn[[i]]))
    col=data.frame(str_split_fixed(colnames( listIn[[i]] ), "-", 2))
    colAvailall=subset(col, X1==ID)
    dfOutL[[i]]=data.frame(it=i, name=names(listIn)[i], ncomp=NA)

    if(NROW(colAvailall)>0){
       count=subset(listIn[[i]], select=paste0(colAvailall$X1, "-", colAvailall$X2))
       na_count <-sapply(count, function(y) sum(length(which(!is.na(y)))))
       dfOutL[[i]]=data.frame(it=i, name=names(listIn)[i], ncomp=sum(na_count))
    }
  }
  
   dfOut=do.call(rbind, dfOutL)
   selDat=subset(dfOut, ncomp ==  max(dfOut$ncomp, na.rm=T))[1,]
   dfR=unique(selDat$name[1])

   print(paste0("Selected dataset: ", selDat$name, ".csv (n complete across all FUs: ", selDat$ncomp))
   return(dfR)
}




# =================== recode Pheno ===============================


checkCoding=function(ID){
    print("Extract variable info and save on cluster")
    dataDic=fread(paste0(HOME, "/data/Data_Dictionary_Showcase.csv"))
    head(dataDic)
    dataDicSel=subset(dataDic, FieldID %in% ID, select=c(FieldID, Field, Coding, Participants, ValueType, Instances, Array, Notes, Link))

    print("Read in file containing coding info of the variables")
    codingDic=read.csv(paste0(HOME, "/data/Codings.csv"))
    
    print("Check coding of the variables")
    
    for ( i in 1:NROW(dataDicSel) ) {
      coding=dataDicSel$Coding[i]
      print(paste0("Check ", dataDicSel$Field[i]))
      print(paste0("Variable type ", dataDicSel$ValueType[i]))
      varType=c("Categorical multiple", "Categorical single")

      if( any(varType %in% dataDicSel$ValueType[i]) ){
          print("Get label for categorical variable")
          codingSel=subset(codingDic, Coding==coding)
          meaningPasted=paste0(codingSel$Value, ": ", codingSel$Meaning, " | ")
          dataDicSel$Label[i]=paste0(meaningPasted, collapse="")
      } else {
            dataDicSel$Label[i]=NA
      }
   print("Save output")
 write.csv(dataDicSel,
              file= paste0(HOME, "/data/dataDicExtracted.csv"),
              row.names = FALSE,
              col.names=T,
              quote=F)
} 
}



recodeDate=function(time, var, df){
    varOut=as.Date(df[[paste0(var, "_", time)]], origin = "1900-01-01")
    print(head(varOut, n=3))
    return(varOut)
}



recodeVar=function(df, var, varInfo){
  print("=============================================")
  print(paste0("=========== ", var, " ======================"))
  print("=============================================")

  varRaw=subset(varInfo, label==var)$extract
  time=length(which(!is.na(match(colnames(df), paste0(var, "_", 0:10)))))
  vaDate=paste0("date_", var, "_",0:(time-1))

  if(varRaw!="yes"){
     funcApply <- eval(parse(text = paste(subset(varInfo, label==var)$recode)))

     dfComp=as.data.frame(extractDate(df=df, var, time=time, return="date", func=funcApply))
     timeR=length(which(!is.na(match(colnames(dfComp), paste0("date_",var, "_", 0:10)))))

     dfFirst=subset(dfComp, select=paste0("first_", var, "_",0:(timeR-1))) # 
     dfDate=subset(dfComp, select=paste0("date_", var, "_",0:(timeR-1))) #
     dfAssess=subset(dfComp, select=paste0(var, "_",0:(timeR-1), "_assessment")) # f
     time=timeR
  } else {
    dfDate=subset(df, select=vaDate)
    dfFirst=subset(df, select=paste0("first_", var, "_",0:(time-1))) # first time assessed
    dfAssess=subset(df, select=paste0(var, "_",0:(time-1), "_assessment")) # info on assessment type
  }

  print(paste0("Number of time points available: ", time))
  funcApply <- eval(parse(text = paste(subset(varInfo, label==var)$recode)))
  dfPheno=as.data.frame(recodeTime(df=df, var=var, time=time, func=funcApply)) # phenot data
  
  dfOut=cbind(dfPheno,  dfDate, dfAssess, dfFirst) # add date
  dfOut$eid=df$eid

  print(head(dfOut, 3))
  return(dfOut)
}


keepRaw=function(time, var, df){
    print("No recoding required")
    varOut=df[[paste0(var, "_", time)]]
    return(varOut)
}


removeZero=function(time, var, df){
  variableOut=df[[paste0(var, "_", time)]]
  variableOutR=ifelse(variableOut==0, NA, variableOut)
  return(variableOutR)
}

recodeKeep=function(time, var, df){
  variableOut=df[[paste0(var, "_", time)]]
  return(variableOut)
}


recodeKeep=function(time, var, df){
  varOut=df[[paste0(var, "_", time)]]
  return(varOut)
}


recodeHealthProblems=function(time, var, df){
   # -1: Do not know | -3: Prefer not to answer | -7: None of the above | 1: Diabetes related eye disease | 2: Glaucoma | 3: Injury or trauma resulting in loss of vision | 4: Cataract | 5: Macular degeneration | 6: Other serious eye condition |   variableOut=as.numeric(as.character(
   # -3: Prefer not to answer | -7: None of the above | 1: Mouth ulcers | 2: Painful gums | 3: Bleeding gums | 4: Loose teeth | 5: Toothache | 6: Dentures | 
    variableOut=as.numeric(as.character(
        revalue( as.factor(df[[paste0(var, "_", time)]]), 
                      c("-1"=NA,  # -1 : Do not know 
                        "-3"= NA, # -3 Prefer not to answer
                        "-7" = 0, #  None of the above
                        "1" = 1, # 
                        "1" = 2, # 
                        "1" = 3, #
                        "1" = 4, #
                        "1" = 5) # 
                      ) ))
  print(table(variableOut))
  return(variableOut)
}



recodeHBMD=function(time, var, df, return="pheno"){
  dfSel=subset(df, select=c(paste0("heel_bone_density_", time), paste0("heel_bone_density_left_", time)))
  variableOut <- rowMeans(dfSel)

  if(return=="date"){
      variableOut=recFirstAss(df=df, varP="srt_hearing_", time=time)
  }

  return(variableOut)
}



recodeSRT=function(time, var, df, return="pheno"){
  dfSel=subset(df, select=c(paste0("srt_hearing_", time), paste0("srt_hearing_right_", time)))
  variableOut <- rowMeans(dfSel)

  if(return=="date"){
      variableOut=recFirstAss(df=df, varP="srt_hearing_", time=time)
  }

  return(variableOut)
}


 
recodeGrip=function(time, var, df, return="pheno"){
  dfSel=subset(df, select=c(paste0("grip_strength_", time), paste0("grip_strength_right_", time)))
  variableOutR <- rowMeans(dfSel)
  variableOut=ifelse(variableOutR==0, NA, variableOutR)

  if(return=="date"){
      variableOut=recFirstAss(df=df, varP="grip_strength_", time=time)
  }

  return(variableOut)
}



recFirstAss=function(df, varP, time){
      selVar=c(paste0("date_", varP, "_", time), paste0("first_", varP, "_", time), paste0( varP, "_", time, "_assessment"))
      variableOut=subset(df, select=selVar)
      
      return(variableOut)
}


recodeHearingSRT=function(time, var, df=df, return="pheno"){
  dfSel=subset(df, select=c(paste0("SRT_right_", time), paste0("SRT_left_", time)))
  variableOut <- rowMeans(dfSel)

   if(return=="date"){
    variableOut=recFirstAss(df=df, varP="SRT_right", time=time)
      #colName=c(paste0("date_", var, "_", time), paste0("first_", var, "_", time), paste0( var, "_", time, "_assessment"))
      #colnames(variableOut)=colName
  }
  return(variableOut)
}




hearingProb=function(time, var, df){
# [-1: Do not know | -3: Prefer not to answer | 0: No | 1: Yes | 99: I am completely deaf | ]
  variableOut=as.numeric(as.character(
              revalue( as.factor(df[[paste0(var, "_", time)]]), 
                      c("-1"=NA,  # -1 : Do not know 
                        "-3"= NA, # -3 Prefer not to answer
                        "0" = 0,
                        "1" = 1,
                        "99" = 2)
                      ) ))
  print(table(variableOut))
  return(variableOut)
}




recodeHealth=function(time, var, df){
    print(paste0("Recode ", var))
      #  #-1: Do not know | -3: Prefer not to answer | 1: Excellent | 2: Good | 3: Fair | 4: Poor | (\"In general how would you rate your overall health?\)
  variableOut=as.numeric(as.character(
              revalue( as.factor(df[[paste0(var,"_", time)]]), 
               c("-1"=NA, 
                                         "-3"=NA,
                                         "1" = "4",
                                         "2" = "3",
                                         "3" = "2",
                                         "4" = "1") )))
  print(table(variableOut))
  return(variableOut)
}


hearingProb=function(time, var, df){
# [-1: Do not know | -3: Prefer not to answer | 0: No | 1: Yes | 99: I am completely deaf | ]
  variableOut=as.numeric(as.character(
              revalue( as.factor(df[[paste0(var, "_", time)]]), 
                      c("-1"=NA,  # -1 : Do not know 
                        "-3"= NA, # -3 Prefer not to answer
                        "0" = 0,
                        "1" = 1,
                        "99" = 1)
                      ) ))
  print(table(variableOut))
  return(variableOut)
}





recodeTimeOutdoors=function(time, var, df=df, return="pheno"){
  winterT=ifelse(df[[paste0("time_outdoors_winter_", time)]]<0, NA, df[[paste0("time_outdoors_winter_", time)]])
  summerT=ifelse(df[[paste0("time_outdoors_summer_", time)]]<0, NA, df[[paste0("time_outdoors_summer_", time)]])
  dfT <- data.frame(winterT, summerT)
  variableOut=rowMeans(dfT)

  if(return=="date"){
     variableOut=recFirstAss(df=df, varP="time_outdoors_summer", time=time)
  }

  return(variableOut)
}


recodeEthnicity=function(time, var, df=df){
  variableOut=as.numeric(as.character(
              revalue( as.factor(df[[paste0(var,"_", time)]]), 
                      c("-1"=NA,  # -1: Do not know 
                        "-3"= NA, # -3 Prefer not to answer
                        "1" = "1",# 1: White 
                        "1001" = "1", # 1001: British 
                        "1002" = "1", # 1002: Irish
                        "1003" = "1", # 1003: Any other white background 
                        "2" = "2",    # 2: Mixed
                        "2001" = "2", # 2001: White and Black Caribbean
                        "2002" = "2", # 2002: White and Black African 
                        "2003" = "2", # 2003: White and Asian
                        "2004" = "2", # 2004: Any other mixed background
                        "3" = "2",    # 3: Asian or Asian British 
                        "3001" = "2", # 3001: Indian 
                        "3002" = "2", # 3002: Pakistani
                        "3003" = "2", # 3003: Bangladeshi
                        "3004" = "2", # 3004: Any other Asian background
                        "4" = "2", # 4: Black or Black British
                        "4001" = "2", # 4001: Caribbean 
                        "4002" = "2", # 4002: African
                        "4003" = "2", # 4003: Any other Black background
                        "5" = "2", # 5: Chinese 
                        "6" = "2") # 6: Other ethnic group
                      ) ))
  print(table(variableOut))
  return(variableOut)
}


SRTrec=function(time, var, df){
   variableOut=as.numeric(df[[paste0(var,"_", time)]])
     return(variableOut)
}

recodeFoodIntake=function(time, var,df=df){
  variableOut=as.numeric(as.character(
              revalue( as.factor(df[[paste0(var, "_",time)]]), 
                      c("-1"=NA,  # -1 : Do not know 
                        "-3"= NA,
                        "-10"="0") # -3 Prefer not to answer
                      ) ))
  print(table(variableOut))
  return(variableOut)
}

recodeUrbanisation=function(time, var, df=df){
  variableOut=as.numeric(as.character(
              revalue( as.factor(df[[paste0(var,"_", time)]]), 
                                                    c("1"= 7, # 1: England/Wales - Urban - sparse 
                                                     "2" = 5,  # 2: England/Wales - Town and Fringe - sparse 
                                                     "3" = 3, # 3: England/Wales - Village - sparse
                                                     "4" = 1, # 4: England/Wales - Hamlet and Isolated dwelling - sparse 
                                                     "5" = 8, # 5: England/Wales - Urban - less sparse 
                                                     "6" = 6, # 6: England/Wales - Town and Fringe - less sparse 
                                                     "7" = 4, # 7: England/Wales - Village - less sparse 
                                                     "8" = 2, # 8: England/Wales - Hamlet and Isolated Dwelling - less sparse 
                                                     "9" = NA, # 9: Postcode not linkable 
                                                     "11" = 8, #  11: Scotland - Large Urban Area 
                                                     "12" = 7, # 12: Scotland - Other Urban Area
                                                     "13" = 6, # 13: Scotland - Accessible Small Town 
                                                     "14" = 5, # 14: Scotland - Remote Small Town
                                                     "15" = 4, # 15: Scotland - Very Remote Small Town
                                                     "16" = 3, # 16: Scotland - Accessible Rural
                                                     "17" = 2, # 17: Scotland - Remote Rural
                                                     "18" = 1 #  18: Scotland - Very Remote Rural 
                                                    ))))
  print(table(variableOut))
  return(variableOut)
}


recodeProspectiveMem=function(time, var, df){
    print(paste0("Recode ", var))
    # prospective_memory
# [0: Instruction not recalled, either skipped or incorrect | 1: Correct recall on first attempt | 2: Correct recall on second attempt | ]
  variableOut=as.numeric(as.character(
              revalue( as.factor(df[[paste0(var,"_", time)]]), 
               c("-1"=NA, 
                                         "-3"=NA,
                                         "1" = "2",
                                         "2" = "1",
                                         "0" = "0") )))
  print(table(variableOut))
  return(variableOut)
}

recodeFriendship=function(time, var,df=df){
  # [ -1: Do not know | -3: Prefer not to answer | 1: Extremely happy | 2: Very happy | 3: Moderately happy | 4: Moderately unhappy | 5: Very unhappy | 6: Extremely unhappy | ]

  variableOut=as.numeric(as.character(
              revalue( as.factor(df[[paste0(var,"_", time)]]), 
                      c("-3"=NA,
                      "-1"=NA, 
                      "6"="0", # 6: Extremely unhappy
                      "5" = "1", # 5: Very unhappy
                      "4" = "2", # 4: Moderately unhappy
                       "3" = "3", # 3: Moderately happy
                      "2" = "4", # 2: Very happy
                      "1" = "5") # 1: Extremely happy
                      ) ))
  print(table(variableOut))
  return(variableOut)
}

recodeAlcFreq=function(time, var,df=df){
  variableOut=as.numeric(as.character(
              revalue( as.factor(df[[paste0(var,"_", time)]]), 
                      c("-3"=NA,
                      "6"="0", # never
                      "5" = "1", # Special occasions only 
                      "4" = "2", #One to three times a month 
                       "3" = "3", # Once or twice a week 
                      "2" = "4", #Three or four times a week 
                      "1" = "5") # daily
                      ) ))
  print(table(variableOut))
  return(variableOut)
}

hearingProb=function(time, var, df){
# [-1: Do not know | -3: Prefer not to answer | 0: No | 1: Yes | 99: I am completely deaf | ]
  variableOut=as.numeric(as.character(
              revalue( as.factor(df[[paste0(var, "_", time)]]), 
                      c("-1"=NA,  # -1 : Do not know 
                        "-3"= NA, # -3 Prefer not to answer
                        "0" = 0,
                        "1" = 1,
                        "99" = 1)
                      ) ))
  print(table(variableOut))
  return(variableOut)
}


#recodePhysicalAct=function(){
 # physical_activity_2 # [0: None | 10: Under 10 minutes | 1030: 10-30 minutes | 12: 1-2 hours | 24: 2-4 hours | 3060: 30-60 minutes | 46: 4-6 hours | 600: 6+ hours | ]
#}


recodeSmokeFreq=function(time, var, df=df, return="pheno"){
  variableOut=rep(NA, NROW(df))
  variableOut=ifelse(df[[paste0("smoking_status_", time)]]=="0", "0", variableOut)
  variableOut=ifelse(df[[paste0("smoking_status_", time)]]=="1", removeNegative(time=time, var="smoking_frequency_former", df=df), variableOut)
  variableOut=ifelse(df[[paste0("smoking_status_", time)]]=="2", removeNegative(time=time, var="smoking_frequency", df=df), variableOut)
  variableOut=as.numeric(variableOut)
  variableOut=ifelse(variableOut<0, NA, variableOut)

  if(return=="date"){
      variableOut=recFirstAss(df=df, varP="smoking_status", time=time)
  }

  return(variableOut)
}

recodeSmokeAge=function(time, var,df=df, return="pheno"){
  variableOut=ifelse( is.na(df[[paste0("age_smoke_current_", time)]])==F, df[[paste0("age_smoke_current_", time)]],  df[[paste0("age_smoke_former_", time)]])
  variableOut=as.numeric(variableOut)

    if(return=="date"){
      variableOut=recFirstAss(df=df, varP="age_smoke_current", time=time)
  }

  return(variableOut)
}

recodeChrono=function(time, var, df=df){
  variableOut=as.numeric(as.character(
              revalue( as.factor(df[[paste0(var, "_", time)]]), 
                      c("1"="2",  # 1 Definitely a 'morning' person
                        "2"="1",  # 2 More a 'morning' than 'evening' person
                        "-1"="0", # -1 Do not know
                        "3"="-1", # 3 More an 'evening' than a 'morning' person
                        "4"="-2", # 4 Definitely an 'evening' person
                        "-3"= NA) # -3 Prefer not to answer
                      ) ))
  print(table(variableOut))
  return(variableOut)
}

recodeAsthma=function(time, var, df=df){
  # -3: Prefer not to answer | -7: None of the above | 5: Blood clot in the leg (DVT) | 6: Emphysema/chronic bronchitis | 7: Blood clot in the lung | 8: Asthma | 9: Hayfever, allergic rhinitis or eczema | 
  variableOut=as.numeric(as.character(
              revalue( as.factor(df[[paste0(var, "_", time)]]), 
                      c("-3"= NA, # -3 Prefer not to answer
                        "-7"="0", # -7: None of the above
                        "5"="0",  # 5: Blood clot in the leg (DVT)
                        "6"="0",  # 6: Emphysema/chronic bronchitis
                        "7"="0",  # 7: Blood clot in the lung
                        "8"="1",  # 8: Asthma
                        "9"="0") # 9: Hayfever, allergic rhinitis or eczema
                      ) ))
  print(table(variableOut))
  return(variableOut)
}

recodeSubstanceStatus=function(time, var, df){

  variableOut=as.numeric(as.character(
               revalue( as.factor(df[[paste0(var, "_", time)]]), 
                      c("-1"=NA, # -1: Do not know 
                        "-3"=NA, # -3: Prefer not to answer 
                        "0"=0, # 0: Never
                        "1"=0, # 1: Previous 
                        "2"=1 # 2: Current
                        )))) 

 print(table(variableOut))
  return(variableOut)
}



recodeBinary=function(time, var, df=df){
  variableOut=as.numeric(as.character(
              revalue( as.factor(df[[paste0(var,"_", time)]]), 
                      c("-1"=NA,  # -1 : Do not know 
                        "-3"= NA, # -3 Prefer not to answer
                        "-818"=NA,
                        "-121"=NA,
                        "1" = "1",# yes
                        "0" = "0") # no
                      ) ))
  print(table(variableOut))
  return(variableOut)
}

removeNegative=function(time, var, df){
  varE=df[[paste0(var, "_", time)]]
  varOut=ifelse(varE<0, NA, varE)
  return(varOut)
}


recodeMedsNo=function(time, var, df){
  varE=df[[paste0(var, "_", time)]]
  varOut=ifelse(varE<0, NA, varE)
  varOutF=ifelse(varOut>=8, 8, varOut)
  return(varOutF)
}


recodeOverallHealth=function(time, var, df){
    varE=df[[paste0(var, "_", time)]]
    # -1: Do not know | -3: Prefer not to answer | 1: Excellent | 2: Good | 3: Fair | 4: Poor | 
    variable_rec=revalue(as.character(varE), c("-1"=NA, 
                                         "-3"=NA,
                                         "1" = "3",
                                         "2" = "3",
                                         "3" = "2",
                                         "4" = "1") )
  return(as.numeric(variable_rec))
}

recodeTime=function(df, var, time, func){
  dfOut= as.data.frame( matrix(NA, nrow=NROW(df), ncol=time))
  print(paste0("Recode ", var, " (", time, " time points)"))
  for(i in 0:(time-1)){
      dfOut[,(i+1)]=func(time=i, var=var, df=df)
  }
  colnames(dfOut)=paste0(var, "_", seq(0,time-1, 1))
return(dfOut)
}

extractDate=function(df, var, time, func, return){
  listOut=list()
for (i in 0:10) {
  tryCatch({
    dfRet <- func(time = i, var = var, df = df, return = return)
    # Append the result to the list only if there is no error
    colnames(dfRet)=c(paste0("date_", var, "_", i), paste0("first_", var, "_", i), paste0(var, "_", i, "_assessment"))
    listOut[[i +1]]=dfRet
    # Continue with the rest of your loop logic
  }, error = function(e) {
    # Handle the error (optional)
    #cat("Error occurred in iteration", i, ": ", conditionMessage(e), "\n")
    # No break needed in this case; continue to the next iteration
  })
}
dfOut=do.call(cbind, listOut)
return(dfOut)
}


# =========================================================================
# ========================= AGE EFFECTS ===================================
# =========================================================================



  getAgeWithin=function(df, var){ # Longitudinal age effects
      dfS=na.omit(subset(df, select=c(pheno, age, ID)))
      n=length(levels(droplevels(as.factor(dfS$ID))))
      mW <-lmer( pheno ~ age + (1|ID), data=dfS)
      mWDF=as.data.frame(coef(summary(mW)))
      dfOut=data.frame(var=var, type=c("within_age"), beta=mWDF$Estimate[2], se=mWDF[["Std. Error"]][2],  n= n)
      return(dfOut)
  }

 getAgeWithinQ=function(df, var){ # Quadratic age effects (longitudinal)
      n=length(levels(droplevels(as.factor( df$ID))))
      mW <-lmer( pheno ~ age +  I(age^2) + (1|ID), data=df)
      mWDF=as.data.frame(coef(summary(mW)))
      dfOut=data.frame(var=var, type=c("non_linear_age"), beta=mWDF$Estimate[3], se=mWDF[["Std. Error"]][3], n= n)
      return(dfOut)
  }


lmeModel=function(df, var, model="cohort"){ # different mixed effects models

  if(model=="cohort"){
      dfS=na.omit(subset(df, select=c(pheno, age_raw, cohort, ID)))
      lC<-lmer( pheno ~ age_raw + cohort + (1|ID), data=dfS)
      m=model
  }
  if(model=="snp_cohort"){
       dfS=na.omit(subset(df, select=c("pheno", "age_raw", "cohort", "ID", "snp", "batch", paste0("PC", 1:20)) ))
       f=as.formula(paste0( "pheno ~ age_raw + cohort + snp + snp:age_raw + snp:cohort + as.factor(batch) + ",  paste0("PC", 1:20, collapse=" + "), " + (1|ID)"))
       lC<-lmer(f, data=dfS)
       m="cohort:snp"
  }
  
  if(model=="snp_age_lme"){
       dfS=na.omit(subset(df, select=c("pheno", "age_raw", "ID", "snp", "batch", paste0("PC", 1:20)) ))
       f=as.formula(paste0( "pheno ~ age_raw + snp + snp:age_raw + as.factor(batch) + ",  paste0("PC", 1:20, collapse=" + "), " + (1|ID)"))
       lC<-lmer(f, data=dfS)
       m="age_raw:snp"
  }
  
    if(model=="snp_age2_lme"){
       dfS=na.omit(subset(df, select=c("pheno", "age_raw", "ID", "snp", "batch", paste0("PC", 1:20)) ))
       f=as.formula(paste0( "pheno ~ age_raw + I(age_raw^2) + snp + snp:age_raw + snp:I(age_raw^2) + as.factor(batch) + ",  paste0("PC", 1:20, collapse=" + "), " + (1|ID)"))
       lC<-lmer(f, data=dfS)
       m="I(age_raw^2):snp"
  }

    n=length(levels(droplevels(as.factor(dfS$ID))))
    lcDF=as.data.frame(coef(summary(lC)))
    lcDF$label=rownames(lcDF)
    lcS=subset(lcDF, label==m)
    dfOut=data.frame(var=var, type=model, beta=lcS[["Estimate"]], se=lcS[["Std. Error"]], n=n)
    return(dfOut)
}




getAge=function(df,  var){
  lmAge = lm(pheno ~ age, data = df) # linear regression (cross-sectional)
  m=as.data.frame(coef(summary(lmAge)))
  dfOut=data.frame(var=var, type=c("between_age"), beta=m$Estimate[2], se=m[["Std. Error"]][2], n=NROW(na.omit(subset(df, select=c(pheno, age)))))
  return(dfOut)
}



ageModelW=function(df, var, model="age"){ # Inverse-probability weighted age effects

      if(model=="snp"){
          df1=data.frame(pheno=df$pheno, age=df$age, weight=df$w, snp=df$snp, sex=df$SEX, batch=df$batch)
          df2=subset(df, select=c(paste0("PC", 1:20)))
          datIN=na.omit(cbind(df1, df2))
          f=as.formula(paste0("pheno ~ snp + snp:age + sex + age + as.factor(batch) + ",  paste0("PC", 1:20, collapse=" + ")))
          label="age_snp_weighted"
          extract="snp:age"
      }
      if(model=="age"){
         datIN=na.omit(data.frame(pheno=df$pheno, age=df$age, weight=df$w))
         f <- as.formula(paste0("as.numeric(pheno) ~ as.numeric(age)")) # I(as.numeric(age)^2) + 
         label="age_weighted"
         extract="as.numeric(age)"
      }
      datIN$wn=datIN$w/mean(datIN$w)
      nEff=round((sum(datIN$wn)^2)/sum(datIN$wn^2),0)
      design.ps <- svydesign(ids=~1,   weights = ~wn, data=datIN)
      
      regglmSTD=svyglm(f, design=design.ps)
      coef=as.data.frame(   coef(summary(regglmSTD)))
      coef$label=rownames(coef)
      coefS=subset(coef, label==extract)
      dfOut=data.frame(var=var, type=label, beta=coefS[["Estimate"]][1], se=coefS$`Std. Error`[1], n=nEff)
      return(dfOut)
    }



ageModel=function(var, return=NULL, data, datAge=NULL, snp=NULL, gene=NULL){ # obtain age effects across all models
  dfSum=NA
  print(paste0("+++++++++++ ", var, " +++++++++++"))
  data$age_raw=data[["age"]]
  data$m=NULL
  head(data)
  print("Prepare longitudinal data")
  data <- data %>% dplyr::group_by(ID) %>% dplyr::mutate(age = age_raw - mean(age_raw)) %>% dplyr::ungroup()
  #data <- data %>% dplyr::group_by(ID) %>% dplyr::mutate(time = age_raw - min(age_raw)) %>% dplyr::ungroup()

  baselineM=mean(datAge[[paste0(var, "_0")]], na.rm=TRUE) #subset(means, pheno==var)$mean
  baselineSD=sd(datAge[[paste0(var, "_0")]], na.rm=TRUE) #subset(means, pheno==var)$sd
  data$phenoWithin=data[[var]] 
  data$pheno = (  data$phenoWithin - baselineM) / baselineSD
  data$phenoWithin=NULL;data[[var]]=NULL
  data$cohort=as.numeric(data$year)-data$age_raw
  data2FU <- as.data.frame(data %>% dplyr::group_by(ID) %>% dplyr::filter(dplyr::n() == 3) %>% dplyr::ungroup())
  dataW=merge(data, data.frame(ID=datAge$eid, w=datAge$w), by="ID", all.x=TRUE) # add sampling weights to longitudinal data

  print("Prepare cross-sectional data")
  datAgeA=datAge
  datAge = data.frame(eid=datAgeA$eid, pheno=(datAgeA[[paste0(var, "_0")]] - baselineM) / baselineSD, pheno_t1=datAgeA[[paste0(var, "_1")]], age=datAgeA[[paste0("age_",var, "_0")]], w=datAgeA$w)


  if(NROW(colnames(snp))==0){
    print("Get within person ages effects")
    ageWitin=getAgeWithin(df=data, var=var)
  
    print("Get within person quadratic age effects")
    ageWitin2=getAgeWithinQ(df=data2FU, var=var)

    print("Get between person ages effects (in all participants)")
    ageBetween= getAge(df=datAge, var=var)

    print("Get between person ages effects (in FU participants)")
    datFU=subset(datAge, is.na(pheno)==F & is.na(pheno_t1)==F )
    ageBetweenFU= getAge(df=datFU, var=var)
    ageBetweenFU$type="between_age_FU"

    print('Get weighted age effects (cross-sectional)')
    ageBetweenW= ageModelW(df=datAge, var=var)

    print('Get cohort effects')
    ageCohort = lmeModel(df=data, var=var, model="cohort")

    print('Combine all results')
    dfSum=rbind(ageWitin, ageWitin2, ageBetween, ageBetweenFU, ageBetweenW, ageCohort ) #ageCohort2
    dfSum$P <- 2 * pnorm(-abs(dfSum$beta/dfSum$se))
    return(dfSum)
    }

    if(NROW(snp)>=1){

    snpVar=levels(as.factor(snp$SNP))
    listSNP=list()
    for(i in 1:NROW(snpVar)){
      s=snpVar[i]
      print(paste0('Get results for ',s, ' on ', var, " (iteration ",i, " out of ", NROW(snpVar), ")" ))
      regenie=readRDS(paste0(HOME,"/output/rds/clump.rds"))
      regenieS=subset(regenie, pheno==var & SNP==s)
      colSNP=regenieS[1,]

      aWithinG=subset(regenieS, type=="interaction_within")
      ageGeneWithin=data.frame(var=var, type="regenie_age_snp_within", beta=aWithinG$BETA, se=aWithinG$SE, n=aWithinG$N)
      aBetweenG=subset(regenieS, type=="interaction_between")
      ageGeneBetween=data.frame(var=var, type="regenie_age_snp_between", beta=aBetweenG$BETA, se=aBetweenG$SE, n=aBetweenG$N)
      
      gDF=subset(gene, select=c("eid", s,"SEX", "batch", paste0("PC", 1:20)))
      names(gDF)[names(gDF) == s] <- 'snp'
      dataG=merge(data, gDF, by.x="ID", by.y="eid", all.x=TRUE ) # longitudinal + genetics
      datAgeG=merge(datAge, gDF, by="eid") # cross-sectional + genetics

      print("Test for interactions with SNPs: cross-sectional")
      ageGeneBetween_lm=effectGene(dfC=datAgeG, s=s, type="interaction_between", pheno=var)

      print("Test for interactions with SNPs: cross-sectional in follow-up sample")
      ageGeneBetweenFU_lm=effectGene(dfC=datAgeG, dfL=dataG, s=s, type="interaction_between_FU", pheno=var)

      print("Test for interactions with SNPs: longitudinal")
      ageGeneWithin_lm=effectGene(dfL=dataG, s=s, type="interaction_within", pheno=var) # from lm
      snpAgeLme = lmeModel(df=dataG, var=var, model="snp_age_lme") # from lme

      print("Test for interactions with SNPs: Cohort")
      snpCohort = lmeModel(df=dataG, var=var, model="snp_cohort")


      print("Test for interactions with SNPs: non-linear age effects")
      dataG2FU=subset(dataG, ID %in% levels(as.factor(data2FU$ID)))
      snpAge2Lme = lmeModel(df=dataG2FU, var=var, model="snp_age2_lme") # from lme (age is)

      print('Get weighted age effects')
      ageGeneW= ageModelW(df=datAgeG, var=var, model="snp") 
    
      print('Combine all results')
      dfSum=rbind(ageGeneWithin, ageGeneBetween, ageGeneBetween_lm, ageGeneBetweenFU_lm, ageGeneWithin_lm, snpCohort, ageGeneW, snpAgeLme, snpAge2Lme ) #snpCohort2
      dfSum$SNP=s;dfSum$CHR=colSNP$CHR;dfSum$POS=colSNP$POS;dfSum$A1=colSNP$A1;dfSum$A2=colSNP$A2
      dfSum$P <- 2 * pnorm(-abs(dfSum$beta/dfSum$se))
      listSNP[[i]]=dfSum
  }
  dfSum=do.call(rbind, listSNP)
  return(dfSum)
  }


 }
  


deriveFU=function(df, var, varDF, sub="all"){
  print(paste0("============== ", var, " ==================="))
  varInfo=subset(varDF, label==var)

  tSum=length(which(!is.na(match(colnames(df), paste0(var, "_", 0:10)))))-1
  t=c(paste0("date_",var, "_", 0:tSum))
  p=c(paste0(var, "_", 0:tSum))


  varLongAll=data.frame(eid=rep(df$eid, length(t)),
                        time=do.call(rbind, lapply(0:tSum, function(x) data.frame(time=rep(x, NROW(df) )))),
                        pheno=do.call(rbind, lapply(p, function(x) data.frame(pheno=df[[x]]))),
                        date=do.call(rbind, lapply(t, function(x) data.frame(date=df[[x]]))))
 varLong=subset(varLongAll, is.na(pheno)==F)

 varLongO=varLong %>%
    dplyr::group_by(eid) %>%
    dplyr::arrange(date, .by_group=TRUE)

 print(paste0('Missing data for dates (should be 0, as only complete data for phenotypes are included): ', NROW(subset(varLongO, is.na(date)==T))) )

 varLongID=as.data.frame(varLongO %>%
                            dplyr::group_by(eid) %>%
                            dplyr::mutate(id=dplyr::row_number()))

  tAvail=as.numeric(levels(as.factor(varLongID$id)))
  print(paste0('Maximum number of follow ups: ', max(varLongID$id)))
  dfL=lapply(tAvail, function(x) selCol(df=varLongID, idSel=x, var=var))
  dfC=Reduce(function(x,y) merge(x = x, y = y, by = "eid", all=T), dfL)
  dfC$eid=as.character(dfC$eid)
  # Filter the data between 2006-03-13 and 2010-10-01 (first assessment center period)
  dfC$dateBaseline=dfC[[paste0("date_", var, "_0")]]
  dfF <- subset(dfC, dateBaseline >= as.Date("2006-03-13") & dateBaseline <= as.Date("2010-10-01"))
  dfF$dateBaseline
  rm(dfC, dfL)

  print("Add age at each follow up per phenotype")
  dfA=merge(subset(df, select=c("eid", paste0( "date_age_", 0), paste0( "age_", 0), "w")), dfF, by="eid", all.x=T)
  ageFUPL=lapply(tAvail-1, function(x) getAgeDiff(df=dfA, var=var, time=x))
  ageFUP=do.call(cbind, ageFUPL)
  dfF=cbind(dfA, ageFUP)
  dfF[[paste0( "date_age_", 0)]]=NULL
  dfF[[paste0( "age_", 0)]]=NULL


  if(sub=="COVID"){
    print('Recode post-COVID assessments as NA')
  print("Check baseline date criteria")
  dfF$dateR0=as.Date(dfF[[paste0("date_",var, "_0")]])
  dfF[[paste0(var, "_0")]]=ifelse(  dfF$dateR0 > as.Date("2010-10-01"), NA,   dfF[[paste0(var, "_0")]])
  dfF[[paste0("age_",var, "_0")]]=ifelse(is.na( dfF[[paste0(var, "_0")]])==TRUE, NA,   dfF[[paste0("age_",var, "_0")]])
  dfF[[paste0("date_",var, "_0")]]=as.Date(ifelse(is.na( dfF[[paste0(var, "_0")]])==TRUE, NA,   dfF[[paste0("date_",var, "_0")]]))

  print("Check t1 date criteria")
  dfF$dateR1=as.Date(dfF[[paste0("date_",var, "_1")]])
  dfF[[paste0(var, "_1")]]=ifelse(  dfF$dateR1 > as.Date("2019-12-31"), NA,   dfF[[paste0(var, "_1")]])
  dfF[[paste0("age_",var, "_1")]]=ifelse(is.na( dfF[[paste0(var, "_1")]])==TRUE, NA,   dfF[[paste0("age_",var, "_1")]])
  dfF[[paste0("date_",var, "_1")]]=as.Date(ifelse(is.na( dfF[[paste0(var, "_1")]])==TRUE, NA,   dfF[[paste0("date_",var, "_1")]]))

  print("Check t2 date criteria")
  dfF$dateR2=as.Date(dfF[[paste0("date_",var, "_2")]])
  dfF[[paste0(var, "_2")]]=ifelse(  dfF$dateR2 > as.Date("2019-12-31"), NA,   dfF[[paste0(var, "_2")]])
  dfF[[paste0("age_",var, "_2")]]=ifelse(is.na( dfF[[paste0(var, "_2")]])==TRUE, NA,   dfF[[paste0("age_",var, "_2")]])
  dfF[[paste0("date_",var, "_2")]]=as.Date(ifelse(is.na( dfF[[paste0(var, "_2")]])==TRUE, NA,   dfF[[paste0("date_",var, "_2")]]))

  if(NROW(na.omit(dfF[[paste0("date_",var, "_3")]]))>0){
    print("Check t3 date criteria")
    dfF$dateR3=as.Date(dfF[[paste0("date_",var, "_3")]])
    dfF[[paste0(var, "_3")]]=ifelse(  dfF$dateR3 > as.Date("2019-12-31"), NA,   dfF[[paste0(var, "_3")]])
    dfF[[paste0("age_",var, "_3")]]=ifelse(is.na( dfF[[paste0(var, "_3")]])==TRUE, NA,   dfF[[paste0("age_",var, "_3")]])
    dfF[[paste0("date_",var, "_3")]]=as.Date(ifelse(is.na( dfF[[paste0(var, "_3")]])==TRUE, NA,   dfF[[paste0("date_",var, "_3")]]))
  } 
  if(NROW(na.omit(dfF[[paste0("date_",var, "_4")]]))>0){
    print("Check t4 date criteria")
    dfF$dateR4=as.Date(dfF[[paste0("date_",var, "_4")]])
    dfF[[paste0(var, "_4")]]=ifelse(  dfF$dateR4 > as.Date("2019-12-31"), NA,   dfF[[paste0(var, "_4")]])
    dfF[[paste0("age_",var, "_4")]]=ifelse(is.na( dfF[[paste0(var, "_4")]])==TRUE, NA,   dfF[[paste0("age_",var, "_4")]])
    dfF[[paste0("date_",var, "_4")]]=as.Date(ifelse(is.na( dfF[[paste0(var, "_4")]])==TRUE, NA,   dfF[[paste0("date_",var, "_4")]]))
  } 
    dfF$dateR0=NULL; dfF$dateR1=NULL; dfF$dateR2=NULL; dfF$dateR3=NULL; dfF$dateR4=NULL
    dfF$eid=as.character(dfF$eid)
    print("Check if dates format is correct:")
    print(  head(    dfF[[paste0("date_",var, "_0")]], n=4))
  }

  return(dfF)
}



selCol=function(df, idSel, var){
    dfSel=subset(df, id==idSel, select=c("eid",  "pheno", "date"))
    colnames(dfSel)=c("eid", paste0(var, "_", idSel-1), paste0("date_",var, "_", idSel-1))
    return(dfSel)
}


getAgeDiff=function(dateBL, ageBL, time, df, var){
    dateFU=df[[paste0("date_",var, "_", time)]]
    dateBL=df[[paste0( "date_age_", 0)]]
    ageBL=df[[paste0( "age_", 0)]]

    diffDate= round(lubridate::time_length(difftime(as.Date(dateFU, origin = "1900-01-01"), as.Date(dateBL, origin = "1900-01-01") ), "years"), 0)
    ageFU=ageBL+diffDate
    print(paste0("Years after initial UKBB assessment for time point ", time))
    print(table(diffDate))
    dfOut=data.frame(ageFU)
    colnames(dfOut)=c(paste0("age_", var, "_",time))
    return(dfOut)
}


getBL=function(var, dfFU, dfBL){

  if(paste0(var, "_0") %in% colnames(dfFU)){
    print(paste0(var, " in FU dataset"))
    df=dfFU
  } else {
    df=dfBL
  }

  print(paste0("Get baseline data for ", var))
   dfOut=subset(df, select=c("eid", paste0(var, "_0") ))
   colnames(dfOut)=c("eid", paste0(var, "_baseline"))  
   dfOut[[paste0(var, "_0")]]=dfOut[[paste0(var, "_baseline")]]
   return(dfOut)
} 



prepareAge=function(df, var, ID){
    print(var)
    timeP=length(which(!is.na(match(colnames(df), paste0(var, "_", 0:10)))))-1
    #timeP=1 # select first follow up only
    #timeP=2
    longP=getLongFE(var=var, time=timeP, ID=ID, df=df)
    longO <- longP[order(longP$ID ),]
    return(longO)
}


getLongFE=function(df, var, time, ID){
    print(paste0("Get long format for ", var, ", including ", time, " FU time points"))
      aP=c(paste0("age_",var, "_", 0:time))
      dP=c(paste0("date_",var, "_", 0:time))
      pP=c(paste0(var, "_", 0:time))

      longP=data.frame(ID=rep(df[[ID]], length(aP)),
                       #w=rep(df[["w"]], length(aP)),
                        do.call(rbind, lapply(pP, function(x) data.frame(pheno=df[[x]]))),
                        do.call(rbind, lapply(dP, function(x) data.frame(date=df[[x]]))),
                        do.call(rbind, lapply(aP, function(x) data.frame(age=df[[x]]))))
      longPc=subset(longP, is.na(pheno)==F)
      longPc$m=paste0(longPc$ID, "_", longPc$age)
      
      repM <- unique(longPc$ID[duplicated(longPc$ID) | duplicated(longPc$ID, fromLast = TRUE)])
      print(paste0(NROW(repM), " individuals with multiple measurements"))
      longPs=subset(longPc, longPc$ID %in% repM)

      names(longPs)[names(longPs) == 'pheno'] <- var
      names(longPs)[names(longPs) == ID ] <- "ID"
      longPs$year=format(as.Date(longPs$date), "%Y")
      longPs$date=NULL
      return(longPs)
  } 


###################################################################################################
# =================================== Get longitudinal data =======================================
###################################################################################################

apcI=function(df, var){ # get info on age, cohort period characteristics in UKB
    datL=combineFU(df=df, var=var)
    
    ageFU=datL[[paste0("age_",var,"_0")]] + datL$fuTime
    fuDat=data.frame(eid=datL$eid, dateFU=datL[[paste0("date_",var,"_1")]], ageFU=ageFU)
    blDat=data.frame(eid=df$eid, dateBL=df[[paste0("date_",var,"_0")]], ageBL=df[[paste0("age_",var,"_0")]])
    dfM=merge(blDat, fuDat, all.x=TRUE)

    yearBaselineDF=as.numeric(format(as.Date(dfM$dateBL), "%Y"))
    yearBSUM=data.frame(type="year_baseline", info=c("mean", "min", "max"), est=c(mean(yearBaselineDF, na.rm=TRUE), min(yearBaselineDF, na.rm=TRUE), max(yearBaselineDF, na.rm=TRUE)))
    
    yearFUDF=as.numeric(format(as.Date(dfM$dateFU), "%Y"))
    yearFUSUM=data.frame(type="year_fu", info=c("mean", "min", "max"), est=c(mean(yearFUDF, na.rm=TRUE), min(yearFUDF, na.rm=TRUE), max(yearFUDF, na.rm=TRUE)))

    fuTimeDF=yearFUDF-yearBaselineDF
    fuTimeSUM=data.frame(type="fu_time", info=c("mean", "min", "max"), est=c(mean(fuTimeDF, na.rm=TRUE), min(fuTimeDF, na.rm=TRUE), max(fuTimeDF, na.rm=TRUE)))

    ageBaselineDF=dfM$ageBL
    ageBaselineSUM=data.frame(type="age_baseline", info=c("mean", "min", "max"), est=c(mean(ageBaselineDF, na.rm=TRUE), min(ageBaselineDF, na.rm=TRUE), max(ageBaselineDF, na.rm=TRUE)))

    ageFUDF=dfM$ageFU
    ageFUSUM=data.frame(type="age_fu", info=c("mean", "min", "max"), est=c(mean(ageFUDF, na.rm=TRUE), min(ageFUDF, na.rm=TRUE), max(ageFUDF, na.rm=TRUE)))

    cohortDF=yearBaselineDF-ageBaselineDF
    cohortSUM=data.frame(type="cohort", info=c("mean", "min", "max"), est=c(mean(cohortDF, na.rm=TRUE), min(cohortDF, na.rm=TRUE), max(cohortDF, na.rm=TRUE)))

    dfSum=rbind(yearBSUM, yearFUSUM, fuTimeSUM, ageBaselineSUM, ageFUSUM, cohortSUM)
    dfSum$estR=round(dfSum$est, 0)
    dfSum$n_BL=NROW(na.omit(ageBaselineDF))
    dfSum$n_FU=NROW(na.omit(ageFUDF))
    dfSum$var=var
    return(dfSum)
}


combineFU=function(df, var){
  pheno1=df[[paste0(var, "_", 1)]]
  date1=df[[paste0("date_",var, "_", 1)]]
  fuPheno=df[[paste0(var, "_", 1)]]

  if(NROW(na.omit(df[[paste0(var, "_", 2)]]))>0){
    print("Add FU with longer FU time")
    date1=data.table::fifelse(is.na(df[[paste0(var, "_2")]])==F, df[[paste0("date_",var, "_2")]], date1)
    fuPheno=data.table::fifelse(is.na(df[[paste0(var, "_2")]])==F, df[[paste0(var, "_2")]], pheno1)

  if(NROW(na.omit(df[[paste0(var, "_", 3)]]))>0){
    date1=data.table::fifelse(is.na(df[[paste0(var, "_3")]])==F, df[[paste0("date_",var, "_3")]], date1)
    fuPheno=data.table::fifelse(is.na(df[[paste0(var, "_3")]])==F, df[[paste0(var, "_3")]], pheno1)
   }
   }
 
  df$fuTime= round(lubridate::time_length(difftime(as.Date(date1, origin = "1900-01-01"), as.Date(df[[paste0("date_",var, "_0")]], origin = "1900-01-01") ), "years"), 0)

  dfSel=subset(df, select=c("eid", paste0("age_",var, "_0") , "fuTime"))
  dfSel[[paste0(var, "_0")]]=df[[paste0(var, "_0")]]
  dfSel[[paste0(var, "_1")]]=fuPheno
  dfSel[[paste0("date_", var, "_0")]]= as.Date(df[[paste0("date_",var, "_", 0)]], origin = "1900-01-01")
  dfSel[[paste0("date_", var, "_1")]]=as.Date(date1, origin = "1900-01-01")

  print(paste0("Min baseline date: ", min(  dfSel[[paste0("date_", var, "_0")]], na.rm=TRUE), " - max baseline date: ",  max(  dfSel[[paste0("date_", var, "_0")]], na.rm=TRUE)))
  print(paste0("Min follow up date: ", min(  date1, na.rm=TRUE), " - max baseline date: ",  max(  date1, na.rm=TRUE)))

  dfSelC <- dfSel[complete.cases(dfSel), ]
  print(paste0("Number of complete cases: ", NROW(dfSelC)))
  return(dfSelC)
}



changeFunc=function(var, df, return="data", col=""){
  print(paste0("========= PROCESS ", toupper(var), " ==================="))
  varSel=combineFU(df=df, var=var)
  varQ=varSel
  print(paste0("Max follow up time: ", max(varSel$fuTime)))
  print(paste0("Min follow up time: ", min(varSel$fuTime)))

  print("Remove outliers")
  t0=varSel[[paste0(var, "_", 0)]]
  t1=varSel[[paste0(var, "_", 1)]]
  m <- abs(mean(t0, na.rm=T))
  sd <- sd(t0, na.rm=T)
  remVal <- m + 10 * sd
  t0c=ifelse(abs(t0) > remVal, NA, t0)
  t1c=ifelse(abs(t1) > remVal, NA, t1)
  nREM=NROW(na.omit(c(t0, t1)))-NROW(na.omit(c(t0c, t1c)))
  print(paste0("Number of outliers removed: ", nREM))

  varSel[[paste0(var, "_", 0)]]=t0c
  varSel[[paste0(var, "_", 1)]]=t1c

  baselineM=mean(varSel[[paste0(var, "_", 0)]], na.rm=TRUE) 
  baselineSD=sd( varSel[[paste0(var, "_", 0)]], na.rm=TRUE)
  t0Scaled=(   varSel[[paste0(var, "_", 0)]] - baselineM) / baselineSD
  t1Scaled=(   varSel[[paste0(var, "_", 1)]] - baselineM) / baselineSD
 
  print("Get slope of change")
  slope=(t1Scaled - t0Scaled) / varSel$fuTime
  slope[is.infinite(slope)] <- NA
  varSel$change=slope
  names(varSel)[names(varSel) == 'change'] <- paste0(var, "_change")
  varSelM=merge(varSel, subset(df, select=c(eid)), by="eid", all.x=T)
  varSelC=subset(varSelM, fuTime>2)
  print(paste0("N observations removed because of too short FU (<3 year): ", NROW(varSelM)-NROW(varSelC)))

 if(return=="data"){
    varSelR=subset(varSelC, select=c("eid", paste0(var, "_0" ), paste0("age_", var, "_0" ), paste0( var, "_change" )))
    return(varSelR)
 } 
 if(return=="timepoints"){
    age0=varSel[[paste0("age_", var, "_0")]]
    age1=age0+varSel$fuTime
    dfTimeP=data.frame(eid=varSel$eid, t0=t0c, t1=t1c, age0=age0, age1=age1)
    colnames(dfTimeP)=c("eid", paste0(var, "_0"), paste0(var, "_1"), "age_0", "age_1")
   return(dfTimeP)
 } 

}
  

getlinearSlope=function(df, s){
  dfS=na.omit(subset(df, eid==s))
  dfS$phenoS=( dfS$pheno - baselineM) / baselineSD
  mWDF=as.data.frame(coef(summary(  lm(phenoS ~ age, data=dfS))))
  dfOut=data.frame(eid=s, slope=mWDF$Estimate[2])
  return(dfOut)
}


modelSlope=function(data, var){
  print(var)
  dfIn=data[[var]]
  df=subset(dfIn, is.na(dfIn[[paste0(var, "_1")]])==FALSE)
  tSum=length(which(!is.na(match(colnames(df), paste0(var, "_", 0:10)))))-1
  t=c(paste0("date_",var, "_", 0:tSum))
  a=c(paste0("age_",var, "_", 0:tSum))
  p=c(paste0(var, "_", 0:tSum))

  baselineM=mean(df[[paste0(var, "_", 0)]], na.rm=TRUE) #subset(means, pheno==var)$mean
  baselineSD=sd( df[[paste0(var, "_", 0)]], na.rm=TRUE) #subset(means, pheno==var)$sd

  varLong=data.frame(eid=rep(df$eid, length(t)),
                        time=do.call(rbind, lapply(0:tSum, function(x) data.frame(time=rep(x, NROW(df) )))),
                        pheno=do.call(rbind, lapply(p, function(x) data.frame(pheno=df[[x]]))),
                        age=do.call(rbind, lapply(a, function(x) data.frame(age=df[[x]]))),
                        date=do.call(rbind, lapply(t, function(x) data.frame(date=df[[x]]))))

  slopeL=lapply(levels(as.factor(varLong$eid)), function(x) getlinearSlope(df=varLong, s=x))
  slope=do.call(rbind, slopeL)
  colnames(slope)=c("eid", paste0(var, "_change"))
  dfOut=merge(subset(dfIn, select=c("eid",paste0(var, "_0"), paste0("age_",var, "_0"))), slope, by="eid", all.y=TRUE)
  return(dfOut)
}


checkCoding=function(varInc){
    print("Extract variable info and save on cluster")
    dataDic=fread(paste0(HOME, "/data/Data_Dictionary_Showcase.csv"))
    dataDicSel=subset(dataDic, FieldID %in% varInc$ID, select=c(FieldID, Field, Coding, Participants, ValueType, Instances, Notes, Link))
    print("Read in file containing coding info of the variables")
    codingDic=read.csv(paste0(HOME, "/data/Codings.csv"))
    
    print("Check coding of the variables")
    
    for ( i in 1:NROW(dataDicSel) ) {
      coding=dataDicSel$Coding[i]
      print(paste0("Check ", dataDicSel$Field[i]))
      print(paste0("Variable type ", dataDicSel$ValueType[i]))
      varType=c("Categorical multiple", "Categorical single")

      if( any(varType %in% dataDicSel$ValueType[i]) ){
          print("Get label for categorical variable")
          codingSel=subset(codingDic, Coding==coding)
          meaningPasted=paste0(codingSel$Value, ": ", codingSel$Meaning, " | ")
          dataDicSel$Label[i]=paste0(meaningPasted, collapse="")
      } else {
            dataDicSel$Label[i]=NA
      }
   print("Save output")
 write.csv(dataDicSel,
              file= paste0(HOME, "/data/dataDicExtracted.csv"),
              row.names = FALSE,
              col.names=T,
              quote=F)
} 
}





prepGWAdat=function(df, var, saveLabel, sample="", save="yes", relatives="yes"){
    
    library(readr)
    print("Prepare data for genome-wide analyses")
    print("Import PCA (genetic) data")
    sample_qc <- as.data.frame(data.table::fread(paste0(UKBB, "/phenotypes/ukb_sqc_v2.txt"), header = F, select = c(3:68)))
    colnames(sample_qc) <- c("array", "batch", "plate", "well", "cluster_CR", "dQC", "dna_concentration", "submitted_gender", "inferred_gender", "X_intensity", "Y_intensity", "submitted_plate", "submitted_well", "missing_rate", "heterozygosity", "heterozygosity_pc_corrected", "heterozygosity_missing_outlier", "PSCA", "in_kinship", "excluded_kinship_inference", "excess_relatives", "white_british", "pca_calculation", paste0("PC", seq(1,40)), "phasing_autosome", "phasing_X", "phasing_Y")
    # Sample eid; This corresponds to a list of all sample eids (with sex information), retrieved from a ".fam" file available from the UKBB portal
    sample_eid <- as.data.frame(data.table::fread(paste0(UKBB, "/genotypes/plink/_001_ukb_cal_chr1_v2.fam"), header = F, select = c(1,2,5), col.names = c("eid", "FID","SEX")))
    cov <- cbind(sample_eid, sample_qc)

    print("Merge with phenotype data")
    df$eid=as.character(df$eid)
    cov$eid=as.character(cov$eid)
    DFgen=merge(df, cov, by="eid", all.y=TRUE)

    #remove related samples (pca_calculation = 1)
    DFgen$removeRelated=ifelse(DFgen$pca_calculation==1, "include", NA) # used in PCA calculations
    DFgen$removeRelated=ifelse(DFgen$pca_calculation==0, "exclude", DFgen$removeRelated)

    print("Remove gender mismatch")
    DFgen$genderMismatch=ifelse(DFgen$submitted_gender==DFgen$inferred_gender, "include", NA) # 
    DFgen$genderMismatch=ifelse(DFgen$submitted_gender!=DFgen$inferred_gender, "exclude", DFgen$genderMismatch)
    
    print("remove non-european ancestry")
    DFgen$removeAncestry=ifelse(DFgen$white_british==1, "include", NA)
    DFgen$removeAncestry=ifelse(DFgen$white_british==0, "exclude", DFgen$removeAncestry)
    
    print("remove outliers in heterozygosity")
    DFgen$removeHetmiss=ifelse(DFgen$heterozygosity_missing_outlier==0, "include", NA)
    DFgen$removeHetmiss=ifelse(DFgen$heterozygosity_missing_outlier==1, "exclude", DFgen$removeHetmiss) #  (0/1)	(no/yes) Indicates samples identified as outliers in heterozygosity and missing rates, which indicates poor-quality genotypes 
    
    if(relatives=="yes"){
      DFgen$includeGeno=ifelse(DFgen$genderMismatch=="include" &
                        DFgen$removeAncestry=="include" &
                        DFgen$removeHetmiss=="include", "include", "exclude")
    }
    if(relatives=="no"){
      print("Remove relatives")
      DFgen$includeGeno=ifelse(DFgen$genderMismatch=="include" &
                        DFgen$removeAncestry=="include" &
                        DFgen$removeRelated=="include" &
                        DFgen$removeHetmiss=="include", "include", "exclude")
    }
             
    print("Recode batch for REGENIE")
    DFgen$batch=as.factor(DFgen$batch)
    DFgen$IID=DFgen$eid
    DFgenInc=subset(DFgen, includeGeno=="include")

    print(paste0("Number of individuals excluded after applying QC filter: ", NROW(DFgen)-NROW(DFgenInc)))
    print(paste0("Number of individuals included in genome-wide analyses: ", NROW(DFgenInc)))
    
    # Phenotype format
    #https://rgcgithub.github.io/regenie/options/#covariate-file-format
    #FID IID Y1 Y2
    #1 1 1.64818554321186 2.2765234736685
    #2 2 -2.67352013711554 -1.53680421614647
    #3 3 0.217542851471485 0.437289912695016
     

     print(paste0("Save phenotype file for ", saveLabel))
     GWAsave=subset(DFgenInc, select=c("FID", "IID",  paste0(var, "_", saveLabel )), IID > 0)

     if(saveLabel=="0"){
      print("Scale baseline phenotype")
      scaled=scale(GWAsave[[paste0(var, "_", saveLabel )]])
      GWAsave[[paste0(var, "_", saveLabel )]]=scaled
     }
     colnames(GWAsave)=c("FID", "IID",  paste0(var, "_", saveLabel , sample))

     if(save=="yes"){
      write.table(GWAsave, # AGE => removed age
              file= paste0(OUT, "/output/regenie/phenofiles/",  paste0(var, "_", saveLabel, sample )),
              sep="\t",
              row.names = FALSE,
              col.names=T,
              quote=F)


      print("Write covariate file")
      DFgen$batch=as.factor(as.numeric(DFgen$batch))
      DFgen$age=DFgen[[paste0("age_", var, "_0" )]] # age 
     # DFgen$age2=DFgen[[paste0("age_", var, "_0" )]]^2 # age^2

      if(saveLabel=="change"){
          print("Age not included as covariate")
          dfGWA=subset(DFgen, select=c("FID","IID","SEX", paste0("PC", 1:20), "batch"))
      }
      if(saveLabel=="0"){
          print("Age included as covariate")
          dfGWA=subset(DFgen, select=c("FID","IID","SEX", "age", paste0("PC", 1:20),"batch"))      
      }
  
  
      write.table(dfGWA, 
              file= paste0(OUT, "/output/regenie/phenofiles/",  paste0(var, "_", saveLabel, sample ), "_Covar"),
              sep="\t",
              row.names = FALSE,
              col.names=T,
              quote=F)
    } else {
      
      GWAsaveF=merge(GWAsave, cov, by="FID", all.x=TRUE)
      return(GWAsaveF)
    }

}



###################################################################################################
# ==================================== Process GWA results ========================================
###################################################################################################

readGWA=function(model){
       print("===================================")
       modelRaw=model
       model <- gsub("_interaction", "", modelRaw)
       print(paste0("Read model ", modelRaw))
       chr=22
       paths = paste0(OUT, "/output/regenie/genofiles/chr",1:chr, "_", model, "_", model, ".regenie") 
       timeCreated=file.info(paths)$ctime
       timeToday=Sys.time()
       timeDiff=round(difftime(timeToday, timeCreated, units = "days"),0)
       print(paste0("Days since GWAS finished: between ", min(timeDiff) ," and " ,max(timeDiff), " days" ))

       GWAList=lapply(paths, function(x) data.table::fread(x ))
       GWADF=do.call(rbind, GWAList)
       GWADF$MAF=ifelse(GWADF$A1FREQ>0.5, 1-GWADF$A1FREQ, GWADF$A1FREQ)
       GWADF=subset(GWADF, MAF>0.01)
       GWADF$P <- 10^(-GWADF$LOG10P)

       if (modelRaw==paste0(model, "_interaction")) {
            print("extract interaction based on between-individual age effects")
            GWADF=subset(GWADF, TEST=="ADD-INT_SNPxage") # test for the main effect of the interaction variable
       }
        if (modelRaw==paste0(model)) {
            print("extract marginal effects")
            GWADF=subset(GWADF, TEST=="ADD") # marginal test where the interacting variable has not been added as a covariate
       }

       GWADFsel=subset(GWADF, select=c(CHROM, GENPOS, ID, ALLELE1, ALLELE0, BETA, SE, P, INFO, N, A1FREQ, MAF))
       colnames(GWADFsel)=c("CHR", "POS", "SNP", "A1", "A2", "BETA", "SE", "P", "INFO", "N", "EAF", "MAF")
       print(head(GWADFsel, n=2))
       print(paste0("SNPs included in sumstats: ",  NROW(GWADF)))
       clump_input=subset(GWADFsel, P<5e-8)
       print(paste0(NROW(clump_input), " n gwas significant SNPs"))
      
       return(GWADFsel)
}


clumpData=function(clump_input, name){
  print(paste0("Clump data for ",name))
  print("Save SNP file")
  print(paste0("Number of included SNPs: ", NROW(clump_input)))

  write.table(clump_input,
              file= paste0(HOME, "/output/clump/",name,".pvalsfile"),
              sep="\t",
              row.names = FALSE,
              col.names=T,
              quote=F)

  system(paste0(PROGS, "/plink --bfile ", DATA, "/clump/g1000_eur --clump ", HOME, "/output/clump/",name, ".pvalsfile --clump-snp-field SNP --clump-field P --out ",  HOME, "/output/clump/" ,name, ".pvalsfile --clump-kb 250 --clump-r2 0.1 --clump-p1 1"))
  print("Remove SNP file")
  system(paste0("rm ", HOME,"/output/clump/",  name,".pvalsfile" ))

  if(file.exists(paste0(HOME, "/output/clump/", name,".pvalsfile.clumped"))==TRUE){
     clumped=read.table(paste0(HOME, "/output/clump/", name,".pvalsfile.clumped"), header = TRUE)
     clumpedSNPsMerged=subset(clump_input, SNP %in% unique(  clumped$SNP))
     clumpedSNPsMerged$pheno=name
     print(paste0("Numbr of SNPs included after clumping ", NROW(clumpedSNPsMerged) ) )
     system(paste0("rm ", HOME, "/output/clump/", name, ".pvalsfile.clumped" ))
    return(clumpedSNPsMerged)
  } else{
    return(data.frame(SNP=0))
  }
 
}


mungeData=function(name, path){
  paste0("start munging for: ", name)
  file.remove(paste0(path, name,".sumstats.gz"))
  setwd(paste0(OUT,"/output/ldsc"))

  munge(files = paste0(path, name),
        hm3 = paste0(DATA, "/ldsc/eur_w_ld_chr/w_hm3.snplist"),
        trait.names=name,
        info.filter = 0.9,
        maf.filter = 0.01)

    paste0("Munging done for: ", name)
}

saveGWA=function(df, name){
      print("Save as text file")
       data.table::fwrite(df, 
       paste0(OUT, "/output/regenie/genofiles/processed/",name),
       col.names = TRUE,
       row.names = FALSE, 
       quote = F, 
       sep = "\t")

}



###################################################################################################
# ==================================== Summarize GWA results ======================================
###################################################################################################


compareSNP=function(var, s=""){
    print(paste0("Combine clumped SNPs for ", var))
    betweenInteracAll=readRDS(paste0(HOME,"/output/clump/clump_",var,"_0_interaction.rds"))
    within=readRDS(paste0(HOME,"/output/clump/clump_",var,"_change.rds"))
    marginal=readRDS(paste0(HOME,"/output/clump/clump_",var,"_0.rds"))
    
    print("Select LD-independent SNPs from all identified SNPs")
    clumpDF=rbind(betweenInteracAll, within, marginal)
    clumpSel=clumpData(clump_input=clumpDF, name=paste0(var, "_eff"))

    SNPs=unique(clumpSel$SNP)
    keff=NROW(SNPs)

    if(NROW(SNPs)>0){
        betweenAllDF=data.table::fread(paste0(OUT, "/output/regenie/genofiles/processed/", var,"_0_interaction"))
        betweenAllS=subset(betweenAllDF, SNP %in% SNPs)
        betweenAllS$type="interaction_between"

        withinDF=data.table::fread(paste0(OUT, "/output/regenie/genofiles/processed/", var,"_change", s))
        withinS=subset(withinDF, SNP %in% SNPs)
        withinS$type=paste0("interaction_within", s)

        marginalDF=data.table::fread(paste0(OUT, "/output/regenie/genofiles/processed/", var,"_0"))
        marginalS=subset(marginalDF, SNP %in% SNPs)
        marginalS$type="marginal"

        dfOut=rbind(betweenAllS, withinS, marginalS)
        dfOut$p_fdr <- do.call(rbind, lapply(dfOut$P, function(x) p.adjust(x, "bonferroni",   keff  ) ))
        dfOut$keff=keff
        dfOut$pheno=var
        print(head(dfOut))
        return(dfOut)
    }
}



periodUKB=function(df, var, y){
    print(paste0("Get UKB long format for " , var))
    tSum=1
    dfSel=subset(df, select=c(paste0(var, "_",0:tSum ), paste0("date_", var, "_",0:tSum ), paste0("age_", var, "_",0:tSum )))

    listY=list()
    for (i in 0:tSum) {
        dfOut=na.omit(subset(df, select=c(paste0(var, "_", i ), paste0("date_", var, "_", i), paste0("age_", var, "_", i ))))
        colnames(dfOut)=c("pheno", "date", "age")
        dfOut$year=format(as.Date(dfOut$date), "%Y")
        dfOut$cohort=as.numeric(dfOut$year)-dfOut$age
        dfOut$date=NULL
        listY[[i+1]]=dfOut
    }
    dfL=do.call(rbind, listY)
    names(dfL)[names(dfL) == 'pheno'] <- var
    dfY=subset(dfL, year %in% y)
    print(paste0(NROW(dfY), " participants out of ",NROW(dfL)," included after applying year restriction"))
    return(dfY)
    }


getMean=function(df, a, type){
    df$var=df[[type]]
    dfA=subset(df, var==a)
    dfout=data.frame(mean=mean(dfA$pheno), n=NROW(dfA), bin=a, type=type)
    return(dfout)
}


getPeriod=function(df, var, y){
    print(y)
    dfS=na.omit(subset(df, select=c(var, "age", "year"), year==y))
    dfS$pheno=dfS[[var]]
    ageBin=40:120
    dfL=lapply(ageBin, function(x) getMean(df=dfS, a=x, type="age"))
    dfr=na.omit(do.call(rbind, dfL))
    if(NROW(dfr)==0){
     dfr= data.frame(mean=NA, n=NA , bin=NA,  type=NA)
    }
    dfr$year=y
    dfr$pheno=var
    dfr$type_x="period"
    dfr$meanTot=mean(dfr$mean, na.rm=TRUE)
    dfr$yearBin_n=length(levels(as.factor(dfr$bin)))

    return(dfr)
}


getCohort=function(df, var, y){
    print(y)
    dfS=na.omit(subset(df, select=c(var, "age", "cohort", "year"), cohort==y))
    dfS$pheno=dfS[[var]]
    yearBin=paste0("20", sprintf("%02d", 6:19))
    dfL=lapply(yearBin, function(x) getMean(df=dfS, a=x, type="year"))
    dfr=na.omit(do.call(rbind, dfL))
    if(NROW(dfr)==0){
     dfr= data.frame(mean=NA, n=NA , bin=NA,  type=NA)
    }
    dfr$bin_x=y
    dfr$pheno=var
    dfr$type_x="cohort"
    dfr$meanTot=mean(dfr$mean, na.rm=TRUE)
    return(dfr)
}


getAgebin=function(df, var, y){
    print(y)
    dfS=na.omit(subset(df, select=c(var, "age", "cohort"), age==y))
    dfS$pheno=dfS[[var]]
    cohortBin=1920:1980
    dfL=lapply(cohortBin, function(x) getMean(df=dfS, a=x, type="cohort"))
    dfr=na.omit(do.call(rbind, dfL))
    if(NROW(dfr)==0){
     dfr= data.frame(mean=NA, n=NA , bin=NA,  type=NA)
    }
    dfr$bin_x=y
    dfr$pheno=var
    dfr$type_x="age"
    dfr$meanTot=mean(dfr$mean, na.rm=TRUE)
    return(dfr)
}



readSNP=function(df, snp){
   print(snp)
   dfIn=subset(df, SNP==snp)[1,]
   chr=dfIn$CHR
   varID=paste0(dfIn$CHR, "_", dfIn$POS, "_", dfIn$A1, "_", dfIn$A2)
   bgen_file <- file.path(paste0(UKBB, "/genotypes/imp/v3/ukb22828_c",chr,"_b0_v3.bgen"))
   bgi_file <- file.path(paste0(UKBB, "/genotypes/imp/v3/ukb22828_c",chr,"_b0_v3.bgen.bgi"))
   sample_file <- file.path(UKBB, "genotypes/imp/v3/ukb22828_c1_22_b0_v3.sample")

   rds <- bigsnpr::snp_readBGEN(bgenfiles = bgen_file,
                                list_snp_id = list( varID), #  Each SNP ID should be in the form "<chr>_<pos>_<a1>_<a2>" (e.g. "1_88169_C_T" or "01_88169_C_T")
                                bgi_dir = dirname(bgen_file),
                                backingfile = tempfile())
   ukbb <- bigsnpr::snp_attach(rds)
   sample=subset(as.data.frame(data.table::fread(sample_file)), ID_1!=0)
   sample$allele=as.matrix(round(ukbb$genotypes[]))
   dfOut=subset(sample, select=c(ID_1, allele))

   dfOut <- dfOut %>% dplyr::mutate(allele = case_when(allele == 2 ~ 0, allele == 1 ~ 1, allele == 0 ~ 2 ))
   # Compute MAF from allele count table
   allele_counts <- table(dfOut$allele)  # Get allele frequency table
   n0 <- as.numeric(allele_counts["0"])  # Homozygous reference
   n1 <- as.numeric(allele_counts["1"])  # Heterozygous
   n2 <- as.numeric(allele_counts["2"])  # Homozygous alternate
   total_alleles <- 2 * sum(allele_counts)  
   effect_allele_count <- (n1 + 2 * n2)
   eaf <- effect_allele_count / total_alleles

   print(paste0("Estimated MAF versus REGENIE MAF: ", round(eaf, 1), " versus ",round(dfIn$EAF, 1)))
   colnames(dfOut)=c("eid", dfIn$SNP)
   return(as.data.frame(dfOut))
}








effectGene=function(dfC=NULL, dfL=NULL, s, pheno, type){
      
      if(  NROW(dfL) > 0){
      p=paste0(pheno, "_change")
      pDF=data.table::fread(paste0(OUT, "/output/regenie/phenofiles/", p))
      cDF=data.table::fread(paste0(OUT, "/output/regenie/phenofiles/", p, "_Covar"))
      dfS <- dfL %>% dplyr::group_by(ID) %>% dplyr::slice(1) %>% dplyr::ungroup() # take first observation
      gDF=na.omit(data.frame(IID=as.character(dfS$ID), gene=dfS$snp)) # unrelated individuals only
      pDF$FID=NULL; cDF$FID=NULL
      cDF$IID=as.character(cDF$IID); pDF$IID=as.character(pDF$IID)
      fDF=na.omit(Reduce(function(x,y) merge(x = x, y = y, by = "IID", all=TRUE ),  list(pDF, cDF, gDF)))
      n=NROW(fDF)
      }
    if(type=="interaction_between"){
        f=as.formula(paste0("pheno ~ snp + snp:age + SEX + age + as.factor(batch) + ",  paste0("PC", 1:20, collapse=" + ")))
        extract="snp:age"
        dfFUS=subset(dfC, select=c("pheno", "snp", "age", "SEX", "batch", paste0("PC", 1:20) ))
        n=NROW(dfFUS)
        fDF=dfC
    }
    if(type=="interaction_between_FU"){
        dfFU=subset(dfC, eid %in% fDF$IID)
        fDF=dfFU
        f=as.formula(paste0("pheno ~ snp + snp:age + SEX + age + as.factor(batch) + ",  paste0("PC", 1:20, collapse=" + ")))
        extract="snp:age"
    }
    if(type=="interaction_within"){
          f=as.formula(paste0(p, " ~ gene + SEX + as.factor(batch) + ",  paste0("PC", 1:20, collapse=" + ")))
          extract="gene"
    }

    m=lm(f, data=fDF)
    sumDF=as.data.frame(coef(summary(m)))
    sumDF$label=rownames(sumDF)
    sumS=subset(sumDF, label==extract)

    dfOut=data.frame(var=pheno, type=paste0("snp_",type, "_lm"), beta=sumS$Estimate, se=sumS[["Std. Error"]], n=n)
    return(dfOut)
}


getBetaGene=function(df, f, timepoint, type){
         dfIn=subset(df, t==timepoint)
         m=lm(f, data=dfIn)
         sumDF=as.data.frame(coef(summary(m)))
         sumDF$label=rownames(sumDF)
         sumS=subset(sumDF, label=="scale(gene)")
         medianAge=median(dfIn$age, na.rm=TRUE)
         dfOut=data.frame(BETA=sumS$Estimate, SE=sumS[["Std. Error"]], timepoint, type, medianAge=medianAge)
         return(dfOut)
    }




directionI=function(df, var, gene, snpList, type, snpOut=list()){
    for(i in 1:NROW(snpList)){
        snp=snpList[i]
        print(paste0("Check ", snp, " on ", var, " (",type,")"))
        f=as.formula(paste0("scale(",var, ") ~ scale(gene) + SEX + age + as.factor(batch) + ",  paste0("PC", 1:20, collapse=" + ")))

    if(type=="interaction_between"){
        split=median(df[[paste0("age_",var,"_0")]], na.rm=TRUE)
        p=paste0(var, "_0")
        dat=subset(df, select=c("eid", p, paste0("age_",p)))
        colnames(dat)=c("IID", var, "age")
        dat0=subset(dat, age<=split)
        dat0$t=0
        dat1=subset(dat, age>split)
        dat1$t=1
        datC=rbind(dat0, dat1)
    }
    if(type=="interaction_within"){
        p=paste0(var, "_change")
        dat=changeFunc(var=var, df=df, return="timepoints")
        dat0=subset(dat, select=c("eid", paste0(var, "_0"), "age_0"))
        colnames(dat0)=c("IID", var, "age")
        dat0$t=0
        dat1=subset(dat, select=c("eid", paste0(var, "_1"), "age_1"))
        colnames(dat1)=c("IID", var, "age")
        dat1$t=1
        datC=rbind(dat0, dat1)
    }

    pDF=data.table::fread(paste0(OUT, "/output/regenie/phenofiles/", p))
    cDF=data.table::fread(paste0(OUT, "/output/regenie/phenofiles/", p, "_Covar"))
    cDF$age=NULL
    gDF=data.frame(IID=as.character(gene$eid), gene=gene[[snp]])
    pDF$FID=NULL; cDF$FID=NULL
    cDF$IID=as.character(cDF$IID); pDF$IID=as.character(pDF$IID)
    fDF=na.omit(Reduce(function(x,y) merge(x = x, y = y, by = "IID", all=TRUE ),  list(datC, cDF, gDF)))
    fDFs=subset(fDF, IID %in% pDF$IID)

    betaL=lapply(c(0,1), function(x) getBetaGene(df=fDFs, timepoint=x, f=f, type=type ))
    betaDF=do.call(rbind, betaL)
    betaDF$direction=ifelse(abs(betaDF$BETA[1]>betaDF$BETA[2]), "intensification", "attenuation")
    betaDF$direction=ifelse(sign(betaDF$BETA[1])!=sign(betaDF$BETA[2]), "crossover",  betaDF$direction)
    betaDF$SNP=snp
    betaDF$pheno=var
    snpOut[[i]]=betaDF
    }
    dfOut=do.call(rbind, snpOut)
    return(dfOut)
}


normPlot=function(df, var){
  print(var)
  df$pheno=df[[paste0(var, "_change")]]

  if(NROW(levels(as.factor( df[[paste0(var, "_0")]])))<=2){
   print("++++binary variable++++")
   }
  if(NROW(levels(as.factor( df[[paste0(var, "_0")]])))>2){
  p1= ggplot(df, aes(x = pheno)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    stat_function(fun = dnorm, args = list(mean = mean(df$pheno), sd = sd(df$pheno)),
                  color = "red", size = 1) +
    labs(title = "", x = "", y = "Density") +
    theme_minimal()
  
  p2= ggplot(df, aes(sample = pheno)) +
    stat_qq() +
    stat_qq_line(color = "red") +
    labs(title = "", x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal()
  
  p=ggpubr::ggarrange(p1, p2) 
  p <- annotate_figure(p, top = text_grob(var, face = "bold", size = 14))
  return(p)
    }
}

