# install.packages("here")
library("here")
here()


#SELECT OPTIONS ======================
args = commandArgs(TRUE)
datdir=as.character(args[1]) #SCRIPT ANALYZES INDIVIDUAL-LEVEL DATA THAT IS RESTRICTED TO COMPLY WITH DATA USAGE AGREEMENTS (see data availability section in manuscript)
datdir="/Users/jamesroe/LCBC/Users/jamesroe/PHD_project/Paper3/data/" #NB! DEL
saveres=0 #1/0
# brainvar="area" #select metric
brainvar="thickness"
#====================================#


# LOAD PACKAGES -----------------------------------------------------------
tmp.packages = c("tidyverse","readr","devtools","magrittr","viridis","corrplot","gamm4","DescTools")
tmpnew.packages = tmp.packages[!(tmp.packages %in% installed.packages()[,"Package"])]
if(length(tmpnew.packages)) {
  install.packages(tmpnew.packages)
}
sapply(tmp.packages, require, character.only = T)
rm(list=ls(pattern="tmp*"))
# devtools::install_github("LKremer/ggpointdensity")
library("ggpointdensity")


#---data
loadData = function(datdir, brainvar) {
  
  DF1 = read.csv(paste0(datdir, "PopAsym_LCBC_", brainvar,".csv"), stringsAsFactors = F, header = T, sep=",")
  DF1 %<>% filter(Age > 18)
  DF2 = read.csv(paste0(datdir, "PopAsym_UKB_", brainvar,"_structuralCovar.csv"), stringsAsFactors = F, header = T, sep="\t")
  DF3 = read.csv(paste0(datdir, "PopAsym_HCP_", brainvar,"_structuralCovar.csv"), stringsAsFactors = F, header = T, sep="\t")
  DF1$Subject = DF1$fsid_base
  
  
  #harmonize vars across cohorts
  DF1$eAge = DF1$Age
  DF2$Age = DF2$eAge
  if (brainvar == "thickness") {
    DF2$cluster.thickness.6p.merge111.label=DF2$cluster.thickness.6P.merge111.label
    DF2$cluster.thickness.6P.merge111.label=NULL
  }
  
  #add ICV
  ep="/Users/jamesroe/LCBC/Users/jamesroe/PHD_project/Paper2/AgeSym/HDD_encrypted"
  B = read.csv(file.path(ep,paste0("Statfiles/aseg.long.stats.LCBC.txt")), stringsAsFactors = F, header = T, sep =" ")
  B$Folder=substr(B$Measure.volume,1,19)
  DF1 = merge(DF1, B)
  FolderLink=DF1$Folder
  
  DF1$ICV=DF1$EstimatedTotalIntraCranialVol
  DF2$EstimatedTotalIntraCranialVol=DF2$ICV
  DF3$Age=DF3$eAge
  DF2$Folder=NA
  DF3$Folder=NA
  # # length(
  #   intersect(names(DF1),names(DF2))
  #   # )
  # # length(
  #   intersect(names(DF1),names(DF3))
  #   # )
  # # length(
  #   intersect(names(DF2),names(DF3))
  #   # )
    
  sitenamelink=DF1$Site_Name
  DF1$Sex
  DF2$Sex = ifelse(DF2$sex == 0, "female", "male")
  DF3$Sex = ifelse(DF3$Sex == "F", "female", "male")
  DF1 %<>% select(one_of(intersect(names(DF1),names(DF2))))
  DF2 %<>% select(one_of(intersect(names(DF1),names(DF2))))
  DF3 %<>% select(one_of(intersect(names(DF1),names(DF2))))
  
  #harmonized
  DF = rbind(DF1 %>% mutate(cohort="LCBC", Site_Name = sitenamelink),
             DF2 %>% mutate(cohort="UKB", Site_Name="UKB",Folder=NA),
             DF3 %>% mutate(cohort="HCP", Site_Name="HCP",Folder=NA))
  return(DF)
}
DF = loadData(datdir, brainvar)


#---clusters
mergelabs = names(DF %>% select(contains("merge")))
names(DF)[(names(DF)) %in% mergelabs] = toupper(mergelabs)
mergelabs = names(DF %>% select(contains("MERGE")))
labs = names(DF %>% select(contains(".label")))


#calc mean across left and right dominant rois
nlabs=labs[grepl("6n",labs)]
if (brainvar=="thickness"){
  plabs=labs[grepl("6p",labs)] [1:11]
} else {
  plabs=labs[grepl("6p",labs)]
}
DF$meanN=rowMeans(DF[,c(nlabs)])
DF$meanP=rowMeans(DF[,c(plabs)])


#weighted average by vtx size
if (brainvar=="area") {
  vtxsizeN=c(7525,6414,2775,1516,1757,599,325)
  vtxsizeP=c(11032,4679,626,1050,1218,607,350)
  DF$wMeanN = as.vector((
    (DF$cluster.area.6n.0001.label * vtxsizeN[1]) +
      (DF$cluster.area.6n.0002.label * vtxsizeN[2]) +
      (DF$cluster.area.6n.0003.label * vtxsizeN[3]) +
      (DF$cluster.area.6n.0004.label * vtxsizeN[4]) +
      (DF$cluster.area.6n.0005.label * vtxsizeN[5]) +
      (DF$cluster.area.6n.0006.label * vtxsizeN[6]) +
      (DF$cluster.area.6n.0007.label * vtxsizeN[6])
  ) / sum(vtxsizeN))
  DF$wMeanP = as.vector((
    (DF$cluster.area.6p.0001.label * vtxsizeP[1]) +
      (DF$cluster.area.6p.0002.label * vtxsizeP[2]) +
      (DF$cluster.area.6p.0003.label * vtxsizeP[3]) +
      (DF$cluster.area.6p.0004.label* vtxsizeP[4]) +
      (DF$cluster.area.6p.0005.label * vtxsizeP[5]) +
      (DF$cluster.area.6p.0006.label * vtxsizeP[6]) +
      (DF$cluster.area.6p.0007.label * vtxsizeP[7])
  ) / sum(vtxsizeP))
  
} else {
  vtxsizeN=c(8379,5117,2290,2218,1762,865,1173,916,653)
  vtxsizeP=c(4110,1581,2441,1485,1135,725,537,595,425,416,779)
  DF$wMeanN = as.vector((
    (DF$cluster.thickness.6n.0001.label * vtxsizeN[1]) +
      (DF$cluster.thickness.6n.0002.label * vtxsizeN[2]) +
      (DF$cluster.thickness.6n.0003.label * vtxsizeN[3]) +
      (DF$cluster.thickness.6n.0004.label * vtxsizeN[4]) +
      (DF$cluster.thickness.6n.0005.label * vtxsizeN[5]) +
      (DF$cluster.thickness.6n.0006.label * vtxsizeN[6]) +
      (DF$cluster.thickness.6n.0007.label * vtxsizeN[7]) +
      (DF$cluster.thickness.6n.0007.label * vtxsizeN[8]) +
      (DF$cluster.thickness.6n.0007.label * vtxsizeN[9])
  ) / sum(vtxsizeN))
  DF$wMeanP = as.vector((
    (DF$cluster.thickness.6p.0001.label * vtxsizeP[1]) +
      (DF$cluster.thickness.6p.0002.label * vtxsizeP[2]) +
      (DF$cluster.thickness.6p.0003.label * vtxsizeP[3]) +
      (DF$cluster.thickness.6p.0004.label* vtxsizeP[4]) +
      (DF$cluster.thickness.6p.0005.label * vtxsizeP[5]) +
      (DF$cluster.thickness.6p.0006.label * vtxsizeP[6]) +
      (DF$cluster.thickness.6p.0007.label * vtxsizeP[7]) +
      (DF$cluster.thickness.6p.0008.label * vtxsizeP[8]) +
      (DF$cluster.thickness.6p.0009.label * vtxsizeP[9]) +
      (DF$cluster.thickness.6p.0010.label * vtxsizeP[10]) +
      (DF$cluster.thickness.6p.0011.label * vtxsizeP[11])
  ) / sum(vtxsizeP))
}
meanlabs = names(DF %>% select(contains("mean")))


folderlink=DF$Folder[DF$hemi==1]
sublink=DF$Subject[DF$hemi==1]
cohortlink=DF$cohort[DF$hemi==1]
sitenamelink=DF$Site_Name[DF$hemi==1]
sexlink=DF$Sex[DF$hemi==1]
LL=DF %>% filter(hemi==1) %>% select(c(
  Subject,
  all_of(nlabs),
  all_of(plabs),
  all_of(meanlabs),
  all_of(mergelabs)
))
RR=DF %>% filter(hemi==0) %>% select(c(
  Subject,
  all_of(nlabs),
  all_of(plabs),
  all_of(meanlabs),
  all_of(mergelabs)
))
identical(LL$Subject,RR$Subject)
LL$Subject = RR$Subject = NULL


#CALC ASYM + INVERSE AI's IN RIGHTWARD ROIS
ASY = (LL - RR) / ((LL + RR) / 2)
ASY[,c(nlabs)]=ASY[,c(nlabs)]*-1


ASY$cohort=cohortlink
ASY$Subject=sublink
ASY$Folder=folderlink
# ASY = rbind(ASY %>% filter(cohort=="LCBC"),
#             ASY %>% filter(cohort=="UKB"),
#             ASY %>% filter(cohort=="HCP"))
# ggplot(ASY,aes(wMeanP,wMeanN*-1,col=as.factor(cohort))) + geom_point(alpha=0.3) + geom_smooth(method="lm")
# plot(ASY$wMeanP)
# plot(ASY$wMeanN)
# identify(ASY$meanN,ASY$meanP,col="pink")


#check outliers based on weighted means across all leftward and rightward clusters
outlierthresh = 6
QL = scale(c(ASY$wMeanP))
QR = scale(c(ASY$wMeanN))
db.out = data.frame(Subject=sublink,QL,QR,
                    QQL=NA,QQR=NA)
db.out$QQL [(abs(QL)>outlierthresh)] = QL[(abs(QL)>outlierthresh)]
db.out$QQR [(abs(QR)>outlierthresh)] = QR[(abs(QR)>outlierthresh)]
outliersL = which(!is.na(db.out$QQL))
outliersR = which(!is.na(db.out$QQR))
outliers = c(db.out[outliersL,1], db.out[outliersR[!(outliersR %in% outliersL)],1]) #outliersL + unique outliersR
cat("\noutlier thresh =", outlierthresh, "removing", length(outliers), "cases\n\n")
print(outliers)
db.out[outliersL,]
db.out[outliersR,]


#remove 2 outliers in HCP for all subsequent analyses
# if (cohort == "HCP") {
  outliershcp = c(894067, 132118)
  ASY = ASY[-which(sublink %in% outliershcp),]
  DF = DF[-which(DF$Subject %in% outliershcp),] #MIGHT HAVE CHANGED SOMETHING
  sitenamelink = sitenamelink[-which(sublink %in% outliershcp)]
  sexlink = sexlink[-which(sublink %in% outliershcp)]
  sublink = sublink[-which(sublink %in% outliershcp)]
# }


ASY %<>% select(Subject,everything())


#all now correct
ASY %>% filter(cohort=="LCBC") %>% select(Subject) %>% unique() %>% nrow()
ASY %>% filter(cohort=="UKB") %>% select(Subject) %>% unique() %>% nrow()
ASY %>% filter(cohort=="HCP") %>% select(Subject) %>% unique() %>% nrow()
ASY %>% filter(cohort=="LCBC") %>% lm(wMeanP ~ wMeanN*-1,. ) %>% summary()
ASY %>% filter(cohort=="HCP") %>% lm(wMeanP ~ wMeanN*-1,. ) %>% summary()


#set up to remove scanner, sex and age effects
extradat=DF %>% filter(hemi == 1) %>% 
  select(Subject, cohort, Sex, Site_Name, eAge, EstimatedTotalIntraCranialVol, Folder) %>%  mutate(ICV = scale(EstimatedTotalIntraCranialVol, T, T))
  
#make cohort ids for merging data
extradat$Folder[extradat$cohort == "UKB"] = paste0(extradat$Subject[extradat$cohort == "UKB"], "_", "UKB")
extradat$Folder[extradat$cohort == "HCP"] = paste0(extradat$Subject[extradat$cohort == "HCP"], "_", "HCP")
ASY$Folder[ASY$cohort == "UKB"] = paste0(ASY$Subject[ASY$cohort == "UKB"], "_", "UKB")
ASY$Folder[ASY$cohort == "HCP"] = paste0(ASY$Subject[ASY$cohort == "HCP"], "_", "HCP")
names(extradat)[1:2]=c("SUB","COH")
ASY = left_join(ASY, extradat, by = "Folder")



# CORRECT FOR AGE SEX SCANNER (ICV) ---------------------------------------------
#lmm
correctCovars = function(ASY, icvcor) {
  c = 1
  ASYcor = ASY
  
  print(paste("returning labels corrected for sex and scanner and ICV =", icvcor))
  
  for (l in c(nlabs, plabs, mergelabs, meanlabs)) {
    c = c + 1
    ASY$Y = ASY[[l]]
    
    if (icvcor == F) {
      
      mmr = summary(lmer(Y ~ eAge + Sex + Site_Name  + (1 |
                                                        Subject), data = ASY))
    
      # ASYcor[, c] = resid(lmer(Y ~ eAge + Sex + Site_Name  + (1 |
      #                                             Subject), data = ASY))
      
      #manually residualise by fixed effects
      coo = mmr$coefficients[, 1]
      ASYcor[, c] = ASY$Y - (ASY$eAge * coo[2]) #correct age
      ASYcor[, c] = ifelse(ASY$Sex == "male", ASYcor[, c] - (1 * coo[3]), ASYcor[, c]) #correct sex
    
      ASYcor[, c] = ifelse(ASY$Site_Name == "ousAvanto", ASYcor[, c] - (1 * coo[4]),
                           ifelse(ASY$Site_Name == "ousPrisma", ASYcor[, c] - (1 * coo[5]),
                                  ifelse(ASY$Site_Name == "ousSkyra", ASYcor[, c] - (1 * coo[6]),
                                         ifelse(ASY$Site_Name == "UKB", ASYcor[, c] - (1 * coo[7]), ASYcor[, c])
                 ))
      ) #correct scanner
    } else if (icvcor == T) {
      
      mmr = summary(lmer(Y ~ eAge + Sex + Site_Name  + ICV + (1 |
                                                                Subject), data = ASY))
      
      # ASYcor[, c] = resid(lmer(Y ~ eAge + Sex + Site_Name  + (1 |
      #                                             Subject), data = ASY))
      
      #manually residualise by fixed effects
      coo = mmr$coefficients[, 1]
      ASYcor[, c] = ASY$Y - (ASY$eAge * coo[2]) #correct age
      ASYcor[, c] = ifelse(ASY$Sex == "male", ASYcor[, c] - (1 * coo[3]), ASYcor[, c]) #correct sex
      
      ASYcor[, c] = ifelse(ASY$Site_Name == "ousAvanto", ASYcor[, c] - (1 * coo[4]),
                              ifelse(ASY$Site_Name == "ousPrisma", ASYcor[, c] - (1 * coo[5]),
                                     ifelse(ASY$Site_Name == "ousSkyra", ASYcor[, c] - (1 * coo[6]),
                                            ifelse(ASY$Site_Name == "UKB", ASYcor[, c] - (1 * coo[7]), ASYcor[, c])
                                     ))
      ) #correct scanner
      ASYcor[, c] = ASYcor[, c] - (ASYcor$ICV * coo[8])
    }
  }
  return(ASYcor)
}
ASYcor = correctCovars(ASY, icvcor = F)
ASYcorICV = correctCovars(ASY, icvcor = T)


# CORRELATION MATRICES WITH AND WITHOUT ICV-ASSOCIATED VARIANCE REMOVED ----------------------------------------------------
ASYcor1 = ASYcor %>% filter(cohort == "LCBC")
ASYcor2 = ASYcor %>% filter(cohort == "UKB")
ASYcor3 = ASYcor %>% filter(cohort == "HCP")

cormat1 = ASYcor1 %>% select(nlabs,plabs)
cormat1 = cor(cormat1,method="pearson")
cormat2 = ASYcor2 %>% select(nlabs,plabs)
cormat2 = cor(cormat2,method="pearson")
cormat3 = ASYcor3 %>% select(nlabs,plabs)
cormat3 = cor(cormat3,method="pearson")


if (saveres == 1) {
  if (! dir.exists("results/structuralCovar")) {
    dir.create("results/structuralCovar")
  }
  #non-annotated
  png(
    filename = here("results/structuralCovar",
                    paste("cormat", "LCBC", brainvar, "correct-publish.png", sep = ".")
    ),
    width = 10,
    height = 10,
    units = "cm",
    res = 300
  )
  corrplot(cormat1,method="color",type="full",col=colorRampPalette(c("blue","white","orange2"))(25),tl.pos="n",tl.col = "blue",cl.cex = 0.8,cl.pos = "r",tl.cex = 0.8,number.cex = 0.8)
  dev.off()
  #non-annotated
  png(
    filename = here("results/structuralCovar",
                    paste("cormat", "UKB", brainvar, "correct-publish.png", sep = ".")
    ),
    width = 10,
    height = 10,
    units = "cm",
    res = 300
  )
  corrplot(cormat2,method="color",type="full",col=colorRampPalette(c("blue","white","orange2"))(25),tl.pos="n",tl.col = "blue",cl.cex = 0.8,cl.pos = "r",tl.cex = 0.8,number.cex = 0.8)
  dev.off()
  #non-annotated
  png(
    filename = here("results/structuralCovar",
      paste("cormat", "HCP", brainvar, "correct-publish.png", sep = ".")
    ),
    width = 10,
    height = 10,
    units = "cm",
    res = 300
  )
  corrplot(cormat3,method="color",type="full",col=colorRampPalette(c("blue","white","orange2"))(25),tl.pos="n",tl.col = "blue",cl.cex = 0.8,cl.pos = "r",tl.cex = 0.8,number.cex = 0.8)
  dev.off()
  save('cormat1','cormat2','cormat3',file = here("results/structuralCovar",paste0("cormat.crossCohorts.",brainvar,".Rda")))
  
  
  #annotated
  png(filename = here("results/structuralCovar",
      paste("cormat",brainvar, "LCBC", "correctaddcoef-publish.png", sep = ".")
    ),
    width = 20,
    height = 20,
    units = "cm",
    res = 300
  )
  corrplot(cormat1,method="color",type="full",addCoef.col = "black",col=colorRampPalette(c("blue","white","orange2"))(25),tl.pos="n",tl.col = "blue",cl.cex = 0.8,cl.pos = "r",tl.cex = 0.8,number.cex = 0.8)
  dev.off()
  #annotated
  png(filename = here("results/structuralCovar",
                      paste("cormat",brainvar, "UKB", "correctaddcoef-publish.png", sep = ".")
  ),
  width = 20,
  height = 20,
  units = "cm",
  res = 300
  )
  corrplot(cormat2,method="color",type="full",addCoef.col = "black",col=colorRampPalette(c("blue","white","orange2"))(25),tl.pos="n",tl.col = "blue",cl.cex = 0.8,cl.pos = "r",tl.cex = 0.8,number.cex = 0.8)
  dev.off()
  #annotated
  png(filename = here("results/structuralCovar",
                      paste("cormat",brainvar, "HCP", "correctaddcoef-publish.png", sep = ".")
  ),
  width = 20,
  height = 20,
  units = "cm",
  res = 300
  )
  corrplot(cormat3,method="color",type="full",addCoef.col = "black",col=colorRampPalette(c("blue","white","orange2"))(25),tl.pos="n",tl.col = "blue",cl.cex = 0.8,cl.pos = "r",tl.cex = 0.8,number.cex = 0.8)
  dev.off()
}


if (brainvar == "area") {
  print(
    paste0(
      "leftward asymmetry in SMG/perisylvian (#1L) was related to higher rightward asymmetry in inferior parietal cortex (#2R; r = ",
      cormat1[8, 2],
      "[LCBC])"
    ))
  print(
    paste0(
      "leftward anterior cingulate asymmetry (ACC; #3L) was related to higher rightward asymmetry in mPFC (#6R, r = ",cormat1[10, 6], ")"
    ))
  print(
    paste0(
      "and leftward asymmetry in a superior frontal cluster (#7L) was related to rightward asymmetry in the cingulate (#3R, r = ", cormat1[14, 3], ")"
    ))
}


#---additionaly remove ICV-associated variance
# ASYcorICV = correctCovars(ASY, icvcor = T)


ASYcorICV1 = ASYcorICV %>% filter(cohort == "LCBC")
ASYcorICV2 = ASYcorICV %>% filter(cohort == "UKB")
ASYcorICV3 = ASYcorICV %>% filter(cohort == "HCP")

cormatICV1 = ASYcorICV1 %>% select(nlabs,plabs)
cormatICV1 = cor(cormatICV1,method="pearson")
cormatICV2 = ASYcorICV2 %>% select(nlabs,plabs)
cormatICV2 = cor(cormatICV2,method="pearson")
cormatICV3 = ASYcorICV3 %>% select(nlabs,plabs)
cormatICV3 = cor(cormatICV3,method="pearson")

max(abs(
  c(
    as.vector(cormat1) - as.vector(cormatICV1),
    as.vector(cormat2) - as.vector(cormatICV2),
    as.vector(cormat3) - as.vector(cormatICV3)
  )
))
        

#---save results
if (saveres == 1) {
  save('cormatICV1','cormatICV2','cormatICV3',file = here("results/structuralCovar",paste0("cormat.crossCohorts.ICVcor",brainvar,".Rda")))
}



# GEODESIC DISTANCE -------------------------------------------------------
#load surfdist matrix
surfdist = read.csv(
  here("results/structuralCovar", paste0("geodesicMeanDistmat_",brainvar,".txt")),stringsAsFactors=F, header=T,row.names=1)

#subset matrices based on metric
if (brainvar=="area"){
  #right ipsilateral
  ripsidist = surfdist[1:7, 1:7] %>% unname()
  
  #left ipsilateral
  lipsidist = surfdist[8:(nrow(surfdist)), 8:(nrow(surfdist))] %>% unname()
  
  #contralateral
  contradist = surfdist[8:(nrow(surfdist)), 1:7] %>% unname()
  
  #subsets of cormat
  contracor1 = cormat1[8:(nrow(surfdist)), 1:7]
  lipsicor1 = cormat1[8:(nrow(surfdist)), 8:(nrow(surfdist))]
  ripsicor1 = cormat1[1:7, 1:7]
  
  contracor2 = cormat2[8:(nrow(surfdist)), 1:7]
  lipsicor2 = cormat2[8:(nrow(surfdist)), 8:(nrow(surfdist))]
  ripsicor2 = cormat2[1:7, 1:7]
  
  contracor3 = cormat3[8:(nrow(surfdist)), 1:7]
  lipsicor3 = cormat3[8:(nrow(surfdist)), 8:(nrow(surfdist))]
  ripsicor3 = cormat3[1:7, 1:7]
} else {
  #right ipsilateral
  ripsidist = surfdist[1:9, 1:9] %>% unname()
  
  #left ipsilateral
  lipsidist = surfdist[10:(nrow(surfdist)), 10:(nrow(surfdist))] %>% unname()
  
  #contralateral
  contradist = surfdist[10:(nrow(surfdist)), 1:9] %>% unname()
  
  #subsets of cormat
  contracor1 = cormat1[10:(nrow(surfdist)), 1:9]
  lipsicor1 = cormat1[10:(nrow(surfdist)), 10:(nrow(surfdist))]
  ripsicor1 = cormat1[1:9, 1:9]
  
  contracor2 = cormat2[10:(nrow(surfdist)), 1:9]
  lipsicor2 = cormat2[10:(nrow(surfdist)), 10:(nrow(surfdist))]
  ripsicor2 = cormat2[1:9, 1:9]
  
  contracor3 = cormat3[10:(nrow(surfdist)), 1:9]
  lipsicor3 = cormat3[10:(nrow(surfdist)), 10:(nrow(surfdist))]
  ripsicor3 = cormat3[1:9, 1:9]
}

contracor1 = FisherZ(rho=contracor1) #fishers R to Z transform
contracor2 = FisherZ(rho=contracor2)
contracor3 = FisherZ(rho=contracor3)
row.names(contradist)=NULL; row.names(lipsidist)=NULL; row.names(ripsidist)=NULL


#correlate opposite direction asymmetries with surfdist
Y = as.vector(contracor1)
X = as.vector(as.matrix(contradist))
if (brainvar == "area") {
  print("geodesic distance was lower between cluster-pairs that were more correlated")
  #LCBC
  print(cor.test(X,Y,method="spearman"))
  #UKB
  Y = as.vector(contracor2)
  print(cor.test(X,Y,method="spearman"))
  #HCP
  Y = as.vector(as.matrix(contracor3))
  print(cor.test(X,Y,method="spearman"))
}


#correlate same direction asymmetries with surfdist
#left
if (brainvar == "area") {
  print("same-direction SA asymmetries were not more correlated if closer in cortex (leftward [all p > .5]; rightward [all p > .5])")
  #LCBC
  lipsidistZ1 = FisherZ(rho=lipsicor1[lower.tri(lipsicor1)]) #fishers R to Z transform
  ripsidistZ1 = FisherZ(rho=ripsicor1[lower.tri(ripsicor1)])
  #UKB
  lipsidistZ2 = FisherZ(rho=lipsicor2[lower.tri(lipsicor1)])
  ripsidistZ2 = FisherZ(rho=ripsicor2[lower.tri(ripsicor1)])
  #HCP
  lipsidistZ3 = FisherZ(rho=lipsicor3[lower.tri(lipsicor1)])
  ripsidistZ3 = FisherZ(rho=ripsicor3[lower.tri(ripsicor1)])

  cor.test(lipsidistZ1, lipsidist[lower.tri(lipsidist)])
  cor.test(lipsidistZ2, lipsidist[lower.tri(lipsidist)])
  cor.test(lipsidistZ3, lipsidist[lower.tri(lipsidist)])

  cor.test(ripsidistZ1, ripsidist[lower.tri(ripsidist)])
  cor.test(ripsidistZ2, ripsidist[lower.tri(ripsidist)])
  cor.test(ripsidistZ3, ripsidist[lower.tri(ripsidist)])
}


mytheme = theme(
  plot.background = element_rect(fill = "white"),
  panel.background = element_rect(fill = "white"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  text = element_text(color = "black", size = 14, family = "Helvetica Neue Light"),
  plot.title = element_text(hjust = 0.5),
  axis.title.y = element_text(color = "black", size = 16, vjust =-1, margin = margin(0,20,0,0)),
  axis.title.x = element_text(color = "black", size = 16, vjust = -2, margin = margin(0,20,10,0)),
  axis.text = element_text(color = "black", size = 14),
  legend.text = element_text(color = "white", size = 14),
  legend.position = "none")


#plot opposite hemisphere geodesic in UKB
Y = as.vector(contracor2)
ggd=data.frame(Y,X)
distplot=ggplot(data=ggd,aes(x=X,Y,col=Y)) +
  geom_point(size=5) +
  geom_smooth(method="lm",col="black") +
  theme_classic() +
  scale_y_continuous(breaks = seq(-.2, .8, by = .3)) +
  xlab("Geodesic distance\n(cross-hemisphere)") +
  ylab("r (Z)") +
  mytheme

if (saveres == 1 && brainvar == "area") {
  ggsave(filename = here("results/structuralCovar",
                       paste("surfdist", "UKB", brainvar, "2Fishers.png", sep = ".")), plot = distplot, width = 8, height=8, dpi=600, units="cm")
}


# POST HOC PCA FOR THICKNESS IN UKB DATA ----------------------------------
if (brainvar == "thickness") {
  
  cohort="UKB"
  if (cohort == "LCBC") { dat=ASYcor1 }
  if (cohort == "UKB") { dat=ASYcor2 }
  if (cohort == "HCP") { dat=ASYcor3 }
  
  #main PCA - all leftward v all rightward
  pcaASYcor=prcomp(dat[,c(nlabs,plabs)],center=T,scale.=T)
  spcaASYcor=summary(pcaASYcor)
  spcaASYcor
  
  #scree plot
  pcaASY.var=spcaASYcor$sdev^2
  pcaASY.var.per=round(pcaASY.var/sum(pcaASY.var)*100,1)
  
  #first 10 components
  gdat=data.frame(var=pcaASY.var.per,
                pc=seq(1,20))[-c(11:20),]
  

  scree = ggplot(data=gdat,aes(x=pc,var)) +
    geom_point(col="#00886e",size=3) +
    geom_bar(stat = "identity",fill="#00886e",alpha=0.5) +
    geom_line(col="#00886e", size=1) +
    theme_classic() +
    ggtitle(cohort) +
    xlab("Principal component") +
    ylab("% explained variance") +
    coord_cartesian(ylim = c(4,23)) +
    mytheme +
    theme(
      axis.ticks = element_blank(),
      axis.text.x = element_blank()
    )
  
  if (cohort != "UKB") {
    scree = scree + coord_cartesian(ylim = c(4,9))
  }
  
  saveres=0
  if (saveres == 1 && brainvar == "thickness") {
    ggsave(filename = here("results/structuralCovar",
                           paste("scree", cohort, brainvar, ".png", sep = ".")), plot = scree, width = 8, height=8, dpi=300, units="cm")
  }

  
  #UKB plot - thickness
  #mean extracted spatially averaging across vertices in merged labels
  plot(ASYcor2$CLUSTER.THICKNESS.6N.MERGE19.LABEL*-1,ASYcor2$CLUSTER.THICKNESS.6P.MERGE111.LABEL,col="pink")
  cor.test(ASYcor2$CLUSTER.THICKNESS.6N.MERGE19.LABEL*-1,ASYcor2$CLUSTER.THICKNESS.6P.MERGE111.LABEL)
  
  
  #weighted mean
  plot(ASYcor1$wMeanN*-1,ASYcor1$wMeanP)
  abline(lm(ASYcor1$wMeanN*-1 ~ASYcor1$wMeanP))
  plot(ASYcor2$wMeanN*-1,ASYcor2$wMeanP)
  abline(lm(ASYcor2$wMeanN*-1 ~ASYcor2$wMeanP))
  plot(ASYcor3$wMeanN*-1,ASYcor3$wMeanP)
  abline(lm(ASYcor3$wMeanN*-1 ~ASYcor3$wMeanP))
  
  cor.test(ASYcor1$wMeanP,ASYcor1$wMeanN*-1)
  cor.test(ASYcor2$wMeanP,ASYcor2$wMeanN*-1)
  cor.test(ASYcor3$wMeanP,ASYcor3$wMeanN*-1)
  
  
  #non weighted mean
  cor.test(ASYcor1$meanN*-1,ASYcor1$meanP)
  cor.test(ASYcor2$meanN*-1,ASYcor2$meanP)
  cor.test(ASYcor3$meanN*-1,ASYcor3$meanP)
  
  
  #PCA by hemi to triple check this relationship
  #UKB
  pcaASYcorL = prcomp(ASYcor2[, c(plabs)], center = T, scale. = T)
  pcaASYcorR = prcomp(ASYcor2[, c(nlabs)], center = T, scale. = T)
  plot(pcaASYcorL$x[,1],pcaASYcorR$x[,1])
  abline(lm(pcaASYcorR$x[,1]~pcaASYcorL$x[,1]), col="blue")
  cor.test(pcaASYcorL$x[,1],pcaASYcorR$x[,1])
  
  #LCBC
  pcaASYcorL = prcomp(ASYcor1[, c(plabs)], center = T, scale. = T)
  pcaASYcorR = prcomp(ASYcor1[, c(nlabs)], center = T, scale. = T)
  plot(pcaASYcorL$x[,1],pcaASYcorR$x[,1])
  abline(lm(pcaASYcorR$x[,1]~pcaASYcorL$x[,1]), col="blue")
  cor.test(pcaASYcorL$x[,1],pcaASYcorR$x[,1])
  
  pcaASYcorL = prcomp(ASYcor3[, c(plabs)], center = T, scale. = T)
  pcaASYcorR = prcomp(ASYcor3[, c(nlabs)], center = T, scale. = T)
  plot(pcaASYcorL$x[,1],pcaASYcorR$x[,1])
  abline(lm(pcaASYcorR$x[,1]~pcaASYcorL$x[,1]), col="blue")
  cor.test(pcaASYcorL$x[,1],pcaASYcorR$x[,1])
  

  #UKB plot
  (pmeanthick = ggplot(data=ASYcor2,aes(x=wMeanN*-1, y=wMeanP)) +
    geom_pointdensity(size=0.75) + #col="dark green") +
    scale_color_viridis(option = "G") +
    geom_vline(xintercept = 0,col="grey",linetype=5,size=0.5) +
    geom_hline(yintercept = 0,col="grey",linetype=5,size=0.5) +
    geom_smooth(method="lm",col="black",se = T,size=0.25) +
    theme_classic() +
    mytheme +
    theme(
  axis.ticks = element_blank()
        )
  )
  #change axis lim for LCBC and HCP
  if (cohort != "UKB") {
    pmeanthick = pmeanthick + theme(axis.text = element_text(size = 11)) #+ coord_cartesian(ylim = c(-0.1,0.15), xlim = c(-0.12,0.12))
  }
  

  if (saveres == 1 && brainvar == "thickness") {
    ggsave(filename = here("results/structuralCovar",
                           paste("weightaverageplotUKB",brainvar, "png", sep = ".")), plot = pmeanthick, width = 8, height=8, dpi=600, units="cm")
  }
}
