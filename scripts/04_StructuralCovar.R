# install.packages("here")
library("here")
here()


#SELECT OPTIONS ======================
args = commandArgs(TRUE)
datdir=as.character(args[1]) #SCRIPT ANALYZES INDIVIDUAL-LEVEL DATA THAT IS RESTRICTED TO COMPLY WITH DATA USAGE AGREEMENTS (see data availability section in manuscript)
datdir="/Users/jamesroe/LCBC/Users/jamesroe/PHD_project/Paper3/data/" #NB! DEL
saveres=0 #1/0
# cohort="LCBC" #select cohort
# cohort="HCP"
cohort="UKB"
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
# DF = read.csv(paste0(datdir, "PopAsym_",cohort,"_", brainvar,"_structuralCovar.csv"), stringsAsFactors = F, header = T, sep="\t")
#DEL
#LCBC
DF = read.csv(paste0(datdir, "PopAsym_",cohort,"_", brainvar,"_structuralCovar.csv"), stringsAsFactors = F, header = T, sep=",")
DF$Subject=DF$fsid_base
DF %<>% filter(Age >=55)

#DEL
DF %<>% filter(eAge <=55)

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


sublink=DF$Subject[DF$hemi==1]
LL=DF %>% filter(hemi==1) %>% select(c(
  all_of(nlabs),
  all_of(plabs),
  all_of(meanlabs),
  all_of(mergelabs)
))
RR=DF %>% filter(hemi==0) %>% select(c(
  all_of(nlabs),
  all_of(plabs),
  all_of(meanlabs),
  all_of(mergelabs)
))


#CALC ASYM + INVERSE AI's IN RIGHTWARD ROIS
ASY = (LL - RR) / ((LL + RR) / 2)
ASY[,c(nlabs)]=ASY[,c(nlabs)]*-1
#DEL
# ASY[,c(nlabs)]=ASY[,c(nlabs)]

#two large outliers in HCP thickness data detected (SI Fig 8)
plot(ASY$wMeanP,ASY$wMeanN*-1,col="pink")
abline(lm(ASY$wMeanN*-1~ASY$wMeanP))
plot(ASY$wMeanP)
plot(ASY$wMeanN)
# identify(ASY$meanN,ASY$meanP,col="pink")


#outliers based on weighted means across all leftward and rightward clusters
#both detected
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


#remove them for all subsequent analyses
if (cohort == "HCP") {
  # outliersvisual = c(894067, 132118)
  ASY = ASY[-which(sublink %in% outliers),]
  DF = DF[-which(DF$Subject %in% outliers),]
  sublink = sublink[-which(sublink %in% outliers)]
}


ASY$Subject = sublink
ASY %<>% select(Subject,everything())


#set up to remove scanner, sex and age effects
if (cohort=="UKB") {
  extradat = DF %>% filter(hemi == 1) %>% 
    select(Subject, eAge, sex, ICV) %>% mutate(ICV = scale(ICV, T, T))
  
} else if (cohort == "LCBC") {
  #DEL
  ep="/Users/jamesroe/LCBC/Users/jamesroe/PHD_project/Paper2/AgeSym/HDD_encrypted"
  B = read.csv(file.path(ep,paste0("Statfiles/aseg.long.stats.",cohort,".txt")), stringsAsFactors = F, header = T, sep =" ")
  B$Folder=substr(B$Measure.volume,1,19)
  DF = merge(DF, B)
  
  extradat = DF %>% filter(hemi == 1) %>% 
    mutate(ICV = scale(EstimatedTotalIntraCranialVol, T, T)) %>% 
    select(Subject, Age, Sex, Site_Name, fsid_base, ICV) 
  #harmonize some vars across cohorts
  extradat$eAge=extradat$Age; extradat$Age=NULL
  extradat$sex=extradat$Sex; extradat$Sex=NULL

} else if (cohort == "HCP") {
  extradat = DF %>% filter(hemi == 1) %>% select(Subject, eAge, sex, ICV) %>% mutate(ICV = scale(ICV, T, T))
}
ASY=merge(ASY,extradat,by="Subject")



# CORRECT FOR AGE SEX SCANNER ---------------------------------------------
#lmm if LCBC
#lm if not
c = 1
ASYcor = ASY
if (cohort != "LCBC") {
  for (l in c(nlabs, plabs, mergelabs, meanlabs)) {
    c = c + 1 #set counter at 2 as first col is subject
    ASY$Y = ASY[[l]]
    
    #manually residualise by fixed effects
    mmr = summary(lm(Y ~ eAge + sex, data = ASY))
    coo = mmr$coefficients[, 1]
    ASYcor[, c] = ASY$Y - (ASY$eAge * coo[2]) #correct age
    ASYcor[, c] = ifelse(ASY$sex == 1, ASYcor[, c] - (1 * coo[3]), ASYcor[, c]) #correct sex
    
  }
} else {
  for (l in c(nlabs, plabs, mergelabs, meanlabs)) {
    c = c + 1
    ASY$Y = ASY[[l]]
    mmr = summary(lmer(Y ~ eAge + sex + Site_Name  + (1 |
                                                        fsid_base), data = ASY))
    
    #manually residualise by fixed effects
    coo = mmr$coefficients[, 1]
    ASYcor[, c] = ASY$Y - (ASY$eAge * coo[2]) #correct age
    ASYcor[, c] = ifelse(ASY$sex == "male", ASYcor[, c] - (1 * coo[3]), ASYcor[, c]) #correct sex
    ASYcor[, c] = ifelse(
      ASY$Site_Name == "ousPrisma",
      ASYcor[, c] - (1 * coo[4]),
    
      ifelse(ASY$Site_Name == "ousSkyra", ASYcor[, c] -
               (1 * coo[5]),
             ASYcor[, c])
    ) #correct scanner
  }
}


# CORRELATION MATRICES WITH AND WITHOUT ICV-ASSOCIATED VARIANCE REMOVED ----------------------------------------------------
cormat = ASYcor %>% select(nlabs,plabs)
cormat = cor(cormat,method="pearson")
if (saveres == 1) {
  if (! dir.exists("results/structuralCovar")) {
    dir.create("results/structuralCovar")
  }
  #non-annotated
  png(
    filename = here("results/structuralCovar",
      paste("cormat", cohort, brainvar, "publish.png", sep = ".")
    ),
    width = 10,
    height = 10,
    units = "cm",
    res = 300
  )
  corrplot(cormat,method="color",type="full",col=colorRampPalette(c("blue","white","orange2"))(25),tl.pos="n",tl.col = "blue",cl.cex = 0.8,cl.pos = "r",tl.cex = 0.8,number.cex = 0.8)
  dev.off()
  #annotated
  png(filename = here("results/structuralCovar",
      paste("cormat", cohort, brainvar, "addcoef-publish.png", sep = ".")
    ),
    width = 20,
    height = 20,
    units = "cm",
    res = 300
  )
  corrplot(cormat,method="color",type="full",addCoef.col = "black",col=colorRampPalette(c("blue","white","orange2"))(25),tl.pos="n",tl.col = "blue",cl.cex = 0.8,cl.pos = "r",tl.cex = 0.8,number.cex = 0.8)
  dev.off()
}

if (cohort == "LCBC" && brainvar == "area") {
  print(
    paste0(
      "leftward asymmetry in SMG/perisylvian (#1L) was related to higher rightward asymmetry in inferior parietal cortex (#2R; r = ",
      cormat[8, 2],
      "[LCBC])"
    ))
  print(
    paste0(
      "leftward anterior cingulate asymmetry (ACC; #3L) was related to higher rightward asymmetry in mPFC (#6R, r = ",cormat[10, 6], ")"
    ))
  print(
    paste0(
      "and leftward asymmetry in a superior frontal cluster (#7L) was related to rightward asymmetry in the cingulate (#3R, r = ", cormat[14, 3], ")"
    ))
}


#---additionaly remove ICV-associated variance
c = 1
ASYcorICV = ASY
if (cohort != "LCBC") {
  for (l in c(nlabs, plabs, mergelabs, meanlabs)) {
    c = c + 1 #set counter at 2 as first col is subject
    ASY$Y = ASY[[l]]
    ASYcorICV[, c] = resid(lm(ASY$Y ~ ASY$eAge + ASY$sex + ASY$ICV))
  }
} else {
  for (l in c(nlabs, plabs, mergelabs, meanlabs)) {
    c = c + 1
    ASY$Y = ASY[[l]]
    mmr = summary(lmer(Y ~ eAge + sex + Site_Name + ICV  + (1 |
                                                              fsid_base), data = ASY))
    #manually residualise by fixed effects
    coo = mmr$coefficients[, 1]
    ASYcorICV[, c] = ASY$Y - (ASY$eAge * coo[2]) #correct age
    ASYcorICV[, c] = ifelse(ASY$sex == "male", ASYcorICV[, c] - (1 * coo[3]), ASYcorICV[, c]) #correct sex
    ASYcorICV[, c] = ifelse( 
      ASY$Site_Name == "ousPrisma",
      ASYcorICV[, c] - (1 * coo[4]),
      
      ifelse(
        ASY$Site_Name == "ousSkyra",
        ASYcorICV[, c] - (1 * coo[5]),
        ASYcorICV[, c]
      ) #correct scanner
    )
    ASYcorICV[, c] = ASYcorICV[, c] - (ASY$ICV * coo[6]) #correct ICV
  }
}

cormatICV = ASYcorICV %>% select(nlabs,plabs)
cormatICV = cor(cormatICV,method="pearson")
max(abs(cormat-cormatICV))
    

#---save results
if (saveres == 1) {
  save('cormat','cormatICV',file = here("results/structuralCovar",paste0("cormat.",cohort,"life.",brainvar,"-publish.Rda")))
}
cor(as.vector(cormat), as.vector(cormatICV))



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
  contracor = cormat[8:(nrow(surfdist)), 1:7]
  lipsicor = cormat[8:(nrow(surfdist)), 8:(nrow(surfdist))]
  ripsicor = cormat[1:7, 1:7]
} else {
  #right ipsilateral
  ripsidist = surfdist[1:9, 1:9] %>% unname()
  
  #left ipsilateral
  lipsidist = surfdist[10:(nrow(surfdist)), 10:(nrow(surfdist))] %>% unname()
  
  #contralateral
  contradist = surfdist[10:(nrow(surfdist)), 1:9] %>% unname()
  
  #subsets of cormat
  contracor = cormat[10:(nrow(surfdist)), 1:9]
  lipsicor = cormat[10:(nrow(surfdist)), 10:(nrow(surfdist))]
  ripsicor = cormat[1:9, 1:9]
}

contracor = FisherZ(rho=contracor) #fishers R to Z transform
row.names(contradist)=NULL; row.names(lipsidist)=NULL; row.names(ripsidist)=NULL


#correlate opposite direction asymmetries with surfdist
Y = as.vector(contracor)
X = as.vector(as.matrix(contradist))
if (brainvar == "area") {
  print("geodesic distance was lower between cluster-pairs that were more correlated")
} else if (brainvar == "thickness") {
  print("Opposite-direction CT asymmetries that were closer in cortex were more negatively correlated in LCBC (rho = .29, p = .003) but not HCP (p = .32) or UKB (p = .84)")
}
cor.test(X,Y,method="spearman")


#correlate same direction asymmetries with surfdist
#left
if (brainvar == "area") {
  print("same-direction SA asymmetries were not more correlated if closer in cortex (leftward [all p > .5]; rightward [all p > .5])")
} else if (brainvar == "thickness") {
  print(paste("whereas CT asymmetry in left-asymmetric (rho = -.40, p = .003 [LCBC]; rho = -.44, p = 8.1-4 [UKB], rho = -.28, p = .04 [HCP])",
              "and right-asymmetric (rho = -.34, p = .04 [LCBC]; rho = -.48, p = .003 [UKB], rho = -.58, p = 2.0-4 [HCP])", 
              "regions was more positively correlated"))
}
#left
lipsidistZ = FisherZ(rho=lipsicor[lower.tri(lipsicor)]) #fishers R to Z transform
ripsidistZ = FisherZ(rho=ripsicor[lower.tri(ripsicor)]) #fishers R to Z transform
# cor.test(lipsidist[lower.tri(lipsidist)], 
#          lipsicor[lower.tri(lipsicor)],method="spearman")
cor.test(lipsidistZ, lipsidist[lower.tri(lipsidist)])
#right
# cor.test(ripsidist[lower.tri(ripsidist)],
#          ripsicor[lower.tri(ripsicor)],method="spearman")
cor.test(ripsidistZ, ripsidist[lower.tri(ripsidist)])
#FINDSTRING

  
plot(lipsidist[lower.tri(lipsidist)],
     lipsicor[lower.tri(lipsicor)])
abline(lm(lipsicor[lower.tri(lipsicor)]~
            lipsidist[lower.tri(lipsidist)]))
plot(ripsidist[lower.tri(ripsidist)],
     ripsicor[lower.tri(ripsicor)])
abline(lm(ripsicor[lower.tri(ripsicor)]~
            ripsidist[lower.tri(ripsidist)]))


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


ggd=data.frame(Y,X)
distplot=ggplot(data=ggd,aes(x=X,Y,col=Y)) +
  geom_point(size=5) +
  geom_smooth(method="lm",col="black") +
  theme_classic() +
  scale_y_continuous(breaks = seq(-.2, .8, by = .3)) +
  xlab("Geodesic distance\n(cross-hemisphere)") +
  ylab("r (Z)") +
  mytheme

if (saveres == 1 && cohort == "UKB" && brainvar == "area") {
  ggsave(filename = here("results/structuralCovar",
                       paste("surfdist", cohort, brainvar, "2Fishers.png", sep = ".")), plot = distplot, width = 8, height=8, dpi=600, units="cm")
}



# POST HOC PCA FOR THICKNESS IN UKB DATA ----------------------------------
if (brainvar == "thickness") {
  
  #main PCA - all leftward v all rightward
  pcaASYcor=prcomp(ASYcor[,c(nlabs,plabs)],center=T,scale.=T)
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
      axis.ticks = element_blank()
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
  plot(ASYcor$CLUSTER.THICKNESS.6N.MERGE19.LABEL*-1,ASYcor$CLUSTER.THICKNESS.6P.MERGE111.LABEL,col="pink")
  cor.test(ASYcor$CLUSTER.THICKNESS.6N.MERGE19.LABEL*-1,ASYcor$CLUSTER.THICKNESS.6P.MERGE111.LABEL)
  if (cohort == "UKB") {
    print("r = -0.56, p = 2.2e-16 (UKB)")
  } else if (cohort == "LCBC") {
    print("r = -0.05, p = .04 (LCBC)")
  } else if (cohort == "HCP") {
    print("r = -0.05, p = .072 (HCP)")
  }
  
  #non weighted mean
  plot(ASYcor$meanN*-1,ASYcor$meanP)
  cor.test(ASYcor$meanN*-1,ASYcor$meanP)
  if (cohort == "UKB") {
    print("r = -0.56, p = 2.2e-16 (UKB)")
  } else if (cohort == "LCBC") {
    print("r = -0.05, p = .03 (LCBC)")
  } else if (cohort == "HCP") {
    print("r = -0.04, p = .15 (HCP)")
  }
  
  #weighted mean
  plot(ASYcor$wMeanN*-1,ASYcor$wMeanP)
  cor.test(ASYcor$wMeanP,ASYcor$wMeanN*-1)
  if (cohort == "UKB") { 
    print("r = -0.61, p = 2.2e-16 (UKB)")
  } else if (cohort == "LCBC") {
    print("r = -0.10, p = .00014 (LCBC)")
  } else if (cohort == "HCP") {
    print("r = -0.112, p = .00016 (HCP)")
  }
  
  #PCA by hemi to triple check this relationship
  pcaASYcorL = prcomp(ASYcor[, c(plabs)], center = T, scale. = T)
  pcaASYcorR = prcomp(ASYcor[, c(nlabs)], center = T, scale. = T)
  plot(pcaASYcorL$x[,1],pcaASYcorR$x[,1])
  cor.test(pcaASYcorL$x[,1],pcaASYcorR$x[,1])
  if (cohort == "UKB") {
    print("r = 0.66, p = 2.2e-16 (UKB)")
  } else if (cohort == "LCBC") {
    print("r = -0.17, p = 2.00e-11 (LCBC)")
  } else if (cohort == "HCP") {
    print("r = -0.13, p = 2.31e-5 (HCP)")
  }
  
  ASYcor$pcaASYcorL = pcaASYcorL$x[,1]
  ASYcor$pcaASYcorR = pcaASYcorR$x[,1]

  #UKB plot
  (pmeanthick = ggplot(data=ASYcor,aes(x=wMeanN*-1, y=wMeanP)) +
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
    pmeanthick = pmeanthick + coord_cartesian(ylim = c(-0.05,0.15), xlim = c(-0.05,0.12))
  }
  
  saveres=1
  if (saveres == 1 && brainvar == "thickness") {
    ggsave(filename = here("results/structuralCovar",
                           paste("weightaverageplot", cohort, "55minus",brainvar, "png", sep = ".")), plot = pmeanthick, width = 8, height=8, dpi=600, units="cm")
  }
}
