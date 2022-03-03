rm(list=ls())
library("here")
here()


#SELECT OPTIONS ======================
# args = commandArgs(TRUE)
# datdir=as.character(args[1]) #SCRIPT ANALYZES INDIVIDUAL-LEVEL DATA THAT IS RESTRICTED TO COMPLY WITH DATA USAGE AGREEMENTS (see data availability section in manuscript)
saveres=0 #1/0
cohort="UKB"
# brainvar="area" #select metric
brainvar="thickness"
#====================================#


#---load packages
tmp.packages = c("tidyverse","readr","devtools","itsadug","magrittr","ggridges","gamm4","viridis") #, "missMDA", "naniar", "VIM")
tmpnew.packages = tmp.packages[!(tmp.packages %in% installed.packages()[,"Package"])]
if(length(tmpnew.packages)) {
  install.packages(tmpnew.packages)
}
sapply(tmp.packages, require, character.only = T)
rm(list=ls(pattern="tmp*"))

# devtools::install_github("LKremer/ggpointdensity")
library("ggpointdensity")


#---data
DF = read.csv(paste0(datdir,"/","PopAsym_",cohort,"_", brainvar,"_structuralCovar.csv"), stringsAsFactors = F, header = T, sep="\t")


#fix global labs
mergelabs = names(DF %>% select(contains("merge")))
names(DF)[(names(DF)) %in% mergelabs] = toupper(mergelabs)
mergelabs = names(DF %>% select(contains("MERGE")))
labs = names(DF %>% select(contains(".label")))
nlabs=labs[grepl("6n",labs)]
if (brainvar == "thickness") {
  plabs=labs[grepl("6p",labs)][1:11]
} else {
  plabs=labs[grepl("6p",labs)]
}
labs = c(nlabs,plabs)


LL = DF %>% filter(hemi==1) %>% select(c(nlabs,plabs))
RR = DF %>% filter(hemi==0) %>% select(c(nlabs,plabs))
ASY = (LL-RR)/( (LL+RR)/2 )
L = DF %>% filter(hemi==1) #add base data
ASY = data.frame(L[,1:13],
               ASY)


#inverse rightward AIs - easier interpretation for lm
ASY[,c(nlabs)]=ASY[,c(nlabs)]*-1


#remove NAs on handedness & ambidexterous
unique(ASY$handedness)
length(ASY$handedness)
ind1 = which(is.na(ASY$handedness)) #NA
ind2 = which(ASY$handedness == -3) #missing
ind3 = which(ASY$handedness == 3) #ambidex
ind = c(ind1, ind2, ind3)
ASY = ASY[-ind, ]
nrow(ASY)
unique(ASY$handedness)
ASY$handedness = ifelse(ASY$handedness == 2, 1, 0) #recode left handers = 1, right handers = 0


# HANDEDNESS, SEX, ICV REGRESSIONS ----------------------------------------
for (i in 1:length(labs)) {
  if (i==1) {
    asyout=pp=pr=sp=m1out=m3out=list()
    phand=psex=picv=0
    ehand=esex=eicv=0
    cihand=cisex=ciicv=0
  }
  
  c=labs[i]
  print(c)
  ASY$Y=ASY[[c]]
  
  #unstandardized Y for plotting
  m1 = lm(
    ASY$Y ~ ASY$handedness +
    + scale(ASY$eAge, center = T, scale = F)
    + scale(ASY$sex, center=T, scale = F) +
      scale(ASY$ICV, center = T, scale = F)
  )
  m1out[[i]]=summary(m1)$coefficients
  
  
  #standardized scores
  m3 = lm(
    scale(ASY$Y, center=T,scale=T) ~ ASY$handedness
    + scale(ASY$eAge, center = T, scale = F)
    + ASY$sex +
      scale(ASY$ICV, center = T, scale = T)
  )
  m3out[[i]]=summary(m3)$coefficients
  
  
  #take standardized stats
  #pvals
  phand[i]=m3out[[i]][2,4]
  psex[i]=m3out[[i]][4,4]
  picv[i]=m3out[[i]][5,4]
  
  #effects
  ehand[i]=m3out[[i]][2,1]
  esex[i]=m3out[[i]][4,1]
  eicv[i]=m3out[[i]][5,1]
  
  #ci's
  cihand[i]=m3out[[i]][2,2]*1.96
  cisex[i]=m3out[[i]][4,2]*1.96
  ciicv[i]=m3out[[i]][5,2]*1.96
  
  ASY$pred1=predict(m1, newdata = data.frame(
    Y=ASY$Y,
    handedness=ASY$handedness,
    eAge=0, #mean age
    sex=0, #mean sex
    ICV = 0)) #mean ICV
  
  
  ASY$partial_resid=ASY$pred1+residuals(m1)
  pr[[i]]=ggplot(ASY) +
    geom_jitter(aes(x=handedness,y=partial_resid,col=as.factor(handedness)),alpha=0.07) +
    geom_violin(aes(x=handedness,y=partial_resid,fill=as.factor(handedness)),alpha=0) +
    geom_boxplot(aes(x=handedness,y=partial_resid,col=as.factor(handedness)),width=.1) +
    theme_classic()
  
}


# CREATE RESULTS TABLE ----------------------------------------------------
if (brainvar=="area") {
  res_phand=as.data.frame(phand) %>% mutate(meas="handedness",
                                         n=seq(labs),
                                         effect=ehand,
                                         ci=cihand,
                                         labs=labs,
                                         dir=c(rep("rightward",7),rep("leftward",7)))
  res_psex=as.data.frame(psex) %>% mutate(meas="sex",
                                       n=seq(labs),
                                       effect=esex,
                                       ci=cisex,
                                       labs=labs,
                                       dir=c(rep("rightward",7),rep("leftward",7)))
  res_picv=as.data.frame(picv) %>% mutate(meas="icv",
                                       n=seq(labs),
                                       effect=eicv,
                                       ci=ciicv,
                                       labs=labs,
                                       dir=c(rep("rightward",7),rep("leftward",7)))
} else {
  res_phand=as.data.frame(phand) %>% mutate(meas="handedness",
                                         n=seq(labs),
                                         effect=ehand,
                                         ci=cihand,
                                         labs=labs,
                                         dir=c(rep("rightward",9),rep("leftward",11)))
  res_psex=as.data.frame(psex) %>% mutate(meas="sex",
                                       n=seq(labs),
                                       effect=esex,
                                       ci=cisex,
                                       labs=labs,
                                       dir=c(rep("rightward",9),rep("leftward",11)))
  res_picv=as.data.frame(picv) %>% mutate(meas="icv",
                                       n=seq(labs),
                                       effect=eicv,
                                       ci=ciicv,
                                       labs=labs,
                                       dir=c(rep("rightward",9),rep("leftward",11)))
  
}
names(res_psex)[1]="p"
names(res_phand)[1]="p"
names(res_picv)[1]="p"
res1 = rbind(res_phand,res_psex,res_picv)
res1$N = as.integer(row.names(res1))


#multiple comparison log
#area = 14 rois * 4 assoc (hand, sex, icv, cog)
#cth = 20 rois * 4 assoc (hand, sex, icv, cog)
bonfarea = 14 * 4
bonfcth = 20 * 4
bonf = bonfarea + bonfcth

#quick overview
tmp = res1[ (-log10(res1$p)) >= (-log10(0.01/bonf)),]


# ASSOCIATION PLOT --------------------------------------------------------
mytheme = theme(
  plot.background = element_rect(fill = "white"),
  panel.background = element_rect(fill = "white"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  text = element_text(color = "black", size = 12, family = "Helvetica Neue Light"),
  plot.title = element_text(hjust = 0.5),
  axis.ticks = element_blank(),
  axis.line.x=element_blank(),
  axis.text = element_text(color = "black", size = 12),
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  legend.text = element_text(color = "white", size = 12))

#isolate only surviving
tmp = res1[ (-log10(res1$p)) >= (-log10(0.01/bonf)),]
p_pval = ggplot(res1) +
  # geom_point(aes(x=N,y=-log10(p),col=factor(meas)),size=2) +
  geom_point(data=tmp,aes(x=N,y=-log10(p),col=factor(meas)),size=2) +
  geom_bar(aes(x=N,y=-log10(p),col=factor(meas)),stat="identity", width = 0.01) +
  geom_hline(yintercept = -log10(0.01)) +
  geom_hline(yintercept = -log10(0.01/bonf),col="black",linetype=2) +
  scale_color_manual(values=wesanderson::wes_palette(n=5, name="Cavalcanti1")[c(1:4,5)]) +
  scale_fill_manual(values=wesanderson::wes_palette(n=5, name="Cavalcanti1")[c(1:4,5)]) +
  theme_classic() +
  mytheme


#effects
eff_survive = res1[ which(-log10(res1$p) > -log10(0.01/bonf)) ,]

#save only effect plot for thickness. Areal one saved later after adding cog effect
if (brainvar=="thickness") {
  p_effect = ggplot(eff_survive) +
    geom_point(aes(x=N,y=effect,col=factor(meas)),size=2) +
    geom_bar(aes(x=N,y=effect,col=factor(meas)),stat="identity", width = 0.01,alpha=0.05) +
    geom_errorbar(aes(x=N,y=effect,ymin=effect-ci,ymax=effect+ci,col=factor(meas)),width=1) +
    geom_hline(yintercept = 0,col="black",linetype=2) +
    scale_color_manual(values=wesanderson::wes_palette(n=5, name="Cavalcanti1")[c(1:4,5)]) +
    scale_color_manual(values=wesanderson::wes_palette(n=5, name="Cavalcanti1")[c(1:4,5)]) +
    coord_cartesian(xlim = c(0,max(eff_survive$N))) +
    ylab("Standardized β") +
    theme_classic() +
    mytheme +
    theme(legend.position = "none")
  
  p_both = cowplot::plot_grid(p_pval,p_effect,ncol=2)
  if (saveres == 1) {
    if (!dir.exists(here("results/indDiffs"))) {
      dir.create(here("results/indDiffs"))
    }
    ggsave(filename = here(paste0("results/indDiffs/p_effects_",brainvar,".png")), plot = p_both, width = 18, height=6.5, dpi=600, units="cm")
  }
}


# BUILD OUTPUT TABLE ------------------------------------------------------
tmp1 = res1 %>% filter(meas=="handedness") %>% select(labs,effect,p)
tmp2 = res1 %>% filter(meas == "sex") %>% select(labs, effect, p)
tmp3 = res1 %>% filter(meas == "icv") %>% select(labs, effect, p)
TABLE1 = data.frame(tmp1, tmp2, tmp3)
TABLE1 %<>% select(-4, -7)
names(TABLE1)[2] = "hand.effect"
names(TABLE1)[4] = "sex.effect"
names(TABLE1)[6] = "icv.effect"
TABLE1$hand.effect = round(TABLE1$hand.effect, 3)
TABLE1$sex.effect = round(TABLE1$sex.effect, 3)
TABLE1$icv.effect = round(TABLE1$icv.effect, 3)

if (brainvar=="thickness"){
  TABLE1$roi = substr(TABLE1$labs,24,25)
  roiname=c("R_Superior-temporal-sulcus",
                   "R_Lateral-occipital",
                   "R_Lingual-gyrus",
                   "R_Posterior-insula_Sylvian",
                   "R_Entorhinal",
                   "R_Superior-insula",
                   "R_Planum-temporale",
                   "R_Posterior-cingulate",
                   "R_Anterior-insula",
                   "L_Cingulate",
                   "L_Superior-frontal-gyrus",
                   "L_Postcentral",
                   "L_Precentral",
                   " L_Supplementary-motor-cortex",
                   "L_Collateral-sulcus",
                   "L_Anterior-transverse-collateral-sulcus",
                   " L_Caudal-middle-frontal",
                   "L_Caudal-superior-frontal",
                   "L_Rostral-superior-frontal",
                   "L_Calcarine-sulcus ")
  TABLE1$roiname = roiname
  #rearrange table rows
  TABLE1 = TABLE1[c(10:20,1:9),]
} else {
  TABLE1$roi = substr(TABLE1$labs,19,20)
  roiname=c("R_Parieto-occipital_sulcus",
                   "R_Inferior-parietal_Lateral-occipital",
                   "R_Cingulate",
                   "R_Middle-frontal-gyrus",
                   "R_Superior-temporal-sulcus",
                   "R_Superior-frontal-gyrus",
                   "R_Gyrus-rectus",
                   "L_Postcentral-gyrus_Supramarginal",
                   "L_Anterior-temporal_Parahippocampal-gyrus",
                   "L_Anterior-cingulate_Subcallosal",
                   "L_Anterior-insula",
                   "L_Retrosplenial-cortex",
                   "L_Temporal-pole_Inferior-temporal-gyrus",
                   "L_Superior-frontal-gyrus")
  TABLE1$roiname = roiname
  #rearrange table rows
  TABLE1 = TABLE1[c(8:14,1:7),]
}
TABLE1$p = formatC(TABLE1$p, format = "e", digits = 2)
TABLE1$p.1 = formatC(TABLE1$p.1, format = "e", digits = 2)
TABLE1$p.2 = formatC(TABLE1$p.2, format = "e", digits = 2)
TABLE1 %<>% select(labs,roi,roiname,everything()) 




# COGNITION ---------------------------------------------------------------
#load UKB cognitive
UKBCOG = read.csv(file.path(datdir,"UKB_cognition.tsv"), stringsAsFactors = F, header = T, sep=" ")
UKBCOG %<>% select(-2,-3,-4,-5,-6,-7,-8,-9,-"eAge",-"ebirthDate")
UKBCOG = merge(ASY,UKBCOG, by = "eid")


#list cogntive vars
cogvars = names(UKBCOG %>% select(contains("_2")))

#rename / transform / recode vars
UKBCOG %<>% mutate(WordPairs = assoc_c20197_2_0, #
                   IncorrectMatches = log(mean_rtmatch_c20023_2_0), #
                   ProspMem = prospmem_c20018_2_0, #
                   FluidIQ = fluidiq_c20016_2_0, #
                   SymbolDigit = symbol_c23324_2_0, #
                   PairMatch = log(pairsincorrectinv_c399_2_1), #
                   TMTA = log(tmta_c6348_2_0),
                   TMTB = log(tmtb_c6350_2_0),
                   NumericMem = nummem_c4282_2_0,
                   LondonTower = tower_c21004_2_0,
                   Matrices = matrix_c6373_2_0,
                   ProspMem = ifelse(ProspMem==1,1,
                                     ifelse(is.na(ProspMem),NA,0)), 
                   TMTA =  ifelse(TMTA==0,NA,TMTA),
                   TMTB =  ifelse(TMTB==0,NA,TMTB),
                   NumericMem = ifelse(NumericMem==-1,NA,NumericMem),
                   PairMatch = ifelse(PairMatch==-Inf,NA,PairMatch),
                   TMTA = ifelse(TMTA==-Inf,NA,TMTA),
                   TMTB = ifelse(TMTB==-Inf,NA,TMTB)
)

# explore cogvars before and after transformation
explorecogvars = function(cogvars) {
  pg=pc=list()
  for (i in 1:length(cogvars)) {
    pg[[i]] = ggplot(UKBCOG, aes_string(cogvars[i], bins = 100)) +
      geom_histogram() +
      theme_minimal()
  }
  library("patchwork")
  cogvarplot = pg[[1]]+pg[[2]]+pg[3]+pg[[4]]+pg[[5]]+pg[[6]]+pg[[7]]+pg[[8]]+pg[[9]]+pg[[10]]+pg[[11]]
  return(cogvarplot)
}
explorecogvars(cogvars)
cogvars=c("WordPairs","IncorrectMatches","ProspMem",
          "FluidIQ","SymbolDigit","PairMatch",
          "TMTA","TMTB","NumericMem",
          "LondonTower","Matrices")
explorecogvars(cogvars)



# IMPUTE & PCA OF COGNITION --------------------------------------------------------
#-1 full cog imputation
#-2 no imputation (no cog NAs)
for (jj in c(1,2)) {
  cognition_only = UKBCOG[, c("eid", cogvars), ]
  
  #n missing
  sapply(cognition_only, function(x) sum(is.na(x)))
  sapply(cognition_only, function(x) sum(is.infinite(x)))
  
  #2371
  #remove N with all cogtests missing
  nocogvars = cognition_only %>% filter(
    is.na(WordPairs) &
      is.na(IncorrectMatches) &
      is.na(ProspMem) &
      is.na(FluidIQ) &
      is.na(SymbolDigit) &
      is.na(PairMatch) &
      is.na(TMTA) &
      is.na(TMTB) &
      is.na(NumericMem) &
      is.na(LondonTower) &
      is.na(Matrices
      )
  )
  
  if (jj==1) {
    print("full cog imputation")
    cognition_only %<>% filter(!eid %in% nocogvars$eid)
    nrow(cognition_only)
  
  } else if (jj == 2) {
    print("checking against nonimputed. Removing all missing cogvars")
    cognition_only = cognition_only %>% na.omit()
    # tmp2 = tmp[, 2:length(tmp)]
    # nrow(tmp2)
  }
  sublink = cognition_only$eid
  
  if (jj == 1) {
    
    #IMPUTE COGNITION #----
    print("imputing")
    imputeData = imputePCA(cognition_only[, 2:length(cognition_only)])
    imputed = imputeData$completeObs

    #check imputation
    # library("naniar"); library("VIM")
    # gg_miss_var(as.data.frame(cognition_only))
    # res <- summary(aggr(cognition_only, sortVar = TRUE))$combinations
    # (res[rev(order(res[,2])),])
    # matrixplot(cognition_only, sortby = 2)
    
    
    #tentatively optimized using general cross validation
    print("optimizing")
    nb <- estim_ncpPCA(cognition_only[, 2:length(cognition_only)],
                       method.cv = "gcv", verbose = FALSE)
    plot(0:5, nb$criterion, xlab = "nb dim", ylab = "MSEP")
    
    
    #impute again
    imputeData2 <- imputePCA(
      cognition_only[, 2:length(cognition_only)], ncp = nb$ncp) # iterativePCA algorithm
    imputed2 = imputeData2$completeObs
    
    #default
    pca1 = prcomp(imputed, center = T, scale. = T)
    summary(pca1)
    pcadf=data.frame(eid=sublink,
                     pc1=pca1$x[,1],
                     pc2=pca1$x[,2])
    #tentatively optimized
    pca2 = prcomp(imputed2, center = T, scale. = T)
    summary(pca2)
    pcadf2=data.frame(eid=sublink,
                      pc1=pca2$x[,1],
                      pc2=pca2$x[,2])
    imputed2 = as.data.frame(imputed2)
    names(imputed2)=paste0("fullimp",names(imputed2))
    imputed2$eid=sublink
    
  } else {
    
    #default PCA if nonimputed
    pca2 = prcomp(cognition_only[, 2:length(cognition_only)], center = T, scale. = T)
    summary(pca2)
    pcadf2=data.frame(eid=sublink,
                      pc1=pca2$x[,1],
                      pc2=pca2$x[,2])
  }
  
  if (jj == 1) {
    if (saveres == 1) {
      if (!dir.exists(here("results/indDiffs"))) {
        dir.create(here("results/indDiffs"))
      }
      write.table(round(pca2$rotation,3),here("results/indDiffs/cogweightsPCA.csv"),quote=F)
      xx=summary(pca2)
      write.table(round(xx$importance[2:3,],3),here("results/indDiffs/cogimportancePCA.csv"),quote=F)
    }
  }

  MERGEDAT = left_join(pcadf2, UKBCOG, by = "eid")
  
  
  
  #inverse PC to correlate negatively with age (if needed)
  if (sign(cor(MERGEDAT$eAge,MERGEDAT$pc1)) == 1) {
    print("inversing PC1")
    MERGEDAT[["pc1"]]= MERGEDAT[["pc1"]]*-1
    cor.test(MERGEDAT$pc1,MERGEDAT$eAge)
  }
  ggplot(MERGEDAT, aes(eAge,pc1)) + geom_point(alpha=0.05)
  
  
  
  # COGNITION EFFECTS -------------------------------------------------------
  for (i in 1:length(labs)) {
    if (i==1) {
      mod_cog=gc=list()
      pcog=0
      ecog=cicog=0
      CC=list()
    }
    
    c=labs[i]
    print(c)
    MERGEDAT$Y=MERGEDAT[[c]]
    
    
    #standardized scores
    mcog = lm(
      scale(MERGEDAT$Y, center=T,scale=T) ~ scale(MERGEDAT$pc1, center = T, scale = T)
      + scale(MERGEDAT$eAge, center = T, scale = F)
      + MERGEDAT$sex +
        scale(MERGEDAT$ICV, center = T, scale = T)
    )
    mod_cog[[i]]=summary(mcog)$coefficients
    
    
    #pvals
    pcog[i]=mod_cog[[i]][2,4]
    
    #effect
    ecog[i]=mod_cog[[i]][2,1]
    
    #effect
    cicog[i]=mod_cog[[i]][2,2]*1.96
    
    MERGEDAT$predcog = predict(mcog, newdata = data.frame(
      Y = MERGEDAT$Y,
      pc1 = MERGEDAT$pc1,
      eAge = 0, #mean age
      sex = 0, #mean sex
      ICV = 0)) #mean ICV
    MERGEDAT$partial_residcog = MERGEDAT$predcog + residuals(mcog)
    
    #save confidence intervals
    CC[[i]] = confint(mcog)[2,]
    
    #plot
    gc[[i]]=ggplot(MERGEDAT,aes(x=pc1,y=partial_residcog,fill=pc1)) +
      geom_point() +
      geom_smooth(method="lm",col="blue")
  }
  
  if (jj==1) {
  res_cog = data.frame(
    p = pcog,
    effect = ecog,
    ci = cicog,
    N = seq(labs),
    labs,
    meas = "cognition")
  } else {
    res_cog_nonimpute = data.frame(
      p = pcog,
      effect = ecog,
      ci = cicog,
      N = seq(labs),
      labs,
      meas = "cognition")
  }
  
  if (jj==1) {
    #save significant cognition effect
    cog_survive = res_cog[-log10(res_cog$p) >= -log10(0.01/bonf),]
    CHECKED = data.frame(res_cog[8 ,], CI_lwrupr = unlist(CC[[8]])) %>% mutate(imputedtype = "non",
                                                                               Nobs = nrow(MERGEDAT))
    
    # ASSOCIATION PLOT (cognition) -------------------------------------------
    #isolate only surviving
    tmp = res_cog[-log10(res_cog$p) >= -log10(0.01/bonf),]
    
    p_cog = ggplot(res_cog) +
      geom_point(data=tmp,aes(x=N,y=-log10(p)),col="#81A88D",size=2) +
      geom_bar(aes(x=N,y=-log10(p)),col="#81A88D",stat="identity", width = 0.01) +
      geom_hline(yintercept = -log10(0.01)) +
      geom_hline(yintercept = -log10(0.01/bonf),col="black",linetype=2) +
      coord_cartesian(ylim=c(0,max(-log10(res1$p)))) +
      theme_classic() +
      mytheme
    
    
    if (saveres == 1) {
      filename = here(paste0("results/indDiffs/p_cog_",brainvar,".png"))
      print(paste("saving",filename))
      
      if (!dir.exists(here("results/indDiffs"))) {
        dir.create(here("results/indDiffs"))
      }
      ggsave(filename = filename, plot = p_cog, width = 3, height=6.5, dpi=600, units="cm")
    }
    
    
    
    #plot surviving cognition effect on areal plot
    if (brainvar=="area") {
      addboth = rbind(eff_survive,
                      cog_survive %>% mutate(dir="leftward",
                                       n=N) %>% select(names(eff_survive))) %>% arrange(N)
      addboth$N[1]=4 #positioning
      
      p_effects = ggplot(addboth) +
        geom_point(aes(x=N,y=effect,col=factor(meas)),size=2) +
        geom_bar(aes(x=N,y=effect,col=factor(meas)),stat="identity", width = 0.01,alpha=0.05) +
        geom_errorbar(aes(x=N,y=effect,ymin=effect-ci,ymax=effect+ci,col=factor(meas)),width=0.75) +
        geom_hline(yintercept = 0,col="black",linetype=2) +
        scale_color_manual(values=c("#81A88D",
                                    wesanderson::wes_palette(n=5, name="Cavalcanti1")[c(1:4,5)])) +
        scale_color_manual(values=c("#81A88D",
                                    wesanderson::wes_palette(n=5, name="Cavalcanti1")[c(1:4,5)])) +
        coord_cartesian(xlim = c(0,max(eff_survive$N))) +
        ylab("Standardized β") +
        theme_classic() +
        mytheme +
        theme(legend.position = "none")  
      
      #append plots
      p_both = cowplot::plot_grid(p_pval,p_effects,ncol=2)
      
      if (saveres == 1) {
        filename = here(paste0("results/indDiffs/p_effects_",brainvar,".png"))
        print(paste("saving",filename))
        ggsave(filename = filename, plot = p_both, width = 18, height=6.5, dpi=600, units="cm")
      }
    }

    
    #COMPILE OUTPUT TABLE
    TABLE2 = res_cog %>% select(labs,effect,p)
    names(TABLE2)[2]="cognition.effect"
    names(TABLE2)[3] = "p.cog"
    TABLE2$cognition.effect=round(TABLE2$cognition.effect,3)
    
    if (brainvar=="thickness"){
      TABLE2$roi=substr(TABLE2$labs,24,25) 
      TABLE2 = rbind(TABLE2[10:nrow(TABLE2),],
                     TABLE2[1:9,])
    } else {
      TABLE2$roi=substr(TABLE2$labs,19,20)
      TABLE2 = rbind(TABLE2[8:nrow(TABLE2),],
                     TABLE2[1:7,])
    }
    TABLE2$p.cog=formatC(TABLE2$p,format="e",digits = 2)
    TABLE2$roiname = TABLE1$roiname
    TABLE2 %<>% select(labs,roi,roiname,everything()) 
      
    
    TABLEOUT = cbind(TABLE1,TABLE2[,4:5])
    TABLEOUT %<>% select(1:3,cognition.effect,p.cog,everything())
    

  } else {

  # NON IMPUTED EFFECT ------------------------------------------------------
    CHECKED = rbind(CHECKED,
                             data.frame(res_cog[8 ,], CI_lwrupr = unlist(CC[[8]])) %>% mutate(imputedtype = "non",
                                                                                              Nobs = nrow(MERGEDAT)))
  }
}
CHECKED
res_cog[which(res_cog==min(res_cog$p)),]
res_cog_nonimpute[which(res_cog_nonimpute==min(res_cog_nonimpute$p)),] #check lowest pval with no imputation

if (saveres == 1) {
  write.table(TABLEOUT,here(paste0("results/indDiffs/effects_",brainvar,"_fullimputation.tsv")),quote=F,row.names=F,col.names = T)
}


