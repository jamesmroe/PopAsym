# install.packages("here")
library("here")


#SELECT OPTIONS ======================
args = commandArgs(TRUE)
datdir=as.character(args[1])
saveplots=1
brainvar="thickness"
# brainvar="area"
#====================================#


#---load packages
tmp.packages = c("tidyverse","readr","devtools","itsadug","magrittr","ggridges","gamm4","here")
tmpnew.packages = tmp.packages[!(tmp.packages %in% installed.packages()[,"Package"])]
if(length(tmpnew.packages)) {
  install.packages(tmpnew.packages)
}
sapply(tmp.packages, require, character.only = T)
rm(list=ls(pattern="tmp*"))


#---data
D = read.csv(file.path(datdir,paste0("PopAsym_LCBC_", brainvar, ".csv")), stringsAsFactors = F, header = T, sep =",")


#---clusters
mergelabs = names(D %>% select(contains("merge")))
names(D)[(names(D)) %in% mergelabs] = toupper(mergelabs)
labs = names(D %>% select(contains(".label")))


#---COMPUTE TRAJECTORIES-------------
pb = txtProgressBar(min=1, max=length(labs), style=3)
span="Life"
for (s in span){
  DF = D

  for (i in 1:length(labs)) {
    if (i==1) {
      kk=6 #knots
      pp=list() #plots
      RR=gammtable=list() #results
      omegasq=0
      outperlab=0 #track outliers per cluster
      hemi_trajectories_knotlevels = list()
      diff_trajectories_knotlevels = list()
      clusters = list()
    }
    setTxtProgressBar(pb,i)
    
    whichlabel=labs[i]
    DF$brainvar = DF[[whichlabel]]
    
    LH = DF %>% filter(hemi==1)
    RH = DF %>% filter(hemi==0)
    
    #====================== GAMM ====================== 
    g = gamm4(brainvar ~ s(Age, by = as.factor(hemi), k = kk) + as.factor(hemi) + Sex + Site_Name, data = DF, random = ~ (1 | fsid_base))
    g.sum = summary(g$gam)
    
    #--- 1) predictions
    predR = itsadug::get_predictions(g$gam,
                                       cond = list(hemi = 0, 
                                                   Sex = "female",
                                                   Site_Name = "ousAvanto",
                                                   Age = DF$Age[DF$hemi == 0], se=T))
    predL = itsadug::get_predictions(g$gam,
                                       cond = list(hemi = 1, 
                                                   Sex = "female",
                                                   Site_Name = "ousAvanto",
                                                   Age = DF$Age[DF$hemi == 1], se=T))
    residualsL <- residuals(g$mer)[DF$hemi == 1]
    residualsR <- residuals(g$mer)[DF$hemi == 0]
    
    
    #quick plots
    l=nrow(DF)/2
    pdat = data.frame("Age" = seq(from = min(DF$Age),to = max(DF$Age),length.out = l),
                      "hemi" = rep(1,l),
                      "Sex" = rep("female",l),
                      "Site_Name" = rep("ousAvanto",l))
    pdat = rbind(pdat,pdat)
    pdat$hemi[l:nrow(DF)] = 0
    Xp = predict(g$gam, newdata = pdat, type = "lpmatrix")
    plot(pdat$Age,Xp %*% coef(g$gam))
    
  
    #add predictions to data
    LH <- LH %>% 
      mutate(
        predictionsL = predL$fit,
        partial_residualsL = predL$fit + residualsL
      )
    RH <- RH %>% 
      mutate(
        predictionsR = predR$fit,
        partial_residualsR = predR$fit + residualsR
      )
    
    
    #set colours based on variable code
    tmpsplit = substr( paste0(strsplit(whichlabel,"\\.")[[1]][2:4]) [2],2,2)
    if (tmpsplit == "n" || tmpsplit == "N")  {
      col = c("palevioletred4","palevioletred1","palevioletred") 
    } else  {
      col = c("lightgoldenrod4","lightgoldenrod1","lightgoldenrod3") 
    } 
    
    
    #--- 2) set threshold for outliers
    Opt.age = predL$Age
    outlierthresh = 6
    
    
    #--- 3) detect and remove outliers
    N.outliers = c()
    QL = scale(g$gam$residuals[DF$hemi == 1])
    QR = scale(g$gam$residuals[DF$hemi == 0])
    db.out = data.frame(LH$Folder,QL,QR)
    db.out$QQL = NA; db.out$QQR = NA
    db.out$QQL [(abs(QL)>outlierthresh)] = QL[(abs(QL)>outlierthresh)] #place outliers in new col
    db.out$QQR [(abs(QR)>outlierthresh)] = QR[(abs(QR)>outlierthresh)]
    outliersL = which(!is.na(db.out$QQL))
    outliersR = which(!is.na(db.out$QQR))
    outliers = c(outliersL, outliersR[!(outliersR %in% outliersL)]) #outliersL + unique outliersR 
    cat("\noutlier thresh =", outlierthresh, "SD from GAMM model: removing", length(outliers), "cases\n\n")
    outperlab[i] = length(outliers)
    
    if (outperlab[i]>0) {

      outobs = LH[["Folder"]][outliers]
      LHc = LH %>% filter(!(Folder %in% outobs)) %>% select(-predictionsL, -partial_residualsL)
      RHc = RH %>% filter(!(Folder %in% outobs)) %>% select(-predictionsR, -partial_residualsR)
      DFc = rbind(LHc, RHc)
      
    }  else {
      
      LHc=LH
      RHc=RH
      DFc=DF
    
    }
    
    #====================== GAMM REDO ====================== 
    #--- 4) redo gamm analysis
    gamm.trajectories = gamm4(brainvar ~ s(Age, by = as.factor(hemi), k = kk) + as.factor(hemi) + Sex + Site_Name, data = DFc, random = ~ (1 |fsid_base))
    gamm.sum = summary(gamm.trajectories$gam)
    g = gamm.trajectories$gam

    #--- 5) redo predictions
    predR = itsadug::get_predictions(gamm.trajectories$gam,
                                     cond = list(hemi = 0, 
                                                 Sex = "female",
                                                 Site_Name = "ousAvanto",
                                                 Age = DFc$Age[DFc$hemi == 0], se=T))
    predL = itsadug::get_predictions(gamm.trajectories$gam,
                                     cond = list(hemi = 1, 
                                                 Sex = "female",
                                                 Site_Name = "ousAvanto",
                                                 Age = DFc$Age[DFc$hemi == 1], se=T))
    residualsL <- residuals(gamm.trajectories$mer)[DFc$hemi == 1]
    residualsR <- residuals(gamm.trajectories$mer)[DFc$hemi == 0]
    
    #add predictions to data
    LHc %<>% 
      mutate(
        predictionsL = predL$fit,
        partial_residualsL = predL$fit + residualsL
      )
    RHc %<>% 
      mutate(
        predictionsR = predR$fit,
        partial_residualsR = predR$fit + residualsR
      )
    
    
    #LIFESPAN PLOTS
    mytheme = theme(
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      text = element_text(color = "black", size = 18, family = "Helvetica Neue Light"),
      plot.title = element_text(hjust = 0.5),
      axis.ticks = element_blank(),
      axis.title.y = element_text(color = "black", size = 22, vjust =-1, margin = margin(0,20,0,0)),
      axis.title.x = element_text(color = "black", size = 22, vjust = -2, margin = margin(0,20,0,0)),
      axis.text = element_text(color = "black", size = 22),
      legend.text = element_text(color = "white", size = 15),
      legend.key.size = unit(1,"cm"))
    
    pp[[i]] =
      ggplot() +
      geom_line(data=LHc,aes(x=Age,partial_residualsL,group=fsid_base),color="grey",alpha=0) +  #keep scaling by making transparent
      geom_point(data=LHc,aes(x=Age,partial_residualsL,group=fsid_base),color=col[1],stat="identity",alpha=0) +  #keep scaling by making transparent
      geom_line(data=RHc,aes(x=Age,partial_residualsR,group=fsid_base),color="grey",alpha=0) +  #keep scaling by making transparent
      geom_point(data=RHc,aes(x=Age,partial_residualsR,group=fsid_base),color=col[3],stat="identity",alpha=0) + #keep scaling by making transparent
      geom_ribbon(data=predL,aes(x=Age,ymin=fit-CI,ymax=fit+CI),alpha=.4,show.legend=F,fill=col[1]) +
      geom_line(data=LHc,aes(x=Age,predictionsL),color=col[1],size=1.5) +
      geom_ribbon(data=predR,aes(x=Age,ymin=fit-CI,ymax=fit+CI),alpha=.4,show.legend=F,fill=col[2]) +
      geom_line(data=RHc,aes(x=Age,predictionsR),color=col[2],size=1.5) +
      geom_line(data=predR,aes(x=Age,fit),color=col[2],size=1.5) +
      geom_line(data=predL,aes(x=Age,y=fit-CI),alpha=1,show.legend=F,color="black",size=0.25) +
      geom_line(data=predL,aes(x=Age,y=fit+CI),alpha=1,show.legend=F,color="black",size=0.25) +
      geom_line(data=predR,aes(x=Age,y=fit-CI),alpha=1,show.legend=F,color="black",size=0.25) +
      geom_line(data=predR,aes(x=Age,y=fit+CI),alpha=1,show.legend=F,color="black",size=0.25) +
      ggtitle(whichlabel) +
      scale_x_continuous(c(0,20,40,60,80)) +
      labs(x = "Age", y = paste(brainvar, "(residualised)")) +
      mytheme
    
    #add border to merge plots
    if (tmpsplit == "N" || tmpsplit == "P")  {
      pp[[i]] = pp[[i]] + theme(
        panel.border = element_rect(colour = "grey", fill=NA, size=4))
    }
    
    #save plot
    if (saveplots == 1) {
      if (! dir.exists(here("results/lifespan"))) {
        dir.create(here("results/lifespan"))
      }
      filename = file.path(here("results/lifespan"), paste0(s,"GammCI_", whichlabel,kk,"-publish.jpg"))
      ggsave(filename = filename, plot = pp[[i]], width = 9, height=13, dpi=600, units="cm")
    }
    
    
    #====================== ORDERED FACTOR GAMM ====================== 
    #order factor to get difference trajectories
    DFc %<>% mutate(ohemi=ifelse(hemi==1,"left","right")) %>% mutate(ohemi=as.ordered(ohemi))
    
    ogc = gamm4(brainvar ~ as.factor(hemi) + s(Age,k = kk) + s(Age, by = ohemi, k = kk) + Sex + Site_Name, data = DFc, random = ~ (1 | fsid_base))
    # plot.gam(ogc$gam)
    osum = summary(ogc$gam)
    hemieffect = osum$p.coeff[2]
    
    
    #make 100 obs age model
    nn=100
    Optage=seq(min(DFc$Age), max(DFc$Age), length.out = nn) 
    fake.frameX = data.frame("Age" = Optage,
                             "hemi" =rep(1,nn),
                             "Sex" = rep("female",nn),
                             "Site_Name" = rep("ousAvanto",nn))
    pdat = rbind(fake.frameX,fake.frameX)
    pdat$hemi[1:100] = 0 #prediction data
    
    
    #simulate from posterior - Xp matrix
    #muliply matrix by model coefficients and sum rows to get predicted values
    Xp = predict(g, newdata = pdat, type = "lpmatrix") 
    # plot(Xp %*% coef(g))
    
    #which cols of Xp relate to which hemi
    c1 = grepl("hemi\\)1", colnames(Xp)) #L
    c2 = grepl("hemi\\)0", colnames(Xp)) #R
    
    #which rows of Xp
    r1 = with(pdat, hemi == 1) #L
    r2 = with(pdat, hemi == 0) #R
    
    
    #differences between smooths
    #cols of differenced Xp matrix that aren't involved in comparison set to zero
    X = Xp[r1, ] - Xp[r2, ] #L-R
    X[, !(c1 | c2)] = 0
    X[, !grepl("s\\(", colnames(Xp))] = 0
    
    
    #get predicted vals of zero-centered hemispheric difference at zero of other covs
    #SE of the difference using variance-covariance matrix of estimated model coefficients
    dif = X %*% coef(g)
    se = sqrt(rowSums((X %*% vcov(g, unconditional =T)) * X))
    
    
    #critical value of the t distribution for confint
    crit95 = qt(.975, df.residual(g))
    upr = dif + (crit95 * se)
    lwr = dif - (crit95* se)
    uprse=dif+se
    lwrse=dif-se
    
    comp = data.frame("OptAge" = Optage, 
                      dif, se, lwr, upr)
    hemieffectcorr = gamm.sum$p.coeff[2]
    
    
    direct=substr(strsplit(whichlabel,"\\.")[[1]][3],2,2)
    if (direct == "n"){
      if (brainvar=="area") {
        rescale = c(-0.16,0) #plot in same space
      } else if (brainvar == "thickness") {
        rescale = c(-0.25,0.05) #plot in same space
      }
    } else if (direct == "p") {
      if (brainvar=="area") {
        rescale = c(0,0.16) #plot in same space
      } else if (brainvar == "thickness") {
        rescale = c(-0.05,0.25) #plot in same space
      }
    }  
    
    #difference
    dpp = ggplot(comp) +
      geom_line(aes(Optage,lwr+hemieffectcorr),col=col[1],size=1,alpha=0.7) +
      geom_line(aes(Optage,upr+hemieffectcorr),col=col[1],size=1,alpha=0.7) +
      geom_ribbon(aes(Optage,ymin=lwr+hemieffectcorr,ymax = upr+hemieffectcorr), alpha = 0.3,fill=col[1]) +
      geom_hline(aes(yintercept=0),linetype=5,col="black") +
      labs(x="Age",y="LH-RH") +
      theme_classic() +
      coord_cartesian(ylim=rescale) +
      theme(text=element_text(size=18),
            plot.title = element_text(hjust = 0.5),
            axis.ticks = element_blank(),
            axis.title.y = element_text(color = "black", size = 22, vjust =-1, margin = margin(0,20,0,0)),
            axis.title.x = element_text(color = "black", size = 22, vjust = -2, margin = margin(0,20,0,0)))
    
    #save plot
    if (saveplots == 1) {
      filename = here("results/lifespan", paste0(s,"DiffgammCI_", whichlabel,kk,"-publish.jpg"))
      ggsave(filename = filename, plot = dpp, width = 7, height=5, dpi=600, units="cm") #width 6
    }
    
    
    #gammstats
    RR$label[i] = whichlabel
    RR$edf[i] = osum$edf[2]
    RR$Fval[i] = osum$s.table[[2,3]]
    RR$p.spl[i] = -log10(osum$s.pv[2])
    RR$pp.spl[i] = osum$s.pv[2]
    RR$rsqad[i] = osum$r.sq
    RR$nobs[i] = dim(DFc)[1]
    RR$hemieffect[i] = osum$p.table[2,1]
    RR$hemip[i] = osum$p.table[2,4]
    RR$sexeffect[i] = osum$p.table[3,1]
    RR$sexp[i] = osum$p.table[3,4]
    RR$dif[[i]] = dif
    # RR$scaneffect[[c]] = ogamm.sum$p.table[4,4]
    
    RR$fsq[i] =( RR$edf[i]*(RR$Fval[i]-1) ) / RR$nobs[i]
    omegasq[i] = RR$fsq[i] / (1+RR$fsq[i])
    gammtable[[i]]=data.frame("cluster" = whichlabel,
                              outperlab[i],
                              "edf" = round(RR$edf[i],1),
                              "modelrsqad" = RR$rsqad[i],
                              "F" = round(RR$Fval[i],1),
                              "p" = formatC(RR$pp.spl[i],format="e",digits=2),
                              "omegasq"=omegasq[i],
                              "hemif" = round(RR$hemieffect[i],2),
                              "hemip"=  formatC(RR$hemip[i],format="e",digits=2),
                              "sexeffect" = round(RR$sexeffect[i],2),
                              "p.sex" = formatC(RR$sexeffect[i],format="e",digits=2),
                              nobs=RR$nobs[i])
  }
}
o = bind_rows(gammtable)
o$modelrsqad=round(o$modelrsqad,3)
o$omegasq=round(o$omegasq,3)
if (brainvar=="area") {
  o$roiname=c("R_Parieto-occipital_sulcus",
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
              "L_Superior-frontal-gyrus",
              "M1","M2")
  
  #reorder table to L R MERGE
  o=rbind(o[8:14,],
          o[1:7,],
          o[16:15,])
} else {
  o$roiname=c("R_Superior-temporal-sulcus",
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
              "L_Supplementary-motor-cortex",
              "L_Collateral-sulcus",
              "L_Anterior-transverse-collateral-sulcus",
              "L_Caudal-middle-frontal",
              "L_Caudal-superior-frontal",
              "L_Rostral-superior-frontal",
              "L_Calcarine-sulcus",
              "M1","M2")
  
  #reorder table to L R MERGE
  o=rbind(o[10:(nrow(o)-2),],
          o[1:9,],
          o[c(22, 21),])
}
o %<>% select(cluster,roiname,edf,"F",p,omegasq,hemif,hemip,sexeffect,p.sex,modelrsqad,nobs,everything()) %>% 
  mutate(Outliers_excluded = outperlab.i.) %>% select(-outperlab.i.)

if (saveplots == 1) {
  write.table(o, here("results/lifespan",paste("gammtable",brainvar,"csv",sep=".")), row.names = F, col.names = T, sep = " ", quote = F)
}
