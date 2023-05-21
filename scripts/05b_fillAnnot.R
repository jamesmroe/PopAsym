rm(list=ls())
library("tidyverse")
library("magrittr")
library("here")
here()

#fill in annotation with heritability values

#SELECT OPTIONS ======================
#SCRIPT TAKES SUMMARY-LEVEL DATA AS INPUT AND IS THUS FULLY REPRODUCIBLE WITH THE DATA PROVIDED HERE
saveres=0 #1=yes, 0=no
typeHerit = "snp" #snp/twin
#====================================#



#FILL RESULTS OF SNP-BASED ANALYSIS
if (typeHerit == "snp") {
  
  rdir=here("results/heritability")
  sdir=here("annotInfo")
  
  RESAREA = read.csv(here("results/heritability/cortexH2SNParea.csv"), sep=" ")
  RESCTH = read.csv(here("results/heritability/cortexH2SNPthickness.csv"), sep=" ")
  RES = rbind(RESAREA, RESCTH)
  
  
  #---FDR
  all.equal(RES$FDR[1:500],p.adjust(RES$sig[1:500], "BH"))
  all.equal(RES$FDR[501:1000],p.adjust(RES$sig[501:1000], "BH"))
  all.equal(RES$bin, ifelse(RES$FDR<.05, TRUE, FALSE))
  sgof::BH(RESAREA$sig,alpha = .05)
  sum(RES$FDR[1:500]<.05)
  sgof::BH(RESCTH$sig,alpha = .05)
  sum(RES$FDR[501:1000]<.05)
  
  
  #cortex-wide SNP heritability for thickness was signifcantly lower
  summary(lm((scale(h2)) ~ type, data = RES))
  ggplot(RES, aes(x=type, y=h2)) + geom_jitter()
  t.test(RES$h2[RES$type=="thickness"], RES$h2[RES$type=="area"])
  
  
  #percentage of significant SNP h2 after FDR per type
  sum(RES$bin[RES$type=="area"])/length(RES$type[RES$type=="area"]) #53% area
  sum(RES$bin[RES$type=="thickness"])/length(RES$type[RES$type=="thickness"]) #11% cth
  
  
  #open vtx annot file and replace codes with vals
  V=read.csv(here("annotInfo/lh.vtxs.lh.Schaefer2018_1000Parcels_17Networks_order.40tvs8ref_sym_20.txt"),stringsAsFactors = F) #vertex code file
  A=read.csv(here("annotInfo/lh.annot.lh.Schaefer2018_1000Parcels_17Networks_order.40tvs8ref_sym_20.txt"),stringsAsFactors = F) #annot code file
  D=RES$name
  
  
  #UKBSNP results labels FDR
  res1snp=RESAREA$name[RESAREA$FDR<.05]
  RESAREA$FDR[RESAREA$FDR<.05]
  RESAREA$name[RESAREA$FDR==min(RESAREA$FDR)]
  RESAREA[RESAREA$h2==max(RESAREA$h2),]
  RESCTH[RESCTH$h2==max(RESCTH$h2),]
  
  res2snp=RESCTH$name[RESCTH$FDR<.05]
  RESCTH$FDR[RESCTH$FDR<.05]
  RESCTH$name[RESCTH$FDR==min(RESCTH$FDR)]
  
  
  C=A$code
  D=lapply(D, function(x) strsplit(x,".mgh")[[1]][1]) %>% unlist()
  S=paste0("lh.",A$struct_names)
  VVA=SSA=VVT=SST=V
  for (i in 1:length(C)) {
    print(i)
    co=C[i] #annot code
    struct=S[i]
    if (i==1) { #if medial wall
      VVA[which(VVA==co),]=0
      VVT[which(VVT==co),]=0
      SSA[which(SSA==co),]=0
      SST[which(SST==co),]=0
    } else {
      VVA[which(VVA==co),]= RESAREA$h2[which(RESAREA$name==struct)]
      SSA[which(SSA==co),]= -log10(RESAREA$FDR[which(RESAREA$name==struct)])
      
      VVT[which(VVT==co),]= RESCTH$h2[which(RESCTH$name==struct)]
      SST[which(SST==co),]= -log10(RESCTH$FDR[which(RESCTH$name==struct)])
    }
  } 
  
  if (saveres == 1) {
    print("writing maps to csv")
    write.table(VVA,here("results/heritability/mapH2UKBSNP.area.csv"),row.names=F,col.names = F,quote=F)
    write.table(SSA,here("results/heritability/mapSigFDRUKBSNP.area.csv"),row.names=F,col.names = F,quote=F)
  
    write.table(VVT,here("results/heritability/mapH2UKBSNP.thickness.csv"),row.names=F,col.names = F,quote=F)
    write.table(SST,here("results/heritability/mapSigFDRUKBSNP.thickness.csv"),row.names=F,col.names = F,quote=F)
    
    print("writing FDR-corrected significant label names to file")
    writeLines(RESAREA$name[RESAREA$FDR<.05],file.path(rdir,"resAreaUKBSNPSigFDR.txt"))
    writeLines(RESCTH$name[RESCTH$FDR<.05],file.path(rdir,"resThicknessUKBSNPSigFDR.txt"))
  }
  
  
  
  
  if (saveres == 1) {
    #bash code to find labels easily (copy to tsd)
    # sdir="/cluster/projects/p23/projects/11_MemP/James/UKBSNP/shared/annotInfo"
    # rdir="/cluster/projects/p23/projects/11_MemP/James/UKBSNP/results.SNPHerit/results.schaf"
    # for i in $(cat resAreaUKBSNPSigFDR.txt); do cp $sdir/$i.mgh.label $rdir/resAreaUKBSNPSigFDR.$i.mgh.label; done
    # for i in $(cat resThicknessUKBSNPSigFDR.txt); do cp $sdir/$i.mgh.label $rdir/resThicknessUKBSNPSigFDR.$i.mgh.label; done
    
    #copy FDR-corrected significant labels
    if (saveres == 1) {
      if (!dir.exists(file.path(rdir,"labels"))) { dir.create(file.path(rdir,"labels")) }
      print("running system command to copy signficant labels to results/labels")
      cmd = paste(paste0('sdir=',sdir),
                  paste0('rdir=',rdir),
                  paste0('for i in $(cat $rdir/resAreaUKBSNPSigFDR.txt); do cp $sdir/labels/$i.mgh.label $rdir/labels/resAreaUKBSNPSigFDR.$i.mgh.label; done'),
                  paste0('for i in $(cat $rdir/resThicknessUKBSNPSigFDR.txt); do cp $sdir/labels/$i.mgh.label $rdir/labels/resThicknessUKBSNPSigFDR.$i.mgh.label; done'),
                  sep='\n')
      system(cmd)
    }
  }
}


#FILL RESULTS OF TWIN ANALYSIS
if (typeHerit == "twin") {

  # b="/cluster/projects/p23/projects/11_MemP/James/UKBSNP"
  # rdir=file.path(b,"results.SNPHerit/results.schaf")
  # sdir=file.path(b,"shared/annotInfo")

  rdir=here("results/heritability")
  sdir=here("annotInfo")
  
  h2area=read.csv(file.path(rdir,"cortexTwinAEarea_revision.csv"),stringsAsFactors = F,sep=" ",header=T)
  h2cth=read.csv(file.path(rdir,"cortexTwinAEthickness_revision.csv"),stringsAsFactors = F,sep=" ",header=T)
  h2area$type="area"
  h2cth$type="thickness"

  
  #open vtx annot file and replace codes with vals
  V=read.csv(file.path(sdir,"lh.vtxs.lh.Schaefer2018_1000Parcels_17Networks_order.40tvs8ref_sym_20.txt"),stringsAsFactors = F)
  A=read.csv(file.path(sdir,"lh.annot.lh.Schaefer2018_1000Parcels_17Networks_order.40tvs8ref_sym_20.txt"),stringsAsFactors = F)
  identical(h2area$roi,h2cth$roi)
  D=h2area$roi

  RES = rbind(h2area,
              h2cth)

  C=A$code
  D=lapply(D, function(x) strsplit(x,".mgh")[[1]][1]) %>% unlist()
  D=lapply(D, function(x) strsplit(x,"AI.")[[1]][2]) %>% unlist()
  RES$name=rep(D,2)
  S=paste0("lh.",A$struct_names)
  VVA=SSA=VVT=SST=V
  
  
  RESAREA=RES[RES$type=="area",]
  RESCTH=RES[RES$type=="thickness",]
  
  #original ACE model
  # RESAREA$FDR=p.adjust(RESAREA$pA,method="BH")
  # sgof::BH(RESAREA$pA,alpha = .05)
  RESAREA$FDR=p.adjust(RESAREA$pArevise,method="BH")
  sgof::BH(RESAREA$pArevise,alpha = .05)
  
  #revised AE model  
  # RESCTH$FDR=p.adjust(RESCTH$pA,method="BH")
  # sgof::BH(RESCTH$pA,alpha = .05)
  RESCTH$FDR=p.adjust(RESCTH$pArevise,method="BH")
  sgof::BH(RESCTH$pArevise,alpha = .05)
  
  #order by h2
  tmp=RESAREA[
    with(RESAREA, order(a2revise)),
  ]
  
  RESAREA[RESAREA$a2revise==max(RESAREA$a2revise),]
  RESCTH[RESCTH$a2revise==max(RESCTH$a2revise),]
  
  #twin results labels FDR
  res1twin=RESAREA$name[RESAREA$FDR<.05]
  RESAREA$pArevise[RESAREA$pArevise<.05]
  RESAREA$name[RESAREA$pArevise<.05]
  
  res2twin=RESCTH$name[RESCTH$FDR<.05]
  RESCTH$pArevise[RESCTH$pArevise<.05]
  RESCTH$name[RESCTH$pArevise<.05]
  
  RESAREA[RESAREA$FDR<.05,] %>% dim()
  RESCTH[RESCTH$FDR<.05,] %>% dim()
  
  
  
  for (i in 1:length(C)) {
    print(i)
    co=C[i]
    struct=S[i]
    if (i==1) {
      #medial wall
      VVA[which(VVA==co),]=0
      VVT[which(VVT==co),]=0
      SSA[which(SSA==co),]=0
      SST[which(SST==co),]=0
    } else {
      VVA[which(VVA==co),]= RESAREA$a2revise[which(RESAREA$name==struct)]
      SSA[which(SSA==co),]= -log10(RESAREA$FDR[which(RESAREA$name==struct)])
      
      VVT[which(VVT==co),]= RESCTH$a2revise[which(RESCTH$name==struct)]
      SST[which(SST==co),]= -log10(RESCTH$FDR[which(RESCTH$name==struct)])
    }
  }  
  
  if (saveres == 1) {
    print("writing maps to csv")
    write.table(VVA,file.path(rdir,"mapH2TwinAEarea_revision.csv"),row.names=F,col.names = F,quote=F)
    write.table(SSA,file.path(rdir,"mapSigFDRTwinAEarea_revision.csv"),row.names=F,col.names = F,quote=F)
    
    write.table(VVT,file.path(rdir,"mapH2TwinAEthickness_revision.csv"),row.names=F,col.names = F,quote=F)
    write.table(SST,file.path(rdir,"mapSigFDRTwinAEthickness_revision.csv"),row.names=F,col.names = F,quote=F)
    
    print("writing FDR-corrected significant label names to file")
    writeLines(RESAREA$name[RESAREA$FDR<.05],file.path(rdir,"resAreaTwinSigFDR_revision.txt"))
    writeLines(RESCTH$name[RESCTH$FDR<.05],file.path(rdir,"resThicknessTwinSigFDR_revision.txt"))
    
    #p > 05 for comparing twin results
    writeLines(RESAREA$name[RESAREA$pArevise<.05],file.path(rdir,"resAreaTwinSig05_revision.txt"))
    writeLines(RESCTH$name[RESCTH$pArevise<.05],file.path(rdir,"resThicknessTwinSig05_revision.txt"))
  }
  
  
  
  if (saveres == 1) {
    #bash code to find labels easily (copy to tsd)
    # sdir="/cluster/projects/p23/projects/11_MemP/James/UKBSNP/shared/annotInfo"
    # rdir="/cluster/projects/p23/projects/11_MemP/James/UKBSNP/results.SNPHerit/results.schaf"
    # for i in $(cat $sdir/labels/resAreaTwinSig_revison.txt); do cp $sdir/labels/$i.mgh.label $rdir/FINALresAreaTwinSig05.$i.mgh.label; done
    # for i in $(cat resThicknessTwinSig_revision.txt); do cp $sdir/$i.mgh.label $rdir/FINALresThicknessTwinSig05.$i.mgh.label; done
    
    #copy FDR-corrected significant labels
    if (saveres == 1) {
      if (!dir.exists(file.path(rdir,"labels"))) { dir.create(file.path(rdir,"labels")) }
      print("running system command to copy signficant labels to results/labels")
      cmd = paste(paste0('sdir=',sdir),
                  paste0('rdir=',rdir),
                  paste0('for i in $(cat $rdir/resAreaTwinSigFDR_revision.txt); do cp $sdir/labels/$i.mgh.label $rdir/labels/resAreaTwinSigFDR.$i.mgh.label; done'),
                  paste0('for i in $(cat $rdir/resThicknessTwinSigFDR_revision.txt); do cp $sdir/labels/$i.mgh.label $rdir/labels/resThicknessTwinSigFDR.$i.mgh.label; done'),
                  sep='\n')
      system(cmd)
      
      #sig 05 for checking
      cmd = paste(paste0('sdir=',sdir),
                  paste0('rdir=',rdir),
                  paste0('for i in $(cat $rdir/resAreaTwinSig05_revision.txt); do cp $sdir/labels/$i.mgh.label $rdir/labels/resAreaTwinSig05.$i.mgh.label; done'),
                  paste0('for i in $(cat $rdir/resThicknessTwinSig05_revision.txt); do cp $sdir/labels/$i.mgh.label $rdir/labels/resThicknessTwinSig05.$i.mgh.label; done'),
                  sep='\n')
      system(cmd)
    }
  }
  
  #spatial correlations between twin and SNP results
  #remove medial wall
  spatcor1area = VVA [which(V!=65793),]
  spatcor1thick = VVT [which(V!=65793),]
  
  #load in UKB H2
  spatcor2area = read.csv(here("results/heritability/mapH2UKBSNP.area.csv"), header=F) 
  spatcor2thick = read.csv(here("results/heritability/mapH2UKBSNP.thickness.csv"), header=F)
  #remove medial wall
  spatcor2area = spatcor2area[which(V!=65793),]
  spatcor2thick = spatcor2thick[which(V!=65793),]
  
  #area
  # plot(spatcor1area,spatcor2area)
  cor.test(spatcor1area,spatcor2area)
  
  
  #thick
  # plot(spatcor1thick,spatcor2thick)
  cor.test(spatcor1thick,spatcor2thick)
  
  
  #cortex-wide twin heritability for thickness was signifcantly lower
  TMP = rbind(RESAREA,RESCTH)
  ggplot(TMP, aes(x=type, y=a2revise)) + geom_jitter()
  summary(lm((scale(a2revise)) ~ type, data = TMP))
  t.test(TMP$a2revise[TMP$type=="thickness"], TMP$a2revise[TMP$type=="area"])
  
  
}

#how many sig in HCP were also sig in UKB?
sum(res1twin %in% res1snp)/length(res1twin)
sum(res2twin %in% res2snp)/length(res2twin)



#difference between area and thickness for cluster-wise analysis - SNP-based results
A=read.csv("/Users/jamesroe/Dropbox/PHD_project/Paper3/Manuscript/Tables/results.SNPHerit/H2.clust.txt",header=F, sep="\t")
A$type="area"
A$type[15:nrow(A)]="thickness"
summary(lm((scale(V2)) ~ type, data = A))


#difference between area and thickness for cluster-wise analysis - twin-based results
TMP1 = read.csv(here("results/heritability",paste0("clustersTwinAEarea_revision.csv")), sep=" ")
TMP2 = read.csv(here("results/heritability",paste0("clustersTwinAEthickness_revision.csv")), sep=" ")
TMP1$a2revise = lapply(TMP1$a2revise, function(x) strsplit(x,"_")[[1]][1]) %>% unlist()
TMP2$a2revise = lapply(TMP2$a2revise, function(x) strsplit(x,"_")[[1]][1]) %>% unlist()
TMP3 = rbind(TMP1 %>% mutate(type="area"),TMP2 %>% mutate(type="thickness"))
ggplot(TMP3, aes(x=type, y=a2revise)) + geom_point()
#only non-siignificant difference
summary(lm(a2revise~type,TMP3))

