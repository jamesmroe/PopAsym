rm(list=ls())
library("tidyverse")
library("magrittr")
library("here")
here()

#fill in annotation with heritability values

#SELECT OPTIONS ======================
#SCRIPT TAKES SUMMARY-LEVEL DATA AS INPUT AND IS THUS FULLY REPRODUCIBLE WITH THE DATA PROVIDED HERE
saveres=1 #1=yes, 0=no
typeHerit = "snp" #snp/twin
#====================================#


#UKBSNP
rm(list=ls())
b="/cluster/projects/p23/projects/11_MemP/James/UKBSNP"
rdir=here("results/heritability")
sdir=here("annotInfo")


RESAREA = read.csv(here("results/heritability/cortexH2SNParea.csv"), sep=" ")
RESCTH = read.csv(here("results/heritability/cortexH2SNPthickness.csv"), sep=" ")
RES = rbind(RESAREA, RESCTH)


#---FDR
all.equal(RES$FDR[1:500],p.adjust(RES$sig[1:500], "BH"))
all.equal(RES$FDR[501:1000],p.adjust(RES$sig[501:1000], "BH"))
all.equal(RES$bin, ifelse(RES$FDR<.05, TRUE, FALSE))


#SNP heritability for thickness was substantially lower
summary(lm((scale(h2)) ~ type, data = RES))


#percentage of significant SNP h2 after FDR per type
sum(RES$bin[RES$type=="area"])/length(RES$type[RES$type=="area"]) #53% area
sum(RES$bin[RES$type=="thickness"])/length(RES$type[RES$type=="thickness"]) #11% cth



#open vtx annot file and replace codes with corvals
V=read.csv(here("annotInfo/lh.vtxs.lh.Schaefer2018_1000Parcels_17Networks_order.40tvs8ref_sym_20.txt"),stringsAsFactors = F) #vertex code file
A=read.csv(here("annotInfo/lh.annot.lh.Schaefer2018_1000Parcels_17Networks_order.40tvs8ref_sym_20.txt"),stringsAsFactors = F) #annot code file
D=RES$name


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
  write.table(VVA,here("results/heritability/mapH2UKBSNP.area.csv"),row.names=F,col.names = F,quote=F)
  write.table(SSA,here("results/heritability/mapSigFDRUKBSNP.area.csv"),row.names=F,col.names = F,quote=F)

  write.table(VVT,here("results/heritability/mapH2UKBSNP.thickness.csv"),row.names=F,col.names = F,quote=F)
  write.table(SST,here("results/heritability/mapSigFDRUKBSNP.thickness.csv"),row.names=F,col.names = F,quote=F)
}


#UKBSNP results labels FDR
RESAREA$name[RESAREA$FDR<.05]
RESAREA$FDR[RESAREA$FDR<.05]
RESAREA$name[RESAREA$FDR==min(RESAREA$FDR)]
RESAREA[RESAREA$h2==max(RESAREA$h2),]
RESCTH[RESCTH$h2==max(RESCTH$h2),]

RESCTH$name[RESCTH$FDR<.05]
RESCTH$FDR[RESCTH$FDR<.05]
RESCTH$name[RESCTH$FDR==min(RESCTH$FDR)]

#write to file copy labels to AREATWINSIG05
# writeLines(RESAREA$name[RESAREA$FDR<.05],file.path(rdir,"resAreaUKBSNPSigFDR.txt"))
# writeLines(RESCTH$name[RESCTH$FDR<.05],file.path(rdir,"resThicknessUKBSNPSigFDR.txt"))


#bash code
# sdir="/cluster/projects/p23/projects/11_MemP/James/UKBSNP/shared/annotInfo"
# rdir="/cluster/projects/p23/projects/11_MemP/James/UKBSNP/results.SNPHerit/results.schaf"
# for i in $(cat resAreaUKBSNPSigFDR.txt); do cp $sdir/$i.mgh.label $rdir/resAreaUKBSNPSigFDR.$i.mgh.label; done
# for i in $(cat resThicknessUKBSNPSigFDR.txt); do cp $sdir/$i.mgh.label $rdir/resThicknessUKBSNPSigFDR.$i.mgh.label; done





#FILL RESULTS OF TWIN ANALYSIS RAN LOCALLY ON MAC
#open vtx annot file and replace codes with vals
rm(list=ls())
b="/cluster/projects/p23/projects/11_MemP/James/UKBSNP"
rdir=file.path(b,"results.SNPHerit/results.schaf")
sdir=file.path(b,"shared/annotInfo")

h2area=read.csv(file.path(rdir,"cortexATwinCEarea.csv"),stringsAsFactors = F,sep=" ",header=T)
h2cth=read.csv(file.path(rdir,"cortexTwinACEthickness.csv"),stringsAsFactors = F,sep=" ",header=T)
h2area$type="area"
h2cth$type="thickness"

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
VVA=V
SSA=V
VVT=V
SST=V
#vertex file
#code names
#res
#colnames


RESAREA=RES[RES$type=="area",]
RESCTH=RES[RES$type=="thickness",]

RESAREA$FDR=p.adjust(RESAREA$pA,method="BH")
RESCTH$FDR=p.adjust(RESCTH$pA,method="BH")

#order by h2
tmp=RESAREA[
  with(RESAREA, order(a2)),
]

RESAREA[RESAREA$a2==max(RESAREA$a2),]
RESCTH[RESCTH$a2==max(RESCTH$a2),]

for (i in 1:length(C)) {
  co=C[i]
  struct=S[i]
  if (i==1) {
    VVA[which(VVA==co),]=0
    VVT[which(VVT==co),]=0
    SSA[which(SSA==co),]=0
    SST[which(SST==co),]=0
  } else {
    VVA[which(VVA==co),]= RESAREA$a2[which(RESAREA$name==struct)]
    SSA[which(SSA==co),]= -log10(RESAREA$FDR[which(RESAREA$name==struct)])
    
    VVT[which(VVT==co),]= RESCTH$a2[which(RESCTH$name==struct)]
    SST[which(SST==co),]= -log10(RESCTH$FDR[which(RESCTH$name==struct)])
  }
}  

write.table(VVA,file.path(rdir,"mapH2cortexACEarea.csv"),row.names=F,col.names = F,quote=F)
write.table(SSA,file.path(rdir,"mapSigcortexACEarea.csv"),row.names=F,col.names = F,quote=F)

write.table(VVT,file.path(rdir,"mapH2cortexACEthickness.csv"),row.names=F,col.names = F,quote=F)
write.table(SST,file.path(rdir,"mapSigcortexACEthickness.csv"),row.names=F,col.names = F,quote=F)



#twin results labels sig before FDR
RESAREA$name[RESAREA$FDR<.05]
RESAREA$pA[RESAREA$pA<.05]
RESAREA$name[RESAREA$pA<.05]

RESCTH$name[RESCTH$FDR<.05]
RESCTH$pA[RESCTH$pA<.05]
RESCTH$name[RESCTH$pA<.05]


#write to file copy labels to AREATWINSIG05
writeLines(RESAREA$name[RESAREA$pA<.05],file.path(rdir,"resAreaTwinSig05-cortexACE.txt"))
writeLines(RESCTH$name[RESCTH$pA<.05],file.path(rdir,"resThicknessTwinSig05-cortexACE.txt"))


#bash code to find labels easily
# sdir="/cluster/projects/p23/projects/11_MemP/James/UKBSNP/shared/annotInfo"
# rdir="/cluster/projects/p23/projects/11_MemP/James/UKBSNP/results.SNPHerit/results.schaf"
# for i in $(cat resAreaTwinSig05-cortexACE.txt); do cp $sdir/$i.mgh.label $rdir/FINALresAreaTwinSig05.$i.mgh.label; done
# for i in $(cat resThicknessTwinSig05-cortexACE.txt); do cp $sdir/$i.mgh.label $rdir/FINALresThicknessTwinSig05.$i.mgh.label; done
