# setwd("/Users/jamesroe/Dropbox/GitHub/PopAsym")
library("here")
library("corrplot")
library("tidyverse")
library("magrittr")
library("viridis")
here()

#SELECT OPTIONS ======================
saveres = 1
#====================================#

RR = read.csv(here("results/heritability/clustersGENCORR.csv"), stringsAsFactors = F, sep = " ")
AREA = RR[RR$meas == "area",]
THICK = RR[RR$meas == "thickness",]


#area = 78 tests
#thick = 48 tests
#pairwise tests only when both traits achieving p<.05 uncorrected SNP-heritability
nrow(AREA)
nrow(THICK)


#confirm FDR-correction
all.equal(AREA$fdr, p.adjust(AREA$pval,method="BH"))
all.equal(THICK$fdr, p.adjust(THICK$pval,method="BH"))


# RETRIEVE AREA GENETIC CORRELATIONS -----------------------------------------------
#grab FDR-corrected results for area
pp1=which(p.adjust(AREA$pval,method="BH")<0.05)
AREA[pp1,]

#create corr matrix to fill
MATarea = matrix(nrow = 14, ncol = 14)
MATarea = lower.tri(MATarea, diag = T)
diag(MATarea)=1
rMATarea=porigMATarea=MATarea #create matrices to store rG, origPvals, fdrPvals
for (i in pp1) {
  t1 = as.integer(AREA[i,]$t1)
  t2 = as.integer(AREA[i,]$t2)
  gP = (AREA[i,]$fdr)
  gPorig = (AREA[i,]$pval)
  rG = (AREA[i,]$rG)
  MATarea[t1,t2]=-log10(gP)
  rMATarea[t1,t2]=rG
  porigMATarea[t1,t2]=-log10(gPorig)
}

#duplicate to lower tri
rMATarea[lower.tri(rMATarea)] <- t(rMATarea)[lower.tri(rMATarea)]
MATarea[lower.tri(MATarea)] <- t(MATarea)[lower.tri(MATarea)]
porigMATarea[lower.tri(porigMATarea)] <- t(porigMATarea)[lower.tri(porigMATarea)]


#VIP - opposite hemisphere genetic correlations need to be sign inversed
#as AI's for rightward clusters are negative (more asymmetry)
#hence negative genetic correlations between leftward and rightward clusters represent asymmetry-asymmetry relationships
corrplot(rMATarea,method="color",col=viridis(50),type="full",addCoef.col = "black",tl.col = "blue",is.corr = F,tl.pos="n")
#opposite correlations in lower left quadrant / upper right quadrant
rMATarea[(8:14), (1:7)] = rMATarea[(8:14), (1:7)] * -1
rMATarea[(1:7), (8:14)] = rMATarea[(1:7), (8:14)] * -1


#rG
corrplot(rMATarea,method="color",col=viridis(50),type="full",addCoef.col = "black",tl.col = "blue",is.corr = F,tl.pos="n")
#-log10(p) [uncorrected]
corrplot(porigMATarea,method="color",col=viridis(50),type="full",addCoef.col = "black",tl.col = "blue",is.corr = F,tl.pos="n")
#-log10(p) [FDR-corrected]
corrplot(MATarea,method="color",col=viridis(50),type="full",addCoef.col = "black",tl.col = "blue",is.corr = F,tl.pos="n")
#set 0s (not tested) to NA
rMATarea[which(rMATarea==0)] = NA


#save fig
if (saveres == 1) {
  filename = here("results/heritability/gencormat.area.png")
  print(paste("saving", filename))
  png(filename = filename, width = 15, height = 15, units = "cm", res = 300)
  corrplot(rMATarea,method="circle",type="lower",
           addCoef.col = "black",col=colorRampPalette(c("blue","white","orange2"))(25),tl.pos="n",na.label=" ",tl.col = "blue",cl.cex = 0.8,cl.pos = "r",tl.cex = 0.8,number.cex = 0.8)
  dev.off()
}

#store min rG for AREA to get same colour/size code for THICKNESS plot below
minarea = min(rMATarea,na.rm = T)



# RETRIEVE THICKNESS GENETIC CORRELATIONS -----------------------------------------------
#grab FDR-corrected results for thickness
pp2 = which(p.adjust(THICK$pval,method="BH")<0.05)
THICK[pp2,]

#create corr matrix to fill
MATthick = matrix(nrow = 20, ncol = 20)
MATthick = lower.tri(MATthick, diag = T)
diag(MATthick)=1
rMATthick=porigMATthick=MATthick #create matrices to store rG, origPvals, fdrPvals
for (i in pp2) {
  t1 = as.integer(THICK[i,]$t1)
  t2 = as.integer(THICK[i,]$t2)
  gP = (THICK[i,]$fdr)
  gPorig = (THICK[i,]$pval)
  rG = (THICK[i,]$rG)
  MATthick[t1,t2]=-log10(gP)
  rMATthick[t1,t2]=rG
  porigMATthick[t1,t2]=-log10(gPorig)
}

#duplicate to lower tri
rMATthick[lower.tri(rMATthick)] <- t(rMATthick)[lower.tri(rMATthick)]
MATthick[lower.tri(MATthick)] <- t(MATthick)[lower.tri(MATthick)]
porigMATthick[lower.tri(porigMATthick)] <- t(porigMATthick)[lower.tri(porigMATthick)]

#one surviving result is not an opposite asymmetry correlation
THICK[pp2,]
#rG
corrplot(rMATthick,method="color",col=viridis(50),type="full",addCoef.col = "black",tl.col = "blue",is.corr = F,tl.pos="n")
#-log10(p) [uncorrected]
corrplot(porigMATthick,method="color",col=viridis(50),type="full",addCoef.col = "black",tl.col = "blue",is.corr = F,tl.pos="n")
#-log10(p) [FDR-corrected]
corrplot(MATthick,method="color",col=viridis(50),type="full",addCoef.col = "black",tl.col = "blue",is.corr = F,tl.pos="n")
#set 0s (not tested) to NA
rMATthick[which(rMATthick==0)] = NA


#add in min rG for area to lower tri (not visualized) to get same colour/size bar as areal results 
rMATthick[20,1] = minarea
corrplot(rMATthick,method="circle",col=viridis(50),type="lower",
         addCoef.col = "black",tl.col = "blue",is.corr = F,tl.pos="n",na.label=" ")
corrplot(rMATthick,method="circle",col=viridis(50),type="upper",tl.col = "blue",is.corr = F,tl.pos="n",na.label=" ")


if (saveres == 1) {
  filename = here("results/heritability/gencormat.thick.png")
  print(paste("saving", filename))
  png(filename = filename, width = 15, height = 15, units = "cm", res = 300)
  corrplot(rMATthick,method="circle",type="upper",
           # addCoef.col = "black",
           col=colorRampPalette(c("blue","white","orange2"))(25),tl.pos="n",na.label=" ",tl.col = "blue",cl.cex = 0.8,cl.pos = "r",tl.cex = 0.8,number.cex = 0.8) #,order = "hclust",addrect = 3,rect.col = "black")
  dev.off()
}

