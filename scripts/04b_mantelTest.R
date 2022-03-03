# restart R
# install.packages("here")
library("here")
here()


#SELECT METRIC=========
#SCRIPT TAKES SUMMARY-LEVEL DATA AS INPUT AND IS THUS FULLY REPRODUCIBLE WITH THE DATA PROVIDED HERE
# brainvar="area"
brainvar="thickness"
#=====================#


# LOAD PACKAGES -----------------------------------------------------------
tmp.packages = c("corrplot","viridis","ade4")
tmpnew.packages = tmp.packages[!(tmp.packages %in% installed.packages()[,"Package"])]
if(length(tmpnew.packages)) {
  install.packages(tmpnew.packages)
}
sapply(tmp.packages, require, character.only = T)
rm(list=ls(pattern="tmp*"))



# LOAD IN CORRELATION MATRICES FOR MANTEL TEST ---------------------------
#LCBC
load(here(
  "results/structuralCovar",
  paste0("cormat.LCBC.", brainvar, "-publish.Rda")
))
M1 = cormat
M1icv = cormatICV
corrplot(M1,method="color",type="full",addCoef.col = "black",col=colorRampPalette(c("blue","white","orange2"))(25),tl.pos="n",tl.col = "blue",cl.cex = 0.8,cl.pos = "r",tl.cex = 0.8,number.cex = 0.8)
rm(cormat, cormatICV)

#UKB
load(here(
  "results/structuralCovar",
  paste0("cormat.UKB.", brainvar, "-publish.Rda")
))
M2 = cormat
M2icv = cormatICV
corrplot(M2,method="color",type="full",addCoef.col = "black",col=colorRampPalette(c("blue","white","orange2"))(25),tl.pos="n",tl.col = "blue",cl.cex = 0.8,cl.pos = "r",tl.cex = 0.8,number.cex = 0.8)
rm(cormat, cormatICV)

#HCP
load(here(
  "results/structuralCovar",
  paste0("cormat.HCP.", brainvar, "-publish.Rda")
))
M3 = cormat
M3icv = cormatICV
corrplot(M3,method="color",type="full",addCoef.col = "black",col=colorRampPalette(c("blue","white","orange2"))(25),tl.pos="n",tl.col = "blue",cl.cex = 0.8,cl.pos = "r",tl.cex = 0.8,number.cex = 0.8)
rm(cormat, cormatICV)


# TEST WHETHER ADDITIONALLY CORRECTING FOR ICV CHANGES INTERRELATIONS ---------------------------
diff1=M1-M1icv
diff2=M2-M2icv
diff3=M3-M3icv
corrplot(diff1,method="color",col=viridis(50,option="C"),type="full",addCoef.col = "black",tl.pos="n",tl.col = "blue",cl.cex = 0.8,cl.pos = "r",tl.cex = 0.8,number.cex = 0.8)
corrplot(diff2,method="color",col=viridis(50,option="C"),type="full",addCoef.col = "black",tl.pos="n",tl.col = "blue",cl.cex = 0.8,cl.pos = "r",tl.cex = 0.8,number.cex = 0.8)
corrplot(diff3,method="color",col=viridis(50,option="C"),type="full",addCoef.col = "black",tl.pos="n",tl.col = "blue",cl.cex = 0.8,cl.pos = "r",tl.cex = 0.8,number.cex = 0.8)
cor(as.vector(M1),as.vector(M1icv))
cor(as.vector(M2),as.vector(M2icv))
cor(as.vector(M3),as.vector(M3icv))
max(abs(diff1))
max(abs(diff2))
max(abs(diff3))
print(paste("None of the relationships could be explained by brain size", 
      "as additionally removing the effect of intracranial volume (ICV) from cluster AI’s had a negligible effect on their interrelations"))

print(paste("max correlation change =",
max(c(max(abs(diff1)),
    max(abs(diff2)),
    max(abs(diff3)))
    )))



# RUN PAIRWISE MANTEL TESTS -----------------------------------------------
for (i in 1:3) {
  if (i == 1) {
    RR=list()
    distmantel=rdistmantel=list()
    print("LCBC v UKB")
    y=M1
    x=M2
  } else if (i == 2) {
    print("LCBC v HCP")
    y=M1
    x=M3
  } else if (i==3) {
    print("UKB v HCP")
    y=M2
    x=M3
  }
  pmat1=y
  pmat2=x
  mat1=as.matrix(pmat1)
  mat2=as.matrix(pmat2)
  #run mantel on distance matrix
  distmat1 <- as.dist(mat1)
  distmat2 <- as.dist(mat2)
  qdistmat1 <- quasieuclid(as.dist(distmat1))
  qdistmat2 <- quasieuclid(as.dist(distmat2))
  
  RR[[i]] = mantel.rtest(qdistmat1,qdistmat2,nrepet=10000)
  print(RR[[i]])
  # plot(r1 <- mantel.rtest(qdistmat1,qdistmat2), main = "Mantel's test")
  # r1
  # plot(qdistmat1,qdistmat2)
  # cor.test(qdistmat1,qdistmat2)
}

