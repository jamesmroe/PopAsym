library("tidyverse")
library("magrittr")
library("here")
here()


#SELECT OPTIONS ======================
args = commandArgs(TRUE)
datdir=as.character(args[1]) #SCRIPT ANALYZES INDIVIDUAL AND FAMILY-LEVEL DATA THAT IS RESTRICTED TO COMPLY WITH HCP DATA USAGE AGREEMENTS (see data availability section in manuscript)
datdir="/Users/jamesroe/LCBC/Users/jamesroe/PHD_project/Paper3/data"
saveres = 0 #1 = save, 0 = nosave
schaefer = 0 #1 = cortex-wide, 0 = cluster-wise

brainvar="area" #select metric
# brainvar="thickness"
#====================================#
runIndividuallevelAnalysis = 1

if (runIndividuallevelAnalysis == 1) {
  print("running indivdual level analysis on restricted HCP data")
  
  #read in family structure
  #TWINS
  twinpairs=read.csv(file.path(datdir,"twinpairs_hcp.csv"),header=T,stringsAsFactors = F,sep=" ")
  twins1sib=read.csv(file.path(datdir,"twins1sib_hcp.csv"),header=T,stringsAsFactors = F,sep=" ")
  twins2sibs=read.csv(file.path(datdir,"twins2sibs_hcp.csv"),header=T,stringsAsFactors = F,sep=" ")
  twins3sibs=read.csv(file.path(datdir,"twins3sibs_hcp.csv"),header=T,stringsAsFactors = F,sep=" ")
  twins4sibs=read.csv(file.path(datdir,"twins4sibs_hcp.csv"),header=T,stringsAsFactors = F,sep=" ")
  
  #SIBS
  sibpairs=read.csv(file.path(datdir,"sibpairs_hcp.csv"),header=T,stringsAsFactors = F,sep=" ")
  sibtrips=read.csv(file.path(datdir,"sibtrips_hcp.csv"),header=T,stringsAsFactors = F,sep=" ")
  sibquads=read.csv(file.path(datdir,"sibquads_hcp.csv"),header=T,stringsAsFactors = F,sep=" ")
  sibquints=read.csv(file.path(datdir,"sibquints_hcp.csv"),header=T,stringsAsFactors = F,sep=" ")
  
  
  
  #ID log
  I=c(twinpairs$twin1,twinpairs$twin2,
      twins1sib$twin1,twins1sib$twin2,twins1sib$Sib1,
      twins2sibs$twin1,twins2sibs$twin2,twins2sibs$Sib1,twins2sibs$Sib2,
      twins3sibs$twin1,twins3sibs$twin2,twins3sibs$Sib1,twins3sibs$Sib2,twins3sibs$Sib3,
      twins4sibs$twin1,twins4sibs$twin2,twins4sibs$Sib1,twins4sibs$Sib2,twins4sibs$Sib3,twins4sibs$Sib4,
      sibpairs$Sib1,sibpairs$Sib2,
      sibtrips$Sib1,sibtrips$Sib2,sibtrips$Sib3,
      sibquads$Sib1,sibquads$Sib2,sibquads$Sib3,sibquads$Sib4,
      sibquints$Sib1,sibquints$Sib2,sibquints$Sib3,sibquints$Sib4,sibquints$Sib5)
  length(I)
  
  #read in demographics
  DF = read.csv(file.path(datdir,"demog_hcp.csv"), stringsAsFactors = F, header = T, sep="")
  
  
  #load in base data if cluster-wise. If cortex-wide does this within loop below
  if (schaefer != 1) {
    DD = read.csv(file.path(datdir,paste0(brainvar,"measures_hcp.csv")),header=T,stringsAsFactors = F,sep=" ")
    identical(DF$Subject, DD$Subject)
    
    if (brainvar == "global") {
      ll = 8
    } else {
      nlabs=names(DD)[grepl("6n",names(DD), ignore.case = F)]
      plabs=names(DD)[grepl("6p",names(DD), ignore.case = F)]
      labs = c(nlabs, plabs)
      ll = length(labs)
    }
  } else {
    labs=c(list.files(file.path(datdir,paste0("split.data.schaefer.",brainvar)),pattern="AI"))
    ll = length(labs)
  }
  
  if (brainvar =="area") {
    roinames = c(
      "R_Parieto-occipital_sulcus",
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
      "L_Superior-frontal-gyrus"
    )
  } else if (brainvar == "thickness") {
    roinames=c(
      "R_Superior-temporal-sulcus",
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
      "L_Calcarine-sulcus"
    )
  }
  print(paste("running twin bivariate models for", brainvar))
  

  c = 0
  for (j in 1:ll) {
    for (jj in ((j+1):(ll))) {
      if ((j+1)==(ll+1)) {
        print("finished")
        break
      }
      if (j == jj) {
        print("skipping")
        next
      }
      
      c=c+1
      if (j==1 & jj==2) {
        trait1=c()
        trait2=c()
        rG=rGp=c()
        rGlwr=c()
        rGupr=c()
        CIout=modelComp=list()
        TT1=c()
        TT2=c()
        pRG=c()
        h2_1=h2_2=c()
      }
      trait1[c]=sprintf("%02d", j)
      trait2[c]=sprintf("%02d", jj)
      roi1=roinames[j]; roi2=roinames[jj]
      print(paste("trait",trait1[c], "v", trait2[c]))
      
      
      #LOAD PRECORRECTED TRAITS
      #reordered by family structure
      idir="/Users/jamesroe/LCBC/Users/jamesroe/PHD_project/Paper3/data/correctedLabs4Gencor"
      pattern1 = paste0(roi1, ".*", brainvar)
      pattern2 = paste0(roi2, ".*", brainvar)
      
      T1 = read.csv(list.files(idir, pattern = pattern1, full.names = T),header=T,stringsAsFactors = F,sep=" ")
      T2 = read.csv(list.files(idir, pattern = pattern2, full.names = T),header=T,stringsAsFactors = F,sep=" ")
      
      # #inverse if rightward #already done now
      # if (substr(roi1,1,1) == "R") {
      #   print(paste("inversing rightward ROI trait1. Mean =", mean(T1$Ycor)))
      #   T1$Ycor = T1$Ycor * -1
      # }
      # if (substr(roi2,1,1) == "R") {
      #   print(paste("inversing rightward ROI trait2. Mean =", mean(T2$Ycor)))
      #   T2$Ycor = T2$Ycor * -1
      # }
      
      #prep twins
      #add missing col
      print("prepping data")
      twinpairs$newZyg2=twinpairs$newZyg1 #add new col to match up cols across family structs
      twinpairs %<>% select(Family_ID,twin1,twin2,newZyg1,newZyg2,everything())
      
      fieldnames=names(twins4sibs)
      nfield=length(fieldnames)
      
      #same nCols in all
      addCols = function(df,X) {
        df=cbind(df,t(rep(NA,X)))
        names(df) [length(names(df))-X+1:X] = fieldnames[length(names(df))-X+1:X]
        return(df)
      }
      
      X=nfield-length(names(twinpairs))
      twinpairsX = addCols(twinpairs,X)
      
      X=nfield-length(names(twins1sib))
      twins1sibX = addCols(twins1sib,X)
      
      X=nfield-length(names(twins2sibs))
      twins2sibsX = addCols(twins2sibs,X)
      
      X=nfield-length(names(twins3sibs))
      twins3sibsX = addCols(twins3sibs,X)
      
      
      #prep Sibs
      fieldnames=names(sibquints)
      nfield=length(fieldnames)
      
      X=nfield-length(names(sibpairs))
      sibpairsX = addCols(sibpairs,X)
      
      X=nfield-length(names(sibtrips))
      sibtripsX = addCols(sibtrips,X)
      
      X=nfield-length(names(sibquads))
      sibquadsX = addCols(sibquads,X)
      
      
      #Put it all together
      sibs = rbind(sibpairsX, sibtripsX, sibquadsX, sibquints)
      twins = rbind(twinpairsX, twins1sibX, twins2sibsX, twins3sibsX, twins4sibs)
      twinIDs = twins %>% select(twin1,twin2,Sib1,Sib2,Sib3,Sib4)
      
      #split by zygosity
      MZdat = twins %>% filter(newZyg1=="MZ")
      DZdat = twins %>% filter(newZyg1=="DZ")
      MZids = MZdat %>% select(twin1,twin2,Sib1,Sib2,Sib3,Sib4) %>% unname()
      DZids = DZdat %>% select(twin1,twin2,Sib1,Sib2,Sib3,Sib4) %>% unname()
      
      sibIDs = sibs %>% select(Sib1,Sib2,Sib3,Sib4,Sib5)
      sibIDs$Sib0 = rep(NA,nrow(sibIDs)) #col 1 reserved for a twin
      sibIDs %<>% select(Sib0,everything())
      names(MZids) = c(1:6)
      
      #grab 1:row num of mz
      lenMZ=dim(MZids)[1]
      lenDZ=dim(DZids)[1]
      names(DZids) = c(1:6)
      names(sibIDs) = c(1:6)
      mtr=rbind(MZids,DZids,sibIDs)
      mtr=as.matrix(mtr %>% unname())
      
      
      #populate matrix of IDs with phenotypic vals
      populate_matrix <- function(ids, values, matrix_template) {
        # Stop if ids and values are not same length.
        if(length(ids) != length(values)) 
          stop("ids and values are not the same length", call. = FALSE)
        # match ids to positions in the matrix
        idx <- match(ids, matrix_template)
        # Create a clean copy of the matrix
        val_mtr <- matrix(NA, 
                          nrow = nrow(matrix_template), 
                          ncol = ncol(matrix_template))
        # populate index locations with values
        val_mtr[idx] <- values
        return(val_mtr)
      }
      mtrX1 <- populate_matrix(T1$ID, T1$Ycor, as.matrix(mtr))
      mtrX2 <- populate_matrix(T2$ID, T2$Ycor, as.matrix(mtr))
      
      #check worked
      # head(mtr)
      # T1[which(T1$ID==139637),]
      # T1[which(T1$ID==196750),]
      # head(mtrX1)
      # tail(mtr)
      # T1[which(T1$ID==959574),]
      # tail(mtrX1)
      
      
      DataMZplus1=mtrX1[1:lenMZ,]
      DataDZplus1=mtrX1[(lenMZ+1):nrow(mtrX1),]
      DataMZplus2=mtrX2[1:lenMZ,]
      DataDZplus2=mtrX2[(lenMZ+1):nrow(mtrX2),]
      
      mtrMZ=mtr[1:lenMZ,]
      mtrDZ=mtr[(lenMZ+1):nrow(mtr),]
      
      
      selVars <- c('v1','v2','v3','v4','v5','v6')
      dimnames(DataMZplus1) <- list(NULL,selVars)
      dimnames(DataDZplus1) <- list(NULL,selVars)
      dimnames(DataMZplus2) <- list(NULL,selVars)
      dimnames(DataDZplus2) <- list(NULL,selVars)
      
      
      #EXCLUDE 2 HCP OUTLIERS
      #identify subject exclusions in mtr and delete from mtrx
      mtrX1=mtrX1[-which(mtr==132118, arr.ind = TRUE)[1],]
      mtrX1=mtrX1[-which(mtr==894067, arr.ind = TRUE)[1],]
      DataMZplus1=DataMZplus1[-which(mtrMZ==132118, arr.ind = TRUE)[1],]
      DataDZplus1=DataDZplus1[-which(mtrDZ==894067, arr.ind = TRUE)[1],]
      
      mtrX2=mtrX2[-which(mtr==132118, arr.ind = TRUE)[1],]
      mtrX2=mtrX2[-which(mtr==894067, arr.ind = TRUE)[1],]
      DataMZplus2=DataMZplus2[-which(mtrMZ==132118, arr.ind = TRUE)[1],]
      DataDZplus2=DataDZplus2[-which(mtrDZ==894067, arr.ind = TRUE)[1],]
      # mtr[73,]
      # mtr[256,]
      # length(mtrX[which(!is.na(mtrX), arr.ind = TRUE)])
      # dim(mtrX)
      
      
      
      #RESCALE DUE TO LOW VARIANCE
      DataMZplus1_sc = scale(DataMZplus1,
                            rep(mean(unlist(DataMZplus1), na.rm = T), ncol(DataMZplus1)),
                            rep(sd(unlist(DataMZplus1), na.rm = T), ncol(DataMZplus1)))
      
      DataDZplus1_sc = scale(DataDZplus1,
                            rep(mean(unlist(DataDZplus1), na.rm = T), ncol(DataDZplus1)),
                            rep(sd(unlist(DataDZplus1), na.rm = T), ncol(DataDZplus1)))
      
      DataMZplus2_sc = scale(DataMZplus2,
                             rep(mean(unlist(DataMZplus2), na.rm = T), ncol(DataMZplus2)),
                             rep(sd(unlist(DataMZplus2), na.rm = T), ncol(DataMZplus2)))
      
      DataDZplus2_sc = scale(DataDZplus2,
                             rep(mean(unlist(DataDZplus2), na.rm = T), ncol(DataDZplus2)),
                             rep(sd(unlist(DataDZplus2), na.rm = T), ncol(DataDZplus2)))
      
      
      
      MZconcat = data.frame(DataMZplus1_sc,DataMZplus2_sc)
      DZconcat = data.frame(DataDZplus1_sc,DataDZplus2_sc)
      names(MZconcat) = c(paste0("trait1_", 1:6), paste0("trait2_", 1:6))
      names(DZconcat) = c(paste0("trait1_", 1:6), paste0("trait2_", 1:6))
      library("OpenMx")
      ff = mxFitFunctionML()
      print("running AE")
      twinAEbivarModel <- mxModel("twinAE", 
                              mxMatrix("Full", nrow=1, ncol=12, free=TRUE, values=0, label=c(rep("mean1", 6), rep("mean2", 6)), name="expMean"),
                              # Matrix expMean for expected mean 
                              # vector for MZ and DZ twins    
                              # -------------------------------------
                              
                              mxMatrix("Lower", nrow=2, ncol=2, free=c(T, T, T), values = sqrt(c(.1, 0.01, .1)), label= c("a11", "a21", "a22"), name="X"), #all parameters specified
                              mxMatrix("Lower", nrow=2, ncol=2, free=F, values=0 * diag(2), label= c("c11", "c21", "c22"), name="Y"), #all parameters specified
                              mxMatrix("Lower", nrow=2, ncol=2, free=c(T, T, T), values = sqrt(c(.9, 0, .9)), label= c("e11", "e21", "e22"), name="Z"), #all parameters specified
                              # Matrices X, Y, and Z to store the 
                              # a, c, and e path coefficients
                              # -------------------------------------
                              
                              mxAlgebra(X %*% t(X), name="A"),
                              mxAlgebra(Y %*% t(Y), name="C"),
                              mxAlgebra(Z %*% t(Z), name="E"),	
                              # Matrixes A, C, and E to compute 
                              # A, C, and E variance components
                              # -------------------------------------
                              mxAlgebra(rbind(cbind(1, 1, .5, .5, .5, .5),
                                              cbind(1, 1, .5, .5, .5, .5),
                                              cbind(.5, .5, 1, .5, .5, .5),
                                              cbind(.5, .5, .5, 1, .5, .5),
                                              cbind(.5, .5, .5, .5, 1, .5),
                                              cbind(.5, .5, .5, .5, .5, 1)),
                                        name="rgMZ"),
                              mxMatrix("Id", 6, 6, F, name = "rE"),
                              
                              mxAlgebra(A %x% rgMZ + E %x% rE, name="expCovMZ"),
                              
                              # Matrix expCOVMZ for expected 
                              # covariance matrix for MZ twins
                              # -------------------------------------
                              
                              mxAlgebra(rbind(rbind(cbind(1, .5, .5, .5, .5, .5),
                                                    cbind(.5, 1, .5, .5, .5, .5),
                                                    cbind(.5, .5, 1, .5, .5, .5),
                                                    cbind(.5, .5, .5, 1, .5, .5),
                                                    cbind(.5, .5, .5, .5, 1, .5),
                                                    cbind(.5, .5, .5, .5, .5, 1))), 
                                        name="rgDZ"),
                              mxAlgebra(A %x% rgDZ + E %x% rE, name="expCovDZ"),
                              # Matrix expCOVMZ for expected 
                              # covariance matrix for DZ twins
                              # -------------------------------------
                              
                              # MZ Model
                              # -------------------------------------
                              mxModel("MZ",
                                      mxData(MZconcat, type="raw"),
                                      mxExpectationNormal("twinAE.expCovMZ", "twinAE.expMean",names(MZconcat)),
                                      mxFitFunctionML()),
                              
                              # DZ Model
                              # -------------------------------------
                              mxModel("DZ", 
                                      mxData(DZconcat, type="raw"),
                                      mxExpectationNormal("twinAE.expCovDZ", "twinAE.expMean",names(DZconcat)),
                                      mxFitFunctionML()),
                              
                              
                              # # Just easier
                              mxFitFunctionMultigroup(c("MZ", "DZ")),
                              # mxAlgebra(MZ.objective + DZ.objective, name="twin"),
                              # mxFitFunctionAlgebra("twin")
                              mxAlgebra(A[1, 1] / (A[1, 1] + E[1, 1]), name = "h2_trait1"),
                              mxAlgebra(A[2, 2] / (A[2, 2] + E[2, 2]), name = "h2_trait2"),
                              mxAlgebra(A[1, 2] / sqrt(A[1, 1] * A[2, 2]), name = "rg"),
                             # 
                              mxCI(c("h2_trait1", "h2_trait2", "rg"))
                              
      )
      twinAEFitRG <- mxTryHard(twinAEbivarModel, intervals = T)
      summary(twinAEFitRG)
      SS=summary(twinAEFitRG)
      
      #update model to test signifcance
      twinAE = omxSetParameters(twinAEbivarModel, labels = "a21", free = F, values = 0, name = "noRG")
      twinAEFit <- mxTryHard(twinAE, intervals = T)
      
      modelComp[[c]] = mxCompare(twinAEFitRG, twinAEFit)
      anova(twinAEFitRG, twinAEFit)
      rGp[c] = modelComp[[c]]$p[2]
        
      estA = mxEval(A, twinAEFitRG)
      estAnew = mxEval(A, twinAEFit)
      mxEval(A %x% rgMZ, twinAEFitRG)
      estE = mxEval(E, twinAEFitRG)
      estP = estA + estE
      estA  / estP
      cov2cor(estA)
      cov2cor(estE)
      
      h2_1[c]=mxEval(h2_trait1, twinAEFitRG)
      h2_2[c]=mxEval(h2_trait2, twinAEFitRG)
      
      #results out
      rG[c]=mxEval(rg, twinAEFitRG)
      rGlwr[c]=SS$CI[3,1]
      rGupr[c]=SS$CI[3,3]
      CIout[[c]]=SS$CI
      TT1[c] = roi1
      TT2[c] = roi2
      
    }
  }
  RRtwin = data.frame(rG,rGp,rGlwr,rGupr,trait1,trait2,TT1,TT2,h2_1=round(h2_1,3),h2_2=round(h2_2,3))
  
  if (brainvar == "area") {
    #remove one cluster not testing for area (non sig SNP-based h2; SI File2 Table 1H)
    RRtwinOut = RRtwin %>% filter(trait1 != "04" & trait2 != "04")
  } else if (brainvar == "thickness") {
    #remove all clusters not testing for thickness (non sig SNP-based h2; SI File2 Table 1L)
    RRtwinOut = RRtwin %>% filter(trait1 != "03" & trait2 != "03",
                                  trait1 != "08" & trait2 != "08",
                                  trait1 != "11" & trait2 != "11",
                                  trait1 != "12" & trait2 != "12",
                                  trait1 != "13" & trait2 != "13",
                                  trait1 != "14" & trait2 != "14",
                                  trait1 != "17" & trait2 != "17",
                                  trait1 != "18" & trait2 != "18",
                                  trait1 != "19" & trait2 != "19"
                                  )
  }
  
  #ensure ROI-pairs match with SNP-based
  # TMP = RRtwinOut
  # names(TMP)[5:8]=c("t1","t2","trait1","trait2")
  # names(TMP)[1]="rG_twin"
  # RR = read.csv(here("results/heritability/clustersGENCORR.csv"), stringsAsFactors = F, sep = " ")
  # AREA = RR[RR$meas == "area",]
  # THICK = RR[RR$meas == "thickness",]
  # THICK$t1 = sprintf("%02d", THICK$t1)
  # THICK$t2 = sprintf("%02d", THICK$t2)
  # TMPTMP = merge(TMP, THICK)
  
  #2 FDR corrected results
  RRtwinOut[which(p.adjust(RRtwinOut$rGp,method="BH")<0.05),]
  RRtwinOut$fdr = p.adjust(RRtwinOut$rGp,method="BH")
  
  if (saveres == 1) {
    print("saving output")
    write.table(RRtwinOut, here("results/heritability",paste0("clustersGENCORRtwin",brainvar,".csv")),quote=F,row.names = F,col.names = T)
  }
  RESERVE = RRtwinOut

} else {
  
  print("reproducing area plot")
  
  # From here takes file and reproduces plot --------------------------------
  RRtwinOut = read.csv(here("results/heritability",paste0("clustersGENCORRtwin",brainvar,".csv")), sep = " ", stringsAsFactors = F)
  
  #rows with FDRsig results
  # pp1 = RRtwinOut[which(p.adjust(RRtwinOut$rGp,method="BH")<0.05),]
  #rows with sig results
  pp1 = RRtwinOut[which(RRtwinOut$rGp<0.05),]
  
  
  #create corr matrix to fill
  if (brainvar == "area") {
    MAT = matrix(nrow = 14, ncol = 14)
  } else if (brainvar == "thickness") {
    MAT = matrix(nrow = 20, ncol = 20)
  }
  MAT = lower.tri(MAT, diag = T)
  diag(MAT)=1
  #create matrices to store rG, origPvals, fdrPvals
  rMAT=porigMAT=cilwrMAT=ciuprMAT=MAT
  # pp1=which(p.adjust(RR$pval,method="BH")<0.05)
  # RRtwinsig  = RRtwin[sign(RRtwin$rGlwr)==sign(RRtwin$rGupr),]
  for (i in 1:nrow(pp1)) {
    t1 = as.integer(pp1[i,]$trait1)
    t2 = as.integer(pp1[i,]$trait2)
    # gP = (pp1[i,]$rGp)
    gPorig = (pp1[i,]$rGp)
    RG = (pp1[i,]$rG)
    cilwrMAT[t1,t2] = (pp1[i,]$rGlwr)
    ciuprMAT[t1,t2] = (pp1[i,]$rGupr)
    MAT[t1,t2]=-log10(gPorig)
    rMAT[t1,t2]=RG
    porigMAT[t1,t2]=-log10(gPorig)
  }
  
  #duplicate to lower tri
  rMAT[lower.tri(rMAT)] <- t(rMAT)[lower.tri(rMAT)]
  MAT[lower.tri(MAT)] <- t(MAT)[lower.tri(MAT)]
  porigMAT[lower.tri(porigMAT)] <- t(porigMAT)[lower.tri(porigMAT)]
  cilwrMAT[lower.tri(cilwrMAT)] <- t(cilwrMAT)[lower.tri(cilwrMAT)]
  ciuprMAT[lower.tri(ciuprMAT)] <- t(ciuprMAT)[lower.tri(ciuprMAT)]
  
  library("corrplot")
  library("viridis")
  #set 0s to NA
  rMAT[which(rMAT==0)] = NA
  #set NA's on upr and lwr to max 
  #(i.e. where rG is estimated to be practically equal to 1 or -1 (e.g. -0.999999999934736), CIupr or CIlwr is NA)
  #so set these to 1 or -1
  cilwrMAT[apply(rMAT, c(1, 2), function(x) all.equal(as.numeric(x),-1) == TRUE)] = -1
  ciuprMAT[apply(rMAT, c(1, 2), function(x) all.equal(as.numeric(x),1) == TRUE)] = 1
  #then set 0s to NA
  cilwrMAT[which(cilwrMAT==0)] = NA
  ciuprMAT[which(ciuprMAT==0)] = NA
  #check worked
  corrplot(rMAT,method="circle",type="lower",
           addCoef.col = "black",col=colorRampPalette(c("blue","white","orange2"))(25),tl.pos="n",na.label=" ",tl.col = "blue",cl.cex = 0.8,cl.pos = "r",tl.cex = 0.8,number.cex = 0.8)
  corrplot(cilwrMAT,method="circle",type="lower",
           addCoef.col = "black",col=colorRampPalette(c("blue","white","orange2"))(25),tl.pos="n",na.label=" ",tl.col = "blue",cl.cex = 0.8,cl.pos = "r",tl.cex = 0.8,number.cex = 0.8)
  corrplot(ciuprMAT,method="circle",type="lower",
           addCoef.col = "black",col=colorRampPalette(c("blue","white","orange2"))(25),tl.pos="n",na.label=" ",tl.col = "blue",cl.cex = 0.8,cl.pos = "r",tl.cex = 0.8,number.cex = 0.8)
  
  #get CI labels
  conf <- paste0("[", round(cilwrMAT, digits=2), ":", round(ciuprMAT, digits=2), "]")
  #mapping
  # xs <- lower.tri(row(cilwrMAT))
  xs <- row(cilwrMAT)
  ys <- (ncol(cilwrMAT)+1) - col(cilwrMAT)
  # ys <- (12) - col(res1[[1]])
  
  corrplot(rMAT,method="circle",type="full",
           addCoef.col = "black",col=colorRampPalette(c("blue","white","orange2"))(25),tl.pos="n",na.label=" ",tl.col = "blue",cl.cex = 0.8,cl.pos = "r",tl.cex = 0.8,number.cex = 0.8)
  conf[conf=="[NA:NA]"]=""
  #remove 1:1 which are only on the diagonal (area)
  if (brainvar == "area") {
    conf[conf=="[1:1]"]=""
  } else if (brainvar == "thickness") {
    conf[conf=="[1:1]"]=""
    conf[7]="[1:1]"
    conf[121]="[1:1]"
  }
  text(xs, ys, conf, pos=1, cex=0.5)
  #works for full matrix, there is one remaining NA because the upr CI estimate was not produced by OpenMx (area)
  
  if (saveres == 1) {
    filename = here(paste0("results/heritability/gencormattwin.",brainvar,".png"))
    print(paste("saving", filename))
    png(filename = filename, width = 15, height = 15, units = "cm", res = 300)
    if (brainvar == "area") {
      corrplot(rMAT,method="circle",type="full",
               addCoef.col = "black",col=colorRampPalette(c("blue","white","orange2"))(25),tl.pos="n",na.label=" ",tl.col = "blue",cl.cex = 0.8,cl.pos = "r",tl.cex = 0.8,number.cex = 0.8)
      text(xs, ys, conf, pos=1, cex=0.5)
    } else if (brainvar == "thickness") {
      corrplot(rMAT,method="circle",type="full",
               addCoef.col = "black",col=colorRampPalette(c("blue","white","orange2"))(25),tl.pos="n",na.label=" ",tl.col = "blue",cl.cex = 0.8,cl.pos = "r",tl.cex = 0.8,number.cex = 0.8)
      text(xs, ys, conf, pos=1, cex=0.5)
    }
    dev.off()
  }
  
  if (brainvar == "area") {
    #two results were FDR corrected sig
    print("two results FDR corrected")
    pp1[pp1$fdr<.05,]
    #pvals
    # corrplot(porigMATarea,method="circle",type="lower",
    #          addCoef.col = "black",col=colorRampPalette(c("blue","white","orange2"))(25),tl.pos="n",na.label=" ",tl.col = "blue",cl.cex = 0.8,cl.pos = "r",tl.cex = 0.8,number.cex = 0.8, is.corr = F)
  } else if (brainvar == "thickness") {
    print("no results FDR corrected")
    pp1[pp1$fdr<.05,]
  }
}

  
  