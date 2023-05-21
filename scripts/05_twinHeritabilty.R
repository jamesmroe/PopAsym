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
  # brainvar="global" #cluster-wise only
  #====================================#
  
  
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
      labs=names(DD)[grepl("label",names(DD))]
      ll = length(labs)
    }
  } else {
    labs=c(list.files(file.path(datdir,paste0("split.data.schaefer.",brainvar)),pattern="AI"))
    ll = length(labs)
  }
  
  
  print(paste("running ACE for", brainvar, ", schaefer =", schaefer),)
  for (j in 1:ll) {
    if (j==1) {
      sig = list(list())
      #original ACE model
      a2 = c2 = e2 = V = pA = pC = LL_ACE = LL_CE = LL_AE = 0
      a2cilwr = a2ciest = a2ciupr = 0
      pp = list()
      #post review AE model
      a2revise = e2revise = Vrevise = a2revisecilwr = a2reviseciest = a2reviseciupr = pArevise = LL_E = 0
      rMZ = rDZsibPairs = 0
    }
    print(j)
    
    
    #if cortex-wide
    if (schaefer == 1) {
      lab=labs[j]
      print(lab)
      D=read.csv(file.path(datdir,paste0("split.data.schaefer.",brainvar),lab),sep=",",stringsAsFactors = F,header=F)
      D=as.data.frame(t(D))
      D=data.frame(D,DF$Subject)
    } else {
      if (brainvar != "global") {
        #if clusterwise
        lab=labs[j]
        D=data.frame(DD[[lab]],DD$Subject)
      } else {
        #if wholehemimeas
        D=data.frame(DD[j],DF$Subject)
      } 
    }
    
    
    #EXTRACT FAM DATA FROM FULL HCP
    #reorders by family structure
    dat=data.frame()
    c=0
    print("extracting family data from mainframe")
    for (ii in I) {
      c=c+1
      dat[c,1] = D[which(DF$Subject==ii),2]
      dat[c,2] = D[which(DF$Subject==ii),1]
      dat[c,3] = DF[which(DF$Subject==ii),"Age_in_Yrs"] #add age
      dat[c,4] = DF[which(DF$Subject==ii),"Sex"] #add sex
    }
    names(dat)=c("ID","lab","Age","Sex")
    
    
    #correct for age and sex
    print("correcting for age and sex")
    co=coef(lm(lab ~ Age + Sex, dat))
    dat %<>% mutate(Ycor = lab - (Age*co[2]) )
    dat %<>% mutate(Ycor = ifelse(Sex=="M",Ycor-co[3],Ycor))
    
    
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
    mtrX <- populate_matrix(dat$ID, dat$Ycor, as.matrix(mtr))
  
    #check worked
    # head(mtr)
    # dat[which(dat$ID==139637),]
    # dat[which(dat$ID==196750),]
    # head(mtrX)
    # tail(mtr)
    # dat[which(dat$ID==959574),]
    # tail(mtrX)
    
    
    DataMZplus=mtrX[1:lenMZ,]
    DataDZplus=mtrX[(lenMZ+1):nrow(mtrX),]
    mtrMZ=mtr[1:lenMZ,]
    mtrDZ=mtr[(lenMZ+1):nrow(mtr),]
    
    
    selVars <- c('v1','v2','v3','v4','v5','v6')
    dimnames(DataMZplus) <- list(NULL,selVars)
    dimnames(DataDZplus) <- list(NULL,selVars)
    
    
    #EXCLUDE 2 HCP OUTLIERS
    #identify subject exclusions in mtr and delete from mtrx
    mtrX=mtrX[-which(mtr==132118, arr.ind = TRUE)[1],]
    mtrX=mtrX[-which(mtr==894067, arr.ind = TRUE)[1],]
    DataMZplus=DataMZplus[-which(mtrMZ==132118, arr.ind = TRUE)[1],]
    DataDZplus=DataDZplus[-which(mtrDZ==894067, arr.ind = TRUE)[1],]
    # length(mtrX[which(!is.na(mtrX), arr.ind = TRUE)])
    # dim(mtrX)
     
    
    
    #RESCALE DUE TO LOW VARIANCE
    DataMZplus_sc = scale(DataMZplus,
                          rep(mean(unlist(DataMZplus), na.rm = T), ncol(DataMZplus)),
                          rep(sd(unlist(DataMZplus), na.rm = T), ncol(DataMZplus)))
    
    DataDZplus_sc = scale(DataDZplus,
                          rep(mean(unlist(DataDZplus), na.rm = T), ncol(DataDZplus)),
                          rep(sd(unlist(DataDZplus), na.rm = T), ncol(DataDZplus)))
    
    if (j==1) {
      print("expected covariance matrix MZ")
      print(rbind(cbind("A+C+E", "A+C", ".5%x%A+C", ".5%x%A+C", ".5%x%A+C", ".5%x%A+C"),
            cbind("A+C", "A+C+E", ".5%x%A+C", ".5%x%A+C", ".5%x%A+C", ".5%x%A+C"),
            cbind(".5%x%A+C", ".5%x%A+C", "A+C+E", ".5%x%A+C", ".5%x%A+C", ".5%x%A+C"),
            cbind(".5%x%A+C", ".5%x%A+C", ".5%x%A+C", "A+C+E", ".5%x%A+C", ".5%x%A+C"),
            cbind(".5%x%A+C", ".5%x%A+C", ".5%x%A+C", ".5%x%A+C", "A+C+E", ".5%x%A+C"),
            cbind(".5%x%A+C", ".5%x%A+C", ".5%x%A+C", ".5%x%A+C", ".5%x%A+C", "A+C+E")))
      
      print("expected covariance matrix DZ")
      print(rbind(cbind("A+C+E", ".5%x%A+C", ".5%x%A+C", ".5%x%A+C", ".5%x%A+C", ".5%x%A+C"),
            cbind(".5%x%A+C", "A+C+E", ".5%x%A+C", ".5%x%A+C", ".5%x%A+C", ".5%x%A+C"),
            cbind(".5%x%A+C", ".5%x%A+C", "A+C+E", ".5%x%A+C", ".5%x%A+C", ".5%x%A+C"),
            cbind(".5%x%A+C", ".5%x%A+C", ".5%x%A+C", "A+C+E", ".5%x%A+C", ".5%x%A+C"),
            cbind(".5%x%A+C", ".5%x%A+C", ".5%x%A+C", ".5%x%A+C", "A+C+E", ".5%x%A+C"),
            cbind(".5%x%A+C", ".5%x%A+C", ".5%x%A+C", ".5%x%A+C", ".5%x%A+C", "A+C+E")))
    }
    
    library("OpenMx")
    ff = mxFitFunctionML()
    print("running ACE")
    twinACEModel <- mxModel("twinACE", 
                            mxMatrix("Full", nrow=1, ncol=6, free=TRUE, values=0, label="mean", name="expMean"),
                            # Matrix expMean for expected mean 
                            # vector for MZ and DZ twins    
                            # -------------------------------------
                            
                            mxMatrix("Lower", nrow=1, ncol=1, free=TRUE, values=.6, label="a", name="X"), #all parameters specified
                            mxMatrix("Lower", nrow=1, ncol=1, free=TRUE, values=.6, label="c", name="Y"), #all parameters specified
                            mxMatrix("Lower", nrow=1, ncol=1, free=TRUE, values=.6, label="e", name="Z"), #all parameters specified
                            # Matrices X, Y, and Z to store the 
                            # a, c, and e path coefficients
                            # -------------------------------------
                            
                            mxAlgebra(X * t(X), name="A"),
                            mxAlgebra(Y * t(Y), name="C"),
                            mxAlgebra(Z * t(Z), name="E"),	
                            # Matrixes A, C, and E to compute 
                            # A, C, and E variance components
                            # -------------------------------------
                            
                            
                            mxAlgebra(rbind(cbind(A+C+E, A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C),
                                            cbind(A+C, A+C+E, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C),
                                            cbind(.5%x%A+C, .5%x%A+C, A+C+E, .5%x%A+C, .5%x%A+C, .5%x%A+C),
                                            cbind(.5%x%A+C, .5%x%A+C, .5%x%A+C, A+C+E, .5%x%A+C, .5%x%A+C),
                                            cbind(.5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, A+C+E, .5%x%A+C),
                                            cbind(.5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, A+C+E)), 
                                      name="expCovMZ"),
                            
                            # Matrix expCOVMZ for expected 
                            # covariance matrix for MZ twins
                            # -------------------------------------
                            
                            mxAlgebra(rbind(cbind(A+C+E, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C),
                                            cbind(.5%x%A+C, A+C+E, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C),
                                            cbind(.5%x%A+C, .5%x%A+C, A+C+E, .5%x%A+C, .5%x%A+C, .5%x%A+C),
                                            cbind(.5%x%A+C, .5%x%A+C, .5%x%A+C, A+C+E, .5%x%A+C, .5%x%A+C),
                                            cbind(.5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, A+C+E, .5%x%A+C),
                                            cbind(.5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, A+C+E)), 
                                      name="expCovDZ"),
                            # Matrix expCOVMZ for expected 
                            # covariance matrix for DZ twins
                            # -------------------------------------
                            
                            # MZ Model
                            # -------------------------------------
                            mxModel("MZ",
                                    mxData(DataMZplus_sc, type="raw"),
                                    mxExpectationNormal("twinACE.expCovMZ", "twinACE.expMean",selVars),
                                    mxFitFunctionML()),
                            
                            # DZ Model
                            # -------------------------------------
                            mxModel("DZ", 
                                    mxData(DataDZplus_sc, type="raw"),
                                    mxExpectationNormal("twinACE.expCovDZ", "twinACE.expMean",selVars),
                                    mxFitFunctionML()),
                            
                            
                            # # Just easier
                            mxFitFunctionMultigroup(c("MZ", "DZ")),
                            # mxAlgebra(MZ.objective + DZ.objective, name="twin"),
                            # mxFitFunctionAlgebra("twin")
                            
                            mxAlgebra(a^2 / (a^2 + c^2 + e^2), name = "h2"),
                            mxCI(c("a", "c", "e", "h2"))
                            
    )
    twinACEFit <- mxRun(twinACEModel, suppressWarnings=FALSE, intervals = TRUE)
    summary(twinACEFit)
    SS=summary(twinACEFit)
    
    # Specify and Run ACE Model with RawData and Matrix-style Input
    # -----------------------------------------------------------------------
    LL_ACE[j] <- mxEval(objective, twinACEFit)
    MZc <- mxEval(expCovMZ, twinACEFit)
    DZc <- mxEval(expCovDZ, twinACEFit)
    M   <- mxEval(expMean, twinACEFit)
    # Retrieve expected mean vector and 
    # expected covariance matrices
    # -------------------------------------
    
    A <- mxEval(A, twinACEFit)
    C <- mxEval(C, twinACEFit)
    E <- mxEval(E, twinACEFit)
    # Retrieve the A, C, and E 
    # variance components
    # -------------------------------------
    
    V[j] <- (A+C+E)
    a2[j] <- A/V[j]
    mxEval(h2, twinACEFit)
    c2[j] <- C/V[j]
    e2[j] <- E/V[j]
    
    a2cilwr[j]=SS$CI[4,1]
    a2ciest[j]=SS$CI[4,2]
    a2ciupr[j]=SS$CI[4,3]
    
    # Fit CE model and compare for sig
    # -----------------------------------------------------------------------
    twinCEModel <- mxModel("twinCE", 
                           mxMatrix("Full", nrow=1, ncol=6, free=TRUE, values=0, label="mean", name="expMean"),
                           # Matrix expMean for expected mean 
                           # vector for MZ and DZ twins    
                           # -------------------------------------
                           
                           mxMatrix("Lower", nrow=1, ncol=1, free=FALSE, values=0, label="a", name="X"), #drop genetic parameter 
                           mxMatrix("Lower", nrow=1, ncol=1, free=TRUE, values=.6, label="c", name="Y"),
                           mxMatrix("Lower", nrow=1, ncol=1, free=TRUE, values=.6, label="e", name="Z"),
                           # Matrices X, Y, and Z to store the 
                           # a, c, and e path coefficients
                           # -------------------------------------
                           
                           mxAlgebra(X * t(X), name="A"),
                           mxAlgebra(Y * t(Y), name="C"),
                           mxAlgebra(Z * t(Z), name="E"),	
                           # Matrixes A, C, and E to compute 
                           # A, C, and E variance components
                           # -------------------------------------
                           
                           mxAlgebra(rbind(cbind(A+C+E, A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C),
                                           cbind(A+C, A+C+E, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C),
                                           cbind(.5%x%A+C, .5%x%A+C, A+C+E, .5%x%A+C, .5%x%A+C, .5%x%A+C),
                                           cbind(.5%x%A+C, .5%x%A+C, .5%x%A+C, A+C+E, .5%x%A+C, .5%x%A+C),
                                           cbind(.5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, A+C+E, .5%x%A+C),
                                           cbind(.5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, A+C+E)), 
                                     name="expCovMZ"),
                           
                           # Matrix expCOVMZ for expected 
                           # covariance matrix for MZ twins
                           # -------------------------------------
                           
                           mxAlgebra(rbind(cbind(A+C+E, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C),
                                           cbind(.5%x%A+C, A+C+E, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C),
                                           cbind(.5%x%A+C, .5%x%A+C, A+C+E, .5%x%A+C, .5%x%A+C, .5%x%A+C),
                                           cbind(.5%x%A+C, .5%x%A+C, .5%x%A+C, A+C+E, .5%x%A+C, .5%x%A+C),
                                           cbind(.5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, A+C+E, .5%x%A+C),
                                           cbind(.5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, A+C+E)), 
                                     name="expCovDZ"),
                           # Matrix expCOVMZ for expected 
                           # covariance matrix for DZ twins
                           # -------------------------------------
                           
                           # MZ Model
                           # -------------------------------------
                           mxModel("MZ",
                                   mxData(DataMZplus_sc, type="raw"),
                                   mxExpectationNormal("twinCE.expCovMZ", "twinCE.expMean",selVars),
                                   mxFitFunctionML()),
                           
                           # DZ Model
                           # -------------------------------------
                           mxModel("DZ", 
                                   mxData(DataDZplus_sc, type="raw"),
                                   mxExpectationNormal("twinCE.expCovDZ", "twinCE.expMean",selVars),
                                   mxFitFunctionML()),
                           
                           
                           # # Just easier
                           mxFitFunctionMultigroup(c("MZ", "DZ"))
                           # mxAlgebra(MZ.objective + DZ.objective, name="twin"),
                           # mxFitFunctionAlgebra("twin")
    )
    twinCEFit <- mxRun(twinCEModel, suppressWarnings=FALSE, intervals=TRUE)
    # LL_CE[j] <- mxEval(objective, twinCEFit)
    
    sig[[j]]=list()
    sig[[j]]$CEcomp = mxCompare(twinACEFit,twinCEFit)
    sig[[j]]$CEsig = mxCompare(twinACEFit,twinCEFit)$p[2]
    pA[j]=mxCompare(twinACEFit,twinCEFit)$p[2]
    LL_ACE[j]=mxCompare(twinACEFit,twinCEFit)$minus2LL[1]
    LL_CE[j]=mxCompare(twinACEFit,twinCEFit)$minus2LL[2]
    
    
    # Fit AE model and compare for sig
    # -----------------------------------------------------------------------
    twinAEModel <- mxModel("twinAE", 
                           mxMatrix("Full", nrow=1, ncol=6, free=TRUE, values=0, label="mean", name="expMean"),
                           # Matrix expMean for expected mean 
                           # vector for MZ and DZ twins    
                           # -------------------------------------
                           
                           mxMatrix("Lower", nrow=1, ncol=1, free=TRUE, values=.6, label="a", name="X"),
                           mxMatrix("Lower", nrow=1, ncol=1, free=FALSE, values=0, label="c", name="Y"), #drop shared environment parameter 
                           mxMatrix("Lower", nrow=1, ncol=1, free=TRUE, values=.6, label="e", name="Z"),
                           # Matrices X, Y, and Z to store the 
                           # a, c, and e path coefficients
                           # -------------------------------------
                           
                           mxAlgebra(X * t(X), name="A"),
                           mxAlgebra(Y * t(Y), name="C"),
                           mxAlgebra(Z * t(Z), name="E"),	
                           # Matrixes A, C, and E to compute 
                           # A, C, and E variance components
                           # -------------------------------------
                           
                           mxAlgebra(rbind(cbind(A+C+E, A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C),
                                           cbind(A+C, A+C+E, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C),
                                           cbind(.5%x%A+C, .5%x%A+C, A+C+E, .5%x%A+C, .5%x%A+C, .5%x%A+C),
                                           cbind(.5%x%A+C, .5%x%A+C, .5%x%A+C, A+C+E, .5%x%A+C, .5%x%A+C),
                                           cbind(.5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, A+C+E, .5%x%A+C),
                                           cbind(.5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, A+C+E)), 
                                     name="expCovMZ"),
                           
                           # Matrix expCOVMZ for expected 
                           # covariance matrix for MZ twins
                           # -------------------------------------
                           
                           mxAlgebra(rbind(cbind(A+C+E, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C),
                                           cbind(.5%x%A+C, A+C+E, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C),
                                           cbind(.5%x%A+C, .5%x%A+C, A+C+E, .5%x%A+C, .5%x%A+C, .5%x%A+C),
                                           cbind(.5%x%A+C, .5%x%A+C, .5%x%A+C, A+C+E, .5%x%A+C, .5%x%A+C),
                                           cbind(.5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, A+C+E, .5%x%A+C),
                                           cbind(.5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, A+C+E)), 
                                     name="expCovDZ"),
                           # Matrix expCOVMZ for expected 
                           # covariance matrix for DZ twins
                           # -------------------------------------
                           
                           # MZ Model
                           # -------------------------------------
                           mxModel("MZ",
                                   mxData(DataMZplus_sc, type="raw"),
                                   mxExpectationNormal("twinAE.expCovMZ", "twinAE.expMean",selVars),
                                   mxFitFunctionML()),
                           
                           # DZ Model
                           # -------------------------------------
                           mxModel("DZ", 
                                   mxData(DataDZplus_sc, type="raw"),
                                   mxExpectationNormal("twinAE.expCovDZ", "twinAE.expMean",selVars),
                                   mxFitFunctionML()),
                           
                           
                           # # Just easier
                           mxFitFunctionMultigroup(c("MZ", "DZ")),
                           # mxAlgebra(MZ.objective + DZ.objective, name="twin"),
                           # mxFitFunctionAlgebra("twin")
                           
                           mxAlgebra(a^2 / (a^2 + c^2 + e^2), name = "h2"),
                           mxCI(c("a", "c", "e", "h2"))
    )
    twinAEFit <- mxRun(twinAEModel, suppressWarnings=FALSE, intervals=TRUE)
    # LL_AE[j] <- mxEval(objective, twinAEFit)
    SS2 = summary(twinAEFit)
    
    sig[[j]]$AEcomp = mxCompare(twinACEFit,twinAEFit)
    sig[[j]]$AEsig = mxCompare(twinACEFit,twinAEFit)$p[2]
    pC[j]=mxCompare(twinACEFit,twinAEFit)$p[2]
    LL_AE[j]=mxCompare(twinACEFit,twinAEFit)$minus2LL[2]
    
    
    MZcAE <- mxEval(expCovMZ, twinAEFit)
    DZcAE <- mxEval(expCovDZ, twinAEFit)
    MAE   <- mxEval(expMean, twinAEFit)
    
    Arevise <- mxEval(A, twinAEFit)
    Crevise <- mxEval(C, twinAEFit)
    Erevise <- mxEval(E, twinAEFit)
    
    Vrevise[j] <- (Arevise+Crevise+Erevise)
    a2revise[j] <- Arevise/Vrevise[j]
    mxEval(h2, twinAEFit)
    e2revise[j] <- Erevise/Vrevise[j]
    
    a2revisecilwr[j]=SS2$CI[4,1]
    a2reviseciest[j]=SS2$CI[4,2]
    a2reviseciupr[j]=SS2$CI[4,3]
    
    
    ##REVISION##
    # compare to fiitted AE model for sig
    # -----------------------------------------------------------------------
    twinEModel <- mxModel("twinE", 
                           mxMatrix("Full", nrow=1, ncol=6, free=TRUE, values=0, label="mean", name="expMean"),
                           # Matrix expMean for expected mean 
                           # vector for MZ and DZ twins    
                           # -------------------------------------
                           
                           mxMatrix("Lower", nrow=1, ncol=1, free=FALSE, values=0, label="a", name="X"), #drop A
                           mxMatrix("Lower", nrow=1, ncol=1, free=FALSE, values=0, label="c", name="Y"), #drop shared environment parameter 
                           mxMatrix("Lower", nrow=1, ncol=1, free=TRUE, values=.6, label="e", name="Z"),
                           # Matrices X, Y, and Z to store the 
                           # a, c, and e path coefficients
                           # -------------------------------------
                           
                           mxAlgebra(X * t(X), name="A"),
                           mxAlgebra(Y * t(Y), name="C"),
                           mxAlgebra(Z * t(Z), name="E"),	
                           # Matrixes A, C, and E to compute 
                           # A, C, and E variance components
                           # -------------------------------------
                           
                           mxAlgebra(rbind(cbind(A+C+E, A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C),
                                           cbind(A+C, A+C+E, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C),
                                           cbind(.5%x%A+C, .5%x%A+C, A+C+E, .5%x%A+C, .5%x%A+C, .5%x%A+C),
                                           cbind(.5%x%A+C, .5%x%A+C, .5%x%A+C, A+C+E, .5%x%A+C, .5%x%A+C),
                                           cbind(.5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, A+C+E, .5%x%A+C),
                                           cbind(.5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, A+C+E)), 
                                     name="expCovMZ"),
                           
                           # Matrix expCOVMZ for expected 
                           # covariance matrix for MZ twins
                           # -------------------------------------
                           
                           mxAlgebra(rbind(cbind(A+C+E, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C),
                                           cbind(.5%x%A+C, A+C+E, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C),
                                           cbind(.5%x%A+C, .5%x%A+C, A+C+E, .5%x%A+C, .5%x%A+C, .5%x%A+C),
                                           cbind(.5%x%A+C, .5%x%A+C, .5%x%A+C, A+C+E, .5%x%A+C, .5%x%A+C),
                                           cbind(.5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, A+C+E, .5%x%A+C),
                                           cbind(.5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, .5%x%A+C, A+C+E)), 
                                     name="expCovDZ"),
                           # Matrix expCOVMZ for expected 
                           # covariance matrix for DZ twins
                           # -------------------------------------
                           
                           # MZ Model
                           # -------------------------------------
                           mxModel("MZ",
                                   mxData(DataMZplus_sc, type="raw"),
                                   mxExpectationNormal("twinE.expCovMZ", "twinE.expMean",selVars),
                                   mxFitFunctionML()),
                           
                           # DZ Model
                           # -------------------------------------
                           mxModel("DZ", 
                                   mxData(DataDZplus_sc, type="raw"),
                                   mxExpectationNormal("twinE.expCovDZ", "twinE.expMean",selVars),
                                   mxFitFunctionML()),
                           
                           
                           # # Just easier
                           mxFitFunctionMultigroup(c("MZ", "DZ"))
                           # mxAlgebra(MZ.objective + DZ.objective, name="twin"),
                           # mxFitFunctionAlgebra("twin")
    )
    twinEFit <- mxRun(twinEModel, suppressWarnings=FALSE, intervals=TRUE)
    LL_E[j] <- mxEval(objective, twinEFit)
    
    
    sig[[j]]$Ecomp = mxCompare(twinAEFit,twinEFit)
    sig[[j]]$Esig = mxCompare(twinAEFit,twinEFit)$p[2]
    pArevise[j]=mxCompare(twinAEFit,twinEFit)$p[2]
    LL_E[j]=mxCompare(twinAEFit,twinEFit)$minus2LL[2]
    
  
    # estimate correlations in MZ DZ ----
    # View(DataMZplus_sc)
    # View(DataDZplus_sc)
    # take entire first 2 cols for MZ correlations (i.e. all MZ obs)
    rMZ[j] = cor(DataMZplus_sc[,1],DataMZplus_sc[,2])
    
    # for DZ/sibs need to take maximum N across two cols possible
    
    # first instance of a sib pair in MZ data, one of which is an MZ with another twin not taken here
    firstSibinMZdat = min(which(!is.na(DataMZplus_sc[,3])))
    TMP1 = DataMZplus_sc[firstSibinMZdat:nrow(DataMZplus_sc),2:3]
    
    # from first instance of a sibpair in DZ data, shift the data into the first two cols
    firstSibsOnlyObsinDZdat = min(which(is.na(DataDZplus_sc[,1])))
    
    # now have all complete MZ pairs
    #+ all complete DZ pairs
    #+ all complete sibpairs
    # correlation for "DZ" computed across the latter two
    
    TMPDZ = DataDZplus_sc
    TMPDZ[firstSibsOnlyObsinDZdat:nrow(TMPDZ),1:2] =  TMPDZ[firstSibsOnlyObsinDZdat:nrow(TMPDZ),2:3]
    #NB! data now duplicate in cols 2 and 3 of TMPDZ, wipe other cols
    TMPDZ=TMPDZ[,-c(3:6)]
    
    TMP2 = rbind(TMP1,TMPDZ)
    
    rDZsibPairs[j] = cor(TMP2[,1],TMP2[,2])
    
  }
  
  
  if (brainvar == "global") {
    rois = names(DD)[1:ll]
  } else {
    rois = labs
  }
  
  if (schaefer != 1) {
    aceout=data.frame(roi=rois,
                      a2=round(a2,3),
                      a2cilwr=round(a2cilwr,5),
                      a2ciupr=round(a2ciupr,3),
                      c2=round(c2,3),
                      e2=round(e2,3),
                      pA=pA,
                      pC=pC,
                      LL_ACE=round(LL_ACE,1),
                      LL_CE=round(LL_CE,1),
                      LL_AE=round(LL_AE,1),
                      LL_E=round(LL_E,1),
                      
                      a2revise=round(a2revise,3),
                      a2revisecilwr=round(a2revisecilwr,5),
                      a2reviseciupr=round(a2reviseciupr,3),
                      e2revise=round(e2revise,3),
                      pArevise=pArevise,
                      rMZ = round(rMZ,3), 
                      rDZsibPairs = round(rDZsibPairs,3))
  } else {
    aceout=data.frame(#roi=labs,
                      a2=round(a2,3),
                      a2cilwr=(a2cilwr),
                      a2ciupr=(a2ciupr),
                      c2=round(c2,3),
                      e2=round(e2,3),
                      pA=pA,
                      pC=pC,
                      LL_ACE=round(LL_ACE,1),
                      LL_CE=round(LL_CE,1),
                      LL_AE=round(LL_AE,1),
                      LL_E=round(LL_E,1),
                      
                      a2revise=round(a2revise,3),
                      a2revisecilwr=round(a2revisecilwr,5),
                      a2reviseciupr=round(a2reviseciupr,3),
                      e2revise=round(e2revise,3),
                      pArevise=pArevise,
                      rMZ = round(rMZ,3), 
                      rDZsibPairs = round(rDZsibPairs,3))
                      
  }
  
  if (schaefer != 1) {
    if (brainvar == "global") {
      roinames = rois
    } else if (brainvar =="area") {
      roinames = c(
      "L_Postcentral-gyrus_Supramarginal",
      "L_Anterior-temporal_Parahippocampal-gyrus",
      "L_Anterior-cingulate_Subcallosal",
      "L_Anterior-insula",
      "L_Retrosplenial-cortex",
      "L_Temporal-pole_Inferior-temporal-gyrus",
      "L_Superior-frontal-gyrus",
      "R_Parieto-occipital_sulcus",
      "R_Inferior-parietal_Lateral-occipital",
      "R_Cingulate",
      "R_Middle-frontal-gyrus",
      "R_Superior-temporal-sulcus",
      "R_Superior-frontal-gyrus",
      "R_Gyrus-rectus")
    } else if (brainvar == "thickness") {
      roinames=c(
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
         "R_Superior-temporal-sulcus",
         "R_Lateral-occipital",
         "R_Lingual-gyrus",
         "R_Posterior-insula_Sylvian",
         "R_Entorhinal",
         "R_Superior-insula",
         "R_Planum-temporale",
         "R_Posterior-cingulate",
         "R_Anterior-insula")
    }
  } else {
    roinames = labs   
  }
  
  
  if (schaefer != 1) {
    #remove mergelabs & reorder labs to left first
    if (brainvar == "area") {
      aceout = aceout[c(9:15,1:7),]
      direct = c(rep("Leftward",7),rep("Rightward",7))
    } else if (brainvar == "thickness") {
      aceout = aceout[c(10:20,1:9),]
      direct = c(rep("Leftward",11),rep("Rightward",9))
    } else if (brainvar == "global") {
      # aceout = aceout[c(1:2,5,7,3:4,6,8),]
      direct = rep(0, nrow(aceout))
    }
  } else {
    direct = rep(0, nrow(aceout))
  }
  
  #Fix NAs on lower CIs
  if (length(is.na(aceout$a2cilwr))) {
    aceout$a2cilwr[is.na(aceout$a2cilwr)] = 0
    aceout$a2revisecilwr[is.na(aceout$a2revisecilwr)] = 0
  }
  
  (TABLE = data.frame(direct, 
                     ROI = seq(1,nrow(aceout)),
                     Cluster = roinames, 
                     
                     #original
                     a2 = paste0(aceout$a2,"_(",aceout$a2cilwr,"_-_",aceout$a2ciupr,")"),
                     pA = aceout$pA,
                     c2 = aceout$c2,
                     e2 = aceout$e2,
                     LLACE = aceout$LL_ACE,
                     LLCE = aceout$LL_CE,
                     
                     #revision
                     a2revise = paste0(aceout$a2revise,"_(",round(aceout$a2revisecilwr,3),"_-_",round(aceout$a2reviseciupr,3),")"),
                     pArevise = aceout$pArevise,
                     e2revise = aceout$e2,
                     LLAE = aceout$LL_AE,
                     LLE = aceout$LL_E,
                     rMZ = aceout$rMZ, 
                     rDZsibPairs = aceout$rDZsibPairs)
                     )
  
  # ggplot(aceout) +
  #   geom_errorbar(aes(y=factor(roi),x=a2,xmin=a2cilwr,xmax=a2ciupr,col=factor(roi)),width=0.1) +
  #   geom_point(aes(y=factor(roi),x=a2)) +
  #   geom_vline(xintercept = 0, linetype = 2)
  # ggplot(aceout) +
  #   geom_errorbar(aes(y=factor(roi),x=a2revise,xmin=a2revisecilwr,xmax=a2reviseciupr,col=factor(roi)),width=0.1) +
  #   geom_point(aes(y=factor(roi),x=a2revise)) +
  #   geom_vline(xintercept = 0, linetype = 2)
  
  if (!dir.exists(here("results/heritability"))) {
    dir.create(here("results/heritability"))
  }
  if (saveres == 1) {
    if (schaefer != 1) {
      print("saving output")
      write.table(TABLE, here("results/heritability",paste0("clustersTwinAE",brainvar,"_revision.csv")),quote=F,row.names = F,col.names = T)
    } else {
      print("saving output")
      write.table(aceout, here("results/heritability",paste0("cortexTwinAE",brainvar,"_revision.csv")),quote=F,row.names = F,col.names = T)
    }
  }

  