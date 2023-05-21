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

# brainvar="area" #select metric
brainvar="thickness"
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


for (j in 1:ll) {
  print(j)
  
  roi1 = roinames[j]
  
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
  
  if (brainvar == "area") {
    jj = j
  } else if (brainvar == "thickness") {
    jj = j+14
  }
  head(dat)
  
  #inverse if rightward
  if (substr(roi1,1,1) == "R") {
    print(paste("inversing rightward ROI trait1. Mean =", mean(dat$lab)))
    dat$lab = dat$lab * -1
  }
  
  
  #correct for age and sex
  print("correcting for age and sex")
  co=coef(lm(lab ~ Age + Sex, dat))
  dat %<>% mutate(Ycor = lab - (Age*co[2]) )
  dat %<>% mutate(Ycor = ifelse(Sex=="M",Ycor-co[3],Ycor),
                  trait = sprintf("%02d", jj),
                  roi = roinames[j])
  if (saveres == 1) {
    filename = paste0("/Users/jamesroe/LCBC/Users/jamesroe/PHD_project/Paper3/data/correctedLabs4Gencor/","invresidtwindatcor_",sprintf("%02d", jj),"_",roinames[j],"_",brainvar,".csv")
    print(paste("saving",filename))
    write.table(dat,filename,quote=F,row.names = F,col.names = T)
  }
}
