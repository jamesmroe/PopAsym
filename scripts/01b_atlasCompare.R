# install.packages("here")
library("here")
library("tidyverse")
here()


#SELECT OPTIONS ======================
saveres=0 #save plots 1/0
#SCRIPT TAKES SUMMARY-LEVEL DATA AS INPUT AND IS THUS FULLY REPRODUCIBLE HERE
#====================================#


cohorts=c("LCBC_55minus","UKB_sample","HCP")
ffs=c("area", "thickness")
pp=logs=mod1=ppc=ROUT=RRR=list()
cc=0
for (f in ffs) {
  for (c in cohorts){
    # c="LCBC_55minus"
    # f="area"
    cc=cc+1
    
    RRR=rbind(read.csv(here("results/atlasCompare",paste0(f,"-",c,"-","valsAparcCorrected.txt")),sep=" ",header=F),
              read.csv(here("results/atlasCompare",paste0(f,"-",c,"-","valsAsymLabelsCorrected.txt")),sep=" ",header=F))
    names(RRR)=c("roi","mmc","mm","ss","ci","nvtx","type")
    # plot(factor(RRR$type),RRR$mmc)
    ROUT[[cc]] = RRR
    
    
    logs[[cc]]=paste0(f,"_",c)
    mod1[[cc]]=summary(lm(scale(RRR$mm)~RRR$type + RRR$type*RRR$nvtx))  
    
    print(paste(f,c,"aparc =", round(mean(RRR$mm[RRR$type=="aparc"]),2),
                "±", round(sd(RRR$mm[RRR$type=="aparc"]),2)))
    
    print(paste(f,c,"asym =", round(mean(RRR$mm[RRR$type=="asym"]),2),
                "±", round(sd(RRR$mm[RRR$type=="asym"]),2)))
    
    #make long format
    RRR_long = RRR %>%
      pivot_longer(
        cols = starts_with("mm"),
        names_to = "m",
        values_to = "mm",
        values_drop_na = TRUE
      ) %>% mutate(desc=paste0(type,"_",m),
                   desc=factor(desc,ordered=T,
                               levels = c("aparc_mm","asym_mm","aparc_mmc","asym_mmc"))) #mm is raw means, mmc is corrected means
    ppc[[cc]] =
      ggplot(RRR_long,aes(desc,mm,col=type)) +
      geom_boxplot(aes(fill=type),alpha=0.3) +
      geom_jitter(width=0.1,size=3) +
      scale_fill_manual(values = c("#37454B", "#F2C500")) +
      scale_color_manual(values = c("#37454B", "#F2C500")) +
      theme_classic() +
      theme(legend.position = "top",
            text = element_text(color="black",size=18,family="Helvetica Neue Light"),
            axis.ticks.x=element_blank()) +
      labs(y="Vertex-mm R")
    
    if (saveres == 1) {
      ggsave(filename = here("results/atlasCompare",paste0("p-",f,"-",c,"-",".png")), plot = ppc[[cc]], width = 8, height=12, dpi=400, units="cm")
    }
  }
}

#area 
#LCBC-v-UKB
cor(ROUT[[1]]$mm,ROUT[[2]]$mm)
#LCBC-v-HCP
cor(ROUT[[1]]$mm,ROUT[[3]]$mm)
#HCP-v-UKB
cor(ROUT[[2]]$mm,ROUT[[3]]$mm)

#area 
#LCBC-v-UKB
cor(ROUT[[4]]$mm,ROUT[[5]]$mm)
#LCBC-v-HCP
cor(ROUT[[4]]$mm,ROUT[[6]]$mm)
#HCP-v-UKB
cor(ROUT[[5]]$mm,ROUT[[6]]$mm)

OUT = data.frame(Dataset=logs %>% unlist(),
                 Metric = c(rep("Area",3),rep("Thickness",3)),
                 B=rbind(round(coef(mod1[[1]])[2,1],2),
                         round(coef(mod1[[2]])[2,1],2),
                         round(coef(mod1[[3]])[2,1],2),
                         round(coef(mod1[[4]])[2,1],2),
                         round(coef(mod1[[5]])[2,1],2),
                         round(coef(mod1[[6]])[2,1],2)
                 ),
                 P=rbind(round(coef(mod1[[1]])[2,4],3),
                         round(coef(mod1[[2]])[2,4],3),
                         round(coef(mod1[[3]])[2,4],3),
                         round(coef(mod1[[4]])[2,4],3),
                         round(coef(mod1[[5]])[2,4],3),
                         round(coef(mod1[[6]])[2,4],3)
                 ),
                 B2=rbind(formatC(coef(mod1[[1]])[3,1],format="e",digits=2),
                          formatC(coef(mod1[[2]])[3,1],format="e",digits=2),
                          formatC(coef(mod1[[3]])[3,1],format="e",digits=2),
                          formatC(coef(mod1[[4]])[3,1],format="e",digits=2),
                          formatC(coef(mod1[[5]])[3,1],format="e",digits=2),
                          formatC(coef(mod1[[6]])[3,1],format="e",digits=2)
                 ),
                 P2=rbind(formatC(coef(mod1[[1]])[3,4],format="e",digits=2),
                          formatC(coef(mod1[[2]])[3,4],format="e",digits=2),
                          formatC(coef(mod1[[3]])[3,4],format="e",digits=2),
                          formatC(coef(mod1[[4]])[3,4],format="e",digits=2),
                          formatC(coef(mod1[[5]])[3,4],format="e",digits=2),
                          formatC(coef(mod1[[6]])[3,4],format="e",digits=2)
                 ),
                 B3=rbind(formatC(coef(mod1[[1]])[4,1],format="e",digits=2),
                          formatC(coef(mod1[[2]])[4,1],format="e",digits=2),
                          formatC(coef(mod1[[3]])[4,1],format="e",digits=2),
                          formatC(coef(mod1[[4]])[4,1],format="e",digits=2),
                          formatC(coef(mod1[[5]])[4,1],format="e",digits=2),
                          formatC(coef(mod1[[6]])[4,1],format="e",digits=2)
                 ),
                 P2=rbind(round(coef(mod1[[1]])[4,4],2),
                          round(coef(mod1[[2]])[4,4],2),
                          round(coef(mod1[[3]])[4,4],2),
                          round(coef(mod1[[4]])[4,4],2),
                          round(coef(mod1[[5]])[4,4],2),
                          round(coef(mod1[[6]])[4,4],2)),
                 r2=rbind(round(mod1[[1]]$adj.r.squared,2),
                          round(mod1[[2]]$adj.r.squared,2),
                          round(mod1[[3]]$adj.r.squared,2),
                          round(mod1[[4]]$adj.r.squared,2),
                          round(mod1[[5]]$adj.r.squared,2),
                          round(mod1[[6]]$adj.r.squared,2)),
                 df=rbind(mod1[[1]]$df[2],
                          mod1[[2]]$df[2],
                          mod1[[3]]$df[2],
                          mod1[[4]]$df[2],
                          mod1[[5]]$df[2],
                          mod1[[6]]$df[2])
                 
)
OUT
