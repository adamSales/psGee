library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(sandwich)
library(lmtest)
library(forcats)
library(arm)
library(randomForest)
library(missForest)
library(estimatr)
library(splines)
library(rstan)
library(kableExtra)
library(texreg)
library(tableone)
library(xtable)
library(purrr)
library(tibble)
library(coefplot)
library(splines)
library(tikzDevice)
select <- dplyr::select




### load in data from the website (accessed 7/3/2025)
print(load(url("https://causeweb.org/tshs/datasets/OPT_Study_PersonLevel_Data.RData")))


############################################################
###
### process the data
###
#########################################################

###
opt <- opt%>%
  mutate(OFIBRIN1=as.numeric(as.character(OFIBRIN1)),
         ETXU_CAT1=as.numeric(as.character(ETXU_CAT1)))

levels(opt$Group)=c(T="Treatment",C="Control")[levels(opt$Group)]

cca <- opt%>%filter(!is.na(OFIBRIN1),!is.na(ETXU_CAT1),!is.na(V5..BOP))%>%
    #select(Group, OFIBRIN1,ETXU_CAT1,V5..BOP,Tx.comp.,Completed.EDC,BL..BOP)%>%
    mutate(S=Tx.comp.=='Yes',
         Y=V5..BOP/100,
         Z=Group=='Treatment')


#### table from the appendix
tab1 <- map(list(`Full Data`=opt,`Complete Cases`=cca),
            function(dat)
                dat%>%
                mutate(Data=if(any(is.na(V5..BOP))) "Full Data" else "Complete Cases")%>%
    #rename("Trt. Completed"=Tx.comp.)%>%
  group_by(Data,Group)%>%
  summarize(n=n(),across(c(OFIBRIN1,ETXU_CAT1,V5..BOP),~paste0(round(mean(.,na.rm=TRUE),1)," (",round(sd(.,na.rm=TRUE),1),")")),
    across(Tx.comp.,#`Trt. Completed`,
           list(
               No =~if(Group[1]=="Control") "-" else paste0(sum(.=="No ",na.rm=TRUE)," (",round(mean(.=="No ",na.rm=TRUE)*100,1),"%)"),
               Und=~if(Group[1]=="Control") "-" else paste0(sum(.=="Und",na.rm=TRUE)," (",round(mean(.=="Und",na.rm=TRUE)*100,1),"%)"),
               Yes=~if(Group[1]=="Control") "-" else paste0(sum(.=="Yes",na.rm=TRUE)," (",round(mean(.=="Yes",na.rm=TRUE)*100,1),"%)")),
           .names="Trt. Completed: {.fn}"))%>%
    select(everything(),V5..BOP)%>%
    rename(Endotoxin=ETXU_CAT1,Fibrinogen=OFIBRIN1,"% Sites Bleeding"=V5..BOP)%>%
    t())%>%
    do.call("cbind",.)

sink("results/optsTable1.tex")
kbl(tab1[-c(1:2),],format="latex",booktabs=TRUE,col.names=tab1[2,],caption="Descriptive statistics---mean and standard deviation or count and percent---for study variables in the full OPT dataset and in the analysis sample (i.e. complete cases)",label="optTab1" )%>%
    add_header_above(c(" ","Full Data"=2,"Complete Cases"=2))
sink()

##### complete case analysis

## "included 640 participants with nonmissing values for the covariates and
## outcome, of whom 314 were assigned to the treatment arm and 326 to the
## control arm"
table(cca$Group)

## "Of those in the treatment arm, 50% were nonadherent"
xtabs(~Group+Tx.comp.,data=cca,addNA=TRUE)
mean(cca$Tx.comp.[cca$Group=='T']=='Yes')

xtabs(~Group+Completed.EDC,data=cca,addNA=TRUE)
mean(cca$Completed.EDC[cca$Group=='T']=='Yes')

cca <- cca%>%
  mutate(S=Tx.comp.=='Yes',
         Y=V5..BOP/100,
         Z=Group=='Treatment')#%>%
#  select(-Tx.comp.)#,-Completed.EDC)

############################################################
###
#### calucuate ATE
###
#########################################################



### "the ITT analysis ... -.24 (95% CI: -0.27,-0.21)"
ATE <- t.test(Y~Group,data=cca)
ATE$coef <- coef(lm(Y~Group,data=cca))[2]


############################################################
###
#### GEEPERS
###
#########################################################

geepers <- est(cca,covFormU=~ETXU_CAT1+OFIBRIN1)


##### check principal score model
psMod <- geepers$psMod
arm::binnedplot(predict(psMod,type='response'),resid(psMod,type='response'))
summary(psMod)

print(aucMod(psMod))

### standardized coefficients:
print(glm(S~scale(ETXU_CAT1)+scale(OFIBRIN1),data=mf,family=binomial))


############################################################
###
#### Bayesian mixture model
###
#########################################################

## put the data in the proper format
Xout <- model.matrix(geepers$outMod)
Xout <- Xout[,-c(1,which(colnames(Xout)%in%c('Z','Sp','Z:Sp')))]

sdat <- list(
      YctlY=cca$Y[cca$Z==0],
      YtrtY=cca$Y[cca$Z==1],
      XctlU=model.matrix(formula(geepers$psMod)[-2],data=subset(cca,Z==0))[,-1],
      XtrtU=model.matrix(geepers$psMod)[,-1],
      XctlY=Xout[cca$Z==0,],
      XtrtY=Xout[cca$Z==1,],
      bottomOuter=cca$S[cca$Z==1],
      nc=sum(1-cca$Z),
      nt=sum(cca$Z)
    )
    sdat$ncovU=ncol(sdat$XctlU)
    sdat$ncovY=ncol(Xout)


psStan <- stan(file = 'code/fh2t/psModSimp.stan',data=sdat)
save(psStan,file='code/OPTS/psStanSimp.RData')
load('code/OPTS/psStanSimp.RData')

ppp <- psStan@model_pars
ppp <- ppp[-grep("pi|__",ppp)]
#ppp[10:11] <- paste0(ppp[6:7],'[2]')
#ppp[6:7] <- paste0(ppp[6:7],'[1]')

traceplot(psStan,par=ppp)

print(psStan,par=ppp)

### mixture effect estimates
mixEffs <- summary(psStan)$summary
mixEffs <- mixEffs[grep("ATE",rownames(mixEffs)),]
rownames(mixEffs) <- c(notbottomOuterATE="eff0",
                       bottomOuterATE="eff1",
                       ATEdiff="diff")[rownames(mixEffs)]


############################################################
###
#### PSW
###
#########################################################


PSW <- psw(cca,psMod)
print(PSW)


############################################################
###
#### Plot estimates for paper
###
#########################################################

ateDF <-
    data.frame(EFF=.25,
               estimates=ATE$coef,
               ymin=-ATE$conf.int[2],
               ymax=-ATE$conf.int[1])

### figure
effects <-
bind_rows(
### geepers
    geepers%>%
    effsFromFit()%>%
    data.frame()%>%
    rownames_to_column("eff")%>%
    mutate(method="GEEPERS"),
### PSW
    PSW$coef%>%
    data.frame()%>%
    rownames_to_column("eff")%>%
    rename(estimates=est,SE=se)%>%
    mutate(method="PSW"),
### mixture model
    mixEffs%>%
    data.frame()%>%
    rownames_to_column("eff")%>%
    transmute(eff,estimates=mean,SE=sd,method="Mixture"),
### BSIV (took numbers from the paper)
    data.frame(
        eff=c("eff0","eff1","diff"),
        estimates=c(-.34,-.14,.19),
        SE=c((.44-.26)/4,(.22-0.07)/4,(.36-.05)/4),
        method="BSIV")
)

effects%>%pivot_wider(id_cols=eff,names_from=method,values_from=c(estimates,SE))

effects%>%pivot_wider(id_cols=method,names_from=eff,values_from=c(estimates,SE))%>%
    select(method,ends_with("eff0"),ends_with("eff1"),ends_with("diff"))




optPlot <- effects%>%filter(eff!="diff")%>%
    mutate(
        EFF=ifelse(eff=='eff0', ## Positioning
            ifelse(method=="GEEPERS",-.105,
            ifelse(method=="BSIV",-.035,
            ifelse(method=='Mixture',0.035,.105))),
            ifelse(method=="GEEPERS",.395,
            ifelse(method=="BSIV",.465,
            ifelse(method=='Mixture',.535,.605)))),
        effect=factor(
            c(eff1="Complete Treatment",eff0="Partial Treatment",diff="Difference")[eff],
            levels=c("Partial Treatment","Complete Treatment","Difference")),
        ymin=estimates-2*SE,
        ymax=estimates+2*SE,
        meth=c(Mixture='\\textsc{pmm}',GEEPERS='\\textsc{geepers}',PSW='\\textsc{psw}',BSIV="\\textsc{bsiv}")[method],
        meth=factor(meth,levels=c("\\textsc{geepers}","\\textsc{bsiv}","\\textsc{pmm}","\\textsc{psw}")))%>%
    ggplot(aes(EFF,estimates,color=meth,ymin=ymin,ymax=ymax))+
    geom_point()+
    geom_errorbar(width=0,linewidth=2)+
    geom_hline(yintercept=0)+
    geom_hline(yintercept=ATE$coef,linetype="dashed")+
    geom_point(data=ateDF,aes(EFF,estimates),color='black',inherit.aes=FALSE)+
    geom_errorbar(data=ateDF,aes(EFF,ymin=ymin,ymax=ymax),width=0,color='black',inherit.aes=FALSE)+
    scale_x_continuous("Principal Stratum",breaks=c(0,.25,.5),minor_breaks=NULL,
                       labels=c("Partial\nTreatment","ATE","Complete\nTreatment"),limits=c(-.15,.65))+
        scale_color_manual(values=palette)+
  labs(color='Method',y='Principal Effect')+
  theme(legend.position="top")#%>%print()

tikz("figure/opts.tex",width=5,height=3,standAlone=TRUE)
print(optPlot)
dev.off()


setwd("figure")
system("lualatex opts.tex")
setwd("..")
