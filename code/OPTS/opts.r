library(rstan)
library(dplyr)
library(tibble)
library(ggplot2)
library(sandwich)
library(lmtest)
library(tableone)
library(purrr)
library(kableExtra)
library(tidyr)

rstan_options(auto_write = FALSE)


source('code/regression.r')
options(mc.cores = 2)

auc <- function(x,y)
  if(length(unique(y))==2 & length(y)==length(x)){
    wilcox.test(x~y)$statistic/(sum(y)*(sum(1-y)))
  }else
    wilcox.test(x,y)$statistic/(legnth(x)*length(y))

aucMod <- function(mod)
  auc(mod$linear,1-mod$y)


print(load('code/OPTS/OPT_Study_PersonLevel_Data.RData'))
## downloaded from https://www.causeweb.org/tshs/obstetrics-and-periodontal-therapy/



opt <- opt%>%
  mutate(OFIBRIN1=as.numeric(as.character(OFIBRIN1)),
         ETXU_CAT1=as.numeric(as.character(ETXU_CAT1)))

levels(opt$Group)=c(T="Treatment",C="Control")[levels(opt$Group)]

xtabs(~is.na(OFIBRIN1)+is.na(ETXU_CAT1)+is.na(V5..BOP),data=opt)

sum(is.na(opt$V5..BOP))
sum((is.na(opt$ETXU_CAT1)|is.na(opt$OFIBRIN1))&!is.na(opt$V5..BOP))

cca <- opt%>%filter(!is.na(OFIBRIN1),!is.na(ETXU_CAT1),!is.na(V5..BOP))%>%
    #select(Group, OFIBRIN1,ETXU_CAT1,V5..BOP,Tx.comp.,Completed.EDC,BL..BOP)%>%
    mutate(S=Tx.comp.=='Yes',
         Y=V5..BOP/100,
         Z=Group=='Treatment')

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


kbl(tab1[-c(1:2),],format="latex",booktabs=TRUE,col.names=tab1[2,],caption="Descriptive statistics---mean and standard deviation or count and percent---for study variables in the full OPT dataset and in the analysis sample (i.e. complete cases)",label="optTab1" )%>%
    add_header_above(c(" ","Full Data"=2,"Complete Cases"=2))



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


##### table for appendix



### "the ITT analysis ... -.24 (95% CI: -0.27,-0.21)"
ATE <- t.test(Y~Group,data=cca)
ATE$coef <- coef(lm(Y~Group,data=cca))[2]

### geepers
geepers <- est(cca,~ETXU_CAT1+OFIBRIN1)


bs=matrix(nrow=1000,ncol=3)
for(i in 1:1000){
  pe=pointEst(cca[sample(1:nrow(cca),nrow(cca),replace=TRUE),],~ETXU_CAT1+OFIBRIN1)
  coefs <- coef(pe$outMod)
  bs[i,] <- c(coefs['ZTRUE'],coefs['ZTRUE']+coefs['ZTRUE:Sp'],coefs['ZTRUE:Sp'])
}
apply(bs,2,sd)

psMod <- geepers$psMod
arm::binnedplot(predict(psMod,type='response'),resid(psMod,type='response'))
summary(psMod)

aucMod(psMod)

## bs <- replicate(5000,
##                 coef(
##                   pointEst(
##                     cca[sample(1:nrow(cca),nrow(cca),replace=TRUE),],
##                     ~ETXU_CAT1+OFIBRIN1)$outMod
##                 )['ZTRUE']
##                 )

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


  psStan <- stan(file = 'code/fh2t/psModSimp.stan')#,data=sdat)
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


PSW <- psw(cca,psMod)
print(PSW)

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

effects%>%filter(eff!="diff")%>%
    mutate(
        EFF=ifelse(eff=='eff0',
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
        ymax=estimates+2*SE)%>%
    ggplot(aes(EFF,estimates,color=method,ymin=ymin,ymax=ymax))+
    geom_point()+
    geom_errorbar(width=0)+
    geom_hline(yintercept=ATE$coef)+
    geom_point(data=ateDF,aes(EFF,estimates),color='black',inherit.aes=FALSE)+
    geom_errorbar(data=ateDF,aes(EFF,ymin=ymin,ymax=ymax),width=0,color='black',inherit.aes=FALSE)+
    scale_x_continuous("Principal Stratum",breaks=c(0,.25,.5),minor_breaks=NULL,
                     labels=c("Partial\nTreatment","ATE","Complete\nTreatment"),limits=c(-.15,.65))+
  labs(color='Method',y='Principal Effect')+
  theme(legend.position="top")
ggsave("figure/opts.pdf",width=5,height=3,units="in")







cca2 <- cca
cca2$Y <- qlogis(cca$Y)
cca2 <- subset(cca2,is.finite(Y))
nrow(cca2)

est(cca2,~ETXU_CAT1+OFIBRIN1)


### other analyses: try different baseline variables
### there are 36, so choose(36,2)=630 combinations of 2
