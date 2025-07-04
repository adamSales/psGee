library(tidyverse)
library(broom)
library(kableExtra)
library(gridExtra)
library(ggpubr)
library(tikzDevice)


options(
tikzLatexPackages = c(
    "\\usepackage{tikz}",
    "\\usepackage{bm}",
    "\\usepackage{amsmath}",
"\\usepackage[active,tightpage]{preview}",
"\\PreviewEnvironment{pgfpicture}",
"\\setlength\\PreviewBorder{0pt}"
),
tikzXelatexPackages = c(
"\\usepackage{tikz}\n",
"\\usepackage[active,tightpage,xetex]{preview}\n",
"\\usepackage{fontspec,xunicode}\n",
"\\PreviewEnvironment{pgfpicture}\n",
"\\setlength\\PreviewBorder{0pt}\n"
),
tikzLualatexPackages = c(
    "\\usepackage{tikz}\n",
    "\\usepackage{bm}",
        "\\usepackage{amsmath}",
"\\usepackage[active,tightpage,psfixbb]{preview}\n",
"\\usepackage{fontspec,xunicode}\n",
"\\PreviewEnvironment{pgfpicture}\n",
"\\setlength\\PreviewBorder{0pt}\n"
)
)

source('code/simulation/readSimFuncs.r')

#### after simulation results have been loaded and pre-processed...
## source('code/simulation/readSim.r')
load('simResults/fullResults.RData')
load('simResults/resultsNs_mu01is0.RData')
load('simResults/resultsB1s_mu01is0.RData')


###############################################################################
### results across ns
###############################################################################
 pdns <- resultsNs %>%
  group_by(b1)%>%
  mutate(AUCf=mean(auc,na.rm=TRUE))%>%
  ungroup()%>%
  mutate(
    errP = est - pop,
    B1 = factor(b1,levels=c(0,0.2,0.5),
                labels=c(bquote(alpha==0),bquote(alpha==0.2),bquote(alpha==0.5))),
    AUCff=paste0('AUC=',round(AUCf,1)),
    N=paste0('n=',n),
    M1=factor(mu01,levels=c("0","0.3"),
              labels=c(bquote(mu[c]^1-mu[c]^0==0),bquote(mu[c]^1-mu[c]^0==0.3))),
    PE=paste('Stratum',eff),
    dist=paste(c(mix='Mixture',unif='Unform',norm='Normal')[errDist],'Errors'),
    estimator=c(bayes='Mixture',mest='GEEPERs',psw='PSW')[estimator],
    interactionZ=ifelse(intZ,"Z interaction","No\nZ interaction"),
    interactionS=ifelse(intS,"S interaction","No\nS interaction"),
    N=paste0("n=",n)
  ) %>%
     filter(rhat<1.1)


### any statistically significant bias?
pdns%>%filter(eff=="1")%>%
    group_by(estimator,n)%>%
    summarize(bias=mean(errP),pval=t.test(errP)$p.value)%>%
    mutate(p.adj=p.adjust(pval,method="BH"),
           p.holm=p.adjust(pval,method="holm"))%>%
    print(n=Inf)



#### figure in the paper
biasN <- pdns%>%
  filter(rhat<1.1)%>%
  group_by(n,eff,estimator)%>%
    summarize(bias=mean(est-pop,na.rm=TRUE),
            ttest=tidy(t.test(est-pop)))%>%
  ungroup()%>%
  bind_cols(.$ttest)%>%
  select(-ttest)


seN <- pdns%>%
  filter(rhat<1.1)%>%
  group_by(n,eff,estimator)%>%
  summarize(se=sd(errP,na.rm=TRUE))%>%ungroup()

full_join(seN,biasN)%>%filter(eff==1,estimator!="PSW")%>%summarize(across(c(se,bias),max),coef=round(se/bias))

plotByN=bind_rows(
  seN%>%mutate(se=se/3,meas="SE")%>%rename(what=se),
  biasN%>%mutate(meas="Bias")%>%rename(what=bias))%>%
  filter(eff==1,estimator!="PSW")%>%
  ggplot(aes(n,what,color=estimator,group=paste0(estimator,meas),linetype=meas,fill=estimator,shape=meas))+
  geom_point()+geom_line()+geom_hline(yintercept=0)+#,linetype="dotted",size=2)+
  scale_y_continuous(name="Bias",sec.axis=sec_axis(trans=~.*3,name='Standard Error'))+
  scale_shape_manual(values=c(0,16))+
    annotate('text',600,.11,label=list(bquote(paste(#atop(
                                "Normal Resid., No Interactions, ",alpha==0.5))),#", "~mu[C]^1-mu[C]^0==0.3))),
             parse=TRUE)+
  scale_x_continuous(name="Sample Size Per Group (n)",breaks=c(seq(100,700,200),1000))+
  #ggtitle("Bias and Standard Error by n")+
  theme(legend.title=element_blank())
ggsave("simFigs/biasSEbyN.jpg",plot=plotByN,width=6,height=3)




###############################################################################
### results across B1
###############################################################################

 pdb1 <- resultsB1s %>%
  group_by(b1)%>%
  mutate(AUCf=mean(auc,na.rm=TRUE))%>%
  ungroup()%>%
  mutate(
    errP = est - pop,
    B1 = factor(b1,levels=seq(0,1,.1),
                labels=lapply(seq(0,1,.1), function(x) bquote(alpha==.(x)))),
    AUCff=paste0('AUC=',round(AUCf,1)),
    N=paste0('n=',n),
    M1=factor(mu01,levels=c("0","0.3"),
              labels=c(bquote(mu[c]^1-mu[c]^0==0),bquote(mu[c]^1-mu[c]^0==0.3))),
    PE=paste('Stratum',eff),
    dist=paste(c(mix='Mixture',unif='Unform',norm='Normal')[errDist],'Errors'),
    estimator=c(bayes='Mixture',mest='GEEPERs',psw='PSW')[estimator],
    interactionZ=ifelse(intZ,"Z interaction","No\nZ interaction"),
    interactionS=ifelse(intS,"S interaction","No\nS interaction"),
    N=paste0("n=",n)
  ) %>%
  filter(rhat<1.1)

### any statistically significant bias?
pdb1%>%filter(eff=="1",rhat<1.1,estimator!="PSW")%>%
    group_by(estimator,b1)%>%
    summarize(bias=mean(errP),pval=round(t.test(errP)$p.value,4))%>%
    mutate(p.adj=p.adjust(pval,method="BH"),
           p.holm=p.adjust(pval,method="holm"))%>%
    print(n=Inf)



#### figure in the paper

biasB1 <- pdb1%>%
  filter(rhat<1.1)%>%
  group_by(b1,eff,estimator)%>%
  summarize(bias=mean(errP,na.rm=TRUE))%>%ungroup()

seB1 <- pdb1%>%
  filter(rhat<1.1)%>%
  group_by(b1,eff,estimator)%>%
  summarize(se=sd(errP,na.rm=TRUE))%>%ungroup()

full_join(seB1,biasB1)%>%filter(eff==1,estimator!="PSW",b1>0)%>%summarize(across(c(se,bias),max),coef=round(se/bias))

plotByAlpha=bind_rows(
  seB1%>%mutate(se=se/5,meas="SE")%>%rename(what=se),
  biasB1%>%mutate(meas="Bias")%>%rename(what=bias))%>%
  filter(b1>0.1,eff==1,estimator!="PSW")%>%
  ggplot(aes(b1,what,color=estimator,group=paste0(estimator,meas),linetype=meas,fill=estimator,shape=meas))+
  geom_point()+geom_line()+geom_hline(yintercept=0)+#,linetype="dashed")+
    scale_y_continuous(name="Bias",sec.axis=sec_axis(trans=~.*5,name='Standard Error'))+
    scale_x_continuous(breaks=seq(0,1,0.2),minor_breaks=seq(0.1,0.9,0.2),labels=sprintf("%0.1f",seq(0,1,0.2)))+
  scale_shape_manual(values=c(0,16))+
    annotate('text',.7,.08,label=list(bquote(
                               paste(#atop(
                                   "Normal Resid., No Interactions, ",n==500))),# "~mu[C]^1-mu[C]^0==0.3))),
             parse=TRUE)+
    xlab(bquote(alpha))#+ggtitle(bquote("Bias and Standard Error by "~alpha))

ggsave("simFigs/biasSEbyB1.jpg",plot=plotByAlpha,width=6,height=3)


#### combine results by n and B1
ggarrange(plotByN,plotByAlpha,ncol=1,common.legend = TRUE, legend = "bottom")
ggsave("simFigs/biasSEbyB1n.jpg",width=5,height=4)




######################################################################
### main results

pd <- results %>%
  group_by(b1)%>%
  mutate(AUCf=mean(auc,na.rm=TRUE))%>%
  ungroup()%>%
  mutate(
    errP = est - pop,
    B1 = factor(b1,levels=c(0,0.2,0.5),
                labels=c(bquote(alpha==0),bquote(alpha==0.2),bquote(alpha==0.5))),
    AUCff=paste0('AUC=',round(AUCf,1)),
    N=paste0('n=',n),
    M1=factor(mu01,levels=c("0","0.3"),
              ##labels=c(bquote(mu[c]^1-mu[c]^0==0),bquote(mu[c]^1-mu[c]^0==0.3))),
              labels=c(bquote(mu[c]^1==0),bquote(mu[c]^1==0.3))),
    PE=paste('Stratum',eff),
    dist=paste(c(mix='Mixture',unif='Unform',norm='Normal')[errDist],'Errors'),
    #estimator=c(bayes='Mixture',mest='GEEPERs',psw='PSW')[estimator],
    estimator=c(bayes='\\textsc{pmm}',mest='\\textsc{geepers}',psw='\\textsc{psw}')[estimator],
    interactionZ=ifelse(intZ,"Z\ninteraction","No Z\ninteraction"),
    interactionS=ifelse(intS,"S\ninteraction","No S\ninteraction"),
    intAll=ifelse(intZ,ifelse(intS,"$\\bm{x}\\text{:}S_T$\\&$\\bm{x}\\text{:}Z$","$\\bm{x}\\text{:}Z$"),ifelse(intS,"$\\bm{x}\\text{:}S_T$","No inter.")),
    intAll=factor(intAll,levels=c("No inter.","$\\bm{x}\\text{:}Z$", "$\\bm{x}\\text{:}S_T$","$\\bm{x}\\text{:}S_T$\\&$\\bm{x}\\text{:}Z$"))
  ) %>%
  filter(estimator=='PSW'|rhat<1.1)

#######################################################
### figure for paper
#######################################################


norm<-bp(pd,
   subset=b1 > 0& errDist=='norm'& eff==1&mu01==0&n==500&PE=='Stratum 1',
   title="Normal Residuals",ylim=c(-1.5,1.5),#c(-1.5,1.5),
   facet=B1~intAll,#interactionZ+interactionS,
   Labeller=labeller(B1="none",#label_parsed,
                     intAll=label_value),labSize=2.5
   )+
    labs(y=bquote("Estimation Error for "~tau^1),#subtitle="No Interactions",
         x=NULL)+
  #theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          strip.text.y=element_blank(),strip.background.y=element_blank(),
                                        #plot.subtitle = element_text(size=10),
        legend.position="none")+
  scale_fill_manual(values=c('#1b9e77','#d95f02','#7570b3'))+
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'))

unif<-bp(pd,
   subset=b1 > 0& errDist=='unif'& eff==1&mu01==0&n==500&PE=='Stratum 1',
   title="Uniform Residuals",ylim=c(-1.5,1.5),#c(-1.5,1.5),
   facet=B1~intAll,#interactionZ+interactionS,
   Labeller=labeller(B1=label_parsed,intAll=label_value),labSize=2.
   )+
    labs(y=NULL,#bquote("Estimation Error for "~tau^1),#subtitle="No Interactions",
         x=NULL)+
  #theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                                        #plot.subtitle = element_text(size=10),
        legend.position="none")+
  scale_fill_manual(values=c('#1b9e77','#d95f02','#7570b3'))+
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'))

#pdf(
tikz("simFigs/boxplotsNew.tex",#pdf",
     width=6.5,height=4,standAlone=TRUE)
grid.arrange(norm,unif,nrow=1)
dev.off()

setwd("simFigs")
system("lualatex boxplotsNew.tex")
setwd("..")


#######################################################
### "boxplots" figure with n=1000 for appendix
#######################################################
norm1000<-bp(pd,
   subset=b1 > 0& errDist=='norm'& eff==1&mu01==0&n==1000&PE=='Stratum 1',
   title="Normal Residuals",ylim=c(-1.5,1.5),#c(-1.5,1.5),
   facet=B1~intAll,#interactionZ+interactionS,
   Labeller=labeller(B1="none",#label_parsed,
                     intAll=label_value),labSize=2.5
   )+
    labs(y=bquote("Estimation Error for "~tau^1),#subtitle="No Interactions",
         x=NULL)+
  #theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          strip.text.y=element_blank(),strip.background.y=element_blank(),
                                        #plot.subtitle = element_text(size=10),
        legend.position="none")+
  scale_fill_manual(values=c('#1b9e77','#d95f02','#7570b3'))+
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'))

unif1000<-bp(pd,
   subset=b1 > 0& errDist=='unif'& eff==1&mu01==0&n==1000&PE=='Stratum 1',
   title="Uniform Residuals",ylim=c(-1.5,1.5),#c(-1.5,1.5),
   facet=B1~intAll,#interactionZ+interactionS,
   Labeller=labeller(B1=label_parsed,intAll=label_value),labSize=2.
   )+
    labs(y=NULL,#bquote("Estimation Error for "~tau^1),#subtitle="No Interactions",
         x=NULL)+
  #theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                                        #plot.subtitle = element_text(size=10),
        legend.position="none")+
  scale_fill_manual(values=c('#1b9e77','#d95f02','#7570b3'))+
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'))

#pdf(
tikz("simFigs/boxplotsAppendixN1000.tex",#pdf",
     width=6.5,height=4,standAlone=TRUE)
grid.arrange(norm1000,unif1000,nrow=1)
dev.off()

setwd("simFigs")
system("lualatex boxplotsAppendixN1000.tex")
setwd("..")


#######################################################
#### rmse appendix
#######################################################

rmse <- pd%>%
  filter(rhat<1.1|estimator=='PSW')%>%
  group_by(n,mu01,errDist,b1,intS,intZ,interactionS,interactionZ,eff,estimator)%>%
   summarize(rmse=sqrt(mean(errP^2,na.rm=TRUE)))%>%ungroup()



sink('writeUps/rmseTabAppendix500.tex')
cbind(
rmse%>%filter(errDist!='mix',n==500,eff!="Diff",b1==.2,mu01==0)%>%
transmute(
    #Parameter=ifelse(eff=="0","$\\eff0$","$\\eff1$"),
    `Residual\nDist.`=c(norm='Normal',unif='Uniform')[errDist],
         `$\\bm{x}:Z$\nInt.?`=ifelse(intZ,'Yes','No'),
         `$\\bm{x}:S_T$\nInt.?`=ifelse(intS,'Yes','No'),
         estimator,#=c(Mixture="\\pmm",GEEPERs="\\geepers",PSW="\\psw")[estimator],
#         `$\\beta_1$`=mu01,
         Parameter=paste0('$\\tau^',eff,'$'),
         rmse)%>%
pivot_wider(names_from=estimator,values_from=rmse),
rmse%>%filter(errDist!='mix',n==500,mu01==0,eff!="Diff",b1==.5)%>%
  transmute(`Res.\nDist.`=c(norm='Norm.',unif='Unif.')[errDist],
         `$\\bm{x}:Z$\nInt.?`=ifelse(intZ,'Yes','No'),
         `$\\bm{x}:S_T$\nInt.?`=ifelse(intS,'Yes','No'),
         estimator,#=c(Mixture="\\pmm",GEEPERs="\\geepers",PSW="\\psw")[estimator],
         n,
         #`$\\beta_1$`=mu01,
         Parameter=paste0('$\\tau^',eff,'$'),
         rmse,
         xx=rep(1:(n()/3),each=3))%>%
select(xx,estimator,rmse)%>%
pivot_wider(names_from=estimator,values_from=rmse)%>%
select(-xx))%>%#GEEPERs,Mix.,PSW))%>%
    kbl('latex',booktabs=TRUE,col.names=linebreak(names(.)),escape=FALSE,digits=2)%>%
    add_header_above(c(" " = 4,  "$\\\\alpha=0.2$" = 3, "$\\\\alpha=0.5$" = 3),escape=FALSE)%>%
     add_header_above(c(" " = 4,  "$n=500$"=6),escape=FALSE)%>%
  collapse_rows(columns=1,latex_hline="major",valign="middle")%>%print()
sink()


sink('writeUps/rmseTabAppendix1000.tex')
cbind(
rmse%>%filter(mu01==0,errDist!='mix',n==1000,eff!="Diff",b1==.2)%>%
transmute(
    #Parameter=ifelse(eff=="0","$\\eff0$","$\\eff1$"),
    `Residual\nDist.`=c(norm='Normal',unif='Uniform')[errDist],
         `$\\bm{x}:Z$\nInt.?`=ifelse(intZ,'Yes','No'),
    `$\\bm{x}:S_T$\nInt.?`=ifelse(intS,'Yes','No'),
     Parameter=ifelse(eff=="0","$\\eff0$","$\\eff1$"),
         estimator,#=c(Mixture="\\pmm",GEEPERs="\\geepers",PSW="\\psw")[estimator],
         #`$\\beta_1$`=mu01,
         #`Prin.\nEff`=paste0('$\\tau^',eff,'$'),
         rmse)%>%
pivot_wider(names_from=estimator,values_from=rmse),
rmse%>%filter(errDist!='mix',n==1000,mu01==0,eff!="Diff",b1==.5)%>%
  transmute(`Res.\nDist.`=c(norm='Norm.',unif='Unif.')[errDist],
         `$\\bm{x}:Z$\nInt.?`=ifelse(intZ,'Yes','No'),
         `$\\bm{x}:S_T$\nInt.?`=ifelse(intS,'Yes','No'),
         estimator,#=c(Mixture="\\pmm",GEEPERs="\\geepers",PSW="\\psw")[estimator],
         n,
         #`$\\beta_1$`=mu01,
         Parameter=paste0('$\\tau^',eff,'$'),
         rmse,
         xx=rep(1:(n()/3),each=3))%>%
select(xx,estimator,rmse)%>%
pivot_wider(names_from=estimator,values_from=rmse)%>%
select(-xx))%>%#GEEPERs,Mix.,PSW))%>%
    kbl('latex',booktabs=TRUE,col.names=linebreak(names(.)),escape=FALSE,digits=2)%>%
  add_header_above(c(" " = 4,  "$\\\\alpha=0.2$" = 3, "$\\\\alpha=0.5$" = 3),escape=FALSE)%>%
add_header_above(c(" " = 4,  "$n=1000$"=6),escape=FALSE)%>%
    collapse_rows(columns=1,latex_hline="major",valign="middle")%>%print()
sink()



#######################################################
#### coverage
#######################################################
 coverage <- pd%>%
  filter(rhat<1.1)%>%
  group_by(n,mu01,errDist,b1,interactionS,interactionZ,intS,intZ,eff,estimator)%>%
   summarize(coverage=mean(CInormL<=pop & CInormU>=pop,na.rm=TRUE))%>%ungroup()

redCov <- function(coverage) paste0("\\rd{",sprintf("%.2f",coverage),"}")
condRed <- function(coverage,intZ,intS,estimator,errDist,b1=1)
    ifelse(intZ|intS,redCov(coverage),
    ifelse(errDist=="unif"&estimator=="\\textsc{pmm}",redCov(coverage),
    ifelse(estimator=="\\textsc{geepers}"&b1==0,redCov(coverage),sprintf("%.2f",coverage))))


sink('writeUps/coverageTab.tex')
cbind(
     coverage%>%
      mutate(coverage=condRed(coverage,intZ,intS,estimator,errDist,b1))%>%
      filter(n==500,eff==1,errDist!='mix',mu01==0,b1==0,estimator!="\\textsc{psw}")%>%
     transmute(
             Parameter=ifelse(eff=="0","$\\eff0$",ifelse(eff=="1","$\\eff1$","$\\eff1-\\eff0$")),
         `Residual\nDist.`=c(norm='Normal',unif='Uniform')[errDist],
         `$\\bm{x}:Z$\nInt.?`=ifelse(intZ,'Yes','No'),
         `$\\bm{x}:S_T$\nInt.?`=ifelse(intS,'Yes','No'),
         estimator,#=c(Mixture="\\pmm",GEEPERs="\\geepers")[estimator],
         coverage)%>%
    pivot_wider(names_from=estimator,values_from=coverage),
   coverage%>%filter(n==500,eff==1,errDist!='mix',mu01==0,b1==.2)%>%
      mutate(coverage=condRed(coverage,intZ,intS,estimator,errDist))%>%
    pivot_wider(names_from=estimator,values_from=coverage)%>%
    select(`\\textsc{geepers}`,`\\textsc{pmm}`),
    coverage%>%filter(n==500,eff==1,errDist!='mix',mu01==0,b1==.5)%>%
      mutate(coverage=condRed(coverage,intZ,intS,estimator,errDist))%>%
    pivot_wider(names_from=estimator,values_from=coverage)%>%
    select(`\\textsc{geepers}`,`\\textsc{pmm}`)#%>%
#    rename("\\pmm"="Mixture","\\geepers"="GEEPERs")
  )%>%
    kbl('latex',booktabs=TRUE,col.names=linebreak(names(.)),escape=FALSE,digits=2)%>%
    add_header_above(c(" " = 3, "$\\\\alpha=0$" = 2,
                       "$\\\\alpha=0.2$" = 2, "$\\\\alpha=0.5$" = 2),escape=FALSE)%>%
    collapse_rows(columns=1,latex_hline="major",valign="middle")%>%
    footnote(general=c("\\\\footnotesize Based on 500 replications. $n=500$. Simulation standard error $\\\\approx 1$ percentage point. Estimates colored \\\\rd{red} indicate cases where the assumptions of the model are not met."),escape=FALSE,footnote_as_chunk = TRUE,threeparttable=TRUE)%>%print()
  sink()

########################
## coverage tab for appendix: n=500
########################

sink('writeUps/coverageTabAppendix500.tex')

cbind(
    coverage%>%filter(mu01==0,errDist!='mix',eff!="Diff",b1==0,estimator!='\\textsc{psw}',n==500,
                      )%>%
    transmute(
            #Parameter=ifelse(eff=="0","\\eff0","\\eff1"),
        `Residual\nDist.`=c(norm='Normal',unif='Uniform')[errDist],
         `$\\bm{x}:Z$\nInt.?`=ifelse(intZ,'Yes','No'),
         `$\\bm{x}:S_T$\nInt.?`=ifelse(intS,'Yes','No'),
         estimator,#=ifelse(estimator=='Mixture','\\pmm',"\\geepers"),
         #`$\\beta_1$`=mu01,
         Parameter=paste0('$\\tau^',eff,'$'),
         coverage=sprintf("%.2f", coverage),#round(coverage,2),
         coverage=ifelse(intZ|intS|(errDist=="unif"&estimator=="\\textsc{pmm}"),paste0("\\rd{",coverage,"}"),coverage))%>%
pivot_wider(names_from=estimator,values_from=coverage),
coverage%>%filter(mu01==0,errDist!='mix',eff!="Diff",b1==.2,estimator!='PSW',n==500)%>%
  transmute(`Residual\nDist.`=c(norm='Normal',unif='Uniform')[errDist],
         `$\\bm{x}:Z$\nInt.?`=ifelse(intZ,'Yes','No'),
         `$\\bm{x}:S_T$\nInt.?`=ifelse(intS,'Yes','No'),
         estimator,#=ifelse(estimator=='Mixture','\\pmm',"\\geepers"),
         #`$\\beta_1$`=mu01,
         Parameter=paste0('$\\tau^',eff,'$'),
         coverage=sprintf("%.2f", coverage),#round(coverage,2),
         coverage=ifelse(intZ|intS|(errDist=="unif"&estimator=="\\textsc{pmm}"),paste0("\\rd{",coverage,"}"),coverage))%>%
pivot_wider(names_from=estimator,values_from=coverage)%>%
select(`\\textsc{geepers}`,`\\textsc{pmm}`),
coverage%>%filter(mu01==0,errDist!='mix',n==500,eff!="Diff",b1==.5,estimator!='\\textsc{psw}')%>%
  transmute(`Residual\nDist.`=c(norm='Normal',unif='Uniform')[errDist],
         `$\\bm{x}:Z$\nInt.?`=ifelse(intZ,'Yes','No'),
         `$\\bm{x}:S_T$\nInt.?`=ifelse(intS,'Yes','No'),
         estimator,#=ifelse(estimator=='Mixture','\\pmm',"\\geepers"),
         n,
         #`$\\beta_1$`=mu01,
         Parameter=paste0('$\\tau^',eff,'$'),
         coverage=sprintf("%.2f", coverage),#round(coverage,2),
         coverage=ifelse(intZ|intS|(errDist=="unif"&estimator=="\\textsc{pmm}"),paste0("\\rd{",coverage,"}"),coverage))%>%
pivot_wider(names_from=estimator,values_from=coverage)%>%
select(`\\textsc{geepers}`,`\\textsc{pmm}`))%>%
kbl('latex',booktabs=TRUE,col.names=linebreak(names(.)),escape=FALSE,digits=2)%>%
    add_header_above(c(" " = 4, "$\\\\alpha=0$" = 2, "$\\\\alpha=0.2$" = 2, "$\\\\alpha=0.5$" = 2),escape=FALSE)%>%
    add_header_above(c(" " = 4, "$n=500$"=6),escape=FALSE)%>%
  collapse_rows(columns=1,latex_hline="major",valign="middle")%>%print()
sink()

########################
## coverage tab for appendix: n=1000
########################
sink('writeUps/coverageTabAppendix1000.tex')
cbind(
    coverage%>%filter(mu01==0,errDist!='mix',eff!="Diff",b1==0,estimator!='\\textsc{psw}',n==1000
                      )%>%
  transmute(`Residual\nDist.`=c(norm='Normal',unif='Uniform')[errDist],
         `$\\bm{x}:Z$\nInt.?`=ifelse(intZ,'Yes','No'),
         `$\\bm{x}:S_T$\nInt.?`=ifelse(intS,'Yes','No'),
         estimator,#=ifelse(estimator=='Mixture','\\pmm',"\\geepers"),
         #`$\\beta_1$`=mu01,
         Parameter=paste0('$\\tau^',eff,'$'),
         coverage=sprintf("%.2f", coverage),#round(coverage,2),
         coverage=ifelse(intZ|intS|(errDist=="unif"&estimator=="\\textsc{pmm}"),paste0("\\rd{",coverage,"}"),coverage))%>%
pivot_wider(names_from=estimator,values_from=coverage),
coverage%>%filter(mu01==0,errDist!='mix',eff!="Diff",b1==.2,estimator!='PSW',n==1000)%>%
  transmute(`Residual\nDist.`=c(norm='Normal',unif='Uniform')[errDist],
         `$\\bm{x}:Z$\nInt.?`=ifelse(intZ,'Yes','No'),
         `$\\bm{x}:S_T$\nInt.?`=ifelse(intS,'Yes','No'),
         estimator,#=ifelse(estimator=='Mixture','\\pmm',"\\geepers"),
         #`$\\beta_1$`=mu01,
         Parameter=paste0('$\\tau^',eff,'$'),
         coverage=sprintf("%.2f", coverage),#round(coverage,2),
         coverage=ifelse(intZ|intS|(errDist=="unif"&estimator=="\\textsc{pmm}"),paste0("\\rd{",coverage,"}"),coverage))%>%
pivot_wider(names_from=estimator,values_from=coverage)%>%
select(`\\textsc{geepers}`,`\\textsc{pmm}`),
coverage%>%filter(mu01==0,errDist!='mix',n==1000,eff!="Diff",b1==.5,estimator!='\\textsc{psw}')%>%
  transmute(`Residual\nDist.`=c(norm='Normal',unif='Uniform')[errDist],
         `$\\bm{x}:Z$\nInt.?`=ifelse(intZ,'Yes','No'),
         `$\\bm{x}:S_T$\nInt.?`=ifelse(intS,'Yes','No'),
         estimator,#=ifelse(estimator=='Mixture','\\pmm',"\\geepers"),
         n,
         #`$\\beta_1$`=mu01,
         Parameter=paste0('$\\tau^',eff,'$'),
         coverage=sprintf("%.2f", coverage),#round(coverage,2),
         coverage=ifelse(intZ|intS|(errDist=="unif"&estimator=="\\textsc{pmm}"),paste0("\\rd{",coverage,"}"),coverage))%>%
pivot_wider(names_from=estimator,values_from=coverage)%>%
select(`\\textsc{geepers}`,`\\textsc{pmm}`))%>%
kbl('latex',booktabs=TRUE,col.names=linebreak(names(.)),escape=FALSE,digits=2)%>%
    add_header_above(c(" " = 4, "$\\\\alpha=0$" = 2, "$\\\\alpha=0.2$" = 2, "$\\\\alpha=0.5$" = 2),escape=FALSE)%>%
    add_header_above(c(" " = 4, "$n=1000$"=6),escape=FALSE)%>%
  collapse_rows(columns=1,latex_hline="major",valign="middle")%>%print()
sink()




### how does AUC very w b1 and n?
resultsB1s%>%
  group_by(b1)%>%
  mutate(meanAUC=mean(auc))%>%
  ggplot(aes(as.factor(b1),auc))+geom_boxplot()+geom_point(aes(y=meanAUC))+geom_smooth(se=FALSE)+ylab('AUC')+xlab(bquote(alpha))
  ggsave('simFigs/alphaAUC.pdf',width=5,height=4)
                                        #  scale_y_continuous('Avg. AUC',seq(.5,1,.1))
