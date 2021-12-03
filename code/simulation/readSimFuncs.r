

library(tidyverse)
library(rstan)
library(parallel)

bias <- function(sss){

    sss1 <- sss[,1:6]
    sss2 <- sss[,7:11]

    colMeans(cbind(sss1[,3:6]-sss1[,c(1,2,1,2)],sss2))
}

rmse <- function(sss){

    sss1 <- sss[,1:6]
    sss2 <- sss[,7:11]

    c(
        sqrt(colMeans((sss1[,3:6]-sss1[,c(1,2,1,2)])^2)),
        colMeans(sss2)
    )
}


summ <- function(sss){
    sss1 <- sss[,1:6]
    sss2 <- sss[,7:11]
    rbind(colMeans(sss1),
          apply(sss1,2,sd),
          sqrt(row1Means((sss-sss[rep(1:2,3),])^2))
          )
    }

muEff <- function(facs,eff)
    with(as.list(facs),
         ifelse(eff==0,mu10,
         ifelse(eff==1,mu11-mu01,mu11-mu01-mu10))
         )


bayesProc1 <- function(eff, res1) {
  if (rownames(res1$mest)[3] == 'diff')
    rownames(res1$mest)[3] <- 'effDiff'
  with(res1,
       map_dfr(
         list(bayes, mest),
         ~
           tibble(
             eff = eff,
             pop = muEff(facs, eff),
             samp = ifelse(eff == 'Diff', true['S1'] -true['S0'], true[paste0('S', eff)]),
             estimator = ifelse('mean' %in% colnames(.), 'bayes', 'mest'),
             est = .[paste0('eff', eff), ifelse('mean' %in% colnames(.), 'mean', 'estimates')],
             se = .[paste0('eff', eff), ifelse('sd' %in%colnames(.), 'sd', 'SE')],
             CInormL = est - 2 *se,
             CInormU = est + 2 *se,
             CIpercL = ifelse('2.5%' %in% colnames(.), .[paste0('eff', eff), '2.5%'], CInormL),
             CIpercU = ifelse('97.5%' %in% colnames(.), .[paste0('eff', eff), '97.5%'], CInormU),
             rhat = bayes[paste0('eff', eff), 'Rhat'],
             auc = attr(mest, 'auc')
           )
       ))
}



bayesProc <- function(res1,facs){
    if(inherits(res1,'try-error')) return(as_tibble(facs))
    facs%>%
        rbind()%>%
        as_tibble()%>%
        bind_cols(
            map_dfr(c(0,1,'Diff'),bayesProc1,res1=res1)
        )
}



loadRes <- function(ext1='',ext2=''){
    load(paste0('simResults',ext1,'/cases',ext2,'.RData'))

    results <- list()
#    summ <- NULL
    #fn <- list.files('./simResults','sim[[:alnum:]]+.RData')
    for(i in 1:nrow(cases)){#length(fn)){
        #if(i %% 10==0)
            cat(round(i/nrow(cases)*100),'%',sep='')#length(fn)*100), '% ')
        load(paste0('simResults',ext1,'/sim',i,ext2,'.RData'))#fn[i]))
        #stopifnot(identical(facs,cases[i,]))
        resT <- map_dfr(res,bayesProc,facs=facs)
        resT$run <- i
        results[[i]] <- resT
        rm(res,facs)
    }

    reduce(results,bind_rows)
}

bp <- function(pd,subset,facet,title=deparse(substitute(subset))){
  
  r <- if (missing(subset)) 
    rep_len(TRUE, nrow(pd))
  else {
    e <- substitute(subset)
    r <- eval(e, pd, parent.frame())
    if (!is.logical(r)) 
      stop("'subset' must be logical")
    r & !is.na(r)
  }
  
  pd <- pd[r, , drop = TRUE]
  
  ggplot(pd,
         aes(
           x = estimator,
           y = errP,
           fill = estimator,
           color = estimator
         )) +
    geom_jitter(alpha = 0.1) +
    geom_violin(#position = "dodge2",
      color = 'black',
      #outlier.shape = NA
      draw_quantiles = 0.5) +
    geom_hline(yintercept = 0) +
    coord_cartesian(ylim = quantile(pd$errP, c(0.01, 0.99))) +
    facet_grid(facet , scales = "free") +
    labs(title = title,
         x = NULL, y = 'Estimation Error') + theme(legend.pos = 'none')
  
}


load1 <- function(i,ext1='',ext2=''){
  load(paste0('simResults',ext1,'/sim',i,ext2,'.RData'))#fn[i]))
  list(res=res,facs=facs)
}

justLoad <- function(ext1='',ext2=''){
  load(paste0('simResults',ext1,'/cases',ext2,'.RData'))
  
  resList <- lapply(1:nrow(cases),load1)
  
  list(resList=resList,cases=cases)
}

justProc <- function(resList)
  map(1:length(resList),
      function(i){
        cat(i,'/',length(resList),' ')
        resT <- try(map_dfr(resList[[i]]$res,bayesProc,facs=resList[[i]]$facs))
        if(inherits(resT,'try-error')) resT <- as_tibble(resList[[i]]$facs)
        resT$run <- i
        resT
      })

