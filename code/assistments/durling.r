load('assistmentsData/dataWcontrasts.RData')

### find experiments where thea durling made 2 hints, one with vid one w/o
durlingVid <- which(dat$contains_video&grepl('durling',dat$tutoring_creator))

durlingTxt <- dat$assigned_tsid[!dat$contains_video&grepl('durling',dat$tutoring_creator)]

durlingTsid <- c(durlingVid,durlingTxt)

durlingVidVSdurlingTxt <- 
  which(
    (dat$alternative_tsid_1[durlingVid]%in%durlingTxt)|
      (dat$alternative_tsid_2[durlingVid]%in%durlingTxt)|
      (dat$alternative_tsid_3[durlingVid]%in%durlingTxt)|
      (dat$alternative_tsid_4[durlingVid]%in%durlingTxt))

durling1 <- which(dat$alternative_tsid_1%in%durlingTxt)
durling2 <- which(dat$alternative_tsid_2%in%durlingTxt)
durling3 <- which(dat$alternative_tsid_3%in%durlingTxt)
durling4 <- which(dat$alternative_tsid_4%in%durlingTxt)



durlingAlts <- unique(dat$alts[intersect(durlingVid,c(durling1,durling2,durling3,durling4))])

durling <- dat%>%
  ungroup()%>%
  filter(alts%in%durlingAlts,grepl('durling',tutoring_creator))


