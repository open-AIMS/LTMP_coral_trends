source('CoralTrends_functions.R')
CoralTrends_checkPackages()

## ---- processManta_getData
load(file=paste0(DATA_PATH, 'primary/manta.RData'))
pcode.mod = read.csv(paste0(DATA_PATH, 'primary/P_CODE_MOD.csv'), strip.white=TRUE)
manta = manta %>% left_join(pcode.mod) %>% droplevels
## ----end


#####################################################################
## Data are collected per Manta Tow.                               ##
## The most appropriate unit for these analyses is the reef level. ##
## - Only include data collected after 1985                        ##
## - Convert data into percent cover                               ##
## - Summarize the data to the reef/year level                     ##
#####################################################################
## ---- processManta_towlevel
manta.tow <- manta %>%
  filter(REPORT_YEAR>1985) %>%
  mutate(Cover=CoralTrends_calcPercent(LIVE_CORAL)) %>%
  mutate(Latitude=REEF_LAT, Longitude=REEF_LONG) %>%
  CoralTrends_calc3ZoneLocation() %>%
  filter(!is.na(Region)) %>% droplevels %>%
  mutate(Region=factor(Region, levels=c('Northern GBR', 'Central GBR', 'Southern GBR')),
         Zone=Region) %>%
  mutate(Year=factor(REPORT_YEAR)) %>%
  ## mutate(Cover=ifelse(Cover==0,0.0001,Cover)) %>%
  mutate(Cover=ifelse(Cover==0,0.0001,Cover)) %>%
  mutate(P_CODE.mod=factor(ifelse(is.na(P_CODE.mod),'Other',P_CODE.mod))) %>%
  as.data.frame
## ----end
##
## ---- processManta_aggregateToReef
manta.sum <- manta %>%
    filter(REPORT_YEAR>1985) %>%
    mutate(Cover=CoralTrends_calcPercent(LIVE_CORAL)) %>%
    group_by(P_CODE.mod,A_SECTOR,SHELF,REEF_NAME,REEF_ID,REPORT_YEAR) %>%
    summarise(Cover=mean(Cover, na.rm=TRUE),
              CoverCat.median=median_cover_cat(LIVE_CORAL),
              Tows=length(unique(TOW_SEQ_NO)),
              Latitude=mean(REEF_LAT, na.rm=TRUE),
              Longitude=mean(REEF_LONG, na.rm=TRUE)) %>%
        mutate(Cover_from_cat=CoralTrends_calcPercent(CoverCat.median)) %>%
    ungroup %>%
    group_by(REEF_NAME) %>%
    ## mutate(VisitFlag = n()<4) %>%  # apply a flag that indicates TRUE when the
    ##                                # reef has fewer than 4 visits
    ungroup %>%
    CoralTrends_calc3ZoneLocation() %>%
    filter(!is.na(Region)) %>% droplevels %>%
    mutate(Region=factor(Region, levels=c('Northern GBR', 'Central GBR', 'Southern GBR')),
           Zone=Region) %>%
    as.data.frame
## ----end

## ---- Compare tow and reef level
bfun <- function(x) {
    if (length(x$Cover)>3 & !all(x$Cover==x$Cover[1])) {
        betareg::betareg(Cover ~ 1, data=x) %>% coef() %>% `[[`(1) %>% plogis()
    } else {
        x$Cover[1]
    }
}

manta.stats <- manta.tow %>%
    mutate(RN=REEF_NAME,YR=Year) %>%
    group_by(REEF_NAME, Year) %>%
    nest() %>%
    mutate(Mean=map_dbl(data, ~mean(.x$Cover)),
           Median=map_dbl(data, ~median(.x$Cover)),
           Logis=map_dbl(data, ~.x$Cover %>% logit %>% mean %>% plogis),
           Beta=map_dbl(data, ~bfun(.x))
           ) %>%
    dplyr::select(-data) %>%
    unnest(cols=c()) %>%
    relocate(Mean, Beta,.after=last_col()) %>%
    arrange(Year) %>%
    ungroup() %>%
    as.data.frame()


manta.sum %>% dplyr::select(REEF_NAME, REPORT_YEAR, Cover) %>%
    left_join(manta.stats %>% mutate(REPORT_YEAR=as.integer(as.character(Year))) %>% dplyr::select(REEF_NAME, REPORT_YEAR, Mean, Beta)) %>%
    head
## ----end

writeLines(paste0("NumberOfReefs=",length(unique(manta.sum$REEF_NAME)),
                  "\nNumberOfSurveys=",nrow(manta.sum)),
           con='../data/processed/Manta.properties')
save(manta.sum, file='../data/processed/manta.sum.RData')
save(manta.tow, file='../data/processed/manta.tow.RData')



#
if (1==2) {
  ## Data are collected per Manta Tow.
  ## The most appropriate unit for these analyses is the reef level.
  ## - Only include data collected after 1985
  ## - Convert data into percent cover
  ## - Summarize the data to the reef/year level
  manta.sum = manta %>%
    filter(REPORT_YEAR>1985, REPORT_YEAR< (finalYear + 1)) %>%
    mutate(Cover=CoralTrends_calcPercent(LIVE_CORAL)) %>%
    group_by(P_CODE.mod,A_SECTOR,SHELF,REEF_NAME,REEF_ID,REPORT_YEAR) %>%
    summarise(Cover=mean(Cover, na.rm=TRUE), Tows=length(unique(TOW_SEQ_NO)), Latitude=mean(REEF_LAT, na.rm=TRUE), Longitude=mean(REEF_LONG, na.rm=TRUE)) %>%
    ungroup()


  ## Include only those reefs that have had more than 4 visits
  wch <- table(manta.sum$REEF_NAME)
  wch<-names(wch[wch>4])
  manta.sum = manta.sum %>% filter(REEF_NAME %in% wch)


  ## Put the reefs into Locations (Zones)

                                        #manta.sum = manta.sum %>% mutate(Location=CoralTrends_calc3ZoneLocations(Latitude))
  manta.sum = CoralTrends_calc3ZoneLocations(manta.sum)
  ## which reefs are excluded...
                                        #as.character(unique(filter(manta.sum, Location=='Outside')$REEF_NAME))
  ## Reefs exclude (Outside)
  exc=manta.sum %>% filter(Location=='Outside')
  manta.sum = manta.sum %>% filter(Location!='Outside') %>% droplevels %>% mutate(Zone=factor(Location, levels=c('Northern','Central','Southern')))
  ## number of reefs used
  length(unique(manta.sum$REEF_NAME))
  ## number of reef surveys
  nrow(manta.sum)
  writeLines(paste0("NumberOfReefs=",length(unique(manta.sum$REEF_NAME)),
                    "\nNumberOfSurveys=",nrow(manta.sum)),
             con='data/processed/Manta.properties')
  save(manta.sum, file='data/processed/manta.sum.RData')



  ##For comparison, we can also attempt to replicate what Glenn had done...
  manta.sum.G = manta %>%
    filter(REPORT_YEAR>1985,REPORT_YEAR<2013) %>%
    mutate(Cover=CoralTrends_calcPercent(LIVE_CORAL)) %>%
    group_by(A_SECTOR,SHELF,REEF_NAME,REEF_ID,REPORT_YEAR) %>%
    summarise(Cover=mean(Cover, na.rm=TRUE), Tows=length(unique(TOW_SEQ_NO)), Latitude=mean(REEF_LAT, na.rm=TRUE), Longitude=mean(REEF_LONG, na.rm=TRUE)) %>%
    ungroup()

  wch <- table(manta.sum.G$REEF_NAME)
  wch<-names(wch[wch>4])
                                        #length(wch)  # number of reefs
  manta.sum.G = manta.sum.G %>% filter(REEF_NAME %in% wch)
                                        #dim(manta.sum) # number of reef surveys
  manta.sum.G = manta.sum.G %>%
    CoralTrends_calc3ZoneLocations %>%
    mutate(Location=ifelse(Latitude>-11.8, 'Outside',as.character(Location)))
  ## Reefs exclude (Outside)
  exc.G=manta.sum.G %>% filter(Location=='Outside')

  manta.sum.G = manta.sum.G %>% filter(Location!='Outside')
  nrow(manta.sum.G)
  length(unique(manta.sum.G$REEF_NAME))
}
