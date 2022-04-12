source('CoralTrends_functions.R') 
CoralTrends_checkPackages()

load('../data/processed/manta.sum.RData')

lookup.reefs <- manta.sum %>%
    dplyr::select(REEF_ID,REEF_NAME,Latitude,Longitude,Zone) %>%
    distinct

## Read in all the bleaching data.  As of 2022, Mike is supplying pre-processed data!

## 0 = 0
## 1 = >0 - 0.1
## 2 = >0.1 - 0.3
## 3 = >0.3 - 0.6
## 4 = >0.6 - 0.9
## 5 = >0.9
## severe 4 and 5 (or 3, 4, 5)
## - all, filter LTMP dec to march
## - filter out 2022

## ---- Bleaching from Mike (LTMP + COE)
bleaching <- read_csv('../data/primary/220405 combined bleaching scores.csv') %>%
    ## The following forces the listed reefs to be in central (which is what their lat/long would do)
    ## the data do not have lat/long, Mike already put them into Regions - however his method was not
    ## consistent with De'ath
    mutate(Region = ifelse(Reef %in% c('BOULDER REEF','15077S','EDGELL REEFS (NO 5)',
                                       '15047S', 'EGRET REEF', 'IRENE REEF', 'LENA REEF',
                                       'RIBBON NO 1 REEF', 'RIBBON NO 3 REEF',
                                       'ROSSER REEF'),
                           'Central', Region)) %>%
    mutate(Region = ifelse(Full.Reef.ID %in% c('20113S'), 'Central', Region)) %>%
    mutate(Zone = Region,
           Location = Region,
           Severe = ifelse(mid.point >= 0.6, 1, 0),
           BleachingCAT = COE_COTScategories(mid.point)) %>%
    filter(!is.na(mid.point),
           REPORT_YEAR < (finalYear + 1)) %>% droplevels() %>%
    group_by(Zone, REPORT_YEAR, Reef, Full.Reef.ID) %>%
    mutate(n = n()) %>%
    mutate(flag = ifelse(any(Project == 'COE'), 1,0),
           flag2 = ifelse(Project == 'LTMP', 1, 0),
           flag3 = ifelse(n>1, 1, 0),
           flag4 = ifelse(flag & flag2 & flag3, 1, 0),
           Date = as.Date(SAMPLE_DATE, '%d-%b-%y'),
           Season = ifelse(months(Date) %in% c('December','January','February','March'), 'Bleaching', 'Non-bleaching')) %>%
    filter(flag4==0) %>%
    dplyr::select(-starts_with('flag'))
## ----end
## ---- Bleaching from Ray B
## Read in the 1998 and 2002 (RB) data
bleaching.1998 = read.xls('../data/primary/2016 Bleaching Response AIMS Site Selection_1998 2002 Repeated Aerial Su....xlsx', sheet=4, header=TRUE)
bleaching.1998 = bleaching.1998 %>% dplyr:::select(REEF_ID, X2002.GBR.Name, BLEACHING_1998, BLEACHING_2002,Latitude=Y_COORD,Longitude=X_COORD) %>% 
    filter(!is.na(REEF_ID)) %>%
    filter(as.character(REEF_ID)!='') %>%
    mutate(BLEACHING_1998 = -1*(BLEACHING_1998 - 5),
           BLEACHING_2002 = -1*(BLEACHING_2002 - 5)) %>%
    gather(key=REPORT_YEAR, value=BLEACHINGcat,BLEACHING_1998,BLEACHING_2002) %>%
    mutate(REPORT_YEAR = gsub('BLEACHING_','',REPORT_YEAR)) %>%
    filter(!is.na(Latitude)) %>% droplevels %>%
    mutate(Severe = ifelse(BLEACHINGcat > 2, 1, 0)) 
## Put the reefs into Locations (Zones)
#bleaching.1998 = CoralTrends_calc3ZoneLocations.df(bleaching.1998)
## bleaching.1998 = CoralTrends_calc3ZoneLocations(bleaching.1998)
bleaching.1998 = CoralTrends_calc3ZoneLocation(bleaching.1998) %>%
    mutate(Region = gsub(' GBR', '', Region))
bleaching.1998 = bleaching.1998 %>% filter(!is.na(BLEACHINGcat)) %>%
    mutate(REEF_ID=stringr::str_sub(as.character(REEF_ID), start=1, end=5),
           BLEACHINGcat=as.numeric(BLEACHINGcat),
           REPORT_YEAR=as.numeric(as.character(REPORT_YEAR))) %>%
    group_by(Region, REEF_ID,REPORT_YEAR) %>%
    summarize(RAYcat=max(BLEACHINGcat)) %>%
    filter(REEF_ID %in% unique(lookup.reefs$REEF_ID)) %>%
    mutate(Zone = Region, Location = Region) %>%
    left_join(lookup.reefs %>%
              mutate(REEF_ID = as.character(REEF_ID)) %>%
              dplyr::select(REEF_ID, Reef = REEF_NAME))

## ----end
## ---- Combine all
bleaching <-
    bleaching %>%
    mutate(REEF_ID=stringr::str_sub(as.character(Full.Reef.ID), start=1, end=5)) %>% 
    full_join(bleaching.1998) %>%   
    mutate(BleachingCAT=ifelse(!is.na(BleachingCAT),BleachingCAT,
                        ifelse(!is.na(RAYcat),RAYcat,NA))) %>%
    mutate(Project = ifelse(is.na(Project), 'RB', Project)) %>%
    mutate(Reef = ifelse(is.na(Reef), REEF_ID, Reef))
## ----end


## ---- Bleaching full
bleaching.full_3Zone <- bleaching %>%
    ungroup() %>%
    dplyr::select(REEF_ID, # = Full.Reef.ID,
                  REPORT_YEAR,
                  REEF_NAME = Reef,
                  BLEACHINGcat = BleachingCAT,
                  Zone,
                  Date,
                  Season,
                  Project)
save(bleaching.full_3Zone, file='../data/modelled/bleaching.full_3Zone.RData')
## ----end

## ---- Bleaching Sum
bleaching.sum.all <-
    bleaching %>%
    group_by(Zone, REPORT_YEAR, BleachingCAT) %>%
    summarise(BLEACHING = n(),
              BLEACHING.type2 = sum((Project == 'LTMP' & Season == 'Bleaching') | Project %in% c('COE', 'RB')), na.rm=TRUE,
              BLEACHING.type2 = ifelse(is.na(BLEACHING.type2), 0, BLEACHING.type2)) %>%
    group_by(Zone, REPORT_YEAR) %>%
    mutate(N = sum(BLEACHING),
           BLEACHING.p = BLEACHING/N,
           N.type2 = sum(BLEACHING.type2),
           BLEACHING.type2.p = BLEACHING.type2/N.type2,
           BLEACHING.type2.p = ifelse(is.na(BLEACHING.type2.p), 0, BLEACHING.type2.p))%>%
    mutate(Zone = factor(Zone, levels = c('Northern', 'Central', 'Southern'),
                         labels = c('Northern GBR', 'Central GBR', 'Southern GBR'))) %>%
    arrange(Zone, REPORT_YEAR, BleachingCAT) %>%
    dplyr::select(Zone,REPORT_YEAR, BLEACHINGcat = BleachingCAT,
                  BLEACHING, N, BLEACHING.p,
                  BLEACHING.type2, N.type2, BLEACHING.type2.p) %>%
    mutate(Location = Zone) %>%
    mutate(BLEACHINGcat = factor(BLEACHINGcat),
           BLEACHING.type2 = factor(BLEACHING.type2))
    ## filter(Severe == 1)

## bleaching.sum.all %>% as.data.frame %>% head 

## bleaching %>% as.data.frame %>% head

save(bleaching.sum.all, file='../data/modelled/bleaching.sum.all_3Zone.RData')    
## load(file='../../new/data/modelled/bleaching.sum.all_3Zone.RData')    
## bleaching.sum.all
## ----end


## ---- tests
bleaching.full_3Zone %>%
    filter(REPORT_YEAR == 1998,
           Zone=='Northern') %>%
    as.data.frame() %>%
    filter(Project == 'RB') %>%
    arrange(REEF_ID)


bleaching.sum.all %>%
    filter(REPORT_YEAR == 1995) %>%
    as.data.frame
## ----end



## ## Read in the 2016 (CoE) data
## #bleaching.2016 = read.xls('data/primary/170627 CoECRS_2016AerialScores_LTMP.xlsx', sheet=1, header=TRUE)
## ## bleaching.2016_2017 = read.xls('data/primary/CoECRS_2016_2017AerialScores_LTMP.xlsx', sheet=1, header=TRUE)
## bleaching.2016_2017_2020 = read.xls('data/primary/CoECRS_2016_2017_2020AerialScores_LTMP.xlsx', sheet=1, header=TRUE)

## ## We could attempt to match up Reefs to LTMP reefs, but this might just restrict
## ## Instead, lets use all reefs and just assign them a Zone.

## #bleaching.2016 = bleaching.2016 %>% filter(Bleaching.cat!='#N/A') %>% droplevels %>%
## #  dplyr:::select(Full.Reef.ID,REEF_NAME=Reef, BLEACHINGcat=Bleaching.cat,Latitude=Reef.Lat, Longitude=Reef.Long) %>% mutate(REPORT_YEAR=2016) 

## bleaching.2016 = bleaching.2016_2017_2020 %>%
##   gather(key=Year, value=Bleaching.cat, Bleaching.cat.2016, Bleaching.cat.2017, Bleaching.cat.2020) %>%
##   separate(Year,c('A','REPORT_YEAR'), sep='Bleaching.cat.') %>%
##   mutate(REPORT_YEAR=as.numeric(as.character(REPORT_YEAR))) %>%
##   dplyr::select(-A) %>% 
##   filter(Bleaching.cat!='#N/A') %>% droplevels %>%
##   dplyr:::select(Full.Reef.ID,REEF_NAME=Reef, BLEACHINGcat=Bleaching.cat,Latitude=Reef.Lat, Longitude=Reef.Long,REPORT_YEAR)


## ## Put the reefs into Locations (Zones)
## bleaching.2016 = CoralTrends_calc3ZoneLocations.df(bleaching.2016)
## save(bleaching.2016, file='data/processed/bleaching.2016_3Zone.RData')

## ## Read in the 1998 and 2002 (RB) data
## bleaching.1998 = read.xls('data/primary/2016 Bleaching Response AIMS Site Selection_1998 2002 Repeated Aerial Su....xlsx', sheet=4, header=TRUE)
## bleaching.1998 = bleaching.1998 %>% dplyr:::select(REEF_ID, X2002.GBR.Name, BLEACHING_1998, BLEACHING_2002,Latitude=Y_COORD,Longitude=X_COORD) %>% 
##     filter(!is.na(REEF_ID)) %>% filter(as.character(REEF_ID)!='') %>%
##     mutate(BLEACHING_1998 = -1*(BLEACHING_1998 - 5),BLEACHING_2002 = -1*(BLEACHING_2002 - 5)) %>%
##     gather(key=REPORT_YEAR, value=BLEACHINGcat,BLEACHING_1998,BLEACHING_2002) %>%
##     mutate(REPORT_YEAR = gsub('BLEACHING_','',REPORT_YEAR)) %>%
##     filter(!is.na(Latitude)) %>% droplevels
## ## Put the reefs into Locations (Zones)
## #bleaching.1998 = CoralTrends_calc3ZoneLocations.df(bleaching.1998)
## bleaching.1998 = CoralTrends_calc3ZoneLocations(bleaching.1998)
## bleaching.1998 = bleaching.1998 %>% filter(!is.na(BLEACHINGcat))


## ##Create a merged version for frequency analysis
## ## It turns out that although Rays data has a field called REEF_ID, this is
## ## probably FULL_REEF_ID as it contains more than 5 digits
## ## We need to pair back these ID's to just 5 digits

## ## Start by creating a master lookp of reefs to include from each source.
## lookup.reefs = manta.sum %>% dplyr::select(REEF_ID,REEF_NAME,Latitude,Longitude,Zone) %>% distinct

## #load(file='data/processed/all.reefs.cyclones.RData')
## #lookup.reefs = all.reefs %>% dplyr::select(REEF_ID,REEF_NAME,Latitude,Longitude,Zone) %>% distinct


## ## Now filter each of the sources to the reefs in this lookup.
## ## In so doing, we must also only use REEF_ID (not FULLREEF.ID)
## bleaching.allyears = read.csv('data/primary/bleaching manta aesthetics.csv', strip.white=TRUE)
## bleaching.allyears = bleaching.allyears %>%
##     mutate(Tmp.score = ifelse(MaxOfBLEACHING_PERHC=='', as.character(MaxOfBLEACHING),
##                        ifelse(MaxOfBLEACHING_PERHC=='0', 0,
##                        ifelse(MaxOfBLEACHING_PERHC=='0+', 0,
##                        ifelse(MaxOfBLEACHING_PERHC=='-1',1,
##                        ifelse(MaxOfBLEACHING_PERHC=='1',1,as.character(MaxOfBLEACHING_PERHC))))))) %>%
##     mutate(aerial.score = ifelse(is.na(aerial.score),
##                                  as.numeric(as.character(Tmp.score)),
##                                  aerial.score)) %>%
##     dplyr::select(-Tmp.score)
## bleaching.allyears = bleaching.allyears %>%
##     #mutate(score = as.numeric(as.character(ifelse(MaxOfBLEACHING_PERHC=='', MaxOfBLEACHING,
##     #                                       ifelse(MaxOfBLEACHING_PERHC=='0', 0,
##     #                                       ifelse(MaxOfBLEACHING_PERHC=='0+', 0,
##     #                                       ifelse(MaxOfBLEACHING_PERHC=='-1',1,
##     #                                       ifelse(MaxOfBLEACHING_PERHC=='1',1,MaxOfBLEACHING_PERHC)))))))) %>%
##     mutate(REEF_ID=stringr::str_sub(as.character(FULLREEF_ID), start=1, end=5)) %>%
##     group_by(REEF_ID,REPORT_YEAR) %>%
##     summarize(REEF_NAME=first(REEF_NAME), MANTAcat=max(aerial.score)) %>% #, MANTAcat1=max(score)) %>%
##     filter(REEF_ID %in% unique(lookup.reefs$REEF_ID))
## dim(bleaching.allyears)

## ## Now same for Terry's data
## bleaching.2016 = bleaching.2016 %>%
##     mutate(REEF_ID=stringr::str_sub(as.character(Full.Reef.ID), start=1, end=5), BLEACHINGcat=as.numeric(as.character(BLEACHINGcat))) %>%
##     group_by(REEF_ID,REPORT_YEAR) %>%
##     summarize(REEF_NAME=first(REEF_NAME), TERRYcat=max(BLEACHINGcat)) %>%
##     filter(REEF_ID %in% unique(lookup.reefs$REEF_ID))
## dim(bleaching.2016)

## ## Now same for Ray's data
## bleaching.1998 = bleaching.1998 %>%
##     mutate(REEF_ID=stringr::str_sub(as.character(REEF_ID), start=1, end=5), BLEACHINGcat=as.numeric(BLEACHINGcat), REPORT_YEAR=as.numeric(as.character(REPORT_YEAR))) %>%
##     group_by(REEF_ID,REPORT_YEAR) %>%
##     summarize(RAYcat=max(BLEACHINGcat)) %>%
##     filter(REEF_ID %in% unique(lookup.reefs$REEF_ID))
## bleaching.1998 %>% dim
## dim(bleaching.1998)

## ## Merge them all together
## bleaching.merge=bleaching.allyears %>% full_join(bleaching.2016) %>% full_join(bleaching.1998) %>%
##     mutate(BLEACHINGcat=ifelse(!is.na(TERRYcat),TERRYcat, ifelse(!is.na(RAYcat),RAYcat,MANTAcat))) %>%
##     dplyr::select(-TERRYcat,-MANTAcat,-RAYcat) %>%
##     group_by(REEF_ID) %>% arrange(REPORT_YEAR) %>% ungroup %>%
##     full_join(lookup.reefs %>% mutate(REEF_ID=as.character(REEF_ID)))
## dim(bleaching.merge)
## save(bleaching.merge, file='data/modelled/bleaching.merge_3Zone.RData')    


## ## Collate the number of reefs of each category/year/Location
## bleaching.sum = bleaching.merge %>%
##     group_by(Zone,REPORT_YEAR) %>%
##     mutate(N=n()) %>%
##     group_by(Zone,REPORT_YEAR,BLEACHINGcat) %>%
##     summarize(BLEACHING=n(), N=mean(N), BLEACHING.p=100*BLEACHING/N) %>%
##     mutate(BLEACHING=ifelse(is.na(BLEACHINGcat),NA,BLEACHING), N=ifelse(is.na(BLEACHINGcat),NA,N),
##            BLEACHING.p=ifelse(is.na(BLEACHINGcat),0,BLEACHING.p))

## bleaching.all = bleaching.merge %>% ungroup %>% droplevels %>%
##     group_by(REPORT_YEAR) %>%
##     mutate(N=n()) %>%
##     group_by(REPORT_YEAR,BLEACHINGcat) %>%
##     summarize(BLEACHING=n(), N=mean(N), BLEACHING.p=100*BLEACHING/N) %>%
##     mutate(BLEACHING=ifelse(is.na(BLEACHINGcat),NA,BLEACHING), N=ifelse(is.na(BLEACHINGcat),NA,N),
##            BLEACHING.p=ifelse(is.na(BLEACHINGcat),0,BLEACHING.p)) %>% mutate(Zone='All')

## bleaching.sum.all = rbind(bleaching.sum, bleaching.all)
## bleaching.sum.all = bleaching.sum.all %>% ungroup %>% mutate(REPORT_YEAR=as.integer(REPORT_YEAR), BLEACHINGcat = factor(BLEACHINGcat))
## bleaching.sum.all$Location <- factor(bleaching.sum.all$Zone, levels=c('All',"Northern","Central","Southern"),
##                                 labels=c('Great Barrier Reef',"Northern","Central","Southern"))

## bleaching.lookup = expand.grid(Location=c('Great Barrier Reef',"Northern","Central", "Southern"),
##                                 REPORT_YEAR=seq.int(1985,2020,by=1))
## #bleaching.sum.all = bleaching.sum.all %>% right_join(bleaching.lookup)
## bleaching.sum.all = bleaching.sum.all %>% full_join(bleaching.lookup)                                
## bleaching.sum.all = bleaching.sum.all %>% mutate(Zone=Location, BLEACHINGcat = ifelse(is.na(BLEACHINGcat), 0, as.character(BLEACHINGcat)), BLEACHING.p = ifelse(is.na(BLEACHING.p),0,BLEACHING.p)) %>%
##     droplevels
## save(bleaching.sum.all, file='data/modelled/bleaching.sum.all_3Zone.RData')    



## ## stop here - dont run any further!!

