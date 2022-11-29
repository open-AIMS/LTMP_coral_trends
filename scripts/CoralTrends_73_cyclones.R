source('CoralTrends_functions.R') 
CoralTrends_checkPackages()

#load('data/manta.sum.newzones.RData')
load('../data/processed/manta.sum.RData')
## Generate a list of reefs we are using to help accumulate other sources of associated data
all.reefs = manta.sum %>%
    dplyr:::select(REEF_NAME,REEF_ID,Latitude,Longitude,Zone) %>%
    group_by(REEF_NAME,REEF_ID,Zone) %>%
    summarize_at(vars(Latitude,Longitude), funs(mean)) %>%
    as.data.frame

## Read in the 2016 (CoE) data
#cyclones = read.table('data/primary/S_by_reef.txt', header=TRUE, sep='\t')
#cyclones.all = read.table('data/primary/Disturb_cyclones.txt', header=TRUE, sep='\t')
#cyclones.all = read.csv('data/primary/170905 Cyclone data from Marji.csv', strip.white=TRUE)
#cyclones.all = read.csv('data/primary/170909 Cyclone wave data from Marji.csv', strip.white=TRUE)
## cyclones.all = read.csv('data/primary/20200330 Cyclone wave data from Marji.csv', strip.white=TRUE) %>% dplyr::select(-max)
cyclones.all = read.csv('../data/primary/20221114 Cyclone wave data from Marji.csv', strip.white=TRUE) %>%
    dplyr::select(-max)

cyclones = cyclones.all %>%
    dplyr:::rename(REEF_NAME=Reef,A_SECTOR=TMP_sector,SHELF=Shelf,REEF_LAT=lat, REEF_LONG=long_) %>%
    filter(!REEF_NAME=='') %>%                                                   #remove cases with missing reef name
    dplyr:::select(-Project,-gbrmpa_sector,-full_reef_id,-gazetted_name,-GBRMPA_ID) %>% #remote new and unnecessary fields
    gather(key=REPORT_YEAR, value=S, -A_SECTOR:-REEF_LONG) %>%                    #melt the data by year columns
    mutate(REPORT_YEAR=as.numeric(as.character(gsub('X','',REPORT_YEAR)))) %>%   #generate a REPORT_YEAR field that is a numeric year
    mutate(CYCLONEcat = S)                                                       #generate a CYCLONEcat field that has the cyclone category


##missing reefs
all.reefs  %>% full_join(cyclones %>% dplyr:::select(-A_SECTOR,-SHELF)) %>%
    filter(is.na(CYCLONEcat)) %>% dplyr:::select(REEF_NAME) %>% distinct
save(all.reefs, file='../data/processed/all.reefs.cyclones.RData')

cyclones = all.reefs  %>%
    left_join(cyclones %>% dplyr:::select(-A_SECTOR,-SHELF)) %>%
    filter(!is.na(CYCLONEcat))

cyclones.full = cyclones %>% dplyr:::select(REEF_ID, REEF_NAME,REPORT_YEAR,CYCLONEcat)
save(cyclones.full, file='../data/modelled/cyclones.full_3Zone.RData')

cyclones$Location <- factor(cyclones$Zone, levels=c("Northern GBR","Central GBR","Southern GBR"),
                                labels=c("Northern GBR","Central GBR","Southern GBR"))

## Fill in the gaps
cyclones.lookup = expand.grid(Location=c("Northern GBR","Central GBR","Southern GBR"),
                              REPORT_YEAR=seq.int(1985,2022,by=1))
cyclones %>% full_join(cyclones.lookup) %>% filter(is.na(Zone)) %>% head

cyclones = cyclones %>% full_join(cyclones.lookup) %>% arrange(REEF_NAME,REPORT_YEAR)  

cyclones = cyclones %>% mutate(Zone=ifelse(Location=='Northern GBR', 'Northern GBR',
                                    ifelse(Location=='Central GBR', 'Central GBR',
                                    ifelse(Location=='Southern GBR', 'Southern GBR', 'Great Barrier Reef'))))

## 2022 ======================================================
## Mike would like calculations of the duration between cyclones per decade

## The intervals could be calculated either by:
## - looking forwards from the focal decade or
## - looking backwards from teh focal decade

## The first iteration of this analysis looked forward.  The problem
## with this approach is that outcomes in later years are going to be
## biased towards having few reefs with intervals On the other hand,
## looking backwards will have a similar effect on early records.

## ---- Look forward
cyclones.interval <- cyclones %>%
    #filter(REEF_ID == 15005) %>%
    filter(CYCLONEcat > 0) %>%
    group_by(REEF_NAME) %>% 
    mutate(lead = lead(REPORT_YEAR),
           lag = lag(REPORT_YEAR),
           Interval = lead - REPORT_YEAR,
           Interval.lag = REPORT_YEAR - lag,
           Decade = factor(case_when(between(REPORT_YEAR, 1981,1990) ~ 1980,
                              between(REPORT_YEAR, 1991,2000) ~ 1990,
                              between(REPORT_YEAR, 2001,2010) ~ 2010,
                              between(REPORT_YEAR, 2011,2020) ~ 2020,
                              between(REPORT_YEAR, 2021,2030) ~ 2030)),
           ) %>%
#    group_by(Location, Decade, REEF_NAME) %>%
#    summarise(Interval = mean(Diff, na.rm = TRUE)) %>%
    ungroup()

ggplot(cyclones.interval) +
    geom_point(aes(y=Interval, x=REEF_NAME)) +
    facet_grid(Decade~ Location, scales='free_x')
ggplot(cyclones.interval) +
    geom_point(aes(y=Interval.lag, x=REEF_NAME)) +
    facet_grid(Decade~ Location, scales='free_x')

cyclones.interval %>%
   group_by(Location, Decade, REEF_NAME) %>%
   summarise(Interval = mean(Interval, na.rm = TRUE)) %>%
    group_by(Decade, Location) %>%
    summarise(Interval = mean(Interval, na.rm=TRUE))
cyclones.interval %>%
   group_by(Location, Decade, REEF_NAME) %>%
   summarise(Interval = mean(Interval.lag, na.rm = TRUE)) %>%
    group_by(Decade, Location) %>%
    summarise(Interval = mean(Interval, na.rm=TRUE))

cyclones.interval <- cyclones.interval %>%
    mutate(Decade_Location = interaction(Decade, Location)) %>%
    filter(!is.na(Interval))
    ## filter(!is.na(Interval.lag)) %>%
    ## dplyr::select(-Interval) %>%  
    ## rename(Interval = Interval.lag)


library(brms)
cyclones.interval.brm<-brm(Interval~0 + Decade_Location+(Decade_Location|REEF_NAME),
                                        #family=poisson(link='log'),
                           family=negbinomial(link = "log", link_shape = "log"),
                           data=cyclones.interval,
                           chains=3,cores = 3,iter=2000,warmup=500,thin=2, seed = 1,
                           backend = 'cmdstanr')   
save(cyclones.interval.brm, file = "../data/modelled/cyclones.interval.brm")
load(file = "../data/modelled/cyclones.interval.brm")
summary(cyclones.interval.brm)
emmeans(cyclones.interval.brm, ~Decade_Location, type='response')
em <- emmeans(cyclones.interval.brm, ~Decade_Location, type='response') %>% as.data.frame()
em %>% filter(!Decade_Location %in% c('2030.Northern GBR', '2030.Central GBR', '2030.Southern GBR')) %>%
    ggplot(aes(y = prob, x = Decade_Location)) +
    geom_pointrange(aes(ymin = lower.HPD, ymax = upper.HPD))
## ggemmeans(cyclones.interval.brm, ~Decade_Location) %>% plot()


cyclones.interval.brm1<-brm(Interval~0 + Decade+(1|REEF_NAME),
                                        #family=poisson(link='log'),
                           family=negbinomial(link = "log", link_shape = "log"),
                           data=cyclones.interval,# %>% filter(Decade != 2030),
                           chains=3,cores = 3,iter=2000,warmup=500,thin=2, seed = 1,
                           backend = 'cmdstanr')   
save(cyclones.interval.brm1, file = "../data/modelled/cyclones.interval.brm1")
load(file = "../data/modelled/cyclones.interval.brm1")
summary(cyclones.interval.brm1)
emmeans(cyclones.interval.brm1, ~Decade, type='response') %>% as.data.frame() %>%
    mutate(nDecade = as.numeric(as.character(Decade))) %>% 
    arrange(nDecade)

## ----end
## ---- Look backwards
cyclones.interval <- cyclones %>%
    #filter(REEF_ID == 15005) %>%
    filter(CYCLONEcat > 0) %>%
    group_by(REEF_NAME) %>% 
    mutate(lead = lead(REPORT_YEAR),
           lag = lag(REPORT_YEAR),
           Interval = lead - REPORT_YEAR,
           Interval.lag = REPORT_YEAR - lag,
           Decade = factor(case_when(between(REPORT_YEAR, 1981,1990) ~ 1980,
                              between(REPORT_YEAR, 1991,2000) ~ 1990,
                              between(REPORT_YEAR, 2001,2010) ~ 2010,
                              between(REPORT_YEAR, 2011,2020) ~ 2020,
                              between(REPORT_YEAR, 2021,2030) ~ 2030)),
           ) %>%
#    group_by(Location, Decade, REEF_NAME) %>%
#    summarise(Interval = mean(Diff, na.rm = TRUE)) %>%
    ungroup()

ggplot(cyclones.interval) +
    geom_point(aes(y=Interval, x=REEF_NAME)) +
    facet_grid(Decade~ Location, scales='free_x')
ggplot(cyclones.interval) +
    geom_point(aes(y=Interval.lag, x=REEF_NAME)) +
    facet_grid(Decade~ Location, scales='free_x')

cyclones.interval %>%
   group_by(Location, Decade, REEF_NAME) %>%
   summarise(Interval = mean(Interval, na.rm = TRUE)) %>%
    group_by(Decade, Location) %>%
    summarise(Interval = mean(Interval, na.rm=TRUE))
cyclones.interval %>%
   group_by(Location, Decade, REEF_NAME) %>%
   summarise(Interval = mean(Interval.lag, na.rm = TRUE)) %>%
    group_by(Decade, Location) %>%
    summarise(Interval = mean(Interval, na.rm=TRUE))

cyclones.interval <- cyclones.interval %>%
    mutate(Decade_Location = interaction(Decade, Location)) %>%
    ## filter(!is.na(Interval))
    filter(!is.na(Interval.lag)) %>%
    dplyr::select(-Interval) %>%  
    rename(Interval = Interval.lag)


library(brms)
cyclones.interval.brm<-brm(Interval~0 + Decade_Location+(Decade_Location|REEF_NAME),
                                        #family=poisson(link='log'),
                           family=negbinomial(link = "log", link_shape = "log"),
                           data=cyclones.interval,
                           chains=3,cores = 3,iter=2000,warmup=500,thin=2, seed = 1,
                           backend = 'cmdstanr')   
save(cyclones.interval.brm, file = "../data/modelled/cyclones.interval.back.brm")
load(file = "../data/modelled/cyclones.interval.back.brm")
summary(cyclones.interval.brm)
emmeans(cyclones.interval.brm, ~Decade_Location, type='response')
em <- emmeans(cyclones.interval.brm, ~Decade_Location, type='response') %>% as.data.frame()
em %>% filter(!Decade_Location %in% c('2030.Northern GBR', '2030.Central GBR', '2030.Southern GBR')) %>%
    ggplot(aes(y = prob, x = Decade_Location)) +
    geom_pointrange(aes(ymin = lower.HPD, ymax = upper.HPD))
## ggemmeans(cyclones.interval.brm, ~Decade_Location) %>% plot()


cyclones.interval.brm1<-brm(Interval~0 + Decade+(1|REEF_NAME),
                                        #family=poisson(link='log'),
                           family=negbinomial(link = "log", link_shape = "log"),
                           data=cyclones.interval,# %>% filter(Decade != 2030),
                           chains=3,cores = 3,iter=2000,warmup=500,thin=2, seed = 1,
                           backend = 'cmdstanr')   
save(cyclones.interval.brm1, file = "../data/modelled/cyclones.interval.back.brm1")
load(file = "../data/modelled/cyclones.interval.back.brm1")
summary(cyclones.interval.brm1)
emmeans(cyclones.interval.brm1, ~Decade, type='response') %>% as.data.frame() %>%
    mutate(nDecade = as.numeric(as.character(Decade))) %>% 
    arrange(nDecade)

## ----end


## ===========================================================

## Generate a summary that calculates the number and percentage of reefs in each zone per year
## that are impacted by each level of severity category
cyclones.sum = cyclones %>% group_by(REPORT_YEAR,Zone) %>%
    mutate(N=n()) %>%
    group_by(Zone,REPORT_YEAR,CYCLONEcat) %>%
    summarize(CYCLONE=n(), N=mean(N), CYCLONE.p=100*CYCLONE/N) %>%
    mutate(CYCLONE=ifelse(is.na(CYCLONEcat),NA,CYCLONE), N=ifelse(is.na(CYCLONEcat),NA,N),
           CYCLONE.p=ifelse(is.na(CYCLONEcat),0,CYCLONE.p))
## as above but only for all the whole GBR
cyclones.all = cyclones %>% ungroup %>% droplevels %>%
    group_by(REPORT_YEAR) %>%
    mutate(N=n()) %>%
    group_by(REPORT_YEAR,CYCLONEcat) %>%
    summarize(CYCLONE=n(), N=mean(N), CYCLONE.p=100*CYCLONE/N) %>%
    mutate(CYCLONE=ifelse(is.na(CYCLONEcat),NA,CYCLONE), N=ifelse(is.na(CYCLONEcat),NA,N),
           CYCLONE.p=ifelse(is.na(CYCLONEcat),0,CYCLONE.p)) %>% mutate(Zone='All')

cyclones.sum.all = rbind(cyclones.sum, cyclones.all)
cyclones.sum.all = cyclones.sum.all %>% ungroup %>% mutate(REPORT_YEAR=as.integer(REPORT_YEAR), CYCLONEcat = factor(CYCLONEcat))
cyclones.sum.all$Location <- factor(cyclones.sum.all$Zone, levels=c('All',"Northern GBR","Central GBR","Southern GBR"),
                                labels=c('Great Barrier Reef',"Northern GBR","Central GBR","Southern GBR"))

## Fill in the gaps
cyclones.sum.all = cyclones.sum.all %>% mutate(CYCLONEcat = ifelse(is.na(CYCLONEcat), '0', as.character(CYCLONEcat)), CYCLONE.p = ifelse(is.na(CYCLONE.p),0,CYCLONE.p)) %>%
    droplevels
save(cyclones.sum.all, file='../data/modelled/cyclones.sum.all_3Zone.RData')
load(file='../data/modelled/cyclones.sum.all_3Zone.RData')

labs = levels(cyclones.sum.all$Location)
#save(cyclones.sum.all, file='data/modelled/cods.sum.all.RData')
gcyclones = ggplot(cyclones.sum.all %>% filter(REPORT_YEAR >1985), aes(y=CYCLONE.p, x=REPORT_YEAR)) +
    geom_bar(stat='identity',position='stack',aes(alpha=CYCLONEcat), fill='red',show.legend=FALSE) +
    #geom_line(stat='identity',position='stack',aes(alpha=CYCLONEcat), fill='red',show.legend=FALSE) +
    #facet_grid(~Location) +
    #facet_wrap(~Location, nrow=1, scales='fixed', labeller=labeller(Location=setNames(paste0("", labs, "\n"), labs)))+
    facet_wrap(~Location, nrow=1, strip.position='bottom',scales='fixed', labeller=labeller(Location=setNames(paste0("", labs, "\n"), labs)))+
    scale_alpha_manual(breaks=c('0','1','2','3','4','5'), values=c(0,seq(0.4,1,len=5))) +
    scale_y_reverse(expression(Reefs~impacted~('%')),expand=c(0,0),lim=c(90,0)) +
    scale_x_continuous('',breaks=seq(1985,2015,by=5),limits=c(1985,2020),position = "top")+
    theme_classic() +
    theme(strip.background=element_rect(fill='#00009930', color='black', size=0.5),
          panel.border=element_rect(fill=NA,color='black'),
          axis.title.y=element_text(margin=margin(r=1,unit='lines')),
          panel.grid.minor=element_line(size=0.5,color='grey40',linetype='dotted'),
          panel.grid.major=element_line(color='grey40',linetype='dotted'),
          panel.grid.minor.x=element_line(size=0,color='white',linetype=NULL),
          panel.grid.major.x=element_line(size=0,color='white',linetype=NULL),
          plot.margin=unit(c(0,0,0,0),'pt'))
gcyclones
ggsave(file='../output/figures/cyclones_3Zone.pdf', gcyclones, width=15, height=3)
