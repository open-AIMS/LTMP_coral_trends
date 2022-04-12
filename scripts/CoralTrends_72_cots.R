library(tidyverse)

load('../data/processed/manta.sum.RData')

writeLines("
SELECT V_RM_SAMPLE.P_CODE, V_RM_SAMPLE.A_SECTOR, V_RM_SAMPLE.SHELF, V_RM_SAMPLE.REEF_ID, V_RM_SAMPLE.REEF_NAME, V_RM_SAMPLE.REEF_LAT,
V_RM_SAMPLE.REEF_LONG, V_RM_SAMPLE.OPENORCLOSED, V_RM_SAMPLE.OPENORCLOSED_AFTER2004, V_RM_SAMPLE.REPORT_YEAR, V_RM_SAMPLE.VISIT_NO,
Avg(RM_MEDIAN.MEAN_LIVE) AS AvgOfMEAN_LIVE, Avg(RM_MEDIAN.MEAN_COTS) AS AvgOfMEAN_COTS
FROM V_RM_SAMPLE INNER JOIN RM_MEDIAN ON V_RM_SAMPLE.SAMPLE_ID = RM_MEDIAN.SAMPLE_ID
GROUP BY V_RM_SAMPLE.P_CODE, V_RM_SAMPLE.A_SECTOR, V_RM_SAMPLE.SHELF, V_RM_SAMPLE.REEF_ID, V_RM_SAMPLE.REEF_NAME, V_RM_SAMPLE.REEF_LAT,
V_RM_SAMPLE.REEF_LONG, V_RM_SAMPLE.OPENORCLOSED, V_RM_SAMPLE.OPENORCLOSED_AFTER2004, V_RM_SAMPLE.REPORT_YEAR, V_RM_SAMPLE.VISIT_NO
HAVING (((V_RM_SAMPLE.P_CODE) Not Like 'TS' And (V_RM_SAMPLE.P_CODE) Not Like 'WA_NI'))
ORDER BY V_RM_SAMPLE.REEF_ID, V_RM_SAMPLE.REPORT_YEAR","../data/primary/cots.sql")

if (goto_database_manta) system("java -jar dbExport.jar ../data/primary/cots.sql ../data/primary/cots.csv reef reefmon")

cots <- read.csv('../data/primary/cots.csv',strip.white=TRUE) %>%
    filter(REPORT_YEAR < (finalYear+1)) %>% droplevels()
summary(cots)   




#blank = data.frame(COTScat = factor(c('Zero','NO','IO','AO'), levels=c('Zero','NO','IO','AO'),N=0))

cots.sum = cots %>% filter(REEF_NAME %in% as.character(unique(manta.sum$REEF_NAME))) %>%
    mutate(COTScat = COTScategories(AVGOFMEAN_COTS), COTScat = factor(COTScat, levels=c('Zero','NO','IO','AO'))) %>%
    left_join(manta.sum %>% dplyr:::select(REEF_NAME,Zone)) %>%
    group_by(Zone,REPORT_YEAR) %>%
    mutate(N=n()) %>%
    #group_by(A_SECTOR,SHELF,REEF_NAME,REEF_ID,REPORT_YEAR,COTScat) %>%
    group_by(Zone,REPORT_YEAR,COTScat) %>%
    summarize(COTS=n(), N=mean(N), COTS.p=100*COTS/N)

cots.full = cots %>% filter(REEF_NAME %in% as.character(unique(manta.sum$REEF_NAME))) %>%
    mutate(COTScat = COTScategories(AVGOFMEAN_COTS), COTScat = factor(COTScat, levels=c('Zero','NO','IO','AO'))) %>%
    left_join(manta.sum %>% dplyr:::select(REEF_NAME,Zone)) %>%
    dplyr:::select(REEF_ID, REEF_NAME,Zone, REPORT_YEAR,COTScat)
save(cots.full, file='../data/modelled/cots.full_3Zone.RData')

cots.all = cots %>% filter(REEF_NAME %in% as.character(unique(manta.sum$REEF_NAME))) %>%
    mutate(COTScat = COTScategories(AVGOFMEAN_COTS), COTScat = factor(COTScat, levels=c('Zero','NO','IO','AO'))) %>%
    left_join(manta.sum %>% dplyr:::select(REEF_NAME,Zone)) %>%
    group_by(REPORT_YEAR) %>%
    mutate(N=n()) %>%
    #group_by(A_SECTOR,SHELF,REEF_NAME,REEF_ID,REPORT_YEAR,COTScat) %>%
    group_by(REPORT_YEAR,COTScat) %>%
    summarize(COTS=n(), N=mean(N), COTS.p=100*COTS/N) %>% mutate(Zone='All')

cots.sum.all = rbind(cots.sum, cots.all)
cots.sum.all$Location <- factor(cots.sum.all$Zone, levels=c('All',"Northern GBR","Central GBR","Southern GBR"),
                                labels=c('Great Barrier Reef',"Northern GBR","Central GBR","Southern GBR"))

labs = levels(cots.sum.all$Location)
save(cots.sum.all, file='../data/modelled/cots.sum.all_3Zone.RData')
gcots = ggplot(cots.sum.all %>% filter(REPORT_YEAR >1985, COTScat %in% c('AO','IO')), aes(y=COTS.p, x=REPORT_YEAR)) +
    geom_bar(stat='identity',position='stack',aes(alpha=COTScat), fill='blue',show.legend=FALSE) +
    #geom_line(stat='identity',position='stack',aes(alpha=COTScat), fill='red',show.legend=FALSE) +
    #facet_grid(~Location) +
    #facet_wrap(~Location, nrow=1, scales='fixed', labeller=labeller(Location=setNames(paste0("", labs, "\n"), labs)))+
    facet_wrap(~Location, nrow=1, strip.position='bottom',scales='fixed', labeller=labeller(Location=setNames(paste0("", labs, "\n"), labs)))+
    scale_alpha_manual(breaks=c('IO','AO'), values=c(0.4,1)) +
    scale_y_reverse(expression(Reefs~impacted~('%')),expand=c(0,0),lim=c(90,0)) +
    scale_x_continuous('',breaks=seq(1985,2015,by=5),position = "top")+
    theme_classic() +
    theme(strip.background=element_rect(fill='#00009930', color='black', size=0.5),
          panel.border=element_rect(fill=NA,color='black'),
          axis.title.y=element_text(margin=margin(r=1,unit='lines')),
          panel.grid.minor=element_line(size=0.5,color='grey40',linetype='dotted'),
          panel.grid.major=element_line(color='grey40',linetype='dotted'),
          panel.grid.minor.x=element_line(size=0,color='white',linetype=NULL),
          panel.grid.major.x=element_line(size=0,color='white',linetype=NULL),
          plot.margin=unit(c(0,0,0,0),'pt'))
gcots
## ggsave(file='output/figures/cots_3Zone.pdf', gcots, height=5, width=15, units='in')

gcots = ggplot(cots.sum.all %>% filter(REPORT_YEAR >1985, COTScat %in% c('AO','IO')), aes(y=COTS.p, x=REPORT_YEAR)) +
    geom_col(position='stack',aes(x=REPORT_YEAR-0.4,alpha=COTScat), width=0.2,fill='red',show.legend=FALSE) +
    #geom_col(bleaching.sum.all %>% filter(REPORT_YEAR >1985, Bleachingcat %in% c('1','2', '3', '4')), position='stack',aes(x=REPORT_YEAR-0,alpha=COTScat), width=0.2,fill='blue',show.legend=FALSE) +
    facet_grid(~Location) +
    #facet_wrap(~Location, nrow=1, scales='fixed', labeller=labeller(Location=setNames(paste0("", levels(newdata.newzones$Location), "\n"), levels(newdata.newzones$Location))))+
    scale_alpha_manual(breaks=c('IO','AO'), values=c(0.1,1)) +
    scale_y_reverse(expression(Reefs~impacted~('%')),expand=c(0,0),lim=c(90,0)) +
    scale_x_continuous('',breaks=seq(1985,2015,by=5),position = "top")+
    theme_classic() +
    theme(strip.background=element_rect(fill='#00009930', color='black', size=0.5),
          panel.border=element_rect(fill=NA,color='black'),
          axis.title.y=element_text(margin=margin(r=1,unit='lines')),
          panel.grid.minor=element_line(size=0.5,color='grey40',linetype='dotted'),
          panel.grid.major=element_line(color='grey40',linetype='dotted'),
          panel.grid.minor.x=element_line(size=0,color='white',linetype=NULL),
          panel.grid.major.x=element_line(size=0,color='white',linetype=NULL),
          plot.margin=unit(c(0,0,0,0),'pt'),
          strip.text=element_text(size=15,lineheight=1.0, face='bold',hjust=0.95))
gcots
## ggsave(file='output/figures/cots1_3Zone.pdf', gcots, height=5, width=15, units='in')


