#setwd("C:/Users/memslie/OneDrive - Australian Institute of Marine Science/Desktop/working folder/GBR coral recovery 2017/Ecology/for submission/2022 data added")

library(plyr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(INLA)
library(brms)
library(broom)
library(rstan)
library(emmeans)
library(tidybayes)
library(DHARMa)
library(grid)
library(gridExtra)


###data wrangling
manual<-read.csv(file='../data/primary/manual recovery indiv manta reefs_2022.csv',strip.white=T)


manual<-manual %>% mutate(abs.decline.bench=ifelse(Disturb.year<=Bench.year,NA,bench.cover-Disturb.cover),
                          rel.decline.bench=ifelse(Disturb.year<=Bench.year,NA,abs.decline.bench/bench.cover*100),
                          abs.decline.prior=prior.cover-Disturb.cover,
                          rel.decline.prior=abs.decline.prior/prior.cover*100,
                          rel.perc.benchmark.recovery=ifelse(recovery.year==Bench.year,NA,recovery.cover/bench.cover*100),
                          rel.perc.prior.recovery=ifelse(recovery.year==prior.year,NA,recovery.cover/prior.cover*100),
                          years.benchmark.recovery=(recovery.year-Bench.year),
                          years.prior.recovery=recovery.year-prior.year,
                          recovery.interval=recovery.year-Disturb.year,
                          recovery.rate=(recovery.cover-Disturb.cover)/(recovery.year-Disturb.year),
                          diff.2.baseline=recovery.cover-bench.cover) %>% 
  mutate(Region=factor(Region,levels=c("Northern","Central","Southern")),
         Reef_name=factor(Reef_name),
         Disturbance=factor(Disturbance,levels=c("Bleaching","COTS","Cyclone","Multiple","Unknown")),
         Disturbance2=factor(Disturbance2,levels=c("Bleaching","COTS","Cumulative","Cyclone","Unknown")),
         Dist1=factor(Dist1,levels=c("Bleaching","COTS","Cyclone","Unknown","Multiple")))

manual<-manual %>% dplyr::mutate(recovered.prior=ifelse(rel.perc.prior.recovery>=100,'yes','no'),
                          recovered.bench=ifelse(rel.perc.benchmark.recovery>=100,'yes','no'))

head(manual)
str(manual)

#########################################################################################
#### Figure 2 ------------------------------------------------------

#########
####relative decline

## sensible priors
manual %>% group_by(Region, Disturbance) %>%
  summarise(median(log(rel.decline.prior)),
            mad(log(rel.decline.prior)))
priors<-prior("normal(4, 0.5)", class = 'Intercept') +
  prior("normal(0,2)", class = 'b')


###################################################################################
#### have disturbances become more severe through time -----------------------

plot3<-ggplot(manual,aes(x=Disturb.year,y=Disturb.cover))+
  geom_point(shape=16,alpha=0.6,stat="summary", fun.y = "mean")+
  geom_errorbar(stat = "summary", fun.data = mean_se,width=0.1)+
  scale_y_continuous("Post disturbance coral cover (%)",breaks=c(0,0.1,0.2,0.3,0.4,0.5),labels = c(0,10,20,30,40,50))+
  scale_x_continuous('Year of disturbance')+
  scale_size_area(name="Pre-disturbance\n cover")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=0))
#annotate("text",x=15,y=420,label='25% of reefs recovered disturbances not benchmark')

plot3

ggsave(plot3,file="disturbance cover.pdf",height=6,width=10)


####
##models


##brm model - did try a polynomial fit which didn't improve residuals
post.dist.cover.brm<-brm(bf(Disturb.cover~Disturb.year+(1|Disturbance)+(1|Reef_name),
                            #shape~Disturb.year+(1|Disturbance)+(1|Reef_name),  
                            family=Gamma(link='log')),
                         data=manual,
                         chains=3,iter=2000,warmup=500,thin=2)

save(post.dist.cover.brm,file='post.dist.cover.brm.RData')
summary(post.dist.cover.brm)

load(file='post.dist.cover.brm.RData')

###sampling diagnostics

stan_trace(post.dist.cover.brm$fit)
stan_dens(post.dist.cover.brm$fit)

##### residuals
preds <- posterior_predict(post.dist.cover.brm,  ndraws=250,  summary=FALSE)
mod.resids <- createDHARMa(simulatedResponse = t(preds),
                           observedResponse = manual$Disturb.cover,
                           fittedPredictedResponse = apply(preds, 2, median),
                           integerResponse = FALSE)
plot(mod.resids)

marginal_effects(post.dist.cover.brm)

save(post.dist.cover.brm,file='post.dist.cover.brm.RData')
load(file='post.dist.cover.brm.RData')


###extracting for linear
mcmc=post.dist.cover.brm %>% as.matrix()
post.dist.cover.slope<-median_hdci(exp(mcmc[,2]))

### calculating 'slope' (percentage change +/- CIs)
(1-post.dist.cover.slope[1])*100
(1-post.dist.cover.slope[2])*100
(1-post.dist.cover.slope[3])*100
### exceedence probabilty
sum(mcmc[,2]<0)/length(mcmc[,2])



## extracting trend for plot

predgrid<-with(manual,seq(min(Disturb.year),max(Disturb.year),length=1000))

post.cover.fitted.trend<-emmeans(post.dist.cover.brm,~Disturb.year,
                                 at=list(Disturb.year=predgrid),type='response') %>%  as.data.frame()
plot3<-plot3+
  geom_ribbon(data=post.cover.fitted.trend,aes(y=response,x=Disturb.year,ymin=lower.HPD,ymax=upper.HPD),alpha=0.2)+
  geom_line(data=post.cover.fitted.trend,aes(y=response,x=Disturb.year),colour='blue')+
  ggtitle('a)')

plot3


########################
#relative coral loss

rel.coral.loss.thru.time.plot<-ggplot(manual,aes(x=Disturb.year,y=rel.decline.prior))+
  geom_point(shape=16,alpha=0.6,stat="summary", fun.y = "mean")+
  geom_errorbar(stat = "summary", fun.data = mean_se,width=0.1)+
  scale_x_continuous('Year of disturbance')+
  scale_y_continuous('Relative coral loss (%)')+
  ggtitle('b)')+
  theme_classic()

rel.coral.loss.thru.time.plot

### models

#############
#### no support for polynomial terms higher than 1 degree therefore use linear.

rel.coral.loss.brm<-brm(rel.decline.prior~Disturb.year+(1|Disturbance)+(1|Reef_name), #tried poly(Disturb.year,degree=3) - no support
                        data=manual,
                        family=Gamma(link='log'),
                        chains=3,iter=2000,warmup=500,thin=2,
                        save_pars=save_pars(all = TRUE))

save(rel.coral.loss.brm,file='rel.coral.loss.brm.RData')

load(file='rel.coral.loss.brm.RData')

summary(rel.coral.loss.brm)

###extracting for linear
mcmc=rel.coral.loss.brm %>% as.matrix()
rel.coral.loss.brm.slope<-median_hdci(exp(mcmc[,2]))

### calculating 'slope' (percentage change +/- CIs)
(1-rel.coral.loss.brm.slope[1])*100
(1-rel.coral.loss.brm.slope[2])*100
(1-rel.coral.loss.brm.slope[3])*100
### exceedence probabilty
sum(mcmc[,2]<0)/length(mcmc[,2])


###extract trend for plot

predgrid<-with(manual,seq(min(Disturb.year),max(Disturb.year),length=1000))

rel.coral.loss.fitted.trend<-emmeans(rel.coral.loss.brm,~Disturb.year,
                                     at=list(Disturb.year=predgrid),type='response') %>%  as.data.frame()
rel.coral.loss.thru.time.plot<-rel.coral.loss.thru.time.plot+
  geom_ribbon(data=rel.coral.loss.fitted.trend,aes(y=response,x=Disturb.year,ymin=lower.HPD,ymax=upper.HPD),alpha=0.2)+
  geom_line(data=rel.coral.loss.fitted.trend,aes(y=response,x=Disturb.year),colour='blue')+
  ggtitle('b)')
# theme(plot.title = element_text(vjust = - 10))

rel.coral.loss.thru.time.plot




##########################################################################################
#### disturbance number and interval per decade ------------------------------


manual$Dist.Decade<-factor(manual$Dist.Decade,levels=c('86 to 90','91 to 00','01 to 10','11 to 20'))
# 
############################################################################################################
#### number of reefs with disturbances

reefs.with.disturbance.dec<-#length(unique(manual$Reef_name))
  manual %>% group_by(Dist.Decade) %>% dplyr::summarise(length(unique(Reef_name)))

colnames(reefs.with.disturbance.dec)[2]<-c("Reefs.with.disturbance")  

reefs.with.disturbance.region<-#length(unique(manual$Reef_name))
  manual %>% group_by(Region,Dist.Decade) %>% dplyr::summarise(length(unique(Reef_name)))

colnames(reefs.with.disturbance.dec)[2]<-c("Reefs.with.disturbance")
colnames(reefs.with.disturbance.region)[3]<-c("Reefs.with.disturbance")

reefs.with.disturbance.region<-reefs.with.disturbance.region[-5,]


###number of reefs surveyed
manta<-read.csv(file='manta tow by reef 2021.csv',strip.white=T)
head(manta)

manta<-manta %>% filter(REPORT_YEAR>1985) %>% 
  filter(SECTOR!='TS') %>% 
  mutate(cREPORT_YEAR=factor(REPORT_YEAR)) %>% 
  mutate(SECTOR=factor(SECTOR)) %>% 
  mutate(Dist.Decade=factor(Dist.Decade,levels=c("86 to 90","91 to 00","01 to 10","11 to 20"))) %>% 
  mutate(Region=factor(Region))


reefs.surveyed.dec<-manta %>% group_by(Dist.Decade) %>% dplyr::summarise(length(unique(REEF_NAME))) %>%
  as.data.frame() #%>% rename(reefs.surveyed=length(unique(REEF_NAME)))



reefs.with.disturbance.dec<-left_join(reefs.with.disturbance.dec,reefs.surveyed.dec) %>% as.data.frame() 
colnames(reefs.with.disturbance.dec)[3]<-c("Reefs.surveyed")

reefs.with.disturbance.dec<-reefs.with.disturbance.dec %>% 
  mutate(Dist.Decade=factor(Dist.Decade,levels=c("86 to 90","91 to 00","01 to 10","11 to 20"))) %>% 
  mutate(perc.reefs.disturb=Reefs.with.disturbance/Reefs.surveyed*100)

perc.reefs.disturb<-ggplot(filter(reefs.with.disturbance.dec,Dist.Decade!='NA'),aes(x=Dist.Decade,y=perc.reefs.disturb))+
  geom_bar(label='n',stat='identity',colour='black')+
  scale_y_continuous("Percent survey reefs  with disturbances")+
  scale_x_discrete("Decade of survey")+
  ggtitle('c)')+
  theme_classic()+
  theme(legend.position = c(0.275,0.95),
        legend.text = element_text(size=7.5),
        legend.key.size =unit(0.5,'cm'))

perc.reefs.disturb

library(grid)
library(gridExtra)



################################################################################
#### bleaching interval using COE and LTMP

interval<-read.csv(file='220509 bleaching intervals.csv',strip.white = T)
head(interval)

interval$Dist_decade<-factor(interval$Dist_decade,levels=c('91 to 00','01 to 10','11 to 20'))


### model


bleach.interval.brm<-brm(Bleach_int~Dist_decade+(Dist_decade|Reef),
                         #family=poisson(link='log'),
                         family=negbinomial(link = "log", link_shape = "log"),
                         data=interval,
                         chains=3,iter=2000,warmup=500,thin=2)



save(bleach.interval.brm,file='bleach.interval.brm.RData')
load(file='bleach.interval.brm.RData')

summary(bleach.interval.brm)

bleach.interval.est<-emmeans(bleach.interval.brm,~Dist_decade,type='response') %>% as.data.frame() %>% 
  rename(Decade=Dist_decade)



## 80s data is missing so add dummy variables

Decade<-"85 to 90"
prob<-0
lower.HPD<-0
upper.HPD<-0

eighties<-cbind(Decade,prob,lower.HPD,upper.HPD) %>% as.data.frame() %>% 
  mutate(Decade=factor(Decade),
         prob=as.numeric(prob),
         lower.HPD=as.numeric(lower.HPD),
         upper.HPD=as.numeric(upper.HPD))

bleach.interval.est<-rbind(eighties,bleach.interval.est)
str(bleach.interval.est)

## rename "95 to 00" as "91 to 00"

bleach.interval.est<-bleach.interval.est %>% 
  mutate(Decade=factor(Decade,levels=c("85 to 90","95 to 00","01 to 10","11 to 20"),
                       labels=c("85 to 90","91 to 00","01 to 10","11 to 20")))

### exceedence probabilty
mcmc=bleach.interval.brm %>% as.matrix()

sum(mcmc[,2]<0)/length(mcmc[,2])


modelled.bleach.interval.plot<-ggplot(bleach.interval.est,aes(x=Decade,y=prob))+
  geom_bar(stat='identity',fill='firebrick3')+
  geom_errorbar(stat = "identity",aes(ymin=lower.HPD,ymax=upper.HPD),width=0)+
  ggtitle('c)')+
  scale_y_continuous("Years between disturbance")+
  scale_x_discrete("Decade of first bleaching",limits=c("85 to 90","91 to 00","01 to 10","11 to 20"),
                   labels=c("85 to 90*","91 to 00**","01 to 10","11 to 20"))+
  #facet_wrap(~Region)+
  theme_classic()
modelled.bleach.interval.plot  

###############################################################
#### cyclones - 2022 update

load(file='cyclones.interval.RData')
head(cyclones.interval)

cyclones.interval.brm1<-brm(Interval~0 + Decade+(1|REEF_NAME),
                            #family=poisson(link='log'),
                            family=negbinomial(link = "log", link_shape = "log"),
                            data=cyclones.interval,# %>% filter(Decade != 2030),
                            chains=3,cores = 3,iter=2000,warmup=500,thin=2, seed = 1,
                            backend = 'cmdstanr')   
summary(cyclones.interval.brm1)

save(cyclones.interval.brm1,file='cyclones.interval.brm1.RData')
load(file='cyclones.interval.brm1.RData')


cyclones.summary<-emmeans(cyclones.interval.brm1,~Decade,type='response') %>% as.data.frame() 

cyclones.summary<-cyclones.summary %>% mutate(Decade = factor(Decade, levels = sort(as.numeric(as.character(unique(Decade)))),
                                                              labels=c("85 to 90","91 to 00","01 to 10","11 to 20"))) %>% 
  arrange(Decade) 


modelled.cyclone.interval.plot<-ggplot(cyclones.summary,aes(x=Decade,y=prob))+
  geom_bar(stat='identity',fill='blue3')+
  geom_errorbar(stat = "identity",aes(ymin=lower.HPD,ymax=upper.HPD),width=0)+
  ggtitle('d)')+
  scale_y_continuous("Years between disturbance",breaks=c(0,2,4,6,8,10))+
  scale_x_discrete("Decade of first cyclone",breaks=c("85 to 90","91 to 00","01 to 10","11 to 20"),
                   labels=c("85 to 90","91 to 00","01 to 10","11 to 20"))+
  #facet_wrap(~Region)+
  theme_classic()
modelled.cyclone.interval.plot  

str(cyclones)

###############################################################
####  cots  
##########################
##   time between outbreaks


cots.interval.outbreak<-read.csv(file='cots.interval.outbreak.csv',strip.white = T)

cots.interval.outbreak<-cots.interval.outbreak %>% 
  mutate(REEF_NAME=factor(REEF_NAME),
         Zone=factor(Zone),
         Decade=factor(Decade))


cots.interval.out.brm<-brm(Interval~0 + Decade+(Decade|REEF_NAME),
                           #family=poisson(link='log'),
                           family=negbinomial(link = "log", link_shape = "log"),
                           data=cots.interval.outbreak,
                           chains=3,iter=2000,warmup=500,thin=2, seed = 1)  

save(cots.interval.out.brm,file='cots.interval.out.brm.RData')
load(file='cots.interval.out.brm.RData')

summary(cots.interval.out.brm)

ave.cots.interval<-emmeans(cots.interval.out.brm, ~Decade, type='response') %>% as.data.frame() %>% 
  mutate(Decade=factor(Decade,levels=c(1980,1990,2000,2010)))


str(ave.cots.interval)

cots2010<-c(2010,0,0,0)
ave.cots.interval<-rbind(ave.cots.interval,cots2010) %>% 
  mutate(Decade=factor(Decade,levels=c(1980,1990,2000,2010),labels=c("85 to 90","91 to 00","01 to 10","11 to 20")))


modelled.cots.interval.plot<-ggplot(ave.cots.interval,aes(x=Decade,y=prob))+
  geom_bar(stat='identity',fill='seagreen')+
  geom_errorbar(stat = "identity",aes(ymin=lower.HPD,ymax=upper.HPD),width=0)+
  ggtitle('d)')+
  scale_y_continuous("Disturbance interval (years)",breaks=c(0,2,4,6,8,10))+
  scale_x_discrete("Decade of first cots outbreak")+
  #scale_x_discrete("Decade of first cots outbreak",labels=c("85 to 90","91 to 00","01 to 10","11 to 20"))+
  #facet_wrap(~Region)+
  theme_classic()
modelled.cots.interval.plot 


###############################################################################################
####  combined plot

###############################
bleach.interval.est<-bleach.interval.est %>% mutate(Disturbance=c('Bleaching','Bleaching','Bleaching','Bleaching'))
ave.cots.interval<-ave.cots.interval %>% mutate(Disturbance=c('COTS','COTS','COTS','COTS')) %>% mutate(Decade=factor(Decade))
# levels(ave.cots.interval$Decade)<-c('85 to 90','91 to 00','01 to 10','11 to 20')
cyclones.summary<-cyclones.summary %>% mutate(Disturbance=c('Cyclone','Cyclone','Cyclone','Cyclone'))

dist_int<-rbind(bleach.interval.est,ave.cots.interval,cyclones.summary) %>%  as.data.frame() 


modelled.combined.interval.plot<-ggplot(dist_int,aes(x=Decade,y=prob,fill=Disturbance))+
  geom_bar(stat='identity',width=0.9,position=position_dodge(width=0.97))+
  geom_errorbar(stat = "identity",aes(ymin=lower.HPD,ymax=upper.HPD),width=0,position=position_dodge(width=0.97))+
  ggtitle('d)')+
  scale_y_continuous("Disturbance interval (years)")+#,breaks=c(0,2,4,6,8,10))+
  scale_x_discrete("Decade of first disturbance")+#,labels=c("85 to 90","91 to 00","01 to 10","11 to 20"))+
  scale_fill_manual(values=c('firebrick3','seagreen','blue3'))+
  theme_classic()+
  theme(legend.position=c(0.8,0.9),
        legend.text = element_text(size=7),
        legend.title = element_text(size=8),
        legend.key.width=unit(0.25, "cm"),
        legend.key.height=unit(0.25, "cm"))
modelled.combined.interval.plot 

####################################################################################
#### combining severity and intervals

#### multi panel figure

png(file='221024 have disturbances become more severe and frequent.png',height=8,width=6,units='in',res=72)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))

pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
print(plot3, newpage = FALSE)
popViewport(1)

pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 1))
print(rel.coral.loss.thru.time.plot, newpage = FALSE)
popViewport(1)

pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2))
print(perc.reefs.disturb, newpage = FALSE)
popViewport(1)


pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 2))
print(modelled.combined.interval.plot, newpage = FALSE)
popViewport(1)

dev.off()


######################################################################################################
# ## Figure 5 ----------------------------------------

manual<-manual %>% mutate(abs.ann.prior.decline=abs.decline.prior/Disturb.time,
                          abs.ann.bench.decline=abs.decline.bench/Disturb.time) %>% 
  mutate(#years.dist.to.recovery=recovery.year-Disturb.year,
    recovered.prior=ifelse(rel.perc.prior.recovery>99.99999,'yes','no'),
    recovered.bench=ifelse(rel.perc.benchmark.recovery>99.99999,'yes','no'))

head(manual)

prior.total<-manual %>% dplyr::count(recovered.prior)
benchmark.total<-manual %>% dplyr::count(recovered.bench)


prop.recovered.bench<-manual %>% dplyr::filter(recovered.bench!='NA') %>% group_by(Region,recovery.year) %>% dplyr::count(recovered.bench) %>% 
  mutate(recovered.bench=factor(recovered.bench)) %>% 
  mutate(data4plot=ifelse(recovered.bench=='no',n*-1,n),
         Region=factor(Region,levels=c('Northern','Central','Southern'))) %>% ungroup()

prop.recovered.prior<-manual %>% dplyr::filter(recovered.prior!='NA') %>% 
  group_by(Region,recovery.year) %>% 
  dplyr::count(recovered.prior) %>% 
  mutate(recovered.prior=factor(recovered.prior)) %>% 
  mutate(data4plot=ifelse(recovered.prior=='no',n*-1,n)) %>% ungroup()


prop.all.bench<-manual %>% dplyr::filter(recovered.bench!='NA') %>%  dplyr::count(recovered.bench)

prop.all.prior<-manual %>% dplyr::filter(recovered.bench!='NA') %>%  dplyr::count(recovered.prior)

Years<-expand.grid(recovery.year=seq(from=1985,to=2022,by=1),
                   Region=levels(prop.recovered.bench$Region)) %>%
  mutate(recovery.year=factor(recovery.year))
Years %>% pull(recovery.year) %>% levels

###########################################
# plots ------------------------------------------


########################################################################
# histogram of recovery ---------------------------------------------------

bench.recovery.hist<-ggplot(manual,aes(x=rel.perc.benchmark.recovery))+
  geom_histogram(fill='grey',colour='black',binwidth=10,boundary=10)+
  geom_vline(xintercept=100,linetype='dashed',colour='red')+
  geom_vline(xintercept=50,linetype='dashed',colour='blue')+
  ggtitle('d)')+
  scale_x_continuous("Relative percent recovery")+
  scale_y_continuous("Count")+
  facet_wrap(~Region)+#,scales = 'free_x')+
  theme_classic()+
  theme(axis.text = element_text(size=8))

bench.recovery.hist


prior.recovery.hist<-ggplot(manual,aes(x=rel.perc.prior.recovery))+
  geom_histogram(fill='grey',colour='black',binwidth=10,boundary=10)+
  geom_vline(xintercept=100,linetype='dashed',colour='red')+
  geom_vline(xintercept=50,linetype='dashed',colour='blue')+
  ggtitle('c)')+
  scale_x_continuous("Relative percent recovery")+
  scale_y_continuous("Count")+
  facet_wrap(~Region)+#,scales = 'free_x')+
  theme_classic()+
  theme(axis.text = element_text(size=8))

prior.recovery.hist


############alter proportion plots to fit 2x2 grid

prop.bench.plot<-ggplot(prop.recovered.bench %>% filter(!is.na(Region)),# %>% filter(!is.na(Disturbance)),
                        aes(x=recovery.year,y=data4plot,fill=recovered.bench))+
  geom_bar(stat='identity',colour='black',show.legend = F)+
  scale_x_continuous('Year of recovery',breaks=c(1990,2000,2010,2020))+
  scale_y_continuous('Number of recoveries',breaks=c(-20,-15,-10,-5,0,5,10,15,20),labels=c(20,15,10,5,0,5,10,15,20),limits=c(-18,17))+
  scale_fill_manual("",values=c('grey','white'),breaks=c('no','yes'),labels=c('Partial recovery','Full recovery'))+
  geom_hline(yintercept=0)+
  ggtitle('b) Recovery to historical benchmark')+
  facet_grid(~Region)+
  #facet_grid(Disturbance~Region,scales = 'free_y')+
  theme_classic()+
  theme(legend.position=c(0.9,0.8),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        axis.text = element_text(size=8)) 

prop.bench.plot

prop.prior.plot<-ggplot(prop.recovered.prior %>% filter(!is.na(Region)),aes(x=recovery.year,y=data4plot,fill=recovered.prior))+
  geom_bar(stat='identity',colour='black')+
  scale_x_continuous('Year of recovery',breaks=c(1990,2000,2010,2020))+
  scale_y_continuous('Number of recoveries',breaks=c(-20,-15,-10,-5,0,5,10,15,20),labels=c(20,15,10,5,0,5,10,15,20),limits=c(-17,17))+
  scale_fill_manual("",values=c('grey','white'),breaks=c('no','yes'),labels=c('Partial recovery','Full recovery'))+
  geom_hline(yintercept=0)+
  ggtitle('a) Recovery to prior benchmark')+
  #facet_grid(~Region)+
  facet_grid(~Region)+
  theme_classic()+
  #guides(shape = guide_legend(override.aes = list(size = 0.1)))+
  theme(legend.position=c(0.175,0.9),
        legend.title = element_text(size=7),
        legend.text = element_text(size=8),
        legend.key.width=unit(0.25, "cm"),
        legend.key.height=unit(0.25, "cm"),
        axis.text = element_text(size=8)) 

prop.prior.plot


pdf(file='221021 proportion of reefs recovering plus histograms.pdf',height=8,width=10)

grid.newpage()
pushViewport(viewport(layout = grid.layout(4, 2)))

pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1:2))
print(prop.prior.plot, newpage = FALSE)
popViewport(1)

pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 3:4))
print(prior.recovery.hist, newpage = FALSE)
popViewport(1)

pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 1:2))
print(prop.bench.plot, newpage = FALSE)
popViewport(1)

pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 3:4))
print(bench.recovery.hist, newpage = FALSE)
popViewport(1)


dev.off()



############################################################################################################
#### Figure 6 --------------------------------------------------------


library(mgcv)
library(mgcv.helper)


##################################################################
##historic benchmark

#model1<-gam(rel.perc.benchmark.recovery~s(log(Disturb.cover))+
model1<-gam(rel.perc.benchmark.recovery~s(log(Disturb.cover))+
              s(Reef_name,bs="re")+
              s(bench.cover)+
              s(years.benchmark.recovery),
              #s(recovery.year),
            data=manual %>% filter(!is.na(rel.perc.benchmark.recovery)) %>% 
              mutate(Reef_name=factor(Reef_name),
                     recovery.year=as.numeric(recovery.year)),REML=T,family='scat')


plot(model1,pages=1,scale=0,seWithMean=T,shift=coef(model1)[1])
gam.check(model1)
mgcv::concurvity(model1)
summary(model1)

######## pretty plots

#dist.cover
gam.pred.distcover.bench<-emmeans(model1, specs = "log(Disturb.cover)",
                                  type='response'
                                  ,rg.limit=13000,
                                  at=list(`log(Disturb.cover)`=seq(min(log(manual$Disturb.cover)),
                                                            max(log(manual$Disturb.cover)),
                                                            length=100))) %>% 
  as.data.frame() %>%
  mutate(Disturb.cover = exp(`log(Disturb.cover)`))
  
## plot  
  
dist.cover.bench.gamplot<-ggplot(gam.pred.distcover.bench,aes(x=Disturb.cover,y=emmean))+
  geom_line()+
  geom_line(aes(y=lower.CL),linetype='dashed')+
  geom_line(aes(y=upper.CL),linetype='dashed')+
  scale_x_log10('Post-disturbance cover (%)',labels=function(x)x*100)+
  scale_y_continuous('')+
  ggtitle('a)')+
  theme_classic()

dist.cover.bench.gamplot

#bench.cover

gam.pred.benchcover.bench<-emmeans(model1,~bench.cover,type='response',rg.limit=13000,at=list(bench.cover=seq(min(manual$bench.cover),
                                                                                               max(manual$bench.cover),
                                                                                               length=100))) %>% 
  as.data.frame() 

bench.cover.bench.gamplot<-ggplot(gam.pred.benchcover.bench,aes(x=bench.cover,y=emmean))+
  geom_line()+
  geom_line(aes(y=lower.CL),linetype='dashed')+
  geom_line(aes(y=upper.CL),linetype='dashed')+
  scale_x_continuous('Historic benchmark cover (%)',labels=function(x)x*100)+
  scale_y_continuous('Recovery to historical benchmark (%)')+
  ggtitle('c)')+
  theme_classic()

bench.cover.bench.gamplot

newdata = manual %>% filter(!is.na(rel.perc.benchmark.recovery)) %>% 
  mutate(years.benchmark.recovery=as.numeric(years.benchmark.recovery)) %>% 
  tidyr::expand(years.benchmark.recovery)

gam.pred.recoveryyear.bench<-emmeans(model1,~years.benchmark.recovery,type='response',at=newdata) %>% 
  as.data.frame() 

recovery.year.bench.gamplot<-ggplot(gam.pred.recoveryyear.bench,aes(x=years.benchmark.recovery,y=emmean))+
  geom_line()+
  geom_line(aes(y=lower.CL),linetype='dashed')+
  geom_line(aes(y=upper.CL),linetype='dashed')+
  scale_x_continuous('Years for recovery')+#,labels=function(x)x*1)+
  scale_y_continuous('')+
  ggtitle('e)')+
  theme_classic()

recovery.year.bench.gamplot



##################################################################
##prior benchmark


model2<-gam(rel.perc.prior.recovery~s(log(Disturb.cover))+
              s(Reef_name,bs="re")+
              s(prior.cover)+
              s(years.prior.recovery),
            data=manual %>% filter(!is.na(rel.perc.prior.recovery)) %>% 
              mutate(Reef_name=factor(Reef_name),
                     recovery.year=as.numeric(recovery.year)),REML=T,family='scat')

summary(model2)
plot(model2,pages=1,scale=0,seWithMean=T,shift=coef(model2)[1])
gam.check(model2)
mgcv::concurvity(model2)


gam.pred.distcover.prior<-emmeans(model2, specs = "log(Disturb.cover)",
                                  type='response'
                                  ,rg.limit=13000,
                                  at=list(`log(Disturb.cover)`=seq(min(log(manual$Disturb.cover)),
                                                                   max(log(manual$Disturb.cover)),
                                                                   length=100))) %>% 
  as.data.frame() %>%
  mutate(Disturb.cover = exp(`log(Disturb.cover)`))



dist.cover.prior.gamplot<-ggplot(gam.pred.distcover.prior,aes(x=Disturb.cover,y=emmean))+
  geom_line()+
  geom_line(aes(y=lower.CL),linetype='dashed')+
  geom_line(aes(y=upper.CL),linetype='dashed')+
  scale_x_log10('Post-disturbance cover (%)',labels=function(x)x*100)+
  scale_y_continuous('')+
  ggtitle('b)')+
  theme_classic()

dist.cover.prior.gamplot



gam.pred.benchcover.prior<-emmeans(model2,~prior.cover,type='response',rg.limit=13000,at=list(prior.cover=seq(min(manual$prior.cover),
                                                                                               max(manual$prior.cover),
                                                                                               length=100))) %>% 
  as.data.frame() 

prior.cover.bench.gamplot<-ggplot(gam.pred.benchcover.prior,aes(x=prior.cover,y=emmean))+
  geom_line()+
  geom_line(aes(y=lower.CL),linetype='dashed')+
  geom_line(aes(y=upper.CL),linetype='dashed')+
  scale_x_continuous('Prior benchmark cover (%)',labels=function(x)x*100)+
  scale_y_continuous('Recovery to prior benchmark (%)')+
  ggtitle('d)')+
  theme_classic()

prior.cover.bench.gamplot

newdata1 = manual %>% filter(!is.na(rel.perc.prior.recovery)) %>% 
  mutate(years.prior.recovery=as.numeric(years.prior.recovery)) %>% 
  tidyr::expand(years.prior.recovery)

gam.pred.recoveryyear.prior<-emmeans(model2,~years.prior.recovery,type='response',at=newdata1) %>% 
  as.data.frame() 

recovery.year.prior.gamplot<-ggplot(gam.pred.recoveryyear.prior,aes(x=years.prior.recovery,y=emmean))+
  geom_line()+
  geom_line(aes(y=lower.CL),linetype='dashed')+
  geom_line(aes(y=upper.CL),linetype='dashed')+
  scale_x_continuous('Years for recovery')+#,labels=function(x)x*100)+
  scale_y_continuous('')+
  ggtitle('f)')+
  theme_classic()

recovery.year.prior.gamplot


pdf(file='230421 gam_plots.pdf',height=8,width=6)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 2)))

pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
print(dist.cover.bench.gamplot, newpage = FALSE)
popViewport(1)

pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2))
print(bench.cover.bench.gamplot, newpage = FALSE)
popViewport(1)

pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 3))
print(recovery.year.bench.gamplot, newpage = FALSE)
popViewport(1)

pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 1))
print(dist.cover.prior.gamplot, newpage = FALSE)
popViewport(1)

pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 2))
print(prior.cover.bench.gamplot, newpage = FALSE)
popViewport(1)

pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 3))
print(recovery.year.prior.gamplot, newpage = FALSE)
popViewport(1)

dev.off()



