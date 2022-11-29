source('CoralTrends_functions.R') 
CoralTrends_checkPackages()

INCLUDE_GBR = FALSE

## ---- load data

load(file='../data/processed/manta.sum.RData')

load(file='../data/modelled/mod.gbr_glmmTMB.beta.disp.RData')
load(file='../data/modelled/manta.tow.gbr.RData')
dat.gbr <- dat.gbr_glmmTMB.beta.disp %>% 
    mutate(Location = 'Great Barrier Reef',
           N = manta.tow.gbr %>% pull(REEF_NAME) %>% unique %>% length) %>%
    dplyr::rename(
               mean = response,
               lower = lower.CL,
               upper = upper.CL)

load(file='../data/modelled/mod.northern_glmmTMB.beta.disp.RData')
load(file='../data/modelled/manta.tow.northern.RData')
dat.northern <- dat.northern_glmmTMB.beta.disp %>% 
    mutate(Location = 'Northern GBR',
           N = manta.tow.northern %>% pull(REEF_NAME) %>% unique %>% length) %>%
    dplyr::rename(
               mean = response,
               lower = lower.CL,
               upper = upper.CL)

load(file='../data/modelled/mod.central_glmmTMB.beta.disp.RData')
load(file='../data/modelled/manta.tow.central.RData')
dat.central <- dat.central_glmmTMB.beta.disp %>% 
    mutate(Location = 'Central GBR',
           N = manta.tow.central %>% pull(REEF_NAME) %>% unique %>% length) %>%
    dplyr::rename(
               mean = response,
               lower = lower.CL,
               upper = upper.CL)
load(file='../data/modelled/mod.southern_glmmTMB.beta.disp.RData')
load(file='../data/modelled/manta.tow.southern.RData')
dat.southern <- dat.southern_glmmTMB.beta.disp %>% 
    mutate(Location = 'Southern GBR',
           N = manta.tow.southern %>% pull(REEF_NAME) %>% unique %>% length) %>%
    dplyr::rename(
               mean = response,
               lower = lower.CL,
               upper = upper.CL)

load(file='../data/spatial/spatial_3Zone.RData')

load(file='../data/modelled/cots.sum.all_3Zone.RData')
load(file='../data/modelled/bleaching.sum.all_3Zone.RData')
load(file='../data/modelled/cyclones.sum.all_3Zone.RData')

hues <- RColorBrewer::brewer.pal(4, "Blues")
## ----end

## ---- number of reefs per selected years
number_of_reefs <- manta.sum %>% mutate(Bins = case_when(
                         REPORT_YEAR <= 1990 ~ 'G1',
                         REPORT_YEAR <=2000 ~ 'G2',
                         REPORT_YEAR <=2010 ~ 'G3',
                         REPORT_YEAR >=2011 ~ 'G4')) %>%
    dplyr::select(Bins, REEF_NAME) %>%
    distinct() %>%
    group_by(Bins) %>%
    count()
## ----end

## ---- Generate the banner
a=oz:::ozRegion(sections=c(3,11:13))
a=oz:::ozRegion()
cc=rbind(xy2df(a$lines[[3]]),
    xy2df(a$lines[[13]]),
    xy2df(a$lines[[12]])[nrow(xy2df(a$lines[[12]])):1,],
    xy2df(a$lines[[11]]))
aa.ps<-SpatialPolygons(list(Polygons(list(Polygon(cc)),ID="QLD")))

gt = ggplot(fortify(aa.ps), aes(y=lat, x=long, group=group)) +
    geom_blank(aes(x=190,y=-20))+coord_map() +
    #geom_blank(aes(x=220,y=-20))+coord_map() +
    #geom_blank(aes(x=270,y=-20))+coord_map() +
    geom_polygon(data=fortify(spatial_3Zone), aes(y=lat, x=long),fill=NA,color='black', size=0.2) +
    coord_equal() +
    theme_classic() +theme(panel.background=element_rect(fill=NA),
                       axis.text.y=element_blank(),
                       axis.text.x=element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank(),
                       axis.ticks=element_blank(),
                       axis.line=element_blank(),
                       plot.background=element_blank(),
                       panel.spacing=unit(0,'pt'),
                       plot.margin=unit(c(0,0,0,0),'pt'))  
gt1=gt+geom_polygon(data=fortify(spatial_3Zone) , aes(y=lat, x=long),fill=hues[4],color=NA)+
geom_polygon(data=fortify(spatial_3Zone), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
    geom_polygon(fill='white', color=hues[4])
gt0=gt1

## ----end

## We have decided to bring the GBR wide figure out on its own
## ---- Figures for Lyndon
newdata=dat.gbr %>% mutate(Location=factor(Location,labels=c('Great Barrier\n\nReef')))
write.csv(newdata %>% arrange(Year) %>% dplyr::select(Year,mean,lower,upper), file='../data/data4Lyndon.csv', quote=FALSE, row.names=FALSE)
max_year = max(as.numeric(as.character(newdata$Year)))
nd = newdata %>%
    group_by(Location) %>% summarize(Year=mean(range(as.numeric(as.character(Year)))),N=paste0('(N=',unique(N),")"))
hues <- RColorBrewer::brewer.pal(4, "Blues")
g1<-ggplot(newdata, aes(y=mean*100, x=as.numeric(as.character(Year))))+
    geom_blank(aes(y=10,x=1995))+geom_blank(aes(y=35,x=1995))+
    annotate(geom='rect', xmin=2012, xmax=max_year, ymin=-Inf,ymax=Inf, fill='grey70',alpha=0.4) +
                                        #facet_wrap(~Location, nrow=1, scales='fixed')+
    facet_wrap(~Location, nrow=1, scales='fixed', labeller=labeller(Location=setNames(paste0("\n", levels(newdata$Location),"\n"), levels(newdata$Location))))+
    geom_blank()+
    #geom_ribbon(aes(ymin=lower*100, ymax=upper*100),fill=hues[2], alpha=0.75)+
    geom_pointrange(aes(ymin=lower*100, ymax=upper*100))+
    geom_line(aes(x=as.numeric(as.character(Year))), color='blue') +
    geom_text(data=nd, aes(y=Inf,x=Year, label=N), vjust=1.2) +
    
    #scale_y_continuous(expression(atop(Coral,cover~('%'))),expand=c(0,0),limits=c(0,35)) +
    scale_y_continuous(expression(Coral~cover~('%')),expand=c(0,0),limits=c(0,35)) +
    scale_x_continuous('',breaks=seq(1985,2020,by=5), limits=c(1985,2021))+
                theme_classic()+
                    theme(strip.background=element_rect(fill=hues[2], color='black', size=0.5),
                          panel.background=element_rect(color='black'),
                          axis.title.y=element_text(size=rel(1.5), margin=margin(r=1,unit='lines')),
                          axis.text.x=element_text(size=rel(1.0)), axis.title.x=element_blank(),
                          axis.text.y=element_text(size=rel(1.2)),
                          panel.grid.minor=element_line(size=0.1,color=NA),
                          panel.grid.major=element_line(size=0.1,color='gray70'),
                          panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
                          panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
                          strip.text=element_text(margin=margin(t=2, b=2,unit='lines'),size=20,lineheight=0.5, face='bold',hjust=0.9,vjust=-1),
                          plot.margin=unit(c(0,0,2,0),'pt'))

g1l <- ggplot(newdata, aes(y=mean*100, x=as.numeric(as.character(Year))))+
    geom_blank(aes(y=10,x=1995))+geom_blank(aes(y=35,x=1995))+
    facet_wrap(~Location, nrow=1, scales='fixed', labeller=labeller(Location=setNames(paste0("\n", levels(newdata$Location),"\n"), levels(newdata$Location))))+
    geom_blank()+
    #geom_ribbon(aes(ymin=lower*100, ymax=upper*100),fill=hues[2], alpha=0.75)+
    geom_pointrange(aes(ymin=lower*100, ymax=upper*100))+
    geom_line(aes(x=as.numeric(as.character(Year))), color='blue') +
    geom_text(data=nd, aes(y=Inf,x=Year, label=N), vjust=1.2) +
    geom_smooth(method='lm', se=FALSE) + 
    #scale_y_continuous(expression(atop(Coral,cover~('%'))),expand=c(0,0),limits=c(0,35)) +
    scale_y_continuous(expression(Coral~cover~('%')),expand=c(0,0),limits=c(0,35)) +
    scale_x_continuous('',breaks=seq(1985,2020,by=5), limits=c(1985,2021))+
                theme_classic()+
                    theme(strip.background=element_blank(),
                          panel.background=element_rect(color='black'),
                          axis.title.y=element_text(size=rel(1.5), margin=margin(r=1,unit='lines')),
                          axis.text.x=element_text(size=rel(1.0)), axis.title.x=element_blank(),
                          axis.text.y=element_text(size=rel(1.2)),
                          panel.grid.minor=element_line(size=0.1,color=NA),
                          panel.grid.major=element_line(size=0.1,color='gray70'),
                          panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
                          panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
                          strip.text=element_blank(),
                          plot.margin=unit(c(0,0,2,0),'pt'))
ggsave(filename='../output/figures/figure4Lyndon.pdf', g1l, width=5, height=3)

##Put on the banner
gT <- ggplot_gtable(ggplot_build(g1))
facets <- grep("strip-t-1-1", gT$layout$name)
gg.gbr <- with(gT$layout[facets,],
           gtable_add_grob(gT, ggplotGrob(gt1),t=t, l=5, b=b, r=5, name="pic_predator"))

grid.draw(gg.gbr)
save(gg.gbr, file='../data/spatial/gg.gbr_3Zone.RData')
ggsave(filename='../output/figures/GBRStan_3Zone.pdf', gg.gbr, width=4, height=4, units='in',dpi=300) 
ggsave(filename='../output/figures/GBRStan_3Zone.png', gg.gbr, width=4, height=4, units='in',dpi=300) 
png(filename='../output/figures/GBRStan_3Zone.png', width=4, height=4, units='in', res=300)
grid.draw(gg.gbr)
dev.off()

## ----end

## ---- GBR
newdata = rbind(dat.gbr,dat.northern,dat.central,dat.southern)
newdata$Location <- factor(newdata$Location, levels=unique(newdata$Location),
                           labels=c('Great Barrier\n\nReef','Northern\n', 'Central\n', 'Southern\n'))
newdata = newdata %>% filter(Location=='Great Barrier\n\nReef') %>% droplevels
nd = newdata %>%
    group_by(Location) %>% summarize(Year=mean(range(as.numeric(as.character(Year)))),N=paste0('(N=',unique(N),")"))
hues <- RColorBrewer::brewer.pal(4, "Blues")
g1<-ggplot(newdata, aes(y=mean*100, x=as.numeric(as.character(Year))))+
    geom_blank(aes(y=10,x=1995))+geom_blank(aes(y=35,x=1995))+
                                        #facet_wrap(~Location, nrow=1, scales='fixed')+
    facet_wrap(~Location, nrow=1, scales='fixed', labeller=labeller(Location=setNames(paste0("\n", levels(newdata$Location),"\n"), levels(newdata$Location))))+
    geom_blank()+
    #geom_ribbon(aes(ymin=lower*100, ymax=upper*100),fill=hues[2], alpha=0.75)+
    geom_pointrange(aes(ymin=lower*100, ymax=upper*100))+
    geom_line(aes(x=as.numeric(as.character(Year))), color='blue') +
    #geom_text(data=nd, aes(y=Inf,x=Year, label=N), vjust=1.2) +
    #scale_y_continuous(expression(atop(Coral,cover~('%'))),expand=c(0,0),limits=c(0,60)) +
    scale_y_continuous(expression(Coral~cover~('%')),expand=c(0,0),limits=c(0,60)) +
    scale_x_continuous('',breaks=seq(1985,2020,by=5), limits=c(1985,2021))+
    theme_classic()+
    theme(strip.background=element_rect(fill=hues[2], color='black', size=0.5),
          panel.background=element_rect(color='black'),
                                        #axis.title.y=element_text(size=rel(1.5), margin=margin(r=1,unit='lines')),
          axis.title.y=element_text(size=rel(1.5), margin=margin(r=0.5,unit='lines')),
          axis.text.x=element_text(size=rel(1.0)), axis.title.x=element_blank(),
          axis.text.y=element_text(size=rel(1.2)),
          panel.grid.minor=element_line(size=0.1,color=NA),
          panel.grid.major=element_line(size=0.1,color='gray70'),
          panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
          panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
                                        #strip.text=element_text(margin=margin(t=2, b=2,unit='lines'),size=20,lineheight=0.5, face='bold',hjust=0.9,vjust=-1),
          strip.text=element_text(margin=margin(t=0.5, b=0.5,unit='lines'),size=20,lineheight=0.5, face='bold',hjust=0.9,vjust=-1),
          plot.margin=unit(c(0,0,2,0),'pt'))

## Add banner map
gT <- ggplot_gtable(ggplot_build(g1))
facets <- grep("strip-t-1-1", gT$layout$name)
gg <- with(gT$layout[facets,],
           gtable_add_grob(gT, ggplotGrob(gt1),t=t, l=5, b=b, r=5, name="pic_predator"))
grid.draw(gg)
save(gg, file='../data/spatial/gg_3Zone_GBR.RData')
ggsave(file='../output/figures/3ZonesGBR.pdf', gg, width=5, height=5, units='in',dpi=300) 
ggsave(file='../output/figures/3ZonesGBR.png', gg, width=5, height=5, units='in',dpi=300)
png(file='../output/figures/3ZonesGBR.png', width=5, height=5, units='in',res=300)
grid.draw(gg)
dev.off()

## ----end

## ---- multipanel 
{
    ## ---- make banners
    gt1=gt+geom_polygon(data=fortify(spatial_3Zone) , aes(y=lat, x=long),fill=hues[4],color=NA)+
        geom_polygon(data=fortify(spatial_3Zone), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
        geom_polygon(fill='white', color=hues[4]) +
                                        #annotate(geom='text', x=Inf, y=Inf, label='Far\nNorthern', vjust=1,hjust=1)
        geom_blank(aes(x=250,y=-20))

    gt2=gt+geom_polygon(data=fortify(spatial_3Zone[1]), aes(y=lat, x=long),fill=hues[4],color=NA)+
        geom_polygon(data=fortify(spatial_3Zone[1]), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
        geom_polygon(fill='white', color=hues[4]) +
        geom_blank(aes(x=250,y=-20))

    gt3=gt+geom_polygon(data=fortify(spatial_3Zone[2]), aes(y=lat, x=long),fill=hues[4],color=NA)+
        geom_polygon(data=fortify(spatial_3Zone[2]), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
        geom_polygon(fill='white', color=hues[4]) +
        geom_blank(aes(x=250,y=-20))

    gt4=gt+geom_polygon(data=fortify(spatial_3Zone[3]), aes(y=lat, x=long),fill=hues[4],color=NA)+
        geom_polygon(data=fortify(spatial_3Zone[3]), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
        geom_polygon(fill='white', color=hues[4]) +
        geom_blank(aes(x=250,y=-20))


    if (!INCLUDE_GBR) gt1=gt2; gt2=gt3; gt3=gt4; #gt4=gt5;

    save(gt1, gt2, gt3, file='../data/spatial/gts.RData')
    ## ----end
    ## ---- data
    newdata = rbind(dat.gbr,dat.northern,dat.central,dat.southern)
    newdata$Location <- factor(newdata$Location, levels=unique(newdata$Location),
                               labels=c('Great Barrier\n\nReef','Northern\n', 'Central\n', 'Southern\n'))
    save(newdata, file='../data/modelled/cellmeans.RData')
    if (!INCLUDE_GBR) newdata = newdata %>% filter(Location!='Great Barrier\n\nReef') %>% droplevels
                                        #newdata = newdata %>% left_join(dat.all %>% select(Location, N) %>% distinct)
    nd = newdata %>%
        group_by(Location) %>% summarize(Year=mean(range(as.numeric(as.character(Year)))),N=paste0('(N=',unique(N),")"))
    hues <- RColorBrewer::brewer.pal(4, "Blues")
    ## ----end
    ## ---- base plot
    g1<-ggplot(newdata, aes(y=mean*100, x=as.numeric(as.character(Year))))+
        geom_blank(aes(y=10,x=1995))+geom_blank(aes(y=35,x=1995))+
                                        #facet_wrap(~Location, nrow=1, scales='fixed')+
        facet_wrap(~Location, nrow=1, scales='fixed', labeller=labeller(Location=setNames(paste0("\n", levels(newdata$Location),"\n"), levels(newdata$Location))))+
        geom_blank()+
                                        #geom_ribbon(aes(ymin=lower*100, ymax=upper*100),fill=hues[2], alpha=0.75)+
        geom_pointrange(aes(ymin=lower*100, ymax=upper*100))+
        geom_line(aes(x=as.numeric(as.character(Year))), color='blue') +
                                        #geom_text(data=nd, aes(y=Inf,x=Year, label=N), vjust=1.2) +
                                        #scale_y_continuous(expression(atop(Coral,cover~('%'))),expand=c(0,0),limits=c(0,60)) +
        scale_y_continuous(expression(Coral~cover~('%')),expand=c(0,0),limits=c(0,60)) +
        scale_x_continuous('',breaks=seq(1985,2020,by=5), limits=c(1985,2021))+
        theme_classic()+
        theme(strip.background=element_rect(fill=hues[2], color='black', size=0.5),
              panel.background=element_rect(color='black'),
                                        #axis.title.y=element_text(size=rel(1.5), margin=margin(r=1,unit='lines')),
              axis.title.y=element_text(size=rel(1.5), margin=margin(r=0.5,unit='lines')),
              axis.text.x=element_text(size=rel(1.0)), axis.title.x=element_blank(),
              axis.text.y=element_text(size=rel(1.2)),
              panel.grid.minor=element_line(size=0.1,color=NA),
              panel.grid.major=element_line(size=0.1,color='gray70'),
              panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
              panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
                                        #strip.text=element_text(margin=margin(t=2, b=2,unit='lines'),size=20,lineheight=0.5, face='bold',hjust=0.9,vjust=-1),
              strip.text=element_text(margin=margin(t=0.5, b=0.5,unit='lines'),size=20,lineheight=0.5, face='bold',hjust=0.9,vjust=-1),
              plot.margin=unit(c(0,0,2,0),'pt'))
    save(g1,file='../data/modelled/g1_trends.Rmd')
    ## ----end
    ## ---- final plot
    gT <- ggplot_gtable(ggplot_build(g1))
    facets <- grep("strip-t-1-1", gT$layout$name)
    gg <- with(gT$layout[facets,],
               gtable_add_grob(gT, ggplotGrob(gt1),t=t, l=5, b=b, r=5, name="pic_predator"))
    facets <- grep("strip-t-2-1", gT$layout$name)
    gg <- with(gg$layout[facets,],
               gtable_add_grob(gg, ggplotGrob(gt2),t=t, l=9, b=b, r=9, name="pic_predator"))
    facets <- grep("strip-t-3-1", gT$layout$name)
    gg <- with(gg$layout[facets,],
               gtable_add_grob(gg, ggplotGrob(gt3),t=t, l=13, b=b, r=13, name="pic_predator"))
    if (INCLUDE_GBR) {
        facets <- grep("strip-t-5-1", gT$layout$name)
        gg <- with(gg$layout[facets,],
                   gtable_add_grob(gg, ggplotGrob(gt5),t=t, l=17, b=b, r=20, name="pic_predator"))
    }

    grid.draw(gg)
    save(gg, file='../data/spatial/gg_3Zone.RData')
    ggsave(file='../output/figures/3Zones.pdf', gg, width=10, height=3, units='in',dpi=300) 
    ggsave(file='../output/figures/3Zones.png', gg, width=10, height=3, units='in',dpi=300) 
    png(file='../output/figures/3Zones.png', width=10, height=3, units='in',res=300) 
    grid.draw(gg)
    dev.off()


    ## ----end
    ## ---- COTS plot
    if (!INCLUDE_GBR) cots.sum.all = cots.sum.all %>% filter(Location!='Great Barrier Reef') %>% droplevels
    labs = levels(cots.sum.all$Location)
                                        #Version with NO
    gcots <- ggplot(cots.sum.all %>% filter(REPORT_YEAR >1985, COTScat %in% c('NO','AO','IO')), aes(y=COTS.p, x=REPORT_YEAR)) +
        geom_bar(stat='identity',position='stack',aes(alpha=COTScat), fill='blue',show.legend=FALSE) +
        facet_wrap(~Location, nrow=1, strip.position='bottom',scales='fixed', labeller=labeller(Location=setNames(paste0("", labs, "\n"), labs)))+
        scale_alpha_manual(breaks=c('NO','IO','AO'), values=c(0.2,0.4,1)) +
        scale_y_reverse(expression(atop(Reefs,impacted~('%'))),expand=c(0,0),lim=c(90,0)) +
        scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top") + #limits=c(1985,(finalYear+1)))+
        theme_classic() +
        coord_cartesian(xlim=c(1985, finalYear)) +
        theme(
                                        #strip.background=element_rect(fill='#00009930', color='black', size=0.5),
            panel.border=element_rect(fill=NA,color='black'),
            axis.title.y=element_text(size=rel(1.5),margin=margin(r=1,unit='lines')),
                                        #axis.text.x=element_text(size=rel(1.2)),
            axis.text.x=element_blank(), axis.title.x=element_blank(),
            axis.text.y=element_text(size=rel(1.2)),
            panel.grid.minor=element_line(size=0.1,color=NA),
            panel.grid.major=element_line(size=0.1,color='gray70'),
            panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
            panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
            plot.margin=unit(c(0,0,0,0),'pt'),
            strip.text.x=element_blank(), strip.background=element_blank())
                                        #Version without NO
    gcots = ggplot(cots.sum.all %>% filter(REPORT_YEAR >1985, COTScat %in% c('AO','IO')), aes(y=COTS.p, x=REPORT_YEAR)) +
        geom_bar(stat='identity',position='stack',aes(alpha=COTScat), fill='blue',show.legend=FALSE) +
                                        #geom_line(stat='identity',position='stack',aes(alpha=COTScat), fill='red',show.legend=FALSE) +
                                        #facet_grid(~Location) +
                                        #facet_wrap(~Location, nrow=1, scales='fixed', labeller=labeller(Location=setNames(paste0("", labs, "\n"), labs)))+
        facet_wrap(~Location, nrow=1, strip.position='bottom',scales='fixed', labeller=labeller(Location=setNames(paste0("", labs, "\n"), labs)))+
        scale_alpha_manual(breaks=c('IO','AO'), values=c(0.4,1)) +
        scale_y_reverse(expression(atop(Reefs,impacted~('%'))),expand=c(0,0),lim=c(90,0)) +
        scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top") + # limits=c(1985,(finalYear+1)))+
        theme_classic() +
        coord_cartesian(xlim=c(1985, finalYear)) +
        theme(
                                        #strip.background=element_rect(fill='#00009930', color='black', size=0.5),
            panel.border=element_rect(fill=NA,color='black'),
            axis.title.y=element_text(size=rel(1.5),margin=margin(r=1,unit='lines')),
                                        #axis.text.x=element_text(size=rel(1.2)),
            axis.text.x=element_blank(), axis.title.x=element_blank(),
            axis.text.y=element_text(size=rel(1.2)),
            panel.grid.minor=element_line(size=0.1,color=NA),
            panel.grid.major=element_line(size=0.1,color='gray70'),
            panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
            panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
            plot.margin=unit(c(0,0,0,0),'pt'),
            strip.text.x=element_blank(), strip.background=element_blank())
    gcots

    gcots <- ggplot_gtable(ggplot_build(gcots))
    gcots$widths = gg$widths

    grid.arrange(gg,gcots, nrow=2,heights=c(2,1),padding=unit(-1,'lines'))

    library(viridis)
    library(RColorBrewer)
    labs = levels(cots.sum.all$Location)
    cots.dat = cots.sum.all %>% filter(REPORT_YEAR >1985, COTScat %in% c('IO','AO'))
    gcots = ggplot(cots.dat, aes(y=COTS.p, x=REPORT_YEAR-0.3)) +
        geom_bar(stat='identity',position='stack',aes(fill=COTScat),width=0.1,show.legend=FALSE)+
        geom_bar(data=cots.dat %>% filter(COTScat %in% c('IO','AO')),stat='identity',position='stack',aes(fill=COTScat),width=0.4,show.legend=FALSE)+
        geom_bar(data=cots.dat %>% filter(COTScat %in% c('AO')),stat='identity',position='stack',aes(fill=COTScat),width=0.4,show.legend=FALSE)+
                                        #geom_col(position='stack',aes(fill=COTScat),width=0.1,show.legend=TRUE)+
                                        #geom_area(position='stack',aes(fill=COTScat),show.legend=TRUE)+
        scale_fill_manual('', breaks=c('IO','AO'), labels=c('IO','AO'), values=scales:::brewer_pal(palette='Greens')(3)[-1]) +
                                        #geom_bar(stat='identity',position='stack',aes(alpha=COTScat), fill=brewer.pal(3,'Dark2')[1],show.legend=FALSE,width=0.3) +
                                        #scale_alpha_manual(breaks=c('NO','IO','AO'), values=c(0.2,0.4,1)) +
        facet_wrap(~Location, nrow=1, strip.position='bottom',scales='fixed', labeller=labeller(Location=setNames(paste0("", labs, "\n"), labs)))+
        scale_y_reverse(expression(atop(Reefs,impacted~('%'))),expand=c(0,0),lim=c(100,0)) +
        scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top") + #, limits=c(1985,2021))+
        coord_cartesian(xlim=c(1985, finalYear)) +
        theme_classic() +
        theme(
                                        #strip.background=element_rect(fill='#00009930', color='black', size=0.5),
            panel.border=element_rect(fill=NA,color='black'),
            axis.title.y=element_text(size=rel(1.5),margin=margin(r=1,unit='lines')),
                                        #axis.text.x=element_text(size=rel(1.2)),
            axis.text.x=element_blank(), axis.title.x=element_blank(),
            axis.text.y=element_text(size=rel(1.2)),
            panel.grid.minor=element_line(size=0.1,color=NA),
            panel.grid.major=element_line(size=0.1,color='gray70'),
            panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
            panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
            plot.margin=unit(c(2,0,5,0),'pt'),
            strip.text.x=element_blank(), strip.background=element_blank())
    gcots

    ## ----end
    ## ---- Bleaching plot
    if (!INCLUDE_GBR) bleaching.sum.all = bleaching.sum.all %>% filter(Location!='Great Barrier Reef') %>% droplevels
    bleaching.dat <- bleaching.sum.all %>% filter(REPORT_YEAR >1985, REPORT_YEAR < (finalYear+1), BLEACHINGcat!='0')
    gbleaching <- ggplot(bleaching.dat, aes(y=BLEACHING.p*100, x=REPORT_YEAR+0.3)) +
        geom_bar(stat='identity',position='stack',aes(fill=BLEACHINGcat),width=0.4,show.legend=FALSE)+
        ## geom_bar(data=bleaching.dat %>% filter(BLEACHINGcat %in% c('2','3','4')),stat='identity',position='stack',aes(fill=BLEACHINGcat),width=0.4,show.legend=FALSE)+
        ## geom_bar(data=bleaching.dat %>% filter(BLEACHINGcat %in% c('3','4')),stat='identity',position='stack',aes(fill=BLEACHINGcat),width=0.4,show.legend=FALSE)+
        ## geom_bar(data=bleaching.dat %>% filter(BLEACHINGcat %in% c('4')),stat='identity',position='stack',aes(fill=BLEACHINGcat),width=0.4,show.legend=FALSE)+
        
                                        #geom_bar(stat='identity',position='stack',aes(fill=BLEACHINGcat),show.legend=FALSE,width=0.3)+
        scale_fill_manual('', breaks=c(1,2,3,4,5), labels=c(1,2,3,4,5), values=c(scales:::brewer_pal(palette='Reds')(6)[-1])) +
                                        #geom_bar(stat='identity',position='stack',aes(alpha=BLEACHINGcat), fill=brewer.pal(3,'Dark2')[2],show.legend=FALSE,width=0.3) +
        facet_wrap(~Location, nrow=1, strip.position='bottom',scales='fixed', labeller=labeller(Location=setNames(paste0("", labs, "\n"), labs)))+
                                        #scale_alpha_manual(breaks=c(0,1,2,3,4), values=c(0,0.25,0.5,0.75,1)) +
        scale_y_reverse(expression(atop(Reefs,impacted~('%'))),expand=c(0,0),lim=c(100,0)) +
        scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top") + #, limits=c(1985,2021))+
        coord_cartesian(xlim=c(1985, finalYear)) +
        theme_classic() +
        theme(
                                        #strip.background=element_rect(fill='#00009930', color='black', size=0.5),
            panel.border=element_rect(fill=NA,color='black'),
            axis.title.y=element_text(size=rel(1.5),margin=margin(r=1,unit='lines')),
                                        #axis.text.x=element_text(size=rel(1.2)),
            axis.text.x=element_blank(), axis.title.x=element_blank(),
            axis.text.y=element_text(size=rel(1.2)),
            panel.grid.minor=element_line(size=0.1,color=NA),
            panel.grid.major=element_line(size=0.1,color='gray70'),
            panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
            panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
            plot.margin=unit(c(2,0,5,0),'pt'),
            strip.text.x=element_blank(), strip.background=element_blank())
    gbleaching

    ## ----end
    ## ---- Cyclones plot
    if (!INCLUDE_GBR) cyclones.sum.all = cyclones.sum.all %>% filter(Location!='Great Barrier Reef') %>% droplevels
    cyclones.dat = cyclones.sum.all %>% filter(REPORT_YEAR >1985, CYCLONEcat!=0)
    gcyclones = ggplot(cyclones.dat, aes(y=CYCLONE.p, x=REPORT_YEAR+0)) +
        geom_bar(stat='identity',position='stack',aes(fill=CYCLONEcat),width=0.4,show.legend=FALSE)+
        geom_bar(data=cyclones.dat %>% filter(CYCLONEcat %in% c('2','3')),stat='identity',position='stack',aes(fill=CYCLONEcat),width=0.4,show.legend=FALSE)+
        geom_bar(data=cyclones.dat %>% filter(CYCLONEcat %in% c('3')),stat='identity',position='stack',aes(fill=CYCLONEcat),width=0.4,show.legend=FALSE)+

                                        #geom_bar(stat='identity',position='stack',aes(alpha=CYCLONEcat), fill=brewer.pal(3,'Dark2')[3],show.legend=FALSE,width=0.3) +
                                        #geom_bar(stat='identity',position='stack',aes(fill=CYCLONEcat),show.legend=FALSE,width=0.3)+
    scale_fill_manual('', breaks=c(1,2,3), labels=c(1,2,3), values=c(scales:::brewer_pal(palette='Blues')(4)[-1])) +
    facet_wrap(~Location, nrow=1, strip.position='bottom',scales='fixed', labeller=labeller(Location=setNames(paste0("", labs, "\n"), labs)))+
                                        #scale_alpha_manual(breaks=c(0,1,2,3,4,5), values=c(0,1/5,2/5,3/5,4/5,5/5)) +
    scale_y_reverse(expression(atop(Reefs,impacted~('%'))),expand=c(0,0),lim=c(100,0)) +
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top") + #, limits=c(1985,2021))+
    coord_cartesian(xlim=c(1985, finalYear)) +
    theme_classic() +
    theme(
                                        #strip.background=element_rect(fill='#00009930', color='black', size=0.5),
        panel.border=element_rect(fill=NA,color='black'),
        axis.title.y=element_text(size=rel(1.5),margin=margin(r=1,unit='lines')),
                                        #axis.text.x=element_text(size=rel(1.2)),
        axis.text.x=element_blank(), axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(2,0,5,0),'pt'),
        strip.text.x=element_blank(), strip.background=element_blank())
    gcyclones


    ## ----end
    ## ---- Put them all together
    gbleaching = gbleaching + theme(panel.background=element_blank())
    gcyclones = gcyclones + theme(panel.background=element_blank())
    gcots <- ggplot_gtable(ggplot_build(gcots))
    gbleaching <- ggplot_gtable(ggplot_build(gbleaching))
    gcyclones <- ggplot_gtable(ggplot_build(gcyclones))
    panels <- grepl("panel", gbleaching$layout$name)

    pp <- c(subset(gcots$layout, grepl("panel", gcots$layout$name), se = t:r))

                                        # Overlap panels for second plot on those of the first plot
    gT <- gtable_add_grob(gcots, gbleaching$grobs[grepl("panel", gcots$layout$name)], 
                          pp$t, pp$l, pp$b, pp$l, name='bleaching')
    gT <- gtable_add_grob(gT, gcyclones$grobs[grepl("panel", gcots$layout$name)], 
                          pp$t, pp$l, pp$b, pp$l, name='cyclones')
    gT$widths = gg$widths

    grid.arrange(gg,gT, nrow=2,heights=c(2,1),padding=unit(0,'lines'))
    save(gT, file='../data/spatial/gT_3Zone.RData')
    ggsave(file='../output/figures/3ZonesFigurePts3_new.pdf',
           grid.arrange(gg,gT, nrow=2,heights=c(2,1),padding=unit(0,'lines')),
           width=10, height=5, units='in',dpi=300) 
    ggsave(file='../output/figures/3ZonesFigurePts1_new.png',
           grid.arrange(gg,gT, nrow=2,heights=c(2,1),padding=unit(0,'lines')),
           width=15, height=5, units='in',dpi=300) 
    png(file='../output/figures/3ZonesFigurePts1_new.png', width=15, height=5, units='in',res=300)
    grid.arrange(gg,gT, nrow=2,heights=c(2,1),padding=unit(0,'lines'))
    dev.off()
    ## ----end
    ## ---- Disturbances only
    {
        ## ---- COTS
        library(viridis)
        library(RColorBrewer)
        if (!INCLUDE_GBR) cots.sum.all = cots.sum.all %>% filter(Location!='Great Barrier Reef') %>% droplevels
        labs = levels(cots.sum.all$Location)
        cots.dat = cots.sum.all %>% filter(REPORT_YEAR >1985, COTScat %in% c('IO','AO'))
        gcots = ggplot(cots.dat, aes(y=COTS.p, x=REPORT_YEAR-0.3)) +
            geom_bar(stat='identity',position='stack',aes(fill=COTScat),width=0.1,show.legend=FALSE)+
            geom_bar(data=cots.dat %>% filter(COTScat %in% c('IO','AO')),stat='identity',position='stack',aes(fill=COTScat),width=0.3,show.legend=TRUE)+
            geom_bar(data=cots.dat %>% filter(COTScat %in% c('AO')),stat='identity',position='stack',aes(fill=COTScat),width=0.3,show.legend=FALSE)+
                                        #geom_col(position='stack',aes(fill=COTScat),width=0.1,show.legend=TRUE)+
                                        #geom_area(position='stack',aes(fill=COTScat),show.legend=TRUE)+
            scale_fill_manual('COTS outbreak status', breaks=c('IO','AO'), labels=c('IO','AO'), values=scales:::brewer_pal(palette='Greens')(3)[-1]) +
                                        #geom_bar(stat='identity',position='stack',aes(alpha=COTScat), fill=brewer.pal(3,'Dark2')[1],show.legend=FALSE,width=0.3) +
                                        #scale_alpha_manual(breaks=c('NO','IO','AO'), values=c(0.2,0.4,1)) +
            facet_wrap(~Location, nrow=1, strip.position='top',scales='fixed', labeller=labeller(Location=setNames(paste0("", labs, "\n"), labs)))+
                                        #scale_y_reverse(expression(atop(Reefs,impacted~('%'))),expand=c(0,0),lim=c(100,0)) +
            scale_y_continuous(expression(Reefs~impacted~('%')),expand=c(0,0),lim=c(0,100)) +
            scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "bottom") +#, limits=c(1985,2021))+
            coord_cartesian(xlim=c(1985, finalYear)) +
            theme_classic() +
            theme(strip.background=element_rect(fill=hues[2], color='black', size=0.5),
                                        #strip.background=element_rect(fill='#00009930', color='black', size=0.5),
                  panel.border=element_rect(fill=NA,color='black'),
                  axis.title.y=element_text(size=rel(1.5),margin=margin(r=0.5,unit='lines')),
                                        #axis.text.x=element_text(size=rel(1.2)),
                  axis.text.x=element_text(size=rel(1.0)), axis.title.x=element_blank(),
                  axis.text.y=element_text(size=rel(1.2)),
                  panel.grid.minor=element_line(size=0.1,color=NA),
                  panel.grid.major=element_line(size=0.1,color='gray70'),
                  panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
                  panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
                  plot.margin=unit(c(2,5,5,0),'pt'),
                  legend.position = 'bottom',
                  strip.text=element_text(margin=margin(t=0.5, b=0.5,unit='lines'),size=20,lineheight=0.5, face='bold',hjust=0.9,vjust=-1))
                                        #strip.text.x=element_blank(), strip.background=element_blank())
        gcots
        gT <- ggplot_gtable(ggplot_build(gcots))

        
        facets <- grep("strip-t-1-1", gT$layout$name)
        gg <- with(gT$layout[facets,],
                   gtable_add_grob(gT, ggplotGrob(gt1),t=t, l=5, b=b, r=6, name="pic_predator"))
        facets <- grep("strip-t-2-1", gT$layout$name)
        gg <- with(gg$layout[facets,],
                   gtable_add_grob(gg, ggplotGrob(gt2),t=t, l=9, b=b, r=9, name="pic_predator"))
        facets <- grep("strip-t-3-1", gT$layout$name)
        gg <- with(gg$layout[facets,],
                   gtable_add_grob(gg, ggplotGrob(gt3),t=t, l=13, b=b, r=13, name="pic_predator"))
        grid.draw(gg)
        ## ggsave(file='output/figures/Disturbances_cots_no_label.png',grid.draw(gg),width=9, height=3, dpi=300)
        ## ggsave(file='output/figures/Disturbances_cots_no_label.pdf',grid.draw(gg),width=9, height=3, dpi=300)
        ggsave(file='../output/figures/Disturbances_cots_no_label.png',gg,width=9, height=3, dpi=300)
        ggsave(file='../output/figures/Disturbances_cots_no_label.pdf',gg,width=9, height=3, dpi=300)
        png(file='../output/figures/Disturbances_cots_no_label.png',width=9, height=3, res=300, units='in')
        grid.draw(gg)
        dev.off()
        ## ----end
        ## ---- Bleaching
        if (!INCLUDE_GBR) bleaching.sum.all = bleaching.sum.all %>% filter(Location!='Great Barrier Reef') %>% droplevels
        bleaching.dat <- bleaching.sum.all %>% filter(REPORT_YEAR >1985, REPORT_YEAR < (finalYear+1), BLEACHINGcat!='0')
        gbleaching = ggplot(bleaching.dat, aes(y=BLEACHING.p*100, x=REPORT_YEAR+0.3)) +
            geom_bar(stat='identity',position='stack',aes(fill=BLEACHINGcat),width=0.3,show.legend=TRUE)+
            ## geom_bar(data=bleaching.dat %>% filter(BLEACHINGcat %in% c('2','3','4')),stat='identity',position='stack',aes(fill=BLEACHINGcat),width=0.3,show.legend=TRUE)+
            ## geom_bar(data=bleaching.dat %>% filter(BLEACHINGcat %in% c('3','4')),stat='identity',position='stack',aes(fill=BLEACHINGcat),width=0.3,show.legend=FALSE)+
            ## geom_bar(data=bleaching.dat %>% filter(BLEACHINGcat %in% c('4')),stat='identity',position='stack',aes(fill=BLEACHINGcat),width=0.3,show.legend=FALSE)+
                                        #geom_bar(stat='identity',position='stack',aes(fill=BLEACHINGcat),show.legend=FALSE,width=0.3)+
            scale_fill_manual('Bleaching severity', breaks=c(1,2,3,4,5), labels=c(1,2,3,4,5), values=c(scales:::brewer_pal(palette='Reds')(6)[-1])) +
                                        #geom_bar(stat='identity',position='stack',aes(alpha=BLEACHINGcat), fill=brewer.pal(3,'Dark2')[2],show.legend=FALSE,width=0.3) +
            facet_wrap(~Location, nrow=1, strip.position='top',scales='fixed', labeller=labeller(Location=setNames(paste0("", labs, "\n"), labs)))+
                                        #scale_alpha_manual(breaks=c(0,1,2,3,4), values=c(0,0.25,0.5,0.75,1)) +
                                        #scale_y_reverse(expression(atop(Reefs,impacted~('%'))),expand=c(0,0),lim=c(100,0)) +
            scale_y_continuous(expression(Reefs~impacted~('%')),expand=c(0,0),lim=c(0,100)) +
            scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "bottom") + #, limits=c(1985,2021))+
            coord_cartesian(xlim=c(1985, finalYear)) +
            theme_classic() +
            theme(strip.background=element_rect(fill=hues[2], color='black', size=0.5),
                                        #strip.background=element_rect(fill='#00009930', color='black', size=0.5),
                  panel.border=element_rect(fill=NA,color='black'),
                  axis.title.y=element_text(size=rel(1.5),margin=margin(r=0.5,unit='lines')),
                                        #axis.text.x=element_text(size=rel(1.2)),
                  axis.text.x=element_text(size=rel(1.0)), axis.title.x=element_blank(),
                  axis.text.y=element_text(size=rel(1.2)),
                  panel.grid.minor=element_line(size=0.1,color=NA),
                  panel.grid.major=element_line(size=0.1,color='gray70'),
                  panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
                  panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
                  plot.margin=unit(c(2,5,5,0),'pt'),
                  legend.position = 'bottom',
                  strip.text=element_text(margin=margin(t=0.5, b=0.5,unit='lines'),size=20,lineheight=0.5, face='bold',hjust=0.9,vjust=-1))
                                        #strip.text.x=element_blank(), strip.background=element_blank())
        gbleaching
        ## ----end
        ## ---- Cyclones
        if (!INCLUDE_GBR) cyclones.sum.all = cyclones.sum.all %>% filter(Location!='Great Barrier Reef') %>% droplevels
        cyclones.dat = cyclones.sum.all %>% filter(REPORT_YEAR >1985, CYCLONEcat!=0)
        gcyclones = ggplot(cyclones.dat, aes(y=CYCLONE.p, x=REPORT_YEAR+0)) +
            geom_bar(stat='identity',position='stack',aes(fill=CYCLONEcat),width=0.3,show.legend=FALSE)+
            geom_bar(data=cyclones.dat %>% filter(CYCLONEcat %in% c('2','3')),stat='identity',position='stack',aes(fill=CYCLONEcat),width=0.3,show.legend=TRUE)+
            geom_bar(data=cyclones.dat %>% filter(CYCLONEcat %in% c('3')),stat='identity',position='stack',aes(fill=CYCLONEcat),width=0.3,show.legend=FALSE)+
                                        #geom_bar(stat='identity',position='stack',aes(alpha=CYCLONEcat), fill=brewer.pal(3,'Dark2')[3],show.legend=FALSE,width=0.3) +
                                        #geom_bar(stat='identity',position='stack',aes(fill=CYCLONEcat),show.legend=FALSE,width=0.3)+
            scale_fill_manual('Cyclone severity', breaks=c(1,2,3), labels=c(1,2,3), values=c(scales:::brewer_pal(palette='Blues')(4)[-1])) +
            facet_wrap(~Location, nrow=1, strip.position='top',scales='fixed', labeller=labeller(Location=setNames(paste0("", labs, "\n"), labs)))+
                                        #scale_alpha_manual(breaks=c(0,1,2,3,4,5), values=c(0,1/5,2/5,3/5,4/5,5/5)) +

    scale_y_continuous(expression(Reefs~impacted~('%')),expand=c(0,0),lim=c(0,100)) +
                                        #scale_y_reverse(expression(atop(Reefs,impacted~('%'))),expand=c(0,0),lim=c(100,0)) +
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "bottom") + #, limits=c(1985,2021))+
    coord_cartesian(xlim=c(1985, finalYear)) +
    theme_classic() +
    theme(strip.background=element_rect(fill=hues[2], color='black', size=0.5),
                                        #strip.background=element_rect(fill='#00009930', color='black', size=0.5),
          panel.border=element_rect(fill=NA,color='black'),
          axis.title.y=element_text(size=rel(1.5),margin=margin(r=0.5,unit='lines')),
                                        #axis.text.x=element_text(size=rel(1.2)),
          axis.text.x=element_text(size=rel(1.0)), axis.title.x=element_blank(),
          axis.text.y=element_text(size=rel(1.2)),
          panel.grid.minor=element_line(size=0.1,color=NA),
          panel.grid.major=element_line(size=0.1,color='gray70'),
          panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
          panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
          plot.margin=unit(c(2,5,5,0),'pt'),
          legend.position = 'bottom',
          ##strip.text.x=element_blank(), strip.background=element_blank())
          strip.text=element_text(margin=margin(t=0.5, b=0.5,unit='lines'),size=20,lineheight=0.5, face='bold',hjust=0.9,vjust=-1))
        gcyclones

    }


    ## ----end
        ## ---- Put them all together
        ggbleaching = gbleaching + theme(panel.background=element_blank()) + ggtitle('a)')
        ggcyclones = gcyclones + theme(panel.background=element_blank()) + ggtitle('a)')
        ggcots = gcots
        gtcots <- ggplot_gtable(ggplot_build(ggcots + ggtitle('a)')))
        gtbleaching <- ggplot_gtable(ggplot_build(ggbleaching))
        gtcyclones <- ggplot_gtable(ggplot_build(ggcyclones))
        panels <- grepl("panel", gtbleaching$layout$name)

        pp <- c(subset(gtcots$layout, grepl("panel", gtcots$layout$name), se = t:r))

                                        # Overlap panels for second plot on those of the first plot
        gT <- gtable_add_grob(gtcots, gtbleaching$grobs[grepl("panel", gtcots$layout$name)], 
                              pp$t, pp$l, pp$b, pp$l, name='bleaching')
        gT <- gtable_add_grob(gT, gtcyclones$grobs[grepl("panel", gtcots$layout$name)], 
                              pp$t, pp$l, pp$b, pp$l, name='cyclones')
                                        #gT$widths = gg$widths

        grid.draw(gT)

        facets <- grep("strip-t-1-1", gT$layout$name)
        gg <- with(gT$layout[facets,],
                   gtable_add_grob(gT, ggplotGrob(gt1),t=t, l=5, b=b, r=5, name="pic_predator"))
        facets <- grep("strip-t-2-1", gT$layout$name)
        gg <- with(gg$layout[facets,],
                   gtable_add_grob(gg, ggplotGrob(gt2),t=t, l=9, b=b, r=9, name="pic_predator"))
        facets <- grep("strip-t-3-1", gT$layout$name)
        gg <- with(gg$layout[facets,],
                   gtable_add_grob(gg, ggplotGrob(gt3),t=t, l=13, b=b, r=13, name="pic_predator"))
        grid.draw(gg)
        save(gg, file='../data/processed/disturbanceBars-gg.RData')
                                        #grid.arrange(gg,gT, nrow=2,heights=c(2,1),padding=unit(0,'lines'))

        ## ----end
        ## ---- flip the figure around
        gt_2021 = ggplot(fortify(aa.ps), aes(y=lat, x=long, group=group)) +
            geom_blank(aes(x=154,y=-10))+coord_map() +
            geom_blank(aes(x=154,y=-90))+coord_map() +
                                        #geom_blank(aes(x=190,y=-20))+coord_map() +
                                        #geom_blank(aes(x=220,y=-20))+coord_map() +
                                        #geom_blank(aes(x=270,y=-20))+coord_map() +
            geom_polygon(data=fortify(spatial_3Zone), aes(y=lat, x=long),fill=NA,color='black', size=0.2) +
            coord_equal() +
            theme_classic() +theme(panel.background=element_rect(fill=NA),
                                   axis.text.y=element_blank(),
                                   axis.text.x=element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank(),
                                   axis.ticks=element_blank(),
                                   axis.line=element_blank(),
                                   plot.background=element_blank(),
                                   panel.spacing=unit(0,'pt'),
                                   plot.margin=unit(c(0,0,0,0),'pt'))  
        gt1_2021=gt_2021+geom_polygon(data=fortify(spatial_3Zone[1]) , aes(y=lat, x=long),fill=hues[4],color=NA)+
            geom_polygon(data=fortify(spatial_3Zone[1]), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
            geom_polygon(fill='white', color=hues[4])
        gt1_2021
        gt2_2021=gt_2021+geom_polygon(data=fortify(spatial_3Zone[2]) , aes(y=lat, x=long),fill=hues[4],color=NA)+
            geom_polygon(data=fortify(spatial_3Zone[2]), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
            geom_polygon(fill='white', color=hues[4])
        gt2_2021
        gt3_2021=gt_2021+geom_polygon(data=fortify(spatial_3Zone[3]) , aes(y=lat, x=long),fill=hues[4],color=NA)+
            geom_polygon(data=fortify(spatial_3Zone[3]), aes(y=lat, x=long),fill=NA,color='black', size=0.2)  +
            geom_polygon(fill='white', color=hues[4])
        gt3_2021

        labs.shorter <- gsub(' GBR','',labs)
        ggbleaching = gbleaching + theme(panel.background=element_blank(), legend.justification = 'left', legend.direction = 'horizontal') + ggtitle('a)') +
            facet_grid(Location~., scales='fixed', labeller=labeller(Location=setNames(paste0("", labs.shorter, "\n"), labs))) +
            theme(panel.spacing.y = unit(15, 'pt')) +
            guides(fill = guide_legend(title.position = 'top'))

        ggcyclones = gcyclones + theme(panel.background=element_blank(), legend.justification = c(0.6, 0.5), legend.direction = 'horizontal') + ggtitle('a)')+
            facet_grid(Location~., scales='fixed', labeller=labeller(Location=setNames(paste0("", labs.shorter, "\n"), labs)))+
            theme(panel.spacing.y = unit(15, 'pt')) +
            guides(fill = guide_legend(title.position = 'top'))


        ggcots = gcots+
            facet_grid(Location~., scales='fixed', labeller=labeller(Location=setNames(paste0("", labs.shorter, "\n"), labs)))+
            theme(panel.spacing.y = unit(15, 'pt'), legend.justification = "right", legend.direction = 'horizontal') +
            guides(fill = guide_legend(title.position = 'top'))


        gtcots <- ggplot_gtable(ggplot_build(ggcots + ggtitle('a)')))
        gtbleaching <- ggplot_gtable(ggplot_build(ggbleaching))
        gtcyclones <- ggplot_gtable(ggplot_build(ggcyclones))
        panels <- grepl("panel", gtbleaching$layout$name)

        pp <- c(subset(gtcots$layout, grepl("panel", gtcots$layout$name), se = t:r))
        pl <- c(subset(gtcots$layout, grepl("guide-box", gtcots$layout$name), se = t:r))

                                        # Overlap panels for second plot on those of the first plot
        gT <- gtable_add_grob(gtcots, gtbleaching$grobs[grepl("panel", gtcots$layout$name)], 
                              pp$t, pp$l, pp$b, pp$l, name='bleaching')
        gT <- gtable_add_grob(gT, gtcyclones$grobs[grepl("panel", gtcots$layout$name)], 
                              pp$t, pp$l, pp$b, pp$l, name='cyclones')
                                        #gT$widths = gg$widths

        gT <- gtable_add_grob(gT, gtcyclones$grobs[grepl("guide-box", gtcots$layout$name)], 
                              pl$t, pl$l, pl$b, pl$r, name='cyclones-guide')
        gT <- gtable_add_grob(gT, gtbleaching$grobs[grepl("guide-box", gtcots$layout$name)], 
                              pl$t, pl$l, pl$b, pl$r, name='bleaching-guide')
        
        grid.draw(gT)

        facets <- grep("strip-r-1-1", gT$layout$name)
        gg <- with(gT$layout[facets,],
                   ## gtable_add_grob(gT, ggplotGrob(gt1),t=t, l=5, b=b, r=5, name="pic_predator"))
                   gtable_add_grob(gT, ggplotGrob(gt1_2021),t=5, l=6, b=7, r=7, name="pic_predator"))
        facets <- grep("strip-r-2-1", gT$layout$name)
        gg <- with(gg$layout[facets,],
                   gtable_add_grob(gg, ggplotGrob(gt2_2021),t=9, l=6, b=9, r=7, name="pic_predator"))
        facets <- grep("strip-r-3-1", gT$layout$name)
        gg_2021 <- with(gg$layout[facets,],
                        gtable_add_grob(gg, ggplotGrob(gt3_2021),t=11, l=6, b=11, r=7, name="pic_predator"))
        grid.draw(gg_2021)

        save(gg_2021, file='../data/processed/disturbanceBars-gg-nolabel_transposed.RData')
        ## ggsave(file='output/figures/Disturbances_all_no_label_transposed.pdf',grid.draw(gg_2021))
        pdf(file='../output/figures/Disturbances_all_no_label_transposed.pdf', width = 7, height = 7*1.05)
        grid.draw(gg_2021)
        dev.off()
        ## The following is not working any more 174mm - could try -units PixelsPerCentimeter
                                        #system('convert -resize 174mmx output/figures/Disturbances_all_no_label.pdf output/figures/Disturbances_all_no_label.pdf')
        ## ggsave(file='output/figures/Disturbances_all_no_label_transposed.png',grid.draw(gg),width=9, height=3, dpi=300)
        png(file='../output/figures/Disturbances_all_no_label_transposed.png', width=2400, height=2400, res=300)
        grid.draw(gg_2021)
        dev.off()

        ## ----end
    }
    ## ----end

}    
## ----end


