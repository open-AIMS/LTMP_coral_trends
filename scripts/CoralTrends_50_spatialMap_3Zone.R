source('CoralTrends_functions.R')
CoralTrends_checkPackages()

##############################
## Load the 3 Zone polygons ##
##############################
load(file='../data/spatial/whagbr.RData')
load(file='../data/spatial/whagbr.n.RData')
load(file='../data/spatial/whagbr.c.RData')
load(file='../data/spatial/whagbr.s.RData')
load('../data/spatial/qld.RData')

load('../data/processed/manta.sum.RData')

#######################################################################
## Summarize to spatial data for each site (marginalizing over time) ##
#######################################################################
manta.sites = manta.sum %>% ungroup %>%
    group_by(REEF_NAME)  %>%
    summarize(Latitude=mean(Latitude,na.rm=TRUE),
              Longitude=mean(Longitude,na.rm=TRUE),
              N=n(),
              Tows=mean(Tows,na.rm=TRUE)) %>%
    ungroup

#####################################################################
## Generate a site map bubble plot in which the size of bubbles is ##
## proportional to the number of visits                            ##
#####################################################################
gm1=ggplot(fortify(qld), aes(y=lat, x=long)) +
    geom_polygon(aes(group=group), fill='grey', color='grey40') +
    geom_polygon(data=fortify(whagbr.n),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.c),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.s),aes(group=group), color='black', fill=NA) +
    geom_point(data=manta.sites, aes(y=Latitude,x=Longitude, size=N/5), fill=NA,shape=21,color=NA) +
    geom_point(data=manta.sites, aes(y=Latitude,x=Longitude), size=manta.sites$N/5,alpha=0.3,fill='blue',shape=21,color='black') +
    coord_equal(xlim=bbox(whagbr)[1,],ylim=bbox(whagbr)[2,]) +
    annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf, fill=NA,color='black')+
    scale_x_continuous('') + scale_y_continuous('')+
    scale_fill_brewer('Year',type='qual',palette='Set1')+
    #scale_size_area('Number of visits') +
#    scale_size_identity('Number of visits', breaks=c(1,10,20),labels=c(1,10,20),guide='legend') +
    scale_size_continuous('Number of visits', breaks=1:6, labels=1:6*5, limits=c(0,30)/5)+
    annotate(geom='text',x=145.5,y=-13, label='Northern GBR', hjust=0) +
    annotate(geom='text',x=147.7,y=-17.5, label='Central GBR', hjust=0) +
    annotate(geom='text',x=153.7,y=-20.2, label='Southern\nGBR', hjust=1) +
    theme_classic() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
                            legend.position = c(0.95,0.95),legend.justification=c(1,1))+
    guides(fill=guide_legend(override.aes = list(size=5)),
           size=guide_legend(override.aes=list(alpha=0.3,fill='blue',color='black')))

#####################################################################
## Generate a site map bubble plot in which the size of bubbles is ##
## proportional to the total number of manta tows                  ##
#####################################################################
gm2=ggplot(fortify(qld), aes(y=lat, x=long)) +
    geom_polygon(aes(group=group), fill='grey', color='grey40') +
    geom_polygon(data=fortify(whagbr.n),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.c),aes(group=group), color='black', fill=NA) +
    geom_polygon(data=fortify(whagbr.s),aes(group=group), color='black', fill=NA) +
    geom_point(data=manta.sites, aes(y=Latitude,x=Longitude, size=Tows/20), fill=NA,shape=21,color=NA) +
    geom_point(data=manta.sites, aes(y=Latitude,x=Longitude), size=manta.sites$Tows/20,alpha=0.3,fill='blue',shape=21,color='black') +
    coord_equal(xlim=bbox(whagbr)[1,],ylim=bbox(whagbr)[2,]) +
    annotate(geom='rect', xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf, fill=NA,color='black')+
    scale_x_continuous('') + scale_y_continuous('')+
    scale_fill_brewer('Year',type='qual',palette='Set1')+
    #scale_size_area('Number of visits') +
#    scale_size_identity('Number of visits', breaks=c(1,10,20),labels=c(1,10,20),guide='legend') +
    scale_size_continuous('Number of Tows', breaks=c(1,10,20,40,80,160)/20, labels=c(1,10,20,40,80,160), limits=c(0,160)/20)+
    annotate(geom='text',x=145.5,y=-13, label='Northern GBR', hjust=0) +
    annotate(geom='text',x=147.7,y=-17.5, label='Central GBR', hjust=0) +
    annotate(geom='text',x=153.7,y=-20.2, label='Southern\nGBR', hjust=1) +
    theme_classic() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
                            legend.position = c(0.95,0.95),legend.justification=c(1,1))+
    guides(fill=guide_legend(override.aes = list(size=5)),
           size=guide_legend(override.aes=list(alpha=0.3,fill='blue',color='black')))

ggsave('../output/figures/MapOfSites.png', grid.arrange(gm1,gm2,nrow=1),width=10, height=5, dpi=300)
ggsave('../output/figures/MapOfSites.pdf', grid.arrange(gm1,gm2,nrow=1),width=10, height=5, dpi=300)
