source('CoralTrends_functions.R')
CoralTrends_checkPackages()

load('../data/processed/manta.sum.RData')


manta.heat = manta.sum %>% group_by(REEF_ID, REPORT_YEAR) %>% summarise(Latitude=mean(Latitude,na.rm=TRUE), Freq=length(unique(Latitude)))
manta.heat = manta.heat %>% ungroup %>% arrange(-Latitude,REPORT_YEAR)
manta.heat$Location=factor(CoralTrends_calc3ZoneLocations(as.numeric(as.character(manta.heat$Latitude))), levels=c('Northern','Central','Southern'))
manta.heat$REPORT_YEAR <- factor(as.character(manta.heat$REPORT_YEAR))
manta.heat$REEF_ID <- factor(as.character(manta.heat$REEF_ID), levels=rev(unique(as.character(manta.heat$REEF_ID))))
g=ggplot(manta.heat, aes(y=REEF_ID, x=REPORT_YEAR, fill=as.factor(Freq))) +
    geom_tile() +
    scale_fill_manual('Surveys', breaks=c(1), values=c('orange'),guide=FALSE) +
        scale_y_discrete('Latitude')+
            scale_x_discrete('')+
                    facet_grid(Location~., scales='free',space='free')+
            theme_classic()+
        theme(axis.text.y=element_text(size=5),
              axis.ticks.y=element_blank(),
              panel.background=element_rect(fill=NA,color='black'),
              strip.background=element_blank(),
              strip.text=element_text(size=12))
ggsave('../output/figures/TemporalHeatMap_3Zone.png', g,width=13, height=15, dpi=300)
ggsave('../output/figures/TemporalHeatMap_3Zone.pdf', g,width=13, height=15, dpi=300)
