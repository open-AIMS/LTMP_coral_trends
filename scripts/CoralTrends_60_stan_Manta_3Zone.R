source('CoralTrends_functions.R')
CoralTrends_checkPackages()

load('../data/processed/manta.sum.RData')
load('../data/processed/manta.tow.RData')

## Generate a list of reefs we are using to help accumulate other
## sources of associated data
all.reefs = manta.sum %>%
    dplyr:::select(P_CODE.mod,REEF_NAME,REEF_ID,Latitude,Longitude) %>%
    group_by(REEF_NAME,REEF_ID) %>%
    summarize_at(vars(Latitude,Longitude), funs(mean)) %>%
    as.data.frame
write.csv(all.reefs,file='../data/all.reefs_3Zone.csv', quote=FALSE, row.names=FALSE)

## Genuine stan cannot handle proportional data for binomial families
## (particularly when weights are applied). A work-around is to
## multiple the proportion by the weights and convert this into an integer
dat.all = manta.sum %>%
    mutate(Location=Region) %>%
    dplyr:::select(Cover, REEF_NAME, Tows,P_CODE.mod,Location,REPORT_YEAR) %>%
    mutate(Year=factor(REPORT_YEAR), N=length(unique(REEF_NAME))) %>% ungroup() %>%
    mutate(Cvr1 = as.integer(as.vector(Cover) * Tows), Cvr0 = Tows - Cvr1)

##original, stan_glmer beta, INLA_tow beta scaled, INLA_reef binomial, INLA_reef beta
##BRMS beta vanilla, BRMS beta disp, MGCV beta , MGCV ordinal, CLMM, BRMS_reef beta,
##
models <- c('BRMS beta vanilla', 'BRMS beta disp', 'glmmTMB beta vanilla', 'glmmTMB beta disp', 'BRMS ordinal')
models <- 'glmmTMB_tow beta disp'
## Fit stan models===================================================================

## GBR
## ---- GBR
{
  dat.all.gbr = dat.all %>% droplevels
  save(dat.all.gbr, file='../data/modelled/dat.all.gbr.RData')

  ## In 2021, we used the glmmTMB beta disp model

  ## ---- Gbr.Data
  ## Tow level data
  manta.tow.gbr = manta.tow %>%
    droplevels %>%
    mutate(oLIVE_CORAL=factor(LIVE_CORAL,
                              levels=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U'),
                              ordered=TRUE),
           nREEF_NAME=as.numeric(as.factor(REEF_NAME)),
           nLIVE_CORAL=as.numeric(oLIVE_CORAL),
           REEF_YEAR = interaction(REEF_NAME, Year)
           )
  save(manta.tow.gbr, file='../data/modelled/manta.tow.gbr.RData')
  ## ----end

  ## Raw cells ----------------------------------------------------------------

  ## ---- Raw cells
  dat.all.gbr.cellmeans <- cellMeansRaw(dat.all.gbr)
  rawAdd <- ggproto_Raw(dat.all.gbr.cellmeans)

  manta.tow.gbr.cellmeans <- cellMeansRaw(manta.tow %>%
                                          group_by(P_CODE.mod, REEF_NAME, Year) %>%
                                          summarise(Cover=mean(Cover), Tows=length(unique(TOW_SEQ_NO))) %>%
                                          ungroup)
  ## ----end

  manta.tow.gbr = manta.tow %>%
    mutate(oLIVE_CORAL=factor(LIVE_CORAL,
                              levels=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U'),
                              ordered=TRUE),
           nREEF_NAME=as.numeric(as.factor(REEF_NAME)),
           nLIVE_CORAL=as.numeric(oLIVE_CORAL),
           REEF_YEAR = interaction(REEF_NAME, Year)
           )

  ## ---- GBR.glmmTMB.tow.beta disp
  if ('glmmTMB_tow beta disp' %in% models) {
    mod.gbr_glmmTMB.beta.disp <- glmmTMB(Cover ~ Year + (1|REEF_NAME/REEF_YEAR),
                                         dispformula = ~Year,
                                         data=manta.tow.gbr,
                                         ## weights=dat.all.gbr$Tows,
                                         family=beta_family())
    dat.gbr_glmmTMB.beta.disp = emmeans(mod.gbr_glmmTMB.beta.disp, ~Year, type='response') %>%
      as.data.frame()
    ## DHARMa::simulateResiduals(mod.gbr_glmmTMB.beta.disp, plot=TRUE)
    ## performance::check_model(mod.gbr_glmmTMB.beta.disp)
    save(mod.gbr_glmmTMB.beta.disp, dat.gbr_glmmTMB.beta.disp, file=paste0('../data/modelled/mod.gbr_glmmTMB.beta.disp.RData'))

  }
  ## ----end
  ## ---- summarise model
  load(file=paste0('../data/modelled/mod.gbr_glmmTMB.beta.disp.RData'))
  summary(mod.gbr_glmmTMB.beta.disp)
  head(dat.gbr_glmmTMB.beta.disp)
  g1 <- dat.gbr_glmmTMB.beta.disp %>%
    ggplot(aes(y=response, x=as.numeric(as.character(Year)))) +
    geom_ribbon(aes(ymin=lower.CL, ymax=upper.CL), fill='blue', alpha=0.3) +
    geom_line(color='blue') +
    scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
    scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) +
    scale_color_discrete('Raw data aggregate') +
    rawAdd +
    theme_classic() +
    theme(axis.title.x=element_blank(),
          legend.position=c(0.01,0.01), legend.justification=c(0,0),
          panel.grid.minor=element_line(),
          panel.grid.major=element_line()) +
    ggtitle('GBR') +
    guides(color=guide_legend(nrow=2, byrow=TRUE))
  g1
  dev.off()

  ## ----end
  ## ---- Extract reef-level predictions (glmmTMB.tow.beta disp)
  load(file=paste0('../data/modelled/mod.gbr_glmmTMB.beta.disp.RData'))

  ## coef(mod.gbr_glmmTMB.beta.disp)$cond$REEF_NAME['AGINCOURT REEFS (NO 1)',]
  ## plogis(-1.068388 - 0.3291372)
  ## ndata = data.frame(Year=factor(1986:2021))
  ## Xmat <- model.matrix(~Year, ndata)
  ## coefs <- coef(mod.gbr_glmmTMB.beta.disp)$cond$REEF_NAME['AGINCOURT REEFS (NO 1)',] %>% as.matrix() %>% as.vector()
  ## ndata = cbind(ndata, Fit=as.vector(plogis(coefs %*% t(Xmat))))
  ## g1 + geom_line(data=ndata, aes(y=Fit, x=as.numeric(as.character(Year))), color='green')


  ## coefs <- coef(mod.gbr_glmmTMB.beta.disp)$cond$REEF_YEAR %>% filter(str_detect(rownames(.), 'AGINCOURT REEFS \\(NO 1\\)'))
  ## rownames(coefs)
  ## #pull(`(Intercept)`)
  ## coefs[1,]
  ## plogis(as.matrix(coefs[1,]) %*% t(Xmat))


  ## coef(mod.gbr_glmmTMB.beta.disp)$cond$REEF_YEAR['AGINCOURT REEFS (NO 1).1989:AGINCOURT REEFS (NO 1)',]
  ## plogis(-1.514835 -0.3291372)




  ## plogis((-1.068388 + 0.3291372) + (-1.514835 +0.3291372))
  ## plogis((-1.068388 + 0.3291372) + (-1.514835 +0.3291372))


  ## aa =  fixef(mod.gbr_glmmTMB.beta.disp)[[1]]
  ## ndata = data.frame(Year=factor(1986:2021))
  ## Xmat <- model.matrix(~Year, ndata)

  ## ab1 = coef(mod.gbr_glmmTMB.beta.disp)$cond$REEF_YEAR %>%
  ##                                    rownames_to_column('R') %>%
  ##                                    filter(str_detect(R,'AGINCOURT REEFS \\(NO 1\\)'))
  ## ab1 %>% dim
  ## ab1 %>% head

  ## ab2 = coef(mod.gbr_glmmTMB.beta.disp)$cond$REEF_NAME %>%
  ##                                    rownames_to_column('R') %>%
  ##                                    filter(str_detect(R,'AGINCOURT REEFS \\(NO 1\\)'))
  ## ab2 %>% dim
  ## ab2 %>% head

  ## sweep(
  ##     x=as.matrix(ab1[,-1]),
  ##     MARGIN=2,
  ##     STATS=as.vector(as.matrix(ab2[-1])),
  ##     FUN='+'
  ##     )
  reefs <- manta.tow.gbr %>% pull(REEF_NAME) %>% unique
  fit.all.reefs <- vector('list', length(reefs))
  names(fit.all.reefs) <- reefs
  for (r in reefs) {
    print(r)
    raw.sum<-manta.tow.gbr %>% filter(REEF_NAME==r) %>%
      group_by(Year) %>%
      summarise(Cover=mean(Cover))

                                        #r='AGINCOURT REEFS (NO 1)'
    a0 <- fixef(mod.gbr_glmmTMB.beta.disp)[[1]][1]
    a1 <- fixef(mod.gbr_glmmTMB.beta.disp)[[1]][-1]
    a1 <- data.frame(Slope=a1) %>% rownames_to_column('Year') %>% mutate(Year=gsub('Year','',Year))
    a0 <- cbind(Intercept=a0, a1)

    a<-ranef(mod.gbr_glmmTMB.beta.disp) %>% `[[`(1) %>% `[[`(1) %>%
      as.data.frame() %>%
      rownames_to_column('R') %>%
      mutate(REEF_NAME=gsub('(.*)\\.[0-9]{4}.*', '\\1', R),
             Year=gsub('.*\\.([0-9]{4}).*', '\\1', R)) %>%
      dplyr::select(-R) %>%
      dplyr::rename(rand.Intercept=`(Intercept)`) %>%
      filter(REEF_NAME==r)
    b<-ranef(mod.gbr_glmmTMB.beta.disp) %>% `[[`(1) %>% `[[`(2) %>%
      as.data.frame %>%
      rownames_to_column('REEF_NAME') %>%
      dplyr::rename(rand.slope=`(Intercept)`) %>%
      filter(REEF_NAME==r)
    a1 <- a %>% left_join(b)
    a3 <- a1 %>% left_join(a0) %>%
      mutate(Cover=plogis(Intercept+rand.Intercept+Slope+rand.slope)) %>%
      left_join(raw.sum %>% dplyr::rename(Raw=Cover))
    ## a3 %>% head
    fit.all.reefs[[r]] <- a3

    ## g1 <- ggplot() +
    ##     geom_line(data=a3, aes(y=Cover, x=as.numeric(as.character(Year))), color='red') +
    ##     geom_line(data=raw.sum, aes(y=Cover, x=as.numeric(as.character(Year))), color='blue')
    ## g1
  }
  fit.all.reefs <- do.call('rbind',fit.all.reefs)

  fit.all.reefs %>% filter(REEF_NAME=='BROOMFIELD REEF') %>%
    ggplot() +
    geom_line(aes(y=Cover, x=as.numeric(as.character(Year))),color='blue') +
    geom_line(aes(y=Raw, x=as.numeric(as.character(Year))), color='red')



  ## Using predict function is very slow and seems only to permit one prediction at a time with newdata?
  ## newdata <- data.frame(Year=1989:2021,
  ##                       REEF_NAME='AGINCOURT REEFS (NO 1)') %>%
  ##     mutate(REEF_YEAR=interaction(REEF_NAME,Year))

  ## predict(mod.gbr_glmmTMB.beta.disp, newdata=data.frame(Year=1990:1992, REEF_NAME=NA, REEF_YEAR=NA), re.form=~0, type='response')
  ## predict(mod.gbr_glmmTMB.beta.disp, newdata=newdata, re.form=NULL, type='response')
  ## predict(mod.gbr_glmmTMB.beta.disp, newdata=newdata[1,], type='response', se.fit=TRUE)



  ## predict(mod.gbr_glmmTMB.beta.disp, newdata=data.frame(Year=1990, REEF_NAME=r, REEF_YEAR=interaction(1990,r)))
  ## ----end
}
## ----end

## Northern
## ---- Northern
{
  ## ---- Northern.Data
  dat.all.northern = dat.all %>%
    filter(Location=='Northern GBR') %>%
    droplevels %>%
    mutate(P_CODE.mod=factor(ifelse(is.na(P_CODE.mod),'Other',P_CODE.mod))) %>%
    group_by(REEF_NAME) %>%
    mutate(W=mean(Tows, na.rm=TRUE)) %>%
    ungroup %>%
    mutate(W1=W/sum(W)) %>%
    group_by(Year) %>%
    mutate(W2=Tows/sum(Tows)) %>%
    ungroup
  save(dat.all.northern, file='../data/modelled/dat.all.northern.RData')
  ## Tow level data
  manta.tow.northern = manta.tow %>%
    filter(Region=='Northern GBR') %>%
    droplevels %>%
    mutate(oLIVE_CORAL=factor(LIVE_CORAL,
                              levels=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U'),
                              ordered=TRUE),
           nREEF_NAME=as.numeric(as.factor(REEF_NAME)),
           nLIVE_CORAL=as.numeric(oLIVE_CORAL),
           REEF_YEAR = interaction(REEF_NAME, Year)
           )
  save(manta.tow.northern, file='../data/modelled/manta.tow.northern.RData')
  ## ----end
  ## Raw cells
  ## ---- Northern.Raw cells
  dat.all.northern.cellmeans <- cellMeansRaw(dat.all.northern)
  rawAdd <- ggproto_Raw(dat.all.northern.cellmeans)

  manta.tow.northern.cellmeans <- cellMeansRaw(manta.tow %>%
                                               filter(Region=='Northern GBR') %>%
                                               droplevels %>%
                                               group_by(P_CODE.mod, REEF_NAME, Year) %>%
                                               summarise(Cover=mean(Cover), Tows=length(unique(TOW_SEQ_NO))) %>%
                                               ungroup)
  ## ----end
  ## ---- Northern.BRMS.tow.beta disp
  if ('BRMS beta disp' %in% models) {
    mod.northern_brms.beta.disp <- brm(bf(Cover ~ Year + (1|REEF_NAME), phi~0+Year),
                                       data=manta.tow.northern,
                                       family=Beta(link='logit'),
                                       iter=1e4,
                                       warmup=5e3,
                                       thin=5,
                                       chains=4, cores=4,
                                       prior = prior(normal(0, 3), class = "b") +
                                         prior(normal(0, 3), class = "Intercept") +
                                         prior(gamma(2, 1), class = "sd") #+
                                       ## prior(gamma(2, 1), class = "phi")
                                       )
    dat.northern_brms.beta.disp = emmeans(mod.northern_brms.beta.disp, ~Year, type='response') %>%
      as.data.frame()
    save(mod.northern_brms.beta.disp, dat.northern_brms.beta.disp, file=paste0('../data/modelled/mod.northern_brms.beta.disp.RData'))
    rm(list=c('dat.northern_brms.beta.disp', 'mod.northern_brms.beta.disp'))
    gc()
  }
  ## ----end
  ## ---- Northern.glmmTMB.tow.beta disp
  if ('glmmTMB_tow beta disp' %in% models) {
    mod.northern_glmmTMB.beta.disp <- glmmTMB(Cover ~ Year + (1|REEF_NAME/REEF_YEAR),
                                              dispformula = ~Year,
                                              data=manta.tow.northern,
                                              ## weights=dat.all.northern$Tows,
                                              family=beta_family())
    dat.northern_glmmTMB.beta.disp = emmeans(mod.northern_glmmTMB.beta.disp, ~Year, type='response') %>%
      as.data.frame()
    ## DHARMa::simulateResiduals(mod.northern_glmmTMB.beta.disp, plot=TRUE)
    ## performance::check_model(mod.northern_glmmTMB.beta.disp)
    save(mod.northern_glmmTMB.beta.disp, dat.northern_glmmTMB.beta.disp, file=paste0('../data/modelled/mod.northern_glmmTMB.beta.disp.RData'))
    ## Lets also capture the full posteriors for years

  }
  ## ----end
  ## ---- Northern.glmmTMB.tow.beta disp random.effects
  if ('glmmTMB_tow beta disp rs' %in% models) {
    nt <- parallel::detectCores()
    mod.northern_glmmTMB.beta.disp.rs <- glmmTMB(Cover ~ Year + (Year|REEF_NAME),
                                                 dispformula = ~Year,
                                                 data=manta.tow.northern,
                                                 ## weights=dat.all.northern$Tows,
                                                 family=beta_family(),
                                                 control = glmmTMBControl(parallel=nt))
    dat.northern_glmmTMB.beta.disp.rs = emmeans(mod.northern_glmmTMB.beta.disp.rs, ~Year, type='response') %>%
      as.data.frame()
    DHARMa::simulateResiduals(mod.northern_glmmTMB.beta.disp.rs, plot=TRUE)
    performance::check_model(mod.northern_glmmTMB.beta.disp.rs)
    save(mod.northern_glmmTMB.beta.disp.rs, dat.northern_glmmTMB.beta.disp.rs, file=paste0('../data/modelled/mod.northern_glmmTMB.beta.disp.rs.RData'))
  }
  ## ----end
  ## Compare the models
  ## ---- Northern Compare models
  if(1==2) {

    load(file='../data/modelled/dat.northern.RData')
    dat.northern.original <- dat.northern
    load(file='../data/modelled/mod.northern.RData')
    load(file='../data/modelled/mod.northern_inla_beta.RData')
    load(file='../data/modelled/newdata.northern_inla_beta.RData')
    newdata.northern_beta <- newdata.northern
    load(file='../data/modelled/mod.northern_inla_binomial.RData')
    load(file='../data/modelled/newdata.northern_inla_binomial.RData')
    newdata.northern_binomial <- newdata.northern
    load(file='../data/modelled/dat.northern_inla_tow.RData')
    load(file='../data/modelled/dat.northern_mgcv.RData')
    load(file='../data/modelled/dat.northern_brms.cumulative.RData')
    load(file='../data/modelled/dat.northern_clmm.RData')
    load(file='../data/modelled/dat.northern_brms.beta.RData')
    load(file='../data/modelled/dat.northern_glmmTMB.RData')
    load(file='../data/modelled/dat.northern_mgcv.beta.RData')

    ## original
    g1 <- dat.northern.original %>%
      ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
      geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
      geom_line(color='blue') +
      scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
      scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) +
      scale_color_discrete('Raw data aggregate') +
      rawAdd +
      theme_classic() +
      theme(axis.title.x=element_blank(),
            legend.position=c(0.01,0.01), legend.justification=c(0,0),
            panel.grid.minor=element_line(),
            panel.grid.major=element_line()) +
      ggtitle('Original (stan binomial reef level)') +
      guides(color=guide_legend(nrow=2, byrow=TRUE))
    g1

    ## inla (reef level binomial)
    g2 <- newdata.northern_binomial %>%
      ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
      geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
      geom_line(color='blue') +
      scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
      scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) +
      scale_color_discrete('Raw data aggregate') +
      rawAdd +
      theme_classic() +
      theme(axis.title.x=element_blank(),
            legend.position=c(0.01,0.01), legend.justification=c(0,0),
            panel.grid.minor=element_line(),
            panel.grid.major=element_line()) +
      ggtitle('INLA reef level binomial') +
      guides(color=guide_legend(nrow=2, byrow=TRUE))
    g2

    ## inla (reef level beta)
    g3 <- newdata.northern_beta %>%
      ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
      geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
      geom_line(color='blue') +
      scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
      scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) +
      scale_color_discrete('Raw data aggregate') +
      rawAdd +
      theme_classic() +
      theme(axis.title.x=element_blank(),
            legend.position=c(0.01,0.01), legend.justification=c(0,0),
            panel.grid.minor=element_line(),
            panel.grid.major=element_line()) +
      ggtitle('INLA reef level beta') +
      guides(color=guide_legend(nrow=2, byrow=TRUE))
    g3

    ## inla (tow level beta)
    g4 <- dat.northern %>%
      ggplot(aes(y=mean, x=as.numeric(as.character(Year)))) +
      geom_ribbon(aes(ymin=lower, ymax=upper), fill='blue', alpha=0.3) +
      geom_line(color='blue') +
      scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
      scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) +
      scale_color_discrete('Raw data aggregate') +
      rawAdd +
      theme_classic() +
      theme(axis.title.x=element_blank(),
            legend.position=c(0.01,0.01), legend.justification=c(0,0),
            panel.grid.minor=element_line(),
            panel.grid.major=element_line()) +
      ggtitle('INLA tow level beta') +
      guides(color=guide_legend(nrow=2, byrow=TRUE))
    g4

    ## mgcv ordinal
    g5 <- dat.northern.mgcv %>%
      ggplot(aes(y=Mean, x=as.numeric(as.character(Year)))) +
      geom_line(color='blue') +
      scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
      scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) +
      scale_color_discrete('Raw data aggregate') +
      rawAdd +
      theme_classic() +
      theme(axis.title.x=element_blank(),
            legend.position=c(0.01,0.01), legend.justification=c(0,0),
            panel.grid.minor=element_line(),
            panel.grid.major=element_line()) +
      ggtitle('mgcv ordinal') +
      guides(color=guide_legend(nrow=2, byrow=TRUE))

    ## brms ordinal
    g6 <- dat.northern_brms.cumulative %>%
      ggplot(aes(y=estimate, x=as.numeric(as.character(Year)))) +
      geom_ribbon(aes(ymin=conf.low, ymax=conf.high), fill='blue', alpha=0.3) +
      geom_line(color='blue') +
      scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
      scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) +
      scale_color_discrete('Raw data aggregate') +
      rawAdd +
      theme_classic() +
      theme(axis.title.x=element_blank(),
            legend.position=c(0.01,0.01), legend.justification=c(0,0),
            panel.grid.minor=element_line(),
            panel.grid.major=element_line()) +
      ggtitle('brms ordinal') +
      guides(color=guide_legend(nrow=2, byrow=TRUE))
    g6

    ## clmm ordinal
    g7 <- dat.northern_clmm %>%
      ggplot(aes(y=Mean, x=as.numeric(as.character(Year)))) +
      geom_line(color='blue') +
      scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
      scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) +
      scale_color_discrete('Raw data aggregate') +
      rawAdd +
      theme_classic() +
      theme(axis.title.x=element_blank(),
            legend.position=c(0.01,0.01), legend.justification=c(0,0),
            panel.grid.minor=element_line(),
            panel.grid.major=element_line()) +
      ggtitle('clmm ordinal') +
      guides(color=guide_legend(nrow=2, byrow=TRUE))
    g7

    ## brms beta
    g8 <- dat.northern_brms.beta %>%
      ggplot(aes(y=response, x=as.numeric(as.character(Year)))) +
      geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD), fill='blue', alpha=0.3) +
      geom_line(color='blue') +
      scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
      scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) +
      scale_color_discrete('Raw data aggregate') +
      rawAdd +
      theme_classic() +
      theme(axis.title.x=element_blank(),
            legend.position=c(0.01,0.01), legend.justification=c(0,0),
            panel.grid.minor=element_line(),
            panel.grid.major=element_line()) +
      ggtitle('brms beta') +
      guides(color=guide_legend(nrow=2, byrow=TRUE))
    g8

    library(patchwork)
    g1 + g2 + g3 + g4 + g5 + g6 + g7 + g8


    dat.northern_brms.beta %>%
      ggplot(aes(y=response, x=as.numeric(as.character(Year)))) +
      ## geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD, fill='BRMS beta'), alpha=0.3) +
      geom_line(aes(color='BRMS beta')) +
      ## geom_ribbon(data=dat.northern_brms.cumulative, aes(y=estimate, ymin=conf.low, ymax=conf.high, fill='BRMS ordinal'), alpha=0.3) +
      geom_line(data=dat.northern_brms.cumulative, aes(y=estimate, color='BRMS ordinal')) +
      ## geom_ribbon(data=newdata.northern_beta, aes(y=mean, ymin=lower, ymax=upper, fill='INLA beta'), alpha=0.3) +
      geom_line(data=newdata.northern_beta, aes(y=mean, color='INLA beta')) +
      geom_line(data=dat.northern_clmm, aes(y=Mean, color='clmm ordinal')) +
      geom_line(data=dat.northern.mgcv, aes(y=Mean, color='mgcv ordinal')) +
      geom_line(data=dat.northern_glmmTMB, aes(y=response, color='glmmTMB beta')) +
      geom_line(data=dat.northern_mgcv.beta, aes(y=response, color='mgcv beta')) +
      scale_y_continuous('Coral cover (%)', labels = function(x) x*100) +
      scale_x_continuous('', breaks=function(x) round(seq(from=x[1], to=x[2], by=2),0)) +
      scale_color_brewer('Model', type='qual') +
      scale_fill_discrete('Model') +
      ## rawAdd +
      theme_classic() +
      theme(axis.title.x=element_blank(),
            legend.position=c(0.01,0.01), legend.justification=c(0,0),
            panel.grid.minor=element_line(),
            panel.grid.major=element_line()) +
      ggtitle('brms beta') +
      guides(color=guide_legend(nrow=2, byrow=TRUE))


    rm(list=c('dat.northern','mod.northern','mod.northern_inla_beta','newdata.northern','newdata.northern_beta', 'mod_northern_inla_binomial','newdata.northern_binomial'))
  }
  ## ----end

}
## ----end


## Central
## ---- Central
## ---- Central.Data
dat.all.central = dat.all %>%
    filter(Location=='Central GBR') %>%
    droplevels %>%
    mutate(P_CODE.mod=factor(ifelse(is.na(P_CODE.mod),'Other',P_CODE.mod))) %>%
    group_by(REEF_NAME) %>%
    mutate(W=mean(Tows, na.rm=TRUE)) %>%
    ungroup %>%
    mutate(W1=W/sum(W)) %>%
    group_by(Year) %>%
    mutate(W2=Tows/sum(Tows)) %>%
    ungroup
save(dat.all.central, file='../data/modelled/dat.all.central.RData')

## Tow level data
manta.tow.central = manta.tow %>%
    filter(Region=='Central GBR') %>%
    droplevels %>%
    mutate(oLIVE_CORAL=factor(LIVE_CORAL,
                              levels=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U'),
                              ordered=TRUE),
           nREEF_NAME=as.numeric(as.factor(REEF_NAME)),
           nLIVE_CORAL=as.numeric(oLIVE_CORAL),
           REEF_YEAR = interaction(REEF_NAME, Year)
           )
save(manta.tow.central, file='../data/modelled/manta.tow.central.RData')
## ----end
## ---- Central.Raw cells
dat.all.central.cellmeans <- cellMeansRaw(dat.all.central)
rawAdd <- ggproto_Raw(dat.all.central.cellmeans)

manta.tow.central.cellmeans <- cellMeansRaw(manta.tow %>%
                                             filter(Region=='Central GBR') %>%
                                             droplevels %>%
                                             group_by(P_CODE.mod, REEF_NAME, Year) %>%
                                             summarise(Cover=mean(Cover), Tows=length(unique(TOW_SEQ_NO))) %>%
                                             ungroup)
## ----end
## ---- Central.glmmTMB.tow.beta disp
if ('glmmTMB_tow beta disp' %in% models) {
    mod.central_glmmTMB.beta.disp <- glmmTMB(Cover ~ Year + (1|REEF_NAME/REEF_YEAR),
                                     dispformula = ~Year,
                                     data=manta.tow.central,
                                     ## weights=dat.all.central$Tows,
                                     family=beta_family())
    dat.central_glmmTMB.beta.disp = emmeans(mod.central_glmmTMB.beta.disp, ~Year, type='response') %>%
        as.data.frame()
    ## DHARMa::simulateResiduals(mod.central_glmmTMB.beta.disp, plot=TRUE)
    ## performance::check_model(mod.central_glmmTMB.beta.disp)
    save(mod.central_glmmTMB.beta.disp, dat.central_glmmTMB.beta.disp, file=paste0('../data/modelled/mod.central_glmmTMB.beta.disp.RData'))
}
## ----end

## ----end


## Southern
## ---- Southern
## ---- Southern.Data
dat.all.southern = dat.all %>%
    filter(Location=='Southern GBR') %>%
    droplevels %>%
    mutate(P_CODE.mod=factor(ifelse(is.na(P_CODE.mod),'Other',P_CODE.mod))) %>%
    group_by(REEF_NAME) %>%
    mutate(W=mean(Tows, na.rm=TRUE)) %>%
    ungroup %>%
    mutate(W1=W/sum(W)) %>%
    group_by(Year) %>%
    mutate(W2=Tows/sum(Tows)) %>%
    ungroup
save(dat.all.southern, file='../data/modelled/dat.all.southern.RData')

## Tow level data
manta.tow.southern = manta.tow %>%
    filter(Region=='Southern GBR') %>%
    droplevels %>%
    mutate(oLIVE_CORAL=factor(LIVE_CORAL,
                              levels=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U'),
                              ordered=TRUE),
           nREEF_NAME=as.numeric(as.factor(REEF_NAME)),
           nLIVE_CORAL=as.numeric(oLIVE_CORAL),
           REEF_YEAR = interaction(REEF_NAME, Year)
           )
save(manta.tow.southern, file='../data/modelled/manta.tow.southern.RData')
## ----end
## ---- Southern.Raw cells
dat.all.southern.cellmeans <- cellMeansRaw(dat.all.southern)
rawAdd <- ggproto_Raw(dat.all.southern.cellmeans)

manta.tow.southern.cellmeans <- cellMeansRaw(manta.tow %>%
                                             filter(Region=='Northern GBR') %>%
                                             droplevels %>%
                                             group_by(P_CODE.mod, REEF_NAME, Year) %>%
                                             summarise(Cover=mean(Cover), Tows=length(unique(TOW_SEQ_NO))) %>%
                                             ungroup)
## ----end
## ---- Southern.glmmTMB.tow.beta disp
if ('glmmTMB_tow beta disp' %in% models) {
    mod.southern_glmmTMB.beta.disp <- glmmTMB(Cover ~ Year + (1|REEF_NAME/REEF_YEAR),
                                     dispformula = ~Year,
                                     data=manta.tow.southern,
                                     ## weights=dat.all.southern$Tows,
                                     family=beta_family())
    dat.southern_glmmTMB.beta.disp = emmeans(mod.southern_glmmTMB.beta.disp, ~Year, type='response') %>%
        as.data.frame()
    ## DHARMa::simulateResiduals(mod.southern_glmmTMB.beta.disp, plot=TRUE)
    ## performance::check_model(mod.southern_glmmTMB.beta.disp)
    save(mod.southern_glmmTMB.beta.disp, dat.southern_glmmTMB.beta.disp, file=paste0('../data/modelled/mod.southern_glmmTMB.beta.disp.RData'))
}
## ----end
## ----end
