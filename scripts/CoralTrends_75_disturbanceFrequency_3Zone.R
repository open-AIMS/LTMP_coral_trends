## ---- prelim
library(tidyverse)
library(emmeans)
library(gridExtra)
library(grid)
library(glmmTMB)
library(scales) # for brewer_pal
source('CoralTrends_functions.R')
library(gtable)
## ----end

pre.start=1985
pre.end = 2011

YearCuts = c(1985,2001,2011,2016)
CutLabels = c('Distant Past', 'Past', 'Recent')
CutLabels = c('1986-2001', '2002-2011', '2012-2016')
## Start with table of disturbances
## ---- LoadData
load(file='../data/modelled/bleaching.full_3Zone.RData')
## load(file='../data/modelled/bleaching.merge_3Zone.RData')
load(file='../data/modelled/cots.full_3Zone.RData')
load(file='../data/modelled/cyclones.full_3Zone.RData')
load(file='../data/processed/all.reefs.cyclones.RData')
## ----end
## New analysis
## ---- helperFunctions
{
    ## ---- calcFreqs
    calcFreqs = function(mod) {
                                        #l1=lsmeans(mod, specs='time', by='Zone',type='response', at=list(time=0))
        l1=emmeans(mod, specs='time', by='Zone',type='response', at=list(time=0))
        list(Intercepts=l1)
        ##ideally, it would be good to be able to nominate the family to use here
        ##since the backtransform for slopes would just be exp rather than the
        ##inverse of a logit.  Add 1 to the slopes and intervals....
        ##SLOPES ARE ON A log odd-ratio scale.  WE SHOULD exp them so that they represent
        ## the factor change per year (or if them multiply by 100, the percent change per year)
                                        #lt=lstrends(mod, specs='Zone', var='time')
        lt=emtrends(mod, specs='Zone', var='time')
        l2=test(lt)
        list(Intercepts=as.data.frame(summary(l1)) %>% full_join(test(l1)) %>% mutate(p.value=round(p.value,3)),
             slopes=as.data.frame(summary(lt)) %>% full_join(l2) %>% mutate(p.value=round(p.value,3)))
    }
    ## ----end
    ## ---- calcFreqs.matrix
    calcFreqs.matrix = function(mod) {
                                        #dat=recover.data.glmmTMB(mod)#mod$frame
        dat=recover.data.glmmTMB(mod)#mod$frame
        form=formula(delete.response(terms(mod)))
        ## Slopes - actually rates (change in probability of being impacted per year)
        newdata = data.frame(time=1, Zone=factor(levels(dat$Zone),levels=levels(dat$Zone)))
        Xmat = model.matrix(form, data=newdata)
        Xmat[,c(1,3,4)]=0
        coefs = fixef(mod)[[1]]
        (fit=as.vector(coefs %*% t(Xmat)))
        SE=sqrt(diag(Xmat %*% vcov(mod)[[1]] %*% t(Xmat)))
        q=qnorm(0.975) #asymptotic (z test)
        l1=data.frame(fit=exp(fit), lower=exp(fit-q*SE), upper=exp(fit+q*SE))
        
        ## Intercepts - probabilty of being impacted at time 0 (1985)
        newdata = data.frame(time=0, Zone=factor(levels(dat$Zone),levels=levels(dat$Zone)))
        Xmat = model.matrix(form, data=newdata)
        coefs = fixef(mod)[[1]]
        fit=as.vector(coefs %*% t(Xmat))
        SE=sqrt(diag(Xmat %*% vcov(mod)[[1]] %*% t(Xmat)))
        q=qnorm(0.975) #asymptotic (z test)
        l2=data.frame(fit=binomial()$linkinv(fit), lower=binomial()$linkinv(fit-q*SE), upper=binomial()$linkinv(fit+q*SE))

        list(Intercept=l2, Slope=l1)
    }
    ## ----end
    ## ---- ACF.glmmTMB
    ACF.glmmTMB <- 
        function (object, maxLag, resType = c("pearson", "response", 
                                              "deviance","raw"), re=names(object$modelInfo$reTrms$cond$flist[1]),...) 
    {
        resType <- match.arg(resType)
        res <- resid(object, type = resType)
        res = split(res,object$modelInfo$reTrms$cond$flist[[re]])
        if (missing(maxLag)) {
            maxLag <- min(c(maxL <- max(lengths(res)) - 1, as.integer(10 * 
                                                                      log10(maxL + 1))))
        }
        val <- lapply(res, function(el, maxLag) {
            N <- maxLag + 1L
            tt <- double(N)
            nn <- integer(N)
            N <- min(c(N, n <- length(el)))
            nn[1:N] <- n + 1L - 1:N
            for (i in 1:N) {
                tt[i] <- sum(el[1:(n - i + 1)] * el[i:n])
            }
            array(c(tt, nn), c(length(tt), 2))
        }, maxLag = maxLag)
        val0 <- rowSums(sapply(val, function(x) x[, 2]))
        val1 <- rowSums(sapply(val, function(x) x[, 1]))/val0
        val2 <- val1/val1[1L]
        z <- data.frame(lag = 0:maxLag, ACF = val2)
        attr(z, "n.used") <- val0
        class(z) <- c("ACF", "data.frame")
        z
    }
    ## ----end
    ## ---- recover.data.glmmTMB
    recover.data.glmmTMB <- function(object, ...) {
        fcall <- getCall(object)
        recover_data(fcall,delete.response(terms(object)),
                     attr(model.frame(object),"na.action"), ...)
    }
    ## ----end
    ## ---- lsm.basis.glmmTMB
    lsm.basis.glmmTMB <- function (object, trms, xlev, grid, vcov.,
                                   mode = "asymptotic", component="cond", ...) {
        if (mode != "asymptotic") stop("only asymptotic mode is available")
        if (component != "cond") stop("only tested for conditional component")
        if (missing(vcov.)) 
            V <- as.matrix(vcov(object)[[component]])
        else V <- as.matrix(.my.vcov(object, vcov.))
        dfargs = misc = list()
        if (!is.null(object$modelInfo$family)) {
            fam = object$modelInfo$family$family
            misc$tran = object$modelInfo$family$link
            misc$inv.lbl = "response"
            if (!is.na(pmatch(fam, "binomial"))) 
                misc$inv.lbl = "prob"
            else if (!is.na(pmatch(fam, "poisson"))) 
                misc$inv.lbl = "rate"
        }
                                        #misc = lsmeans:::.std.link.labels(object$modelInfo$family, misc)
        if (mode == "asymptotic") {
            dffun = function(k, dfargs) NA
        }
        ## use this? misc = .std.link.labels(family(object), misc)
        contrasts = attr(model.matrix(object), "contrasts")
        m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
        X = model.matrix(trms, m, contrasts.arg = contrasts)
        bhat = fixef(object)[[component]]
        if (length(bhat) < ncol(X)) {
            kept = match(names(bhat), dimnames(X)[[2]])
            bhat = NA * X[1, ]
            bhat[kept] = fixef(object)[[component]]
            modmat = model.matrix(trms, model.frame(object), contrasts.arg = contrasts)
            nbasis = estimability::nonest.basis(modmat)
        }
        else nbasis = estimability::all.estble
        list(X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun, 
             dfargs = dfargs, misc = misc,...)
    }
    ## ----end
    ## ---- plotEffects
    plotEffects = function(mod, firstYear,individualReefs=NULL,
                           ribbonFillColor='blue', points = TRUE,
                           type = 1, maxYear = (finalYear)) {
        dat=recover.data.glmmTMB(mod)#mod$frame
        tt <- terms(mod)
        Terms <- delete.response(tt)

        newdata = NULL
        for (z in unique(dat$Zone)) {
            dat1 = dat %>% filter(Zone==z)
            pts = with(dat1,
                       expand.grid(Zone=z,
                                   time=seq(min(time), max(time), len=100)))
            newdata = rbind(newdata,pts)
        }
        newdata = newdata %>%
            mutate(Zone=factor(Zone,
                               levels=c('Northern GBR','Central GBR','Southern GBR')))
        m <- model.frame(Terms, newdata)
        Xmat <- model.matrix(Terms, m)
        coefs = fixef(mod)[[1]]
        fit=as.vector(coefs %*% t(Xmat))
        SE=sqrt(diag(Xmat %*% vcov(mod)[[1]] %*% t(Xmat)))
        q=qnorm(0.975) #asymptotic (z test)
        newdata = cbind(newdata, data.frame(fit=binomial()$linkinv(fit),
                                            lower=binomial()$linkinv(fit-q*SE),
                                            upper=binomial()$linkinv(fit+q*SE))) %>%
            mutate(Date=firstYear + time)

        if (!is.null(individualReefs)) {
            individualReefs = individualReefs %>%
                mutate(Date=firstYear + time, Time=factor(time))
            if (type ==1) {  
                g1=ggplot(newdata, aes(y=fit, x=Date)) +
                    geom_blank()
                
                if(points) g1 <- g1 + geom_point(data=individualReefs, aes(y=fit), color='grey')
                g1 <- g1 + 
                    geom_boxplot(data=individualReefs, aes(y=fit, group=Time), outlier.shape=NA) +
                    facet_grid(~Zone) +
                    geom_ribbon(aes(ymin=lower, ymax=upper), fill=ribbonFillColor, alpha=0.5, color=NA) +
                    geom_line(color=ribbonFillColor) +
                    scale_x_continuous('', expand=c(0,0)) + #, limits=c(1984,maxYear))+
                    coord_cartesian(xlim=c(1984, (finalYear+1))) +
                    scale_y_continuous('Pr(impact)', limits=c(0,1.00))+
                    theme_classic()+
                    theme(strip.background = element_blank(), plot.margin=unit(c(0,2,0,1),'lines'),
                          panel.spacing=unit(1,'lines'))
            } else {  ## if need to transpose
                individualReefs <- individualReefs %>% mutate(Col = 1) %>%
                    bind_rows(individualReefs %>% mutate(Col = 2)) %>%
                    bind_rows(individualReefs %>% mutate(Col = 3))
                
                g1=ggplot(newdata, aes(y=fit, x=Date)) +
                    geom_blank()
                if(points) g1 <- g1 + geom_point(data=individualReefs, aes(y=fit), color='grey')
                g1 <- g1 + 
                    geom_boxplot(data=individualReefs, aes(y=fit, group=Time), outlier.shape=NA) +
                    facet_grid(Zone~Col) +
                    geom_ribbon(aes(ymin=lower, ymax=upper), fill=ribbonFillColor, alpha=0.5, color=NA) +
                    geom_line(color=ribbonFillColor) +
                    scale_x_continuous('', expand=c(0,0)) + #, limits=c(1984,2020))+
                    coord_cartesian(xlim=c(1984, (finalYear+1))) +
                    scale_y_continuous('Pr(impact)', limits=c(0,1.00))+
                    theme_classic()+
                    theme(strip.background = element_blank(), plot.margin=unit(c(0,2,0,1),'lines'),
                          panel.spacing=unit(1,'lines'))
            }
        } else {
            g1=ggplot(newdata, aes(y=fit, x=Date)) +
                geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, color=NA) +
                geom_line() +
                facet_grid(~Zone) +
                scale_y_continuous('Pr(impact)')+
                theme_classic()
        }
        g1
    }
    ## ----end
    ## ---- modelEachReef
    modelEachReef = function(mod,form=CYCLONEany~time, dat=mod$frame) {
        df = list()
        for (i in 1:length(unique(dat$REEF_NAME))) {
            d1=dat %>% filter(REEF_NAME==unique(dat$REEF_NAME)[i]) %>% droplevels
            if (nrow(d1)>3) {
                newmod = glm(form,data=d1, family=binomial)
                d1=data.frame(REEF_NAME=unique(d1$REEF_NAME),
                              Zone=unique(d1$Zone),
                              time=seq(min(d1$time),max(d1$time),by=1))
                df[[i]]=cbind(d1, fit=predict(newmod, newdata=d1, type='response'))
            }
        }
        df = do.call('rbind',df)
    }
    ## ----end
    ## ---- modelEachreefar1
    modelEachReefAR1 = function(mod,form=CYCLONEany~time) {
        dat=recover.data.glmmTMB(mod) %>% full_join(mod$frame)
        wch=grep('poly',colnames(dat))
        dat=dat[,-wch]
        df = list()
        for (i in 1:length(unique(dat$REEF_NAME))) {
            d1=dat %>% dplyr::filter(REEF_NAME==unique(dat$REEF_NAME)[i]) %>% droplevels
            newmod = glmmTMB(form,data=d1, family=binomial)
            d1=data.frame(REEF_NAME=unique(d1$REEF_NAME), Zone=unique(d1$Zone),time=seq(min(d1$time),max(d1$time),by=1))
            df[[i]]=cbind(d1, fit=predict(newmod, newdata=d1, type='response'))
        }
        df = do.call('rbind',df)
    }
    ## ----end
}
## ----end

dist.table = list()

## ---- Cyclones
{
    ## ---- cycloneModel
cyclones.full = all.reefs  %>% left_join(cyclones.full) %>% filter(!is.na(CYCLONEcat)) %>% filter(REPORT_YEAR < (finalYear +1))
## spread the data so that for each reef/year there is a binary response for each category
cyclones.binary = cyclones.full %>%
  mutate(CYCLONEcat=as.factor(CYCLONEcat)) %>%
  bind_cols %>% data.frame(model.matrix(~-1+CYCLONEcat, data=.)) %>%
  mutate(CYCLONEany=ifelse(CYCLONEcat0==1, 0,1),
         CYCLONEsevere=ifelse((CYCLONEcat2+CYCLONEcat3)<1,0,1),
         time = REPORT_YEAR-min(REPORT_YEAR),
         Zone=factor(Zone, levels=c('Northern GBR','Central GBR','Southern GBR')))



# Start off by exploring very short-term changes with splines.
library(mgcv)
cyclones.any.gamm = gamm(CYCLONEany ~ s(REPORT_YEAR, by=Zone), random=list(REEF_NAME=~1),
                           data=cyclones.binary, family=binomial())
summary(cyclones.any.gamm$gam)
plot(cyclones.any.gamm$gam, pages=1, shift=fixef(cyclones.any.gamm$lme)[[1]], trans=binomial()$linkinv, ylim=c(-5,5))

## GAM's over fit.  We really want to be able to say whether the long-term frequency of
## various disturbances has changed over time - GAMS are too short term.


library(glmmTMB)
cyclones.any.glmmTMB = glmmTMB(CYCLONEany ~ time*Zone + (1|REEF_NAME),
                               data=cyclones.binary, family=binomial())
#cyclones.any.glmmTMB = glmmTMB(CYCLONEany ~ poly(time,3)*Zone + (1|REEF_NAME),
#                               datq=cyclones.binary, family=binomial())
cyclones.any.glmmTMB1 = glmmTMB(CYCLONEany ~ time*Zone + (Zone|REEF_NAME),
                                data=cyclones.binary, family=binomial())
anova(cyclones.any.glmmTMB,cyclones.any.glmmTMB1)
## random intercept models fine..

plot(ACF(cyclones.any.glmmTMB, resType="pearson"), alpha=0.05)
## no evidence of temporal autocorrelation)

#cyclones.any.glmmTMB2 = glmmTMB(CYCLONEany ~ time*Zone + (1|REEF_NAME) + ar1(-1+time|Zone/REEF_NAME),
#                                data=cyclones.binary, family=binomial())

save(cyclones.any.glmmTMB, file='../data/modelled/cyclones.any.glmmTMB.RData')
## ----end
    ## ---- LoadCycloneModel
load(file='../data/modelled/cyclones.any.glmmTMB.RData')
summary(cyclones.any.glmmTMB)
## ----end
    ## Note, in the following, I would prefer the slopes where 1+slope estimate
    ## The underlying lstrends function backtransforms from logit (to produce odds) rather than
    ## log (to produce odds ratio).
    ## Ratio would be more intuitive on the natural scale as we can then say that
    ## for every one unit increase in time, the probability of experiencing an event
    ## increases by a factor of ...
    ## When backtransformed to odds, it does not have an interpretation.
    ## ---- CyclonesAny
calcFreqs(cyclones.any.glmmTMB)
calcFreqs.matrix(cyclones.any.glmmTMB)
df=modelEachReef(cyclones.any.glmmTMB,form=CYCLONEany~time)
d1=plotEffects(cyclones.any.glmmTMB, firstYear=min(cyclones.binary$REPORT_YEAR), individualReefs=df)
ggsave(file='../output/figures/Disturbances_Cyclones.any.pdf',d1,width=7, height=3)
ggsave(file='../output/figures/Disturbances_Cyclones.any.png',d1,width=7, height=3, dpi=300)
## ----end
    ## ---- CyclonesSevere - This is the one used..
    dist.table[['Cyclones']] = list()
    cyclones.severe.glmmTMB = glmmTMB(CYCLONEsevere ~ time*Zone + (1|REEF_NAME),
                                      data=cyclones.binary, family=binomial())
    plot(ACF(cyclones.severe.glmmTMB, resType="pearson"), alpha=0.05)
    calcFreqs(cyclones.severe.glmmTMB)
    calcFreqs.matrix(cyclones.severe.glmmTMB)
    dist.table[['Cyclones']][['Intercept']] = calcFreqs(cyclones.severe.glmmTMB)[[1]] %>%
        dplyr::select(-time) %>%
        mutate(Disturbance='Cyclones', Stat='Intercept') %>%
        dplyr::select(Disturbance,Stat,everything())
    dist.table[['Cyclones']][['Slope']] = calcFreqs(cyclones.severe.glmmTMB)[[2]] %>%
        mutate(Disturbance='Cyclones', Stat='Slope') %>%
        dplyr::select(Disturbance,Stat,everything())                                
    dist.table[['Cyclones']][['PercentChange']] = ((calcFreqs.matrix(cyclones.severe.glmmTMB)[[2]]-1)*100) %>%
        mutate(Disturbance='Cyclones', Stat='PercentChange') %>%
        dplyr::select(Disturbance,Stat,everything())  
    dist.table[['Cyclones']][['InterceptProb']] = calcFreqs.matrix(cyclones.severe.glmmTMB)[[1]] %>%
        mutate(Disturbance='Cyclones', Stat='InterceptProb') %>%
        dplyr::select(Disturbance,Stat,everything())  

    df=modelEachReef(cyclones.severe.glmmTMB,form=CYCLONEsevere~time)
    d1=plotEffects(cyclones.severe.glmmTMB, firstYear=min(cyclones.binary$REPORT_YEAR), individualReefs=df, ribbonFillColor=brewer_pal(palette='Blues')(4)[4], points=FALSE)
    g.severe.cyclones=d1
    save(g.severe.cyclones, file='../data/modelled/g.severe.cyclones.RData')
    d1_2021=plotEffects(cyclones.severe.glmmTMB, firstYear=min(cyclones.binary$REPORT_YEAR), individualReefs=df, ribbonFillColor=brewer_pal(palette='Blues')(4)[4], points=FALSE, type=2)
    g.severe.cyclones_2021=d1_2021
    save(g.severe.cyclones_2021, file='../data/modelled/g.severe.cyclones_2021.RData')
    ggsave(file='../output/figures/Disturbances_Cyclones.severe.pdf',d1,width=7, height=3)
    ggsave(file='../output/figures/Disturbances_Cyclones.severe.png',d1,width=7, height=3, dpi=300)
    save(cyclones.severe.glmmTMB, file='../data/modelled/cyclones.severe.glmmTMB.RData')
    ## ----end

    ## ---- CyclonesCat3
ggplot(cyclones.binary, aes(y=CYCLONEcat3, x=REPORT_YEAR)) + geom_point() + facet_wrap(~Zone)
cyclones.cat3.glmmTMB = glmmTMB(CYCLONEcat3 ~ time*Zone + (1|REEF_NAME),
                               data=cyclones.binary, family=binomial())
plot(ACF(cyclones.cat3.glmmTMB, resType="pearson"), alpha=0.05)
calcFreqs(cyclones.cat3.glmmTMB)
calcFreqs.matrix(cyclones.cat3.glmmTMB)
df=modelEachReef(cyclones.cat3.glmmTMB,form=CYCLONEcat3~time)
d1=plotEffects(cyclones.cat3.glmmTMB, firstYear=min(cyclones.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Blues')(4)[4])
ggsave(file='../output/figures/Disturbances_Cyclones.cat3.pdf',d1,width=7, height=3)
ggsave(file='../output/figures/Disturbances_Cyclones.cat3.png',d1,width=7, height=3, dpi=300)
save(cyclones.cat3.glmmTMB, file='../data/modelled/cyclones.cat3.glmmTMB.RData')
## ----end

    ## ---- CyclonesCat2
ggplot(cyclones.binary, aes(y=CYCLONEcat2, x=REPORT_YEAR)) + geom_point() + facet_wrap(~Zone)
cyclones.cat2.glmmTMB = glmmTMB(CYCLONEcat2 ~ time*Zone + (1|REEF_NAME),
                               data=cyclones.binary, family=binomial())
plot(ACF(cyclones.cat2.glmmTMB, resType="pearson"), alpha=0.05)
calcFreqs(cyclones.cat2.glmmTMB)
calcFreqs.matrix(cyclones.cat2.glmmTMB)
df=modelEachReef(cyclones.cat2.glmmTMB,form=CYCLONEcat2~time)
d1=plotEffects(cyclones.cat2.glmmTMB, firstYear=min(cyclones.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Blues')(4)[4])
ggsave(file='../output/figures/Disturbances_Cyclones.cat2.pdf',d1,width=7, height=3)
ggsave(file='../output/figures/Disturbances_Cyclones.cat2.png',d1,width=7, height=3, dpi=300)
save(cyclones.cat2.glmmTMB, file='../data/modelled/cyclones.cat2.glmmTMB.RData')
## ----end

    ##polynomials
    ## ---- CyclonesAnyPoly
cyclones.any.glmmTMB3 = glmmTMB(CYCLONEany ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME),
                                data=cyclones.binary, family=binomial())
save(cyclones.any.glmmTMB3, file='../data/modelled/cyclones.any.glmmTMB3.RData')

plot(ACF(cyclones.any.glmmTMB3, resType="pearson"), alpha=0.05)
summary(cyclones.any.glmmTMB3)
df=modelEachReef(cyclones.any.glmmTMB3,form=CYCLONEany~poly(time,3),dat=cyclones.binary)
d1=plotEffects(cyclones.any.glmmTMB3, firstYear=min(cyclones.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Blues')(4)[4])
ggsave(file='../output/figures/Disturbances_Cyclones.any_polygons.pdf',d1,width=7, height=3)
ggsave(file='../output/figures/Disturbances_Cyclones.any_polygons.png',d1,width=7, height=3, dpi=300)
## ----end

    ## ---- CyclonesSeverePoly
cyclones.severe.glmmTMB3 = glmmTMB(CYCLONEsevere ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME),
                                   data=cyclones.binary, family=binomial())
save(cyclones.severe.glmmTMB3, file='../data/modelled/cyclones.severe.glmmTMB3.RData')
plot(ACF(cyclones.severe.glmmTMB3, resType="pearson"), alpha=0.05)
summary(cyclones.severe.glmmTMB3)
df=modelEachReef(cyclones.severe.glmmTMB3,form=CYCLONEsevere~poly(time,3),dat=cyclones.binary)
d1=plotEffects(cyclones.severe.glmmTMB3, firstYear=min(cyclones.binary$REPORT_YEAR), individualReefs=df)
ggsave(file='../output/figures/Disturbances_Cyclones.severe_polygons.pdf',d1,width=7, height=3)
ggsave(file='../output/figures/Disturbances_Cyclones.severe_polygons.png',d1,width=7, height=3, dpi=300)
## ----end
}
## ----end

## ---- COTS
{
    cots.full = all.reefs  %>% left_join(cots.full %>% distinct) %>% filter(!is.na(COTScat)) %>%
        filter(REPORT_YEAR < (finalYear + 1))
    ## spread the data so that for each reef/year there is a binary response for each category
    ##For some reason, some reefs (e.g. '16017S') have two different COTScat in a particular year (1994)
    ## To correct for this, I will give them the max category.
    ## ---- COTSBinary
cots.binary = cots.full %>%
    mutate(COTScat=factor(COTScat, levels=c('Zero','NO','IO','AO'))) %>%
    group_by(REEF_NAME,REEF_ID,Zone,Latitude,Longitude,REPORT_YEAR) %>%
    summarize(COTScat=levels(COTScat)[max(as.numeric(COTScat))]) %>%
    ungroup %>%
    mutate(COTScat=factor(COTScat, levels=c('Zero','NO','IO','AO'))) %>%
    bind_cols %>% data.frame(model.matrix(~-1+COTScat, data=.)) %>%
    mutate(COTSany=ifelse(COTScatZero==1, 0,1),
           COTSsevere=ifelse((COTScatIO+COTScatAO)<1,0,1),
           time = REPORT_YEAR-min(REPORT_YEAR),
           Time=as.factor(time),
           Zone=factor(Zone, levels=c('Northern GBR','Central GBR','Southern GBR')))
save(cots.binary, file='../data/modelled/cots.binary.RData')

g1 = ggplot(cots.binary, aes(y=COTSany, x=REPORT_YEAR)) +
    geom_count(aes(size=..prop.., group=REPORT_YEAR)) + scale_size_area() +
    facet_wrap(~Zone)
ggsave(file='../output/figures/cots_any.pdf', g1,width=10, height=3)
## ----end
                                        #cots.any.glmmTMB = glmmTMB(COTSany ~ time*Zone + (1|REEF_NAME) + ar1(as.factor(time)-1|REEF_NAME),
                                        #                           data=cots.binary, family=binomial())
    ## ---- CotsAny
cots.any.glmmTMB = glmmTMB(COTSany ~ time*Zone + (1|REEF_NAME)+ ar1(-1+factor(time)|REEF_NAME),
                           data=cots.binary, family=binomial())
summary(cots.any.glmmTMB)
plot(ACF(cots.any.glmmTMB, resType="pearson"), alpha=0.05)
calcFreqs(cots.any.glmmTMB)
calcFreqs.matrix(cots.any.glmmTMB)

#df=modelEachReefAR1(cots.any.glmmTMB,form=COTSany~time+ar1(-1+factor(time)|REEF_NAME))
df=modelEachReef(cots.any.glmmTMB,form=COTSany~time,dat=cots.binary)
d1=plotEffects(cots.any.glmmTMB, firstYear=min(cots.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Greens')(3)[3])
ggsave(file='../output/figures/Disturbances_COTS_any.pdf',d1,width=7, height=3)
ggsave(file='../output/figures/Disturbances_COTS_any.png',d1,width=7, height=3, dpi=300)
#ggsave(file='figures/cots.any.glmmTMB.pdf', d1,width=7, height=3)
save(cots.any.glmmTMB, file='../data/modelled/cots.any.glmmTMB.RData')
## ----end

    ##polynomials
    ## ---- CotsAnyPoly
cots.any.glmmTMB3 = glmmTMB(COTSany ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME) + ar1(-1+factor(time)|REEF_NAME),
                                data=cots.binary, family=binomial())
save(cots.any.glmmTMB3, file='../data/modelled/cots.any.glmmTMB3.RData')

plot(ACF(cots.any.glmmTMB3, resType="pearson"), alpha=0.05)
summary(cots.any.glmmTMB3)
#df=modelEachReefAR1(cots.any.glmmTMB3,form=COTSany~poly(time,3) + ar1(-1+factor(time)|REEF_NAME))
#save(df, file='data/modelled/cots.any.glmmTMB3.df.RData')
df=modelEachReef(cots.any.glmmTMB3,form=COTSany~poly(time,3),dat=cots.binary)
d1=plotEffects(cots.any.glmmTMB3, firstYear=min(cots.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Greens')(3)[3])
                                        #ggsave(file='figures/cots.any.glmmTMB3.pdf', g1,width=10, height=3)
ggsave(file='../output/figures/Disturbances_COTS_any_poly.pdf',d1,width=7, height=3)
ggsave(file='../output/figures/Disturbances_COTS_any_poly.png',d1,width=7, height=3, dpi=300)
## ----end
    ## ---- CotsSevere
g1=ggplot(cots.binary, aes(y=COTSsevere, x=REPORT_YEAR)) +
    geom_count(aes(size=..prop.., group=REPORT_YEAR)) + scale_size_area() +
    facet_wrap(~Zone)
ggsave(file='../output/figures/cots_severe.pdf', g1,width=10, height=3)

#cots.severe.glmmTMB = glmmTMB(COTSsevere ~ time*Zone + (1|REEF_NAME) + ar1(-1+factor(time)|REEF_NAME),
#                              data=cots.binary, family=binomial())
cots.severe.glmmTMB = glmmTMB(COTSsevere ~ time*Zone + (1|REEF_NAME),
                               data=cots.binary, family=binomial())
plot(ACF(cots.severe.glmmTMB, resType="pearson"), alpha=0.05)
calcFreqs(cots.severe.glmmTMB)
calcFreqs.matrix(cots.severe.glmmTMB)
## ----
#df=modelEachReefAR1(cots.severe.glmmTMB,form=COTSsevere~time+ar1(-1+factor(time)|REEF_NAME))
df=modelEachReef(cots.severe.glmmTMB,form=COTSsevere~time, dat=cots.binary)
d1=plotEffects(cots.severe.glmmTMB, firstYear=min(cots.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Greens')(3)[3], points=FALSE)
g.severe.cots = d1
d1_2021=plotEffects(cots.severe.glmmTMB, firstYear=min(cots.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Greens')(3)[3], points=FALSE, type=2)
g.severe.cots_2021 = d1_2021
ggsave(file='../output/figures/Disturbances_COTS_severe.pdf',d1,width=7, height=3)
ggsave(file='../output/figures/Disturbances_COTS_severe.png',d1,width=7, height=3, dpi=300)
#ggsave(file='figures/cots.severe.glmmTMB.pdf', g1,width=10, height=3)
save(cots.severe.glmmTMB, file='../data/modelled/cots.severe.glmmTMB.RData')
## ----end

    ##polynomials
    ## ---- CotsSeverePoly
#cots.severe.glmmTMB3 = glmmTMB(COTSsevere ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME) + ar1(-1+factor(time)|REEF_NAME),
#                               data=cots.binary, family=binomial())
cots.severe.glmmTMB3 = glmmTMB(COTSsevere ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME),
                                data=cots.binary, family=binomial())
save(cots.severe.glmmTMB3, file='../data/modelled/cots.severe.glmmTMB3.RData')

plot(ACF(cots.severe.glmmTMB3, resType="pearson"), alpha=0.05)
summary(cots.severe.glmmTMB3)
#df=modelEachReefAR1(cots.severe.glmmTMB3,form=COTSsevere~poly(time,3) + ar1(-1+factor(time)|REEF_NAME))
#save(df, file='data/modelled/cots.severe.glmmTMB3.df.RData')
df=modelEachReef(cots.severe.glmmTMB3,form=COTSsevere~poly(time,3),dat=cots.binary)
d1=plotEffects(cots.severe.glmmTMB3, firstYear=min(cots.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Greens')(3)[3])
ggsave(file='../output/figures/Disturbances_COTS_severe_poly.pdf',d1,width=7, height=3)
ggsave(file='../output/figures/Disturbances_COTS_severe_poly.png',d1,width=7, height=3, dpi=300)
#ggsave(file='figures/cots.severe.glmmTMB3.pdf', g1,width=10, height=3)
## ----end

    ## ---- CotsAO - This is the one used..
    dist.table[['COTS']] = list()
    g1=ggplot(cots.binary, aes(y=COTScatAO, x=REPORT_YEAR)) +
        geom_count(aes(size=..prop.., group=REPORT_YEAR)) + scale_size_area() +
        facet_wrap(~Zone)
    ggsave(file='../figures/cots_catAO.pdf', g1,width=10, height=3)

    cots.catAO.glmmTMB = glmmTMB(COTScatAO ~ time*Zone + (1|REEF_NAME),# + ar1(-1+time|REEF_NAME),
                                 data=cots.binary, family=binomial())
    save(cots.catAO.glmmTMB, file='../data/modelled/cots.catAO.glmmTMB.RData')
    plot(ACF(cots.catAO.glmmTMB, resType="pearson"), alpha=0.05)

    calcFreqs(cots.catAO.glmmTMB)
    calcFreqs.matrix(cots.catAO.glmmTMB)
    dist.table[['COTS']][['Intercept']] = calcFreqs(cots.catAO.glmmTMB)[[1]] %>%
        dplyr::select(-time) %>%
        mutate(Disturbance='COTS', Stat='Intercept') %>%
        dplyr::select(Disturbance,Stat,everything())
    dist.table[['COTS']][['Slope']] = calcFreqs(cots.catAO.glmmTMB)[[2]] %>%
        mutate(Disturbance='COTS', Stat='Slope') %>%
        dplyr::select(Disturbance,Stat,everything())
    dist.table[['COTS']][['PercentChange']] = ((calcFreqs.matrix(cots.severe.glmmTMB)[[2]]-1)*100) %>%
        mutate(Disturbance='COTS', Stat='PercentChange') %>%
        dplyr::select(Disturbance,Stat,everything())  
    dist.table[['COTS']][['InterceptProb']] = calcFreqs.matrix(cots.severe.glmmTMB)[[1]] %>%
        mutate(Disturbance='COTS', Stat='InterceptProb') %>%
        dplyr::select(Disturbance,Stat,everything())  

                                        #df=modelEachReefAR1(cots.catAO.glmmTMB,form=COTScatAO~time+ar1(-1+time|REEF_NAME))
    df=modelEachReef(cots.catAO.glmmTMB,form=COTScatAO~time)
    d1=plotEffects(cots.catAO.glmmTMB, firstYear=min(cots.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Greens')(3)[3], points=FALSE)
    g.AO.cots = d1
    d1_2021=plotEffects(cots.catAO.glmmTMB, firstYear=min(cots.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Greens')(3)[3], points=FALSE, type=2)
    g.AO.cots_2021 = d1_2021
    save(g.AO.cots, file='../data/modelled/g.AO.cots.RData')
    save(g.AO.cots_2021, file='../data/modelled/g.AO.cots_2021.RData')
    ggsave(file='../output/figures/Disturbances_COTS_AO.pdf',d1,width=7, height=3)
    ggsave(file='../output/figures/Disturbances_COTS_AO.png',d1,width=7, height=3, dpi=300)
                                        #ggsave(file='figures/cots.catAO.glmmTMB.pdf', g1,width=10, height=3)
    ## ----end

    ##polynomials
    ## ---- CotsCatAOPoly
cots.catAO.glmmTMB3 = glmmTMB(COTScatAO ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME),
                                data=cots.binary, family=binomial())
save(cots.catAO.glmmTMB3, file='../data/modelled/cots.catAO.glmmTMB3.RData')

plot(ACF(cots.catAO.glmmTMB3, resType="pearson"), alpha=0.05)
summary(cots.catAO.glmmTMB3)
#df=modelEachReefAR1(cots.catAO.glmmTMB3,form=COTScatAO~poly(time,3) + ar1(-1+factor(time)|REEF_NAME))
#save(df, file='data/modelled/cots.catAO.glmmTMB3.df.RData')
df=modelEachReef(cots.catAO.glmmTMB3,form=COTScatAO~poly(time,3),dat=cots.binary)
d1=plotEffects(cots.catAO.glmmTMB3, firstYear=min(cots.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Greens')(3)[3])
ggsave(file='../output/figures/Disturbances_COTS_AO_poly.pdf',d1,width=7, height=3)
ggsave(file='../output/figures/Disturbances_COTS_AO_poly.png',d1,width=7, height=3, dpi=300)
#ggsave(file='figures/cots.catAO.glmmTMB3.pdf', g1,width=10, height=3)
## ----end

    ## ---- COTSJUNK
## ggplot(cots.binary, aes(y=COTScatIO, x=REPORT_YEAR)) + geom_point() + facet_wrap(~Zone)
## cots.catIO.glmmTMB = glmmTMB(COTScatIO ~ time*Zone + (1|REEF_NAME) + ar1(-1+time|REEF_NAME),
##                              data=cots.binary, family=binomial())
## cots.catIO.glmmPQL = glmmPQL(COTScatIO ~ time*Zone, random=~1|REEF_NAME,
##                              data=cots.binary, family=binomial(), correlation=corAR1(form=~time|REEF_NAME))
## cots.catIO.glmmPQL = glmmPQL(COTScatIO ~ poly(time,3)*Zone, random=~1|REEF_NAME,
##                              data=cots.binary, family=binomial(), correlation=corAR1(form=~time|REEF_NAME))
## plot(allEffects(cots.catIO.glmmPQL))
## plot(ACF(cots.catIO.glmmTMB, resType="pearson"), alpha=0.05)
## ## ---- CotsIO
## calcFreqs(cots.catIO.glmmTMB)
## calcFreqs.matrix(cots.catIO.glmmTMB)
## df=modelEachReefAR1(cots.catIO.glmmTMB,form=COTScatIO~poly(time,3)+ar1(-1+time|REEF_NAME))
## plotEffects(cots.catIO.glmmTMB, firstYear=min(cots.binary$REPORT_YEAR), individualReefs=df)
## ## ----
## save(cots.catIO.glmmTMB, file='data/modelled/cots.catIO.glmmTMB.RData')


## ##polynomials
## cyclones.any.glmmTMB3 = glmmTMB(CYCLONEany ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME),
##                                 data=cyclones.binary, family=binomial())
## save(cyclones.any.glmmTMB3, file='data/modelled/cyclones.any.glmmTMB3.RData')

## plot(ACF(cyclones.any.glmmTMB3, resType="pearson"), alpha=0.05)
## summary(cyclones.any.glmmTMB3)
## ## ---- CyclonesAnyPoly
## df=modelEachReef(cyclones.any.glmmTMB3,form=CYCLONEany~poly(time,3),dat=cyclones.binary)
## plotEffects(cyclones.any.glmmTMB3, firstYear=min(cyclones.binary$REPORT_YEAR), individualReefs=df)
## ## ----

## cyclones.severe.glmmTMB3 = glmmTMB(CYCLONEsevere ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME),
##                                    data=cyclones.binary, family=binomial())
## save(cyclones.severe.glmmTMB3, file='data/modelled/cyclones.severe.glmmTMB3.RData')
## plot(ACF(cyclones.severe.glmmTMB3, resType="pearson"), alpha=0.05)
## summary(cyclones.severe.glmmTMB3)
## ## ---- CyclonesSeverePoly
## df=modelEachReef(cyclones.severe.glmmTMB3,form=CYCLONEsevere~poly(time,3),dat=cyclones.binary)
## plotEffects(cyclones.severe.glmmTMB3, firstYear=min(cyclones.binary$REPORT_YEAR), individualReefs=df)
## ## ----

## ----end
}
## ----end

## ---- bleaching
{
    ## mike supplied some bleaching data to fill in the gaps between aerial surveys
    ## load(file='data/modelled/bleaching.merge_3Zone.RData')
    load(file='../data/modelled/bleaching.full_3Zone.RData')
                                        #bleaching.full = all.reefs  %>% left_join(bleaching.merge %>% distinct) %>% filter(!is.na(BLEACHINGcat))
    bleaching.full = bleaching.full_3Zone %>% filter(!is.na(Zone)) %>% droplevels %>%
        filter(REPORT_YEAR < (finalYear + 1)) %>%
        group_by(REEF_ID) %>% arrange(REPORT_YEAR) %>% ungroup %>%
        mutate(Zone=factor(Zone, levels=c('Northern','Central','Southern'),
                           labels=c('Northern GBR','Central GBR','Southern GBR')))

    bleaching.binary = bleaching.full %>% filter(!is.na(BLEACHINGcat)) %>%
        mutate(BLEACHINGcat=as.factor(BLEACHINGcat)) %>%
        bind_cols %>% data.frame(model.matrix(~-1+BLEACHINGcat, data=.)) %>%
        mutate(BLEACHINGany=ifelse(BLEACHINGcat0==1, 0,1),
               BLEACHINGsevere=ifelse((BLEACHINGcat3+BLEACHINGcat4+BLEACHINGcat5)<1,0,1),
               ## BLEACHINGsevere=ifelse((BLEACHINGcat4+BLEACHINGcat5)<1,0,1),
               time = REPORT_YEAR-min(REPORT_YEAR),
               Zone=factor(Zone, levels=c('Northern GBR','Central GBR','Southern GBR')))
    save(bleaching.binary, file='../data/modelled/bleaching.binary.RData')

    ## now a version that exludes LTMP data outside of December - March
    bleaching.binary.2 = bleaching.full %>%
        filter((Project == 'LTMP' & Season == 'Bleaching') | Project == 'COE') %>%
        filter(!is.na(BLEACHINGcat)) %>% droplevels() %>% 
        mutate(BLEACHINGcat=as.factor(BLEACHINGcat)) %>%
        bind_cols %>% data.frame(model.matrix(~-1+BLEACHINGcat, data=.)) %>%
        mutate(BLEACHINGany=ifelse(BLEACHINGcat0==1, 0,1),
               BLEACHINGsevere=ifelse((BLEACHINGcat3+BLEACHINGcat4+BLEACHINGcat5)<1,0,1),
               ## BLEACHINGsevere=ifelse((BLEACHINGcat4+BLEACHINGcat5)<1,0,1),
               time = REPORT_YEAR-min(REPORT_YEAR),
               Zone=factor(Zone, levels=c('Northern GBR','Central GBR','Southern GBR')))
    save(bleaching.binary.2, file='../data/modelled/bleaching.binary.2.RData')

    ## ---- BleachingAny
    g1=ggplot(bleaching.binary, aes(y=BLEACHINGany, x=REPORT_YEAR)) +
        geom_count(aes(size=..prop.., group=REPORT_YEAR)) + scale_size_area() +
        facet_wrap(~Zone)
    ggsave(file='../output/figures/bleaching_any.pdf', g1,width=10, height=3)

    bleaching.any.glmmTMB = glmmTMB(BLEACHINGany ~ time*Zone + (1|REEF_NAME),
                                    data=bleaching.binary, family=binomial())
    save(bleaching.any.glmmTMB, file='../data/modelled/bleaching.any.glmmTMB.RData')

    summary(bleaching.any.glmmTMB)
    plot(ACF(bleaching.any.glmmTMB, resType="pearson"), alpha=0.05)
    calcFreqs(bleaching.any.glmmTMB)
    calcFreqs.matrix(bleaching.any.glmmTMB)
    ## ----
    df=modelEachReef(mod = bleaching.any.glmmTMB,form=BLEACHINGany~time) %>%
        mutate(Zone=factor(Zone, levels=c('Northern GBR','Central GBR','Southern GBR')))
    d1=plotEffects(bleaching.any.glmmTMB, firstYear=min(bleaching.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Reds')(5)[5])
    ggsave(file='../output/figures/Disturbances_Bleaching_any.pdf',d1,width=7, height=3)
    ggsave(file='../output/figures/Disturbances_Bleaching_any.png',d1,width=7, height=3, dpi=300)
                                        #ggsave(file='figures/bleaching.any.glmmTMB.pdf', g1,width=10, height=3)
    ## ----end


    ##polynomials
    ## ---- BleachingAnyPoly
    bleaching.any.glmmTMB3 = glmmTMB(BLEACHINGany ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME),
                                     data=bleaching.binary, family=binomial())
    save(bleaching.any.glmmTMB3, file='../data/modelled/bleaching.any.glmmTMB3.RData')

    plot(ACF(bleaching.any.glmmTMB3, resType="pearson"), alpha=0.05)
    summary(bleaching.any.glmmTMB3)
                                        #df=modelEachReefAR1(bleaching.any.glmmTMB3,form=BLEACHINGany~poly(time,3) + ar1(-1+factor(time)|REEF_NAME))
                                        #save(df, file='data/modelled/bleaching.any.glmmTMB3.df.RData')
    df=modelEachReef(bleaching.any.glmmTMB3,form=BLEACHINGany~poly(time,3),dat=bleaching.binary)
    d1=plotEffects(bleaching.any.glmmTMB3, firstYear=min(bleaching.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Reds')(5)[5])
    ggsave(file='../output/figures/Disturbances_Bleaching_poly.pdf',d1,width=7, height=3)
    ggsave(file='../output/figures/Disturbances_Bleaching_poly.png',d1,width=7, height=3, dpi=300)
                                        #ggsave(file='figures/bleaching.any.glmmTMB3.pdf', g1,width=10, height=3)
    ## ----end

    ## ---- Bleaching severe - This is the one used..
    g1=ggplot(bleaching.binary, aes(y=BLEACHINGsevere, x=REPORT_YEAR)) +
        geom_count(aes(size=..prop.., group=REPORT_YEAR)) + scale_size_area() +
        facet_wrap(~Zone)
    ggsave(file='../output/figures/bleaching_severe.pdf', g1,width=10, height=3)

    dist.table[['Bleaching']] = list()
    bleaching.severe.glmmTMB = glmmTMB(BLEACHINGsevere ~ time*Zone + (1|REEF_NAME),
                                       data=bleaching.binary, family=binomial())
    save(bleaching.severe.glmmTMB, file='../data/modelled/bleaching.severe.glmmTMB.RData')

    summary(bleaching.severe.glmmTMB)
    plot(ACF(bleaching.severe.glmmTMB, resType="pearson"), alpha=0.05)

    ## emmeans(bleaching.severe.glmmTMB, ~ time, by = 'Zone', at=list(time=c(0,1,2,3)))
    ## emmeans(bleaching.severe.glmmTMB, ~ time, by = 'Zone', at=list(time=c(0,1,2,3)), type='response')
    ## emmeans(bleaching.severe.glmmTMB, ~ time, by = 'Zone', at=list(time=c(0,1,2,3))) %>% pairs() %>% regrid()
    ## (0.008-0.00695)/(0.00695)
    ## (0.00920 - 0.00800)/(0.00800)
    ## (0.01058 - 0.00920)/(0.00920)
    
    calcFreqs(bleaching.severe.glmmTMB)
    calcFreqs.matrix(bleaching.severe.glmmTMB)
    dist.table[['Bleaching']][['Intercept']] = calcFreqs(bleaching.severe.glmmTMB)[[1]] %>%
        dplyr::select(-time) %>%
        mutate(Disturbance='Bleaching', Stat='Intercept') %>%
        dplyr::select(Disturbance,Stat,everything())
    dist.table[['Bleaching']][['Slope']] = calcFreqs(bleaching.severe.glmmTMB)[[2]] %>%
        mutate(Disturbance='Bleaching', Stat='Slope') %>%
        dplyr::select(Disturbance,Stat,everything())  
    dist.table[['Bleaching']][['PercentChange']] = ((calcFreqs.matrix(bleaching.severe.glmmTMB)[[2]]-1)*100) %>%
        mutate(Disturbance='Bleaching', Stat='PercentChange') %>%
        dplyr::select(Disturbance,Stat,everything())  
    dist.table[['Bleaching']][['InterceptProb']] = calcFreqs.matrix(bleaching.severe.glmmTMB)[[1]] %>%
        mutate(Disturbance='Bleaching', Stat='InterceptProb') %>%
        dplyr::select(Disturbance,Stat,everything())  
    ## ----
                                        #df=modelEachReefAR1(bleaching.severe.glmmTMB,form=BLEACHINGsevere~time+ar1(-1+time|REEF_NAME))
    df=modelEachReef(bleaching.severe.glmmTMB,form=BLEACHINGsevere~time)
                                        #df=modelEachReef(bleaching.severe.glmmTMB,form=BLEACHINGsevere~time)
    d1=plotEffects(bleaching.severe.glmmTMB, firstYear=min(bleaching.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Reds')(5)[5], points=FALSE)
    g.severe.bleaching = d1
    save(g.severe.bleaching, file='../data/modelled/g.severe.bleaching.RData')
    d1_2021=plotEffects(bleaching.severe.glmmTMB, firstYear=min(bleaching.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Reds')(5)[5], points=FALSE, type=2)
    g.severe.bleaching_2021 = d1_2021
    save(g.severe.bleaching_2021, file='../data/modelled/g.severe.bleaching_2021.RData')
    ggsave(file='../output/figures/Disturbances_Bleaching_severe.pdf',d1,width=7, height=3)
    ggsave(file='../output/figures/Disturbances_Bleaching_severe.png',d1,width=7, height=3, dpi=300)
                                        #ggsave(file='figures/bleaching.severe.glmmTMB.pdf', g1,width=10, height=3)
    ## ----

    ##polynomials - does not converge..
                                        #bleaching.severe.glmmTMB3 = glmmTMB(BLEACHINGsevere ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME) + ar1(-1+factor(time)|REEF_NAME),
                                        #                               data=bleaching.binary, family=binomial())
    bleaching.severe.glmmTMB3 = glmmTMB(BLEACHINGsevere ~ poly(time,2,raw=FALSE)*Zone + (1|REEF_NAME),
                                        data=bleaching.binary, family=binomial())
    save(bleaching.severe.glmmTMB3, file='../data/modelled/bleaching.severe.glmmTMB3.RData')

    plot(ACF(bleaching.severe.glmmTMB3, resType="pearson"), alpha=0.05)
    summary(bleaching.severe.glmmTMB3)
    ## ----end
    ## ---- Bleaching severe version 2 (exclude LTMP not in summer) - This is the one used..
    g1=ggplot(bleaching.binary.2, aes(y=BLEACHINGsevere, x=REPORT_YEAR)) +
        geom_count(aes(size=..prop.., group=REPORT_YEAR)) + scale_size_area() +
        facet_wrap(~Zone)
    ggsave(file='../output/figures/bleaching_severe.2.pdf', g1,width=10, height=3)

    dist.table2 = dist.table
    dist.table2[['Bleaching']] = list()
    bleaching.severe.2.glmmTMB = glmmTMB(BLEACHINGsevere ~ time*Zone + (1|REEF_NAME),
                                         data=bleaching.binary.2, family=binomial())
    save(bleaching.severe.2.glmmTMB, file='../data/modelled/bleaching.severe.2.glmmTMB.RData')

    summary(bleaching.severe.2.glmmTMB)
    plot(ACF(bleaching.severe.2.glmmTMB, resType="pearson"), alpha=0.05)

    calcFreqs(bleaching.severe.2.glmmTMB)
    calcFreqs.matrix(bleaching.severe.2.glmmTMB)
    dist.table2[['Bleaching']][['Intercept']] = calcFreqs(bleaching.severe.2.glmmTMB)[[1]] %>% dplyr::select(-time) %>% mutate(Disturbance='Bleaching', Stat='Intercept') %>% dplyr::select(Disturbance,Stat,everything())
    dist.table2[['Bleaching']][['Slope']] = calcFreqs(bleaching.severe.2.glmmTMB)[[2]] %>% mutate(Disturbance='Bleaching', Stat='Slope') %>% dplyr::select(Disturbance,Stat,everything())  
    ## ----
                                        #df=modelEachReefAR1(bleaching.severe.2.glmmTMB,form=BLEACHINGsevere~time+ar1(-1+time|REEF_NAME))
    df=modelEachReef(bleaching.severe.2.glmmTMB,form=BLEACHINGsevere~time)
                                        #df=modelEachReef(bleaching.severe.2.glmmTMB,form=BLEACHINGsevere~time)
    d1=plotEffects(bleaching.severe.2.glmmTMB, firstYear=min(bleaching.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Reds')(5)[5], points=FALSE)
    g.severe.bleaching.2 = d1
    save(g.severe.bleaching.2, file='../data/modelled/g.severe.bleaching.2.RData')
    d1_2021=plotEffects(bleaching.severe.2.glmmTMB, firstYear=min(bleaching.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Reds')(5)[5], points=FALSE, type=2)
    g.severe.bleaching_2021.2 = d1_2021
    save(g.severe.bleaching_2021, file='../data/modelled/g.severe.bleaching_2021.2.RData')
    ggsave(file='../output/figures/Disturbances_Bleaching_severe.2.pdf',d1,width=7, height=3)
    ggsave(file='../output/figures/Disturbances_Bleaching_severe.2.png',d1,width=7, height=3, dpi=300)
                                        #ggsave(file='figures/bleaching.severe.2.glmmTMB.pdf', g1,width=10, height=3)
    ## ----

    ##polynomials - does not converge..
                                        #bleaching.severe.2.glmmTMB3 = glmmTMB(BLEACHINGsevere ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME) + ar1(-1+factor(time)|REEF_NAME),
                                        #                               data=bleaching.binary, family=binomial())
    bleaching.severe.2.glmmTMB3 = glmmTMB(BLEACHINGsevere ~ poly(time,2,raw=FALSE)*Zone + (1|REEF_NAME),
                                          data=bleaching.binary.2, family=binomial())
    save(bleaching.severe.2.glmmTMB3, file='../data/modelled/bleaching.severe.2.glmmTMB3.RData')

    plot(ACF(bleaching.severe.2.glmmTMB3, resType="pearson"), alpha=0.05)
    summary(bleaching.severe.2.glmmTMB3)
    ## ----end

    ## ---- BleachingSeverePoly
                                        #df=modelEachReefAR1(bleaching.severe.glmmTMB3,form=BLEACHINGsevere~poly(time,3) + ar1(-1+factor(time)|REEF_NAME))
                                        #save(df, file='data/modelled/bleaching.severe.glmmTMB3.df.RData')
    df=modelEachReef(bleaching.severe.glmmTMB3,form=BLEACHINGsevere~poly(time,3),dat=bleaching.binary)
    d1=plotEffects(bleaching.severe.glmmTMB3, firstYear=min(bleaching.binary$REPORT_YEAR), individualReefs=df,brewer_pal(palette='Reds')(5)[5])
    ggsave(file='../output/figures/Disturbances_Bleaching_severe_poly.pdf',d1,width=7, height=3)
    ggsave(file='../output/figures/Disturbances_Bleaching_severe_poly.png',d1,width=7, height=3, dpi=300)
                                        #ggsave(file='figures/bleaching.severe.glmmTMB3.pdf', g1,width=10, height=3)
    ## ----end
}

## ----end

## ---- Collect together
{
    ## ---- Save the disturbance tables
    save(dist.table, file='../data/modelled/dist.table.RData')

    write.csv(do.call('rbind',lapply(dist.table, `[[`, 'Intercept')), file='../data/modelled/dist.table.intercepts.csv', quote=FALSE, row.names=FALSE)
    write.csv(do.call('rbind',lapply(dist.table, `[[`, 'Slope')), file='../data/modelled/dist.table.slopes.csv', quote=FALSE, row.names=FALSE)
    write.csv(do.call('rbind',lapply(dist.table, `[[`, 'PercentChange')), file='../data/modelled/dist.table.percentchange.csv', quote=FALSE, row.names=FALSE)
    write.csv(do.call('rbind',lapply(dist.table, `[[`, 'InterceptProb')), file='../data/modelled/dist.table.interceptProb.csv', quote=FALSE, row.names=FALSE)


    ## Combine all disturbances
    disturb.binary = cyclones.binary %>%  dplyr::select(REEF_ID,REPORT_YEAR,Zone,CYCLONEany,CYCLONEsevere) %>%
        full_join(cots.binary %>%  dplyr::select(REEF_ID,REPORT_YEAR,Zone,COTSany, COTSsevere)) %>%
        full_join(bleaching.binary %>%  dplyr::select(REEF_ID,REPORT_YEAR,Zone,BLEACHINGany, BLEACHINGsevere) %>% mutate(REEF_ID=as.numeric(as.character(REEF_ID)))) %>%
        distinct

    disturb.binary <- disturb.binary %>%
        mutate(across(c(-Zone, -REPORT_YEAR, -REEF_ID), function(x) ifelse(is.na(x), 0, x)))
    disturb.binary = disturb.binary %>% filter(!is.na(CYCLONEany), !is.na(BLEACHINGany), !is.na(COTSany)) %>%
        mutate(DISTURBany=ifelse((CYCLONEany+COTSany+BLEACHINGany)>0,1,0),
               DISTURBsevere = ifelse((CYCLONEsevere+COTSsevere+BLEACHINGsevere)>0,1,0),
               time = REPORT_YEAR-min(REPORT_YEAR),
               Zone=factor(Zone, levels=c('Northern GBR','Central GBR','Southern GBR')),
               REEF_NAME=REEF_ID
               )
    save(disturb.binary, file='../data/modelled/disturb.binary.RData')
    g1=ggplot(disturb.binary, aes(y=DISTURBany, x=REPORT_YEAR)) +
        geom_count(aes(size=..prop.., group=REPORT_YEAR)) + scale_size_area() +
        facet_wrap(~Zone)
    ggsave(file='../output/figures/disturb_any.pdf', g1,width=10, height=3)

    disturb.any.glmmTMB = glmmTMB(DISTURBany ~ time*Zone + (1|REEF_NAME),
                                  data=disturb.binary, family=binomial())
    save(disturb.any.glmmTMB, file='../data/modelled/disturb.any.glmmTMB.RData')
    ## ----end
    ## ---- Save the disturbance tables version 2
    save(dist.table2, file='../data/modelled/dist.table2.RData')

    write.csv(do.call('rbind',lapply(dist.table2, `[[`, 'Intercept')), file='../data/modelled/dist.table.2.intercepts.csv', quote=FALSE, row.names=FALSE)
    write.csv(do.call('rbind',lapply(dist.table2, `[[`, 'Slope')), file='../data/modelled/dist.table.2.slopes.csv', quote=FALSE, row.names=FALSE)

    ## Combine all disturbances
    disturb.binary.2 = cyclones.binary %>%  dplyr::select(REEF_ID,REPORT_YEAR,Zone,CYCLONEany,CYCLONEsevere) %>%
        full_join(cots.binary %>%  dplyr::select(REEF_ID,REPORT_YEAR,Zone,COTSany, COTSsevere)) %>%
        full_join(bleaching.binary.2 %>%  dplyr::select(REEF_ID,REPORT_YEAR,Zone,BLEACHINGany, BLEACHINGsevere) %>% mutate(REEF_ID=as.numeric(as.character(REEF_ID)))) %>%
        distinct

    disturb.binary.2 <- disturb.binary.2 %>%
        mutate(across(c(-Zone, -REPORT_YEAR, -REEF_ID), function(x) ifelse(is.na(x), 0, x)))
    disturb.binary.2 = disturb.binary.2 %>% filter(!is.na(CYCLONEany), !is.na(BLEACHINGany), !is.na(COTSany)) %>%
        mutate(DISTURBany=ifelse((CYCLONEany+COTSany+BLEACHINGany)>0,1,0),
               DISTURBsevere = ifelse((CYCLONEsevere+COTSsevere+BLEACHINGsevere)>0,1,0),
               time = REPORT_YEAR-min(REPORT_YEAR),
               Zone=factor(Zone, levels=c('Northern GBR','Central GBR','Southern GBR')),
               REEF_NAME=REEF_ID
               )
    save(disturb.binary.2, file='../data/modelled/disturb.binary.2.RData')
    g1=ggplot(disturb.binary, aes(y=DISTURBany, x=REPORT_YEAR)) +
        geom_count(aes(size=..prop.., group=REPORT_YEAR)) + scale_size_area() +
        facet_wrap(~Zone)
    ggsave(file='../output/figures/disturb_any.pdf', g1,width=10, height=3)

    disturb.any.2.glmmTMB = glmmTMB(DISTURBany ~ time*Zone + (1|REEF_NAME),
                                    data=disturb.binary.2, family=binomial())
    save(disturb.any.2.glmmTMB, file='../data/modelled/disturb.any.2.glmmTMB.RData')
    ## ----end
}
## ----end

## ---- DisturbAny
{
    ## ---- Disturb any
    summary(disturb.any.glmmTMB)
    plot(ACF(disturb.any.glmmTMB, resType="pearson"), alpha=0.05)
    calcFreqs(disturb.any.glmmTMB)
    calcFreqs.matrix(disturb.any.glmmTMB)
    ## ----
    df=modelEachReef(mod = disturb.any.glmmTMB, form=DISTURBany~time)
    d1=plotEffects(disturb.any.glmmTMB, firstYear=min(disturb.binary$REPORT_YEAR), individualReefs=df)
    ggsave(file='../output/figures/Disturbances_any.pdf',d1,width=7, height=3)
    ggsave(file='../output/figures/Disturbances_any.png',d1,width=7, height=3, dpi=300)
                                        #ggsave(file='figures/disturb.any.glmmTMB.pdf', g1,width=10, height=3)
    ## ----end

    ## ---- DisturbAny version 2
    summary(disturb.any.2.glmmTMB)
    plot(ACF(disturb.any.2.glmmTMB, resType="pearson"), alpha=0.05)
    calcFreqs(disturb.any.2.glmmTMB)
    calcFreqs.matrix(disturb.any.2.glmmTMB)
    ## ----
    df=modelEachReef(disturb.any.2.glmmTMB,form=DISTURBany~time)
    d1=plotEffects(disturb.any.2.glmmTMB, firstYear=min(disturb.binary.2$REPORT_YEAR), individualReefs=df)
    ggsave(file='../output/figures/Disturbances_any.2.pdf',d1,width=7, height=3)
    ggsave(file='../output/figures/Disturbances_any.2.png',d1,width=7, height=3, dpi=300)
                                        #ggsave(file='figures/disturb.any.glmmTMB.pdf', g1,width=10, height=3)
    ## ----end
    ##polynomials
    ## ---- DisturbAnyPoly
disturb.any.glmmTMB3 = glmmTMB(DISTURBany ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME),
                                data=disturb.binary, family=binomial())
save(disturb.any.glmmTMB3, file='../data/modelled/disturb.any.glmmTMB3.RData')
plot(ACF(disturb.any.glmmTMB3, resType="pearson"), alpha=0.05)
summary(disturb.any.glmmTMB3)
df=modelEachReef(disturb.any.glmmTMB3,form=DISTURBany~poly(time,3),dat=disturb.binary)
d1=plotEffects(disturb.any.glmmTMB3, firstYear=min(disturb.binary$REPORT_YEAR), individualReefs=df)
ggsave(file='../output/figures/Disturbances_any_poly.pdf',d1,width=7, height=3)
ggsave(file='../output/figures/Disturbances_any_poly.png',d1,width=7, height=3, dpi=300)
#ggsave(file='figures/disturb.any.glmmTMB3.pdf', g1,width=10, height=3)
## ----end
    ## ---- DisturbAnyPoly version 2
disturb.any.2.glmmTMB3 = glmmTMB(DISTURBany ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME),
                                data=disturb.binary.2, family=binomial())
save(disturb.any.2.glmmTMB3, file='../data/modelled/disturb.any.2.glmmTMB3.RData')
plot(ACF(disturb.any.2.glmmTMB3, resType="pearson"), alpha=0.05)
summary(disturb.any.2.glmmTMB3)
df=modelEachReef(disturb.any.2.glmmTMB3,form=DISTURBany~poly(time,3),dat=disturb.binary.2)
d1=plotEffects(disturb.any.2.glmmTMB3, firstYear=min(disturb.binary.2$REPORT_YEAR), individualReefs=df)
ggsave(file='../output/figures/Disturbances_any.2_poly.pdf',d1,width=7, height=3)
ggsave(file='../output/figures/Disturbances_any.2_poly.png',d1,width=7, height=3, dpi=300)
#ggsave(file='figures/disturb.any.glmmTMB3.pdf', g1,width=10, height=3)
## ----end
}
## ----end

## ---- Severe disturbances
{
    extract_stats <- function(tab) {
        data.frame(tab$Intercept$Zone,
                   Intecept = sprintf("%0.3g (%0.3g, %0.3g)",tab$Intercept$prob, tab$Intercept$lower.CL, tab$Intercept$upper.CL),
                   Slope = sprintf("%0.3g (%0.3g, %0.3g)",tab$Slope$time.trend, tab$Slope$lower.CL, tab$Slope$upper.CL),
                   AnnualOddsChange = sprintf("%0.3g", 100*(exp(tab$Slope$time.trend)-1)),
                   t.ratio = sprintf("%0.3g", tab$Slope$t.ratio),
                   p.value = sprintf("%0.3g", tab$Slope$p.value)
                   )
    }

    do.call('rbind', lapply(dist.table, extract_stats)) %>%
        rownames_to_column(var="Disturbance") %>%
        mutate(Disturbance = gsub('(.*)\\.[0-9]','\\1',Disturbance)) %>%
        dplyr::rename(Zone=tab.Intercept.Zone)


    ## ---- DisturbSevere
    g1=ggplot(disturb.binary, aes(y=DISTURBsevere, x=REPORT_YEAR)) +
        geom_count(aes(size=..prop.., group=REPORT_YEAR)) + scale_size_area() +
        facet_wrap(~Zone)
    ggsave(file='../output/figures/disturb_severe.pdf', g1,width=10, height=3)

    disturb.severe.glmmTMB = glmmTMB(DISTURBsevere ~ time*Zone + (1|REEF_NAME),
                                     data=disturb.binary, family=binomial())
    save(disturb.severe.glmmTMB, file='../data/modelled/disturb.severe.glmmTMB.RData')
    summary(disturb.severe.glmmTMB)
    plot(ACF(disturb.severe.glmmTMB, resType="pearson"), alpha=0.05)
    calcFreqs(disturb.severe.glmmTMB)
    calcFreqs.matrix(disturb.severe.glmmTMB)
    ## ----
    df=modelEachReef(disturb.severe.glmmTMB,form=DISTURBsevere~time)
    d1=plotEffects(disturb.severe.glmmTMB, firstYear=min(disturb.binary$REPORT_YEAR), individualReefs=df, points=FALSE)
    ## d1=plotEffects(mod = disturb.severe.glmmTMB, firstYear=1985, individualReefs=df, points=FALSE)
    g.severe.any = d1
    save(g.severe.any, file='../data/modelled/g.severe.any.RData')
    ggsave(file='../output/figures/Disturbances_severe.pdf',d1,width=7, height=3)
    ggsave(file='../output/figures/Disturbances_severe.png',d1,width=7, height=3, dpi=300)
                                        #ggsave(file='figures/disturb.severe.glmmTMB.pdf', g1,width=10, height=3)
    ## ----end

    ##polynomials
    ## ---- DisturbSeverePoly
    disturb.severe.glmmTMB3 = glmmTMB(DISTURBsevere ~ poly(time,3,raw=FALSE)*Zone + (1|REEF_NAME),
                                      data=disturb.binary, family=binomial())
    save(disturb.severe.glmmTMB3, file='../data/modelled/disturb.severe.glmmTMB3.RData')
    plot(ACF(disturb.severe.glmmTMB3, resType="pearson"), alpha=0.05)
    summary(disturb.severe.glmmTMB3)
    df=modelEachReef(disturb.severe.glmmTMB3,form=DISTURBsevere~poly(time,3),dat=disturb.binary)
    d1=plotEffects(disturb.severe.glmmTMB3, firstYear=min(disturb.binary$REPORT_YEAR), individualReefs=df)
    ggsave(file='../output/figures/Disturbances_severe_poly.pdf',d1,width=7, height=3)
    ggsave(file='../output/figures/Disturbances_severe_poly.png',d1,width=7, height=3, dpi=300)
                                        #ggsave(file='figures/disturb.severe.glmmTMB3.pdf', g1,width=10, height=3)
    ## ----end
}
## ----end

## ---- Compilation Plots
## Now form a single multipanel figure with the following
## 1. All severe disturbances
## 2. COTS severe disturbances (IO and AO)
## 3. Cyclone severe disturbances (Cat 2,3)
## 4. Bleaching severe disturbances (Cat 3,4)
load(file='../data/modelled/g.severe.any.RData')
load(file='../data/modelled/g.AO.cots.RData')
load(file='../data/modelled/g.severe.cyclones.RData')
load(file='../data/modelled/g.severe.bleaching.RData')

## ---- Disturbances_severe_compilation
gg=grid.arrange(g.severe.any + ggtitle('All severe disturbances'),
             g.severe.cots + ggtitle('A. cf. solaris IO and AO'),
             g.severe.cyclones + ggtitle('Cyclones Hs cat 2 and 3'),
             g.severe.bleaching + ggtitle('Bleaching cat 3 and 4'), ncol=1)

## ggsave(file='output/figures/Disturbances_severe_compilation.pdf',grid.draw(gg),width=9, height=12)
ggsave(file='../output/figures/Disturbances_severe_compilation.pdf',gg,width=9, height=12)
## ggsave(file='output/figures/Disturbances_severe_compilation.png',grid.draw(gg),width=9, height=12, dpi=300)
png(file='../output/figures/Disturbances_severe_compilation.png',width=9, height=12, res=300, units='in')
grid.draw(gg)
dev.off()

## Now form a single multipanel figure with the following
## 1. COTS severe disturbances (IO and AO)
## 2. Cyclone severe disturbances (Cat 2,3)
## 3. Bleaching severe disturbances (Cat 3,4)
gg=grid.arrange(
             g.severe.cots + ggtitle('A. cf. solaris IO and AO'),
             g.severe.cyclones + ggtitle('Cyclones Hs cat 2 and 3'),
             g.severe.bleaching + ggtitle('Bleaching cat 3 and 4'), ncol=1)

## ggsave(file='output/figures/Disturbances_severe_compilation_noAll.pdf',grid.draw(gg),width=9, height=9)
ggsave(file='../output/figures/Disturbances_severe_compilation_noAll.pdf',gg,width=9, height=9)
## ggsave(file='output/figures/Disturbances_severe_compilation_noAll.png',grid.draw(gg),width=9, height=9, dpi=300)
png(file='../output/figures/Disturbances_severe_compilation_noAll.png',width=9, height=9, res=300, units='in')
grid.draw(gg)
dev.off()

## ----end

## ---- Distribuances_severe_full_compilation
## Now combine with disturbances from CoralTrend_trend_manta_3Zone.R
load(file='../data/processed/disturbanceBars-gg.RData')
g.severe.any1=g.severe.any +
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top") + #, limits=c(1985,2021))+
    coord_cartesian(xlim=c(1985, finalYear)) +
    #scale_y_continuous(expression(Pr(impact)),lim=c(0,1.00)) +
    scale_y_continuous(expression(Probability~of~being~impacted~phantom("(")),lim=c(0,1.00)) +
    theme(
        panel.border=element_rect(fill=NA,color='black'),
        axis.title.y=element_text(size=rel(1.5),margin=margin(l=0.5,r=0.5,unit='lines')),
                                        #axis.text.x=element_text(size=rel(1.2)),
        axis.text.x=element_blank(), axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(0,5,5,0),'pt'),
        strip.text.x=element_blank(),
        strip.background=element_blank())
gg2 = ggplot_gtable(ggplot_build(g.severe.any1))
## oldgg = gg
## oldgg2 = gg2
gg$widths
gg2$widths

## gg2$widths[6] <- unit(0, 'cm')


#gg3=gg2

#gg2=gg3
gg$widths[3] =  gg2$widths[3] = unit(1,"cm")
gg$widths[4] =  gg2$widths[4]
## gg$widths[7] <- unit(0,'points')
#gg2$widths=gg$widths
gg2$widths[6] =  gg$widths[7]
gg2$widths[8] =  gg$widths[7] 
grid.arrange(gg,gg2)



## ggsave(file='output/figures/Disturbances_severe_full_compilation.pdf',grid.arrange(gg,gg2),width=9, height=6)
ggsave(file='../output/figures/Disturbances_severe_full_compilation.pdf',grid.arrange(gg,gg2),width=12, height=8)
## ggsave(file='output/figures/Disturbances_severe_full_compilation.png',grid.arrange(gg,gg2),width=9, height=6, dpi=300)
png(file='../output/figures/Disturbances_severe_full_compilation.png',width=12, height=8, res=300, units='in')
    grid.arrange(gg,gg2)
dev.off()
## ----end
## ---- Disturbances_severe_compilation
##Finally, a figure that combines the disturbance bars from CoralTrends_trend_manta_3Zone.R with the severe version of each major disturbance.
load(file='../data/processed/disturbanceBars-gg.RData')
g.severe.cots1= g.severe.cots + ggtitle('(c) A. cf. solaris IO and AO') +
    #scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top", limits=c(1985,2020))+
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "bottom") + #, limits=c(1985,2021))+
    coord_cartesian(xlim=c(1985, finalYear)) +
    #scale_y_continuous(expression(Probability~of~being~impacted~phantom("(")),lim=c(0,1.00)) +
    scale_y_continuous(expression(phantom("(")),lim=c(0,1.00)) +
    theme(
        panel.border=element_rect(fill=NA,color='black'),
        #axis.title.y=element_text(size=rel(1.5),margin=margin(r=0.5,unit='lines')),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=rel(1.0)),
        #axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(0,5,5,6),'pt'),
        strip.text.x=element_blank(),
        strip.background=element_blank())
g.severe.cots2 = ggplot_gtable(ggplot_build(g.severe.cots1))

g.severe.cots1= g.AO.cots + ggtitle(expression(paste('(c) ',italic(A.~cf.~solaris), ' AO'))) +
    #scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top", limits=c(1985,2020))+
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "bottom") + # ,limits=c(1985,2021))+
    coord_cartesian(xlim=c(1985, finalYear)) +
    #scale_y_continuous(expression(Probability~of~being~impacted~phantom("(")),lim=c(0,1.00)) +
    scale_y_continuous(expression(phantom("(")),lim=c(0,1.00)) +
    theme(
        panel.border=element_rect(fill=NA,color='black'),
        #axis.title.y=element_text(size=rel(1.5),margin=margin(r=0.5,unit='lines')),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=rel(1.0)),
        #axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(0,5,5,6),'pt'),
        strip.text.x=element_blank(),
        strip.background=element_blank())
g.severe.cots2 = ggplot_gtable(ggplot_build(g.severe.cots1))


g.severe.cyclones1= g.severe.cyclones + ggtitle('(d) Cyclones Hs cat 2 and 3') +
    #scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top", limits=c(1985,2020))+
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "bottom") + #, limits=c(1985,2021))+
    coord_cartesian(xlim=c(1985, finalYear)) +
    scale_y_continuous(expression(Probability~of~being~impacted~phantom("(")),lim=c(0,1.00)) +
    theme(
        panel.border=element_rect(fill=NA,color='black'),
        #axis.title.y=element_text(size=rel(1.5),margin=margin(r=0.5,unit='lines')),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=rel(1.0)),
        #axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(0,5,5,6),'pt'),
        strip.text.x=element_blank(),
        strip.background=element_blank())
g.severe.cyclones2 = ggplot_gtable(ggplot_build(g.severe.cyclones1))

g.severe.bleaching1= g.severe.bleaching + ggtitle('(b) Bleaching cat 3, 4 and 5') +
    #scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top", limits=c(1985,2020))+
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "bottom") + #, limits=c(1985,2021))+
    coord_cartesian(xlim=c(1985, finalYear)) +
    #scale_y_continuous(expression(Probability~of~being~impacted~phantom("(")),lim=c(0,1.00)) +
    scale_y_continuous(expression(phantom("(")),lim=c(0,1.00)) +
    theme(
        panel.border=element_rect(fill=NA,color='black'),
                                        #axis.title.y=element_text(size=rel(1.5),margin=margin(r=0.5,unit='lines')),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=rel(1.0)),
        #axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(0,5,5,6),'pt'),
        strip.text.x=element_blank(),
      strip.background=element_blank())

g.severe.bleaching2 = ggplot_gtable(ggplot_build(g.severe.bleaching1))
gg$widths[4] =  g.severe.bleaching2$widths[4]
g.severe.bleaching2$widths[6] =  gg$widths[7]
g.severe.bleaching2$widths[8] =  gg$widths[7]
g.severe.cots2$widths[6] =  gg$widths[7]
g.severe.cots2$widths[8] =  gg$widths[7]
g.severe.cyclones2$widths[6] =  gg$widths[7]
g.severe.cyclones2$widths[8] =  gg$widths[7]

gg1=grid.arrange(
    g.severe.bleaching2,
    g.severe.cots2,
    g.severe.cyclones2,
    ncol=1,
    left=textGrob(expression(Probability~of~being~impacted~phantom("(")), vjust=.375,rot=90, gp=gpar(fontsize=16)))
    #left=textGrob(expression(Probability~of~being~impacted), rot=90, gp=gpar(fontsize=16)))
    #left=expression(Probability~of~being~impacted~phantom("(")))
                                        #left='Probability of being impacted')

grid.arrange(gg,gg1, heights=c(1,3))
ggsave(file='../output/figures/Disturbances_severe_compilation.pdf',grid.arrange(gg,gg1, heights=c(1,2)),width=9, height=8)
ggsave(file='../output/figures/Disturbances_severe_compilation.png',grid.arrange(gg,gg1, heights=c(1,2)),width=9, height=8, dpi=300)

## Replace the ggtile in the gTable of gg
matches <- grep(pattern = '^title', gg$layout$name)
gg$grobs[[matches]]$children[[1]]$label = "(a)"
grid.arrange(gg,gg1, heights=c(1,3))
ggsave(file='../output/figures/Disturbances_severe_compilationNewLabels.pdf',grid.arrange(gg,gg1, heights=c(1,2)),width=9, height=8)
ggsave(file='../output/figures/Disturbances_severe_compilationNewLabels.png',grid.arrange(gg,gg1, heights=c(1,2)),width=9, height=8, dpi=300)
png(file='../output/figures/Disturbances_severe_compilationNewLabels.png',width=9, height=8, res=300, units='in')
    grid.arrange(gg,gg1, heights=c(1,2))
dev.off()
## ----end

## Split up figure for Cell Press format (2020) ---------------------------------
## ---- Disturbances_severe_compilation_new.pdf
load(file='../data/processed/disturbanceBars-gg.RData')
load(file='../data/spatial/gts.RData')
g.severe.cots1= g.AO.cots +
  geom_text(data=data.frame(x=1985,y=1,Zone=factor('Northern GBR', levels=c('Northern GBR', 'Central GBR', 'Southern GBR')), l = 'B'),
            aes(y=y,x=x, label=l), hjust=0, vjust=1) +   #ggtitle(expression(paste('(c) ',italic(A.~cf.~solaris), ' AO'))) +
    #scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top", limits=c(1985,2020))+
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "bottom") + #, limits=c(1985,2020))+
    coord_cartesian(xlim=c(1985, finalYear)) +
    #scale_y_continuous(expression(Probability~of~being~impacted~phantom("(")),lim=c(0,1.00)) +
    scale_y_continuous(expression(phantom("(")),lim=c(0,1.00)) +
    theme(
        panel.border=element_rect(fill=NA,color='black'),
        #axis.title.y=element_text(size=rel(1.5),margin=margin(r=0.5,unit='lines')),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=rel(1.0)),
        #axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(0,5,5,6),'pt'),
        strip.text.x=element_blank(),
      strip.background=element_blank(),
      plot.title.position = 'panel')
g.severe.cots2 = ggplot_gtable(ggplot_build(g.severe.cots1))

g.severe.cyclones1= g.severe.cyclones +
  geom_text(data=data.frame(x=1985,y=1,Zone=factor('Northern GBR', levels=c('Northern GBR', 'Central GBR', 'Southern GBR')), l = 'C'),
            aes(y=y,x=x, label=l), hjust=0, vjust=1) +   #ggtitle('(d) Cyclones Hs cat 2 and 3') +
    #scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top", limits=c(1985,2020))+
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "bottom") + #, limits=c(1985,2020))+
    coord_cartesian(xlim=c(1985, finalYear)) +
    scale_y_continuous(expression(Probability~of~being~impacted~phantom("(")),lim=c(0,1.00)) +
    theme(
        panel.border=element_rect(fill=NA,color='black'),
        #axis.title.y=element_text(size=rel(1.5),margin=margin(r=0.5,unit='lines')),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=rel(1.0)),
        #axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(0,5,5,6),'pt'),
        strip.text.x=element_blank(),
        strip.background=element_blank())
g.severe.cyclones2 = ggplot_gtable(ggplot_build(g.severe.cyclones1))

g.severe.bleaching1= g.severe.bleaching +
    geom_text(data=data.frame(x=1985,y=1,Zone=factor('Northern GBR', levels=c('Northern GBR', 'Central GBR', 'Southern GBR')), l = 'A'),
            aes(y=y,x=x, label=l), hjust=0, vjust=1) +   #ggtitle('(b) Bleaching cat 3 and 4') +
    #scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top", limits=c(1985,2020))+
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "bottom") +#, limits=c(1985,2020))+
    coord_cartesian(xlim=c(1985, finalYear)) +
    #scale_y_continuous(expression(Probability~of~being~impacted~phantom("(")),lim=c(0,1.00)) +
    scale_y_continuous(expression(phantom("(")),lim=c(0,1.00)) +
    theme(
        panel.border=element_rect(fill=NA,color='black'),
                                        #axis.title.y=element_text(size=rel(1.5),margin=margin(r=0.5,unit='lines')),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=rel(1.0)),
        #axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(0,5,5,6),'pt'),
        strip.text.x=element_blank(),
      strip.background=element_blank())
g.severe.bleaching2 = ggplot_gtable(ggplot_build(g.severe.bleaching1))

gg$widths[4] =  g.severe.bleaching2$widths[4]
g.severe.bleaching2$widths[6] =  gg$widths[7]
g.severe.bleaching2$widths[8] =  gg$widths[7]
g.severe.cots2$widths[6] =  gg$widths[7]
g.severe.cots2$widths[8] =  gg$widths[7]
g.severe.cyclones2$widths[6] =  gg$widths[7]
g.severe.cyclones2$widths[8] =  gg$widths[7]

## Put the big strip on
g.severe.bleaching3 = g.severe.bleaching2
facets <- grep("strip-t-1", g.severe.bleaching2$layout$name)
facets1 <- grep("strip-t-1", gg$layout$name)
g.severe.bleaching3$grobs[facets] = gg$grobs[facets1]
g.severe.bleaching3$grobs[1] = gg$grobs[1]
g.severe.bleaching3 <- with(g.severe.bleaching3$layout[facets,],
           gtable_add_grob(g.severe.bleaching3, ggplotGrob(gt1),t=t, l=5, b=b, r=5, name="pic_predator"))
facets <- grep("strip-t-2", g.severe.bleaching2$layout$name)
facets1 <- grep("strip-t-2", gg$layout$name)
g.severe.bleaching3$grobs[facets] = gg$grobs[facets1]
g.severe.bleaching3 <- with(g.severe.bleaching3$layout[facets,],
           gtable_add_grob(g.severe.bleaching3, ggplotGrob(gt2),t=t, l=7, b=b, r=7, name="pic_predator"))
facets <- grep("strip-t-3", g.severe.bleaching2$layout$name)
facets1 <- grep("strip-t-3", gg$layout$name)
g.severe.bleaching3$grobs[facets] = gg$grobs[facets1]
g.severe.bleaching3 <- with(g.severe.bleaching3$layout[facets,],
           gtable_add_grob(g.severe.bleaching3, ggplotGrob(gt3),t=t, l=9, b=b, r=9, name="pic_predator"))

g.severe.bleaching3$heights[7] = gg$heights[7]
#grid.newpage()
#grid.draw(g.severe.bleaching3)
gg1=grid.arrange(
  g.severe.bleaching3,
  g.severe.cots2,
  g.severe.cyclones2,
  ncol=1,
  heights=c(1.35,1,1),
  ## left=textGrob(expression(Probability~of~being~impacted~phantom("(")), vjust=.375,rot=90, gp=gpar(fontsize=16))
  left=textGrob(expression(Probability~of~being~impacted~phantom("(")), just="centre", rot=90, gp=gpar(fontsize=16))
  )
grid.arrange(gg,gg1, heights=c(1,3))
grid.arrange(gg1)
## ggsave(file='output/figures/Disturbances_severe_compilation_new.pdf',grid.draw(gg1),width=9, height=6)
ggsave(file='../output/figures/Disturbances_severe_compilation_new.pdf',gg1,width=9, height=6)
ggsave(file='../output/figures/Disturbances_severe_compilation_new.png',grid.draw(gg1),width=9, height=6, dpi=300)
png(file='../output/figures/Disturbances_severe_compilation_new.png',width=9, height=6, res=300, units='in')
grid.draw(gg1)
dev.off()


## Or try patchwork
library(patchwork)
gg1a <- patchwork::wrap_plots(textGrob(expression(Probability~of~being~impacted~phantom("(")), just="centre", rot=90, gp=gpar(fontsize=16)),
                      wrap_plots(g.severe.bleaching3, g.severe.cots2, g.severe.cyclones2, ncol=1), widths=c(0.2, 10), ncol=2)
                                        #left=textGrob(expression(Probability~of~being~impacted), rot=90, gp=gpar(fontsize=16)))
    #left=expression(Probability~of~being~impacted~phantom("(")))
                                        #left='Probability of being impacted')

ggsave(file='../output/figures/Disturbances_severe_compilation_newA.pdf',gg1a,width=9, height=6)
ggsave(file='../output/figures/Disturbances_severe_compilation_newA.png',gg1a,width=9, height=6, dpi=300)
png(file='../output/figures/Disturbances_severe_compilation_newA.png',width=9, height=6, res=300, units='in')
gg1a
dev.off()
## ----end


## Transpose version (2021) ---------------------------------
## ---- Disturbances_severe_compilation_new_transposed.pdf
## Need to redo most plots
load(file='../data/processed/disturbanceBars-gg.RData')
load(file='../data/spatial/gts.RData')
hues <- RColorBrewer::brewer.pal(4, "Blues")


g.severe.cots1= g.AO.cots_2021 +
    ## geom_text(data=data.frame(x=1985,y=1,Zone=factor('Northern GBR', levels=c('Northern GBR', 'Central GBR', 'Southern GBR')), l = '(b)'),
    geom_text(data=data.frame(x=1985,y=1,Zone=factor('Northern GBR', levels=c('Northern GBR', 'Central GBR', 'Southern GBR')), l = '(b)'),
            aes(y=y,x=x, label=l), hjust=0, vjust=1) +   #ggtitle(expression(paste('(c) ',italic(A.~cf.~solaris), ' AO'))) +
    #scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top", limits=c(1985,2020))+
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "bottom", limits=c(1984.5,(finalYear+1)))+
    coord_cartesian(xlim=c(1985, finalYear)) +
    #scale_y_continuous(expression(Probability~of~being~impacted~phantom("(")),lim=c(0,1.00)) +
    scale_y_continuous(expression(phantom("(")),lim=c(0,1.00)) +
    facet_grid(Zone ~Col, scales='fixed', labeller=labeller(Zone=setNames(paste0("", labs.shorter, "\n"), labs))) +
    theme(
        panel.border=element_rect(fill=NA,color='black'),
        #axis.title.y=element_text(size=rel(1.5),margin=margin(r=0.5,unit='lines')),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=rel(1.0)),
        #axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(0,5,5,6),'pt'),
        strip.text.x=element_blank(),
        strip.text=element_text(margin=margin(l=0.5, r=0.5,unit='lines'),size=20,lineheight=0.5, face='bold',hjust=0.9,vjust=0),
      ## strip.background=element_blank(),
      strip.background=element_rect(fill=hues[2], color='black', size=0.5),
      plot.title.position = 'panel')
g.severe.cots1

g.severe.cots2 = ggplot_gtable(ggplot_build(g.severe.cots1))
g.severe.cots2 %>% grid.draw()


g.severe.cyclones1= g.severe.cyclones_2021 +
    geom_text(data=data.frame(x=1985,y=1,Zone=factor('Northern GBR', levels=c('Northern GBR', 'Central GBR', 'Southern GBR')), l = '(c)'),
              aes(y=y,x=x, label=l), hjust=0, vjust=1) +   #ggtitle('(d) Cyclones Hs cat 2 and 3') +
                                        #scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top", limits=c(1985,2020))+
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "bottom") + #, limits=c(1985,2021))+
    coord_cartesian(xlim=c(1985, finalYear)) +
    scale_y_continuous(expression(Probability~of~being~impacted~phantom("(")),lim=c(0,1.00)) +
    ## facet_grid(Zone ~Col, scales='fixed') +
    facet_grid(Zone ~Col, scales='fixed', labeller=labeller(Zone=setNames(paste0("", labs.shorter, "\n"), labs))) +
    theme(
        panel.border=element_rect(fill=NA,color='black'),
                                        #axis.title.y=element_text(size=rel(1.5),margin=margin(r=0.5,unit='lines')),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=rel(1.0)),
                                        #axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(0,5,5,6),'pt'),
        ## strip.text.x=element_blank(),
        strip.text=element_text(margin=margin(l=0.5, r=0.5,unit='lines'),size=20,lineheight=0.5, face='bold',hjust=0.9,vjust=0),
        strip.text.x=element_blank(),
        ## strip.background=element_blank(),
        strip.background=element_rect(fill=hues[2], color='black', size=0.5),
        plot.title.position = 'panel')
g.severe.cyclones2 = ggplot_gtable(ggplot_build(g.severe.cyclones1))
g.severe.cyclones2 %>% grid.draw()

g.severe.bleaching1= g.severe.bleaching_2021 +
  geom_text(data=data.frame(x=1985,y=1,Zone=factor('Northern GBR', levels=c('Northern GBR', 'Central GBR', 'Southern GBR')), l = '(a)'),
            aes(y=y,x=x, label=l), hjust=0, vjust=1) +   #ggtitle('(b) Bleaching cat 3 and 4') +
    #scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top", limits=c(1985,2020))+
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "bottom") + #, limits=c(1985,2021))+
    coord_cartesian(xlim=c(1985, finalYear)) +
    #scale_y_continuous(expression(Probability~of~being~impacted~phantom("(")),lim=c(0,1.00)) +
    scale_y_continuous(expression(phantom("(")),lim=c(0,1.00)) +
    ## facet_grid(Zone ~Col, scales='fixed') +
    facet_grid(Zone ~Col, scales='fixed', labeller=labeller(Zone=setNames(paste0("", labs.shorter, "\n"), labs))) +
    theme(
        panel.border=element_rect(fill=NA,color='black'),
                                        #axis.title.y=element_text(size=rel(1.5),margin=margin(r=0.5,unit='lines')),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=rel(1.0)),
                                        #axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(0,5,5,6),'pt'),
        ## strip.text.x=element_blank(),
        strip.text=element_text(margin=margin(l=0.5, r=0.5,unit='lines'),size=20,lineheight=0.5, face='bold',hjust=0.9,vjust=0),
        strip.text.x=element_blank(),
        ## strip.background=element_blank(),
        strip.background=element_rect(fill=hues[2], color='black', size=0.5),
        plot.title.position = 'panel')
g.severe.bleaching2 = ggplot_gtable(ggplot_build(g.severe.bleaching1))
g.severe.bleaching2 %>% grid.draw()

## gg$widths[4] =  g.severe.bleaching2$widths[4]
## g.severe.bleaching2$widths[6] =  gg$widths[7]
## g.severe.bleaching2$widths[8] =  gg$widths[7]
## g.severe.cots2$widths[6] =  gg$widths[7]
## g.severe.cots2$widths[8] =  gg$widths[7]
## g.severe.cyclones2$widths[6] =  gg$widths[7]
## g.severe.cyclones2$widths[8] =  gg$widths[7]

## Put the big strip on
g.severe.bleaching3 = g.severe.bleaching2
panels_src <- grep("panel-1-1", g.severe.cots2$layout$name)
panels_target <- grep("panel-1-2", g.severe.bleaching3$layout$name)
g.severe.bleaching3$grobs[panels_target] <- g.severe.cots2$grobs[panels_src]
panels_src <- grep("panel-1-1", g.severe.cyclones2$layout$name)
panels_target <- grep("panel-1-3", g.severe.bleaching3$layout$name)
g.severe.bleaching3$grobs[panels_target] <- g.severe.cyclones2$grobs[panels_src]

panels_src <- grep("panel-2-1", g.severe.cots2$layout$name)
panels_target <- grep("panel-2-2", g.severe.bleaching3$layout$name)
g.severe.bleaching3$grobs[panels_target] <- g.severe.cots2$grobs[panels_src]
panels_src <- grep("panel-2-1", g.severe.cyclones2$layout$name)
panels_target <- grep("panel-2-3", g.severe.bleaching3$layout$name)
g.severe.bleaching3$grobs[panels_target] <- g.severe.cyclones2$grobs[panels_src]

panels_src <- grep("panel-3-1", g.severe.cots2$layout$name)
panels_target <- grep("panel-3-2", g.severe.bleaching3$layout$name)
g.severe.bleaching3$grobs[panels_target] <- g.severe.cots2$grobs[panels_src]
panels_src <- grep("panel-3-1", g.severe.cyclones2$layout$name)
panels_target <- grep("panel-3-3", g.severe.bleaching3$layout$name)
g.severe.bleaching3$grobs[panels_target] <- g.severe.cyclones2$grobs[panels_src]
g.severe.bleaching3 %>% grid.draw()
facets <- grep("strip-r-1", g.severe.bleaching3$layout$name)
## facets1 <- grep("strip-t-1-1", gg$layout$name)
## g.severe.bleaching3$grobs[facets] = gg$grobs[facets1]
## g.severe.bleaching3$grobs[1] = gg$grobs[1]
g.severe.bleaching3 <- with(g.severe.bleaching3$layout[facets,],
           gtable_add_grob(g.severe.bleaching3, ggplotGrob(gt1_2021),t=6, l=10, b=8, r=10, name="pic_predator"))
facets <- grep("strip-r-2", g.severe.bleaching3$layout$name)
## facets1 <- grep("strip-r-2", gg$layout$name)
## g.severe.bleaching3$grobs[facets] = gg$grobs[facets1]
g.severe.bleaching3 <- with(g.severe.bleaching3$layout[facets,],
           gtable_add_grob(g.severe.bleaching3, ggplotGrob(gt2_2021),t=10, l=10, b=10, r=10, name="pic_predator"))
facets <- grep("strip-r-3", g.severe.bleaching3$layout$name)
## facets1 <- grep("strip-r-3", gg$layout$name)
## g.severe.bleaching3$grobs[facets] = gg$grobs[facets1]
g.severe.bleaching3 <- with(g.severe.bleaching3$layout[facets,],
           gtable_add_grob(g.severe.bleaching3, ggplotGrob(gt3_2021),t=12, l=10, b=12, r=10, name="pic_predator"))
g.severe.bleaching3 %>% grid.draw()

gg1=grid.arrange(
  g.severe.bleaching3,
  ncol=1,
  left=textGrob(expression(Probability~of~being~impacted~phantom("(")), just="centre", rot=90, gp=gpar(fontsize=16))
  )

## Shorten the strip lables
gg1

grid.newpage()
grid.draw(gg1)
## g.severe.bleaching3$heights[7] = gg$heights[7]
## #grid.newpage()
## #grid.draw(g.severe.bleaching3)
## gg1=grid.arrange(
##   g.severe.bleaching3,
##   g.severe.cots2,
##   g.severe.cyclones2,
##   ncol=1,
##   heights=c(1.35,1,1),
##   ## left=textGrob(expression(Probability~of~being~impacted~phantom("(")), vjust=.375,rot=90, gp=gpar(fontsize=16))
##   left=textGrob(expression(Probability~of~being~impacted~phantom("(")), just="centre", rot=90, gp=gpar(fontsize=16))
##   )
## grid.arrange(gg,gg1, heights=c(1,3))
## grid.arrange(gg1)
## ggsave(file='output/figures/Disturbances_severe_compilation_new_2021.pdf',grid.draw(gg1),width=9, height=10)
ggsave(file='../output/figures/Disturbances_severe_compilation_new_2021.pdf',gg1,width=9, height=8)
## ggsave(file='output/figures/Disturbances_severe_compilation_new_2021.png',grid.draw(gg1),width=9, height=10, dpi=300)
png(file='../output/figures/Disturbances_severe_compilation_new_2021.png',width=9, height=8, res=300, units='in')
grid.draw(gg1)
dev.off()

pdf(file='../output/figures/Disturbances_severe_compilation_new_2021.pdf',width=9, height=8)
grid.draw(gg1)
dev.off()
png(file='../output/figures/Disturbances_severe_compilation_new_2021.png',width=9, height=8, units = 'in', res = 300)
grid.draw(gg1)
dev.off()


## Or try patchwork
library(patchwork)
gg1a <- patchwork::wrap_plots(textGrob(expression(Probability~of~being~impacted~phantom("(")), just="centre", rot=90, gp=gpar(fontsize=16)),
                      wrap_plots(g.severe.bleaching3, ncol=1), widths=c(0.2, 10), ncol=2)
                                        #left=textGrob(expression(Probability~of~being~impacted), rot=90, gp=gpar(fontsize=16)))
    #left=expression(Probability~of~being~impacted~phantom("(")))
                                        #left='Probability of being impacted')

ggsave(file='../output/figures/Disturbances_severe_compilation_newA_2021.pdf',gg1a,width=9, height=7)
ggsave(file='../output/figures/Disturbances_severe_compilation_newA_2021.png',gg1a,width=9, height=7, dpi=300)
## ----end



## Graphical Abstract Figure-----------------------------------------------------
## ---- GraphicalAbstract
load(file='../data/modelled/g.severe.any.RData')
load(file='../data/modelled/g.AO.cots.RData')
load(file='../data/modelled/g.severe.cyclones.RData')
load(file='../data/modelled/g.severe.bleaching.RData')
#load(file='data/spatial/gg_3Zone.RData')
load(file='../data/modelled/g1_trends.Rmd')
gg.trend=g1 + theme(strip.text=element_blank(), plot.margin=unit(c(0,5,5,6),'pt')) + ggtitle('(e)')
gg.trend2 = ggplot_gtable(ggplot_build(gg.trend))
gg.trend2$widths=  gg$widths

matches <- grep(pattern = '^title', gg$layout$name)
gg$grobs[[matches]]$children[[1]]$label = "(a)"
matches <- grep(pattern = '^title', g.severe.bleaching2$layout$name)
g.severe.bleaching2$grobs[[matches]]$children[[1]]$label = "(b)"
matches <- grep(pattern = '^title', g.severe.cots2$layout$name)
g.severe.cots2$grobs[[matches]]$children[[1]]$label = "(c)"
matches <- grep(pattern = '^title', g.severe.cyclones2$layout$name)
g.severe.cyclones2$grobs[[matches]]$children[[1]]$label = "(d)"

gg.dist=grid.arrange(
  g.severe.bleaching2, # + ggtitle('Bleaching cat 3 and 4'), ncol=1)
  g.severe.cots2,# + ggtitle('A. cf. solaris IO and AO'),
  g.severe.cyclones2,# + ggtitle('Cyclones Hs cat 2 and 3'),
  left=textGrob(expression(Probability~of~being~impacted~phantom("(")), vjust=.375,rot=90, gp=gpar(fontsize=16)))


grid.arrange(gg,gg.dist,gg.trend2, heights=c(1.3,3,1))
ggsave(file='../output/figures/GraphicalAbstract.pdf',grid.arrange(gg,gg.dist,gg.trend2, heights=c(1.3,3,1)), width=9, height=8)
ggsave(file='../output/figures/GraphicalAbstract.png',grid.arrange(gg,gg.dist,gg.trend2, heights=c(1.3,3,1)), width=9, height=8, dpi=300)
png(file='../output/figures/GraphicalAbstract.png', width=9, height=8, res=300, units='in')
    grid.arrange(gg,gg.dist,gg.trend2, heights=c(1.3,3,1))
dev.off()
## ----end

## ----end

## Zip up some stuff for Mike
zip('output/DisturbanceFrequencyFigures4Mike.zip', list.files('output/figures', pattern='Disturbances.*.p..|^GraphicalAbstract', full.names=TRUE))

load(file='data/processed/disturbanceBars-gg.RData')
g.severe.any1=g.severe.any +
    scale_x_continuous('',breaks=seq(1985,2020,by=5),position = "top", limits=c(1985,2020))+
    #scale_y_continuous(expression(Pr(impact)),lim=c(0,1.00)) +
    scale_y_continuous(expression(Probability~of~being~impacted~phantom("(")),lim=c(0,1.00)) +
    theme(
        panel.border=element_rect(fill=NA,color='black'),
        axis.title.y=element_text(size=rel(1.5),margin=margin(r=0.5,unit='lines')),
                                        #axis.text.x=element_text(size=rel(1.2)),
        axis.text.x=element_blank(), axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.2)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0.1,color=NA,linetype='dashed'),
        panel.grid.major.x=element_line(size=0.1,color='gray70',linetype='dashed'),
        plot.margin=unit(c(0,5,5,0),'pt'),
        strip.text.x=element_blank(),
        strip.background=element_blank())
gg2 = ggplot_gtable(ggplot_build(g.severe.any1))
    



