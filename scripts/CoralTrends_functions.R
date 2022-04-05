###################################################################
## The following function checks to ensure that all the required ##
## packages are available on the system.                         ##
###################################################################
CoralTrends_checkPackages <- function() {
    require(gdata) # load this first as it masks many
    require(tidyverse)
    require(gtable)
    require(grid)
    require(gridExtra)
    require(xtable)
    library(broom)
    require(rgdal)
    require(rgeos)
    require(sp)
    require(oz)
    require(maps)
    require(mapdata)
    require(ggsn)
    require(scales)
    require(mapping) # consider replacing this with a self contained function
    require(maptools)
    require(raster)

    require(INLA)
    require(rstanarm)
    require(coda)
    require(sf)
    require(glmmTMB)
}

######################################################################
## The following function is used when issuing a system call from R ##
## (e.g. running xelatex).  It ensures that warnings are captured.  ##
######################################################################
CoralTrends_system <- function (sys.command) {
    ret.val<- system (sys.command, intern=TRUE)

    if (!is.null(attributes(ret.val)$status)) {
        warning (paste('CoralTrends_WARNING', sys.command, " failed | ",paste(ret.val, collapse='\n')))
    }
    ret.val
}


#########################################################################
## The following function appends a log statement into a log file      ##
## parameters:                                                         ##
##     status:    a string indicating either 'FAILURE', 'SUCCESS' or   ##
##                'WARNING'                                            ##
##     logFile:   a character string representation of the log file    ##
##                name (including path relative to current working     ##
##                directory)                                           ##
##     Category:  a character string with a category to appear         ##
##                verbatim in the log                                  ##
##     msg1:      a character string with a message to appear verbatim ##
##                in the log                                           ##
#########################################################################
CoralTrends_log <- function (status, logFile='data/logs/env.log',Category, msg1) {
    options(digits.secs=2)              ## switch to subsecond display
    ## Check if the log file exists, and if it does not, create it
    d=dirname(logFile)
    files <- list.files(d)
    if(!any(grepl(paste0('^',logFile,'$'),files))) system(paste0('touch ',logFile))
    now <- Sys.time()

    msg <- paste0(now, '|',status, ': ', Category, ' ',msg1)
    if( !is.null(msg)){ write(msg,file=paste0(logFile),append=TRUE)}

}

CoralTrends_tryCatch <- function(expr, logFile,Category, expectedClass=NULL, msg=NULL, return=NULL, showWarnings=FALSE) {
    #msg <- paste0(now, '| ', msg)
    max.warnings<-4
    warnings<-0
    W <- NULL
    w.handler <- function(w){ # warning handler
        m<-w$message
        if ((warnings < max.warnings) && (grepl ('CoralTrends_WARNING', m)>0)) {
            CoralTrends_log('WARNING', logFile,Category, paste(warnings, msg, m))
            warnings<<-warnings+1
        }
        invokeRestart("muffleWarning")
    }
    ret <- list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                    warning = w.handler),warning = W)
    if(!is.atomic(ret$value) && !is.null(ret$value$message)){
        ## An error occurred
        class(ret) <- "try-error"
        CoralTrends_log('WARNING', logFile,Category, paste(msg, ret$value$message))
        if(!is.null(return)) {
            FALSE
        }#else return()
    } else {    #no error check for warning
        CoralTrends_log('INFO', logFile, Category, msg)
        if(!is.null(return)) {
            TRUE
        }
    }
}


CoralTrends_calcPercent = function(x) {
    ifelse(x=='0', 0,
    ifelse(x=='1', 0.05,
    ifelse(x=='1L', 0.025,
    ifelse(x=='1U', 0.075,
    ifelse(x=='2', 0.2,
    ifelse(x=='2L', 0.15,
    ifelse(x=='2U', 0.25,
    ifelse(x=='3', 0.4,
    ifelse(x=='3L', 0.35,
    ifelse(x=='3U', 0.45,
    ifelse(x=='4', 0.625,
    ifelse(x=='4L', 0.5625,
    ifelse(x=='4U', 0.6875,
    ifelse(x=='5', 0.875,
    ifelse(x=='5L',0.8125,0.9375)))))))))))))))
}

COTScategories = function(x) {
    case_when(
        x == 0 ~ 'Zero',
        x > 0 & x < 0.22 ~ 'NO',
        x >= 0.22 & x < 1 ~ 'IO',
        x >= 1 ~ 'AO'
    )
}

ML_gClip <- function(shp, bb){
    if(class(bb) == "matrix") {
        if (identical(dim(bb), c(2L,2L))) {
            b_poly <- as(raster:::extent(as.vector(t(bb))), "SpatialPolygons")
        } else b_poly = bb
    } else if (class(bb) =='SpatialPolygons') {
        b_poly=bb
    } else b_poly <- as(raster:::extent(bb), "SpatialPolygons")
  gIntersection(shp, b_poly, byid = T)
}


######################################################################
## The following function generates Zone categories based on De'ath ##
## 2012's latitudinal divisions.                                    ##
##   parameters:                                                    ##
##      x:     a numeric vector of latitudes                        ##
##   returns:  a categorical vector of Zones                        ##
######################################################################
CoralTrends_calc3ZoneLocations <- function(x) {
    factor(ifelse(x<= -10.68 & x > -15.4, 'Northern',  #glenn's version is -11.8
           ifelse(x<= -15.4 & x > -20.0, 'Central',
           ifelse(x<= -20.0 & x > -23.92, 'Southern','Outside'))))
}
CoralTrends_calc3ZoneLocation <- function(dat) {
    load(paste0(DATA_PATH, 'primary/gbr_3Zone.RData'))
    dat %>%
        st_as_sf(coords = c("Longitude", "Latitude"), crs = st_crs(gbr_3Zone)) %>%
        st_join(gbr_3Zone) %>%
        cbind(Longitude=st_coordinates(.)[,1],Latitude=st_coordinates(.)[,2]) %>%
        st_drop_geometry
}

## The following function converts a cover into a category representing the median of the category
median_cover_cat <- function(dat) {
    n <- length(dat)
    if (n %% 2 == 0) {
        d <- paste(unique(sort(dat)[c((n-1)/2, (n+1)/2)]), collapse='/')
    } else {
        d <- paste(unique(sort(dat)[(n+1)/2]))
    }
    factor(d)
}

################################################################
## The following function converts a list with x and y into a ##
## data.frame                                                 ##
##   parameters:                                              ##
##      xy:    a list with elements x and y                   ##
##   returns:  a data.frame with fields x and y               ##
################################################################
xy2df<-function(xy) {
  data.frame(x=xy$x,y=xy$y)
}


towns12 <-
structure(list(town = c("Cooktown", "Cairns", "Townsville", "Bowen",
"Proserpine", "Mackay", "Rockhampton", "Gladstone", "Port Douglas",
"Ayr", "Innisfail", "Cardwell"), long = c(145.229629516602, 145.740203857422,
146.78092956543, 148.190048217773, 148.696716308594, 149.112060546875,
150.409729003906, 151.094924926758, 145.405410766602, 147.405349731445,
146.016479492188, 146.024261474609), lat = c(-15.4718799591065,
-16.9359798431396, -19.2582454681396, -19.9939575195313, -20.311653137207,
-21.1522636413574, -23.3642120361328, -23.8742027282715, -16.4558181762695,
-19.5710220336914, -17.516508102417, -18.2670269012451)), .Names = c("town",
"long", "lat"), row.names = c("Cooktown", "Cairns", "Townsville",
"Bowen", "Proserpine", "Mackay", "Rockhampton", "Gladstone",
"Port Douglas", "Ayr", "Innisfail", "Cardwell"), class = "data.frame", col = 1, pch = 21, bg = "red", cex = 1, txt.text = c("Cooktown",
"Cairns", "Townsville", "Bowen", "Proserpine", "Mackay", "Rockhampton",
"Gladstone", "Port Douglas", "Ayr", "Innisfail", "Cardwell"), txt.col = 1, txt.cex = 1, txt.pos = 2)



my_wtd_q = function(x, w, prob, n = 4096)
  with(density(x, weights = w/sum(w), n = n),
       x[which.max(cumsum(y*(x[2L] - x[1L])) >= prob)])

cellMeansRaw <- function(dat) {
    dat %>%
        group_by(Year) %>%
        mutate(W2=Tows/sum(Tows)) %>%
        summarise(Mean=mean(Cover,na.rm=TRUE),
                  Median=median(Cover, na.rm=TRUE),
                  Mean.w = weighted.mean(Cover,W2, na.rm=TRUE),
                  Median.w = my_wtd_q(Cover, W2, 0.5)) %>%
        ungroup
}


ggproto_Raw <- function(dat) {
    list(
        geom_line(data=dat, aes(y=Mean, color='Mean')),
        geom_line(data=dat, aes(y=Median, color='Median')),
        geom_line(data=dat,aes(y=Mean.w, color='Mean.w')),
        geom_line(data=dat,aes(y=Median.w, color='Median.w'))
    )
}



ModelOriginal <- function(form, dat, location) {
    print(form)
    mod <- fitOriginal(form, dat)
    newdata <- cellMeansOriginal(mod, dat, location)
    return(list(mod=mod, newdata=newdata))
}

fitOriginal <- function(form, dat) {
    mod <-  stan_glmer(form,
                       data=dat,
                       family=binomial,
                       iter=5000,
                       warmup=2500,
                       chains=3,cores=3)
    return(mod)
}

cellMeansOriginal <- function(mod, dat, location) {
    newdata = data.frame(Location=location,
                         Year=unique(dat$Year),
                         N=length(unique(dat$REEF_NAME)))
    Xmat = model.matrix(~Year, newdata)
    coefs = data.frame(mod) %>% dplyr:::select(matches('^X.Intercept*|^Year.*'))
    Fit=binomial()$linkinv(as.matrix(coefs) %*% t(Xmat))
    newdata = cbind(newdata,
                    plyr:::adply(Fit,2,function(x) {
                        data.frame(mean=mean(x), coda:::HPDinterval(coda:::as.mcmc(x)))
                    })
                    )
    return(newdata)
}
