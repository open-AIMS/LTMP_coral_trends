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

    require(rstanarm)
    require(coda)
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

