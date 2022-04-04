########################################################################
## This module loads a cofiguration file that contains paramaters for ##
## settings to be applied across all aspects of the analyses          ##
## This file is a text file with key=value pairs                      ##
##                                                                    ##
## Specifically, the pairs are:                                       ##
## Size:                    the number of bootstrapp samples          ##
## Timeseries.plot.width    the width of a timeseries plot            ##
## EndDate                  the last day (date) of the current focal  ##
##                            year (must be in YYYY-MM-DD format)     ##
## Timeseries.plot.height   the height of a timeseries plot           ##
## FocalYear                the current report card year              ##
## StartDate                the lower date range cutoff (must be in   ##
##                            YYY-MM-DD format)                       ##
########################################################################


CoralTrends_tryCatch(
{
  PATH <<- '../'
  DATA_PATH <<- paste0(PATH, 'data/')
  OUTPUT_PATH <<- paste0(PATH, 'output/')

  if(!any(grepl('^data$', list.files(PATH)))) system(paste0('mkdir ', DATA_PATH))
  files <- list.files(DATA_PATH)
  if(!any(grepl('^primary$',files))) system(paste0('mkdir ', DATA_PATH, 'primary'))
  if(!any(grepl('^processed$',files))) system(paste0('mkdir ', DATA_PATH, 'processed'))
  if(!any(grepl('^spatial$',files))) system(paste0('mkdir ', DATA_PATH, 'spatial'))

  files <- list.files(PATH)
  if(!any(grepl('^figures$',files))) system(paste0('mkdir ', PATH, 'figures'))

  files <- list.files(PATH)
  if(!any(grepl('^output$',files))) system(paste0('mkdir ', OUTPUT_PATH))
  if(!any(grepl('^output$',files))) system(paste0('mkdir ', OUTPUT_PATH, 'figures'))

  files <- list.files(PATH)
  if(!any(grepl('^logs$',files))) system(paste0('mkdir ', PATH, 'logs'))
  return=NULL
}, '../logs/all.log','--Config--',msg='configure necessary folders', return=NULL)

CoralTrends_tryCatch(
{
  config = readLines(paste0(PATH, 'parameters/CoralTrends.conf'))
  config = gsub('(.*Date)=(.*)','\\1=as.Date(\'\\2\')',config)
  eval(parse(text=config))
  return=NULL
}, '../logs/all.log','--Config--',msg='load general configurations', return=NULL)
