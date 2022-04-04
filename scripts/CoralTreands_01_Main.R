#########################################################################
## Read in the project specific functions                              ##
## This also includes running a function (CoralTrends_checkPackages()) ##
## that assesses whether all the required packages are present         ##
#########################################################################
source('CoralTrends_functions.R') 
CoralTrends_checkPackages()

#########################################################################
## Read and load the configurations set in parameters/CoralTrends.conf ##
#########################################################################
source('CoralTrends_10_config.R')

##############################################################
## Generate two alternative zoning configurations           ##
## 1. a three zone configuration based on De'ath et al 2012 ##
##    Zones are stored in:                                  ##
##    - data/spatial/whagbr.RData                           ##
##    - data/spatial/whagbr.n.RData                         ##
##    - data/spatial/whagbr.c.RData                         ##
##    - data/spatial/whagbr.s.RData                         ##
##    - data/spatial/qld.RData                              ##
## 2. a four zone configuration based on GBRMPA regions     ##
##    Zones are stored in:                                  ##
##    - data/spatial/management.RData                       ##
##    - data/spatial/qld.RData                              ##
##############################################################
source('CoralTrends_20_spatial_3Zone.R')
