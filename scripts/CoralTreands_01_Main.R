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

##################################################################################################################################################
## Extract and read in the Manta-tow data                                                                                                       ##
## and store it under data/primary/manta.csv and data/primary/manta.RData                                                                       ##
## The SQL is as follows                                                                                                                        ##
## SELECT V_RM_SAMPLE.A_SECTOR, V_RM_SAMPLE.SHELF, V_RM_SAMPLE.REEF_NAME,V_RM_SAMPLE.REEF_ID, V_RM_SAMPLE.REEF_LONG, V_RM_SAMPLE.REEF_LAT,      ##
##  V_RM_SAMPLE.REPORT_YEAR, RM_MANTA.TOW_SEQ_NO, RM_MANTA.LIVE_CORAL, V_RM_SAMPLE.SAMPLE_CLASS                                                 ##
## FROM RM_MANTA INNER JOIN V_RM_SAMPLE ON RM_MANTA.SAMPLE_ID = V_RM_SAMPLE.SAMPLE_ID                                                           ##
## WHERE (((V_RM_SAMPLE.SAMPLE_CLASS) In ('K','C','G','Z') Or (V_RM_SAMPLE.SAMPLE_CLASS) Is Null))                                              ##
## ORDER BY V_RM_SAMPLE.A_SECTOR, V_RM_SAMPLE.SHELF, V_RM_SAMPLE.REEF_NAME, V_RM_SAMPLE.REEF_ID, V_RM_SAMPLE.REEF_LAT, V_RM_SAMPLE.REPORT_YEAR, ##
## RM_MANTA.TOW_SEQ_NO","data/manta.sql")                                                                                                       ##
##################################################################################################################################################
source('CoralTrends_30_getData_Manta.R')
source('CoralTrends_31_getData_Manta_Habitat.R')
## source('CoralTrends_32_getData_Manta_Composition.R')

#################################################
## Process the Manta-tow data                  ##
## - Only include data collected after 1985    ##
## - Convert data into percent cover           ##
## - Summarize the data to the reef/year level ##
## - Assign a zone                             ##
## Data are stored in:                         ##
## - data/processed/manta.sum.RData            ##
## - data/processed/manta.sum.newzones.RData   ##
#################################################
source('CoralTrends_40_processData_Manta_3Zone.R')

###################################
## Generate site maps            ##
## Maps are stored in:           ##
## - output/figures/MapOfSites.* ##
###################################
source('CoralTrends_50_spatialMap_3Zone.R')

######################################################
## Generate temporal head maps of data availability ##
## Heat maps are stored in:                         ##
## - output/figures/TemporalHeatMap_3Zone.*        ##
######################################################
source('CoralTrends_51_temporalHeatMap_3Zone.R')

###############################
## Fit STAN (or INLA) models ##
###############################
if (rerun_models) source('CoralTrends_60_stan_Manta_3Zone.R')

##################################
## Disturbance figures and maps ##
##################################
source('CoralTrends_71_bleaching.R')
source('CoralTrends_72_cots.R')
source('CoralTrends_73_cyclones.R')

############################
## Construct trend graphs ##
############################
source('CoralTrends_80_trend_manta_3Zone.R')

## Map
source('CoralTrends_map_manta.R') 

## Disturbance Frequency
source('CoralTrends_75_disturbanceFrequency_3Zone.R')


## Maps of all the disturbances
source('CoralTrend_spatialFigures.R') 


## Zip up some stuff for Mike
zip('output/NewFigures4Mike.zip', files=c('output/figures/3Zones.pdf',
                                          'output/figures/Disturbances_severe_compilation_new.pdf',
                                          'output/figures/Disturbances_all_no_label.pdf',
                                          'output/figures/Fig1mapE_3Zones_new.pdf',
                                          'output/figures/CoralChangeFig.pdf',
                                          'data/modelled/dist.table.intercepts.csv',
                                          'data/modelled/dist.table.slopes.csv'))



