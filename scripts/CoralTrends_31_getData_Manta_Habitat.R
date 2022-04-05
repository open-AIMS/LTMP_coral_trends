library(tidyverse)

## SQL for extracting manta tow data from oracle
writeLines("
SELECT V_RM_SAMPLE.P_CODE, V_RM_SAMPLE.A_SECTOR, V_RM_SAMPLE.SHELF, V_RM_SAMPLE.REEF_NAME, V_RM_SAMPLE.REEF_ID, V_RM_SAMPLE.REEF_LONG, V_RM_SAMPLE.REEF_LAT, V_RM_SAMPLE.REPORT_YEAR, V_MANTA_ZONES.REEF_ZONE_CODE, RM_MANTA.TOW_SEQ_NO, RM_MANTA.LIVE_CORAL
FROM (RM_MANTA INNER JOIN V_RM_SAMPLE ON RM_MANTA.SAMPLE_ID = V_RM_SAMPLE.SAMPLE_ID) INNER JOIN V_MANTA_ZONES ON (V_RM_SAMPLE.SAMPLE_ID = V_MANTA_ZONES.SAMPLE_ID) AND (RM_MANTA.MANTA_SID = V_MANTA_ZONES.MANTA_SID)
ORDER BY V_RM_SAMPLE.A_SECTOR, V_RM_SAMPLE.SHELF, V_RM_SAMPLE.REEF_NAME, V_RM_SAMPLE.REPORT_YEAR, V_MANTA_ZONES.REEF_ZONE_CODE, RM_MANTA.TOW_SEQ_NO", paste0(DATA_PATH, "primary/manta.hab.sql"))

if (goto_database_manta) system(paste0("java -jar dbExport.jar ",
                                       DATA_PATH, "primary/manta.hab.sql ",
                                       DATA_PATH, "primary/manta.hab.csv reef reefmon"))

manta.hab <- read.csv(paste0(DATA_PATH, 'primary/manta.hab.csv'),strip.white=TRUE)

save(manta.hab, file=paste0(DATA_PATH, 'primary/manta.hab.RData'))
