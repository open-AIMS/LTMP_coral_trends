library(tidyverse)

## SQL for extracting manta tow data from oracle
writeLines("
SELECT V_RM_SAMPLE.P_CODE, V_RM_SAMPLE.A_SECTOR, V_RM_SAMPLE.SHELF, V_RM_SAMPLE.REEF_NAME,V_RM_SAMPLE.REEF_ID, V_RM_SAMPLE.REEF_LONG, V_RM_SAMPLE.REEF_LAT,
  V_RM_SAMPLE.REPORT_YEAR, RM_MANTA.TOW_SEQ_NO, RM_MANTA.LIVE_CORAL, V_RM_SAMPLE.SAMPLE_CLASS
FROM RM_MANTA INNER JOIN V_RM_SAMPLE ON RM_MANTA.SAMPLE_ID = V_RM_SAMPLE.SAMPLE_ID
WHERE (((V_RM_SAMPLE.SAMPLE_CLASS) In ('K','C','G','Z') Or (V_RM_SAMPLE.SAMPLE_CLASS) Is Null))
ORDER BY V_RM_SAMPLE.P_CODE, V_RM_SAMPLE.A_SECTOR, V_RM_SAMPLE.SHELF, V_RM_SAMPLE.REEF_NAME, V_RM_SAMPLE.REEF_ID, V_RM_SAMPLE.REEF_LAT, V_RM_SAMPLE.REPORT_YEAR,
RM_MANTA.TOW_SEQ_NO", paste0(DATA_PATH, "primary/manta.sql"))

## if (goto_database_manta) system(paste0("java -jar ", DATA_PATH, "scripts/dbExport.jar ",
if (goto_database_manta) system(paste0("java -jar dbExport.jar ",
                                       DATA_PATH, "primary/manta.sql ",
                                       DATA_PATH, "primary/manta.csv reef reefmon"))

manta <- read.csv(paste0(DATA_PATH, "primary/manta.csv"),strip.white=TRUE)

save(manta, file=paste0(DATA_PATH, 'primary/manta.RData'))
