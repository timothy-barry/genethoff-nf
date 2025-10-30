library(openxlsx)
library(readr)

dir <- paste0(.get_config_path("LOCAL_BAUER_LAB_DATA_DIR"), "guideseq_genethoff/fastq/")
sample_ids <- c("293T_SpRY_Cas9_1620_GSPneg", "293T_SpRY_Cas9_1620_GSPplus", "293T_SpRY_Cas9_55SG3_GSPneg", "293T_SpRY_Cas9_55SG3_GSPplus")
sample_names_prefixes <- c("293T_SpRY_Cas9_1620_GSPneg_S79_S87_L001_", "293T_SpRY_Cas9_1620_GSPplus_S76_S84_L001_",
                           "293T_SpRY_Cas9_55SG3_GSPneg_S80_S88_L001_", "293T_SpRY_Cas9_55SG3_GSPplus_S77_S85_L001_")
sample_suffix <- "_001.fastq.gz"

dt <- lapply(X = seq_along(sample_names_prefixes), FUN = function(i) {
  file_names <- paste0(dir, paste0(sample_names_prefixes[i], c("I1", "I2", "R1", "R2"), sample_suffix))
  m <- t(matrix(c(sample_ids[i], file_names)))
  colnames(m) <- c("sample_id", "I1", "I2", "R1", "R2")
  as.data.frame(m)
}) |> data.table::rbindlist()
write_tsv(x = dt, file = "datasheet.tsv")
