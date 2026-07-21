#!/usr/bin/env Rscript

# obtain the command line arguments
args <- commandArgs(trailingOnly = TRUE)
paired_counts_fp <- args[1] # "paired_df"
r2_counts_fp <- args[2] # "r2_df"
to_save <- args[3] # "293T_SpRY_Cas9_55SG3_GSPneg_count_table.rds"

col_names <- c("chr", "start", "stop", "primer", "umi", "strand", "read_count", "mean_mapq")
process_df <- function(df) df |> dplyr::select(chr, start, umi, strand, read_count, mean_mapq)
paired_counts_df <- read.delim(paired_counts_fp, header = FALSE, col.names = col_names) |> process_df()
r2_counts_df <- read.delim(r2_counts_fp, header = FALSE, col.names = col_names) |> process_df()

df <- rbind(paired_counts_df, r2_counts_df) |>
  dplyr::distinct() |>
  dplyr::group_by(chr, start, strand) |>
  dplyr::summarize(umi_count = dplyr::n(),
                   n_reads_per_umi = I(list(read_count)), # paste0(read_count, collapse = "-") ,
                   total_read_count = sum(read_count),
                   mean_mapq = sum(read_count/total_read_count * mean_mapq)) |>
  dplyr::ungroup()
saveRDS(object = df, file = to_save)