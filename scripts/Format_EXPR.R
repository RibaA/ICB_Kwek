library(data.table)
library(R.utils)

args <- commandArgs(trailingOnly = TRUE)
#input_dir <- args[1]
#output_dir <- args[2]

input_dir <- "data/input"
output_dir <- "data/output"

tpm = read.csv(
  file.path(input_dir, "EXP_TPM.tsv"),
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  sep = '\t'
)

#############################################################################
#############################################################################

case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )
tpm = tpm[ , colnames(tpm) %in% case[ case$expr %in% 1 , ]$patient ]

write.table( tpm , file= file.path(output_dir, "EXPR.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=TRUE )
