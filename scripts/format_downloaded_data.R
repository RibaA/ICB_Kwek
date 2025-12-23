library(GEOquery)
library(data.table)
library(readxl)
library(stringr)
library(AnnotationDbi)
library(org.Hs.eg.db)

#args <- commandArgs(trailingOnly = TRUE)
#work_dir <- args[1]
work_dir <- "data"

# CLIN.txt
gse <- getGEO("GSE39688", GSEMatrix = TRUE)
clin <- pData(gse$GSE39688_series_matrix.txt.gz)
clin$patientid <- clin$geo_accession
clin <- clin[order(clin$patientid), ]

write.table(clin, file=file.path(work_dir, 'CLIN.txt'), sep = "\t" , quote = FALSE , row.names = FALSE)

# EXP_TPM.tsv
gset <- getGEO("GSE39688", GSEMatrix =TRUE)
expr  <- exprs(gset$GSE39688_series_matrix.txt.gz)
annot <- fData(gset$GSE39688_series_matrix.txt.gz)

stopifnot(nrow(expr) == nrow(annot))  # should be TRUE

acc <- sub("\\..*$", "", trimws(annot$GB_ACC))
acc[acc == ""] <- NA

# 1) Try REFSEQ for NM_/NR_/XM_/XR_
is_refseq <- !is.na(acc) & grepl("^(NM|NR|XM|XR)_", acc)
keys_ref  <- unique(acc[is_refseq])

sym_ref <- mapIds(
  org.Hs.eg.db,
  keys = keys_ref,
  keytype = "REFSEQ",
  column = "SYMBOL",
  multiVals = "first"
)

symbol_map <- rep(NA_character_, length(acc))
symbol_map[is_refseq] <- unname(sym_ref[acc[is_refseq]])

# 2) Fill remaining using ACCNUM (covers many BC*/AK*/etc)
need <- is.na(symbol_map) & !is.na(acc)
keys_acc <- unique(acc[need])

sym_acc <- mapIds(
  org.Hs.eg.db,
  keys = keys_acc,
  keytype = "ACCNUM",
  column = "SYMBOL",
  multiVals = "first"
)

symbol_map[need] <- unname(sym_acc[acc[need]])

# Check how many still NA
rownames(expr) <- symbol_map
expr <- expr[!is.na(rownames(expr)), ]
rownames(expr) <- make.unique(rownames(expr))
expr <- expr[, order(colnames(expr))]

write.table(expr, file=file.path(work_dir, 'EXP_TPM.tsv'), sep = "\t" , quote = FALSE , row.names = TRUE, col.names=TRUE)

