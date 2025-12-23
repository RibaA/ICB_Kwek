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
expr  <- exprs(eset)
annot <- fData(eset)

stopifnot(nrow(expr) == nrow(annot))  # should be TRUE

# Clean accessions
refseq <- sub("\\..*$", "", annot$GB_ACC)
refseq <- trimws(refseq)
refseq[refseq == ""] <- NA

# Only map RefSeq-like IDs; keep full length vector
is_refseq <- !is.na(refseq) & grepl("^(NM|NR|XM|XR)_", refseq)

# Map unique keys once
keys_use <- unique(refseq[is_refseq])

symbol_by_refseq <- AnnotationDbi::mapIds(
  x = org.Hs.eg.db,
  keys = keys_use,
  keytype = "REFSEQ",
  column = "SYMBOL",
  multiVals = "first"
)

# Build row-aligned symbol vector (length == nrow(expr))
symbol_map <- rep(NA_character_, nrow(expr))
symbol_map[is_refseq] <- unname(symbol_by_refseq[refseq[is_refseq]])

length(symbol_map) == nrow(expr)  # should be TRUE
rownames(expr) <- symbol_map
expr <- expr[, order(colnames(expr))]

write.table(expr, file=file.path(work_dir, 'EXP_TPM.tsv'), sep = "\t" , quote = FALSE , row.names = TRUE, col.names=TRUE)

