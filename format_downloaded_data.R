library(GEOquery)
library(data.table)
library(readxl)
library(stringr)
library(illuminaHumanv4.db)

# Look at the first rows
head(data)

#args <- commandArgs(trailingOnly = TRUE)
#work_dir <- args[1]
work_dir <- "C:/ICB_Kwek/data"

# CLIN.txt
gse <- getGEO("GSE39688", GSEMatrix = TRUE)
clin <- pData(gse$GSE39688_series_matrix.txt.gz)

write.table(clin, file=file.path(work_dir, 'CLIN.txt'), sep = "\t" , quote = FALSE , row.names = FALSE)

# EXP_TPM.tsv
untar(file.path(work_dir, "GSE39688_RAW.tar"),
      exdir = file.path(work_dir, "GSE39688"))

gset <- getGEO("GSE39688", GSEMatrix =TRUE)
expr <- exprs(gset$GSE39688_series_matrix.txt.gz)

# convert probe ID to gene symbol
prob.gene <- select(illuminaHumanv4.db, 
             keys = rownames(expr), 
             columns=c("SYMBOL"), 
             keytype="PROBEID")
  
prob.gene <- prob.gene[!is.na(prob.gene$SYMBOL), ]
prob.gene <- prob.gene[!duplicated(prob.gene$SYMBOL), ]
prob.gene <- prob.gene[!duplicated(prob.gene$PROBEID), ]
expr <- expr[rownames(expr) %in% prob.gene$PROBEID, ]
prob.gene <- prob.gene[order(prob.gene$PROBEID), ]
expr <- expr[order(row.names(expr)), ]
rownames(expr) <- prob.gene$SYMBOL

load("C:/ICB_Cloughesy/data/output/Gencode.v19.annotation.RData")
expr <- expr[rownames(expr) %in% features_gene$gene_name, ]

write.table(expr, file=file.path(work_dir, 'EXP_TPM.tsv'), sep = "\t" , quote = FALSE , row.names = TRUE, col.names=TRUE)

