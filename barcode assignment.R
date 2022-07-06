# This R script provides an example of how WILD-seq barcodes were assigned to a single cell sequencing run

library(dplyr)
library(tidyr)
library(Seurat)

## This section processes the barcode reads that were extracted from the PCR enrichment sequencing run
DV1_enrich <- read.table("SITTB1_EnrichPCR_table_full.txt", header=T)

aggregate <- DV1_enrich %>% group_by(CB,barcode_name) %>% summarize(UMI_count=length(UMI))
DV1_CB_counts <- aggregate %>% group_by(CB) %>% summarize(total_UMI=sum(UMI_count))
DV1_enrich <- merge(aggregate, DV1_CB_counts, by="CB")
DV1_enrich$fraction <- DV1_enrich$UMI_count/DV1_enrich$total_UMI

DV1_enrich.top <- DV1_enrich %>% group_by(CB) %>% slice_max(n = 1, order_by = fraction)
DV1_enrich.second <- DV1_enrich %>% group_by(CB) %>% arrange(desc(fraction)) %>% slice(2)

DV1_enrich <- merge(DV1_enrich.top, DV1_enrich.second, all = T, by = "CB")
DV1_enrich$ratio <- DV1_enrich$fraction.x / DV1_enrich$fraction.y

DV1_filtered <- subset(DV1_enrich, is.na(DV1_enrich$ratio) == T | DV1_enrich$ratio >= 2) # at least twice as many UMIs for top asssignment compared to the next best
DV1_filtered <- subset(DV1_filtered, DV1_filtered$UMI_count.x >= 10) #At least 10 UMIs supporting a WILD-seq barcode assignment are required 
DV1_filtered <- DV1_filtered[,c("CB", "barcode_name.x", "UMI_count.x", "total_UMI.x", "fraction.x", "fraction.y", "ratio")]
colnames(DV1_filtered) <- c("CBC", "barcode_name", "UMI_count", "total_UMI", "fraction.first", "fraction.second", "ratio")


## This section processes the barcode reads that were extracted from the single cell transcriptomic data
DV1.10x <- read.table("SITTA11_CBC_table_full.txt", header = T)
DV1.10x$CBC <- gsub("CB:Z:", "", DV1.10x$CBC)

aggregate <- DV1.10x %>% group_by_(.dots=c("CBC","barcode_name")) %>% summarize(UMI_count=length(UMI))
CB_counts <- aggregate %>% group_by_("CBC") %>% summarize(total_UMI=sum(UMI_count))
DV1.10x.table <- merge(aggregate, CB_counts, by="CBC")
DV1.10x.table$fraction <- DV1.10x.table$UMI_count/DV1.10x.table$total_UMI
DV1.10x.table <- subset(DV1.10x.table, DV1.10x.table$UMI_count >=2 & DV1.10x.table$fraction > 0.5)# at least 2 UMIs and more that 50% of total UMIs for the cell supporting assignment


## Barcode assignment from the single cell transcriptomic sequencing data and from the PCR enrichment are combined to get final assignment
DV1.merge <- merge(DV1.10x.table, DV1_filtered, by = "CBC", all = T)

DV1.merge$clonex_BC <- ifelse(is.na(DV1.merge$barcode_name.x) == TRUE, DV1.merge$barcode_name.y, ifelse(is.na(DV1.merge$barcode_name.y) == T, DV1.merge$barcode_name.x, ifelse(DV1.merge$barcode_name.x == DV1.merge$barcode_name.y, DV1.merge$barcode_name.x, "non-match")))
DV1.merge$origin <- ifelse(is.na(DV1.merge$barcode_name.x) == TRUE, "pcr-enrichment", ifelse(is.na(DV1.merge$barcode_name.y) == T, "10X", ifelse(DV1.merge$barcode_name.x == DV1.merge$barcode_name.y, "both", "non-match")))
DV1.merge <- DV1.merge[!(DV1.merge$origin == "10X" & DV1.merge$UMI_count.x < 5),] # more stringent filtering if not represented in both PCR enrichment and 10X data
DV1.merge <- DV1.merge[!(DV1.merge$origin == "pcr-enrichment" & DV1.merge$UMI_count.y < 30),] # more stringent filtering if not represented in both PCR enrichment and 10X data

barcode.table <- DV1.merge[,c("CBC","clonex_BC", "origin")]

write.table(barcode.table, "assigned_barcodes_DV1.txt", sep = "\t", row.names = F, quote = F)
