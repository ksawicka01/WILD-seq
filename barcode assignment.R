##This R script provides an example of how WILD-seq barcodes were assigned to the single cell sequencing run for the D2A1 docetaxel experiment described in the paper

library(dplyr)
library(tidyr)
library(Seurat)

#This section processes the barcode reads that were extracted from the PCR enrichment sequencing run
DV1_enrich <- read.table("SITTB1_EnrichPCR_table_full.txt", header=T)
DV2_enrich <- read.table("SITTB2_EnrichPCR_table_full.txt", header=T)
DV3_enrich <- read.table("SITTB3_EnrichPCR_table_full.txt", header=T)
D1_enrich <- read.table("SITTB4_EnrichPCR_table_full.txt", header=T)
D2_enrich <- read.table("SITTB5_EnrichPCR_table_full.txt", header=T)
D3_enrich <- read.table("SITTB6_EnrichPCR_table_full.txt", header=T)

aggregate <- DV1_enrich %>% group_by(CB,barcode_name) %>% summarize(UMI_count=length(UMI))
DV1_CB_counts <- aggregate %>% group_by(CB) %>% summarize(total_UMI=sum(UMI_count))
DV1_enrich <- merge(aggregate, DV1_CB_counts, by="CB")
DV1_enrich$fraction <- DV1_enrich$UMI_count/DV1_enrich$total_UMI

aggregate <- DV2_enrich %>% group_by(CB,barcode_name) %>% summarize(UMI_count=length(UMI))
DV2_CB_counts <- aggregate %>% group_by(CB) %>% summarize(total_UMI=sum(UMI_count))
DV2_enrich <- merge(aggregate, DV2_CB_counts, by="CB")
DV2_enrich$fraction <- DV2_enrich$UMI_count/DV2_enrich$total_UMI

aggregate <- DV3_enrich %>% group_by(CB,barcode_name) %>% summarize(UMI_count=length(UMI))
DV3_CB_counts <- aggregate %>% group_by(CB) %>% summarize(total_UMI=sum(UMI_count))
DV3_enrich <- merge(aggregate, DV3_CB_counts, by="CB")
DV3_enrich$fraction <- DV3_enrich$UMI_count/DV3_enrich$total_UMI

aggregate <- D1_enrich %>% group_by(CB,barcode_name) %>% summarize(UMI_count=length(UMI))
D1_CB_counts <- aggregate %>% group_by(CB) %>% summarize(total_UMI=sum(UMI_count))
D1_enrich <- merge(aggregate, D1_CB_counts, by="CB")
D1_enrich$fraction <- D1_enrich$UMI_count/D1_enrich$total_UMI

aggregate <- D2_enrich %>% group_by(CB,barcode_name) %>% summarize(UMI_count=length(UMI))
D2_CB_counts <- aggregate %>% group_by(CB) %>% summarize(total_UMI=sum(UMI_count))
D2_enrich <- merge(aggregate, D2_CB_counts, by="CB")
D2_enrich$fraction <- D2_enrich$UMI_count/D2_enrich$total_UMI

aggregate <- D3_enrich %>% group_by(CB,barcode_name) %>% summarize(UMI_count=length(UMI))
D3_CB_counts <- aggregate %>% group_by(CB) %>% summarize(total_UMI=sum(UMI_count))
D3_enrich <- merge(aggregate, D3_CB_counts, by="CB")
D3_enrich$fraction <- D3_enrich$UMI_count/D3_enrich$total_UMI


DV1_enrich.top <- DV1_enrich %>% group_by(CB) %>% slice_max(n = 1, order_by = fraction)
DV1_enrich.second <- DV1_enrich %>% group_by(CB) %>% arrange(desc(fraction)) %>% slice(2)

DV1_enrich <- merge(DV1_enrich.top, DV1_enrich.second, all = T, by = "CB")
DV1_enrich$ratio <- DV1_enrich$fraction.x / DV1_enrich$fraction.y

DV1_filtered <- subset(DV1_enrich, is.na(DV1_enrich$ratio) == T | DV1_enrich$ratio >= 2)
DV1_filtered <- subset(DV1_filtered, DV1_filtered$UMI_count.x >= 10)
DV1_filtered <- DV1_filtered[,c("CB", "barcode_name.x", "UMI_count.x", "total_UMI.x", "fraction.x", "fraction.y", "ratio")]
colnames(DV1_filtered) <- c("CBC", "barcode_name", "UMI_count", "total_UMI", "fraction.first", "fraction.second", "ratio")

DV2_enrich.top <- DV2_enrich %>% group_by(CB) %>% slice_max(n = 1, order_by = fraction)
DV2_enrich.second <- DV2_enrich %>% group_by(CB) %>% arrange(desc(fraction)) %>% slice(2)

DV2_enrich <- merge(DV2_enrich.top, DV2_enrich.second, all = T, by = "CB")
DV2_enrich$ratio <- DV2_enrich$fraction.x / DV2_enrich$fraction.y

DV2_filtered <- subset(DV2_enrich, is.na(DV2_enrich$ratio) == T | DV2_enrich$ratio >= 2)
DV2_filtered <- subset(DV2_filtered, DV2_filtered$UMI_count.x >= 10)
DV2_filtered <- DV2_filtered[,c("CB", "barcode_name.x", "UMI_count.x", "total_UMI.x", "fraction.x", "fraction.y", "ratio")]
colnames(DV2_filtered) <- c("CBC", "barcode_name", "UMI_count", "total_UMI", "fraction.first", "fraction.second", "ratio")

DV3_enrich.top <- DV3_enrich %>% group_by(CB) %>% slice_max(n = 1, order_by = fraction)
DV3_enrich.second <- DV3_enrich %>% group_by(CB) %>% arrange(desc(fraction)) %>% slice(2)

DV3_enrich <- merge(DV3_enrich.top, DV3_enrich.second, all = T, by = "CB")
DV3_enrich$ratio <- DV3_enrich$fraction.x / DV3_enrich$fraction.y

DV3_filtered <- subset(DV3_enrich, is.na(DV3_enrich$ratio) == T | DV3_enrich$ratio >= 2)
DV3_filtered <- subset(DV3_filtered, DV3_filtered$UMI_count.x >= 10)
DV3_filtered <- DV3_filtered[,c("CB", "barcode_name.x", "UMI_count.x", "total_UMI.x", "fraction.x", "fraction.y", "ratio")]
colnames(DV3_filtered) <- c("CBC", "barcode_name", "UMI_count", "total_UMI", "fraction.first", "fraction.second", "ratio")


D1_enrich.top <- D1_enrich %>% group_by(CB) %>% slice_max(n = 1, order_by = fraction)
D1_enrich.second <- D1_enrich %>% group_by(CB) %>% arrange(desc(fraction)) %>% slice(2)

D1_enrich <- merge(D1_enrich.top, D1_enrich.second, all = T, by = "CB")
D1_enrich$ratio <- D1_enrich$fraction.x / D1_enrich$fraction.y

D1_filtered <- subset(D1_enrich, is.na(D1_enrich$ratio) == T | D1_enrich$ratio >= 2)
D1_filtered <- subset(D1_filtered, D1_filtered$UMI_count.x >= 5)
D1_filtered <- D1_filtered[,c("CB", "barcode_name.x", "UMI_count.x", "total_UMI.x", "fraction.x", "fraction.y", "ratio")]
colnames(D1_filtered) <- c("CBC", "barcode_name", "UMI_count", "total_UMI", "fraction.first", "fraction.second", "ratio")

D2_enrich.top <- D2_enrich %>% group_by(CB) %>% slice_max(n = 1, order_by = fraction)
D2_enrich.second <- D2_enrich %>% group_by(CB) %>% arrange(desc(fraction)) %>% slice(2)

D2_enrich <- merge(D2_enrich.top, D2_enrich.second, all = T, by = "CB")
D2_enrich$ratio <- D2_enrich$fraction.x / D2_enrich$fraction.y

D2_filtered <- subset(D2_enrich, is.na(D2_enrich$ratio) == T | D2_enrich$ratio >= 2)
D2_filtered <- subset(D2_filtered, D2_filtered$UMI_count.x >= 5)
D2_filtered <- D2_filtered[,c("CB", "barcode_name.x", "UMI_count.x", "total_UMI.x", "fraction.x", "fraction.y", "ratio")]
colnames(D2_filtered) <- c("CBC", "barcode_name", "UMI_count", "total_UMI", "fraction.first", "fraction.second", "ratio")

D3_enrich.top <- D3_enrich %>% group_by(CB) %>% slice_max(n = 1, order_by = fraction)
D3_enrich.second <- D3_enrich %>% group_by(CB) %>% arrange(desc(fraction)) %>% slice(2)

D3_enrich <- merge(D3_enrich.top, D3_enrich.second, all = T, by = "CB")
D3_enrich$ratio <- D3_enrich$fraction.x / D3_enrich$fraction.y

D3_filtered <- subset(D3_enrich, is.na(D3_enrich$ratio) == T | D3_enrich$ratio >= 2)
D3_filtered <- subset(D3_filtered, D3_filtered$UMI_count.x >= 5)
D3_filtered <- D3_filtered[,c("CB", "barcode_name.x", "UMI_count.x", "total_UMI.x", "fraction.x", "fraction.y", "ratio")]
colnames(D3_filtered) <- c("CBC", "barcode_name", "UMI_count", "total_UMI", "fraction.first", "fraction.second", "ratio")

#This section processes the barcode reads that were extracted from the single cell transcriptomic data
DV1.10x <- read.table("SITTA11_CBC_table_full.txt", header = T)
DV2.10x <- read.table("SITTB11_CBC_table_full.txt", header = T)
DV3.10x <- read.table("SITTC11_CBC_table_full.txt", header = T)
D1.10x <- read.table("SITTD11_CBC_table_full.txt", header = T)
D2.10x <- read.table("SITTE11_CBC_table_full.txt", header = T)
D3.10x <- read.table("SITTF11_CBC_table_full.txt", header = T)

DV1.10x$CBC <- gsub("CB:Z:", "", DV1.10x$CBC)
DV2.10x$CBC <- gsub("CB:Z:", "", DV2.10x$CBC)
DV3.10x$CBC <- gsub("CB:Z:", "", DV3.10x$CBC)
D1.10x$CBC <- gsub("CB:Z:", "", D1.10x$CBC)
D2.10x$CBC <- gsub("CB:Z:", "", D2.10x$CBC)
D3.10x$CBC <- gsub("CB:Z:", "", D3.10x$CBC)

aggregate <- DV1.10x %>% group_by_(.dots=c("CBC","barcode_name")) %>% summarize(UMI_count=length(UMI))
CB_counts <- aggregate %>% group_by_("CBC") %>% summarize(total_UMI=sum(UMI_count))
DV1.10x.table <- merge(aggregate, CB_counts, by="CBC")
DV1.10x.table$fraction <- DV1.10x.table$UMI_count/DV1.10x.table$total_UMI
DV1.10x.table <- subset(DV1.10x.table, DV1.10x.table$UMI_count >=2 & DV1.10x.table$fraction > 0.5) #4942

aggregate <- DV2.10x %>% group_by_(.dots=c("CBC","barcode_name")) %>% summarize(UMI_count=length(UMI))
CB_counts <- aggregate %>% group_by_("CBC") %>% summarize(total_UMI=sum(UMI_count))
DV2.10x.table <- merge(aggregate, CB_counts, by="CBC")
DV2.10x.table$fraction <- DV2.10x.table$UMI_count/DV2.10x.table$total_UMI
DV2.10x.table <- subset(DV2.10x.table, DV2.10x.table$UMI_count >=2 & DV2.10x.table$fraction > 0.5) #329

aggregate <- DV3.10x %>% group_by_(.dots=c("CBC","barcode_name")) %>% summarize(UMI_count=length(UMI))
CB_counts <- aggregate %>% group_by_("CBC") %>% summarize(total_UMI=sum(UMI_count))
DV3.10x.table <- merge(aggregate, CB_counts, by="CBC")
DV3.10x.table$fraction <- DV3.10x.table$UMI_count/DV3.10x.table$total_UMI
DV3.10x.table <- subset(DV3.10x.table, DV3.10x.table$UMI_count >=2 & DV3.10x.table$fraction > 0.5) #824

aggregate <- D1.10x %>% group_by_(.dots=c("CBC","barcode_name")) %>% summarize(UMI_count=length(UMI))
CB_counts <- aggregate %>% group_by_("CBC") %>% summarize(total_UMI=sum(UMI_count))
D1.10x.table <- merge(aggregate, CB_counts, by="CBC")
D1.10x.table$fraction <- D1.10x.table$UMI_count/D1.10x.table$total_UMI
D1.10x.table <- subset(D1.10x.table, D1.10x.table$UMI_count >=2 & D1.10x.table$fraction > 0.5) #1284

aggregate <- D2.10x %>% group_by_(.dots=c("CBC","barcode_name")) %>% summarize(UMI_count=length(UMI))
CB_counts <- aggregate %>% group_by_("CBC") %>% summarize(total_UMI=sum(UMI_count))
D2.10x.table <- merge(aggregate, CB_counts, by="CBC")
D2.10x.table$fraction <- D2.10x.table$UMI_count/D2.10x.table$total_UMI
D2.10x.table <- subset(D2.10x.table, D2.10x.table$UMI_count >=2 & D2.10x.table$fraction > 0.5) #1553

aggregate <- D3.10x %>% group_by_(.dots=c("CBC","barcode_name")) %>% summarize(UMI_count=length(UMI))
CB_counts <- aggregate %>% group_by_("CBC") %>% summarize(total_UMI=sum(UMI_count))
D3.10x.table <- merge(aggregate, CB_counts, by="CBC")
D3.10x.table$fraction <- D3.10x.table$UMI_count/D3.10x.table$total_UMI
D3.10x.table <- subset(D3.10x.table, D3.10x.table$UMI_count >=2 & D3.10x.table$fraction > 0.5) #2028


#Barcode assignment from the single cell transcriptomic sequencing data and from the PCR enrichment are combined to get final assignment
DV1.merge <- merge(DV1.10x.table, DV1_filtered, by = "CBC", all = T)
DV2.merge <- merge(DV2.10x.table, DV2_filtered, by = "CBC", all = T)
DV3.merge <- merge(DV3.10x.table, DV3_filtered, by = "CBC", all = T)
D1.merge <- merge(D1.10x.table, D1_filtered, by = "CBC", all = T)
D2.merge <- merge(D2.10x.table, D2_filtered, by = "CBC", all = T)
D3.merge <- merge(D3.10x.table, D3_filtered, by = "CBC", all = T)

DV1.merge$clonex_BC <- ifelse(is.na(DV1.merge$barcode_name.x) == TRUE, DV1.merge$barcode_name.y, ifelse(is.na(DV1.merge$barcode_name.y) == T, DV1.merge$barcode_name.x, ifelse(DV1.merge$barcode_name.x == DV1.merge$barcode_name.y, DV1.merge$barcode_name.x, "non-match")))
DV1.merge$origin <- ifelse(is.na(DV1.merge$barcode_name.x) == TRUE, "pcr-enrichment", ifelse(is.na(DV1.merge$barcode_name.y) == T, "10X", ifelse(DV1.merge$barcode_name.x == DV1.merge$barcode_name.y, "both", "non-match")))

DV2.merge$clonex_BC <- ifelse(is.na(DV2.merge$barcode_name.x) == TRUE, DV2.merge$barcode_name.y, ifelse(is.na(DV2.merge$barcode_name.y) == T, DV2.merge$barcode_name.x, ifelse(DV2.merge$barcode_name.x == DV2.merge$barcode_name.y, DV2.merge$barcode_name.x, "non-match")))
DV2.merge$origin <- ifelse(is.na(DV2.merge$barcode_name.x) == TRUE, "pcr-enrichment", ifelse(is.na(DV2.merge$barcode_name.y) == T, "10X", ifelse(DV2.merge$barcode_name.x == DV2.merge$barcode_name.y, "both", "non-match")))

DV3.merge$clonex_BC <- ifelse(is.na(DV3.merge$barcode_name.x) == TRUE, DV3.merge$barcode_name.y, ifelse(is.na(DV3.merge$barcode_name.y) == T, DV3.merge$barcode_name.x, ifelse(DV3.merge$barcode_name.x == DV3.merge$barcode_name.y, DV3.merge$barcode_name.x, "non-match")))
DV3.merge$origin <- ifelse(is.na(DV3.merge$barcode_name.x) == TRUE, "pcr-enrichment", ifelse(is.na(DV3.merge$barcode_name.y) == T, "10X", ifelse(DV3.merge$barcode_name.x == DV3.merge$barcode_name.y, "both", "non-match")))

D1.merge$clonex_BC <- ifelse(is.na(D1.merge$barcode_name.x) == TRUE, D1.merge$barcode_name.y, ifelse(is.na(D1.merge$barcode_name.y) == T, D1.merge$barcode_name.x, ifelse(D1.merge$barcode_name.x == D1.merge$barcode_name.y, D1.merge$barcode_name.x, "non-match")))
D1.merge$origin <- ifelse(is.na(D1.merge$barcode_name.x) == TRUE, "pcr-enrichment", ifelse(is.na(D1.merge$barcode_name.y) == T, "10X", ifelse(D1.merge$barcode_name.x == D1.merge$barcode_name.y, "both", "non-match")))

D2.merge$clonex_BC <- ifelse(is.na(D2.merge$barcode_name.x) == TRUE, D2.merge$barcode_name.y, ifelse(is.na(D2.merge$barcode_name.y) == T, D2.merge$barcode_name.x, ifelse(D2.merge$barcode_name.x == D2.merge$barcode_name.y, D2.merge$barcode_name.x, "non-match")))
D2.merge$origin <- ifelse(is.na(D2.merge$barcode_name.x) == TRUE, "pcr-enrichment", ifelse(is.na(D2.merge$barcode_name.y) == T, "10X", ifelse(D2.merge$barcode_name.x == D2.merge$barcode_name.y, "both", "non-match")))

D3.merge$clonex_BC <- ifelse(is.na(D3.merge$barcode_name.x) == TRUE, D3.merge$barcode_name.y, ifelse(is.na(D3.merge$barcode_name.y) == T, D3.merge$barcode_name.x, ifelse(D3.merge$barcode_name.x == D3.merge$barcode_name.y, D3.merge$barcode_name.x, "non-match")))
D3.merge$origin <- ifelse(is.na(D3.merge$barcode_name.x) == TRUE, "pcr-enrichment", ifelse(is.na(D3.merge$barcode_name.y) == T, "10X", ifelse(D3.merge$barcode_name.x == D3.merge$barcode_name.y, "both", "non-match")))


DV1.merge <- DV1.merge[!(DV1.merge$origin == "10X" & DV1.merge$UMI_count.x < 5),]
DV2.merge <- DV2.merge[!(DV2.merge$origin == "10X" & DV2.merge$UMI_count.x < 5),]
DV3.merge <- DV3.merge[!(DV3.merge$origin == "10X" & DV3.merge$UMI_count.x < 5),]
D1.merge <- D1.merge[!(D1.merge$origin == "10X" & D1.merge$UMI_count.x < 5),]
D2.merge <- D2.merge[!(D2.merge$origin == "10X" & D2.merge$UMI_count.x < 5),]
D3.merge <- D3.merge[!(D3.merge$origin == "10X" & D3.merge$UMI_count.x < 5),]

DV1.merge <- DV1.merge[!(DV1.merge$origin == "pcr-enrichment" & DV1.merge$UMI_count.y < 30),]
DV2.merge <- DV2.merge[!(DV2.merge$origin == "pcr-enrichment" & DV2.merge$UMI_count.y < 30),]
DV3.merge <- DV3.merge[!(DV3.merge$origin == "pcr-enrichment" & DV3.merge$UMI_count.y < 30),]
D1.merge <- D1.merge[!(D1.merge$origin == "pcr-enrichment" & D1.merge$UMI_count.y < 30),]
D2.merge <- D2.merge[!(D2.merge$origin == "pcr-enrichment" & D2.merge$UMI_count.y < 30),]
D3.merge <- D3.merge[!(D3.merge$origin == "pcr-enrichment" & D3.merge$UMI_count.y < 30),]


DV1.merge$CBC <- paste0(DV1.merge$CBC, "-1_1")
DV2.merge$CBC <- paste0(DV2.merge$CBC, "-1_2")
DV3.merge$CBC <- paste0(DV3.merge$CB, "-1_3")
D1.merge$CBC <- paste0(D1.merge$CBC, "-1_4")
D2.merge$CBC <- paste0(D2.merge$CBC, "-1_5")
D3.merge$CBC <- paste0(D3.merge$CB, "-1_6")

barcode.table <- rbind(DV1.merge, DV2.merge, DV3.merge, D1.merge, D2.merge, D3.merge)

write.table(barcode.table, "assigned_docetaxel_d2a1_clone_barcodes.txt", sep = "\t", row.names = F, quote = F)


##The WILD-seq barcode assignment can then be incorporated into the metadata of the single cell seurat object
seurat <- readRDS("d2a1_DTX_merged.rds")
barcode.table <- read.table("assigned_docetaxel_d2a1_clone_barcodes.txt", header = T)


barcode.table <- barcode.table[,c("CBC","clonex_BC", "origin")]
rownames(barcode.table) <- barcode.table$CBC
barcode.table[,1] <- NULL

seurat <- AddMetaData(object = seurat, metadata = barcode.table, col.name=c("clone", "identification.method"))
seurat@meta.data$clone[is.na(seurat@meta.data$clone)] <- "Unknown"
seurat@meta.data$identification.method[is.na(seurat@meta.data$identification.method)] <- "None"

saveRDS(seurat, "d2a1_DTX_resequenced_with_clonal_assignment.rds")


