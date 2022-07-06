## This R script processes a file containing WILD-seq barcode sequences and their associated UMI counts such as that generated using the pipeline described in Barcode Whitelist fastq processing to generate a whitelist of barcodes present in the sample

library(dplyr)
library(tidyr)
library(stringdist)

data <- read.table("SLX-21864.DNAA001.BC.UMI.counts.txt") #This was the file created using the Barcode Whitelist fastq processing pipeline
colnames(data) <- c("count", "UMI", "barcode")

BC_counts <- data %>% group_by(barcode) %>% summarize(UMI_count=n())
BC_filtered <- subset(BC_counts, UMI_count >= quantile(UMI_count, 0.90))
distance <- stringdistmatrix(BC_filtered$barcode, method="hamming")
REDIST<-as.dist(distance)
sum(is.infinite(REDIST))
hc0 <- hclust(REDIST,method = "complete")
plot(hc0,hang =-1)
ct <- cutree(hc0, h=5)
cl <-rect.hclust(hc0, h=5, border=c(1:length(unique(ct)))+1)

BC_clustered <- cbind(BC_filtered, cluster=0)
row.names(BC_clustered) <- 1:dim(BC_clustered)[1]

for( i in 1:length(cl)){
    BC_clustered[cl[[i]],"cluster"] <- i
  }

BC_whitelist <- BC_clustered %>% group_by(cluster) %>% top_n(1, UMI_count)

write.table(BC_clustered, file="DNAA001_BC_all.txt", sep = "\t", quote = F, row.names = F)
write.table(BC_whitelist, file="DNAA001_BC_whitelist.txt"), sep = "\t", quote = F, row.names = F)

# create a fasta file of barcode sequences to use to generate a bowtie index
output<- character(nrow(DNAA001_BC_whitelist) * 2)
output[c(TRUE, FALSE)] <- paste0(">", DNAA001_BC_whitelist$cluster)
output[c(FALSE, TRUE)] <- as.character(DNAA001_BC_whitelist$barcode)
writeLines(output, "DNAA001_BC_whitelist.fasta")
