# Assignment of WILD-seq barcode extraction from PCR enrichment sequencing data
This pipeline extracts the WILD-seq barcodes and corresponding 10X cell barcode and UMI information from the PCR enrichment

## Prerequisits

**umi_tools** https://umi-tools.readthedocs.io/en/latest/

**cutadapt** https://cutadapt.readthedocs.io/en/stable/#

**bowtie** http://bowtie-bio.sourceforge.net/manual.shtml

## Generate a whitelist of cell barcodes 
A whitelist of cell barcodes is generated based on corresponing 10X transcriptome sequencing run
```
# If required combine together all Read1 fastq files from 10X transcriptomics experiment
cat SITTA11_*_R1*.fastq.gz > SITTA11_R1_mRNA_All.fastq.gz

gunzip SITTA11_R1_mRNA_All.fastq.gz
umi_tools whitelist --stdin SITTA11_R1_mRNA_All.fastq --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN --set-cell-number=6000 --log2stderr > whitelist_SIGAD10_mRNA_All.txt
```

## Extract CBs and UMIs 
Cell barcode and UMI information is extracted from Read1 of PCR run and added into the Read2 names
```
# If required combine all read1 fastq files and all read2 fastq files into a single file each
cat SITTB1*_R1*.fastq > SITTB1_PCR_All_R1.fastq
cat SITTB1*_R2*.fastq > SITTB1_PCR_All_R2.fastq

# Filter cell barcodes in the PCR sequencing for those on the transcriptomics run whitelist
umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN --stdin SITTB1_PCR_All_R1.fastq --read2-stdout --read2-in=SITTB1_PCR_All_R2.fastq --stdout=SITTB1_PCR_All_R2_extracted.fastq --whitelist=whitelist_SIGAD10_mRNA_All.txt --filter-cell-barcode
```

## Extract barcode sequence 
Leading and trailing sequences are removed and only readss which match the expected WILD-seq barcode pattern are retained.
```
cutadapt -g NNNNNNNNNNCAGCCATGCGCTCGTTTACTATACGAT --discard-untrimmed -o SITTB1_PCR_All_R2_extracted_trim5.fastq SITTB1_PCR_All_R2_extracted.fastq
cutadapt -a CGGATAGAACT --discard-untrimmed -o SITTB1_PCR_All_R2_extracted_trim5_trim3.fastq SITTB1_PCR_All_R2_extracted_trim5.fastq
egrep -A 2 -B 1 [ATGC]{12}TGCATCGGTTAACCGATGCA[ATGC]{12} SITTB1_PCR_All_R2_extracted_trim5_trim3.fastq | sed '/^--$/d' > SITTB1_PCR_All_R2_extracted_trim5_trim3_filtered.fastq
```

## Map the WILD-seq barcodes
The WILD-seq barcode sequences are mapped to the previously generated barcode whitelist bowtie index
```
bowtie -v 2 BC_index SITTB1_PCR_All_R2_extracted_trim5_trim3_filtered.fastq > SITTB1_PCR_All.bowtie
```

## Summarise cell barcode, WILD-seq barcode and UMI information
The bowtie file is parsed using the following R script to generate the required table for downstream analysis
```
library(dplyr)
library(tidyr)

bowtie <- read.table("SITTB1_PCR_All.bowtie", sep="\t", header=F)
colnames(bowtie) <- c("read","strand","barcode_name","offset","seq","qualities","X", "mismatches")

bowtie$read <- as.character(bowtie$read)
bowtie <- bowtie %>% separate(read, into=c("read", "index"), sep=" ")
bowtie <- bowtie %>% separate(read, into=c("read", "CB", "UMI"), sep="_")

table <- bowtie[,c("CB", "UMI", "barcode_name")]

table <- distinct(table)
write.table(table, "SITTB1_PCR_All.bowtie_EnrichPCR_table_full.txt", sep="\t", row.names = F, quote = F)
```

