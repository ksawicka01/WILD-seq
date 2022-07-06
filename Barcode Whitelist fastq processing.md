# Processing fastq sequencing files from WILD-seq barcode RT-PCR to generate whitelist of barcodes
This pipeline generates a whitelist of barcodes present in the WILD-seq cell pool. This whitelist is used as a reference for assignment of WILD-seq barcodes in subsequent single cell experiments

## Prerequisits
**cutadapt** Available from https://cutadapt.readthedocs.io/en/stable/#

**umi_tools** Available from https://umi-tools.readthedocs.io/en/latest/

## Input file
RT-PCR is performed as described in Materials and Methods subsection 'Whitelist generation of WILD-seq barcodes' and the fastq files from the sequencing of the PCR products processed using the following pipeline to generate a table of barcodes and associated UMI counts

The following commands are based on starting from a fastq file called SLX-21864.DNAA001.000000000-KF8JL.s_1.r_1.fq.gz which is available in the Data folder.

## Trim 3'end of reads 
Reads are trimmed on the 3'end to remove flanking sequence after the barcode
```
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --discard-untrimmed -o SLX-21864.DNAA001.trimmed.fq.gz  SLX-21864.DNAA001.000000000-KF8JL.s_1.r_1.fq.gz
```

## Filter reads
Reads are filtered for those that contain the barcode sequence including the common linker sequence and the UMI from the RT primer
```
zcat SLX-21864.DNAA001.trimmed.fq.gz | egrep -A 2 -B 1 ^[ATGC]{12}TGCATCGGTTAACCGATGCA[ATGC]{12}CGGATAGAACTTTGAATCGCTTG[ATGC]{8}$ | sed '/^--$/d' > SLX-21864.DNAA001.trimmed.filtered.fq
```

## Extract UMI sequence
```
umi_tools extract --extract-method=string --3prime --bc-pattern=NNNNNNNN --stdin SLX-21864.DNAA001.trimmed.filtered.fq --log=processed.log --stdout SLX-21864.DNAA001.trimmed.filtered.UMI.fq
```

## Trim 5' end of reads
Reads are trimmed on the 5'end to remove flanking sequence after the barcode
```
cutadapt -a CGGATAGAACTTTGAATCGCTTG --discard-untrimmed -o SLX-21864.DNAA001.BC.UMI.fq SLX-21864.DNAA001.trimmed.filtered.UMI.fq
```

## Parse file for downstream analysis
This prepares a table with every detected barcode-UMI combination and the number of reads with this combination. This file is then used for further analysis in R.
```
awk 'BEGIN{RS="@M"}{print $1,"\t",$3}' SLX-21864.DNAA001.BC.UMI.fq | awk -F"\t|_" '{print $2,$3}' | awk -e '$1 ~ /\w/ {print $0}'> SLX-21864.DNAA001.BC.UMI.txt
sort SLX-21864.DNAA001.barcodes.UMI.txt | uniq -c > SLX-21864.DNAA001.barcodes.UMI.counts.txt
```

## Generate barcode whitelist
Use SLX-21864.DNAA001.barcodes.UMI.counts.txt as input for Create_Barcode_Whitelist.R

## Generate bowtie index for barcode whitelist
Use the barcode whitelist to generate a bowtie index against which barcode sequences from single cellexperiments will be mapped 

