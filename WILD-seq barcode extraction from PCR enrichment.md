# Assignment of WILD-seq barcode extraction from PCR enrichment sequencing data
This pipeline extracts the WILD-seq barcodes and corresponding 10X cell barcode and UMI information from the PCR enrichment

## Prerequisits

**umi_tools** https://umi-tools.readthedocs.io/en/latest/

**cutadapt** https://cutadapt.readthedocs.io/en/stable/#

**bowtie** http://bowtie-bio.sourceforge.net/manual.shtml

## Generate a whitelist of cell barcodes 
A whitelist of cell barcodes is generated based on corresponing 10X transcriptome sequencing run
```
