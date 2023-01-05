# WILD-seq : Wholistic Interrogation of Lineage Dynamics by sequencing for clonal transcriptomic analysis and identification of chemoresistance.

This repository contains code related to the publication:
Sophia A Wild Ian G Cannell Ashley Nicholls Katarzyna Kania Dario Bressan CRUK IMAXT Grand Challenge Team Gregory J Hannon Kirsty Sawicka (2022) Clonal transcriptomics identifies mechanisms of chemoresistance and empowers rational design of combination therapies eLife 11:e80981.
https://doi.org/10.7554/eLife.80981

Details for analysis of the following steps in the protocol are provided:

1. Generating a whitelist of WILD-seq barcodes present in a cell pool from the RT-PCR sequencing run
- Process the fastq files to extract the barcode sequences
- Create the barcode whitelist

2. Assign WILD-seq barcodes to cells from a 10X scRNA-seq experiment
- extracting WILD-seq barcodes from the 10X sequencing file
- extracting WILD-seq barcodes from the PCR enrichment sequencing file
- assigning WILD-seq barcodes to single cells 
