# Transcriptome Analysis Pipeline

## Overview
This pipeline automates the analysis of the Human herpesevirus 5, also known as Human cytomegalovirus (HCMV), transcriptomes at **2 and 6 days post-infection (dpi)**. It retrieves sequencing data, processes reads, and identifies differentially expressed genes. This pipeline can be used for other transcriptomes as well.

## Features
- Downloads transcriptomic data from NCBI SRA
- Builds a transcriptome index for **kallisto**
- Quantifies transcript abundances
- Identifies differentially expressed genes using **sleuth (R)**
- Maps reads to the HCMV genome using **Bowtie2**
- Assembles transcriptomes using **SPAdes**
- Performs BLAST searches to compare sequences

---

## **1. Installation & Dependencies**
Ensure you have the following dependencies installed:

### **Required Software**

| **Unix/Linux** | fasterq-dump, kallisto, Bowtie2, SPAdes, samtools, BLAST+ |

| **Python**    | biopython, pandas, os, subprocess | 

| **R** | sleuth, dplyr  |

## **2. Accessing Data**
Human herpesvirus 5 (HCMV) transcriptome data from two donors, both containing transcriptomes 2 days post infection (2dpi) and 6 days post infection (6dpi). Transcriptome data was found in NCBI, using ```{sh}wget``` to download each SRR run paired end read data. Next, ```{sh}fasterq-dump``` was utilized to uncompress the paired end fastq files and stored in a separate directory to be utilized in the ```{py}wrapper.py``` script.

SRR Data: [SRR5660030 (Donor 1 - 2dpi)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR5660030&display=data-access)

SRR Data: [SRR5660033 (Donor 1 - 6dpi)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR5660033&display=data_access)

SRR Data: [SRR5660044 (Donor 3 - 2dpi)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR5660044&display=data-access)

SRR Data: [SRR5660045 (Donor 3 - 6dpi)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR5660045&display=data-access)
