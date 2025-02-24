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

| **Unix/Linux** | kallisto, Bowtie2, SPAdes, samtools, BLAST+ |

| **Python**    | biopython, pandas, os, subprocess | 

| **R** | sleuth, dplyr  |

## **2. Accessing Data**
