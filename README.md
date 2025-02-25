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

| **Unix/Linux** -> fasterq-dump, kallisto, Bowtie2, SPAdes, samtools, BLAST+ |

| **Python** -> biopython, pandas, os, subprocess | 

| **R** -> sleuth, dplyr  |

## **2. Accessing Data**
Human herpesvirus 5 (HCMV) transcriptome data from two donors, both containing transcriptomes 2 days post infection (2dpi) and 6 days post infection (6dpi). Transcriptome data was found in NCBI, using `wget` to download each SRR run paired end read data. Next, `fasterq-dump` was utilized to uncompress the paired end fastq files and stored in a separate directory to be utilized in the `wrapper.py` script.

SRR Data: [SRR5660030 (Donor 1 - 2dpi)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR5660030&display=data-access)

SRR Data: [SRR5660033 (Donor 1 - 6dpi)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR5660033&display=data_access)

SRR Data: [SRR5660044 (Donor 3 - 2dpi)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR5660044&display=data-access)

SRR Data: [SRR5660045 (Donor 3 - 6dpi)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR5660045&display=data-access)

To download and prepare the full HCMV paired end read data, please run the script below:
```
sh downloadData.sh
```
* **MAKE SURE TO CLONE REPOSITORY BEFORE RUNNING SHELL SCRIPT IF UTILIZING FULL DATASET** *

Also, within this repository there is **sample** data ready if you do not want to download the full dataset and want to test out the `wrapper.py` script. Please utilize the sample data held within the **"sampleData"** folder as the input data while testing `wrapper.py`. The **"sampleData"** folder contains shortened reads of the paired end fastaq files from the full dataset.

## **3. Running the Wrapper**

To utilize this pipeline and wrapper, clone this repository, ensure all dependencies are downloaded, and data is properly stored (unless using **sampleData**). Run this script after preparations are met:
```
python wrapper.py -i <INPUT DATA DIRECTORY> -o <OUTPUT DIRECTORY> -e <EMAIL ADDRESS (for NCBI/Entrez)> -a <GENOME ACCESSION NUMBER> -l <LOG FILE OUTPUT NAME> -s <SUBFAMILY NAME>
```
| `-i | --input` -> Input directory where initial read/sequence data is stored. |

| `-o | --output` -> Name the output directory where all subsequnet directories and files will be kept. |

| `-e | --email` -> Email address to access NCBI Entrez searches. |

| `-a | --accession` -> NCBI accession number of transcriptome to build indicies. |

| `-l | --logFileName` -> Name the log file output. |

| `-s | --subfamily` -> Name subfamily whose genome we are making a BLAST+ database of. |

Keep in mind the `-i | --input` flag will change using the full dataset (`data`) or using the sample dataset (`sampleData`)
