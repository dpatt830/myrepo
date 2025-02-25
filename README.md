# Transcriptome Analysis Pipeline

## Overview
This pipeline automates the analysis of the Human herpesevirus 5, also known as Human cytomegalovirus ([HCMV](https://www.ncbi.nlm.nih.gov/nuccore/NC_006273.2)), transcriptomes at **2 and 6 days post-infection (dpi)**. It retrieves sequencing data, processes reads, and identifies differentially expressed genes. Data is derived from [Cheng et al. 2017](https://pubmed.ncbi.nlm.nih.gov/29158406/). This pipeline can be used for other transcriptomes as well.

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

| **Unix/Linux** -> fasterq-dump, kallisto, Bowtie2, SPAdes, NCBI Datasets CLI, BLAST+ |

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

## Examples
If you are utilizing the **sampleData**, here is an example script to run:
```
python wrapper.py -i sampleData -o PipelineProject_Drew_Patterson -e dpatterson3@luc.edu -a NC_006273.2 -l logFile.log -s Betaherpesvirinae
```

If you are utilizing the **full** dataset, here is an example script to run:
```
python wrapper.py -i data -o PipelineProject_Drew_Patterson -e dpatterson3@luc.edu -a NC_006273.2 -l logFile.log -s Betaherpesvirinae
```
**AGAIN** keep in mind before running `wrapper` script, download the full data set by running the shell command listed in **Part 2**

## Extra
Provided below are links to each dependency/library homepage or documentation for more information on how to use as well as install:
* **Unix/Linux** -> [fasterq-dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump), [kallisto](https://github.com/pachterlab/kallisto), [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml), [SPAdes](https://github.com/ablab/spades), [NCBI Datasets CLI](https://github.com/ncbi/datasets), [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)

* **Python** -> [biopython](https://biopython.org), [pandas](https://pandas.pydata.org), [os](https://docs.python.org/3/library/os.html), [subprocess](https://docs.python.org/3/library/subprocess.html)

* **R** -> [sleuth](https://github.com/pachterlab/sleuth), [dplyr](https://dplyr.tidyverse.org)
