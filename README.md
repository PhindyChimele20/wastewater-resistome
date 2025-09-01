#  README

## 1. Data sources

This task is based on on publicly available sequencing data from the study "Metagenomic insights into the wastewater resistome before and after purification at large‐scale wastewater treatment plants in the Moscow city" which performed a metagenomic analysis of the raw wastewater, activated sludge and treated wastewater from two large WWTPs responsible for the treatment of urban wastewater in Moscow, Russia. Metagenomic DNA was sequenced using the Illumina HiSeq2500 platform. For this workflow we selected 1 metagenomic sample.

---

## 2. How to download

The raw data for this project is available in the NCBI Sequence Read Archive (SRA) and are available via the BioProject PRJNA945245.
### Code for downloading

```bash
#getting paired-end reads from NCBI
sratoolkit.3.2.1-win64\bin\fasterq-dump --split-files SRR23892276
#####################

#Pairwise overlapping reads were merged using FLASH v.1.2.11
module load flash-1.2.11
flash SRR23892276_1.fastq  SRR23892276_2.fastq -o SRR23892276
```


## 3. How the workflow works
The workflow files is stored in workflow/ and it is divided into different steps:
The workflow files are stored in `workflow/`.

---

### Step 1 – Quality Check

**Purpose:** The workflow takes each FASTQ.qz file (raw reads), assess the quality of the reads and give the scores and overall stats on the quality of reads.
**Tools:** `fastqc`
**Inputs:** Raw reads FASTQ files (from `data/`)
**Outputs:** quality matrix (html)
**Command:**

```bash
module load fastqc-0.11.7
fastqc SRR23892276.extendedFrags.fastq                                        

```

---

### Step 2 - Reads Cleaning/Trimming

**Purpose:** trimming bases with quality < Q30
**Tools:** 'sickle'
**Inputs:** fastq files
**Outputs:** trimmed fastq files
**Command:**

```bash
module load sickle
sickle se -f SRR23892276.extendedFrags.fastq -t sanger -o SRR23892276_merged_trimmed.fastq -q 30

```
---

### Step 3 – Dereplicate sequences

**Purpose:** Sequences clustering into operational taxonomic units
**Tools:** 'usearch'
**Inputs:** trimmed fastq file
**Outputs:** uniques.fa
**Command:**
```bash
module load usearch
usearch -fastx_uniques SRR23892276_merged_trimmed.fastq -fastaout uniques.fa -sizeout

```
### Step 4 - Chimera filtering

**Purpose:** OTU clustering at 97% identity and removal of singletons
**Tools:** 'usearch'
**Inputs:** .sam file
**Outputs:** .bam file, .bam.bai file
**Command:**

```bash
usearch -cluster_otus uniques.fa -otus otus.fa -relabel OTU


```
---
### Step 5 - Calculate abundance

**Purpose:** Map reads back to OTUs to calculate abundance table
**Tools:** 'usearch'
**Inputs:** trimmed fastq file
**Outputs:** otus.fa file, .txt file
**Command:**

```bash
usearch -usearch_global SRR23892276_merged_trimmed.fastq -db otus.fa -id 0.97 -strand both -otutabout otu_table.txt


```
---
