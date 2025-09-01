#getting paired-end reads from NCBI
sratoolkit.3.2.1-win64\bin\fasterq-dump --split-files SRR23892276
#####################

#Pairwise overlapping reads were merged using FLASH v.1.2.11
module load flash-1.2.11
flash SRR23892276_1.fastq  SRR23892276_2.fastq -o SRR23892276
#############################

#quality check
module load fastqc-0.11.7
fastqc SRR23892276.extendedFrags.fastq
###############################

#trimming bases with quality < Q30
module load sickle
sickle se -f SRR23892276.extendedFrags.fastq -t sanger -o SRR23892276_merged_trimmed.fastq -q 30
###############################

#check quality
fastqc SRR23892276_merged_trimmed.fastq
###############################

#Sequences clustering into operational taxonomic units
#Dereplicate sequences
module load usearch
usearch -fastx_uniques SRR23892276_merged_trimmed.fastq -fastaout uniques.fa -sizeout

#OTU clustering at 97% identity 
#removal of singletons
#chimera filtering
usearch -cluster_otus uniques.fa -otus otus.fa -relabel OTU

#Map reads back to OTUs to calculate abundance table
usearch -usearch_global SRR23892276_merged_trimmed.fastq -db otus.fa -id 0.97 -strand both -otutabout otu_table.txt



