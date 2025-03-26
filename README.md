# ERVscanner
A data analysis pipeline to estimate ERV insertion based on short-read sequence data in a fast and efficient way.
]
ERVscanner is a pipeline to estimate non-reference ERV insertions based on short-read whole genome sequence data.
The pipeline consists of two parallel workflows, one for detecting insertions within annotated repeat regions (MRs) and one for detecting insertions outside MRs. 
The genotyped VCF file of insertions outside of MRs is output in <your_datapath>/vcf and the genotyped VCF file of insertions within MRs is output in <your_datapath>/inMR/vcf. In these VCF files, the flag `MI` tag in the INFO field was added to distinguish these two output. You can merge these files using tools such as `bcftools` or our in-house `merge_vcf.py` script.

ERVscanner consists of seven shell scripts, which should be run concequtively.

## Required Tools and Environment
- Unix-like operationg system
- Python3
- bedtools 2.27.1
- samtools 1.20 or higher
- bwa 

## Input files

Before you start to run the pipeline, you have to prepare the following files for input.

1. BAM or CRAM for each sample
2. Tab-separated DFAM_ERV dictionary text file you are focosing (You can search copy and paste the content at https://dfam.org/browse.). Here is an example.
```
    #Accession Name Classification Clades Description Length
    DF000001785	IAPLTR1a_Mm	ERV2	Mus musculus	Mouse family of LTR retrotransposons	337
    DF000001786	IAPLTR2a2_Mm	ERV2	Mus musculus	Long terminal repeat of ERV2 Endogenous Retrovirus from mouse.	444
```
3. Multi fasta file of repeat sequences. This file could include all non-target repeat sequences such as ALU and LINE. Including non-ERV sequences decrease the false-positive rate.
4. Line-delimited list of all samples (sample list file)
5. BED file of ERV regions obtained from DFAM
6. Reference genome sequence
7. A line-separated list of alternative chromosomes in the reference genome of the organism to be analyzed

## Description of each shell script

1. `mkdir.sh`
   If you run this shell script in your working directory, it will generate nessesary directories. The 
1. filter_read, add_pipe1
  process1 Extract both read pairs in which one of the paired ends is mapped to an ERV region in the reference genome (add_pipe1 uses pairs in which both are mapped to MRs)
  process2 Estimation of insertion position Extract reads that map to non-ERV regions.
  process3 Create fastq file of RERVreads Extract reads mapped to ERV
1. Pipeline2
  process4 Merge and map fastqs Reads attached to ERVs and map them back to the database
  process5 Create list of insert positions
1.  Pipeline3
  process6 Match each sample's insertions with position IDs
1. Pipeline4
  process7 Filter by insertion sequence estimation, create cross table of 01
1. Pipeline5, add_pipe5
  process8 Genotyping
1. Pipeline6
  process9 Creation of VCF

## How to run

First, create the directory where the data will be placed in advance with .
Pipeline1 and add_pipe5 can be run in parallel by dividing samples.
Pipeline 2 is a merging process and should be executed after pipeline 1 has been completed for all samples.
Pipeline 3 can be run in parallel with separate samples.
Pipeline 4 is a merging and filtering process and should be executed after Pipeline 3 has been completed for all samples.
Pipeline 5 and add_pipe5 can be run in parallel with separate samples.
Pipeline 6 is a merging process and should be executed after Pipeline 5 or add_pipe5 has been completed for all samples.




Notes
The threshold argument for filtering is not a percentage, but a ratio.
