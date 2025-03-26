# ERVscanner
A data analysis pipeline to estimate ERV insertion based on short-read sequence data in a fast and efficient way.
]
ERVscanner is a pipeline to estimate non-reference ERV insertions based on short-read whole genome sequence data.
The pipeline consists of two parallel workflows, one for detecting insertions within annotated repeat regions (MRs) and one for detecting insertions outside MRs. 
The genotyped VCF file of insertions outside of MRs is output in `<DATA_PATH>/vcf` and the genotyped VCF file of insertions within MRs is output in `<DATA_PATH>/inMR/vcf`. In these VCF files, the flag `MI` tag in the INFO field was added to distinguish these two output. You can merge these files using tools such as `bcftools` or our in-house `merge_vcf.py` script.

ERVscanner consists of seven shell scripts, which should be run concequtively. Some process could be parallelized to increase speed.

## Required Tools and Environment
- Unix-like operationg system, bash
- Python3
- BEDtools 2.27.1
- SAMtools 1.20 or higher
- BWA 

## Input files

Before you start to run the pipeline, you have to prepare the following files for input.

1. BAM or CRAM for each sample
2. Tab-separated dictionary text file (`<DFAM_ERV>` file) you are focusing (You can search copy and paste the content at https://dfam.org/browse.). Here is an example.
   ```
    #Accession Name Classification Clades Description Length
    DF000001785	IAPLTR1a_Mm	ERV2	Mus musculus	Mouse family of LTR retrotransposons	337
    DF000001786	IAPLTR2a2_Mm	ERV2	Mus musculus	Long terminal repeat of ERV2 Endogenous Retrovirus from mouse.	444
   ```
3. A line-separated list of ERV classes you want to analyze, corresponding to the third column of `<DFAM_ERV` file
3. Multi-fasta file of repeat sequences. This file could include all non-target repeat sequences such as SINE and LINE. Including non-ERV sequences decrease the false-positive rate (`<ALL_REPEAT_FASTA>` file)
4. Multi-fasta file of target repeat sequences you want to identify (`<TARGET_REPEAT_FASTA>` file)
5. Line-delimited list of all samples (`<SAMPLE_LIST>` file)
6. BED file of ERV regions obtained from DFAM (`<QUERY_BED>` file). The shell script assumes this file is stored in `<DATA_PATH>` directory.
7. Reference genome sequence (`<REF_GENOME>` file)
8. A line-separated list of alternative chromosomes in the reference genome of the organism (`<ALT_CHR_LIST>` file)

## Description of each shell script

1. `mkdir.sh`
   If you run this shell script in your working directory, it will generate nessesary directories. 
2. `filter_reads.sh`, add_pipe1
   - Extracting both read pairs in which one of the paired ends is mapped to an ERV region in the reference genome (add_pipe1 uses pairs in which both are mapped to MRs)
   - Estimation of insertion position Extract reads that map to non-ERV regions
   - Creating fastq file of RERVreads Extract reads mapped to ERV
3. `remap_reads.sh`
   - Merge and map fastqs Reads attached to ERVs and map them back to the database
   - Create list of insert positions
4. `identify_loci.sh`
   - Matching each sample's insertions with position IDs
5. `filter_loci.sh`
   - Filtering data by insertion sequence estimation and creating a list of insertions across samples
6. `genotype_ins.sh`, add_pipe5
   - Genotyping insetions using the information of boundary-overlapping reads
7. `make_vcf.sh`
   - Creation of final output in VCF file format

## How to run

You first run `mkdir.sh` to prepare directories nessesary for the analysis. ERVscanner produces a lot of intermediate files for checking purpose, but after finishing all procecces, you can delete all intermediate files if you want. Following is the example of command line.
```
bash mkdir.sh -s <SAMPLE_LIST> -d <DATA_PATH>
```
Here, `<DATA_PATH>` is a directory where all output files are stored.

After making the directories, `<QUERY_BED>` file is moved to `<DATA_PATH>`. In addition, `<TARGET_REPEAT_FATA>` file is moved to `<DATA_PATH>/chech_seq/bwa/subject`.

```
mv <QUERY_BED> <DATA_PATH>
mv <TARGET_REPEAT_FASTA> <DATA_PATH>/chech_seq/bwa/subject/
```

After making directories, run `filter_reads.sh`.
```
bash filter_reads.sh -s <SAMPLE_LIST> -r <REF_GENOME> -i <INPUT_PATH> -t <INPUT_TYPE> -d <DATA_PATH> -n <NCORE> -b <QUERY_BED> -q <QUALITY> -c <CLUSTER_THRESHOLD> -a <ALT_CHR_LIST>
```
This process requires many parameters and takes the longest time if samplesize is big. Make sure all nessesary parameters are given.
- -s: A file name of sample list `<SAMPLE_LIST>`
- -r: Reference genome file `<REF_GENOME>`
- -i: Directory where your bam or cram files are stored `<INPUT_PATH>`
- -t: Input type. bam, BAM, cram, CRAM. Case sensitive. default: bam
- -d : A path to data. `<DATA_PATH>` This path should be the same as the path given in `mkdir.sh`.
- -n: Number of core used. default: 1
- -b: BED file defining masked regions (MRs). `QUERY_BED`.
- -q: Threshold to define uniquely mapped reads. default: 30
- -c: Threshold of the number of uniqly mapped reads to define read culster. default: 5
- -a: A line-separated list of alternative chromosomes in the reference genome of the organism. `<ALT_CHR_LIST>`

`filter_reads.sh` can be parallelized. By dividing `<SAMPLE_LIST>` file, you can do it manually.

Run `remap_reads.sh`. 
```
bash remap_reads.sh -f <TARGET_REPEAT_FASTA> -n <NCORE> -d <DATA_PATH>
```
- -f: Multi-fasta file of target repeat sequences you want to identify. `<TARGET_REPEAT_FASTA>`
- -d: A path to data. `<DATA_PATH>` This path should be the same as the path given in `mkdir.sh`.
- -n: Number of core used. default: 1

Run `identify_loci.sh`. By dividing `<SAMPLE_LIST>` file, you can run the process in parallel manually.

```
bash identify_loci.sh -s <SAMPLE_LIST> -d <DATA_PATH> -n <NCORE>
```
- -s: A file name of sample list `<SAMPLE_LIST>`
- -d: A path to data. `<DATA_PATH>` This path should be the same as the path given in `mkdir.sh`.
- -n: Number of core used. default: 1

Run `identify_loci.sh`. By dividing `<SAMPLE_LIST>` file, you can run the process in parallel manually.

```
bash identify_loci.sh -s <SAMPLE_LIST> -d <DATA_PATH> -n <NCORE>
```
- -s: A file name of sample list `<SAMPLE_LIST>`
- -d: A path to data. `<DATA_PATH>` This path should be the same as the path given in `mkdir.sh`.
- -n: Number of core used. default: 1

Run `filter_loci.sh`. This script process all samples at once.
```
bash filter_loci.sh -d <DATA_PATH> -s <SAMPLE_LIST> -p <IDENTITY_THRESHOLD> -e <ERV_CLASS>
```



First, create the directory where the data will be placed in advance with .
Pipeline1 and add_pipe5 can be run in parallel by dividing samples.
Pipeline 2 is a merging process and should be executed after pipeline 1 has been completed for all samples.
Pipeline 3 can be run in parallel with separate samples.
Pipeline 4 is a merging and filtering process and should be executed after Pipeline 3 has been completed for all samples.
Pipeline 5 and add_pipe5 can be run in parallel with separate samples.
Pipeline 6 is a merging process and should be executed after Pipeline 5 or add_pipe5 has been completed for all samples.




Notes
The threshold argument for filtering is not a percentage, but a ratio.
