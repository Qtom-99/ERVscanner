# ERVscanner
A data analysis pipeline to estimate ERV insertion based on short-read sequence data in a fast and efficient way.
]
ERVscanner is a pipeline designed to estimate non-reference ERV (Endogenous Retrovirus) insertions using short-read whole-genome sequencing data.

The pipeline consists of two parallel workflows:

1. One for detecting insertions within annotated repeat regions (MRs), and

1. One for detecting insertions outside MRs.

The process described below focuses on identifying insertions outside MRs. The final genotyped VCF file of insertions is output to `<DATA_PATH>/vcf`.

To detect insertions within MRs, you must repeat the pipeline with the following modifications:

Replace `filter_reads.sh` with `filter_reads_MR.sh`

Replace `genotype_ins.sh` with `genotype_ins_MR.sh`

Change the `<DATA_PATH>` directory accordingly
Then, rerun the pipeline from the beginning.

Detecting insertions in MRs is more challenging due to the repetitive nature of these regions, leading to a higher false-positive rate. This is why we separated the workflow.

To differentiate between the two output types, `MI` tag is added to the `INFO` field of the VCF files. You can merge the two VCF files using tools such as bcftools or our in-house script, `merge_vcf.py`.

ERVscanner is composed of several shell scripts stored in the `pipeline` directory. These scripts should be executed consecutively. Some steps can be parallelized to improve performance.

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
3. A line-separated list of ERV classes you want to analyze, corresponding to the third column of `<DFAM_ERV>` file. (`<ERV_CLASS>` file)
4. Multi-fasta file of repeat sequences. This file could include all non-target repeat sequences such as SINE and LINE. Including non-ERV sequences decrease the false-positive rate (`<ALL_REPEAT_FASTA>` file)
5. Multi-fasta file of target repeat sequences you want to identify (`<TARGET_REPEAT_FASTA>` file)
6. Line-delimited list of all samples (`<SAMPLE_LIST>` file)
7. BED file of ERV regions obtained from DFAM (`<QUERY_BED>` file). The shell script assumes this file is stored in `<DATA_PATH>` directory.
8. Reference genome sequence (`<REF_GENOME>` file)
9. A line-separated list of alternative chromosomes in the reference genome of the organism (`<ALT_CHR_LIST>` file)

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

After making the directories, the following files are moved to generated directories.
- `<QUERY_BED>` file is moved to `<DATA_PATH>`.
- `<TARGET_REPEAT_FATA>` file is moved to `<DATA_PATH>/chech_seq/bwa/subject`.
- `<DRAM_INFO>` file is moved to `<DATA_PATH>/dfam_info`.
- `<ERV_CLASS>` file is move to `<DATA_PATH>`.

```
mv <QUERY_BED> <DATA_PATH>
mv <TARGET_REPEAT_FASTA> <DATA_PATH>/chech_seq/bwa/subject/
mv <DFAM_INFO> <DATA_PATH>/dfam_info/
mv <ERV_CLASS> <DATA_PATH>/
```

After making directories, run `filter_reads.sh`.
```
bash filter_reads.sh -s <SAMPLE_LIST> -d <DATA_PATH> -r <REF_GENOME> -i <INPUT_PATH> -t <INPUT_TYPE> -n <NCORE> -b <QUERY_BED> -q <QUALITY> -c <CLUSTER_THRESHOLD> -a <ALT_CHR_LIST>
```
This process requires many parameters and takes the longest time if samplesize is big. Make sure all nessesary parameters are given.
- -s: A file name of sample list `<SAMPLE_LIST>`
- -d : A path to data. `<DATA_PATH>` This path should be the same as the path given in `mkdir.sh`.
- -r: Reference genome file `<REF_GENOME>`
- -i: Directory where your bam or cram files are stored `<INPUT_PATH>`
- -t: Input type. bam, BAM, cram, CRAM. Case sensitive. default: bam
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
bash filter_loci.sh -s <SAMPLE_LIST> -d <DATA_PATH> -p <IDENTITY_THRESHOLD> -e <ERV_CLASS>
```
- -s: A file name of sample list `<SAMPLE_LIST>`
- -d: A path to data. `<DATA_PATH>` This path should be the same as the path given in `mkdir.sh`.
- -p: Threthold for filtering loci. The fraction of consistent insertion contents and directions in merged dataset. default: 0.7
- -e: A line-separated list of ERV classes you want to analyze, corresponding to the third column of `<DFAM_ERV>` file.

Run `genotype_ins.sh`. This process can be manually parallelized.

```
bash genotype_ins.sh -s <SAMPLE_LIST> -d <DATA_PATH> -r <REF_GENOME> -i <INPUT_PATH> -t <INPUT_TYPE> -n <NCORE>
```
- -s: A file name of sample list `<SAMPLE_LIST>`
- -d : A path to data. `<DATA_PATH>` This path should be the same as the path given in `mkdir.sh`.
- -r: Reference genome file `<REF_GENOME>`
- -i: Directory where your bam or cram files are stored `<INPUT_PATH>`
- -t: Input type. bam, BAM, cram, CRAM. Case sensitive. default: bam
- -n: Number of core used. default: 1

Run `make_vcf`. 
```
bash make_vcf.sh -s <SAMPLE_LIST> -d <DATA_PATH>
```
The final output for the insertions in MRs is stored in `<DATA_PATH>/vcf`.


