# ERVscanner
A data analysis pipeline to estimate ERV insertion based on short-read sequence data in a fast and efficient way.
---
ERVscanner is a pipeline designed to estimate non-reference ERV (Endogenous Retrovirus) insertions using short-read whole-genome sequencing data.

The pipeline consists of two parallel workflows:

- Detecting insertions within annotated repeat regions in the reference genome (masked region, MRs)
- Detecting insertions outside MRs.

The process described below focuses on identifying insertions outside MRs. The final genotyped VCF file of insertions is output to `<DATA_PATH>/vcf`.

To detect insertions within MRs, you must repeat the pipeline with the following modifications:

- Replace `filter_reads.sh` with `filter_reads_MR.sh`
- Replace `genotype_ins.sh` with `genotype_ins_MR.sh`

Change the `<DATA_PATH>` directory accordingly. Then, rerun the pipeline from the beginning.

Detecting insertions in MRs is more challenging due to the repetitive nature of these regions, leading to a higher false-positive rate. This is why we separated the workflow.

To differentiate between the two output types, `MI` tag is added to the `INFO` field of the VCF files. You can merge the two VCF files using tools such as bcftools or our in-house script, `merge_vcf.py`.

ERVscanner is composed of several shell scripts stored in the `pipeline` directory. These scripts should be executed consecutively. Some steps can be parallelized to improve performance.

## Required Tools and Environment
- Unix-like operationg system, bash
- Python3
- BEDtools 2.27.1
- SAMtools 1.20 or higher
- BWA 0.7.17

Higher version will be fine as long as it works.

Following python librries are also required. You can install those libraries using `pip` or anaconda.

- pandas
- pysam

## Input files

Before you start to run the pipeline, you have to prepare the following files for input.

1. BAM or CRAM for each sample
1. Tab-separated table file (`<DFAM_INFO>` file) of your organism (You can search copy and paste the content at https://dfam.org/browse.). All kinds of repeaences are preferable. Here is an example.
   ```
    DF000001785	IAPLTR1a_Mm	ERV2	Mus musculus	Mouse family of LTR retrotransposons	337
    DF000001786	IAPLTR2a2_Mm	ERV2	Mus musculus	Long terminal repeat of ERV2 Endogenous Retrovirus from mouse.	444
   ```
   The labels of columns are, "Accession", "Name", "Classification", "Clades", "Description", "Length".
1. Line-delimited list of all samples (`<SAMPLE_LIST>` file)
1. Reference genome sequence (`<REF_GENOME>` file)
1. A line-separated list of alternative chromosomes in the reference genome of the organism (`<ALT_CHR_LIST>` file). Here is an example.
   ```
   chr1_GL456210_random
   chr5_JH584299_random
   chr7_GL456219_random
   chrX_GL456233_random
   chrY_JH584300_random
   chrUn_GL456239
   ```
1. Dfam annotation of all and non-redundant repeats. The annotation of Dfam3.8 is downloaded from [Dfam website](https://www.dfam.org/releases/Dfam_3.8/annotations/). There are two types of files in each <ASSEMBLY_VERSION> folder, such as `mm10`. Both `<ASSEMBLY_VERSION>.hits.gz` and `<ASSEMBLY_VERSION>.nrph.hits.gz` should be downloaded. These files correspond to `<ALL_REPEAT>` and `<NRPH_REPEAT>` in `preprocess.sh`.

## Description of shell scripts

1. `preprocess.sh`
   - If you run this shell script in your working directory, it will generate nessesary directories and prepare all necessary files. 
1. `filter_reads.sh`, `filter_reads_MR.sh`
   - Extracting both read pairs in which one of the paired ends is mapped to an ERV region in the reference genome (add_pipe1 uses pairs in which both are mapped to MRs)
   - Estimation of insertion position and extracting reads that map to non-ERV regions
   - Creating fastq file of ERV reads and extracting reads mapped to ERV
1. `remap_reads.sh`
   - Merge and map fastqs Reads attached to ERVs and map them back to the database
   - Create list of insert positions
1. `identify_loci.sh`
   - Matching each sample's insertions with position IDs
1. `filter_loci.sh`
   - Filtering data by insertion sequence estimation and creating a list of insertions across samples
1. `genotype_ins.sh`, `genotype_ins_MR.sh`
   - Genotyping insetions using the information of boundary-overlapping reads
1. `make_vcf.sh`
   - Creation of final output in VCF file format

## How to run

Download all files from github. The command create a directory `ERVscanner`.

```
git clone https://github.com/Qtom-99/ERVscanner.git
```

Fill in the control file, `config.txt`. You can find it in the `pipeline` directory.

The following items should be given.

__SAMPLE__
- A file name of sample list. They should be matched to the basename of CRAM/BAM files. <SAMPLE_LIST>

__DATA_PATH__
- Data directory. All data is stored in this directory. <DATA_PATH>

__INPUT_PATH__
- Directory where your bam or cram files are stored
  
__TYPE__
- Input type. bam, BAM, cram, CRAM. case sensitive
  
__REF_GENOME__
- Reference genome file <REF_GENOME>

__DFAM_INFO__
- A table of Dfam annotation you want to analyze. <DFAM_INFO>

__NRPH_REPEAT__
- Non-redundant repeat annotation in Dfam (<ASSEMBLY>.nrph.hits.gz). <NRPH_REPEAT>

__ALL_REPEAT__
- All repeat annotation in Dfam (<ASSEMBLY>.hits.gz). <ALL_REPEAT>

__ALT_CHR_LIST__
- A line-separated list of alternative chromosomes in the reference genome of the organism. <ALT_CHR_LIST>

__PY_PATH__ 
- Path to the downloaded ERVscanner python directory. ex) `/home/username/ERVscanner/python`

__NCORE__
- Number of core used. default: 1

__MAPQ_THRESHOLD__
- MAPQ threshold to define uniquely mapped reads. default: 30

__CLUSTER_THRESHOLD__
- Threshold of the number of uniqly mapped reads to define read culster. default: 5

__IDENTITY_THRESHOLD__
- Threthold for filtering loci. The fraction of consistent insertion contents and directions in merged dataset. default: 0.7


All fields are required. Except for <SAMPLE_LIST>, the parameters should be consistent across the script.

You first run `preprocess.sh` to prepare directories and files nessesary for the analysis. ERVscanner produces a lot of intermediate files for checking purpose, but after finishing all procecces, you can delete all intermediate files if you want. Following is the example of command line.

```
bash preprocess.sh config.txt
```
The `preprocess.sh` generates a file `<DATA_PATH>/dfam_info/target_class.txt`. This file is a list of class of ERVs you are going to anlayze. You can edit the file if you want to modify targets. You can later change the parameters for threshold setting and re-run the pipeline without running `preprocess.sh`.

After finishing preparation, run `filter_reads.sh`.
```
bash filter_reads.sh config.txt
```
`filter_reads.sh` can be parallelized. By dividing `<SAMPLE_LIST>` file, you can do it manually.

Run `remap_reads.sh`. 
```
bash remap_reads.sh config.txt
```
Run `identify_loci.sh`. By dividing `<SAMPLE_LIST>` file, you can run the process in parallel manually.

```
bash identify_loci.sh config.txt
```
Run `filter_loci.sh`. This script process all samples at once.
```
bash filter_loci.sh config.txt
```
Run `genotype_ins.sh`. This process can be manually parallelized.

```
bash genotype_ins.sh config.txt
```
Run `make_vcf.sh`. 
```
bash make_vcf.sh config.txt
```
The final output for the insertions in MRs is stored in `<DATA_PATH>/vcf`.


