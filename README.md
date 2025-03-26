# ERVscanner
A data analysis pipeline to estimate ERV insertion based on short-read sequence data in a fast and efficient way.
]
ERVscanner is a pipeline designed to estimate non-reference ERV (Endogenous Retrovirus) insertions using short-read whole-genome sequencing data.

The pipeline consists of two parallel workflows:

- Detecting insertions within annotated repeat regions in the reference genome (MRs)
- Detecting insertions outside MRs.

The process described below focuses on identifying insertions outside MRs. The final genotyped VCF file of insertions is output to `<DATA_PATH>/vcf`.

To detect insertions within MRs, you must repeat the pipeline with the following modifications:

- Replace `filter_reads.sh` with `filter_reads_MR.sh`
- Replace `genotype_ins.sh` with `genotype_ins_MR.sh`

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
- BWA 0.7.17

Higher version will be fine as long as it works.

Following python librries are also required. You can install those libraries using `pip` or anaconda.

- pandas



## Input files

Before you start to run the pipeline, you have to prepare the following files for input.

1. BAM or CRAM for each sample
1. Tab-separated dictionary text file (`<DFAM_INFO>` file) you are focusing (You can search copy and paste the content at https://dfam.org/browse.). Here is an example.
   ```
    #Accession Name Classification Clades Description Length
    DF000001785	IAPLTR1a_Mm	ERV2	Mus musculus	Mouse family of LTR retrotransposons	337
    DF000001786	IAPLTR2a2_Mm	ERV2	Mus musculus	Long terminal repeat of ERV2 Endogenous Retrovirus from mouse.	444
   ```
1. A line-separated list of ERV classes you want to analyze, corresponding to the third column of `<DFAM_ERV>` file. (`<ERV_CLASS>` file)
1. Multi-fasta file of repeat sequences. This file could include all non-target repeat sequences such as SINE and LINE. Including non-ERV sequences decrease the false-positive rate (`<ALL_REPEAT_FASTA>` file)
1. Line-delimited list of all samples (`<SAMPLE_LIST>` file)
1. BED file of ERV regions obtained from DFAM (`<QUERY_BED>` file). The shell script assumes this file is stored in `<DATA_PATH>` directory.
1. Reference genome sequence (`<REF_GENOME>` file)
1. A line-separated list of alternative chromosomes in the reference genome of the organism (`<ALT_CHR_LIST>` file)

#### How to prepare <QUERY_BED>

Annotation of repeat regions (`*.hits.gz`) can be downloaded from the website of Dfam. Here is a link to the version 3.8 release ([Dfam 3.8 annotation](https://dfam.org/releases/Dfam_3.8/annotations/)).

First, download `<ASSEMBLY>.hits.gz`, such as `mm10.hits.gz` and run the following command. The python script `wordgrep.py` is in `utils` directory.

<ERV_LIST> is a line-separated non-redundant text file listing `NAME` column of Dfam annotation (see `<DFAM_INFO>`), which you want to find. You can make the file using the following command


```
cut -f 2 <DFAM_INFO> | uniq > <ERV_LIST>
```
The <ERV_LIST> file is used for the input of `wordgrep.py` script.
```
python wordgrep.py <ASSEMBLY>.hits.gz <ERV_LIST> <OUTPUT_FILE.hits.gz>
```

Next, run `make_bed.py` in `utils` directory to convert the file to BED format.

```
python make_bed.py <ASSEMBLY>.hits.gz <QUERY_BED>
```
The output file is `<QUERY_BED>` file.

#### How to prepare <ALL_REPEAT_FASTA>

Download `<ASSEMBLY>.nrph.hits.gz`, such as `mm10.nrph.hits.gz` and run `make_bed.py` in `utils` directory.
```
python make_bed.py <ASSEMBLY>.nrph.hits.gz <OUTPUT.bed>
```
Sort the bed file.
```
sort -V -k1,1 -k2,2 <OUTPUT.bed> | bedtools merge -d 50 - > <MERGED.bed>
```
Then generate multi-fasta file.
```
bedtools getfasta -fi <REF_GENOME> -name+ -bed <MERGED.bed> -s -fo <ALL_REPEAT_FASTA>
```
The fasta file is later move to a different directory.

## Description of each shell script

1. `mkdir.sh`
   If you run this shell script in your working directory, it will generate nessesary directories. 
1. `filter_reads.sh`, `filter_reads_MR.sh`
   - Extracting both read pairs in which one of the paired ends is mapped to an ERV region in the reference genome (add_pipe1 uses pairs in which both are mapped to MRs)
   - Estimation of insertion position Extract reads that map to non-ERV regions
   - Creating fastq file of RERVreads Extract reads mapped to ERV
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

You first run `mkdir.sh` to prepare directories nessesary for the analysis. ERVscanner produces a lot of intermediate files for checking purpose, but after finishing all procecces, you can delete all intermediate files if you want. Following is the example of command line.
```
bash mkdir.sh -s <SAMPLE_LIST> -d <DATA_PATH>
```
Here, `<DATA_PATH>` is a directory where all output files are stored.

After making the directories, the following files are moved to generated directories.
- `<QUERY_BED>` file is moved to `<DATA_PATH>`.
- `<ALL_REPEAT_FATA>` file is moved to `<DATA_PATH>/check_seq/bwa/subject`.
- `<DRAM_INFO>` file is moved to `<DATA_PATH>/dfam_info`.
- `<ERV_CLASS>` file is move to `<DATA_PATH>`.

```
mv <QUERY_BED> <DATA_PATH>
mv <ALL_REPEAT_FASTA> <DATA_PATH>/chech_seq/bwa/subject/
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
bash remap_reads.sh -f <ALL_REPEAT_FASTA> -n <NCORE> -d <DATA_PATH>
```
- -f: Multi-fasta file of target repeat sequences you want to identify. `<ALL_REPEAT_FASTA>`
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

Run `make_vcf.sh`. 
```
bash make_vcf.sh -s <SAMPLE_LIST> -d <DATA_PATH>
```
The final output for the insertions in MRs is stored in `<DATA_PATH>/vcf`.


