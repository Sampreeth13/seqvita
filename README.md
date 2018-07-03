# SeqVItA: A Fast and Efficient Pipeline for Detection and Annotation of Sequence Variants in Next Generation Sequencing Data

Sequence Variants Identification and Annotation Pipeline

## Requirements
Please add the following tools to your PATH variable 
* samtools ([link](https://sourceforge.net/projects/samtools/))
* bedtools (>=2.26) ([link](http://bedtools.readthedocs.org/en/latest/content/installation.html))
* R-packages - dplyr, tidyr 
* Bioconductor packages - VariantAnnotation, rfPred, SNPlocs.Hsapiens.dbSNP144.GRCh37

## Installation
Download the source code from https://github.com/sampreeth13/seqvita, extract the zip file

```
unzip SeqVItA.zip
cd SeqVItA

```
## Usage

Variant Calling can be carried out using any of the following modules:
* germline (single sample/multiple samples)
* somatic (case-control samples)
* population (multiple samples)

```
# VARIANT CALLING FROM ALIGNMENT FILE
./variantCalling -v SNP -ib Test.bam -r hg19.fa -o output [options]
./variantCalling -v INDEL -ib Test.bam -r hg19.fa -o output [options]
./variantCalling -v germline -ib Test.bam -r hg19.fa -o output [options]
./variantCalling -v somatic --normal normal.bam --tumor tumor.bam -r hg19.fa -o output [options]
./variantCalling -v population -ib 1.bam 2.bam 3.bam 4.bam ... -r hg19.fa -o output [options]

# VARIANT CALLING FROM MPILEUP FILE

./variantCalling -v SNP -im Test.mpileup -o output [options] [mpileup should contain all the samples information]
./variantCalling -v INDEL -im Test.mpileup -o output [options] [mpileup should contain all the samples information]
./variantCalling -v germline -im Test.mpileup -o output [options] [mpileup should contain all the samples information]
./variantCalling -v somatic -im Test.mpileup -o output [options] [mpileup should contain normal-tumor information in the same order]
./variantCalling -v population -im Test.mpileup -o output [options] [mpileup should contain all the samples information]

# VARIANT CALLING FROM WES/TS FILES

./variantCalling -v SNP -ib Test.bam -r hg19.fa --bed Coordinate_file -o output
./variantCalling -v INDEL -ib Test.bam -r hg19.fa --bed Coordinate_file -o output 
./variantCalling -v germline -ib Test.bam -r hg19.fa --bed Coordinate_file -o output
./variantCalling -v somatic --normal normal.bam --tumor tumor.bam --bed Coordinate_file -r hg19.fa  -o output
./variantCalling -v population -ib 1.bam 2.bam 3.bam 4.bam ... -r hg19.fa --bed Coordinate_file -o output [options]

# VARIANT ANNOTATION

./annotate -i Test.vcf -d genebased -o Output [Default:Variant based (rs ID) drug mapping]

```
## Options available in SeqVItA

| Key Parameter | Parameter Description | Default Value |
|---|---|---|
| --Mqread | Mapping quality Cut-off (Only when alignment in BAM format is used as input)| 20 |
| --Mqcorr | Mapping quality correction using Samtools (Only when alignment in BAM format is used as input) | 0 |
| --Qbase |	Base quality Cut-off | 15 |
| --RD_th |	If no. of reads at a position > RD_th, the site is considered for variant calling |	10 |
| --VAR_th |	If no. of reads supporting alternate allele at a position > VAR_th | 2 |
| --VAF_th | Variant allele frequency  cut-off	 | 0.20 |
|	--Strand_Bias	|Ignore variants with > 90% support is from the same strand |	1 |
| --p-value	| p-value cut-off for calling variants | 0.01 |
| --VAF_homo	| Variant with VAF > RD_alt, the variant is homozygous, else heterozygous	| 0.75 |
| --somatic-p-value | p-value cut-off for calling somatic and LOH variants (Only used in somatic module)| 0.05 |

