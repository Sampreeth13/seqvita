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
./VariantCalling -v SNP -ib Test.bam -r hg19.fa -o output [options]
./VariantCalling -v INDEL -ib Test.bam -r hg19.fa -o output [options]
./VariantCalling -v germline -ib Test.bam -r hg19.fa -o output [options]
./VariantCalling -v somatic --normal normal.bam --tumor tumor.bam -r hg19.fa -o output [options]
./VariantCalling -v population -ib 1.bam 2.bam 3.bam 4.bam ... -r hg19.fa -o output [options]

# VARIANT CALLING FROM MPILEUP FILE

./VariantCalling -v SNP -im Test.mpileup -o output [options] [mpileup should contain all the samples information]
./VariantCalling -v INDEL -im Test.mpileup -o output [options] [mpileup should contain all the samples information]
./VariantCalling -v germline -im Test.mpileup -o output [options] [mpileup should contain all the samples information]
./VariantCalling -v somatic -im Test.mpileup -o output [options] [mpileup should contain normal-tumor information in the same order]
./VariantCalling -v population -im Test.mpileup -o output [options] [mpileup should contain all the samples information]

# VARIANT CALLING FROM WES/TS FILES

./VariantCalling -v SNP -ib Test.bam -r hg19.fa --bed Coordinate_file -o output
./VariantCalling -v INDEL -ib Test.bam -r hg19.fa --bed Coordinate_file -o output 
./VariantCalling -v germline -ib Test.bam -r hg19.fa --bed Coordinate_file -o output
./VariantCalling -v somatic --normal normal.bam --tumor tumor.bam --bed Coordinate_file -r hg19.fa  -o output
./VariantCalling -v population -ib 1.bam 2.bam 3.bam 4.bam ... -r hg19.fa --bed Coordinate_file -o output [options]

# VARIANT ANNOTATION

./Annotate -i Test.vcf -d genebased -o Output [Default:Variant based (rs ID) drug mapping]

```
## Options available in SeqVItA

| Key Parameter | Parameter Description | Default Value |
|---|---|---|
| --Mqread | Mapping quality Cut-off (Only when alignment in BAM format is used as input)| 20 |
| --Mqcorr | Mapping quality correction using Samtools (Only when alignment in BAM format is used as input)(A value of 50 may be considered in case the BAM file is generated using Bowtie or BWA methods | 0 |
| --Qbase |	Base quality Cut-off | 15 |
| --RD_th |	If no. of reads at a position > RD_th, the site is considered for variant calling |	10 |
| --VAR_th |	If no. of reads supporting alternate allele at a position > VAR_th | 2 |
| --VAF_th | Variant allele frequency  cut-off	 | 0.20 |
|	--Strand_Bias	|Ignore variants with > 90% support is from the same strand |	1 |
| --p-value	| p-value cut-off for calling variants | 0.01 |
| --VAF_homo	| Variant with VAF > RD_alt, the variant is homozygous, else heterozygous	| 0.75 |
| --somatic-p-value | p-value cut-off for calling somatic and LOH variants (Only used in somatic module)| 0.05 |

