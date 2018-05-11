# SeqVItA: A Fast and Efficient Pipeline for Detection and Annotation of Sequence Variants in Next Generation Sequencing Data

Sequence Variants Identification and Annotation Pipeline

## Requirements
Please add the following tools to your PATH variable 
* samtools ([link](https://sourceforge.net/projects/samtools/))
* bedtools (>=2.26) ([link](http://bedtools.readthedocs.org/en/latest/content/installation.html))
* R-packages - VariantAnnotation, rfPred, dplyr, tidyr, SNPlocs.Hsapiens.dbSNP144.GRCh37

## Installation
Download the source code from https://github.com/sampreeth13/seqvita, extract the zip file

```
unzip SeqVItA.zip
cd SeqVItA
make 
```
## Usage

```
# VARIANT CALLING FROM ALIGNMENT FILE
./variantCalling snp -ib Test.bam -r hg19.fa -o output.vcf
./variantCalling indel -i Test.bam ....
./variantCalling germline -i Test.bam
./variantCalling somatic --normal normal.bam --tumor tumor.bam ...

# VARIANT CALLING FROM MPILEUP FILE

./variantCalling snp -im Test.mpileup -o output.vcf
./variantCalling indel -i Test.mpileup ....
./variantCalling germline -i Test.mpileup
./variantCalling somatic --normal normal.mpileup --tumor tumor.mpileup ...

# VARIANT CALLING FROM WES/TS FILES

./variantCalling snp -ib Test.bam -r hg19.fa --bed Coordinate_file -o output.vcf
./variantCalling indel -i Test.bam --bed Coordinate_file ....
./variantCalling germline -i Test.bam --bed Coordinate_file
./variantCalling somatic --normal normal.bam --tumor tumor.bam --bed Coordinate_file ...


# VARIANT ANNOTATION

./annotate -i Test.vcf --geneBasedDrug -o Output [Default:Variant based (rs ID) drug mapping]

```
