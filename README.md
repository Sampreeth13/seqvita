# SeqVItA: A Fast and Efficient Pipeline for Detection and Annotation of Sequence Variants in Next Generation Sequencing Data

Sequence Variants Identification and Annotation Pipeline

## Requirements
Please add the following tools to your PATH variable 
* samtools ([link](https://sourceforge.net/projects/samtools/))
* bedtools (>=2.26) ([link](http://bedtools.readthedocs.org/en/latest/content/installation.html))
* R-packages - VariantAnnotation

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
./calOptBinSize -c <config file> -i <input BAM>

# VARIANT CALLING FROM MPILEUP FILE
./prepareData -m <mappability file> -g <gc content file> --win <desired window size> --genome_file <Genome file> -o <Output file name prefix>

# VARIANT ANNOTATION
./pretreatment -i <input BAM> -z <bed file containing bins> --mapfile < file contaning mappability values> -o <output prefix> --gcfile <file containing GC scores>

```
