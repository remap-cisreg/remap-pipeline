
# Indexes for Bowtie2 and genome in FASTA format

This folder will contain the Bowtie Indexes and genomes in FASTA format. 

We download them from [igenome from Illumina](https://emea.support.illumina.com/sequencing/sequencing_software/igenome.html). 

## H. sapiens, NCBI GRCh38

### Bowtie 2 Genome - H. sapiens, NCBI GRCh38
- http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

```
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz
```

### Need to tar -zxvf AND remove some folders (BWA, etc...) to save som space
```
tar -zxvf Homo_sapiens_UCSC_hg38.tar.gz
```

### Remove annotation folder
```
cd Homo_sapiens/UCSC/hg38/
rm -rf Annotation/
```

### Remove unnecessary Indexes
```
cd Sequence/
rm -rf BowtieIndex/ BWAIndex/ Chromosomes/
cd ../../../../
```

### Remove TAR.GZ file
```
cd /absolute/path/to/remap-pipeline/3.genome/
rm -f Homo_sapiens_UCSC_hg38.tar.gz #-- Clean if you wabt
rm README.txt  #-- Ilumina README
```



# A. thaliana, TAIR10

## Bowtie 2 Genome - Arabidopsis thaliana TAIR10
- http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

```
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Arabidopsis_thaliana/Ensembl/TAIR10/Arabidopsis_thaliana_Ensembl_TAIR10.tar.gz 
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Arabidopsis_thaliana/NCBI/TAIR10/Arabidopsis_thaliana_NCBI_TAIR10.tar.gz
```

### Need to tar -zxvf AND remove some folders (BWA, etc...) to save som space
```
tar -zxvf Arabidopsis_thaliana_NCBI_TAIR10.tar.gz
```

### Remove annotation folder
```
cd Arabidopsis_thaliana/NCBI/TAIR10/
rm -rf Annotation/
```

### Remove unnecessary Indexes
```
cd Sequence/
rm -rf BowtieIndex/ BWAIndex/ Chromosomes/
cd ../../../../
```

### Remove TAR.GZ file
```
rm Arabidopsis_thaliana_Ensembl_TAIR10.tar.gz Arabidopsis_thaliana_NCBI_TAIR10.tar.gz
rm README.txt  #-- Ilumina README
```
