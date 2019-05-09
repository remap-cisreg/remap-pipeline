# Generate data for 1.metadata

## Introduction
This folder contains the meta-data (annotated GSEs) which will be processed by the pipeline

## Formating the meta-data
Use the script `extract_info_download_v2.py` located in 2.scripts/bla.bla/ to parse your annotated file and this generate the necessary files. 

### Input:  
Give a TAB/Space formated (TSV) file whcih contains

### Output:


### Command line:
module purge
module load userspace/all
module load python3/3.6.3

`python3 ../2.scripts/utils/python3/extract_info_download_v3.py -h`


python3 ../2.scripts/utils/python3/extract_info_download_v3.py -wd /gpfs/tagc/home/simler/remap-pipeline/  -f /gpfs/tagc/home/simler/remap-pipeline/1.metadata/remap3_.tsv  2> .err 


python3 ../2.scripts/utils/python3/extract_info_download_v3.py -wd /scratch/bballester/thaliana/remap-thaliana/ -f /scratch/bballester/thaliana/remap-thaliana/1.metadata/remap3_prod_genome_occupancy_ArabidopsisThaliana.tsv  2> remap3_prod_genome_occupancy_ArabidopsisThaliana.err 


