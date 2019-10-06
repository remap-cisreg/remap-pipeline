# Generate data for 1.metadata

## Introduction
This folder contains the meta-data (annotated GSEs) which will be processed by the pipeline

## Formating the meta-data
Use the script `extract_info_download_v3.py` or `extract_info_download_encode_v3.py` located in `2.scripts/utils/python2/` to parse your annotated file and this generate the necessary files.

### Input:  
Give a TAB/Space formated (TSV) file which contains : 
- Experiment ID (GSExxxx)
- Target Name (FOXA1) use HGNC to use offical name convention
- full biotype (MCF-7) use conventianla names
- replicate(s) ID separated by ",",  (GSMxxx IDs) 
- control(s) ID separate bu ",",   (GSMxxx IDs) 
- ENA ID (SRPXXXX) 

This input is the result of your curation and annotation process. Annotation and curation should be done with care, using the right Gene names and Biotype names. We use HGNC for human genes, and Ensembl gene names for other species. Please use the [EMBL-EBI Ontology Lookup Service](https://www.ebi.ac.uk/ols/index) or the [ExPASy Cellosaurus](https://web.expasy.org/cellosaurus/) from the SIB - Swiss Institute of Bioinformatics  for annotating biotypes. 

### Output:
- A TSV by experiment containing : 
- Replicate name
- Short description found on ENA
- Library type
- if control
- Download URL
- md5

### Command line:
```
# This is specific to our cluster. You may alerady have Python loaded. 
module purge
module load userspace/all
module load python3/3.6.3
```

```
python3 ../2.scripts/utils/python3/extract_info_download_v3.py -h

python3 ../2.scripts/utils/python3/extract_info_download_v3.py -wd /gpfs/tagc/home/simler/remap-pipeline/  -f /gpfs/tagc/home/simler/remap-pipeline/1.metadata/remap3_.tsv  2> .err

python3 ../2.scripts/utils/python3/extract_info_download_v3.py -wd /scratch/bballester/thaliana/remap-thaliana/ -f /scratch/bballester/thaliana/remap-thaliana/1.metadata/remap3_prod_genome_occupancy_ArabidopsisThaliana.tsv  2> remap3_prod_genome_occupancy_ArabidopsisThaliana.err

python3 ../2.scripts/utils/python3/extract_info_download_v3.py -wd /scratch/bballester/thaliana/remap-histones-at/  -f /scratch/bballester/thaliana/remap-histones-at/1.metadata/remap3_prod_genome_occupancy_ArabidopsisThaliana_histone.tsv  2>  remap3_prod_genome_occupancy_ArabidopsisThaliana_histone.err`
```
