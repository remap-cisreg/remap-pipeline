# Generate data for 1.metadata

## Introduction
This folder contains the meta-data (annotated GSEs) which will be processed by the pipeline

## Formating the meta-data
Use the script `extract_info_download_v3.py` or `extract_info_download_encode_v3.py` located in `2.scripts/utils/python2/` to parse your annotated file and this generate the necessary files.

### Input:  
Give a TAB/Space formated (TSV) file which contains : 
- Experiment ID (GSExxxx)
- Target Name (FOXA1) use HGNC to use offical name convention
- full biotype+treatment (if needed) (MCF-7) use convention names
- replicate(s) ID separated by ",",  (GSMxxx IDs) 
- control(s) ID separate bu ",",   (GSMxxx IDs) 
- ENA ID (SRPXXXX) 

Just like this:

```
GSE104399	ESR1	Breast-tumor_Male_8	GSM2797097	GSM2797159	SRP119087
GSE104399	ESR1	Breast-tumor_Male_9	GSM2797098	GSM2797159	SRP119087
GSE104399	FOXA1	Breast-tumor_Female_1	GSM2797077	GSM2797157	SRP119087
GSE104399	FOXA1	Breast-tumor_Female_2	GSM2797078	GSM2797157	SRP119087
```

This input is the result of your curation and annotation process. 
Annotation and curation should be done with extreme care, using the correct naming conventions, the right Gene names and Biotype names. Regarding cell line treatment, we tend to summarize the information `HUVEC-C_VEGF_12h`, as it just need to be unique from another experiment. We use underscore `_` as a convention to delimit the main biotype `HUVEC-C_VEGF-12h` to the the treatment.

We use HGNC for human genes, and Ensembl gene names for other species. Please use the [EMBL-EBI Ontology Lookup Service](https://www.ebi.ac.uk/ols/index) or the [ExPASy Cellosaurus](https://web.expasy.org/cellosaurus/) from the SIB - Swiss Institute of Bioinformatics  for annotating biotypes. 

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
python3 /absolute/path/to/remap-pipeline/2.scripts/utils/python3/extract_info_download_v3.py -h

python3 /absolute/path/to/remap-pipeline/2.scripts/utils/python3/extract_info_download_v3.py -wd /absolute/path/to/remap-pipeline/  -f /absolute/path/to/remap-pipeline/1.metadata/remap3_.tsv  2> .err

python3 /absolute/path/to/remap-pipeline/2.scripts/utils/python3/extract_info_download_v3.py -wd /absolute/path/to/remap-pipeline/remap-thaliana/ -f /absolute/path/to/remap-pipeline/remap-thaliana/1.metadata/remap3_prod_genome_occupancy_ArabidopsisThaliana.tsv  2> remap3_prod_genome_occupancy_ArabidopsisThaliana.err

python3 ./absolute/path/to/remap-pipeline/2.scripts/utils/python3/extract_info_download_v3.py -wd /absolute/path/to/remap-pipeline/remap-histones-at/  -f /absolute/path/to/remap-pipeline/remap-histones-at/1.metadata/remap3_prod_genome_occupancy_ArabidopsisThaliana_histone.tsv  2>  remap3_prod_genome_occupancy_ArabidopsisThaliana_histone.err`
```
