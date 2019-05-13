# remap-pipeline

## Introduction

## Quick start guide

### Prepare the meta-data (the datasets)
 - See README.md in 1.metadata/


### Download Genome
 - From Illumina for Human (See README.md in 3.genome/)
 - From Illumina/NCBI/Ensembl for Thaliana (See README.md in 3.genome/)

## Compute clusters specificities

### Conda environement

### Singularity environement


#### Homo sapiens example using Singularity
`snakemake   --use-singularity --singularity-args "-B /scratch/bballester:/scratch/bballester" --snakefile 2.scripts/snakefiles_workflow/Snakefile_remap_v4.py --cores 100  --printshellcmds --cluster-config 2.scripts/snakemake_configuration/Snakefile_config_remap_hsap_meso.json --cluster "sbatch  -A {cluster.project-name}  --job-name {cluster.job-name} -p {cluster.partition}  --ntasks {cluster.ntasks} --cpus-per-task={cluster.thread}  -o {cluster.stdout} -e {cluster.stderr}  --time {cluster.time}  --mem-per-cpu={cluster.memory} " --keep-going --configfile cluster_configuration/cluster_singularity_slurm_hsap.json --rerun-incomplete --resources res=100  `

#### Thaliana example using Singularity
`snakemake   --use-singularity --singularity-args "-B /scratch/bballester:/scratch/bballester" --snakefile 2.scripts/snakefiles_workflow/Snakefile_remap_v4.py --cores 100  --printshellcmds --cluster-config 2.scripts/snakemake_configuration/Snakefile_config_remap_thaliana_meso.json --cluster "sbatch  -A {cluster.project-name}  --job-name {cluster.job-name} -p {cluster.partition}  --ntasks {cluster.ntasks} --cpus-per-task={cluster.thread}  -o {cluster.stdout} -e {cluster.stderr}  --time {cluster.time}  --mem-per-cpu={cluster.memory} " --keep-going --configfile cluster_configuration/cluster_singularity_slurm_thaliana.json --rerun-incomplete --resources res=100  `



### Slurm workload manager

### PBS (Torque) workload manager 

#### Make a BED file from all narrowpeaks
```
cd 6.peakcalling/
ls -1  > ../list_experiments.txt 

#-- Exec script
2.scripts/utils/bash/coloring_tf_bed.sh  |  sort -k1,1n -k2,2n  > at_macs2_nomodel.bed
```



## People involved
- Jeanne Cheneby
- Lionel Spinelli
- Benoit Ballester

## Contact information


## References
