#!/bin/sh

##############################
# chargement des modules
##############################
module purge
module load userspace/all
module load python3/3.6.3


##############################
# Snakemake command OLD - prod1
##############################
#snakemake   --use-singularity --singularity-args "-B /scratch/bballester:/scratch/bballester" --snakefile 2.scripts/workflows/Snakefile_remap_v3.5.py --cores 100  --printshellcmds --cluster-config 2.scripts/workflows/config/mesocentre.json --cluster "sbatch  -A {cluster.project-name}  --job-name {cluster.job-name} -p {cluster.partition}  --ntasks {cluster.ntasks} --cpus-per-task={cluster.thread}  -o {cluster.stdout} -e {cluster.stderr}  --time {cluster.time}    " --keep-going --configfile 2.scripts/workflows/config/Snakefile_config_remap.json --rerun-incomplete --resources res=100


##############################
# Snakemake command - prod 2
##############################
snakemake   --use-singularity --singularity-args "-B /scratch/bballester:/scratch/bballester" --snakefile 2.scripts/workflows/Snakefile_remap_v4.py --cores 100  --printshellcmds --cluster-config 2.scripts/workflows/config/mesocentre_v2.json --cluster "sbatch  -A {cluster.project-name}  --job-name {cluster.job-name} -p {cluster.partition}  --ntasks {cluster.ntasks} --cpus-per-task={cluster.thread}  -o {cluster.stdout} -e {cluster.stderr}  --time {cluster.time}  --mem-per-cpu={cluster.memory} " --keep-going --configfile 2.scripts/workflows/config/Snakefile_config_remap.json --rerun-incomplete --re\sources res=100  




