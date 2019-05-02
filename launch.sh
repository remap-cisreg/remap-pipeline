#!/bin/sh
# You need to make sure you have
# python3 and snakemake correctly installed in your compute environment

############################################
# Snakemake command for 
# CONDA environment
# TORQUE as job manager
############################################
#snakemake   --use-singularity --singularity-args "-B /scratch/bballester:/scratch/bballester" --snakefile 2.scripts/workflows/Snakefile_remap_v3.5.py --cores 100  --printshellcmds --cluster-config 2.scripts/workflows/config/mesocentre.json --cluster "sbatch  -A {cluster.project-name}  --job-name {cluster.job-name} -p {cluster.partition}  --ntasks {cluster.ntasks} --cpus-per-task={cluster.thread}  -o {cluster.stdout} -e {cluster.stderr}  --time {cluster.time}    " --keep-going --configfile 2.scripts/workflows/config/Snakefile_config_remap.json --rerun-incomplete --resources res=100


############################################
# Snakemake command for 
# SINGULARITY environment
# SLURM as job manager
############################################

# Loading environment modules (if needed / cluster specific)
module purge
module load userspace/all
module load python3/3.6.3

#-- Homo sapiens example using Singularity
snakemake   --use-singularity --singularity-args "-B /scratch/bballester:/scratch/bballester" --snakefile 2.scripts/workflows/Snakefile_remap_v4.py --cores 100  --printshellcmds --cluster-config 2.scripts/workflows/config/mesocentre_v2.json --cluster "sbatch  -A {cluster.project-name}  --job-name {cluster.job-name} -p {cluster.partition}  --ntasks {cluster.ntasks} --cpus-per-task={cluster.thread}  -o {cluster.stdout} -e {cluster.stderr}  --time {cluster.time}  --mem-per-cpu={cluster.memory} " --keep-going --configfile 2.scripts/workflows/config/Snakefile_config_remap.json --rerun-incomplete --resources res=100  

#### Thaliana example using Singularity
snakemake   --use-singularity --singularity-args "-B /scratch/bballester:/scratch/bballester" --snakefile 2.scripts/snakefiles_workflow/Snakefile_remap_v4.py --cores 100  --printshellcmds --cluster-config 2.scripts/snakemake_configuration/Snakefile_config_remap_thaliana_meso.json --cluster "sbatch  -A {cluster.project-name}  --job-name {cluster.job-name} -p {cluster.partition}  --ntasks {cluster.ntasks} --cpus-per-task={cluster.thread}  -o {cluster.stdout} -e {cluster.stderr}  --time {cluster.time}  --mem-per-cpu={cluster.memory} " --keep-going --configfile cluster_configuration/cluster_singularity_slurm_thaliana.json --rerun-incomplete --resources res=100 

