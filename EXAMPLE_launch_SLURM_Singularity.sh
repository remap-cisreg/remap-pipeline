#!/bin/sh

# Example using SLURM and Docker with 100 parallele jobs
snakemake  --use-singularity --singularity-args "-B /where/to/mount:/what/to/mount" --snakefile 2.scripts/snakefiles_workflow/Snakefile_workflow.py --cores 100 --cluster-config 2.scripts/cluster_configuration/cluster.json --cluster "sbatch  -A {cluster.project-name}  --job-name {cluster.job-name} -p {cluster.partition}  --ntasks {cluster.ntasks} --cpus-per-task={cluster.thread}  -o {cluster.stdout} -e {cluster.stderr}  --time {cluster.time}  --mem-per-cpu={cluster.memory} " --configfile 2.scripts/snakemake_configuration/Snakemake_config.json --resources res=101
