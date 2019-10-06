#!/bin/sh

# Example using TORQUE ans Conda
snakemake --use-conda--snakefile 2.scripts/snakefiles_workflow/Snakefile_workflow.py --cores 100 --cluster-config 2.scripts/cluster_configuration/cluster.json --cluster "qsub -V -q {cluster.queue} -l nodes={cluster.node}:ppn={cluster.thread} -o {cluster.stdout} -e {cluster.stderr}"  --configfile 2.scripts/snakemake_configuration/Snakemake_config.json --resources res=101
