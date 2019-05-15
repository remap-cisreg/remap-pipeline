"""
Author: Jeanne Chèneby
Affiliation: TAGC
Aim: Workflow ReMap human
Date: 01-12-16
last update: 13-09-2018
Run:snakemake --snakefile 2.scripts/snakefiles_workflow/Snakefile_remap_post_processing.py --printshellcmds --cores 10 --cluster-config config/sacapus.json --cluster "qsub -V -q {cluster.queue} -l nodes={cluster.node}:ppn={cluster.thread} -o {cluster.stdout} -e {cluster.stderr}" --keep-going --configfile config/Snakefile_config_remap_saccapus.json --use-conda

dag: snakemake --snakefile 2.scripts/snakefiles_workflow/Snakefile_remap_post_processing.py --printshellcmds --cores 10 --cluster-config config/sacapus.json --cluster "qsub -V -q {cluster.queue} -l nodes={cluster.node}:ppn={cluster.thread} -o {cluster.stdout} -e {cluster.stderr}" --keep-going --configfile config/Snakefile_config_remap_saccapus.json --dag 2> /dev/null | dot -T svg > dag.svg
rulegraph: snakemake --snakefile Snakefile_remap_v4.py --printshellcmds --cores 10 --cluster-config config/sacapus.json --cluster "qsub -V -q {cluster.queue} -l nodes={cluster.node}:ppn={cluster.thread} -o {cluster.stdout} -e {cluster.stderr}" --keep-going --configfile config/Snakefile_config_remap_saccapus.json --rulegraph 2> /dev/null | dot -T svg > rulegraph.svg


Run (meso): snakemake --snakefile Snakefile_remap_v4.py --printshellcmds --cores 99 --cluster-config config/mesocentre.json --cluster "srun -p {cluster.partition} -N {cluster.node} -n {cluster.thread} -o {cluster.stdout} -e {cluster.stderr}" --keep-going --configfile config/Snakefile_config_remap.json


Latest modification:
  - todo
"""

#================================================================#
#                        Imports/Configuration file              #
#================================================================#

import os
import csv
import pprint



#================================================================#
#     Global variables                                           #
#================================================================#

BASE_DIR = config["working_dir"]["base_dir"]
workdir: BASE_DIR

# Directories
TAB_DIR = config["working_dir"]["tab_dir"]
PREPROCESSING_DIR = config["working_dir"]["preprocessing"]
PEAKCALLING_DIR = config["working_dir"]["outdir"]
QUALITY_DIR = config["working_dir"]["quality"]
BAM_DIR = config["working_dir"]["bam_dir"]
RULE_DIR = config["working_dir"]["rule_dir"]

# Header CSV
CONTROL_HEADER = config["header"]["control"]
FILENAME_HEADER = config["header"]["filename"]
LIBRARY_HEADER = config["header"]["library"]
LIBRARY_URL = config["header"]["url"]
LIBRARY_MD5 = config["header"]["md5"]

# cluster info
# nb_thread_quality_NSC_RSC = cluster_config[ 'quality_NSC_RSC'][ 'thread']


#================================================================#
#                         Includes                               #
#================================================================#


# include: os.path.join(BASE_DIR, RULE_DIR, "merge_bam.rules")
# include: os.path.join(BASE_DIR, RULE_DIR, "quality_NSC_RSC.rules")
# include: os.path.join(BASE_DIR, RULE_DIR, "quality_FRiP.rules")
# include: os.path.join(BASE_DIR, RULE_DIR, "quality_all.rules")

include: os.path.join(BASE_DIR, RULE_DIR, "delete_trim.rules")

#================================================================#
#                     Defining dataset                           #
#================================================================#


# Get all experiments
list_objects_indir = os.listdir( TAB_DIR) # get all files' and folders' names in the current directory
list_exp = []
dict_experiment_chip_filename = {}


for objects_indir in list_objects_indir: # loop through all the files and folders
    # if os.path.isfile( os.path.join( TAB_DIR, objects_indir)): # check whether the current object is a folder or not
    if objects_indir.endswith( "_summary.tab"): # check whether the current object is a folder or not

        experiment_name = objects_indir.rsplit( "_", 1)[ 0]
        list_exp.append( experiment_name)


        with open( os.path.join( TAB_DIR, objects_indir)) as csvfile:


            reader = csv.DictReader( csvfile, delimiter = "\t")
            for row in reader:

                current_filename = row[ FILENAME_HEADER].strip()
                current_experiment_type = row[ CONTROL_HEADER].strip()

                dict_experiment_chip_filename[ experiment_name] = []

                if current_experiment_type == "1":
                    dict_experiment_chip_filename[ experiment_name].append( current_filename)




#================================================================#
#                         Workflow                               #
#================================================================#


rule all:
    input:  expand( os.path.join( QUALITY_DIR, "macs2_{experiment_name}.FRiP_NSC_RSC"), experiment_name = list_exp)



rule merge_bam:
    input:
            bam = lambda wildcards : expand( os.path.join( BAM_DIR, "{chip_name}.bam"), chip_name = dict_experiment_chip_filename[wildcards.experiment_name])
    output:
            temp( os.path.join( PREPROCESSING_DIR,  "merge_bam", "{experiment_name}.bam"))
    singularity:
            config[ "singularity"][ "samtools"]
    conda:
            config[ "conda"][ "samtools"]
    resources:
            res=1
    log:
            os.path.join( PREPROCESSING_DIR,  "merge_bam", "log", "{experiment_name}_merge.log")
    params:
            samtools_merge = config[ "merge_bam"][ "samtools-merge"]

    shell:"""samtools merge {params.samtools_merge} {output} {input.bam}"""



rule quality_NSC_RSC:
    input:
            bam = os.path.join( PREPROCESSING_DIR, "merge_bam", "{experiment_name}.bam")
    output:
            temp( os.path.join( PREPROCESSING_DIR, "NSC_RSC", "{experiment_name}.NSC_RSC"))
    singularity:
            config[ "singularity"][ "phantompeak"]
    conda:
            config[ "conda"][ "phantompeak"]
    resources:
            res=1
    log:
            os.path.join( PREPROCESSING_DIR, "NSC_RSC", "log", "{experiment_name}_NSC_RSC.log")
    params:
            outdir = os.path.join( PREPROCESSING_DIR, "NSC_RSC"),
            run_spp = config[ "quality_NSC_RSC"][ "run_spp"]
    shell:"""
    Rscript /usr/bin/phantompeakqualtools-1.2/run_spp.R -c={input.bam} -savp -out={output} -p=2 -odir={params.outdir} -rf {params.run_spp}
    """



rule quality_FRiP:
    input:
            # bam = lambda wildcards : expand( os.path.join( BAM_DIR, "{chip_name}.bam"), chip_name = dict_experiment_chip_filename[wildcards.experiment_name]["chip"]),
            bam = os.path.join( PREPROCESSING_DIR, "merge_bam", "{experiment_name}.bam"),
            peaks = os.path.join( PEAKCALLING_DIR, "{experiment_name}", "macs2", "{experiment_name}_peaks.narrowPeak"),
            nsc_rsc = os.path.join( PREPROCESSING_DIR, "NSC_RSC", "{experiment_name}.NSC_RSC")
    output:
            temp( os.path.join( PREPROCESSING_DIR, "FRiP", "macs2_{experiment_name}.FRiP"))
    singularity:
            config[ "singularity"][ "bedtools"]
    conda:
            config[ "conda"][ "bedtools"]
    resources:
            res=1
    log:
            os.path.join( PREPROCESSING_DIR, "FRiP", "log", "{experiment_name}_FRiP.log")
    params:
            bedtools_intersect = config[ "quality_FRiP"][ "bedtools-intersect"]

    shell:"""
        total_reads=$(cut -f2 {input.nsc_rsc})

        total_peaks_in_read=$(bedtools intersect {params.bedtools_intersect} -bed -u -a {input.bam} -b {input.peaks} | wc -l)

        FRiP=$(awk -v nr="$total_reads" -v ni="$total_peaks_in_read" 'BEGIN{{print ni/nr * 100}}')
        echo $FRiP > {output}
    """


rule quality_all:
    input:
            nsc_rsc = os.path.join( PREPROCESSING_DIR, "NSC_RSC", "{experiment_name}.NSC_RSC"),
            frip = os.path.join( PREPROCESSING_DIR, "FRiP", "macs2_{experiment_name}.FRiP")
    output:
            os.path.join( QUALITY_DIR, "macs2_{experiment_name}.FRiP_NSC_RSC")
    resources:
            res=1
    log:
            os.path.join( QUALITY_DIR, "log", "{experiment_name}_quality_all.log")
    params:
            outdir = os.path.join( PREPROCESSING_DIR, "NSC_RSC"),
            exp_name = "{experiment_name}"
    shell:"""
        FRiP=$(cat {input.frip})
        NSC=$(cut -f9 {input.nsc_rsc})
        RSC=$(cut -f10 {input.nsc_rsc})

        echo "experiment_name\tFRiP\tNSC\tRSC\n{params.exp_name}\t$FRiP\t$NSC\t$RSC\n" > {output}
    """
