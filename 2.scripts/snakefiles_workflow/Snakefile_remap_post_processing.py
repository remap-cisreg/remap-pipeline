"""
Author: Jeanne Chèneby
Affiliation: TAGC
Aim: Workflow ReMap human
Date: 01-12-16
last update: 13-09-2018
Run : snakemake --snakefile 2.scripts/snakefiles_workflow/Snakefile_remap_post_processing.py --printshellcmds --cores 10 --cluster-config 2.scripts/cluster_configuration/cluster_conda_torque_remap_post_processing.json --cluster "qsub -V -q {cluster.queue} -l nodes={cluster.node}:ppn={cluster.thread} -o {cluster.stdout} -e {cluster.stderr}" --keep-going --configfile 2.scripts/snakemake_configuration/Snakefile_config_remap_sacapus_post_processing.json --use-cond

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

PEAKCALLER = config[ "info"][ "peakcaller"]

# Extension peakcaller
EXTENSION_PEAK = config[ "extention"][ "peak"]

# File names
REMAP_FILENAME_PREFIX = "remap" + \
                        config["info"]["remap_version"]
REMAP_FILENAME_SUFFIX = PEAKCALLER +\
                        "_" + config["info"]["assembly"] + \
                        "v" + config["info"]["update_version"] + \
                        "_" + config["info"]["patch_version"]

REMAP_FULLANME = REMAP_FILENAME_PREFIX + "_" + REMAP_FILENAME_SUFFIX

#================================================================#
#                         Includes                               #
#================================================================#


include: os.path.join(BASE_DIR, RULE_DIR, "merge_bam.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "quality_NSC_RSC.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "quality_FRiP_nbPeaks.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "quality_all_experiment.rules")
#include: os.path.join(BASE_DIR, RULE_DIR, "quality_all.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "filtering_quality_all.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "graph_quality_all.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "creating_remap_catalogue_bed.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "sort_remap_bed.rules")

#include: os.path.join(BASE_DIR, RULE_DIR, "delete_trim.rules")

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
            dict_experiment_chip_filename[ experiment_name] = []
            for row in reader:

                current_filename = row[ FILENAME_HEADER].strip()
                current_experiment_type = row[ CONTROL_HEADER].strip()



                if current_experiment_type == "0":
                    dict_experiment_chip_filename[ experiment_name].append( current_filename)


# print( list_exp)
# pprint.pprint( dict_experiment_chip_filename)
#print( expand( os.path.join( QUALITY_DIR, "macs2_{experiment_name}.quality_all"), experiment_name = list_exp))


#================================================================#
#                         Workflow                               #
#================================================================#


rule all:
    input:  # expand( os.path.join( QUALITY_DIR, "macs2_{experiment_name}.quality_all"), experiment_name = list_exp),
            # os.path.join( QUALITY_DIR, "results", "macs2.quality_all")
            # os.path.join( QUALITY_DIR, "results", "macs2_passed.quality_all"),
            # os.path.join( "remap2020_unsorted.bed"),
            os.path.join( REMAP_FULLANME + ".bed"),
            os.path.join( QUALITY_DIR,  "results", "macs2.quality_all.pdf")

rule quality_all:
    input:
            expand( os.path.join( QUALITY_DIR, "macs2_{experiment_name}.quality_all"), experiment_name = list_exp)
    output:
            os.path.join( QUALITY_DIR,  "results", "macs2.quality_all")
    resources:
            res=1
    log:
            os.path.join( QUALITY_DIR, "log", "quality_all.log")
    params:
               indir = QUALITY_DIR,
               outdir = os.path.join( QUALITY_DIR, "results")
    shell:"""
        mkdir -p {params.outdir}

        cat {input} | awk "NR%2==0" | awk -F" " 'BEGIN {{print "experiment_name\tNSC\tRSC\tFRiP\tnb_peaks\tscore_NSC\tscore_RSC\tscore_FRiP\tscore_total" }}
                                                  {{
                                                    if( $2>=1.10) nsc=2; else if( $2>=1.05) nsc=1; else nsc=0;
                                                    if( $3>=1) rsc=2; else if( $3>=0.8) rsc=1; else rsc=0;
                                                    if( $4>=1) frip=1; else frip=0; score=nsc+rsc+frip; print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"nsc"\t"rsc"\t"frip"\t"score
                                                  }}' >> {output}
    """
