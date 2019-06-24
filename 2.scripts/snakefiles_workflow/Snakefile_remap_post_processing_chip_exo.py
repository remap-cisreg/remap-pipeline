"""
Author: Jeanne Chèneby
Affiliation: TAGC
Aim: Workflow ReMap human
Date: 01-12-16
last update: 13-09-2018
Run:snakemake --snakefile 2.scripts/snakefiles_workflow/Snakefile_remap_chip_exo.py --printshellcmds --cores 10 --cluster-config .scripts/cluster_configuration/cluster_conda_toque.json --cluster "qsub -V -q {cluster.queue} -l nodes={cluster.node}:ppn={cluster.thread} -o {cluster.stdout} -e {cluster.stderr}" --keep-going --configfile 2.snakemake_configuration/Snakefile_config_remap_torque_saccapus_chip_exo.json --use-conda

dag: snakemake --snakefile .scripts/snakefiles_workflow/Snakefile_remap_chip_exo.py --printshellcmds --cores 10 --cluster-config config/sacapus.json --cluster "qsub -V -q {cluster.queue} -l nodes={cluster.node}:ppn={cluster.thread} -o {cluster.stdout} -e {cluster.stderr}" --keep-going --configfile config/Snakefile_config_remap_saccapus.json --dag 2> /dev/null | dot -T svg > dag.svg
rulegraph: snakemake --snakefile.scripts/snakefiles_workflow/Snakefile_remap_chip_exo.py --printshellcmds --cores 10 --cluster-config config/sacapus.json --cluster "qsub -V -q {cluster.queue} -l nodes={cluster.node}:ppn={cluster.thread} -o {cluster.stdout} -e {cluster.stderr}" --keep-going --configfile config/Snakefile_config_remap_saccapus.json --rulegraph 2> /dev/null | dot -T svg > rulegraph.svg


Run (meso): snakemake --snakefile Snakefile_remap_chip_exo.py --printshellcmds --cores 99 --cluster-config config/mesocentre.json --cluster "srun -p {cluster.partition} -N {cluster.node} -n {cluster.thread} -o {cluster.stdout} -e {cluster.stderr}" --keep-going --configfile config/Snakefile_config_remap.json


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

PEAKCALLER = config["info"]["peakcaller"]
EXTENSION_PEAK_INPUT = config["extention"]["peak_input"]
EXTENSION_PEAK_OUTPUT = config["extention"]["peak_output"]


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

include: os.path.join(BASE_DIR, RULE_DIR, "quality_chip_exo_intersect_reads_bed.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "quality_chip_exo_nb_reads_bam.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "quality_chip_exo_filtering_peaks.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "quality_chip_exo_filtering_peaks_getting_kept_blanced_peaks.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "quality_all_chip_exo_nb_peak.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "creating_remap_catalogue_bed.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "sort_remap_bed.rules")



#================================================================#
#                     Defining dataset                           #
#================================================================#


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



#================================================================#
#                         Workflow                               #
#================================================================#


rule all:
    # input:  expand( os.path.join( QUALITY_DIR, "macs2_{experiment_name}.FRiP_NSC_RSC"), experiment_name = list_exp),
    input:
            # expand( os.path.join( PEAKCALLING_DIR, "{experiment_name}", "macs2", "{experiment_name}_peaks.narrowPeak"), experiment_name = list_exp),
            # expand( os.path.join( PEAKCALLING_DIR, "{experiment_name}", "macs2", "{experiment_name}_filtered_peaks.narrowPeak" ), experiment_name = list_exp),
            REMAP_FULLANME + ".bed"
