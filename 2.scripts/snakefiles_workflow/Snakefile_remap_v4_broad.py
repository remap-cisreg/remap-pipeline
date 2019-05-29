"""
Author: Jeanne Chèneby
Affiliation: TAGC
Aim: Workflow ReMap human
Date: 01-12-16
last update: 13-09-2018
Run:snakemake --snakefile Snakefile_remap_v4.py --printshellcmds --cores 10 --cluster-config config/sacapus.json --cluster "qsub -V -q {cluster.queue} -l nodes={cluster.node}:ppn={cluster.thread} -o {cluster.stdout} -e {cluster.stderr}" --keep-going --configfile config/Snakefile_config_remap_saccapus.json --use-conda

dag: snakemake --snakefile Snakefile_remap_v4.py --printshellcmds --cores 10 --cluster-config config/sacapus.json --cluster "qsub -V -q {cluster.queue} -l nodes={cluster.node}:ppn={cluster.thread} -o {cluster.stdout} -e {cluster.stderr}" --keep-going --configfile config/Snakefile_config_remap_saccapus.json --dag 2> /dev/null | dot -T svg > dag.svg
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

include: os.path.join(BASE_DIR, RULE_DIR, "aria2.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "trim_galore.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "bowtie2.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "sam_to_bam.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "remove_mismatch.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "samtools_sort.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "samtools_remove_pcr_duplicate.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "macs2_broad.rules")
#include: os.path.join(BASE_DIR, RULE_DIR, "merge_bam.rules")
#include: os.path.join(BASE_DIR, RULE_DIR, "quality_NSC_RSC.rules")
#include: os.path.join(BASE_DIR, RULE_DIR, "quality_FRiP.rules")
#include: os.path.join(BASE_DIR, RULE_DIR, "quality_all.rules")

include: os.path.join(BASE_DIR, RULE_DIR, "delete_trim.rules")

#================================================================#
#                     Defining dataset                           #
#================================================================#


# Get all experiments
list_objects_indir = os.listdir( TAB_DIR) # get all files' and folders' names in the current directory
list_exp = []
list_replicats = []
list_forward_trim = []
list_reverse_trim = []
list_paired_trim_file = []

dict_experiment_chip_filename = {}
dict_experiment_chip_filename_alt = {}
dict_experiment_chip_filename_alt[ 'all'] = {}
dict_replicat_trim_filename = {}
dict_fastq_info = {}
dict_trim_fastq_filename = {}
dict_fastq_info = {}

for objects_indir in list_objects_indir: # loop through all the files and folders
    # if os.path.isfile( os.path.join( TAB_DIR, objects_indir)): # check whether the current object is a folder or not
    if objects_indir.endswith( "_summary.tab"): # check whether the current object is a folder or not

        experiment_name = objects_indir.rsplit( "_", 1)[ 0]
        list_exp.append( experiment_name)
        dict_experiment_chip_filename[ experiment_name] = {}

        dict_type = {}
        dict_type[ 'chip'] = []
        dict_type[ 'control'] = []

        dict_all_files = {}



        with open( os.path.join( TAB_DIR, objects_indir)) as csvfile:


            reader = csv.DictReader( csvfile, delimiter = "\t")
            for row in reader:

                current_filename = row[ FILENAME_HEADER].strip()
                current_url = row[ LIBRARY_URL].strip()
                current_md5 = row[ LIBRARY_MD5].strip()
                dict_exp_info = {}
                list_replicats.append( current_filename)


                dict_replicat_trim_filename[ current_filename] = {}
                # populating dictionary associating replicat and fastq file name


                # If replicat is single
                if row[ LIBRARY_HEADER] == 'SINGLE':
                    # dict_exp_info[ 'library_type'] = 'SE'
                    dict_replicat_trim_filename[ current_filename][ 'run_type'] = 'SE'



                    # Dealing with trim-galore annoying waus to name output
                    trim_filename = current_filename + "_trimmed"

                    dict_replicat_trim_filename[ current_filename][ 'trim_filename_forward'] = [ trim_filename]
                    dict_replicat_trim_filename[ current_filename][ 'trim_filename_reverse'] = [ ]
                    dict_replicat_trim_filename[ current_filename][ 'trim_filename_list'] = [ trim_filename]


                    dict_replicat_trim_filename[ current_filename][ 'trim_basename'] = [ current_filename, "ok"]
                    dict_replicat_trim_filename[ current_filename][ 'fastq_filename'] = [ current_filename]
                    dict_replicat_trim_filename[ current_filename][ 'fastq_filename_forward'] = [ current_filename]
                    dict_replicat_trim_filename[ current_filename][ 'fastq_filename_reverse'] = []

                    dict_fastq_info[ current_filename] = {}
                    dict_fastq_info[ current_filename][ 'fastq_filename'] = [ current_filename]
                    dict_fastq_info[ current_filename][ 'url'] = [ current_url]
                    dict_fastq_info[ current_filename][ 'md5'] = [ current_md5]
                    dict_fastq_info[ current_filename][ 'experimentName'] = [ objects_indir]

                    dict_trim_fastq_filename[ trim_filename] = {}
                    dict_trim_fastq_filename[ trim_filename][ 'fastq_filename'] = current_filename
                    dict_trim_fastq_filename[ trim_filename][ 'run_type'] = 'SE'


                # If replicat is paired
                elif row[ LIBRARY_HEADER] == 'PAIRED':
                    # dict_exp_info[ 'library_type'] = 'PE'
                    dict_replicat_trim_filename[ current_filename][ 'run_type'] = 'PE'

                    list_current_filename = current_filename.split( ".", 1)
                    list_ID_fastq_files = list_current_filename[ 0].split( "-")

                    forward_filename = ".".join( [ list_ID_fastq_files[ 0], list_current_filename[ 1]])
                    reverse_filename = ".".join( [ list_ID_fastq_files[ 1], list_current_filename[ 1]])

                    forward_url, reverse_url = current_url.split( ";")
                    forward_md5, reverse_md5 = current_md5.split( ";")
                    forward_trim_filename = forward_filename + '_val_1'
                    reverse_trim_filename = reverse_filename + '_val_2'

                    list_forward_trim.append( forward_trim_filename)
                    list_reverse_trim.append( reverse_trim_filename)
                    list_paired_trim_file.append( current_filename)


                    # FOrmating needed info for aria2
                    dict_fastq_info[ forward_filename] = {}
                    dict_fastq_info[ forward_filename][ 'url'] = [ forward_url]
                    dict_fastq_info[ forward_filename][ 'md5'] = [ forward_md5]
                    dict_fastq_info[ forward_filename][ 'experimentName'] = [ objects_indir]

                    dict_fastq_info[ reverse_filename] = {}
                    dict_fastq_info[ reverse_filename][ 'url'] = [ reverse_url]
                    dict_fastq_info[ reverse_filename][ 'md5'] = [ reverse_md5]
                    dict_fastq_info[ reverse_filename][ 'experimentName'] = [ objects_indir]


                    # Dealing with trim-galore annoying waus to name output

                    # list_trim_filename = []
                    # list_trim_filename = [ forward_trim_filename, reverse_trim_filename]
                    # list_trim_filename = [ forward_trim_filename]
                    dict_replicat_trim_filename[ current_filename][ 'trim_filename_forward'] = [ forward_trim_filename]
                    dict_replicat_trim_filename[ current_filename][ 'trim_filename_reverse'] = [ reverse_trim_filename]
                    dict_replicat_trim_filename[ current_filename][ 'trim_filename_list'] = [ forward_trim_filename, reverse_trim_filename]

                    dict_trim_fastq_filename[ forward_trim_filename] = {}
                    dict_trim_fastq_filename[ forward_trim_filename][ 'fastq_filename'] = forward_filename
                    dict_trim_fastq_filename[ forward_trim_filename][ 'run_type'] = 'PE'

                    dict_trim_fastq_filename[ reverse_trim_filename] = {}
                    dict_trim_fastq_filename[ reverse_trim_filename][ 'run_type'] = 'PE'


                    dict_replicat_trim_filename[ current_filename][ 'trim_basename'] = [ forward_filename, reverse_filename]
                    dict_replicat_trim_filename[ current_filename][ 'fastq_filename'] = [ forward_filename, reverse_filename]
                    dict_replicat_trim_filename[ current_filename][ 'fastq_filename_forward'] = [ forward_filename]
                    dict_replicat_trim_filename[ current_filename][ 'fastq_filename_reverse'] = [ reverse_filename]





                else:
                    raise ValueError(  'abnormal library type value in ', objects_indir, current_filename)



                dict_all_files[ current_filename] = dict_exp_info





                # Is file a control
                if row[ CONTROL_HEADER] == "0":
                    dict_type[ 'chip'].append( current_filename)

                elif row[ CONTROL_HEADER] == "1":
                    dict_type[ 'control'].append( current_filename)

                else:
                    raise ValueError(  'abnormal control type value in ', objects_indir, current_filename)





        dict_type[ 'all'] = list( dict_all_files.keys())

        dict_experiment_chip_filename[ experiment_name] = dict_type

        dict_experiment_chip_filename_alt[ experiment_name] = {}
        dict_experiment_chip_filename_alt[ experiment_name][ 'chip'] = {}
        dict_experiment_chip_filename_alt[ experiment_name][ 'chip'] = [ experiment_name]# = dict_type[ 'chip']
        dict_experiment_chip_filename_alt[ 'all'][ experiment_name] = dict_type[ 'chip']

        dict_experiment_chip_filename_alt[ experiment_name][ 'control'] = [ ]

        if len( dict_type[ 'control']) > 0:
            dict_experiment_chip_filename_alt[ experiment_name][ 'control'] = {}

            ID_input = sorted( dict_type[ 'control'])[ 0].split( ".", 1)[ 0]
            biosample_input = ".".join( dict_type[ 'control'][ 0].split( ".")[ 1:3])
            uniq_input_name_temp = ".".join( [ ID_input, biosample_input])

            if uniq_input_name_temp in dict_experiment_chip_filename_alt[ experiment_name][ 'control']:
                for i in range( len( dict_type[ 'control'])):

                    ID_input = sorted( dict_type[ 'control'])[ i].split( ".", 1)[ 0]
                    biosample_input = ".".join( dict_type[ 'control'][ 0].split( ".")[ 1:3])
                    if sorted( dict_type[ 'control']) != sorted( dict_experiment_chip_filename_alt[ experiment_name][ 'control'][ uniq_input_name_temp]):
                        break
                dict_experiment_chip_filename_alt[ experiment_name][ 'control'][ uniq_input_name_temp] = dict_type[ 'control']
                dict_experiment_chip_filename_alt[ 'all'][ uniq_input_name_temp] = dict_type[ 'control']

            else:
                dict_experiment_chip_filename_alt[ experiment_name][ 'control'] = [ uniq_input_name_temp]
                dict_experiment_chip_filename_alt[ 'all'][ uniq_input_name_temp] = dict_type[ 'control']




#================================================================#
#                         Workflow                               #
#================================================================#


rule all:
    # input:  expand( os.path.join( QUALITY_DIR, "macs2_{experiment_name}.FRiP_NSC_RSC"), experiment_name = list_exp),
	input:	expand( os.path.join( PEAKCALLING_DIR, "{experiment_name}", "macs2", "{experiment_name}_peaks.broadPeak"), experiment_name = list_exp),
            expand( os.path.join( PREPROCESSING_DIR, "trim_fastq", "{replicat_name_paired}", "{replicat_name_paired}_del.ok"), zip, replicat_name_paired = list_paired_trim_file)
