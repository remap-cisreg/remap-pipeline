"""
Author: Jeanne Chèneby
Affiliation: TAGC
Aim: Workflow ReMap human
Date: 01-12-16
last update: 13-09-2018
Run : snakemake --snakefile 2.scripts/snakefiles_workflow/Snakefile_remap_non_redundant.py --printshellcmds --cores 10 --cluster-config 2.scripts/cluster_configuration/cluster_conda_torque_remap_non_redundant.json --cluster "qsub -V -q {cluster.queue} -l nodes={cluster.node}:ppn={cluster.thread} -o {cluster.stdout} -e {cluster.stderr}" --keep-going --configfile 2.scripts/snakemake_configuration/Snakefile_config_remap_non_redundant.json --use-conda

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

BASE_DIR = config["working_dir"]["base"]
workdir: BASE_DIR

# Directories
PREPROCESSING_DIR = config["working_dir"]["preprocessing"]
TAB_DIR = config["working_dir"]["tab"]
PEAKCALLING_DIR = config["working_dir"]["outdir"]
QUALITY_DIR = config["working_dir"]["quality"]
BAM_DIR = config["working_dir"]["bam"]
RULE_DIR = config["working_dir"]["rule"]
MODIF_DIR = config["working_dir"]["modif"]
if not os.path.exists( MODIF_DIR):
    os.makedirs( MODIF_DIR)
if not os.path.exists( os.path.join( MODIF_DIR, "log")):
    os.makedirs( os.path.join( MODIF_DIR, "log"))

# Files
PATH_FILE_ALL_TSV = config["file"]["tsv_all_metadata"]
PATH_FILE_MODIF = config["file"]["tsv_modification"]
PATH_FILE_ALL_BED = config["file"]["bed_all"]

FILE_ALL_TSV_NAME = os.path.basename( PATH_FILE_ALL_TSV)


#================================================================#
#                         Includes                               #
#================================================================#


# include: os.path.join(BASE_DIR, RULE_DIR, "tf_allPeaks.rules")



#================================================================#
#                     Defining dataset                           #
#================================================================#

# Create a dictionary where the key is the old name and the value the new one
dict_renaming = {}

file_modif = open( PATH_FILE_MODIF, 'r')
for line in file_modif:
    split_line =  line.strip().split( "\t")
    if len( split_line) > 1:
        if split_line[ 0] not in dict_renaming:
            dict_renaming[ split_line[ 0].strip()] = split_line[ 1].strip()


file_modif.close()

# pprint.pprint( dict_renaming)



# Get all experiments

dict_exp_modif = {}
dict_bam_modif = {}

path_file_log_renaming_exp = os.path.join( MODIF_DIR, "log", "list_exp_modification.tab")
path_file_log_renaming_bam = os.path.join( MODIF_DIR, "log", "list_bam_modification.tab")
file_log_renaming_exp = open( path_file_log_renaming_exp, 'w')
file_log_renaming_bam = open( path_file_log_renaming_bam, 'w')

list_objects_indir = os.listdir( TAB_DIR) # get all files' and folders' names in the current directory

for objects_indir in list_objects_indir: # loop through all the files and folders
    if objects_indir.endswith( "_summary.tab"): # check whether the current object is a folder or not
        experiment_name = objects_indir.rsplit( "_", 1)[ 0]

        if any( current_key in experiment_name for current_key in dict_renaming):

            old_experiment_name = experiment_name
            new_experiment_name = experiment_name

            list_key_modif = []
            # dict_bam_modif_reverse = {}

            for current_key in dict_renaming:
                if experiment_name.find( current_key) >= 0:
                    new_experiment_name = new_experiment_name.replace( current_key, dict_renaming[ current_key])


                    list_key_modif.append( current_key)


            dict_exp_modif[ new_experiment_name] = {}
            dict_exp_modif[ new_experiment_name][ "old_name"] = old_experiment_name
            dict_exp_modif[ new_experiment_name][ "old_modif"] = list_key_modif


            with open( os.path.join( TAB_DIR, objects_indir), 'r') as file_tab:
                next( file_tab)
                for line in file_tab:
                    new_bam_name = line
                    for current_key in dict_renaming:
                        if experiment_name.find( current_key) >= 0:
                            new_bam_name = line.split( "\t")[ 0].replace( current_key, dict_renaming[ current_key])

                    dict_bam_modif[ new_bam_name] = line.split( "\t")[ 0]
                    file_log_renaming_bam.write( "\t".join( [ line.split( "\t")[ 0], new_bam_name]) + "\n")

                    # dict_bam_modif_reverse[ line.split( "\t")[ 0]] = new_bam_name

            # dict_exp_modif[ new_experiment_name][ "dict_bam_modif_reverse"] = dict_bam_modif_reverse


            file_log_renaming_exp.write( "\t".join( [ old_experiment_name, new_experiment_name]) + "\n")


file_log_renaming_exp.close()
file_log_renaming_bam.close()
# pprint.pprint( dict_exp_modif)
# pprint.pprint( dict_bam_modif)





#================================================================#
#                         Workflow                               #
#================================================================#

rule all:
    input:
        os.path.join( MODIF_DIR, "1.metadata", FILE_ALL_TSV_NAME),
        expand( os.path.join( MODIF_DIR, TAB_DIR, "{experiment_name}_summary.tab"), experiment_name = dict_exp_modif),
        expand( os.path.join( MODIF_DIR, BAM_DIR, "{new_bam_name}.bam"), new_bam_name = dict_bam_modif),
        # expand( os.path.join( MODIF_DIR, PREPROCESSING_DIR, "trim_fastq", "{new_bam_name}", "{new_bam_name}_del.ok") , new_bam_name = dict_bam_modif),
        # expand( os.path.join( MODIF_DIR, PEAKCALLING_DIR, "{experiment_name}", "macs2", "{experiment_name}_model.r"), experiment_name = dict_exp_modif),
        expand( os.path.join( MODIF_DIR, PEAKCALLING_DIR, "{experiment_name}", "macs2", "{experiment_name}_peaks.narrowPeak"), experiment_name = dict_exp_modif),
        expand( os.path.join( MODIF_DIR, PEAKCALLING_DIR, "{experiment_name}", "macs2", "{experiment_name}_peaks.xls"), experiment_name = dict_exp_modif),
        # expand( os.path.join( MODIF_DIR, PEAKCALLING_DIR, "{experiment_name}", "macs2", "{experiment_name}_summit.bed"), experiment_name = dict_exp_modif),
        expand( os.path.join( MODIF_DIR, PEAKCALLING_DIR, "{experiment_name}", "macs2", "log", "{experiment_name}.log"), experiment_name = dict_exp_modif),
        os.path.join( MODIF_DIR, QUALITY_DIR, "results", "macs2_passed.quality_all"),
        os.path.join( MODIF_DIR, QUALITY_DIR, "results", "macs2_failed.quality_all"),
        os.path.join( MODIF_DIR, QUALITY_DIR, "results", "macs2.quality_all"),
        os.path.join( MODIF_DIR, PATH_FILE_ALL_BED)




rule renaming_big_TSV:
    input:
            file = PATH_FILE_ALL_TSV
    output:
            file = os.path.join( MODIF_DIR, "1.metadata", FILE_ALL_TSV_NAME)
    resources:
            res=1
    run:
        file_tsv_input = open( input.file, 'r')
        file_tsv_output = open( output.file, 'w')

        # parse all line
        for line in file_tsv_input:

            # only try to replace line containing at list one value to change
            if any(current_key in line for current_key in dict_renaming):
                new_line = line.strip()

                for current_key in dict_renaming:
                    if new_line.find( current_key) >= 0:

                        new_line = new_line.replace( current_key, dict_renaming[ current_key])
                file_tsv_output.write( new_line + "\t" + "1" + "\n")



            else:
                file_tsv_output.write( line.strip() + "\t" + "0" + "\n")

        file_tsv_input.close()
        file_tsv_output.close()



rule renaming_tab:
    input:
            file = lambda wildcards : expand( os.path.join( TAB_DIR,  "{tab_filename}_summary.tab"), tab_filename = dict_exp_modif[ wildcards.experiment_name][ "old_name"])
    output:
            file = os.path.join( MODIF_DIR, TAB_DIR, "{experiment_name}_summary.tab")
    resources:
            res=1
    run:
        file_tab_input = open( str( input.file), 'r')
        file_tab_output = open( output.file, 'w')

        for line in file_tab_input:
            # only try to replace line containing at list one value to change
            if any(current_key in line for current_key in dict_renaming):
                new_line = line.strip()

                for current_key in dict_renaming:
                    if new_line.find( current_key) >= 0:

                        new_line = new_line.replace( current_key, dict_renaming[ current_key])
                file_tab_output.write( new_line + "\n")



            else:
                file_tab_output.write( line.strip() + "\n")


        file_tab_input.close()
        file_tab_output.close()








rule renaming_bam:
    input:
            lambda wildcards : expand( os.path.join( BAM_DIR, "{old_bam_name}.bam"), old_bam_name = dict_bam_modif[ wildcards.new_bam_name])
    output:
            os.path.join( MODIF_DIR, BAM_DIR, "{new_bam_name}.bam")
    resources:
            res=1
    shell:
            """cp {input} {output}"""


rule renaming_trim:
    input:
            lambda wildcards : expand( os.path.join( PREPROCESSING_DIR, "trim_fastq", "{old_bam_name}", "{old_bam_name}_del.ok"), old_bam_name = dict_bam_modif[ wildcards.new_bam_name])
    output:
            os.path.join( MODIF_DIR, PREPROCESSING_DIR, "trim_fastq", "{new_bam_name}", "{new_bam_name}_del.ok")
    resources:
            res=1
    shell:
            """cp {input} {output}"""







rule renaming_macs2_narrowPeak:
    input:
            file = lambda wildcards : expand( os.path.join( PEAKCALLING_DIR, "{old_experiment_name}", "macs2", "{old_experiment_name}_peaks.narrowPeak"), old_experiment_name = dict_exp_modif[ wildcards.experiment_name][ "old_name"])
    output:
            file = os.path.join( MODIF_DIR, PEAKCALLING_DIR, "{experiment_name}", "macs2", "{experiment_name}_peaks.narrowPeak")
    resources:
            res=1
    run:
        file_tab_input = open( str( input.file), 'r')
        file_tab_output = open( output.file, 'w')

        for line in file_tab_input:
            file_tab_output.write( line.replace( dict_exp_modif[ wildcards.experiment_name][ "old_name"], wildcards.experiment_name))

        file_tab_input.close()
        file_tab_output.close()


rule renaming_macs2_model:
    input:
            file = lambda wildcards : expand( os.path.join( PEAKCALLING_DIR, "{old_experiment_name}", "macs2", "{old_experiment_name}_model.r"), old_experiment_name = dict_exp_modif[ wildcards.experiment_name][ "old_name"])
    output:
            file = os.path.join( MODIF_DIR, PEAKCALLING_DIR, "{experiment_name}", "macs2", "{experiment_name}_model.r")
    resources:
            res=1
    run:
        file_tab_input = open( str( input.file), 'r')
        file_tab_output = open( output.file, 'w')

        for line in file_tab_input:
            file_tab_output.write( line.replace( dict_exp_modif[ wildcards.experiment_name][ "old_name"], wildcards.experiment_name))

        file_tab_input.close()
        file_tab_output.close()


rule renaming_macs2_xls:
    input:
            file = lambda wildcards : expand( os.path.join( PEAKCALLING_DIR, "{old_experiment_name}", "macs2", "{old_experiment_name}_peaks.xls"), old_experiment_name = dict_exp_modif[ wildcards.experiment_name][ "old_name"])
    output:
            file = os.path.join( MODIF_DIR, PEAKCALLING_DIR, "{experiment_name}", "macs2", "{experiment_name}_peaks.xls")
    resources:
            res=1
    run:
        file_tab_input = open( str( input.file), 'r')
        file_tab_output = open( output.file, 'w')

        for line in file_tab_input:
            file_tab_output.write( line.replace( dict_exp_modif[ wildcards.experiment_name][ "old_name"], wildcards.experiment_name))

        file_tab_input.close()
        file_tab_output.close()


rule renaming_macs2_summits:
    input:
            file = lambda wildcards : expand( os.path.join( PEAKCALLING_DIR, "{old_experiment_name}", "macs2", "{old_experiment_name}_summits.bed"), old_experiment_name = dict_exp_modif[ wildcards.experiment_name][ "old_name"])
    output:
            file = os.path.join( MODIF_DIR, PEAKCALLING_DIR, "{experiment_name}", "macs2", "{experiment_name}_summits.bed")
    resources:
            res=1
    run:
        file_tab_input = open( str( input.file), 'r')
        file_tab_output = open( output.file, 'w')

        for line in file_tab_input:
            file_tab_output.write( line.replace( dict_exp_modif[ wildcards.experiment_name][ "old_name"], wildcards.experiment_name))

        file_tab_input.close()
        file_tab_output.close()


rule renaming_macs2_log:
    input:
            file = lambda wildcards : expand( os.path.join( PEAKCALLING_DIR, "{old_experiment_name}", "macs2", "log","{old_experiment_name}.log"), old_experiment_name = dict_exp_modif[ wildcards.experiment_name][ "old_name"])
    output:
            file = os.path.join( MODIF_DIR, PEAKCALLING_DIR, "{experiment_name}", "macs2", "log", "{experiment_name}.log")
    resources:
            res=1
    run:
        file_tab_input = open( str( input.file), 'r')
        file_tab_output = open( output.file, 'w')

        for line in file_tab_input:
            file_tab_output.write( line.replace( dict_exp_modif[ wildcards.experiment_name][ "old_name"], wildcards.experiment_name))

        file_tab_input.close()
        file_tab_output.close()






rule renaming_macs2_quality_all_passed:
    input:
            file = os.path.join( QUALITY_DIR, "results", "macs2_passed.quality_all")
    output:
            file = os.path.join( MODIF_DIR, QUALITY_DIR, "results", "macs2_passed.quality_all")
    resources:
            res=1
    run:
        file_input = open( str( input.file), 'r')
        file_output = open( output.file, 'w')

        # parse all line
        for line in file_input:

            # only try to replace line containing at list one value to change
            if any(current_key in line for current_key in dict_renaming):
                new_line = line.strip()

                for current_key in dict_renaming:
                    if new_line.find( current_key) >= 0:

                        new_line = new_line.replace( current_key, dict_renaming[ current_key])
                file_output.write( new_line + "\n")



            else:
                file_output.write( line.strip() + "\n")

        file_input.close()
        file_output.close()


rule renaming_macs2_quality_all_failed:
    input:
            file = os.path.join( QUALITY_DIR, "results", "macs2_failed.quality_all")
    output:
            file = os.path.join( MODIF_DIR, QUALITY_DIR, "results", "macs2_failed.quality_all")
    resources:
            res=1
    run:
        file_input = open( str( input.file), 'r')
        file_output = open( output.file, 'w')

        # parse all line
        for line in file_input:

            # only try to replace line containing at list one value to change
            if any(current_key in line for current_key in dict_renaming):
                new_line = line.strip()

                for current_key in dict_renaming:
                    if new_line.find( current_key) >= 0:

                        new_line = new_line.replace( current_key, dict_renaming[ current_key])
                file_output.write( new_line + "\n")



            else:
                file_output.write( line.strip() + "\n")

        file_input.close()
        file_output.close()



rule renaming_macs2_quality_all:
    input:
            file = os.path.join( QUALITY_DIR, "results", "macs2.quality_all")
    output:
            file = os.path.join( MODIF_DIR, QUALITY_DIR, "results", "macs2.quality_all")
    resources:
            res=1
    run:
        file_input = open( str( input.file), 'r')
        file_output = open( output.file, 'w')

        # parse all line
        for line in file_input:

            # only try to replace line containing at list one value to change
            if any(current_key in line for current_key in dict_renaming):
                new_line = line.strip()

                for current_key in dict_renaming:
                    if new_line.find( current_key) >= 0:

                        new_line = new_line.replace( current_key, dict_renaming[ current_key])
                file_output.write( new_line + "\n")



            else:
                file_output.write( line.strip() + "\n")

        file_input.close()
        file_output.close()



rule renaming_macs2_bed_all:
    input:
            file = PATH_FILE_ALL_BED
    output:
            file = os.path.join( MODIF_DIR, PATH_FILE_ALL_BED)
    resources:
            res=1
    run:
        file_input = open( str( input.file), 'r')
        file_output = open( output.file, 'w')

        # parse all line
        for line in file_input:

            # only try to replace line containing at list one value to change
            if any(current_key in line for current_key in dict_renaming):
                new_line = line.strip()

                for current_key in dict_renaming:
                    if new_line.find( current_key) >= 0:

                        new_line = new_line.replace( current_key, dict_renaming[ current_key])
                file_output.write( new_line + "\n")



            else:
                file_output.write( line.strip() + "\n")

        file_input.close()
        file_output.close()
