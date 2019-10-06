"""
Author: Jeanne Chèneby
Affiliation: TAGC
Aim: Workflow renaming files after runing core
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
PATH_FILE_LIST_EXP_ALL_BED = config["file"]["list_exp_bed_all"]

FILE_ALL_TSV_NAME = os.path.basename( PATH_FILE_ALL_TSV)


#================================================================#
#                         Includes                               #
#================================================================#


# include: os.path.join(BASE_DIR, RULE_DIR, "tf_allPeaks.rules")



#================================================================#
#                     Defining dataset                           #
#================================================================#

# Create a dictionary where the key is the old name and the value the new one
dict_renaming_tf = {}
dict_renaming_cell = {}

file_modif = open( PATH_FILE_MODIF, 'r')
for line in file_modif:
    split_line =  line.strip().split( "\t")
    if len( split_line) > 1:
        if split_line[ 2] == "tf":
            if split_line[ 0] not in dict_renaming_tf:
                dict_renaming_tf[ split_line[ 0].strip()] = split_line[ 1].strip()
        elif split_line[ 2] == "cell":
            if split_line[ 0] not in dict_renaming_cell:
                dict_renaming_cell[ split_line[ 0].strip()] = split_line[ 1].strip()


file_modif.close()
# Get all experiments

dict_exp_modif = {}
dict_bam_modif = {}
dict_bam_modif_reverse = {}

path_file_log_renaming_exp = os.path.join( MODIF_DIR, "log", "list_exp_modification.tab")
path_file_log_renaming_bam = os.path.join( MODIF_DIR, "log", "list_bam_modification.tab")
file_log_renaming_exp = open( path_file_log_renaming_exp, 'w')
file_log_renaming_bam = open( path_file_log_renaming_bam, 'w')

list_objects_indir = os.listdir( TAB_DIR) # get all files' and folders' names in the current directory
for objects_indir in list_objects_indir: # loop through all the files and folders
    if objects_indir.endswith( "_summary.tab"): # check whether the current object is a folder or not
        experiment_name = objects_indir.rsplit( "_", 1)[ 0]
        new_experiment_name = ""

        old_tf = experiment_name.split( ".", 2)[ 1]
        old_cell = experiment_name.split( ".", 2)[ 2]

        # Dealing with renaming tf
        for current_key in dict_renaming_tf:
            current_regex = "^" + current_key + "$"
            current_regex_modify = "^" + current_key + "_+"

            if re.match( current_regex, old_tf):
                old_tf = dict_renaming_tf[ current_key]

            elif re.match( current_regex_modify, old_tf):
                old_tf = re.sub( current_regex_modify, dict_renaming_tf[ current_key] + "_", old_tf)


        # Dealing with renaming cell
        for current_key in dict_renaming_cell:
            current_regex_pure = "^" + current_key + "$"
            current_regex_modify = "^" + current_key + "_+"

            if re.match( current_regex_pure, old_cell):
                old_cell = dict_renaming_cell[ current_key]


            elif re.match( current_regex_modify, old_cell):
                old_cell = re.sub( current_regex_modify, dict_renaming_cell[ current_key] + "_", old_cell)

        new_experiment_name = ".".join( [ experiment_name.split( ".", 1)[ 0], old_tf, old_cell])

        if new_experiment_name != experiment_name:
            dict_exp_modif[ new_experiment_name] = {}
            dict_exp_modif[ new_experiment_name][ "old_name"] = experiment_name

            file_log_renaming_exp.write( "\t".join( [ experiment_name, new_experiment_name]) + "\n")
            with open( os.path.join( TAB_DIR, objects_indir), 'r') as file_tab:
                next( file_tab)
                for line in file_tab:
                    split_line = line.split( "\t")
                    old_bam_name = split_line[ 0]
                    new_bam_name = ""
                    if split_line[ 3] == "0":
                        if len( old_bam_name.split( ".")) == 4:
                            new_bam_name = ".".join( [ old_bam_name.split( ".", 1)[ 0], old_tf, old_cell, old_bam_name.rsplit( ".", 1)[ -1]])
                        else:
                            new_bam_name = ".".join( [ old_bam_name.split( ".", 1)[ 0], old_tf, old_cell])
                    elif split_line[ 3] == "1":
                        if len( old_bam_name.split( ".")) == 4:
                            new_bam_name = ".".join( old_bam_name.split( ".", 2)[ 0:2] + [ old_cell, old_bam_name.rsplit( ".", 1)[ -1]])
                        else:
                            new_bam_name = ".".join( old_bam_name.split( ".", 2)[ 0:2] + [ old_cell])
                    # print( old_bam_name)
                    # print( new_bam_name)
                    dict_bam_modif[ new_bam_name] = old_bam_name
                    dict_bam_modif_reverse[ old_bam_name] = new_bam_name
                    file_log_renaming_bam.write( "\t".join( [ old_bam_name, new_bam_name]) + "\n")

file_log_renaming_exp.close()
file_log_renaming_bam.close()


file_list_exp_all_bed = open( PATH_FILE_LIST_EXP_ALL_BED, 'r')
dict_exp_modif_reverse = {}
dict_exp_modif_plus = {}

file_list_exp_all_modif = open( os.path.join( MODIF_DIR, "log", "list_exp_modif_from_file.txt"), 'w')


for line in file_list_exp_all_bed:
    experiment_name = line.strip()
    new_experiment_name = ""

    old_tf = experiment_name.split( ".", 2)[ 1]
    old_cell = experiment_name.split( ".", 2)[ 2]

    # Dealing with renaming tf
    for current_key in dict_renaming_tf:
        current_regex = "^" + current_key + "$"
        current_regex_modify = "^" + current_key + "_+"

        if re.match( current_regex, old_tf):
            old_tf = dict_renaming_tf[ current_key]

        elif re.match( current_regex_modify, old_tf):
            old_tf = re.sub( current_regex_modify, dict_renaming_tf[ current_key] + "_", old_tf)


    # Dealing with renaming cell
    for current_key in dict_renaming_cell:
        current_regex_pure = "^" + current_key + "$"
        current_regex_modify = "^" + current_key + "_+"

        if re.match( current_regex_pure, old_cell):
            old_cell = dict_renaming_cell[ current_key]

        elif re.match( current_regex_modify, old_cell):
            old_cell = re.sub( current_regex_modify, dict_renaming_cell[ current_key] + "_", old_cell)

    new_experiment_name = ".".join( [ experiment_name.split( ".", 1)[ 0], old_tf, old_cell])

    if new_experiment_name != experiment_name:

        dict_exp_modif_reverse[ experiment_name.strip()] = new_experiment_name.strip()
        dict_exp_modif_plus[ new_experiment_name.strip()] = experiment_name.strip()

        file_list_exp_all_modif.write( experiment_name.strip() + "\t" + new_experiment_name.strip() + "\n")

file_list_exp_all_bed.close()
file_list_exp_all_modif.close()


#================================================================#
#                         Workflow                               #
#================================================================#

rule all:
    input:
        os.path.join( MODIF_DIR, "1.metadata", FILE_ALL_TSV_NAME),
        expand( os.path.join( MODIF_DIR, TAB_DIR, "{experiment_name}_summary.tab"), experiment_name = dict_exp_modif),
        expand( os.path.join( MODIF_DIR, BAM_DIR, "{new_bam_name}.bam"), new_bam_name = dict_bam_modif),
        expand( os.path.join( MODIF_DIR, PREPROCESSING_DIR, "trim_fastq", "{new_bam_name}", "{new_bam_name}_del.ok") , new_bam_name = dict_bam_modif),
        expand( os.path.join( MODIF_DIR, PEAKCALLING_DIR, "{experiment_name}", "macs2", "{experiment_name}_model.r"), experiment_name = dict_exp_modif),
        expand( os.path.join( MODIF_DIR, PEAKCALLING_DIR, "{experiment_name}", "macs2", "{experiment_name}_peaks.narrowPeak"), experiment_name = dict_exp_modif_plus),
        expand( os.path.join( MODIF_DIR, PEAKCALLING_DIR, "{experiment_name}", "macs2", "{experiment_name}_peaks.xls"), experiment_name = dict_exp_modif_plus),
        expand( os.path.join( MODIF_DIR, PEAKCALLING_DIR, "{experiment_name}", "macs2", "{experiment_name}_summit.bed"), experiment_name = dict_exp_modif),
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
    params:
        column_tf = int( config["info"]["tsv_all_metadata_tf_column"]),
        column_cell = int( config["info"]["tsv_all_metadata_cell_column"])
    run:
        file_tsv_input = open( input.file, 'r')
        file_tsv_output = open( output.file, 'w')

        # parse all line
        for line in file_tsv_input:

            modify = False
            split_line = line.split( "\t")
            old_tf = split_line[ params.column_tf]
            old_cell = split_line[ params.column_cell]

            # Dealing with TF
            for current_key in dict_renaming_tf:

                current_regex_pure = "^" + current_key + "$"
                current_regex_modify = "^" + current_key + "_"

                if re.match( current_regex_pure, old_tf):
                    old_tf = re.sub( current_regex_pure, "\t" + dict_renaming_tf[ current_key], old_tf)
                    modify = True

                elif re.match( current_regex_modify, old_tf):
                    old_tf = re.sub( current_regex_modify, "\t" + dict_renaming_tf[ current_key] + "_", old_tf)
                    modify = True

            # Dealing with renaming cell
            for current_key in dict_renaming_cell:

                current_regex_pure = "^" + current_key + "$"
                current_regex_modify = "^" + current_key + "_"

                if re.match( current_regex_pure, old_cell):
                    old_cell = re.sub( current_regex_pure, dict_renaming_cell[ current_key], old_cell)
                    modify = True

                elif re.match( current_regex_modify, old_cell):
                    old_cell = re.sub( current_regex_modify, dict_renaming_cell[ current_key] + "_", old_cell)
                    modify = True

            if modify:
                new_line = "\t".join( split_line[ 0:params.column_tf] + [ old_tf] + [ old_cell] + split_line[ params.column_cell + 1:-1])
                file_tsv_output.write( new_line.strip() + "\t" + "1\n")
            else:
                file_tsv_output.write( line.strip() + "\t" + "0\n")

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

            split_line = line.split( "\t")
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
            file = lambda wildcards : expand( os.path.join( PEAKCALLING_DIR, "{old_experiment_name}", "macs2", "{old_experiment_name}_peaks.narrowPeak"), old_experiment_name = dict_exp_modif_plus[ wildcards.experiment_name])
    output:
            file = os.path.join( MODIF_DIR, PEAKCALLING_DIR, "{experiment_name}", "macs2", "{experiment_name}_peaks.narrowPeak")
    resources:
            res=1
    run:
        file_tab_input = open( str( input.file), 'r')
        file_tab_output = open( output.file, 'w')

        for line in file_tab_input:
            file_tab_output.write( line.replace( dict_exp_modif_plus[ wildcards.experiment_name], wildcards.experiment_name))

        file_tab_input.close()
        file_tab_output.close()


rule renaming_macs2_model:
    input:
            file = lambda wildcards : expand( os.path.join( PEAKCALLING_DIR, "{old_experiment_name}", "macs2", "{old_experiment_name}_model.r"), old_experiment_name = dict_exp_modif_plus[ wildcards.experiment_name])
    output:
            file = os.path.join( MODIF_DIR, PEAKCALLING_DIR, "{experiment_name}", "macs2", "{experiment_name}_model.r")
    resources:
            res=1
    run:
        file_tab_input = open( str( input.file), 'r')
        file_tab_output = open( output.file, 'w')

        for line in file_tab_input:
            file_tab_output.write( line.replace( dict_exp_modif_plus[ wildcards.experiment_name], wildcards.experiment_name))

        file_tab_input.close()
        file_tab_output.close()


rule renaming_macs2_xls:
    input:
            file = lambda wildcards : expand( os.path.join( PEAKCALLING_DIR, "{old_experiment_name}", "macs2", "{old_experiment_name}_peaks.xls"), old_experiment_name = dict_exp_modif_plus[ wildcards.experiment_name])
    output:
            file = os.path.join( MODIF_DIR, PEAKCALLING_DIR, "{experiment_name}", "macs2", "{experiment_name}_peaks.xls")
    resources:
            res=1
    run:
        file_tab_input = open( str( input.file), 'r')
        file_tab_output = open( output.file, 'w')

        for line in file_tab_input:
            file_tab_output.write( line.replace( dict_exp_modif_plus[ wildcards.experiment_name], wildcards.experiment_name))

        file_tab_input.close()
        file_tab_output.close()


rule renaming_macs2_summits:
    input:
            file = lambda wildcards : expand( os.path.join( PEAKCALLING_DIR, "{old_experiment_name}", "macs2", "{old_experiment_name}_summits.bed"), old_experiment_name = dict_exp_modif_plus[ wildcards.experiment_name])
    output:
            file = os.path.join( MODIF_DIR, PEAKCALLING_DIR, "{experiment_name}", "macs2", "{experiment_name}_summits.bed")
    resources:
            res=1
    run:
        file_tab_input = open( str( input.file), 'r')
        file_tab_output = open( output.file, 'w')

        for line in file_tab_input:
            file_tab_output.write( line.replace( dict_exp_modif_plus[ wildcards.experiment_name], wildcards.experiment_name))

        file_tab_input.close()
        file_tab_output.close()


rule renaming_macs2_log:
    input:
            file = lambda wildcards : expand( os.path.join( PEAKCALLING_DIR, "{old_experiment_name}", "macs2", "log","{old_experiment_name}.log"), old_experiment_name = dict_exp_modif_plus[ wildcards.experiment_name])
    output:
            file = os.path.join( MODIF_DIR, PEAKCALLING_DIR, "{experiment_name}", "macs2", "log", "{experiment_name}.log")
    resources:
            res=1
    run:
        file_tab_input = open( str( input.file), 'r')
        file_tab_output = open( output.file, 'w')

        for line in file_tab_input:
            file_tab_output.write( line.replace( dict_exp_modif_plus[ wildcards.experiment_name], wildcards.experiment_name))

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
            if any(current_key in line for current_key in dict_exp_modif_reverse):
                new_line = line.strip()

                for current_key in dict_exp_modif_reverse:
                    if new_line.find( current_key) >= 0:

                        new_line = new_line.replace( current_key, dict_exp_modif_reverse[ current_key])
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
            if any(current_key in line for current_key in dict_exp_modif_reverse):
                new_line = line.strip()

                for current_key in dict_exp_modif_reverse:
                    if new_line.find( current_key) >= 0:

                        new_line = new_line.replace( current_key, dict_exp_modif_reverse[ current_key])
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
            if any(current_key in line for current_key in dict_exp_modif_reverse):
                new_line = line.strip()

                for current_key in dict_exp_modif_reverse:
                    if new_line.find( current_key) >= 0:

                        new_line = new_line.replace( current_key, dict_exp_modif_reverse[ current_key])
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
    log:
        os.path.join( MODIF_DIR, "log", "list_exp_modification_all_bed_supl.txt")
    run:
        file_input = open( str( input.file), 'r')
        file_output = open( output.file, 'w')
        file_log = open( log[ 0], 'w')

        # parse all line if experiment in tab
        for line in file_input:
            split_line = line.split( "\t")
            experiment_name = split_line[ 3].strip()


            if experiment_name in dict_exp_modif_reverse:
                new_experiment_name = dict_exp_modif_reverse[ experiment_name]
                file_output.write( "\t".join( split_line[ 0:3] + [ new_experiment_name] + split_line[ 4:]))
                file_log.write( experiment_name + "\t" + new_experiment_name + "\n")
            # if experiement not in tab
            else:
                file_output.write( line)

        file_input.close()
        file_output.close()
        file_log.close()
