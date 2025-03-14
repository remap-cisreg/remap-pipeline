#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test initialisation fichier snakemake

Created on Thu Mar  6 09:08:58 2025

@author: souville
"""
import os
# BASE_DIR = config["working_dir"]["base_dir"]
# workdir: BASE_DIR

# # Directories
# TAB_DIR = config["working_dir"]["tab_dir"]
# PREPROCESSING_DIR = config["working_dir"]["preprocessing"]
# PEAKCALLING_DIR = config["working_dir"]["outdir"]
# QUALITY_DIR = config["working_dir"]["quality"]
# BAM_DIR = config["working_dir"]["bam_dir"]
# RULE_DIR = config["working_dir"]["rule_dir"]

# # Header CSV
# CONTROL_HEADER = config["header"]["control"]
# FILENAME_HEADER = config["header"]["filename"]
# LIBRARY_HEADER = config["header"]["library"]
# LIBRARY_URL = config["header"]["url"]
# LIBRARY_MD5 = config["header"]["md5"]

# # Extension peakcaller
# EXTENSION_PEAK = config[ "extention"][ "peak"]

BASE_DIR = "/run/user/1001/gvfs/sftp:host=core.cluster.france-bioinformatique.fr,user=souvillel/shared/projects/remap/remap-pipeline/2025/Minimum_slurm_fonctionnel"
workdir: BASE_DIR

# Directories
TAB_DIR = "1.metadata/tab"
PREPROCESSING_DIR = "4.preprocessing"
PEAKCALLING_DIR = "6.peakcalling"
QUALITY_DIR = "8.quality"
BAM_DIR = "5.bam"
RULE_DIR = config["working_dir"]["rule_dir"]

# Header CSV
CONTROL_HEADER = "isControl"
FILENAME_HEADER = "filename"
LIBRARY_HEADER = "library_type"
LIBRARY_URL = "url"
LIBRARY_MD5 = "md5"

# Extension peakcaller
EXTENSION_PEAK = "_peaks.narrowPeak"


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
    if objects_indir.endswith( "_summary.tab"): # check whether the current object is a folder or not
        print(objects_indir)
        #print(BASE_DIR)
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

                #print(row[ LIBRARY_HEADER])
                # If replicat is single
                if row[ LIBRARY_HEADER] == 'SINGLE':
                    # dict_exp_info[ 'library_type'] = 'SE'
                    dict_replicat_trim_filename[ current_filename][ 'run_type'] = 'SE'<



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
        dict_experiment_chip_filename_alt[ experiment_name][ 'chip'] = [ experiment_name]
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

#print(dict_experiment_chip_filename.values())
print("alt one \n")
print(dict_experiment_chip_filename_alt.values())
print("fastq \n")
print(dict_fastq_info)

