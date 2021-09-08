"""
Author: Jeanne Chèneby
Affiliation: TAGC
Aim: Workflow ReMap post-processing ChIP-seq and DAP-seq
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
PREPROCESSING_DIR = config["working_dir"]["preprocessing"]
BED_DIR = config["working_dir"]["bed_dir"]
FASTA_DIR = config["working_dir"]["fasta_dir"]
RULE_DIR = config["working_dir"]["rule_dir"]

# Files
ALL_PEAKS = config["file"]["all_peaks"]

# Name of files
SPECIES = config["info"]["species"]
ASSEMBLY = config["info"]["assembly"]
PEAKCALLER = config["info"]["peakcaller"]
REMAP_VERSION = config["info"]["remap_version"]
UPDATE_VERSION = config["info"]["update_version"]
PATCH_VERSION = config["info"]["patch_version"]

PREFIX = "remap" + REMAP_VERSION
SUFFIX =  PEAKCALLER + "_" + ASSEMBLY + "_v" + UPDATE_VERSION + "_" + PATCH_VERSION
NB_CELL_INFO = int( config["info"]["nb_cell_info"])
NR_NAME = os.path.join(PREFIX + "_nr_" + SUFFIX + ".bed")

#================================================================#
#                         Includes                               #
#================================================================#


include: os.path.join(BASE_DIR, RULE_DIR, "tf_allPeaks.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "tf_nrPeaks_median_size_peak.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "tf_nrPeaks_merge.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "tf_nrPeaks.rules")

include: os.path.join(BASE_DIR, RULE_DIR, "cell_allPeaks.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "cell_nrPeaks_median_size_peak.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "cell_nrPeaks_merge.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "cell_nrPeaks.rules")


#================================================================#
#                     Defining dataset                           #
#================================================================#

set_tf = set()
set_cell = set()
dict_tf_in_cell = {}
# Get all TF
file_peaks = open( ALL_PEAKS, 'r')

for line in file_peaks:
    info_column = line.strip().split( "\t")[ 3]
    split_info_column = info_column.split( ".", 2)
    tf = split_info_column[ 1].split( "_")[ 0]
    set_tf.add( tf)

    cell = "_".join( split_info_column[ -1].split( "_")[ 0: NB_CELL_INFO])
    set_cell.add( cell)

    if cell not in dict_tf_in_cell:
        dict_tf_in_cell[ cell] = set()
    dict_tf_in_cell[ cell].add( tf)

file_peaks.close()

list_tf = []
list_cell = []

for key_cell in dict_tf_in_cell:
    for value_tf in dict_tf_in_cell[ key_cell]:
        list_tf.append( value_tf)
        list_cell.append( key_cell)

#================================================================#
#                         Workflow                               #
#================================================================#


rule all:
    input:

            expand( os.path.join( BED_DIR, "TF", "{tf}" , PREFIX + "_{tf}_nr_" + SUFFIX + ".bed"), tf = set_tf),
            expand( os.path.join( BED_DIR, "TF", "{tf}" , PREFIX + "_{tf}_all_" + SUFFIX + ".bed"), tf = set_tf),
            expand( os.path.join( FASTA_DIR, "{tf}" , PREFIX + "_{tf}_nr_" + SUFFIX + ".fasta"), tf = set_tf),
            expand( os.path.join( BED_DIR, "CELL","{cell}" , PREFIX + "_{cell}_nr_" + SUFFIX + ".bed"), cell = set_cell),
            expand( os.path.join( BED_DIR, "CELL","{cell}" , PREFIX + "_{cell}_all_" + SUFFIX + ".bed"), cell = set_cell),
            NR_NAME




rule tf_nrPeaks_fasta:
    input:
            os.path.join( BED_DIR, "TF", "{tf}" , PREFIX + "_{tf}_nr_" + SUFFIX + ".bed")
    output:
            os.path.join( FASTA_DIR, "{tf}" , PREFIX + "_{tf}_nr_" + SUFFIX + ".fasta")
    singularity:
            config[ "singularity"][ "bedtools"]
    conda:
            config[ "conda"][ "bedtools"]
    params:
            genome_fasta = config["genome"]["fasta"]
    shell: """bedtools getfasta -fi {params.genome_fasta} -bed {input} -fo {output}"""
        

rule merge_all_nr_bed:
    input:
            expand( os.path.join( BED_DIR, "TF", "{tf}" , PREFIX + "_{tf}_nr_" + SUFFIX + ".bed"), tf = set_tf)
    output:
            os.path.join(PREFIX + "_nr_" + SUFFIX + ".bed")	
    params:
            unsort_bed = temp(os.path.join(PREFIX + "_nr_" + SUFFIX + "_unsort.bed"))
            
    shell: """cat  {input} > {params.unsort_bed}
              sort -k1,1 -k2,2n {params.unsort_bed} > {output}"""
