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





#include: os.path.join(BASE_DIR, RULE_DIR, "delete_trim.rules")

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
    #print( split_info_column)
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

# print( list_tf)
# print( list_cell)
# print( dict_tf_in_cell)
#================================================================#
#                         Workflow                               #
#================================================================#


rule all:
    input:
            # expand( os.path.join( BED_DIR, "{tf}" , PREFIX + "_{tf}_all_" + SUFFIX + ".bed"), tf = set_tf)
            # expand( os.path.join( BED_DIR, "{tf}" , PREFIX + "_{tf}_nr_" + SUFFIX + "_part2.bed"), tf = set_tf)
            expand( os.path.join( BED_DIR, "TF", "{tf}" , PREFIX + "_{tf}_nr_" + SUFFIX + ".bed"), tf = set_tf),
            expand( os.path.join( BED_DIR, "TF", "{tf}" , PREFIX + "_{tf}_all_" + SUFFIX + ".bed"), tf = set_tf),
            expand( os.path.join( FASTA_DIR, "{tf}" , PREFIX + "_{tf}_nr_" + SUFFIX + ".fasta"), tf = set_tf),
            # expand( os.path.join( BED_DIR, "CELL","{cell}" , PREFIX + "_{cell}_all_" + SUFFIX + ".bed"), cell = set_cell)
            expand( os.path.join( BED_DIR, "CELL","{cell}" , PREFIX + "_{cell}_nr_" + SUFFIX + ".bed"), cell = set_cell),
            expand( os.path.join( BED_DIR, "CELL","{cell}" , PREFIX + "_{cell}_all_" + SUFFIX + ".bed"), cell = set_cell),
            # expand( os.path.join( BED_DIR, "CELL","{cell}" , PREFIX + "_{cell}_{tf}_" + SUFFIX + "_temp.bed"), zip, cell = list_cell, tf = list_tf)
            # expand( os.path.join( BED_DIR, "CELL","{cell}" , PREFIX + "_{cell}_{tf}_" + SUFFIX + "_part2.bed"), zip, cell = list_cell, tf = list_tf)
            # expand( os.path.join( BED_DIR, "CELL","{cell}" , PREFIX + "_{cell}_{tf}_" + SUFFIX + "_part3.bed"), zip, cell = list_cell, tf = list_tf),
            # expand( os.path.join( BED_DIR, "CELL", "{cell}", PREFIX + "_{cell}_nr_tf_" + SUFFIX + ".bed"), cell = set_cell)



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


# rule cell_nrPeaks_separate_tf:
#     input:
#             os.path.join( BED_DIR, "CELL", "{cell}", PREFIX + "_{cell}_nr_" + SUFFIX + "_part1.bed")
#     output:
#             temp( os.path.join( BED_DIR, "CELL", "{cell}", PREFIX + "_{cell}_{tf}_" + SUFFIX + "_temp.bed"))
#     shell: "grep -E '\.{wildcards.tf}(_[a-zA-Z0-9]+)*\.' {input} > {output}"
#
#
#
# rule cell_tf_nrPeaks_merge:
#     input:
#             os.path.join( BED_DIR, "CELL", "{cell}", PREFIX + "_{cell}_{tf}_" + SUFFIX + "_temp.bed")
#     output:
#             temp( os.path.join( BED_DIR, "CELL", "{cell}", PREFIX + "_{cell}_{tf}_" + SUFFIX + "_part2.bed"))
#     singularity:
#             config[ "singularity"][ "bedtools"]
#     conda:
#             config[ "conda"][ "bedtools"]
#     params:
#             per_overlap = config["cell_nrPeaks"]["perc_overlap"]
#     shell: """
#
#     # Getting number of nucleotide from median size to reprecent asked percentage overlap
#     MEDIAN_SIZE=$(awk -F"\t" '{{ print $3-$2}}' {input} | sort -n | awk ' {{ a[i++]=$1; }} END {{ print a[int(i/2)]; }}')
#     NT_OVERLAP=$(echo $MEDIAN_SIZE | awk -v percOverlap='{params.per_overlap}' '{{print int($1/100*percOverlap)}}')
#
#     bedtools merge -d -$NT_OVERLAP -c 4,7,9 -o collapse,mean,distinct -i {input} | awk -F"\t" '{{print $1"\t"$2"\t"$3"\t"$4"\t"int($5)"\t"$6}}' > {output}
#     """
#
#
# rule cell_tf_nrPeaks:
#     input:
#             os.path.join( BED_DIR, "CELL", "{cell}", PREFIX + "_{cell}_{tf}_" + SUFFIX + "_part2.bed")
#     output:
#             temp( os.path.join( BED_DIR, "CELL", "{cell}", PREFIX + "_{cell}_{tf}_" + SUFFIX + "_part3.bed"))
#     run:
#         import statistics
#
#         file_merge = open( input[ 0], 'r')
#         file_output = open( output[ 0], 'w')
#
#
#         for line in file_merge:
#             # get all relevant info
#             split_line = line.strip(). split( "\t")
#             chromosome = split_line[ 0]
#             summit = int( split_line[ 4])
#             color = split_line[ 5]
#
#             # Spliting info column by peaks in nr
#             raw_peaks_in_nr =  split_line[ 3].split( ",")
#
#             set_tf = set()
#
#             list_peaks_start = []
#             list_peaks_end = []
#
#             for peaks in  raw_peaks_in_nr:
#
#                 # Spliting info and position
#                 split_peaks_in_nr = peaks.split( ":")
#
#                 list_peaks_start.append ( int( split_peaks_in_nr[ 1].split( "-")[ 0]))
#                 list_peaks_end.append ( int( split_peaks_in_nr[ 1].split( "-")[ 1]))
#
#                 # get raw cell line
#                 raw_tf = split_peaks_in_nr[ 0].split( ".", 3)[ 1].split( "_", 1)[ 0]
#                 set_tf.add( raw_tf)
#
#
#             mean_start =  int( round( statistics.mean( list_peaks_start), 0))
#
#             mean_end =  int( round( statistics.mean( list_peaks_end), 0))
#
#             file_output.write( "\t".join([ chromosome, str( mean_start),  str( mean_end), wildcards.cell + ":" + ",".join( list( set_tf)), str( len( raw_peaks_in_nr)), ".", str( summit), str( summit + 1), color]) +"\n")
#
#
#         file_merge.close()
#         file_output.close()
#
# rule cell_tf_concatenate:
#     input:
#             lambda wildcards : expand( os.path.join( BED_DIR, "CELL", "{cell}", PREFIX + "_{cell}_{tf}_" + SUFFIX + "_part3.bed"), tf = dict_tf_in_cell[ wildcards.cell], cell = wildcards.cell)
#     output:
#             os.path.join( BED_DIR, "CELL", "{cell}", PREFIX + "_{cell}_nr_tf_" + SUFFIX + "_unsort.bed")
#     shell: "cat {input} > {output}"
#
#
# rule cell_tf_sort:
#     input:
#             os.path.join( BED_DIR, "CELL", "{cell}", PREFIX + "_{cell}_nr_tf_" + SUFFIX + "_unsort.bed")
#     output:
#             os.path.join( BED_DIR, "CELL", "{cell}", PREFIX + "_{cell}_nr_tf_" + SUFFIX + ".bed")
#     shell: "sort -k1,1 -k2,2n {input} > {output}"
