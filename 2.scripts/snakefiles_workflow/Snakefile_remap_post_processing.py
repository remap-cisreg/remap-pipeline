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

# cluster info
# nb_thread_quality_NSC_RSC = cluster_config[ 'quality_NSC_RSC'][ 'thread']


#================================================================#
#                         Includes                               #
#================================================================#


# include: os.path.join(BASE_DIR, RULE_DIR, "merge_bam.rules")
# include: os.path.join(BASE_DIR, RULE_DIR, "quality_NSC_RSC.rules")
# include: os.path.join(BASE_DIR, RULE_DIR, "quality_FRiP.rules")
# include: os.path.join(BASE_DIR, RULE_DIR, "quality_all.rules")

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
    input:  #expand( os.path.join( QUALITY_DIR, "macs2_{experiment_name}.quality_all"), experiment_name = list_exp),
            # os.path.join( QUALITY_DIR, "results", "macs2.quality_all")
            # os.path.join( QUALITY_DIR, "results", "macs2_passed.quality_all"),
            # os.path.join( "remap2020_unsorted.bed"),
            os.path.join( "remap2020.bed"),
            # os.path.join( QUALITY_DIR,  "results", "macs2.quality_all.pdf")



rule merge_bam:
    input:
            bam = lambda wildcards : expand( os.path.join( BAM_DIR, "{chip_name}.bam"), chip_name = dict_experiment_chip_filename[wildcards.experiment_name])
    output:
            # temp( os.path.join( PREPROCESSING_DIR,  "merge_bam", "{experiment_name}.bam"))
            os.path.join( PREPROCESSING_DIR,  "merge_bam", "{experiment_name}.bam")
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
            # temp( os.path.join( PREPROCESSING_DIR, "NSC_RSC", "{experiment_name}.NSC_RSC"))
            os.path.join( PREPROCESSING_DIR, "NSC_RSC", "{experiment_name}.NSC_RSC")
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
    run_spp.R -c={input.bam} -savp -out={output} -p=2 -odir={params.outdir} -rf {params.run_spp}
    """



rule quality_FRiP_nbPeaks:
    input:
            # bam = lambda wildcards : expand( os.path.join( BAM_DIR, "{chip_name}.bam"), chip_name = dict_experiment_chip_filename[wildcards.experiment_name]["chip"]),
            bam = os.path.join( PREPROCESSING_DIR, "merge_bam", "{experiment_name}.bam"),
            peaks = os.path.join( PEAKCALLING_DIR, "{experiment_name}", "macs2", "{experiment_name}_peaks.narrowPeak"),
            nsc_rsc = os.path.join( PREPROCESSING_DIR, "NSC_RSC", "{experiment_name}.NSC_RSC")
    output:
            # temp( os.path.join( PREPROCESSING_DIR, "FRiP", "macs2_{experiment_name}.FRiP"))
            os.path.join( PREPROCESSING_DIR, "FRiP", "macs2_{experiment_name}.FRiP")
    singularity:
            config[ "singularity"][ "bedtools"]
    conda:
            config[ "conda"][ "bedtools"]
    resources:
            res=2
    log:
            os.path.join( PREPROCESSING_DIR, "FRiP", "log", "{experiment_name}_FRiP.log")
    params:
            bedtools_intersect = config[ "quality_FRiP"][ "bedtools-intersect"]

    shell:"""
        total_reads=$(cut -f2 {input.nsc_rsc})

        total_peaks_in_read=$(bedtools intersect {params.bedtools_intersect} -bed -u -a {input.bam} -b {input.peaks} | wc -l)

        FRiP=$(awk -v nr="$total_reads" -v ni="$total_peaks_in_read" "BEGIN{{print ni/nr * 100}}")

        total_peaks=$(wc -l < {input.peaks} )

        echo $FRiP\t$total_peaks > {output}
    """


rule quality_all_experiment:
    input:
            nsc_rsc = os.path.join( PREPROCESSING_DIR, "NSC_RSC", "{experiment_name}.NSC_RSC"),
            frip = os.path.join( PREPROCESSING_DIR, "FRiP", "macs2_{experiment_name}.FRiP")
    output:
            os.path.join( QUALITY_DIR, "macs2_{experiment_name}.quality_all")
    resources:
            res=1
    log:
            os.path.join( QUALITY_DIR, "log", "{experiment_name}quality_all_experiment.log")
    params:
            outdir = os.path.join( PREPROCESSING_DIR, "NSC_RSC"),
            exp_name = "{experiment_name}"
    shell:"""
        FRiP=$(cat {input.frip})
        NSC=$(cut -f9 {input.nsc_rsc})
        RSC=$(cut -f10 {input.nsc_rsc})

        echo "experiment_name\tNSC\tRSC\tFRiP\tnb_peaks\n{params.exp_name}\t$NSC\t$RSC\t$FRiP" > {output}
    """


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


rule filtering_quality_all:
    input:
            os.path.join( QUALITY_DIR,  "results", "macs2.quality_all")
    output:
            passed = os.path.join( QUALITY_DIR,  "results", "macs2_passed.quality_all"),
            failed = os.path.join( QUALITY_DIR,  "results", "macs2_failed.quality_all")
    resources:
            res=1
    log:
            os.path.join( QUALITY_DIR, "log", "filtering_quality_all.log")
    shell:"""

         awk '{{
                if( $5<100 || $9<2) {{ print > "{output.failed}"}}

                else {{ print > "{output.passed}"}}

              }}'  {input}

    """


# rule graph_quality_all:
#     input:
#             os.path.join( QUALITY_DIR,  "results", "macs2.quality_all")
#     output:
#             os.path.join( QUALITY_DIR,  "results", "macs2.quality_all.pdf")
#     conda:
#             config[ "conda"][ "R_quality"]
#     resources:
#             res=1
#     log:
#             os.path.join( QUALITY_DIR, "log", "graph_quality_all.log")
#
#     shell:"Rscript 2.scripts/utils/r/graph_quality_all.R --file {input} --outfile {output}"



rule creating_remap_catalogue_bed:
    input:
            os.path.join( QUALITY_DIR,  "results", "macs2_passed.quality_all")
    output:
            os.path.join( "remap2020_unsorted.bed")
    resources:
            res=1
    log:
            os.path.join( QUALITY_DIR, "log", "creating_remap_catalogue_bed.log")
    params:
               path_to_peaks = PEAKCALLING_DIR,
               peakcaller = PEAKCALLER,
               path_file_color_tf = config[ "creating_remap_catalogue_bed"][ "color_file"]
    run:

        # create dicrionary with color by tf
        dict_tf_color = {}
        file_color_tf = open(  params.path_file_color_tf, 'r')

        for line in file_color_tf:
            split_line = line.strip().split( "\t")
            if split_line[ 0] not in dict_tf_color:
                dict_tf_color[ split_line[ 0]] = str( split_line[ 1]) + "," + str( split_line[ 2]) + "," + str( split_line[ 3])
        file_color_tf.close()


        # Puting in a list all experiement passing filtering
        list_experiment_pass = []
        file_experiment_pass = open(  input[ 0], 'r')

        for line in file_experiment_pass:
            if not line.startswith( "experiment_name"):
                list_experiment_pass.append( line.split( "\t")[ 0].strip())
        file_experiment_pass.close()

        # writing results
        file_outfile = open( output[ 0], 'w')

        for experiment in list_experiment_pass:

            # getting correct color for tf
            tf = experiment.split( ".")[ 1]
            color_tf = dict_tf_color[ tf]

            # formating narrowPeak ReMap style
            path_file_experiment_peak =  os.path.join( params.path_to_peaks, experiment, params.peakcaller, experiment + "_peaks.narrowPeak")

            file_experiment_peak = open( path_file_experiment_peak, 'r')

            for line in file_experiment_peak:
                split_line = line.strip().split( "\t")
                summit_position = int( split_line[ -1]) + int( split_line[ 1])
                correct_line = split_line[ 0:4] + split_line[ 8:9] + split_line[ 5:6] + [  str( summit_position)] + [ str( summit_position + 1)] + [ color_tf]
                file_outfile.write( "\t".join(correct_line) + "\n")

        file_outfile.close()



rule sort_remap_bed:
    input:
            os.path.join( "remap2020_unsorted.bed")
    output:
            # temp( os.path.join( PREPROCESSING_DIR,  "merge_bam", "{experiment_name}.bam"))
            os.path.join( "remap2020.bed")
    resources:
            res=1
    log:
            os.path.join( QUALITY_DIR, "log", "sort_remap_bed.log")
    params:
            other = ""

    shell:"""sort -k1,1 -k2,2n {input} > {output}"""
