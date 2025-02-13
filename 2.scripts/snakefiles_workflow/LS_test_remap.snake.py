###LS_test_remap.snake



include: os.path.join(BASE_DIR, RULE_DIR, "aria2.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "trim_galore.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "bowtie2.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "sam_to_bam.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "remove_mismatch.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "samtools_sort.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "samtools_remove_pcr_duplicate.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "macs2.rules")
include: os.path.join(BASE_DIR, RULE_DIR, "delete_trim.rules")



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

# Extension peakcaller
EXTENSION_PEAK = config[ "extention"][ "peak"]

#================================================================#
#                     Defining dataset                           #
#================================================================#

SAMPLE, = glob_wildcards(TAB_DIR  + "{sample}_summary.tab")



##########################
#Script post-basecalling#
#########################

#se lance depuis la racine du projet
# Doit contenir: 
# - dossier basecall avec un fichier bam nommé: {projet}_basecall.bam	<= Prévoir en cas de basecall au format fastQ

# Exemple de commande:
#DATA_DIR="/home/gbim/Documents/DATA/2023-42_set4" FAST_DIR="sample_name/aleatoire/" REF="hg38" snakemake.mt -s /home/gbim/Documents/GIT/GBIM/snake/longread.snake --dryrun

################
# Path dossiers# 		<= Proposition passer en alias ou variable sys 		Négatif: Rend le script moins reproductible en lui même 	Voir avec Camille
################
GIT= "/home/gbim/Documents/GIT/"
SYNO= "/home/gbim/syno"
TOOLS = SYNO + "/Bioinfo_scripts/TOOLS"
GDB = SYNO + "/Bioinfo_scripts/genomic_DB"
SCRIPTS = GIT + "GBIM"
R_SCRIPTS = SCRIPTS + "/R/"

import pathlib
DATA_DIR = os.environ.get("DATA_DIR")
REF = os.environ.get("REF")
# FAST_DIR = os.environ.get("FAST_DIR")

# Extraction nom du projet (Extrait le dernier dossier de DATA_DIR)
DIR_NAME = os.path.basename(os.path.normpath(DATA_DIR))

print ("DATA_DIR :", DATA_DIR)
print ("REF :", REF)
# print ("FAST_DIR :", FAST_DIR)
print("DIR_NAME:", DIR_NAME)

# print("Do you which to analyse methylation YES or NO")
# METH = ""
# while METH not in {"YES","NO"}:
# 	METH = input("	Please enter YES or NO: 	")



#FAST_PATH = DATA_DIR + FAST_DIR 



LOG_DIR = DATA_DIR + "/log/"
QC_DIR = DATA_DIR + "/QC/"
OUTPUTS_DIR = DATA_DIR + "/outputs/"
BASECALL_DIR = DATA_DIR + "/basecall/"
MAP_DIR = DATA_DIR + "/map/"
COUNT_MRNA_DIR = DATA_DIR + "/count_mRNA/"
VAR_DIR = DATA_DIR + "/var/"

# SAMPLE = glob_wildcards(FAST_DIR + "sequencing_summary*")

SAMPLE, = glob_wildcards(BASECALL_DIR  + "{sample}_basecall.bam")
print("SAMPLE : ",SAMPLE)



if REF == "hg38":
	GP = GDB + "/GRCh38/"
	GTF_FILE = GP + "gencode.v34.annotation.gtf"
	# GTF_FILE = GP + "gencode.v34.long_noncoding_RNAs.gtf"
	FASTA = GP + "GRCh38.p13.genome.fa"
	org = "org.Hs.eg.db"

if REF == "mm39":
	GP = GDB + "/GRCm39/"
	GTF_FILE = GP + "Mus_musculus.GRCm39.108.chr.gtf"
	FASTA = GP + "mm39.fa"
	org = "org.Mm.eg.db"

if REF == "rheMac10":
	GP = GDB + "/Mmul_10/"
	GTF_FILE = GP + "Macaca_mulatta.Mmul_10.105.chr.gtf"
	FASTA = GP + "/index_BWA/Macaca_mulatta.Mmul_10.fa"

if REF == "GRCh37":
	GP = GDB + "/GRCh37/"
	GTF_FILE = GP + "gencode.v26lift37.annotation.gtf"
	FASTA = GP + "GRCh37.p13.genome.fa"

if REF == "hg19":
	GP = GDB + "/GRCh37/"
	GTF_FILE = GP + "gencode.v26lift37.annotation.gtf"
	# FASTA = GP + "hs37d5.fa"
	FASTA = "/home/gbim/Documents/DATA/ref/hs37/hs37d5.fa"

	INDEX_DIR = GP + "/index_STAR_hs37d5"

if REF == "hg38" or REF == "mm39" or REF == "rheMac10" or REF == "GRCh37":
	INDEX_DIR = GP + "/index_STAR"
	print ("INDEX_DIR : ", INDEX_DIR)

if REF == "T2T" :
	GP = GDB + "/T2T-CHM13/"
	# GTF_FILE = GP + "GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz"
	# BED_FILE = GP + "T2T_ncbiRefSeqCurated_UCSC.bed"
	#FASTA = GP + "GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"
	FASTA = GP + "hg002v1.0.1.fasta"
	# FASTA = "/home/gbim/Documents/DATA/ref/T2T/hg002/hg002v1.0.1.fasta"
if REF == "custom":
	FASTA = os.environ.get("FASTA")
	BED_FILE = os.environ.get("BED")
	print(FASTA)
	print(BED_FILE)
	# if BED_FILE is None:
	# 	print("working")
	
	# else: 
	# 	print("failed")


# pod_loc=$(find -name "pod5" | grep DIR_NAME)

fasta_name= os.path.basename(os.path.normpath(FASTA))
# Mode de fonctionnement de Pepper-Margin-DV
mode = "ont_r10_q20"


rule final:
	input:
		# POD5
		# FAST_PATH + "/pod5/" ,
		# Basecalling
		# BASECALL_DIR + DIR_NAME + "_basecall.bam",
		# alignement
		# MAP_DIR + DIR_NAME + ".aln.bam",
		expand(MAP_DIR  + "{sample}.aln.sort.bam", sample=SAMPLE),
		# ## VARIANTS
		# #svim
		# VAR_DIR + "SV/" + "svim/{sample}/variants.vcf"
		# # NanoSV		
		# # VAR_DIR + DIR_NAME + "_nanosv_variants.vcf",
		# #cuteSV
		# VAR_DIR + "SV/" + "{sample}._cuteSV_variants.vcf",		
		# #sniffles2
		# VAR_DIR + "SV/" + "{sample}.sniffles2_variants.vcf",
		##SURVIVOR OUTPUT
		expand(VAR_DIR +"SV/"+ "{sample}._SVmerged.vcf", sample=SAMPLE),
		expand(VAR_DIR + "SV/" + "{sample}._cuteSV_variants.vcf", sample=SAMPLE),
		expand(VAR_DIR + "SV/" + "{sample}.sniffles2_variants.vcf", sample=SAMPLE),
		expand(VAR_DIR + "SV/" + "svim/{sample}/variants.vcf", sample=SAMPLE),
		##SNV
		# expand(VAR_DIR + "{sample}/pepper/"  + REF + mode +  "{sample}_PEPPER_Margin_DeepVariant.vcf.gz", sample=SAMPLE),
		
		##QC
		#pycoqc
		# expand(QC_DIR + "pycoQC/" + "{sample}_PycoQC-report.html", sample=SAMPLE),
		#Nanoplot
		expand(QC_DIR + "nanoplot/" + "{sample}_NanoPlot-report.html", sample=SAMPLE),
		##Methy
		expand(DATA_DIR + "/METHYL/" +"{sample}._modkit.bedmethyl", sample=SAMPLE)




# rule Fast5_2_Pod5:
# 	input:
# 		# FAST5_DIR + "*.fast5"
# 		FAST_PATH + "/fast5/"
# 	output:
# 		directory(FAST_PATH + "/pod5/")
# 	# conda:
# 	# 	"/home/gbim/Bureau/Projets_en_cours/conda/R3.yaml"
# 	threads:
# 		30
# 	shell:
# 		"""
# 		mkdir -p {output}
# 		pod5 convert fast5 {input} --output {output}/
# 		"""



## Basecall & Align

# if METH == "Y":
# 	rule dorado_basecalling:
# 		input:
# 			FAST_PATH + "/pod5/"
# 		output:
# 			BASECALL_DIR + "{DIR_NAME}_basecall.bam"
# 		shell:
# 			"""
# 			dorado basecaller sup,5mCG_5hmCG {input} > {output}
# 			dorado summary {output} > basecall_summary.tsv
# 			"""
# else:
# 	rule dorado_basecalling:
# 		input:
# 			FAST_PATH + "/pod5/"
# 		output:
# 			BASECALL_DIR + "{DIR_NAME}_basecall.bam"
# 		shell:
# 			"""

# 			dorado basecaller sup {input} > {output}
# 			dorado summary {output} > basecall_summary.tsv
# 			"""


rule dorado_align:
	input:
		BASECALL_DIR + "{sample}_basecall.bam"
	output:
		# MAP_DIR + DIR_NAME + ".aln.bam",
		BAM_FILE = MAP_DIR  + "{sample}.aln.sort.bam",
		BAM_INDEX =  MAP_DIR  + "{sample}.aln.sort.bambai"
	log:
		LOG_DIR + ".{sample}_dorado_align.log"
	shell:
		"""
		#dorado aligner  -k 15 -w 10 --bandwidth 500,20000 {FASTA} {input} > {output}
		dorado aligner -t 50 --bandwidth 500,20000 {FASTA} {input} | samtools sort -@ 50 -o {output.BAM_FILE}
		samtools index -@ 50 {output.BAM_FILE} {output.BAM_INDEX}
		"""

# rule model_extract:
# 	input:
# 		MAP_DIR + DIR_NAME + ".aln.bam"
# 	output:
# 		DATA_DIR + "model_info.txt"
# 		shell:
# 		"""
# 		samtools view {input} |head |awk -F$"\t" 'NR>3 && NR <5 {print $7}'| cut -d' ' -f1,2 > model_info.txt
# 		""" 


# rule Bamindex :
# 	input:
# 		MAP_DIR + DIR_NAME + ".aln.bam"
# 	output:
# 		BAM_FILE = MAP_DIR +  DIR_NAME + ".aln.sort.bam",
# 		BAM_INDEX = MAP_DIR +  DIR_NAME + ".aln.sort.bam.bai"
# 	log:
# 		LOG_DIR + DIR_NAME + ".Bamindex.log"
# 	shell:
# 		"""
# 		samtools sort {input} > {output.BAM_FILE}
# 		samtools index {output.BAM_FILE} {output.BAM_INDEX} 
# 		"""



##QC

rule fastQC:
	input:
		BASECALL_DIR + "{sample}.bam"
	output:
		QC_DIR + "/FASTQC/" + "{sample}_fastqc.zip"
	# conda:
	# 	"/home/gbim/Bureau/Projets_en_cours/conda/R3.yaml"
	threads:
		4
	shell:
		"""
		fastqc -o {QC_DIR}/FASTQC/ -t {threads} {input}
		"""

rule nanoplot:
	input:
		BAM_FILE = MAP_DIR  + "{sample}.aln.sort.bam",
		# BAM_INDEX = MAP_DIR + DIR_NAME + ".aln.bam.bai"
	output:
		NP_report = QC_DIR + "nanoplot/" + "{sample}_NanoPlot-report.html",
	threads:
		10
	params:
		prefix = "{sample}"
	conda:
		"nanoplot"
	shell:
		"""
		NanoPlot -t {threads} --verbose --prefix {params.prefix}_ -o {QC_DIR}nanoplot --bam {input.BAM_FILE} --only-report
		
		""" 



# rule pycoQC:
# 	input:
# 		BAM_FILE = MAP_DIR  + "{sample}.aln.sort.bam",
# 		# BAM_INDEX = MAP_DIR + DIR_NAME + ".aln.bam.bai"
		

# 	output:
# 		PQC_report = QC_DIR + "pycoQC/" + "{sample}_PycoQC-report.html",
	
# 	log:
# 		LOG_DIR + "{sample}.pycoqc.log"
# 	threads:
# 		10
# 	conda:
# 	 	"pycoQC"
# 	shell:
# 		"""
# 		REPORT=$(find -name sequencing_summary*)
# 		# echo $REPORT
# 		pycoQC --summary_file $REPORT --bam_file {input.BAM_FILE} -o {output}
		
		# """ 



## Variant calling

rule svim_variant_calling:
	input:
		BAM_FILE = MAP_DIR  + "{sample}.aln.sort.bam",
		# BAM_INDEX = MAP_DIR + DIR_NAME + ".aln.bam.bai"

	output:
		# VAR_DIR + DIR_NAME + "_svim_variants.vcf"
		directory(VAR_DIR + "SV/svim/{sample}/")
	params:
		NAME = "{sample}"
	conda:
		"svim_env"
	shell:
		"""
		svim alignment --sample {params.NAME} {output} {input.BAM_FILE} {FASTA}
		""" 

###########
#ATTENTION#
###########
#Ceci est un pansement horrible, il faudrait une meilleure solution
# rule nanosv_variant_calling:
# 	input:
# 		BAM_FILE = MAP_DIR + DIR_NAME + ".aln.sort.bam",
# 		# BAM_INDEX = MAP_DIR + DIR_NAME + ".aln.bam.bai"
# 	output:
# 		VAR_FILE = VAR_DIR + DIR_NAME + "_nanosv_variants.vcf",
# 		NANOSV_BED = temp(VAR_DIR + "random.bed")
# 	conda:
# 		"nanosv"
# 	threads:
# 		20
# 	shell:
# 		"""
# 		bedtools random -l 1 -seed 1251 -g {FASTA}.fai > {VAR_DIR}random.bed
# 		sam=$(which samtools)
# 		echo $sam
# 		NanoSV -t {threads} -o {output.VAR_FILE} -b {output.NANOSV_BED} {input.BAM_FILE} -s $sam
# 		"""

rule sniffles2_variant_calling:
	input:
		BAM_FILE = MAP_DIR  + "{sample}.aln.sort.bam",
		# BAM_INDEX = MAP_DIR + DIR_NAME + ".aln.bam.bai"
	output:
		# VAR_DIR + DIR_NAME + "_svim_variants.vcf"
		VAR_DIR + "SV/" + "{sample}.sniffles2_variants.vcf"
	conda:
		"sniffles2"
	threads:
		8
	shell:
		"""
		sniffles -i {input.BAM_FILE} --threads {threads} -v {output}
		""" 

rule cuteSV_variant_calling:
	input:
		BAM_FILE = MAP_DIR  + "{sample}.aln.sort.bam",
		# BAM_INDEX = MAP_DIR + DIR_NAME + ".aln.bam.bai"
	output:
		# VAR_DIR + DIR_NAME + "_svim_variants.vcf"
		output_dir = directory(VAR_DIR + "SV/cutesv/{sample}/"),
		vcf = VAR_DIR +"SV/" + "{sample}._cuteSV_variants.vcf",
	# params:
	# 	output_dir = VAR_DIR + "SV/{sample}/"
	conda:
		"cutesv"
	shell:
		"""
		mkdir {output.output_dir}
		cuteSV {input.BAM_FILE} {FASTA} {output.vcf} {output.output_dir}

		""" 



# rule SVIM_ASM_variant_calling:
# 	input:
# 		BAM_FILE = MAP_DIR + DIR_NAME + ".aln.sort.bam",
# 		# BAM_INDEX = MAP_DIR + DIR_NAME + ".aln.bam.bai"
# 	output:
# 		# VAR_DIR + DIR_NAME + "_svim_variants.vcf"
# 		VAR_DIR + "sniffles2_variants.vcf"
# 	conda:
# 		"sniffles2"
# 	shell:
# 		"""
# 		sniffles -i {input.BAM_FILE} -v {output}
# 		""" 

# rule NanoVAR_variant_calling:
# 	input:
# 		BAM_FILE = MAP_DIR + DIR_NAME + ".aln.sort.bam",
# 		# BAM_INDEX = MAP_DIR + DIR_NAME + ".aln.bam.bai"
# 	output:
# 		# VAR_DIR + DIR_NAME + "_svim_variants.vcf"
# 		VAR_DIR + "sniffles2_variants.vcf"
# 	conda:
# 		"sniffles2"
# 	shell:
# 		"""
# 		sniffles -i {input.BAM_FILE} -v {output}
# 		""" 


rule VCF_list:
	input:
		expand(VAR_DIR + "SV/" + "{sample}._cuteSV_variants.vcf", sample=SAMPLE),
		expand(VAR_DIR + "SV/" + "{sample}.sniffles2_variants.vcf", sample=SAMPLE),
		expand(VAR_DIR + "SV/" + "svim/{sample}/", sample=SAMPLE),
	output:
		VCF_LIST = VAR_DIR + "SV/" + "{sample}_list_vcf.txt",
	params:
		name="{sample}"
	shell:
		"""
		find -maxdepth 4 -name '{params.name}.*.vcf' > {output}

		"""

rule SURVIVOR:
	input:
		VAR_DIR + "SV/" + "{sample}._cuteSV_variants.vcf",
		VAR_DIR + "SV/" + "{sample}.sniffles2_variants.vcf",
		VAR_DIR + "SV/svim/{sample}/",
		VCF_LIST = VAR_DIR + "SV/" + "{sample}_list_vcf.txt",
	output:
		# VCF_LIST = VAR_DIR + "{sample}_list_vcf.txt",
		MERGED_SV = VAR_DIR +"SV/"+ "{sample}._SVmerged.vcf"
	shell:
		"""

		SURVIVOR merge {input.VCF_LIST} 1000 2 1 1 0 30 {output.MERGED_SV}
		"""

rule methylation_bed_modkit:
	input:
		BAM_FILE = MAP_DIR + "{sample}.aln.sort.bam",
	output:
		DATA_DIR + "/METHYL/" +"{sample}._modkit.bedmethyl"
	threads:
		20
	shell:
		"""
		modkit pileup --ref {FASTA} --threads {threads} --interval-size 100000 {input} {output}
		"""


rule pepper:
	input:
		MAP_DIR + "{sample}.aln.sort.bam",
	output:
		vcf = VAR_DIR + "{sample}/pepper/"  + REF + mode +  "{sample}_PEPPER_Margin_DeepVariant.vcf.gz"
	params:
		inputbam = "/input/map/{sample}.aln.sort.bam",
		refFasta = "/ref/"+ fasta_name,
		threads = "64",
		# outfolder = "/input/var/pepper",
		outfolder = "/input/var/{sample}/",
		outfile = REF + mode +  "{sample}_PEPPER_Margin_DeepVariant"
	log:
		LOG_DIR + "{sample}.pepper.log"
	shell:
		"""
	docker run -v {DATA_DIR}:"/input" -v {GP}:"/ref" kishwars/pepper_deepvariant:r0.8  \
	run_pepper_margin_deepvariant call_variant \
		-b {params.inputbam} \
		-f {params.refFasta} \
		-o {params.outfolder} \
		-p {params.outfile} \
		-t {params.threads} \
		--{mode} 
		"""


rule CNV_spectre:
	input:
		BAM_FILE = MAP_DIR + "{sample}.aln.sort.bam",
		SV = VAR_DIR + "SV/" + "{sample}.sniffles2_variants.vcf",
		vcf = VAR_DIR + "pepper/" + REF + mode +  "{sample}_PEPPER_Margin_DeepVariant.vcf.gz"
	output:
		mosdepth_regions = "mosdepth/{sample}.regions.bed.gz",
		snfj = "CNV/spectre/{sample}.sniffles2_variants.snfj"
	params:
		mosdepth_dir = "mosdepth/{sample}",
		CNV = "CNV/spectre",
		recommended_spectre = "-x -b 1000 -Q 20",
		NAME = "{sample}",
		threads = 8
		
	shell:
		"""
		snf2json {input.SV} {output.snfj}
		mosdepth -t {params.threads} {params.recommended_spectre} {params.mosdepth_dir} {input.BAM_FILE}
		spectre CNVCaller -c {output.mosdepth_regions} --sample-id {params.NAME} --output-dir {params.CNV} --snv {input.vcf}  --reference {fasta}

		"""