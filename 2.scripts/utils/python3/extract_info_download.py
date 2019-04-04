#!/usr/bin/env 


'''
From TF annoation file to necessary files to run the workflow
input: tab file
output: one file by experiment (line) containing info read by workflow
'''
__author__ = "Jeanne Chèneby"
__version__ = "0.1"
__status__ = "Production"

import os
import argparse
import urllib.request

import requests

# Check if correct format for a ID
def check_ID_format( ID_to_test):

	if ID_to_test.startswith( "GSE"):
		return True
	elif ID_to_test.startswith( "ERP"):
		return True
	else:
		return False

# Check if correct SRP format
def check_SRP_format( elm_to_test):

	if elm_to_test.startswith( "SRP"):
		return True
	else:
		return False

# Check if correct format for a ID
def check_sampleID_format( elm_to_test):

	if elm_to_test.startswith( "GSM"):
		return True
	elif elm_to_test.startswith( "ERR"):
		return True
	else:
		return False

# dealing with GSM problems
def numbering_replicat( list_remap_sampleID, list_ena_line):
	# priming biological replicat
	nb_rep_bio = 1
	list_sample_name = []

	# parsing all sample of experiment
	for remap_GSM in list_remap_sampleID:

		# priming technical replicat
		nb_rep_tech = 1

		for ena_line in list_ena_line:

			if remap_GSM in ena_line:

				list_sample_name.append( [ str( nb_rep_bio), str( nb_rep_tech)] + ena_line)

				nb_rep_tech += 1

		
		nb_rep_bio += 1
	return  list_sample_name

if __name__ == "__main__":

	#================================================================#
	#                        Working path_outdir		             #
	#================================================================#

	parser = argparse.ArgumentParser(
        description="Create a summary file per experiment line containing all info necessary to run the remap workflow",
        epilog="Jeanne Cheneby")
	parser.add_argument("-wd", "--working_directory",
                        help="Working derictory where are folders 1.metadata, 2.script, ...",
                        action="store",
                        required=True)
	parser.add_argument("-f", "--annotation_file", help="path to file containing all annotation of ChIP-seq",
                        action="store",
                        required=True)
	args = parser.parse_args()






	# setting project working directory
	# path_working_dir = "/gpfs/tagc/home/cheneby/latest_remap"
	path_working_dir = args.working_directory
	os.chdir( path_working_dir)

	# setting TF file
	# path_file_remap_metadata = "1.metadata/TF_catalogue_run_20_01_17_FINAL.tab"
	path_file_remap_metadata = args.annotation_file

	path_outdir = "4.preprocessing"


	#================================================================#
	#                        Necessary variable			             #
	#================================================================#

	general_error =''





	# setting column of interest
	col_ID_serie = 0
	col_TF = 1
	col_cell = 2
	col_SRP = 6
	col_GSM_chip = 3
	col_GSM_control = 4


	# setting rest info
	url_experiment_srt = "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession="


		# info to get
	url_experiment_temp_end_1 = "&result=read_run&fields="
	url_experiment_temp_end_2 = "&download=txt"


	list_info_ena = [
						"study_alias", 					# GSE
						"secondary_study_accession", 	# SRP	
						"experiment_accession", 		# SRX	
						"run_accession", 				# SRR										
						"sample_alias", 				# GSM
						"library_layout", 				# single/paired
						"fastq_ftp", 					# Download link
						"fastq_md5", 					# server md5 of fastq 
						"sample_title"					# Name of sample (non standart)

					]


	info_ena =  ",".join( list_info_ena)

		# final url end
	url_experiment_end = url_experiment_temp_end_1 + info_ena + url_experiment_temp_end_2



	#================================================================#
	#                        Start of script			             #
	#================================================================#



	# get relevant info
	file_remap_metadata = open( path_file_remap_metadata, 'r')

	for line in file_remap_metadata:
		general_error = []

		split_line = line.strip().split( "\t")

		try:	
		# getting serie's ID
			ID_serie = split_line[ col_ID_serie].strip()

			# Check if correct format for a ID
			if check_ID_format( ID_serie):
				pass
			else:
				ID_serie = "wrong_ID_serie_format"
				general_error.append( "wrong_ID_serie_format")
				
		except IndexError:
			general_error.append( "missing_ID_column")


		try:
		# getting TF's name
			TF = split_line[ col_TF].strip()

			# Check if correct format for a ID
			if TF:
				pass
			else:
				TF = "missing_TF"	
				general_error.append( "missing_TF")
		except IndexError:
			general_error.append( "missing_TF_column")


		try:
		# getting TF's name
			cell = split_line[ col_cell].strip()

			# Check if a cell column is empty
			if cell:
				pass
			else:
				cell = "missing_cell"
				general_error.append( "missing_cell")	
		except IndexError:
			general_error.append( "missing SRP column")


		try:
		# getting serie's SRP
			SRP = split_line[ col_SRP].strip()

			# Check if correct SRP format
			if check_SRP_format( SRP):
				pass
			else:
				SRP = "wrong_SRP_format"
				general_error.append( "wrong_SRP_format")
		except IndexError:
			general_error.append( "missing_SRP_column")


		try:
		# getting ChIP GSM(s)
			list_chip_GSM = split_line[ col_GSM_chip].strip().split(";")

			for GSM in list_chip_GSM:
				# Check if correct format for a ID
				if check_sampleID_format( GSM):
					pass
				else:
					list_chip_GSM = "wrong_chip_GSM_format"	
					general_error.append( "wrong_chip_GSM_format")
		except IndexError:
			general_error.append( "missing_chip_GSM_column")


		try:
		# getting ChIP GSM(s)
			list_control_GSM = split_line[ col_GSM_control].strip().split(";")

			for GSM in list_control_GSM:
				# Check if correct format for a ID
				if check_sampleID_format( GSM):
					pass
				else:
					list_chip_GSM = "wrong_chip_GSM_format"	

					general_error.append( "wrong_control_GSM_format")
		except IndexError:
			general_error.append( "no_control")







		# error if missing column


		# print( ID_serie)
		# print( TF)
		# print( cell)
		# print( SRP)
		# print( list_chip_GSM)
		# print( list_control_GSM)
		# print( general_error)


		

		list_line_chip_ena = []
		list_line_control_ena = []

		# Only continue if no error
		if not general_error or general_error == ['no_control']:
			url_ena = url_experiment_srt +  SRP + url_experiment_end
			#print( url_ena)

			# getting info from ena
			ena_resp = requests.get( url_ena)
			if ena_resp.status_code != 200:
			    # This means something went wrong.
			    raise ApiError('GET /tasks/ {}'.format( ena_resp.status_code))

			# formating info
			list_line_ena_resp = ena_resp.text.split( '\n')

			iter_line_ena_resp = iter(list_line_ena_resp)

			#ignore header
			next(iter_line_ena_resp)
			# parse line
			for line_ena_resp in iter_line_ena_resp:

				list_info_ena = line_ena_resp.split( '\t')
				try:
					# add only relevant line
					if list_info_ena[ 4] in list_chip_GSM:
						list_line_chip_ena.append( list_info_ena)
					if list_info_ena[ 4] in list_control_GSM:
						list_line_control_ena.append( list_info_ena)

				except IndexError:
					general_error.append( "missing_gsm_ena")



			# Creating replicat number for each file
			list_chip_exp_ena = numbering_replicat( list_chip_GSM, list_line_chip_ena) 
			list_control_exp_ena = numbering_replicat( list_control_GSM, list_line_control_ena)



			# path outfile
			proposed_name_experiment = ".".join( [ ID_serie, TF, cell])
			proposed_name_experiment_control = ".".join( [ ID_serie, "input", cell])
			path_outdir_experiment = path_outdir + "/" + proposed_name_experiment
			path_log = path_outdir_experiment + "/log"
			path_outfile_experiment = path_outdir_experiment + "/" + proposed_name_experiment + "_summary.tab"
			# Creating outidir
			if not os.path.exists( path_outdir_experiment):
				os.makedirs( path_outdir_experiment)
			# Creating logdir
			if not os.path.exists( path_log):
				os.makedirs( path_log)

			outfile_experiment = open( path_outfile_experiment, 'w')

			line_header = "\t".join( ["filename", "info", "library_type", "isControl", "url", "md5"]) + "\n"

			outfile_experiment.write( line_header)

			# Creating list containing one file by line
			list_experiment_all_info = []
			for sample in list_chip_exp_ena:
				# dealing with single-end sample

				list_sample_info = []
				# Adding name
				list_sample_info.append( proposed_name_experiment + "." + sample[ 0] + "-" + sample[ 1])
				# Adding experiment title
				list_sample_info.append( sample[ 10])
				# Adding libray type
				list_sample_info.append( sample[ 7])
				# Adding control_type
				list_sample_info.append( "0")
				# # Adding SRR
				# list_sample_info.append( sample[ 5])
				# Adding download link
				list_sample_info.append( sample[ 8])
				# Adding md5
				list_sample_info.append( sample[ 9])

				outfile_experiment.write( "\t".join( list_sample_info) + "\n")




			# Adding control file
			for sample in list_control_exp_ena:
				# dealing with single-end sample
				list_sample_info = []
				# Adding name
				list_sample_info.append( proposed_name_experiment_control + "." + sample[ 0] + "-" + sample[ 1])
				# Adding experiment title
				list_sample_info.append( sample[ 10])
				# Adding libray type
				list_sample_info.append( sample[ 7])
				# Adding control_type
				list_sample_info.append( "1")
				# # Adding SRR
				# list_sample_info.append( sample[ 5])
				# Adding download link
				list_sample_info.append( sample[ 8])
				# Adding md5
				list_sample_info.append( sample[ 9])

				outfile_experiment.write( "\t".join( list_sample_info) + "\n")

	


			outfile_experiment.close()
