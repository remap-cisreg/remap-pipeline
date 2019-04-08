#!/usr/bin/env


'''
From TF annoation file to necessary files to run the workflow
input: tab file
output: one file by experiment (line) containing info read by workflow
'''
__author__ = "Jeanne Chèneby"
__version__ = "0.2"
__status__ = "Production"

import os
import sys
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

# Format and check info from ena into remap tab format
def process_list_experiment( list_exp, name_experiment):
	# Creating list containing one file by line
	list_experiment_all_info = []
	list_all_error = []
	list_all_sample_info = []


	# Check if list of experiment contain chip or control
	if "input" in name_experiment:
		chip_type = "1"
	else:
		chip_type = "0"

	# Parse all sample of experiment
	for sample in list_exp:

		# Initialize info for sample
		list_sample_info = []
		list_error = []



		if sample[ 7] == 'SINGLE':
			final_name = sample[ 5] + "." + name_experiment.split( ".", 1)[1]
		else:
			final_name = sample[ 5] + "r1-" + sample[ 5] + "r2." + name_experiment.split( ".", 1)[1]

		# Adding name
		# final_name = name_experiment + "." + sample[ 0] + "-" + sample[ 1]

		list_sample_info, list_error = check_column( final_name, name_experiment + " missing filename column", list_sample_info, list_error)
		# if final_name:
		# 	list_sample_info.append( final_name)
		# else:
		# 	current_error = name_experiment + " missing filename column"
		# 	list_error.append( current_error)
		# 	sys.stderr.write( current_error + "\n")


		# Adding experiment title
		final_exp_title = sample[ 10]
		list_sample_info.append( final_exp_title)


		# Adding libray type
		final_library = sample[ 7]
		list_sample_info, list_error = check_column( final_library, final_name + " missing libray type column", list_sample_info, list_error)
		# if final_library:
		# 	list_sample_info.append( final_library)
		# else:
		# 	current_error = final_name + " missing libray type column"
		# 	list_error.append( current_error)
		# 	sys.stderr.write( current_error + "\n")

		# Adding control_type
		list_sample_info.append( chip_type)

		# # Adding SRR
		# list_sample_info.append( sample[ 5])

		# Adding download link
		final_link = sample[ 8]
		problem_link_paired = False
		if final_link:
			# Check if too many link (Only works for paired-end)
			if final_library == 'PAIRED':
				# If there is more than 2 link
				if len( final_link.split( ";")) > 2:
					problem_link_paired = True
					list_split_finale_link = final_link.split( ";")
					list_temp_final_link = []
					list_index = []
					# Parse all link
					for i in range( len( list_split_finale_link)):
						# _1 or _2 add link to list
						if "_1" in list_split_finale_link[ i]:
							list_temp_final_link.append( "ftp://" + list_split_finale_link[ i])
							list_index.append( i)
						elif "_2" in list_split_finale_link[ i]:
							list_temp_final_link.append( "ftp://" + list_split_finale_link[ i])
							list_index.append( i)
						# else don't add link in list
						else:
							current_error = final_name + " has too many donwload links: " + final_link + " and md5 is removed"
							list_error.append( current_error)
							sys.stderr.write( current_error + "\n")
					list_sample_info.append( ";".join( list_temp_final_link))
				# If there is only two link
				elif len( final_link.split( ";")) == 2:
					split_link = final_link.split( ";")
					list_sample_info.append( "ftp://" + split_link[0] + ";ftp://" + split_link[ -1])
				else:
					current_error = "Paired sample " + final_name + " is missing donwload links: " + final_link
					list_error.append( current_error)
					sys.stderr.write( current_error + "\n")

			elif final_library == 'SINGLE':
				if len( final_link.split( ";")) == 1:
					list_sample_info.append(  "ftp://" + final_link)
				else:
					current_error = "Single sample " + final_name + " has too many donwload links: " + final_link
					list_error.append( current_error)
					sys.stderr.write( current_error + "\n")

		else:
			current_error = final_name + " missing download link column"
			list_error.append( current_error)
			sys.stderr.write( current_error + "\n")


		# Adding md5
		final_md5 = sample[ 9]
		if final_md5:
			if problem_link_paired:
				list_split_final_md5 = final_md5.split( ";")
				list_temp_md5 = []
				for current_index in list_index:
					list_temp_md5.append( list_split_final_md5[ current_index])
				list_sample_info.append( ";".join( list_temp_md5))
			else:
				list_sample_info.append( final_md5)
		else:
			current_error = final_name + " missing md5 column"
			list_error.append( current_error)
			sys.stderr.write( current_error + "\n")
		list_all_sample_info.append( list_sample_info)
		list_all_error.append( list_error)

	return list_all_sample_info, list_all_error

def check_column( info_column, error_message, list_current_info, list_current_error):

	if info_column:
		list_current_info.append( info_column)
	else:
		list_current_error.append( error_message)
		sys.stderr.write( error_message + "\n")
	return list_current_info, list_current_error

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

	path_outdir = os.path.join( "1.metadata", "tab")
	path_outdir_error = os.path.join( "1.metadata", "tab_error")

	path_outdir_report = "E2.report"
	# Creating logdir
	if not os.path.exists( path_outdir_report):
		os.makedirs( path_outdir_report)

	#================================================================#
	#                        Necessary variable			             #
	#================================================================#

	general_error = []





	# setting column of interest (start from 0)
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
						"sample_title",					# Name of sample (non standart)
						"run_alias"						# Alternative GSM

					]


	info_ena =  ",".join( list_info_ena)

		# final url end
	url_experiment_end = url_experiment_temp_end_1 + info_ena + url_experiment_temp_end_2



	#================================================================#
	#                        Start of script			             #
	#================================================================#



	# get relevant info
	file_remap_metadata = open( path_file_remap_metadata, 'r')
	nb_line = 0
	for line in file_remap_metadata:
		line_error = []
		nb_line += 1

		split_line = line.strip().split( "\t")

		try:
		# getting serie's ID
			ID_serie = split_line[ col_ID_serie].strip()

			# Check if correct format for a ID
			if check_ID_format( ID_serie):
				pass
			else:
				ID_serie = "wrong_ID_serie_format"
				current_error = "line" + str( nb_line) + ":wrong_ID_serie_format"
				line_error.append( current_error)
				sys.stderr.write( current_error + "\n")

		except IndexError:
			current_error = "line" + str( nb_line) + ":missing_ID_column"
			line_error.append( current_error)
			sys.stderr.write( current_error)


		try:
		# getting TF's name
			TF = split_line[ col_TF].strip()

			# Check if correct format for a ID
			if TF:
				pass
			else:
				TF = "missing_TF"
				current_error = "line" + str( nb_line) + ":missing_TF"
				line_error.append( current_error)
				sys.stderr.write( current_error + "\n")
		except IndexError:
			current_error = "line" + str( nb_line) + ":missing_TF_column"
			line_error.append( current_error)
			sys.stderr.write( current_error + "\n")


		try:
		# getting cell's name
			cell = split_line[ col_cell].strip()

			# Check if a cell column is empty
			if cell:
				pass
			else:
				cell = "missing_cell_type"
				current_error = "line" + str( nb_line) + ":missing_cell_type"
				line_error.append( current_error)
				sys.stderr.write( current_error  + "\n")
		except IndexError:
			current_error = "line" + str( nb_line) + ":missing_cell_type column"
			line_error.append( current_error)
			sys.stderr.write( current_error + "\n")


		try:
		# getting serie's SRP
			SRP = split_line[ col_SRP].strip()

			# Check if correct SRP format
			if check_SRP_format( SRP):
				pass
			else:
				SRP = "wrong_SRP_format"
				current_error = "line" + str( nb_line) + ":missing_wrong_SRP_format"
				line_error.append( current_error)
				sys.stderr.write( current_error + "\n")
		except IndexError:
			current_error = "line" + str( nb_line) + ":missing_SRP_column"
			line_error.append( current_error)
			sys.stderr.write( current_error + "\n")

		list_chip_GSM = []
		try:
		# getting ChIP GSM(s)
			list_chip_GSM = split_line[ col_GSM_chip].strip().split(";")

			for GSM in list_chip_GSM:
				# Check if correct format for a ID
				if check_sampleID_format( GSM):
					pass
				else:
					list_chip_GSM = "wrong_chip_GSM_format"
					current_error = "line" + str( nb_line) + ":wrong_chip_GSM_format"
					line_error.append( current_error)
					sys.stderr.write( current_error + "\n")

		except IndexError:
			current_error = "line" + str( nb_line) + ":missing_chip_GSM_column"
			line_error.append( current_error)
			sys.stderr.write( current_error + "\n")

		list_control_GSM = []
		try:
		# getting control GSM(s)
			list_control_GSM = split_line[ col_GSM_control].strip().split(";")

			for GSM in list_control_GSM:
				control_are_you_here = "Boo !"
				# Check if correct format for a ID
				if check_sampleID_format( GSM):
					pass
				elif GSM == "":
					control_are_you_here = "nop"
				else:
					list_chip_GSM = "wrong_chip_GSM_format"

					current_error = "line" + str( nb_line) + ":wrong_control_GSM_format or no control"
					line_error.append( current_error)
					sys.stderr.write( current_error + "\n")
		except IndexError:
			# line_error.append( "no_control")
			control_are_you_here = "nop"

		# proposed name
		proposed_name_experiment = ".".join( [ ID_serie, TF, cell])
		proposed_name_experiment_control = ".".join( [ ID_serie, "input", cell])


		list_line_chip_ena = []
		list_line_control_ena = []
		list_line_chip_ena_alt = []
		list_line_control_ena_alt = []
		list_error_sample_chip = []
		list_error_sample_control = []

		# Only continue if no error
		if  not line_error:
			url_ena = url_experiment_srt +  SRP + url_experiment_end
			# proposed name
			proposed_name_experiment = ".".join( [ ID_serie, TF, cell])
			proposed_name_experiment_control = ".".join( [ ID_serie, "input", cell])


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
			#List of GSM from files found in ENA
			list_found_GSM = []
			list_found_GSM_alt = []
			# parse line
			for line_ena_resp in iter_line_ena_resp:

				list_info_ena = line_ena_resp.split( '\t')

				try:
					# add only relevant line
					if list_info_ena[ 4] in list_chip_GSM:
						list_line_chip_ena.append( list_info_ena)
						list_found_GSM.append( list_info_ena[ 4])
					if list_info_ena[ 4] in list_control_GSM:
						list_line_control_ena.append( list_info_ena)
						list_found_GSM.append( list_info_ena[ 4])

					maybe_gsm = list_info_ena[ 9].split( "_")[ 0].strip()
					if maybe_gsm in list_chip_GSM:
						list_info_ena[ 4] = maybe_gsm
						list_line_chip_ena_alt.append( list_info_ena)
						list_found_GSM_alt.append( maybe_gsm)
					if maybe_gsm in list_control_GSM:
						list_info_ena[ 4] = maybe_gsm
						list_line_control_ena_alt.append( list_info_ena)
						list_found_GSM_alt.append( maybe_gsm)

				except IndexError:
					pass
					# line_error.append( "missing_gsm_ena for line " + str( nb_line))



			# Check if all GSM from files were found in ENA
			list_all_chip = list_chip_GSM + list_control_GSM
			list_all_chip = list( filter( None, list_all_chip))

			check_GSM_ena =  all(elem in list_found_GSM  for elem in list_all_chip)
			if not check_GSM_ena:
				check_GSM_ena =  all(elem in list_found_GSM_alt  for elem in list_all_chip)

			# Only process if ALL GSM are found
			if check_GSM_ena:

				list_chip_exp_ena = []
				list_control_exp_ena = []

				# Only process info if GEO info are found on ENA
				if list_line_chip_ena:

					# Creating replicat number for each file
					list_chip_exp_ena = numbering_replicat( list_chip_GSM, list_line_chip_ena)
					list_control_exp_ena = numbering_replicat( list_control_GSM, list_line_control_ena)






					# Get all info from ena
					list_sample_chip, list_error_sample_chip = process_list_experiment( list_chip_exp_ena, proposed_name_experiment)
					list_sample_control, list_error_sample_control = process_list_experiment( list_control_exp_ena, proposed_name_experiment_control)


					list_all_sample = list_sample_chip + list_sample_control
					list_all_sample_error = list_error_sample_chip +  list_error_sample_control


					error_in_exp = False
					for elm in list_all_sample_error:
						if len( elm) > 0:
							error_in_exp = True

					# If there is an error in the experiment
					if error_in_exp:
						path_outdir_experiment = os.path.join( path_outdir_error)
						# Create outdir
						path_log = os.path.join( path_outdir_experiment, "log")
						# Creating outidir
						if not os.path.exists( path_outdir_experiment):
							os.makedirs( path_outdir_experiment)


						path_outfile_experiment = os.path.join( path_outdir_experiment, proposed_name_experiment + "_ena_parsing.err")
						outfile_experiment = open( path_outfile_experiment, 'w')

						for all_sample_error in list_all_sample_error:
							for sample_error in all_sample_error:
								if sample_error:
									outfile_experiment.write( sample_error + "\n")
						outfile_experiment.close()
					# If there is no error in the experiment
					else:
						path_outdir_experiment = os.path.join( path_outdir)

						# Create outdir
						path_log = os.path.join( path_outdir_experiment, "log")
						# Creating outidir
						if not os.path.exists( path_outdir_experiment):
							os.makedirs( path_outdir_experiment)


						# Writting results
						path_outfile_experiment = os.path.join( path_outdir_experiment, proposed_name_experiment + "_summary.tab")
						outfile_experiment = open( path_outfile_experiment, 'w')

						line_header = "\t".join( ["filename", "info", "library_type", "isControl", "url", "md5"]) + "\n"
						outfile_experiment.write( line_header)

						for sample in list_all_sample:
							outfile_experiment.write( "\t".join( sample) + "\n")
						outfile_experiment.close()

				# trying to match GSM with run_alias
				elif list_line_chip_ena_alt:
					# Creating replicat number for each file
					list_chip_exp_ena = numbering_replicat( list_chip_GSM, list_line_chip_ena_alt)
					list_control_exp_ena = numbering_replicat( list_control_GSM, list_line_control_ena_alt)


					# Get all info from ena
					list_sample_chip, list_error_sample_chip = process_list_experiment( list_chip_exp_ena, proposed_name_experiment)
					list_sample_control, list_error_sample_control = process_list_experiment( list_control_exp_ena, proposed_name_experiment_control)

					list_all_sample = list_sample_chip + list_sample_control
					list_all_sample_error = list_error_sample_chip +  list_error_sample_control


					error_in_exp = False
					for elm in list_all_sample_error:
						if len( elm) > 0:
							error_in_exp = True

					# If there is an error in the experiment
					if error_in_exp:
						path_outdir_experiment = os.path.join( path_outdir_error)

						# Create outdir
						path_log = os.path.join( path_outdir_experiment, "log")
						# Creating outidir
						if not os.path.exists( path_outdir_experiment):
							os.makedirs( path_outdir_experiment)


						path_outfile_experiment = os.path.join( path_outdir_experiment, proposed_name_experiment + "_ena_parsing.err")
						outfile_experiment = open( path_outfile_experiment, 'w')

						for all_sample_error in list_all_sample_error:
							for sample_error in all_sample_error:
								if sample_error:
									outfile_experiment.write( sample_error + "\n")
						outfile_experiment.close()
					# If there is no error in the experiment
					else:
						path_outdir_experiment = os.path.join( path_outdir)

						# Create outdir
						path_log = os.path.join( path_outdir_experiment, "log")
						# Creating outidir
						if not os.path.exists( path_outdir_experiment):
							os.makedirs( path_outdir_experiment)


						# Writting results
						path_outfile_experiment = os.path.join( path_outdir_experiment, proposed_name_experiment + "_summary.tab")
						outfile_experiment = open( path_outfile_experiment, 'w')

						line_header = "\t".join( ["filename", "info", "library_type", "isControl", "url", "md5"]) + "\n"
						outfile_experiment.write( line_header)

						for sample in list_all_sample:
							outfile_experiment.write( "\t".join( sample) + "\n")
						outfile_experiment.close()

					# If there are parsing error (no GSM found)
				else:

					# Create outdir
					path_outdir_experiment = os.path.join( path_outdir_error)
					path_log = os.path.join( path_outdir_experiment, "log")

					# Creating outidir
					if not os.path.exists( path_outdir_experiment):
						os.makedirs( path_outdir_experiment)


					# writting error message
					path_outfile_experiment = os.path.join( path_outdir_experiment, proposed_name_experiment + "no_GSM__ena.err")

					outfile_experiment = open( path_outfile_experiment, 'w')
					outfile_experiment.write( "Can't find GSM in ENA, check: " + url_ena)
					outfile_experiment.close()

			# If not all GSM were found
			else:

				list_missing_gsm = list( set( list_all_chip) - set( list_found_GSM_alt))
				# Create outdir
				path_outdir_experiment = os.path.join( path_outdir_error)
				path_log = os.path.join( path_outdir_experiment, "log")

				# Creating outidir
				if not os.path.exists( path_outdir_experiment):
					os.makedirs( path_outdir_experiment)


				# writting error message
				path_outfile_experiment = os.path.join( path_outdir_experiment, proposed_name_experiment + "_missing_GSM_ena.err")

				outfile_experiment = open( path_outfile_experiment, 'w')
				for missing_gsm in list_missing_gsm:
					outfile_experiment.write( "Missing " + missing_gsm + " in ENA, check: " + url_ena + "\n")
				outfile_experiment.close()

		# If there were error processing annotation file
		else:
			# Create outdir
			path_outdir_experiment = os.path.join( path_outdir_error)
			path_log = os.path.join( path_outdir_experiment, "log")
			# Creating outidir
			if not os.path.exists( path_outdir_experiment):
				os.makedirs( path_outdir_experiment)

			path_outfile_experiment = os.path.join( path_outdir_experiment, proposed_name_experiment + "_file_parsing.err")
			outfile_experiment = open( path_outfile_experiment, 'w')

			for elm in line_error:
				outfile_experiment.write( elm + "\n")
			outfile_experiment.close()
