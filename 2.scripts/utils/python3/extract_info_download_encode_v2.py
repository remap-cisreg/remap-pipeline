#!/usr/bin/env


'''
From ENCODE result search file to necessary files to run the workflow
input: tab file
output: one file by experiment (line) containing info read by workflow
'''
__author__ = "Jeanne Chèneby"
__version__ = "0.2"
__status__ = "Production"

import argparse
import os
import pprint
import json
import urllib.request
import sys

if __name__ == "__main__":

	#================================================================#
	#                        Working path_outdir		             #
	#================================================================#

	parser = argparse.ArgumentParser(
        description="Create a summary file per experiment line containing all info necessary to run the remap workflow from ENCODE data",
        epilog="Jeanne Cheneby")
	parser.add_argument("-wd", "--working_directory",
                        help="Working derictory where are folders 1.metadata, 2.script, ...",
                        action="store",
                        required=True)
	parser.add_argument("-f", "--annotation_file", help="path to ENCODE search results file containing all annotation of ChIP-seq",
                        action="store",
                        required=True)
	args = parser.parse_args()






	# setting project working directory
	# path_working_dir = "/gpfs/tagc/home/cheneby/latest_remap"
	path_working_dir = args.working_directory
	os.chdir( path_working_dir)

	# setting TF file
	# path_file_remap_metadata = "1.metadata/TF_catalogue_run_20_01_17_FINAL.tab"
	path_file_encode_metadata = args.annotation_file

	path_outdir = "4.preprocessing"
	path_outdir_error = "E4.preprocessing"

	path_outdir_report = "E2.report"
	# Creating logdir
	if not os.path.exists( path_outdir_report):
		os.makedirs( path_outdir_report)


	#================================================================#
	#                        Start of script			             #
	#================================================================#

	file_encode_metadata = open( path_file_encode_metadata, 'r')
	# Dictionary where keys are experiments and vlue a list of file info composing it
	dict_experiment = {}
	final_dict_experiement = {}
	nb_line = 0
	list_experiment_ID = []

	iter_line_ena_resp = iter( file_encode_metadata)
	#ignore header
	next( file_encode_metadata)

	for line in file_encode_metadata:
		nb_line += 1



		split_line = line.split( "\t")

		list_experiment_ID.append( split_line[ 3])

	list_experiment_ID = list( set( list_experiment_ID))


	# Parse all experiment
	for experiment_ID in list_experiment_ID:

		# Loop variables
		dict_replicat = {}
		dict_replicat_single = {}
		dict_replicat_paired = {}
		temp_dict_experiment = {}
		list_control = []
		list_error = []


		# getting info on encode
		test_url = "https://www.encodeproject.org/experiments/" + experiment_ID + "/?format=json"
		with urllib.request.urlopen( test_url) as url:
			dict_exp_raw = json.loads(url.read().decode())
		list_info_files = dict_exp_raw[ "files"]



		# get all info from json
		for dict_info_file in list_info_files:
			if dict_info_file[ "file_type"] == "fastq":
				# get necessasy info
				cell_line =  dict_info_file[ "replicate"][ "experiment"][ "biosample_term_name"].replace( " ", "_")
				TF = dict_info_file[ "replicate"][ "experiment"][ "target"][ "label"]
				biological_replicate_number = dict_info_file[ "replicate"][ "biological_replicate_number"]
				technical_replicate_number = dict_info_file[ "replicate"][ "technical_replicate_number"]
				sample_ID = dict_info_file[ "title"]
				url_info = "https://www.encodeproject.org/files/" + sample_ID
				url_download = "https://www.encodeproject.org/" + dict_info_file[ "href"]
				md5 = dict_info_file[ "md5sum"]

				proposed_name_experiment = ".".join( [ experiment_ID, TF, cell_line])
				proposed_name_replicate = proposed_name_experiment + \
											"." + str( biological_replicate_number) + \
											"-" + str( technical_replicate_number)
				try:
					list_control = list_control + dict_info_file[ "controlled_by"]
				except KeyError:
					sys.stderr.write( proposed_name_replicate + " has no control but it's okay, that's life\n")





				# If sample is paired end
				if dict_info_file[ "run_type"] == "paired-ended":
					run_type = "PAIRED"
					paired_with =  dict_info_file[ "paired_with"].split( "/")[ -2]

					list_all_info = [
										sample_ID,
										url_info,
										run_type,
										"0",
										url_download,
										md5,
										biological_replicate_number,
										technical_replicate_number,
										paired_with
									]

					if proposed_name_replicate in dict_replicat_paired:
						dict_replicat_paired[ proposed_name_replicate].append( list_all_info)
					else:
						dict_replicat_paired[ proposed_name_replicate] = [ list_all_info]
				# If sample is single end
				else:
					run_type = "SINGLE"
					paired_with =  None


					list_all_info = [
										sample_ID,
										url_info,
										run_type,
										"0",
										url_download,
										md5,
										biological_replicate_number,
										technical_replicate_number,
										paired_with
									]
					if proposed_name_replicate in dict_replicat_single:
						dict_replicat_single[ proposed_name_replicate].append( list_all_info)
					else:
						dict_replicat_single[ proposed_name_replicate] = [ list_all_info]


		list_control = list( set( list_control))

		# Formating info ReMap style
		list_exp_line = []
		for key in dict_replicat_paired:
			if len( dict_replicat_paired[ key]) == 2:
				temp_line =		[ 	key,
									dict_replicat_paired[ key][ 0][1],
									dict_replicat_paired[ key][ 0][ 2],
									dict_replicat_paired[ key][ 0][ 3],
									dict_replicat_paired[ key][ 0][ 4] + ";" + dict_replicat_paired[ key][ 1][ 4],
									dict_replicat_paired[ key][ 0][ 5] + ";" + dict_replicat_paired[ key][ 1][ 5]
								]
				list_exp_line.append( temp_line)
			else:
				count_paired_file = 0
				while count_paired_file < len( dict_replicat_paired[ key])//2:
					current_sample_ID = dict_replicat_paired[ key][ count_paired_file][ 0]
					current_buddy = dict_replicat_paired[ key][ count_paired_file][ -1]
					for element in  dict_replicat_paired[ key]:
						if element[ 0] == current_buddy:

							temp_line = [ 	key,
											dict_replicat_paired[ key][ count_paired_file][1],
											dict_replicat_paired[ key][ count_paired_file][ 2],
											dict_replicat_paired[ key][ count_paired_file][ 3],
											dict_replicat_paired[ key][ count_paired_file][ 4] + ";" + element[ 4],
											dict_replicat_paired[ key][ count_paired_file][ 5] + ";" + element[ 5]
										]
							list_exp_line.append( temp_line)
							list_error.append( key + " has at least two files with the same replicat number")
							count_paired_file +=1

					# print( current_sample_ID)
					# print( current_buddy)




		for key in dict_replicat_single:
			if len( dict_replicat_single[ key]) == 1:
				temp_line = [
								key,
								dict_replicat_single[ key][ 0][ 1],
								dict_replicat_single[ key][ 0][ 2],
								dict_replicat_single[ key][ 0][ 3],
								dict_replicat_single[ key][ 0][ 4],
								dict_replicat_single[ key][ 0][ 5]
				]
				list_exp_line.append( temp_line)

			else:
				list_error.append( key + " has at least two files with the same replicat number")


		# pp = pprint.PrettyPrinter(indent=4)
		# pp.pprint( list_exp_line)
		# # print( list_control)
		# print( )

		temp_dict_experiment[ 'list_exp'] = list_exp_line
		temp_dict_experiment[ 'control'] = list_control

		dict_experiment[ experiment_ID] = temp_dict_experiment



	# with open( '1.metadata/info_encode.json', 'w') as fp:
	# 	json.dump( dict_experiment, fp)


	# with open( '1.metadata/info_encode.json') as f:
	# 	dict_experiment = json.load(f)



	# Get info for control
	for key in dict_experiment:
		list_control = []
		dict_replicat_paired = {}
		final_outdir = path_outdir

		if len( dict_experiment[ key][ 'control']) > 0:
			for end_url_control in dict_experiment[ key][ 'control']:
				# print( control_ID)

				# get json from encode
				url_control = "https://www.encodeproject.org" + end_url_control + "?format=json"
				url_info = "https://www.encodeproject.org" + end_url_control
				with urllib.request.urlopen( url_control) as url:
					dict_control_file_raw = json.loads(url.read().decode())


				# Get all Info
				replicat_md5sum = dict_control_file_raw["md5sum"]
				run_type = dict_control_file_raw["run_type"]
				end_donwload_link = dict_control_file_raw["href"]
				biological_replicat = dict_control_file_raw['replicate']['biological_replicate_number']
				technical_replicat = dict_control_file_raw['replicate']['technical_replicate_number']
				cell_line = dict_control_file_raw['replicate']['experiment']['biosample_term_name'].replace( " ", "_")
				replicat_number = str( biological_replicat) + "-" + str( technical_replicat)
				proposed_name_replicate = key + ".INPUT." + cell_line + "." +	replicat_number
				download_link = "https://www.encodeproject.org" + end_donwload_link


				# If paired add to dictionary to deal with
				if run_type == 'paired-ended':
					run_type = 'PAIRED'
					list_all_info = [ proposed_name_replicate, url_info, run_type, '1', download_link, replicat_md5sum]

					if proposed_name_replicate in dict_replicat_paired:
						dict_replicat_paired[ proposed_name_replicate].append( list_all_info)
					else:
						list_all_info = [ 	proposed_name_replicate,
											url_info, run_type,
											'1',
											download_link,
											replicat_md5sum]
						dict_replicat_paired[ proposed_name_replicate] = [ list_all_info]
				# If single end just add to list
				else:
					run_type = 'SINGLE'
					list_all_info = [ proposed_name_replicate, url_info, run_type, '1', download_link, replicat_md5sum]
					list_control.append( list_all_info)

			for paired_key in dict_replicat_paired:

				if len( dict_replicat_paired[ paired_key]) == 2:
					list_all_info =		[ 	paired_key,
										dict_replicat_paired[ paired_key][ 0][1],
										dict_replicat_paired[ paired_key][ 0][ 2],
										dict_replicat_paired[ paired_key][ 0][ 3],
										dict_replicat_paired[ paired_key][ 0][ 4] + ";" + dict_replicat_paired[ paired_key][ 1][ 4],
										dict_replicat_paired[ paired_key][ 0][ 5] + ";" + dict_replicat_paired[ paired_key][ 1][ 5]
									]
					list_control.append( list_all_info)


				else:
					final_outdir = path_outdir_error
					sys.stderr.write( paired_key + " paired-end replicat is missing a file\n")


		final_dict_experiement[ key] = {}
		final_dict_experiement[ key][ 'list_exp'] = dict_experiment[ key][ 'list_exp']
		final_dict_experiement[ key][ 'list_control'] = list_control

		# print( key)

		try:
			proposed_name_experiment = dict_experiment[ key][ 'list_exp'][0][0].rsplit( ".", 1)[0]
		except IndexError:
			sys.stderr.write( key + " missing experiment file(s)\n")


		if proposed_name_experiment:

			path_outdir_experiment = os.path.join( final_outdir, proposed_name_experiment)
			path_log = os.path.join( path_outdir_experiment, "log")


			# print( proposed_name_experiment)
			# print( path_outdir_experiment)
			# print()

			# Creating outidir
			if not os.path.exists( path_outdir_experiment):
				os.makedirs( path_outdir_experiment)
			# Creating logdir
			if not os.path.exists( path_log):
				os.makedirs( path_log)

			path_outfile_experiment = os.path.join( path_outdir_experiment, proposed_name_experiment + "_summary.tab")
			outfile_experiment = open( path_outfile_experiment, 'w')

			outfile_experiment.write( "\t".join( [ "filename", "info", "library_type", "isControl", "url", "md5"]) + "\n")
			for line_rep in dict_experiment[ key][ 'list_exp']:
				outfile_experiment.write( "\t".join( line_rep) + "\n")
			for line_rep_control in list_control:
				outfile_experiment.write( "\t".join( line_rep_control) + "\n")

			outfile_experiment.close()




	# pp = pprint.PrettyPrinter(indent=4)
	# pp.pprint( final_dict_experiement)

	# json = json.dumps( final_dict_experiement)
	# outfile_dict_control = open( 'test.json', "w")
	# outfile_dict_control.write(json)
	# outfile_dict_control.close()
