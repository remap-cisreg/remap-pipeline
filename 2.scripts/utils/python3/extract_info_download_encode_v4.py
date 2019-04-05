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





def join_paired_info( list_info_tech_rep, index_info):

	paired_current_info = ";".join(
							[
								list_info_tech_rep[ 0][ index_info],
								list_info_tech_rep[ 1][ index_info]
							]
						   )

	return paired_current_info




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

	path_outdir = os.path.join( "1.metadata", "tab")
	path_outdir_error = os.path.join( "1.metadata", "tab_error")
	path_outdir_duplicate = "D4.preprocessing"

	path_outdir_report = "E2.report"
	# Creating logdir
	if not os.path.exists( path_outdir_report):
		os.makedirs( path_outdir_report)


	#================================================================#
	#                        globale variables			             #
	#================================================================#

	dict_experiment = {}
	column_file_accession = 0
	column_experiment_accession = 3
	column_biosample_term_id = 5 #EFO/BTO/...
	column_biosample_term_name = 6
	column_biosample_treatment = 9
	column_biosample_treatment_amount = 10
	column_biosample_treatment_duration = 11
	column_experiment_target = 12
	column_biological_replicate = 24
	column_technical_replicate = 25
	column_run_type = 28
	column_pair_with = 30
	column_md5 = 34
	column_url_download = 36
	column_controlled_by = 39
	column_file_status = 40
	column_audit_warning = 42
	column_audit_internal_action = 43
	column_audit_not_compliant = 44


	#================================================================#
	#                        Start of script			             #
	#================================================================#


	# Dictionary where keys are experiments and vlue a list of file info composing it
	final_dict_experiment = {}
	dict_control_file = {}


	# Open metadata file
	file_encode_metadata = open( path_file_encode_metadata, 'r')
	iter_line_ena_resp = iter( file_encode_metadata)
	#ignore header
	next( file_encode_metadata)
	for line in file_encode_metadata:
		# print( line)

		split_line = line.split( "\t")

		file_accession = split_line[ column_file_accession]
		experiment_accession = split_line[ column_experiment_accession]



		# Creating a list with biosample and treatment information
		list_biosample_term_name = [ split_line[ column_biosample_term_name]]
		list_biosample_term_name.append( split_line[ column_biosample_treatment])
		list_biosample_term_name.append( split_line[ column_biosample_treatment_amount])
		list_biosample_term_name.append( split_line[ column_biosample_treatment_duration])
		biosample_name = "_".join( filter (None, list_biosample_term_name))

		experiment_target = split_line[ column_experiment_target]

		# Creating experiment name
		experiment_name = ".".join( [ experiment_accession, experiment_target, biosample_name])

		# Creating replicat name (one per line)
		biological_replicate = split_line[ column_biological_replicate]
		technical_replicate = split_line[ column_technical_replicate]
		biological_replicat_name = ".".join([ experiment_name, biological_replicate])

		run_type = split_line[ column_run_type]

		# If there are no information about this experiment
		if experiment_name not in final_dict_experiment:
			final_dict_experiment[ experiment_name] = {}
			final_dict_experiment[ experiment_name][ 'chip'] = {}
			final_dict_experiment[ experiment_name][ 'control'] = []
			final_dict_experiment[ experiment_name][ 'experiment_accession'] = experiment_accession
			final_dict_experiment[ experiment_name][ 'biosample_name'] = biosample_name

		# If this biological replicat was not already in the dictionary
		if biological_replicat_name not in final_dict_experiment[ experiment_name][ 'chip']:
			final_dict_experiment[ experiment_name][ 'chip'][ biological_replicat_name] = {}
			final_dict_experiment[ experiment_name][ 'chip'][ biological_replicat_name][ 'paired-end'] = {}
			final_dict_experiment[ experiment_name][ 'chip'][ biological_replicat_name][ 'single-end'] = []



		# Getting non run type related info for final file
		info = os.path.join( "https://www.encodeproject.org/files/", file_accession)
		isControl = "0"
		url_download = split_line[ column_url_download]
		md5 = split_line[ column_md5]

		# If this replicat is single end
		if run_type == 'single-ended':
			list_info = [ 	info,
							'SINGLE',
							isControl,
							url_download,
							md5
						]

			final_dict_experiment[ experiment_name][ 'chip'][ biological_replicat_name][ 'single-end'].append( list_info)
		# If this replicat is paired end
		else:

			list_info = [
							info,
							'PAIRED',
							isControl,
							url_download,
							md5
						]
			pair_with = split_line[ column_pair_with]

			# Check if buddy is already in the dictionary
			if pair_with in final_dict_experiment[ experiment_name][ 'chip'][ biological_replicat_name][ 'paired-end']:
				final_dict_experiment[ experiment_name][ 'chip'][ biological_replicat_name][ 'paired-end'][ pair_with].append( list_info)
			# If buddy is not in the dictionary add a key with file name
			else:
				final_dict_experiment[ experiment_name][ 'chip'][ biological_replicat_name][ 'paired-end'][ file_accession] = [ list_info]

		# Listing control in dictionary
		controlled_by = split_line[ column_controlled_by]
		if controlled_by:
			list_controlled_by = controlled_by.split( ",")

			list_clean_controlled_by = []

			for url_controlled_by in list_controlled_by:
				clean_url_controlled_by = url_controlled_by.strip()
				list_clean_controlled_by.append( clean_url_controlled_by)

				if clean_url_controlled_by not in dict_control_file:
					dict_control_file[ clean_url_controlled_by] = {}

			final_dict_experiment[ experiment_name][ 'control'] = list( set( final_dict_experiment[ experiment_name][ 'control'] + list_clean_controlled_by))






	# pprint.pprint( final_dict_experiment)
	# pprint.pprint( dict_control)

	# Dealing with control
	dict_control_info = {}

	for url_control in dict_control_file:
		dict_control_info[ url_control] = {}
		# print( url_control)

		# get json from encode
		url_control_encode = "https://www.encodeproject.org" + url_control.strip()
		url_control_encode_json = url_control_encode + "?format=json"
		# print( url_control_encode_json)
		with urllib.request.urlopen( url_control_encode_json) as url:
			dict_control_file_raw = json.loads( url.read().decode())
		# print( "ok")
		# print( )

		# Get all Info

		dict_control_info[ url_control][ 'biological_replicat'] = dict_control_file_raw['replicate']['biological_replicate_number']
		dict_control_info[ url_control][ 'run_type'] = dict_control_file_raw["run_type"]
		dict_control_info[ url_control][ 'url_download']= "https://www.encodeproject.org" + dict_control_file_raw["href"]
		dict_control_info[ url_control][ 'md5'] = dict_control_file_raw["md5sum"]
		dict_control_info[ url_control][ 'url_control_encode'] = url_control_encode


		# If this replicat is single end
		if dict_control_info[ url_control][ 'run_type'] == 'single-ended':
			dict_control_info[ url_control][ 'pair_with'] = dict_control_file_raw[ 'paired_with'] = ''
		else:
			dict_control_info[ url_control][ 'pair_with'] = dict_control_file_raw[ 'paired_with']


	# with open( '1.metadata/info_control_encode.json', 'w') as fp:
	# 	json.dump( dict_control_info, fp)
	#
	# with open( '1.metadata/info_control_encode.json') as f:
	# 	dict_control_info = json.load(f)


	dict_paired_control = {}
	dict_paired_control[ 'SINGLE'] = {}
	dict_paired_control[ 'PAIRED'] = {}
	for url_control in dict_control_info:


		if dict_control_info[ url_control][ 'run_type'] == 'paired-ended':

			if dict_control_info[ url_control][ 'pair_with'] in dict_paired_control[ 'PAIRED']:
				dict_paired_control[ 'PAIRED'][ dict_control_info[ url_control][ 'pair_with']][ '2'] = dict_control_info[ url_control]
			else:
				dict_paired_control[ 'PAIRED'][ url_control] = {}
				dict_paired_control[ 'PAIRED'][ url_control][ '1'] = dict_control_info[ url_control]
				dict_paired_control[ 'PAIRED'][ url_control][ 'already_used'] = False

		elif dict_control_info[ url_control][ 'run_type'] == 'single-ended':
			dict_paired_control[ 'SINGLE'][ url_control] = {}
			dict_paired_control[ 'SINGLE'][ url_control][ 'info'] = dict_control_info[ url_control]
			dict_paired_control[ 'SINGLE'][ url_control][ 'already_used'] = False




	# pprint.pprint( dict_paired_control)
	# Formating data
	dict_used_control = {}

	for experiment_name in final_dict_experiment:
		used_control = False

		list_new_control = final_dict_experiment[ experiment_name][ 'control'].copy()

		# Get list of uniq control
		for url_control in final_dict_experiment[ experiment_name][ 'control']:

			if url_control in dict_paired_control[ 'PAIRED']:
				try:
					list_new_control.remove( dict_paired_control[ 'PAIRED'][ url_control][ '1'][ 'pair_with'])
				except ValueError:
					# sys.stderr.write( dict_paired_control[ 'PAIRED'][ url_control][ '1'][ 'pair_with'])
					# sys.stderr.write( " ,".join( list_new_control))
					# sys.stderr.write( " ,".join( final_dict_experiment[ experiment_name][ 'control']))
					pass
			else:
				# sys.stderr.write( " ".join( ['I hope you have a buddy', url_control, '...\n']))
				pass






		error_experiment = False

		list_line = [
					 [
						"filename",
						"info",
						"library_type",
						"isControl",
						"url",
						"md5"
					 ]
					]

		# Formating control dictionary to match chip dictionary structure
		# For all control of this experiment get info
		# technical_replicate_number = 1
		control_true_name = ".".join( [
									"INPUT",
									final_dict_experiment[ experiment_name][ 'biosample_name']
								 ]
								)

		for url_controlled_by in list_new_control:

			control_name = ".".join( [
										url_controlled_by.split( "/")[ 2],
										"INPUT",
										final_dict_experiment[ experiment_name][ 'biosample_name']
									 ]
									)
			# If control is single-ended
			if url_controlled_by in dict_paired_control[ 'SINGLE']:

				# Get uniq name and name in experiment
				file_id = dict_control_info[ url_controlled_by][ 'url_control_encode'].rsplit( "/", 2)[ 1]

				biological_replicat_name = ".".join( [ control_name, str( dict_paired_control[ 'SINGLE'][ url_controlled_by][ 'info'][ 'biological_replicat'])])
				true_biological_replicat_name = ".".join( [ ".".join( [ file_id, control_true_name]), str( dict_paired_control[ 'SINGLE'][ url_controlled_by][ 'info'][ 'biological_replicat'])])

				# If this biological replicat was not already in the dictionary
				list_info = [
								# "-".join( [ true_biological_replicat_name, str( technical_replicate_number)]),
								true_biological_replicat_name,
								dict_control_info[ url_controlled_by][ 'url_control_encode'],
								'SINGLE',
								"1",
								dict_control_info[ url_controlled_by][ 'url_download'],
								dict_control_info[ url_controlled_by][ 'md5']
							]

				list_line.append( list_info)
				# technical_replicate_number += 1

				# Check if key already used in this batch
				if biological_replicat_name not in dict_used_control:
					dict_used_control[ biological_replicat_name] = experiment_name
					dict_used_control[ experiment_name] = {}
					dict_used_control[ experiment_name][ 'first_unique'] = {}
					dict_used_control[ experiment_name][ 'dict_duplicate'] = {}
					dict_used_control[ experiment_name][ 'first_unique'][ experiment_name] = biological_replicat_name
				else:
					used_control = True
					# sys.stderr.write( biological_replicat_name + " already used in " + dict_used_control[ biological_replicat_name] + "\n")

			# if control ispaired-end
			elif url_controlled_by in dict_paired_control[ 'PAIRED']:

				# pprint.pprint( dict_paired_control[ 'PAIRED'][ url_controlled_by])
				try:
					file_id = "-".join( [ 	dict_paired_control[ 'PAIRED'][ url_controlled_by][ '1'][ 'url_control_encode'].rsplit( "/", 2)[ 1],
											dict_paired_control[ 'PAIRED'][ url_controlled_by][ '2'][ 'url_control_encode'].rsplit( "/", 2)[ 1]
										])

					# Get uniq name and name in experiment
					biological_replicat_name = ".".join( [ control_name, str( dict_paired_control[ 'PAIRED'][ url_controlled_by][ '1'][ 'biological_replicat'])])
					true_biological_replicat_name = ".".join( [ file_id + "." + control_true_name, str( dict_paired_control[ 'PAIRED'][ url_controlled_by][ '1'][ 'biological_replicat'])])



					list_info = [
									# "-".join( [ true_biological_replicat_name, str( technical_replicate_number)]),
									true_biological_replicat_name,
									dict_paired_control[ 'PAIRED'][ url_controlled_by][ '1'][ 'url_control_encode'] + ";" + dict_paired_control[ 'PAIRED'][ url_controlled_by][ '2'][ 'url_control_encode'],
									'PAIRED',
									"1",
									dict_paired_control[ 'PAIRED'][ url_controlled_by][ '1'][ 'url_download'] + ";" + dict_paired_control[ 'PAIRED'][ url_controlled_by][ '2'][ 'url_download'],
									dict_paired_control[ 'PAIRED'][ url_controlled_by][ '1'][ 'md5'] + ";" + dict_paired_control[ 'PAIRED'][ url_controlled_by][ '2'][ 'md5']
								]

					list_line.append( list_info)
					# technical_replicate_number += 1

					# Check if key already used in this batch
					if biological_replicat_name not in dict_used_control:
						dict_used_control[ biological_replicat_name] = experiment_name
					else:
						used_control = True
				except KeyError:
					error_experiment = True
					sys.stderr.write( url_controlled_by + " paired file is missing a buddy\n")


				# pprint.pprint( dict_paired_control[ 'PAIRED'][ url_controlled_by])





		# for chip biological replicat
		for biological_replicat in final_dict_experiment[ experiment_name][ 'chip']:
			# dealing with single end
			if  len( final_dict_experiment[ experiment_name][ 'chip'][ biological_replicat][ 'single-end']):
				for i in range( len( final_dict_experiment[ experiment_name][ 'chip'][ biological_replicat][ 'single-end'])):
					file_id = final_dict_experiment[ experiment_name][ 'chip'][ biological_replicat][ 'single-end'][ i][ 0].rsplit( "/", 1)[ 1]
					new_name = ".".join( [ file_id, biological_replicat.split( ".", 1)[ 1]])

					# list_final_info = [ new_name + "-" + str( i + 1)] + \
					list_final_info = [ new_name ] + \
					 					final_dict_experiment[ experiment_name][ 'chip'][ biological_replicat][ 'single-end'][ i]


					list_line.append( list_final_info)
					# file_experiment.write( "\t".join( list_final_info) + "\n")
					# print( "\t".join( list_final_info) + "\n")

			# dealing with paired end
			if  len( final_dict_experiment[ experiment_name][ 'chip'][ biological_replicat][ 'paired-end']):
				nb_technical_replicat = 0
				# for all technical replicat
				for technical_replicat in final_dict_experiment[ experiment_name][ 'chip'][ biological_replicat][ 'paired-end']:
					nb_technical_replicat += 1
					if len( final_dict_experiment[ experiment_name][ 'chip'][ biological_replicat][ 'paired-end'][ technical_replicat]) == 2:

						paired_info = join_paired_info( final_dict_experiment[ experiment_name][ 'chip'][ biological_replicat][ 'paired-end'][ technical_replicat], 0)
						paired_url_download =  join_paired_info( final_dict_experiment[ experiment_name][ 'chip'][ biological_replicat][ 'paired-end'][ technical_replicat], 3)
						paired_md5 = join_paired_info( final_dict_experiment[ experiment_name][ 'chip'][ biological_replicat][ 'paired-end'][ technical_replicat], 4)

						file_ids = "-".join( [ final_dict_experiment[ experiment_name][ 'chip'][ biological_replicat][ 'paired-end'][ technical_replicat][ 0][ 0].rsplit( "/", 1)[ 1],
											   final_dict_experiment[ experiment_name][ 'chip'][ biological_replicat][ 'paired-end'][ technical_replicat][ 1][ 0].rsplit( "/", 1)[ 1]
											 ])



						new_name = ".".join( [ file_ids, biological_replicat.split( ".", 1)[ 1]])

						# list_final_info =  	[ new_name + "-" + str( nb_technical_replicat),
						list_final_info =  	[ new_name,
											paired_info] + \
											final_dict_experiment[ experiment_name][ 'chip'][ biological_replicat][ 'paired-end'][ technical_replicat][ 0][ 1:3] + \
											[ paired_url_download,
											paired_md5]
						list_line.append( list_final_info)
						# file_experiment.write( "\t".join( list_final_info) + "\n")
						# print( "\t".join( list_final_info) + "\n")

					elif len( final_dict_experiment[ experiment_name][ 'chip'][ biological_replicat][ 'paired-end'][ technical_replicat]) <2:
						error_experiment = True
						sys.stderr.write( biological_replicat + "-" + str( nb_technical_replicat) + " paired end replicat has too few files\n")

						# list_final_info =  	[ biological_replicat + "-" + str( nb_technical_replicat)]
						list_final_info =  	[ biological_replicat]
						list_line.append( list_final_info)
					else:
						sys.stderr.write( biological_replicat + "-" + str( nb_technical_replicat) + " paired end replicat has too many files\n")
						error_experiment = True
						# list_final_info =  	[ biological_replicat + "-" + str( nb_technical_replicat),
						list_final_info =  	[ biological_replicat,
											paired_info]
						list_line.append( list_final_info)




		# If problem with info writting in error dir
		if error_experiment:
			dir_path_file_experiment = os.path.join( path_outdir_error)
		else:
			dir_path_file_experiment = os.path.join( path_outdir)


		# creating output dir
		if not os.path.exists( dir_path_file_experiment):
			os.makedirs( dir_path_file_experiment)

		# writting in file
		path_file_experiment = os.path.join( dir_path_file_experiment, experiment_name + "_summary.tab")
		file_experiment = open( path_file_experiment, 'w')

		for line in list_line:
			file_experiment.write( "\t".join( line) + "\n")

		file_experiment.close()

	# pprint.pprint( dict_used_control)
