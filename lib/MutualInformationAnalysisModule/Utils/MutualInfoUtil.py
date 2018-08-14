import time
import os
import errno
import uuid
import math
import pandas as pd
import numpy as np
import collections
import natsort
import uuid
import shutil
import itertools
import json
from itertools import combinations
#from pandas import DataFrame
import matplotlib
from fba_tools.fba_toolsClient import fba_tools
import matplotlib.pyplot as plt
from collections import OrderedDict
from copy import deepcopy

from Workspace.WorkspaceClient import Workspace as Workspace
from DataFileUtil.DataFileUtilClient import DataFileUtil
from KBaseReport.KBaseReportClient import KBaseReport

def log(message, prefix_newline=False):
	"""Logging function, provides a hook to suppress or redirect log messages."""
	print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


class MutualInfoUtil:
	def __init__(self, config):
		self.ws_url = config["workspace-url"]
		self.callback_url = config['SDK_CALLBACK_URL']
		self.token = config['KB_AUTH_TOKEN']
		self.shock_url = config['shock-url']
		self.dfu = DataFileUtil(self.callback_url)
		self.ws = Workspace(self.ws_url, token=self.token)
		self.scratch = config['scratch']

	def _mkdir_p(self, path):
		"""
		_mkdir_p: make directory for given path
		"""
		if not path:
			return
		try:
			os.makedirs(path)
		except OSError as exc:
			if exc.errno == errno.EEXIST and os.path.isdir(path):
				pass
			else:
				raise

	def test_dfu(self):
		output_directory = self.scratch
		#output_directory = "/kb/module/test1/"
		#os.mkdir(output_directory)
		#self._mkdir_p(output_directory)

		test_file = os.path.join(output_directory, 'index.html')
		with open(test_file, 'w') as file:
			file.write("test!")
		print("OUTPUT DIR")
		print(output_directory)
		print(os.listdir(output_directory))
		print("file_to_shock")
		report_shock_id = self.dfu.file_to_shock(
			{
				'file_path': output_directory,
				'pack': 'targz'
			})
		print(report_shock_id)
		return

	def _validate_run_flux_mutual_information_analysis_params(self, params):
		"""
		_validate_run_flux_mutual_information_analysis_params:
				validates params passed to run_flux_mutual_information_analysis method
		"""

		log('start validating validate_run_flux_mutual_information_analysis params')

		# check for required parameters
		for p in ['fbamodel_id', 'compounds', 'media_id', 'workspace_name']:
			if p not in params:
				raise ValueError('"{}" parameter is required, but missing'.format(p))

	def _get_file_from_ws(self, ref):
		try:
			file_path = self.ws.get_objects2({'objects':[{'ref': ref}]})
			file_path = file_path['data'][0]
		except Exception as e:
			raise ValueError(
				'Unable to get object from workspace: (' +
				ref + ')' + str(e))
		return file_path

	def _make_media_files(self, ws_name, base, compounds):
		"""
		Build and store media objects for each combination of compound added to the base media.
		:param base: The base media file
		:param compounds: the set of compound to test
		:return: A list of media ids and a matrix with each media combination defined
		"""

		ref = ws_name + "/" + base
		if base.find("/") != -1:
			ref = base

		output = self._get_file_from_ws(ref)
		base_media = output['data']
		base = output['info'][1]
		myuuid = str(uuid.uuid4())
		media_ids = [base]
		new_media_list = []
		media_matrix = [[""]+compounds]
		media_matrix.append([[base]+[0]*len(compounds)])
		for n_comp in range(1, len(compounds)+1):
			for combo in combinations(compounds, n_comp):
				new_media_id = base + '_v%s' % len(media_matrix)
				media_ids.append(new_media_id)
				media_matrix.append([new_media_id]+[1 if comp in combo else 0 for comp in compounds])
				new_media = deepcopy(base_media)
				new_media['id'] = new_media_id
				new_media['name'] = new_media_id
				for new_comp in combo:
					new_media['mediacompounds'].append(
						{'compound_ref': '48/1/1/compounds/id/%s' % new_comp.split('_')[0],
						 'concentration': 1.0, 'maxFlux': 1000, 'minFlux': -1000})
				new_media_list.append(new_media)

		print("Made %s Media Files" % (len(media_ids)-1))
		info = self.ws.save_objects(
			{'workspace': ws_name,
			 "objects": [{
				 "hidden": 1,
				 "type": "KBaseBiochem.Media",
				 "data": media,
				 "name": myuuid + "-" + media['name']
			 } for media in new_media_list]
			 })
		#print(info)
		return media_ids, media_matrix, myuuid

	def _run_fba(self, workspace_name, media_id_list, fbamodel_id, myuuid, base_media):
		print('running fba')
		fba_tool_obj = fba_tools(self.callback_url, service_ver='dev')
		new_media_list = []
		for media in media_id_list:
			if media == base_media:
				new_media_list.append(workspace_name + "/" + media)
			else:
				new_media_list.append(workspace_name + "/" + myuuid + "-" + media)


		fba_tool_obj.run_flux_balance_analysis({
			"max_c_uptake" : 6,
			"workspace" : workspace_name,
			"fbamodel_id" : fbamodel_id,
			"fba_output_id" : fbamodel_id + ".mifba",
			"fbamodel_workspace" : workspace_name,
			"media_id_list" : new_media_list,
			"target_reaction" : "bio1",
			"minimize_flux" : 1
			})
		output = self.ws.get_objects2({
			'objects' : [{
				'ref' : workspace_name + "/" + fbamodel_id + '.mifba'
			}]
			})

		#json.dump(output, open(self.scratch+'/fba.json', 'w'))

		fba = output['data'][0]['data']
		biomass_data = "FBAs,Biomass\n"
		secretion_file = ","+','.join(media_id_list)+"\n"
		full_secretion_file = ","+','.join(media_id_list)+"\n"
		full_flux_file = ","+','.join(media_id_list)+"\n"
		flux_file = ","+','.join(media_id_list)+"\n"
		objectives = fba['other_objectives']
		for i in range(0, len(objectives)):
			biomass_data = biomass_data + media_id_list[i] + "," + str(objectives[i]) + "\n"

		flux_vars = fba['FBAReactionVariables']
		for var in flux_vars:
			id = var['modelreaction_ref'].split("/").pop()
			flux_file = flux_file + id
			full_flux_file = full_flux_file + id
			fluxes = var['other_values']
			for i in range(0, len(objectives)):
				if objectives[i] == 0:
					full_flux_file = full_flux_file + ",0"
					flux_file = flux_file + ",0"
				else:
					full_flux_file = full_flux_file + "," + str(fluxes[i])
					if abs(fluxes[i]) < 1e-7:
						flux_file = flux_file + ",0"
					else:
						flux_file = flux_file + ",1"
			flux_file = flux_file + "\n"
			full_flux_file = full_flux_file + "\n"

		secretion_vars = fba['FBACompoundVariables']
		for var in secretion_vars:
			id = var['modelcompound_ref'].split("/").pop()
			secretion_file = secretion_file + id
			full_secretion_file = full_secretion_file + id
			fluxes = var['other_values']
			for i in range(0, len(objectives)):
				if objectives[i] == 0:
					full_secretion_file = full_secretion_file + ",0"
					secretion_file = secretion_file + ",0"
				else:
					full_secretion_file = full_secretion_file + "," + str(fluxes[i])
					if abs(fluxes[i]) < 1e-7:
						secretion_file = secretion_file + ",0"
					elif fluxes[i] < 0:
						secretion_file = secretion_file + ",-1"
					else:
						secretion_file = secretion_file + ",1"
			secretion_file = secretion_file + "\n"
			full_secretion_file = full_secretion_file + "\n"

		output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
		self._mkdir_p(output_directory)
		biomass_path = os.path.join(output_directory, 'biomass.csv')
		secretion_path = os.path.join(output_directory, 'secretion.csv')
		flux_path = os.path.join(output_directory, 'flux.csv')
		full_secretion_path = os.path.join(output_directory, 'full_secretion.csv')
		full_flux_path = os.path.join(output_directory, 'full_flux.csv')

		with open(biomass_path, 'w') as biomass_f:
			biomass_f.write(biomass_data)

		with open(secretion_path, 'w') as secretion_f:
			secretion_f.write(secretion_file)

		with open(flux_path, 'w') as flux_f:
			flux_f.write(flux_file)

		with open(full_secretion_path, 'w') as full_secretion_f:
			full_secretion_f.write(full_secretion_file)

		with open(full_flux_path, 'w') as full_flux_f:
			full_flux_f.write(full_flux_file)

		return [biomass_path,secretion_path,flux_path,full_secretion_path,full_flux_path]

	def _generate_html_report(self, result_directory, mutual_info_dict):

		"""
		_generate_html_report: generate html summary report
		"""
		#scratch, uui, datafileutil, file_to_shock, shockId, extended report

		log('start generating html report')

		html_report = list()

		output_directory = os.path.join(self.scratch, str(uuid.uuid4()))

		self._mkdir_p(output_directory)
		result_file_path = os.path.join(output_directory, 'mutual_information_report.html')

		shutil.copy(os.path.join(result_directory, 'MI_plot.png'),
					os.path.join(output_directory, 'MI_plot.png'))

		overview_content = ''
		overview_content += '<table><tr><th>Mutual Information for various chemical compound combinations'
		overview_content += ' Object</th></td>'
		overview_content += '<tr><th>Input Chemical Compound Combination</th>'
		overview_content += '<th>Mutual Information (in Bits)</th>'
		overview_content += '</tr>'
		for k, v in mutual_info_dict.items():
			overview_content += '<tr><td>{}</td><td>{}</td></tr>'.format(k, v)
		overview_content += '</table>'

		with open(result_file_path, 'w') as result_file:
			with open(os.path.join(os.path.dirname(__file__), 'report_template.html'),
					  'r') as report_template_file:
				report_template = report_template_file.read()
				report_template = report_template.replace('<p>Overview_Content</p>',
														  overview_content)
				result_file.write(report_template)

		report_shock_id = self.dfu.file_to_shock({'file_path': output_directory,
												  'pack': 'zip'})['shock_id']

		html_report.append({'shock_id': report_shock_id,
							'name': os.path.basename(result_file_path),
							'label': os.path.basename(result_file_path),
							'description': 'HTML summary report for Mutual Information App'})

		return html_report

	def _generate_report(self, result_directory, mutual_info_dict, workspace_name):
		"""
		_generate_report: generate summary report
		"""
		print('-->I am here *************')
		print(mutual_info_dict)
		uuidStr = str(uuid.uuid4())

		output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
		# output_dir = os.path.join(result_directory, uuidStr)
		self._mkdir_p(output_directory)
		shutil.copy2(os.path.join(os.path.dirname(__file__), 'data', 'index.html'),
					 output_directory)

		# shutil.copy('/kb/module/data/index.html', result_directory + '/' + uuidStr + '/index.html')
		json.dump(mutual_info_dict, open(os.path.join(output_directory, 'pdata.json'), 'w'))
		#shutil.copy('pdata.json', result_directory + '/' + uuidStr + '/pdata.json')

		# DataFileUtils to shock
		print(output_directory)
		print(os.listdir(output_directory))
		report_shock_result = self.dfu.file_to_shock({'file_path': output_directory,
												  'pack': 'zip'})
		report_shock_id = report_shock_result['shock_id']
		print(report_shock_result)

		report_file = {'name': 'index.html',
					   'description': 'the report',
					   'shock_id': report_shock_id}

		#biomass_file = {'name': 'biomass_file.csv',
		#			   'description': 'biomass_file',
		#			   'path': paths[0]}

		#flux_file = {'name': 'flux_file.csv',
		#			   'description': 'flux_file',
		#			   'path': paths[2]}

		#full_flux_file = {'name': 'full_flux_file.csv',
		#			   'description': 'full_flux_file',
		#			   'path': paths[4]}

		#secretion_file = {'name': 'secretion_file.csv',
		#			   'description': 'secretion_file',
		#			   'path': paths[1]}

		#full_secretion_file = {'name': 'full_secretion_file.csv',
		#			   'description': 'full_secretion_file',
		#			   'path': paths[3]}

		log('creating report')
		#output_html_files = self._generate_html_report(result_directory,
		#											   mutual_info_dict)
		report_params = {'message': '',
						 'workspace_name': workspace_name,
						 'html_links': [report_file],
						 'file_links': [],
						 'direct_html_link_index': 0,
						 'html_window_height': 333,
						 'report_object_name': 'MutualInfomation_report_' + uuidStr}

		kbase_report_client = KBaseReport(self.callback_url)
		output = kbase_report_client.create_extended_report(report_params)

		report_output = {'report_name': output['name'], 'report_ref': output['ref']}

		return report_output

	def _generate_mutual_info(self, media_matrix, fba_file, mi_options):

		#print('this is fba_file')
		#print(fba_file)
		df1 = pd.read_csv(fba_file[0])
		df1.values

		#df1.as_matrix()
		#print('-->printing df1')# **** rm
		#print(df1.to_string()) # **** rm
		#print(type(df1))  # **** rm
		#print('-->printing media_matrix')
		# print(media_matrix)

		df3= pd.DataFrame(columns=media_matrix[0][1:])
		for i in range(1, len(media_matrix)):
			if i==1:
				df3.loc[media_matrix[i][0][0]] = media_matrix[i][0][1:]
			else:
				df3.loc[media_matrix[i][0]] = media_matrix[i][1:]

		#print('-->*************OK')
		#print(df3)

	    #----Input validation of Media/FBAs with Binary Matrix FBAs------
		# 1.0 Number of rows in Media.csv file =  (Number of columns -1)
		#   1.0. If they are different: Through an ERROR saying missed match number of FBAs in media and binary matrix.
		# 1.1 Check whether the elements in Media.csv file contains only binary values (i.e. 0 and 1)
		#   1.1. If the elements are different: Through an ERROR saying not approapriate input values
		# 1.2 Check whether the compounds in Media.csv file match with number of FBAs
		#   1.2. If the compounds are different from number of FBAs: Through an ERROR saying not appropriate input values

		media_matrix = df3
		s_df1 = df1.shape
		s_df2 = media_matrix.shape
		#print(media_matrix,type(media_matrix))

		Temp_df2 = np.array(media_matrix.values)
		#print('-->******')
		#print(Temp_df2)
		# Create matrix with only the elements remove first column and all the rows
		Temp_df2 = Temp_df2[0:,1:]

		Bin_val_check =  np.array_equal(Temp_df2, Temp_df2.astype(bool))
		#num_compounds = (s_df2[1])-1
		num_compounds = s_df2[1]

		if ((s_df1[1]-1) != s_df2[0]) or (Bin_val_check != True) or (int(math.log(s_df2[0],2)) != num_compounds):
			print ('invalid input values')

		#-----All possible combination of the chemical compounds----------------------
		# 2.0 Sperating m0 from rest of the lables

		Temp1_df2 = media_matrix
		#print('-->*************OK')
		#print(Temp1_df2)
		cols = Temp1_df2.columns
		for i in range(0,len(cols)):
			Temp1_df2.loc[Temp1_df2[cols[i]] == 1 , cols[i]] = cols[i]
		#print('-->*************OK')
		#print (Temp1_df2)

		# 2.1 Creating a disctionary for all FBAs except m0
		#print(len(Temp1_df2))
		#print('--->*********')
		#print(Temp1_df2)

		mydict = {}
		for x in range(0, len(Temp1_df2)):
			for i in range(0, s_df2[1]):
				currentvalue = Temp1_df2.iloc[x, i]
				currentid = Temp1_df2.index[x]
				mydict.setdefault(currentid, [])
				if currentvalue != 0:
					mydict[currentid].append(currentvalue)
				# Add the first key as m0
		media_0_name = Temp1_df2.index[0]
		mydict[media_0_name] = ["0"]
		# Sort the keys
		mydict = collections.OrderedDict(natsort.natsorted(mydict.items()))
		#print ('--> ********')
		compoundslist = Temp1_df2.columns.get_values()
		compoundslist.tolist()
		#print(compoundslist)
		#print('all possible combination')
		#print(len(compoundslist))

		# List of Compounds combination in the list
		my_combi_list = []
		for L in range(0, len(compoundslist) + 1):
			for subset in itertools.combinations(compoundslist, L):
				my_combi_list.append(list(subset))

		my_combi_list[0] = [0]
		# print(my_combi_list)
		'''
		for k, v in mydict.iteritems():
			#print('--> ********')
			print(k, v)
		'''

		# Created a dictionary where the keys:
		# list of compounds combination
		# values are corresponding FBAs list in df2
		result_dict = {}
		for element in my_combi_list[1:]:
			for k, v in mydict.iteritems():
				if set(v).issubset(set(map(lambda x:str(x), element))):
					key = ','.join(map(lambda x:str(x), element))
					if result_dict.get(key):
						media_list = result_dict[key]
						media_list.append(k)
						media_list = list(set(media_list))
						result_dict.update({key: media_list})
					else:
						result_dict.update({key: [media_0_name, k]})

		# Sort the keys
		result_dict['0'] = [media_0_name]
		result_dict = collections.OrderedDict(natsort.natsorted(result_dict.items()))
		# print(result_dict)
		#print('-->I am here **** OK')
		#print(result_dict)
		#print (df1)

		# Created a dictionary where the keys are:
		# list of compounds combination
		# values are compounds combination FBAs with df1 vaules
		All_Comp_Combi_dic = {}
		for column, value in result_dict.items():
			All_Comp_Combi_dic.update({column: df1.get(value)})

		# print('-->All_Comp_Combi_dic******')
		# print (All_Comp_Combi_dic)
		# print(result_dict)

		# To print an item from the All_Comp_Combi_dic
		df = (pd.DataFrame(All_Comp_Combi_dic.items()))
		#print('--> printing df')
		#print(df[0].to_string())
		#print(df[1][7])

		######### INTRACELLULAR FLUX MUTUAL INFORMATION CALCULATION #############
		if mi_options == "flux":
				print('Intracellular flux')
				MI_dict = {}
				for k in range(0, len(df[0])):
					drop_rows_df = df[1][k].drop_duplicates(keep="first")
					drop_columns_df = drop_rows_df.T.drop_duplicates(keep="first").T
					remove = []
					removed = {}
					count_values = []
					cols = df[1][k].columns
					for i in range(len(cols)-1):
						duplicated = []
						v = df[1][k][cols[i]].values
						for j in range(i+1,len(cols)):
							if np.array_equal(v,df[1][k][cols[j]].values):
								remove.append(cols[j])
								duplicated.append(cols[j])
						if duplicated and cols[i] not in remove:
							removed.update({cols[i]:duplicated})
						count = {}
						for key, value in removed.items():
							count.update({key: len(value)})

						#print v

						# print drop_columns_df
						count_values = count.values()
						# print count_values
						count_values = map(lambda x: x+1, count_values)
						# print count_values
						d = {x:count_values.count(x) for x in count_values}
					#print('-->count_values')
					#print(count_values)

					#-------Mutual Inforamtion (MI) calculation-------------
					FBAs = len(df[1][k].columns)
					pure_entropy = math.log(FBAs,2)
					#print (pure_entropy) (-->ok rm)

					# If No duplicates exist and list "value" is empty
					if not count_values:
						#print("List is empty")
						No_duplicate_FBAs = len(drop_columns_df.columns)
						conditional_entropy = -1 * (No_duplicate_FBAs*((1/No_duplicate_FBAs)*((1/1)*math.log(1.0/1.0,2))));
						Mutual_Info = pure_entropy - conditional_entropy
						#print('Mutaul Info:', Mutual_Info)

					if count_values:
					# If duplicates exist and list "value" is not empty
						conditional_entropy = 0
						for key in d:
							#print key, d[key]
							Temp = -1 * d[key] * (key/float(FBAs)) * key * (1.0/key) * math.log(1.0/key,2)
							conditional_entropy = Temp + conditional_entropy
						#print "%3f" %Temp
						Mutual_Info = pure_entropy - conditional_entropy

					MI_dict[df[0][k]] = Mutual_Info
					MI_dict['0'] = 0.0

				#Sorted MI_dict
				MI_dict = sorted(MI_dict.items(), key=lambda x: (-len(x[0]), x[0]))
				MI_dict = OrderedDict(MI_dict)
				#print(MI_dict)

				#print('-->rest')
				#print(compoundslist)
				#print(num_compounds)

				x_groups = [[] for x in range(num_compounds)]
				y_groups = [[] for x in range(num_compounds)]
				names = [[] for x in range(num_compounds)]
				Comp_Mapping = [[] for x in range(num_compounds)]

				for key, val in MI_dict.iteritems():
					del_count = key.count(',')
					x_groups[del_count].append(key)
					y_groups[del_count].append(val)

					# for x, y in zip(x_groups, y_groups):
					# data.append(go.Bar(x=x, y=y, name='test'))

				pdata = []
				for i in range(0, len(x_groups)):
					names[i] = str(i + 1) + ' Compound Combination'
					Comp_Mapping = str(i + 1) + '-' + compoundslist[i]

					record = {}
					record["x"] = []
					for e in x_groups[i]:
						record["x"].append("c" + e)
					record["y"] = y_groups[i]
					record["names"] = names[i]
					record["Comp_Mapping"] = Comp_Mapping
					pdata.append(record)

				#print (pdata)
				#json.dump(pdata, open(self.scratch+'/pdata.json', 'w'))
				return pdata
				#return MI_dict

######### INPUT COMPONENTS AND BIOMASS FLUX MUTUAL INFORMATION CALCULATION #############
		if mi_options == "biomass":
				# Load the file contain the information of FBAs(media) along with corresponding Biomass (growth)
				print('biomass flux')
				df2 = pd.read_csv(fba_file[1])
				df2.values
				#print(df)

				MI_dict_biomass = {}
				for r in range(0, len(df[0])):
					reaction_states = df[1][r].head(1000)

					def get_groups(flux_df):
						groups = collections.defaultdict(list)
						unique = flux_df.aggregate(lambda x: hash(str(x.values)))
						for k, v in unique[0:].iteritems():
							groups[v].append(k)
						return dict([(i, g) for i, g in enumerate(groups.values())])

					n_group = collections.defaultdict(int)
					groups = get_groups(reaction_states)

					for group in groups.values():
						n_group[len(group)] += 1

					groups_count = {}
					for key, values in groups.items():
						groups_count[key] = len(values)
						# print groups_count

					# Take first FBA label of every group
					group_id = {}
					for k, v in groups.items():
						group_id.update({k: groups.values()[k][0]})

					# Obtain the Biomass of each Group
					cols_df = group_id.values()
					cols_df2 = df2.columns
					#print (cols_df)

					# Dictionary of first FBA label of every group and its corresponding number of members
					groups_label_count = {}
					for k, v in groups_count.items():
						groups_label_count.update({cols_df[k]: v})
					#print('groups_label_count')
					#print(groups_label_count)

					def get_cond_count(re_group):
						media_cond = 0
						for media in re_group['FBAs']:
							media_cond += groups_label_count[media]
						return media_cond


					# Extract FBA Groups biomass inside df2
					Groups_Biomass = df2[df2['FBAs'].isin(cols_df)]
					#print('-->I am here')
					#print(Groups_Biomass)

					# Regroup based on the biomass values
					re_group = Groups_Biomass.groupby('Biomass')
					biomass_FBAs_groups = re_group.aggregate(get_cond_count)

					biomass_FBAs_label_groups = Groups_Biomass.groupby("Biomass", sort=True).sum()
					#print(biomass_FBAs_label_groups)

					#print (biomass_FBAs_label_groups)

					Summery = pd.merge(left=biomass_FBAs_label_groups, left_index=True, right=biomass_FBAs_groups,
									   right_index=True,
									   how='inner')
					Data_4_CondMI = Summery.groupby('FBAs_y').count()
					Data_4_CondMI = Data_4_CondMI.to_dict(orient='dict')
					for k, v in Data_4_CondMI.items():
						Data_4_CondMI = v

					Num_of_FBAs = Data_4_CondMI.keys()
					Count_Num_of_FBAs = Data_4_CondMI.values()

					# -------Mutual Inforamtion (MI) calculation Stage II (input compounds respect to BIOMASS-------------
					# Pure Entropy
					FBAs = len(df[1][r].columns)
					pure_entropy = math.log(FBAs, 2)

					conditional_entropy = 0.0
					for l in range(0, len(Count_Num_of_FBAs)):
						temp = -1 * Count_Num_of_FBAs[l] * (Num_of_FBAs[l] / float(FBAs)) * Num_of_FBAs[l] * (
							1.0 / float(Num_of_FBAs[l]) * (math.log(1.0 / float(Num_of_FBAs[l]), 2)))
						conditional_entropy += temp

					Mutual_Info_Biomass = pure_entropy - conditional_entropy
					# print('Mutaul Info:', Mutual_Info_Biomass)

					#print(Mutual_Info_Biomass)
					MI_dict_biomass.update({df[0][r]: Mutual_Info_Biomass})


					#print(MI_dict_biomass)


				# Sorted MI_dict_biomass
				MI_dict_biomass = sorted(MI_dict_biomass.items(), key=lambda x: (-len(x[0]), x[0]))
				MI_dict_biomass = OrderedDict(MI_dict_biomass)

				#print(MI_dict_biomass)

				x_groups = [[] for x in range(num_compounds)]
				y_groups = [[] for x in range(num_compounds)]
				names = [[] for x in range(num_compounds)]
				Comp_Mapping = [[] for x in range(num_compounds)]

				for key, val in MI_dict_biomass.iteritems():
					del_count = key.count(',')
					x_groups[del_count].append(key)
					y_groups[del_count].append(val)

				pdata = []
				for i in range(0, len(x_groups)):
					names[i] = str(i + 1) + ' Compound Combination'
					Comp_Mapping = str(i + 1) + '-' + compoundslist[i]

					record = {}
					record["x"] = []
					for e in x_groups[i]:
						record["x"].append("c" + e)
					record["y"] = y_groups[i]
					record["names"] = names[i]
					record["Comp_Mapping"] = Comp_Mapping
					pdata.append(record)
				return pdata

######### INPUT COMPONENTS AND BIOMASS, SECRETION FLUX MUTUAL INFORMATION CALCULATION #############

		if mi_options == "secretion":
			"""
            #Load the file contain the information of FBAs(media) along with corresponding Biomass (growth)
            #print('secretion flux')
            #df4 = pd.read_csv(fba_file[2])
            #df4.values
            #print(df4)
            #df_biomass = pd.read_csv(fba_file[1])
            #df_biomass.values
            #print(df_biomass)

            MI_dict_b_u_s = {}
            for r in range(0, len(df[0])):

                reaction_states = df[1][r]

                def get_groups(flux_df):
                    groups = collections.defaultdict(list)
                    unique = flux_df.aggregate(lambda x: hash(str(x.values)))
                    for k, v in unique[0:].iteritems():
                        groups[v].append(k)
                    return dict([(i, g) for i, g in enumerate(groups.values())])

                n_group = collections.defaultdict(int)
                groups = get_groups(reaction_states)
                for group in groups.values():
                    n_group[len(group)] += 1
                # print n_group
                # print groups

                groups_count = {}
                for key, values in groups.items():
                    groups_count[key] = len(values)
                # print groups_count

                # Take first FBA label of every group
                group_id = {}
                for k, v in groups.items():
                    group_id.update({k: groups.values()[k][0]})

                # Obtain the Biomass of each Group
                cols_df = group_id.values()
                cols_df4 = df4.columns
                # print cols_df

                # Dictionary of first FBA label of every group and its corresponding number of members
                groups_label_count = {}
                for k, v in groups_count.items():
                    groups_label_count.update({cols_df[k]: v})
                print(groups_label_count)

                # Extract FBA Groups biomass inside df4
                Groups_Biomass = df4[df4['FBAs'].isin(cols_df)].reset_index(drop=True)
                # print cols_df
                # print Groups_Biomass

                # Regroup based on the biomass values
                re_group = Groups_Biomass.groupby(
                    ['cpd00254_e0', 'cpd00067_e0', 'cpd00205_e0', 'cpd00034_e0', 'cpd00009_e0', 'cpd00013_e0',
                     'cpd10515_e0', 'cpd00099_e0', 'cpd00030_e0', 'cpd00012_e0', 'cpd00048_e0', 'cpd00149_e0',
                     'cpd00073_e0', 'cpd10516_e0', 'cpd00001_e0', 'cpd00011_e0', 'cpd00007_e0', 'cpd11416_c0',
                     'cpd00058_e0', 'cpd00084_e0', 'cpd00063_e0', 'cpd00027_e0', 'cpd00029_e0', 'cpd00028_e0',
                     'Biomass'])

                my_list = []
                for index, values in re_group:
                    my_list.append(values['FBAs'].values)

                B_U_S_dict = {}
                for media in my_list:
                    if len(media) > 1:
                        media_cond = 0
                        for i in (0, len(media) - 1):
                            media_cond += groups_label_count[media[i]]
                        B_U_S_dict.update({str(media)[1:-1]: media_cond})
                    # final_my_dict.update({tuple(media.tolist()):media_cond})
                    else:
                        B_U_S_dict.update(
                            {str(media)[1:-1]: groups_label_count[str(tuple(media.tolist()))[1:-1][:-1][1:-1]]})
                B_U_S_dict = {eval(k): v for k, v in B_U_S_dict.iteritems()}
                print(B_U_S_dict)
                Summery = pd.DataFrame(B_U_S_dict.items(), columns=['FBAs_x', 'FBAs_y'])

                Data_4_CondMI = Summery.groupby('FBAs_y').count()
                Data_4_CondMI = Data_4_CondMI.to_dict(orient='dict')

                print(Data_4_CondMI)

                for k, v in Data_4_CondMI.items():
                    Data_4_CondMI = v

                Num_of_FBAs = Data_4_CondMI.keys()
                Count_Num_of_FBAs = Data_4_CondMI.values()
                print(Num_of_FBAs)
                print(Count_Num_of_FBAs)

                # -------Mutual Inforamtion (MI) calculation Stage II (input compounds respect to Biomass, Uptake and Secretion-------------
                # Pure Entropy
                FBAs = len(df[1][r].columns)
                pure_entropy = math.log(FBAs, 2)

                conditional_entropy = 0.0
                for l in range(0, len(Count_Num_of_FBAs)):
                    temp = -1 * Count_Num_of_FBAs[l] * (Num_of_FBAs[l] / float(FBAs)) * Num_of_FBAs[l] * (
                        1.0 / float(Num_of_FBAs[l]) * (math.log(1.0 / float(Num_of_FBAs[l]), 2)))
                    conditional_entropy += temp

                Mutual_Info_B_U_S = pure_entropy - conditional_entropy
                # print('Mutaul Info:', Mutual_Info_B_U_S)

                MI_dict_b_u_s.update({df[0][r]: Mutual_Info_B_U_S})

            # Sorted MI_dict_biomass
            MI_dict_b_u_s = sorted(MI_dict_b_u_s.items(), key=lambda x: (-len(x[0]), x[0]))
            MI_dict_b_u_s = OrderedDict(MI_dict_b_u_s)

            print(MI_dict_b_u_s)

            x_groups = [[] for x in range(num_compounds)]
            y_groups = [[] for x in range(num_compounds)]
            names = [[] for x in range(num_compounds)]
            Comp_Mapping = [[] for x in range(num_compounds)]

            for key, val in MI_dict_b_u_s.iteritems():
                del_count = key.count(',')
                x_groups[del_count].append(key)
                y_groups[del_count].append(val)

            # for x, y in zip(x_groups, y_groups):
            # data.append(go.Bar(x=x, y=y, name='test'))

            compound_IDs = ['H2', 'Vitamin K', 'Hematin', 'Glucose', 'Acetate', 'Formate', 'B12']

            pdata = []
            for i in range(0, len(x_groups)):
                names[i] = str(i + 1) + ' Compound Combination'
                Comp_Mapping = str(i + 1) + '-' + compound_IDs[i]

                record = {}
                record["x"] = []
                for e in x_groups[i]:
                    record["x"].append("c" + e)
                record["y"] = y_groups[i]
                record["names"] = names[i]
                record["Comp_Mapping"] = Comp_Mapping
                pdata.append(record)

            print(pdata)
"""