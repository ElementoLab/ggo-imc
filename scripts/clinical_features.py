import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import os
import matplotlib
import anndata
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["axes.grid"] = False
matplotlib.use('Agg')

import yaml
metadata_filename = 'metadata/ggo_config.yml'

with open(metadata_filename, "r") as stream:
	try:
		metadata = yaml.safe_load(stream)
	except yaml.YAMLError as exc:
		print(exc)

adata_dict = dict()
print(f"Reading {metadata['PANEL_G']['AnnData']['phenotyped_file_name']}...")
adata_dict['PANEL_G'] = sc.read(metadata['PANEL_G']['AnnData']['phenotyped_file_name'])
print(f"Reading {metadata['PANEL_H']['AnnData']['phenotyped_file_name']}...")
adata_dict['PANEL_H'] = sc.read(metadata['PANEL_H']['AnnData']['phenotyped_file_name'])

df = pd.read_excel(metadata['clinical_data'])
'''
for panel in adata_dict:
	adata = adata_dict[panel]
	adata.obs['PID'] = adata.obs['sample_y'].str[:-1]

	for feature in ['radio', 'pathology', 'pred', 'Smoking Status', 'Gender', 'Race']:
		print('Patient: ', adata.obs[['PID', feature]].drop_duplicates().groupby(feature).count())
		print('ROI: ', adata.obs[['roi', feature]].drop_duplicates().groupby(feature).count())

	df = adata.obs
	race = df[['PID', 'Race']].drop_duplicates()
	print(race.groupby('Race').count()['PID'])
'''
os.makedirs('figures/demographics', exist_ok = True)
for panel in adata_dict:
	adata = adata_dict[panel]
	adata.obs['PID'] = adata.obs['sample_y'].str[:-1]

	features = ['radio', 'pathology', 'pred', 'Smoking Status', 'Gender', 'Race', 'Age', 'clinical tumor SIZE (cm)', 'Packyrs']
	for group in ['PID', 'roi', 'sample_y']:
		feat = [group] + features

		df = adata.obs[feat].drop_duplicates()
		for feature in features:
			print('Patient: ', df.groupby(feature).count()[group])
			print('ROI: ', df.groupby(feature).count()[group])
		#df['pathology'] = pd.Categorical(df['pathology'])
		df['radio'] = pd.Categorical(df['radio'], categories = ['PNS', 'PS', 'S'])
		df['Smoking Status'] = pd.Categorical(df['Smoking Status'], categories = reversed(['Never', 'Former smokes', 'Current smokes']))

		fig, ax = plt.subplots(1,1,dpi = 300, figsize = (2.5,1.2))
		sns.histplot(df, x = 'Age', ax = ax, color = '#EEEEEE')
		sns.despine()
		plt.savefig(f'figures/demographics/age_{panel}_{group}_hist.pdf', bbox_inches = 'tight')
		plt.close()

		fig, ax = plt.subplots(1,1,dpi = 300, figsize = (2.5,1.2))
		sns.histplot(df, x = 'clinical tumor SIZE (cm)', ax = ax, color = '#EEEEEE', cumulative = True, stat = 'proportion')
		sns.despine()
		plt.xlabel('Tumor size (cm)')
		plt.xlim(0,4)
		plt.savefig(f'figures/demographics/tumor_size_{panel}_{group}_hist.pdf', bbox_inches = 'tight')
		plt.close()

		fig, ax = plt.subplots(1,1,dpi = 300, figsize = (2.5,1.2))
		sns.kdeplot(df, x = 'clinical tumor SIZE (cm)', ax = ax, color = '#222222', cumulative = True)
		sns.despine()
		plt.xlabel('Tumor size (cm)')
		plt.xlim(0,4)
		plt.savefig(f'figures/demographics/tumor_size_{panel}_{group}_kde.pdf', bbox_inches = 'tight')
		plt.close()

		fig, ax = plt.subplots(1,1,dpi = 300, figsize = (3.5,3))
		df.groupby('Smoking Status').count()[group].plot(kind='pie', autopct='%1.1f%%', fontsize=10, ax = ax, colors = metadata['smoking_status_color'])
		plt.ylabel('')
		plt.tight_layout()
		plt.savefig(f'figures/demographics/smoker_{panel}_{group}_ratio.pdf', bbox_inches = 'tight')
		plt.close()

		fig, ax = plt.subplots(1,1,dpi = 300, figsize = (3.5,3))
		df['Race'] = df['Race'].replace('African American', 'Black or African American')
		df['Race'] = pd.Categorical(df['Race'], categories = ['White', 'Asian', 'Black or African American', 'Other', 'Unknown',])
		df.groupby('Race').count()[group].plot(kind='pie', autopct='%1.1f%%', fontsize=10, ax = ax, colors = ['#FFF7D2', '#FDE782', '#934F0A', '#222222', '#A2A2A2'])
		plt.ylabel('')
		plt.tight_layout()
		plt.savefig(f'figures/demographics/race_{panel}_{group}_ratio.pdf', bbox_inches = 'tight')
		plt.close()

		fig, ax = plt.subplots(1,1,dpi = 300, figsize = (3,2))
		x = 'clinical tumor SIZE (cm)'
		y = 'radio'
		sns.violinplot(df[df['radio'] != 'N'], x = x, y = y, hue = y, dodge = False, ax = ax, palette = metadata['radio_color'][1:])
		sns.swarmplot(df[df['radio'] != 'N'], x = x, y = y, color = 'black', s = 4, ax = ax)
		ax.legend().remove()
		plt.ylabel('')
		plt.xlabel('Clinical Tumor Size (cm)')
		sns.despine()
		plt.savefig(f'figures/demographics/clinical_tumor_size (cm)_{panel}_{group}.pdf', bbox_inches = 'tight')
		plt.close()


		fig, ax = plt.subplots(1,1,dpi = 300, figsize = (3,2))
		x = 'Packyrs'
		y = 'Smoking Status'
		df['Smoking Status'] = pd.Categorical(df['Smoking Status'], categories = reversed(['Never', 'Former smokes', 'Current smokes']))
		df['Packyrs'] = df['Packyrs'].fillna(0)
		ax = sns.violinplot(df, x = x, y = y, hue = y, dodge = False, ax = ax, palette=metadata['smoking_status_color'])
		sns.swarmplot(df, x = x, y = y, color = 'black', s = 4, ax = ax)
		plt.ylabel('')
		plt.xlabel('Packyrs')
		ax.legend().remove()
		sns.despine()
		plt.savefig(f'figures/demographics/packyrs_{panel}_{group}.pdf', bbox_inches = 'tight')
		plt.close()

		fig, ax = plt.subplots(1,1, dpi = 300, figsize = (3,1.5))
		df.groupby(['radio','pathology']).count()[group].reset_index().pivot(index = 'radio', columns = 'pathology').plot(kind='barh', stacked=True, width=0.8, ax = ax, edgecolor = 'black', color = metadata['pathology_color'])
		ax.legend().remove()
		ax.set_ylabel('Solidity')
		ax.set_xlabel('Frequency')
		sns.despine()
		plt.savefig(f'figures/demographics/Solidity_vs_Invasiveness_count_{panel}_{group}.pdf', bbox_inches = 'tight')
		plt.close()

# filtered comparison
for panel in adata_dict:
	adata = adata_dict[panel]
	filt = [{'pathology': 'Normal', 'radio': 'N'}, {'pathology': 'AIS', 'radio': 'PNS'}, {'pathology': 'MIA', 'radio': 'PS'}, {'pathology': 'IAC', 'radio': 'S'}]
	
	adata = anndata.concat([
		adata[(adata.obs['pathology'] == s['pathology']) & (adata.obs['radio'] == s['radio'])] for s in filt
	])
	adata.obs['PID'] = adata.obs['sample_y'].str[:-1]
	feat = ['PID'] + features
	df = adata.obs[feat].drop_duplicates()
	df = df[df['radio'] != 'N']



	sns.scatterplot(
		df,
		x = 'Packyrs',
		y = 'clinical tumor SIZE (cm)',
		hue = 'radio'
	)

	plt.ylabel('Clinical Tumor Size (cm)')
	plt.xlabel('Pack Years')
	sns.despine()
	plt.savefig(f'figures/demographics/Packyrs vs tumor size_{panel}_{group}.pdf', bbox_inches = 'tight')
	plt.close()

	fig, ax = plt.subplots(1,1,dpi = 300, figsize = (3,2))
	x = 'clinical tumor SIZE (cm)'
	y = 'radio'
	sns.violinplot(df[df['radio'] != 'N'], x = x, y = y, hue = y, dodge = False, ax = ax, palette = metadata['radio_color'][1:])
	sns.swarmplot(df[df['radio'] != 'N'], x = x, y = y, color = 'black', s = 4, ax = ax)
	ax.legend().remove()
	plt.ylabel('')
	plt.xlabel('Clinical Tumor Size (cm)')
	sns.despine()
	plt.savefig(f'figures/demographics/clinical_tumor_size (cm)_{panel}_{group}_filtered.pdf', bbox_inches = 'tight')
	plt.close()
		# filt = [{'pathology': 'Normal', 'radio': 'N'}, {'pathology': 'AIS', 'radio': 'PNS'}, {'pathology': 'MIA', 'radio': 'PS'}, {'pathology': 'IAC', 'radio': 'S'}]
		
		# density = anndata.concat([
		#     density[(density.obs['pathology'] == s['pathology']) & (density.obs['radio'] == s['radio'])] for s in filt
		# ])