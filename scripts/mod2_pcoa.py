#!/usr/bin/env python

import click
import matplotlib.pyplot as plt
from skbio.stats.ordination import OrdinationResults
import numpy as np
import pandas as pd
from collections import defaultdict

ALPHA = 1.0
LINE_WIDTH = 0.3

@click.group()
def mod1_pcoa():
	pass

@mod1_pcoa.command()
@click.option('--coords', required=True, type=click.Path(resolve_path=True, readable=True, exists=True), help='Coordinates file')
@click.option('--mapping_file', required=True, type=click.Path(resolve_path=True, readable=True, exists=True), help='Mapping file')
def body_site(coords,mapping_file):
	"""Generates as many figures as samples in the coordinates file (Figure 1)"""
	o = OrdinationResults.from_file(coords)
	x, y = o.site[:,0], o.site[:,1]

	# coordinates
	c_df = pd.DataFrame(o.site, o.site_ids)

	# mapping file
	mf = pd.read_csv(mapping_file, '\t', converters=defaultdict(str), index_col='#SampleID')
	mf = mf.loc[o.site_ids]

	# coloring:  
	# 'HMP-FECAL', 'GG-FECAL', 'PGP-FECAL' -- light red
	# 'AGP-FECAL' -- dark red
	# 'HMP-ORAL', 'PGP-ORAL' -- light blue
	# 'AGP-ORAL' -- dark blue
	# 'HMP-SKIN', 'PGP-SKIN' -- light orange
	# 'AGP-SKIN' -- dark orange

	cat_colors = {'HMP-FECAL': '#CB3F34', 'GG-FECAL': '#CB3F34', 'PGP-FECAL': '#CB3F34', 
				  'AGP-FECAL': '#74151C', 
				  'HMP-ORAL': '#64B4C5', 'PGP-ORAL': '#64B4C5',
				  'AGP-ORAL': '#1D3773', 
				  'HMP-SKIN': '#EDEE8E', 'PGP-SKIN': '#EDEE8E',
				  'AGP-SKIN': '#D96B2D'}

	for sample in mf.index:
		for cat, color in cat_colors.iteritems():
			sub_coords = c_df[mf.TITLE_BODY_SITE == cat]

			plt.scatter(sub_coords[0], sub_coords[1], color=color, edgecolor='k', lw=LINE_WIDTH, alpha=ALPHA)

		plt.scatter(c_df.loc[sample][0], c_df.loc[sample][1],
					color=cat_colors[mf.loc[sample]['TITLE_BODY_SITE']],
					s=270, edgecolor='w')
		plt.scatter(c_df.loc[sample][0], c_df.loc[sample][1],
					color=cat_colors[mf.loc[sample]['TITLE_BODY_SITE']],
					s=250, edgecolor='k')
		plt.axis('off')
		my_dpi = 72
		plt.savefig(sample+'.pdf', figsize=(1000/my_dpi, 1000/my_dpi), dpi=my_dpi)
		plt.close()

@mod1_pcoa.command()
@click.option('--coords', required=True, type=click.Path(resolve_path=True, readable=True, exists=True), help='Coordinates file')
@click.option('--mapping_file', required=True, type=click.Path(resolve_path=True, readable=True, exists=True), help='Mapping file')
def country(coords,mapping_file):
	"""Generates as many figures as samples in the coordinates file (Figure 1)"""
	o = OrdinationResults.from_file(coords)
	x, y = o.site[:,0], o.site[:,1]

	# coordinates
	c_df = pd.DataFrame(o.site, o.site_ids)

	# mapping file
	mf = pd.read_csv(mapping_file, '\t', converters=defaultdict(str), index_col='#SampleID')
	mf = mf.loc[o.site_ids]

	# coloring:  
	# Australia
	# Belgium
	# Canada
	# China
	# Finland
	# France
	# Germany
	# Great Britain
	# Ireland
	# Japan
	# Malawi - BLUE
	# Netherlands
	# New Zealand
	# Norway
	# Scotland
	# Spain
	# Switzerland
	# Thailand
	# United Arab Emirates
	# United Kingdom - PINK
	# United States of America
	# Venezuela - YELLOW
	# no_data

	cat_colors = {'Australia': '#FF0000',
				  'Belgium': '#FF0000',
				  'Canada': '#FF0000',
				  'China': '#FF0000',
				  'Finland': '#FF0000',
				  'France': '#FF0000',
				  'Germany': '#FF0000',
				  'Great Britain': '#FF0000',
				  'Ireland': '#FF0000',
				  'Japan': '#FF0000',
				  'Malawi': '#0000FF',
				  'Netherlands': '#FF0000',
				  'New Zealand': '#FF0000',
				  'Norway': '#FF0000',
				  'Scotland': '#FF0000',
				  'Spain': '#FF0000',
				  'Switzerland': '#FF0000',
				  'Thailand': '#FF0000',
				  'United Arab Emirates': '#FF0000',
				  'United Kingdom': '#FF66FF',
				  'United States of America': '#FF0000',
				  'Venezuela': '#FFFF00',
				  'no_data': '#FF0000'}

	for sample in mf.index:
		for cat, color in cat_colors.iteritems():
			sub_coords = c_df[mf.COUNTRY == cat]

			plt.scatter(sub_coords[0], sub_coords[1], color=color, edgecolor='k', lw=LINE_WIDTH, alpha=ALPHA)

		plt.scatter(c_df.loc[sample][0], c_df.loc[sample][1],
					color=cat_colors[mf.loc[sample]['COUNTRY']],
					s=270, edgecolor='w')
		plt.scatter(c_df.loc[sample][0], c_df.loc[sample][1],
					color=cat_colors[mf.loc[sample]['COUNTRY']],
					s=250, edgecolor='k')
		plt.axis('off')
		my_dpi = 72
		plt.savefig(sample+'.pdf', figsize=(1000/my_dpi, 1000/my_dpi), dpi=my_dpi)
		plt.close()

@mod1_pcoa.command()
@click.option('--coords', required=True, type=click.Path(resolve_path=True, readable=True, exists=True), help='Coordinates file')
@click.option('--mapping_file', required=True, type=click.Path(resolve_path=True, readable=True, exists=True), help='Mapping file')
@click.option('--color', required=True, type=str, help='Metadata category to set color by')
def gradient(coords,mapping_file, color):
	"""Generates as many figures as samples in the coordinates file (Figures 2 & 3)"""
	o = OrdinationResults.from_file(coords)
	c_df = pd.DataFrame(o.site, o.site_ids)

	# mapping file
	mf = pd.read_csv(mapping_file, '\t', converters=defaultdict(str), index_col='#SampleID')
	mf = mf.loc[o.site_ids]
	mf[color] = mf[color].convert_objects(convert_numeric=True)

	numeric = mf[~pd.isnull(mf[color])]
	non_numeric = mf[pd.isnull(mf[color])]

	color_array = plt.cm.jet_r(numeric[color]/max(numeric[color]))

	for sample in mf.index:

		# NUMERIC
		ids = numeric.index
		x, y = c_df.loc[ids][0], c_df.loc[ids][1]
		plt.scatter(x, y, c=numeric[color], cmap=plt.get_cmap('jet_r'),
					alpha=ALPHA, lw=LINE_WIDTH)

		#plt.colorbar()

		# NON-NUMERIC
		ids = non_numeric.index
		x, y = c_df.loc[ids][0], c_df.loc[ids][1]
		plt.scatter(x, y, c='0.5', alpha=ALPHA, lw=LINE_WIDTH)

		# INDIVIDUAL BIG DOT
		try:
			color_index = numeric.index.tolist().index(sample)
		except ValueError:
			color_index = None

		if color_index is None:
			_color = '0.5'
		else:
			_color = color_array[color_index]

		plt.scatter(c_df.loc[sample][0], c_df.loc[sample][1],
			color=_color,
			s=270, edgecolor='w')
		plt.scatter(c_df.loc[sample][0], c_df.loc[sample][1],
			color=_color,
			s=250, edgecolor='k')

		plt.axis('off')
		my_dpi = 72
		plt.savefig(sample+'.pdf', figsize=(1000/my_dpi, 1000/my_dpi), dpi=my_dpi)
		plt.close()

@mod1_pcoa.command()
@click.option('--coords', required=True, type=click.Path(resolve_path=True, readable=True, exists=True), help='Coordinates file')
@click.option('--mapping_file', required=True, type=click.Path(resolve_path=True, readable=True, exists=True), help='Mapping file')
@click.option('--color', required=True, type=str, help='Metadata category to set color by')
def double_gradient(coords,mapping_file, color):
	"""Generates as many figures as samples in the coordinates file (Figures 2 & 3)"""
	o = OrdinationResults.from_file(coords)
	c_df = pd.DataFrame(o.site, o.site_ids)

	# mapping file
	mf = pd.read_csv(mapping_file, '\t', converters=defaultdict(str), index_col='#SampleID')
	mf = mf.loc[o.site_ids]
	mf[color] = mf[color].convert_objects(convert_numeric=True)

	numeric = mf[~pd.isnull(mf[color])]
	non_numeric = mf[pd.isnull(mf[color])]

	numeric_UK = numeric[numeric.COUNTRY == 'United Kingdom']
	numeric_Other = numeric[numeric.COUNTRY != 'United Kingdom']

	color_array_UK = plt.cm.spring(numeric_UK[color]/max(numeric_UK[color]))
	color_array_Other = plt.cm.winter(numeric_Other[color]/max(numeric_Other[color]))

	for sample in mf.index:

		# NUMERIC OTHER
		ids = numeric_Other.index
		x, y = c_df.loc[ids][0], c_df.loc[ids][1]
		plt.scatter(x, y, c=numeric_Other[color], cmap=plt.get_cmap('winter'),
					alpha=ALPHA, lw=LINE_WIDTH)

		# NUMERIC UK
		ids = numeric_UK.index
		x, y = c_df.loc[ids][0], c_df.loc[ids][1]
		plt.scatter(x, y, c=numeric_UK[color], cmap=plt.get_cmap('spring'),
					alpha=ALPHA, lw=LINE_WIDTH)

		#plt.colorbar()

		# NON-NUMERIC
		ids = non_numeric.index
		x, y = c_df.loc[ids][0], c_df.loc[ids][1]
		plt.scatter(x, y, c='0.5', alpha=ALPHA, lw=LINE_WIDTH)

		# INDIVIDUAL BIG DOT
		try:
			if mf.loc[sample].COUNTRY == 'United Kingdom':
				color_index = numeric_UK.index.tolist().index(sample)
			else:
				color_index = numeric_Other.index.tolist().index(sample)
		except ValueError:
			color_index = None

		if color_index is None:
			_color = '0.5'
		else:
			if mf.loc[sample].COUNTRY == 'United Kingdom':
				_color = color_array_UK[color_index]
			else:
				_color = color_array_Other[color_index]

		plt.scatter(c_df.loc[sample][0], c_df.loc[sample][1],
			color=_color,
			s=270, edgecolor='w')
		plt.scatter(c_df.loc[sample][0], c_df.loc[sample][1],
			color=_color,
			s=250, edgecolor='k')

		plt.axis('off')
		my_dpi = 72
		plt.savefig(sample+'.pdf', figsize=(1000/my_dpi, 1000/my_dpi), dpi=my_dpi)
		plt.close()

if __name__ == '__main__':
	mod1_pcoa()
