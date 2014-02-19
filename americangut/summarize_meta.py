#!/usr/bin/env python 
# summarize_meta.py

from __future__ import division
from make_phyla_plots import map_to_2D_dict
from build_cat_table import check_category, identify_groups, sort_alphabetically
from numpy import zeros, array, ones


def get_meta_counts(mapping, category):
	"""Determines the number of elements in each group

	INPUTS:
		mapping -- a 2D dictionary of the mapping data for the study

		category -- a string describing the metadata field to be summarized
	OUTPUTS:
		counts -- a dictionary keying each group to the associated counts

	"""

	# Identifies the groups in the metadata
	groups = identify_groups(mapping, category)

	counts = {g:0 for g in list(groups)}
	for (sample, meta) in mapping.iteritems():
		counts[meta[category]] = counts[meta[category]] + 1

	return counts

def convert_counts(counts, ordering=None, na_val=None, ignore_nas=False, 
	normalize=False):
	"""Converts a group counts dictionary to a plotable vector
	INPUTS:
		counts -- a dictionary of the group and group frequency

		ordering -- a string or list describing the way data should be ordered.
					The string can have preset ordering values of 'ALPHA', 
					indicating group sorting from A to Z, 'ALPHA_REV' indicating 
					group sorting from Z to A, 'FREQ' indicating sorting of 
					group frequencies from least common to most common, or 
					'FREQ_REV', sorting from most frequent to least frequent. 
					Ordering can also be a list specifying the custom ordering 
					for the groups. A value of None will not order the data.
					DEFAULT: None
		
		na_vals -- a string describing the way NA values are specified.
					DEFAULT: None

		ignore_nas -- a binary value specifying if NA values should be omitted 
					during the sorting. If ignore_nas is true, than an na_val 
					must be specified. This value is ignored if a custom 
					ordering list is supplied containing the na_val.
					DEFAULT: False

		normalize -- a binary value idicating whether the data should be 
					normalized.
					DEFAULT: False

	OUTPUTS:
		ordered_groups -- a list of groups in the correct order

		ordered_counts -- a list of counts corresponding to the groups
	"""

	# Checks the sanity of the inputs
	if not isinstance(counts, dict):
		raise TypeError('counts must be a dict.')
	# Performs sanity checks on NA handing inputs
	if not isinstance(na_val, (str, int, float)) and not na_val == None:
		raise TypeError('na_val must be a string or number.')
	if not isinstance(ignore_nas, bool):
		raise TypeError('ignore_nas must be a boolian.')
	if na_val == None and ignore_nas == True:
		raise ValueError('An na_val must be supplied to ignore NAs.')

	# Preforms sanity check on the normalize argument
	if not isinstance(normalize, bool):
		raise TypeError('Normalize must be a boolian.')

	# Preforms sanity check ordering
	if not isinstance(ordering, (str, list)) and ordering is not None:
		raise TypeError('Ordering must be a key word string or list.')
	if isinstance(ordering, str) and ordering not in set(['ALPHA', 'ALPHA_REV', 
		'FREQ', 'FREQ_REV']):
		raise ValueError('%s is not a supported ordering mode.' % ordering)
	if isinstance(ordering, list) and \
		not set(ordering).issubset(set(counts.keys())):
		raise ValueError('ordering contains a group not in counts.')

	# Handles group ordering
	groups = list(counts.keys())		
	num_groups = len(groups)

	# Assigns each group its corresponding count
	freq = []
	# Orders the data using the custom list
	if isinstance(ordering, list):
		for item in ordering:
			freq.append(counts[item])
		groups = ordering
	else:			
		new_groups = []
		for group in groups:
			if ignore_nas and group == na_val:
				continue
			new_groups.append(group)
			freq.append(counts[group])
		groups = new_groups

	num_groups = len(groups)
	
	# Sorts frequency data
	if ordering == 'ALPHA':
		(s_groups, order) = sort_alphabetically(list(counts.keys()))
	elif ordering == 'ALPHA_REV':
		(s_groups, order) = sort_alphabetically(list(counts.keys()))
		order = order[::-1]
	elif ordering == 'FREQ':
		order = sorted(range(len(freq)), key=freq.__getitem__)
	elif ordering == 'FREQ_REV':
		order = sorted(range(len(freq)), key=freq.__getitem__)
		order = order[::-1]	
	else:
		order = [i for (i, g) in enumerate(groups)]

	# Outputs the final data
	ordered_groups = []
	ordered_counts = []

	for idx in order:
		ordered_groups.append(groups[idx])
		ordered_counts.append(freq[idx])

	# Normalizes if desired
	if normalize:
		ordered_counts = array([ordered_counts])/sum(ordered_counts)
	else:
		ordered_counts = array([ordered_counts])


	return ordered_groups, ordered_counts

def compare_counts(maps, category, ordering=None, na_val=None, 
	ignore_nas=False, normalize=False):
	"""Combines studies into a numpy array
	INPUTS:
		maps -- a list of mapping dictionaries which contain a common category

		category -- a category in the mapping dictionary

		ordering -- a string or list describing the way data should be ordered.
					The string can have preset ordering values of 'ALPHA', 
					indicating group sorting from A to Z, 'ALPHA_REV' indicating 
					group sorting from Z to A, 'FREQ' indicating sorting of 
					group frequencies from least common to most common, or 
					'FREQ_REV', sorting from most frequent to least frequent. 
					Ordering can also be a list specifying the custom ordering 
					for the groups. A value of None will not order the data.
					DEFAULT: None
		
		na_vals -- a string describing the way NA values are specified.
					DEFAULT: None

		ignore_nas -- a binary value specifying if NA values should be omitted 
					during the sorting. If ignore_nas is true, than an na_val 
					must be specified. This value is ignored if a custom 
					ordering list is supplied containing the na_val.
					DEFAULT: False

		normalize -- a binary value idicating whether the data should be 
					normalized.
					DEFAULT: False

	OUTPUTS:
		ordered_groups -- a list of groups in the correct order

		ordered_table -- a list of counts corresponding to the groups"""
	
	# Checks maps is a list
	if not isinstance(maps, list):
		raise TypeError('maps must be a list')
	# Performs sanity checks on NA handing inputs
	if not isinstance(na_val, (str, int, float)) and not na_val == None:
		raise TypeError('na_val must be a string or number.')
	if not isinstance(ignore_nas, bool):
		raise TypeError('ignore_nas must be a boolian.')
	if na_val == None and ignore_nas == True:
		raise ValueError('An na_val must be supplied to ignore NAs.')

	# Preforms sanity check on the normalize argument
	if not isinstance(normalize, bool):
		raise TypeError('Normalize must be a boolian.')

	# Preforms sanity check ordering
	if not isinstance(ordering, (str, list)) and ordering is not None:
		raise TypeError('Ordering must be a key word string or list.')
	if isinstance(ordering, str) and ordering not in set(['ALPHA', 'ALPHA_REV', 
		'FREQ', 'FREQ_REV']):
		raise ValueError('%s is not a supported ordering mode.' % ordering)

	# Gets counts for the category in each map
	counts = []
	groups = set()
	for mapping in maps:
		map_counts = get_meta_counts(mapping=mapping, category=category)
		for key in map_counts:
			groups.add(key)
		counts.append(map_counts)

	# Removes NA from the list of groups
	if ignore_nas:
		groups.remove(na_val)
	
	# Pre-alocates an output array
	num_sets = len(maps)
	num_groups = len(groups)
	freq_table = zeros((num_groups, num_sets))	
	groups = list(groups)

	# Assigns counts to each group
	new_groups = []
	for x_idx, group in enumerate(groups):
		if group == na_val and ignore_nas:
			print 'SKIP'
			continue
		new_groups.append(group)
		for y_idx, map_counts in enumerate(counts):
			if group in map_counts:
				freq_table[x_idx, y_idx] = map_counts[group]

	# Determines how the data should be sorted
	if ordering == 'ALPHA':
		(s_group, order) = sort_alphabetically(groups)
	elif ordering == 'ALPHA_REV':
		(s_group, order) = sort_alphabetically(groups)
		order = order[::-1]
	elif ordering == 'FREQ':
		total_counts = freq_table.sum(axis=1)
		order = total_counts.argsort()
	elif ordering == 'FREQ':
		total_counts = freq_table.sum(axis=1)
		order = total_counts.argsort(reverse=True)
	elif isinstance(ordering, list):
		order = []
		for item in ordering:
			order.append(groups.index(item))
	else:
		order = [i for (i, g) in enumerate(groups)]

	num_groups = len(order)

	# Outputs the final data
	ordered_table = zeros((num_groups, num_sets))
	ordered_groups = []

	if normalize:
		freq_table = freq_table/(ones((num_sets, 1))*freq_table.sum(axis=0))

	for idx, pos in enumerate(order):
		ordered_groups.append(groups[pos])
		ordered_table[idx,:] = freq_table[pos,:]

	return ordered_groups, ordered_table

