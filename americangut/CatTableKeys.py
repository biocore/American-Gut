from numpy import mean, median, std, ndarray, sum as nsum, array, ones
from inspect import isfunction
from TaxPlot import check_data_array
from build_cat_table import (identify_groups, 
							 build_sub_table, 
							 identify_sample_group)

KEY_PROPERTIES = set(['data_type', 'data_mode', 'category', 'match_id', 
					  'name_type', 'name_val', 'name_disp', 
					  'name_delim', 'name_delim_pos',
					  'first_cat', 'reverse', 'data', 'samples', 'taxa', 
					  'meta'])

POSSIBLE_TYPES = set(['ID', 'POP', 'GROUP'])    

POSSIBLE_MODES = set(['ALL', 'AVR', 'MED', 'STDEV', 'BOOL', 'SUM', 'CUSTOM'])   

REQUIRED_KEYS = set(['TYPE', 'MODE', 'CATEGORY', 'MATCH_ID', 'NAME_FUN',
                      'NAME_DISPLAY', 'ORDER', 'REVERSE'])

POSSIBLE_NAMES = set(['RAW', 'SPLIT', 'SUB', 'CLEAN', 'S&C'])

POSSIBLE_DISPLAYS = set(['ID', 'CAT', 'DESCR'])

def return_original(data):
	"""Returns the original data set"""
	return data

def return_binary(data):
	"""Returns a binary version of the data set"""
	if not isinstance(data, (int, float, ndarray)):
		raise TypeError('Data must be numeric')
	return data > 0

def return_binary_sum(data):
	"""Returns the number of instances of the data set"""
	if not isinstance(data, (int, float, ndarray)):
		raise TypeError('Data must be numeric')
	return nsum(data>0, 1)

def return_mean(data):
	"""Returns the mean of each row in the data set"""
	if not isinstance(data, (int, float, ndarray)):
		raise TypeError('Data must be numeric')
	return mean(data, 1)

def return_median(data):
	"""Returns the mean of each row in the data set"""
	if not isinstance(data, (int, float, ndarray)):
		raise TypeError('Data must be numeric')
	return median(data, 1)

def return_stdev(data):
	"""Returns the standard deviation of each row in the data set"""
	if not isinstance(data, (int, float, ndarray)):
		raise TypeError('Data must be numeric')
	return std(data, 1)

def return_sum(data):
	"""Returns the sum of each row in the data set"""
	if not isinstance(data, (int, float, ndarray)):
		raise TypeError('Data must be numeric')
	return nsum(data, 1)


FUNCTION_LOOKUP = {'ALL': (return_original, return_original),
                   'BIN': (return_binary, return_original),
                   'BIS': (return_binary_sum, 'Counts'),
                   'AVR': (return_mean, 'Mean'),
                   'MED': (return_median, 'Median'),
                   'STD': (return_stdev, 'Standard Deviation'),
                   'SUM': (return_sum, 'Sum')}

def check_meta(meta, samples, category):
	"""Checks the meta data is sane"""
	# Checks the metadata object is a dictionary
	if not isinstance(meta, dict) or not isinstance(meta[meta.keys()[0]], dict):
		raise TypeError('meta must be a 2-D dictionary.')
	# Checks that every sample id is respresented in the metadata
	if samples is not None:
		for id_ in samples:
			if id_ not in meta.keys():
				raise ValueError('Not all samples are represented in the meta'
					' data')

	# Checks the category is in the mapping data
	if category is not None:
		groups = identify_groups(meta, category)

class CatTable:	
	"""DOC STRING

	<MORE INFORMATION>"""

	# Sets up intial values
	data_type = 'POP'
	data_mode = 'ALL'
	category = None
	group = None
	match_id = None
	name_type = 'RAW'
	name_delim = None
	name_delim_pos = None
	name_val = None
	name_disp = 'ID'
	first_cat = None
	reverse = False
	data = None
	samples = None
	taxa = None
	meta = None
	table_out = None
	samples_out = None
	names_out = None
	# Sets up function names

	def __init__(self, **kwargs):
		"""Initializes an instance of the class"""
		# Sets any positional arguments
		self = self.add_attributes(**kwargs)
		# Checks the object is sane
		self = self.check_cat_table()

	def add_attributes(self, **kwargs):
		"""Adds keyword arguments"""
		for (k, v) in kwargs.iteritems():
			if k in KEY_PROPERTIES:
				setattr(self, k, v)
			else:
				raise ValueError('%s is not a CatTable property.' %k)
		return self

	def define_data_object(self):
		"""Creates an output data object"""

		# Gets a list of sample ids to plot
		if self.data_type == 'POP':
			self.samples_out = self.samples
			self.group = 'Population'

		if self.data_type == 'ID':
			self.samples_out = [self.match_id]
			self.group = self.match_id
		
		if self.data_type == 'GROUP':
			# Determines where samples fall in each group
			group_assign = identify_sample_group(sample_ids=self.samples,
												 mapping=self.meta,
												 category=self.category)
			# Determines if we are matching a group or a sample
			if self.match_id in group_assign:
				self.group = self.match_id
			else:
				self.group = self.meta[self.match_id][self.category]

			self.samples_out = group_assign[self.group]

		# Creates a subtable of those sample ids
		sub_table = build_sub_table(data=self.data, 
									sample_ids=self.samples, 
									target_ids=self.samples_out)

		# Summarizes the data using the data mode
		function = FUNCTION_LOOKUP[self.data_mode][0]
		self.table_out = function(sub_table)
		if len(self.table_out.shape) == 1:
			# print 'Flag'
			self.table_out = (array([[1]])*self.table_out).transpose()

		return self

	def define_name_object(self):
		"""Sets up the output names in ID mode"""
		if self.samples_out is None:
			raise ValueError('Samples out cannot be None.')

		# Defines the name function
		naming = {'RAW':   lambda x: x,
				  'SPLIT': lambda x: x.split(self.name_delim)
				  		   	   [self.name_delim_pos],
				  'SUB':   lambda x: self.name_val,
				  'CLEAN': lambda x: x.replace('_', ' ').capitalize(),
				  'S&C':   lambda x: x.split(self.name_delim)
				  		   	   [self.name_delim_pos].capitalize()}
		name_function = naming[self.name_type]

		# For a single ID, the cleaned ID is returned.
		if self.data_type == 'ID':
			self.names_out = [name_function(self.samples_out[0])]
		# Cleans up the IDS if the name mode is ID and multiple samples 
		# have been supplied.
		elif (self.data_mode == 'ALL' or self.data_mode == 'BIN') and \
			(self.name_disp == 'ID' or self.name_disp == 'DESCR'):
			self.names_out = []
			for id_ in self.samples_out:
				self.names_out.append(name_function(id_))
		# Cleans up the group name if the name mode is ID and a single group is 
		# output. For POP data, the group is designated as 'Population'.
		elif self.name_disp == 'ID' and not (self.data_mode == 'ALL' or \
			self.data_mode == 'BIN'):
			self.names_out = [name_function(self.group)]
		# Handles CATEGORICAL data. For population data, the category is 
		# 'Population'
		elif self.name_disp == 'CAT' and self.data_mode == 'ALL' or \
			self.data_mode == 'BIN':
			self.names_out = [name_function(self.category)]*len(self.names_out)
		# Handles categorical data for a single output
		elif self.name_disp == 'CAT' and self.data_type == 'GROUP':
			self.names_out = [name_function(self.category)]
		elif self.name_disp == 'CAT':
			self.names_out = [name_function('Population')]
		# Handles description data
		elif self.name_disp ==  'DESCR':
			self.names_out = [name_function(FUNCTION_LOOKUP[self.data_mode][1])]

		return self

	def check_cat_table(self):
		"""Checks the validity of a cat table object"""

		# Checks the data type
		if self.data_type not in POSSIBLE_TYPES:
			raise ValueError('The data type is not supported.')

		# Checks the data mode
		if self.data_mode not in POSSIBLE_MODES:
			raise ValueError('The data mode is not supported.')

		# Checks the category type is supported
		if not isinstance(self.category, str) and \
			self.category is not None:
			raise TypeError('The category must be a string.')

		# Checks the id is supported
		if not isinstance(self.match_id, str) and \
			self.match_id is not None:
			raise TypeError('The sample to match must be a string.')

		# Checks the name input is sane
		if self.name_type not in POSSIBLE_NAMES:
			raise ValueError('The name type is not supported.')

		# Checks the name display is sane
		if self.name_disp not in POSSIBLE_DISPLAYS:
			raise ValueError('The name display is not supported')
	
		# Checks the first category argument is a supported type
		if not isinstance(self.first_cat, (str, int)) and \
			self.first_cat is not None:
			raise TypeError('The first category must be an integer or a string.')

		# Checks if the sort direction is sane
		if not isinstance(self.reverse, bool):
			raise TypeError('Reverse sorting must be a binary object.')

		# Checks the data is sane
		if not (self.data is None and self.samples is None and \
				self.taxa is None) and (self.data is None or \
				self.samples is None or self.taxa is None):
				raise ValueError('Values must be supplied for the data matrix, '
					'samples and columns.')
		if self.data is not None and self.samples is not None and \
			self.taxa is not None:
			self.data = check_data_array(data = self.data, 
											  row_names = self.taxa, 
											  col_names = self.samples,
											  data_id = 'data', 
											  row_id='taxa',
											  col_id = 'samples')

		if self.meta is not None:
			check_meta(self.meta, self.samples, self.category)


		# Checks that arguments are sane based on the data type.
		if self.data_type == 'ID':
		 	if self.match_id is None:
				raise ValueError('A sample ID must be specified for ID data.')
			if self.samples is not None and \
				self.match_id not in set(self.samples):
				raise ValueError('The sample is not in the dataset.')

		if self.data_type == 'GROUP':
			if self.match_id is None:
				raise ValueError('A sample ID must be specified for GROUP data.')
			if self.category is None:
				raise ValueError('A category must be specified for GROUP data.')
			if self.meta is not None and self.samples is not None:
				possible_groups = identify_groups(mapping = self.meta, 
												  category = self.category)
				# name_set = set(self.samples).union(set())
				if self.match_id not in set(self.samples).union(
					set(possible_groups)):
					raise ValueError('The sample or group to match cannot be '
						'found.')

		# Updates the data object
		if self.data is not None and (self.data_mode is not 'GROUP' or \
			(self.data_mode is 'GROUP' and not self.self==None)):
			self = self.define_data_object()
			self = self.define_name_object()

		return self

	def get_data_type(self):
		"""Gives the type of data to be summarized in the table"""
		return self.data_type

	def get_data_mode(self):
		"""Gives the way the data should be summarized"""		
		return self.data_mode

	def get_category(self):
		"""Gives the metadata category"""
		return self.category

	def get_id(self):
		"""Returns the target id"""
		return self.match_id

	def get_name_format(self):
		"""Returns a description of how names are cleaned"""
		return self.name_type

	def get_name_display(self):
		"""Gives the name display mode"""
		return self.name_disp

	def get_first_category(self):
		"""Returns the first metadata category to be displayed"""
		return self.first_cat

	def get_sort_reverse(self):
		"""Returns the direction of sorting (True = High to Low)"""
		return self.reverse

	def get_data(self):
		"""Returns data object"""
		self = self.check_cat_table()
		return self.table_out, self.names_out, self.taxa

	def get_table_out(self):
		"""Returns the table for plotting"""
		if self.data is None:
			print('No data has been supplied.')
		self = self.check_cat_table()
		return self.table_out, self.names_out, self.taxa

	def set_first_cat(self, category):
		"""Sets the first meta data category to be displayed"""
		self.first_cat = category

	def set_data_type(self, data_type, **kwargs):
		"""Sets the data type"""
		self.data_type = data_type
		self = self.add_attributes(kwargs)
		# Checks the object is sane
		self = self.check_cat_table()
		# Returns the object
		return self

	def set_data_mode(self, data_mode, **kwargs):
		"""Sets the data type"""
		self.data_mode = data_mode
		self = self.add_attributes(kwargs)
		# Checks the object is sane
		self = self.check_cat_table()
		# Returns the object
		return self

	def set_match_id(self, match_id):
		"""Sets the reference"""
		self.match_id = match_id
		# Checks the object is sane
		self = self.check_cat_table()
		# Returns the object
		return self

	def set_category(self, category):
		"""Sets the meta data category"""
		self.category = category
		# Checks the object is sane
		self = self.check_cat_table()
		# Returns the object
		return self

	def set_name_mode(self, name_mode, **kwargs):
		"""Sets up the way the name will be rendered"""
		self.name_type = name_mode
		self = self.add_attributes(kwargs)
		# Checks the object is sane
		self = self.check_cat_table()
		# Returns the object
		return self

	def set_custom_name_mode(self, custom_function):
		"""Allows the user to set a custom function for handling the name"""
		self.name_type = 'CUSTOM'
		self.name_fun = custom_function
		# Checks the object is sane
		self = self.check_cat_table()
		# Returns the object
		return self

	def set_custom_data_mode(self, custom_function):
		"""Allows the user to set a custom function for handling the name"""
		self.data_type = 'CUSTOM'
		self.name_fun = custom_function
		# Checks the object is sane
		self = self.check_cat_table()
		# Returns the object
		return self

	def set_name_display_value(self, display_text):
		"""Sets custom text for name display"""
		self.name_type = 'SUB'
		self.name_vale = display_text
		# Checks the object is sane
		self = self.check_cat_table()
		# Returns the object
		return self

	def set_name_delimiters(self, delimiter, delim_pos):
		"""Sets delimiters for splitting names"""
		self.name_delim = delimiter
		self.name_delim_pos = delim_pos
		# Checks the object is sane
		self = self.check_cat_table()
		# Returns the object
		return self

	def get_data(self):
		"""Returns the metadata associated with the data"""
		return self.data, self.samples, self.taxa		

	def get_metadata(self):
		"""Returns the metadata associated with the data"""
		return self.meta

	def set_data(self, data, samples, taxa):
		"""Adds data to the object"""
		self.data = data
		self.samples = samples
		self.taxa = taxa
		# Checks the object is sane
		self = self.check_cat_table()
		# Returns the object
		return self

	def set_metadta(self, meta):
		"""Adds metadata to the object"""
		self.meta = meta
		# Checks the object is sane
		self = self.check_cat_table()
		# Returns the object
		return self

		
