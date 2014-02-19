#!/usr/bin/env python
from __future__ import division
from unittest import TestCase, main
from summarize_meta import get_meta_counts, convert_counts, compare_counts
from numpy import array


class SummarizeMetaTest(TestCase):

	def setUp(self):
		"""Sets up the data set for testing"""
		# Sets up a mapping file
		self.mapping_1 = {'Harry':    {'SEX':'Male', 'HOUSE':'Gryffendor', 
								       'PARTNER':'Ginny', 'PET':'Hedwig', 
								       'LIVING_PARENTS': 'NA', 'AGE':'17'},
						   'Ron':     {'SEX':'Male', 'HOUSE':'Gryffendor', 
					 		    	   'PARTNER':'Hermione', 'PET':'Scabbers',
					 		     	   'LIVING_PARENTS': 'Arthur and Molly', 
					 		     	   'AGE':'17'},
						   'Hermione':{'SEX':'Female', 'HOUSE':'Gryffendor', 
					 			       'PARTNER':'Ron', 'PET':'Crookshanks',
					 			       'LIVING_PARENTS': 'Helen and Unknown',
					 			       'AGE':'18'}}
		self.mapping_2 = {'Stark':    {'SEX':'Male', 'VERSE':'Marvel', 'AGE':
									   '40', 'POWER':'No', 'HOME':'New York', 
									   'SERIES':'Avengers'}, 
		                  'Romanov':  {'SEX':'Female', 'VERSE':'Marvel', 'AGE':
		                  			   'NA', 'POWER':'Yes', 'HOME':'New York', 
		                  			   'SERIES':'Avengers'}, 
	                      'Allen':    {'SEX':'Male', 'VERSE':'DC', 'AGE':'19', 
	                      			   'POWER':'NA', 'HOME':'Central City', 
	                      			   'SERIES':'Arrow'}}
		self.maps = [self.mapping_1, self.mapping_2]
		self.category = 'SEX'
		self.na_val = 'NA'
		self.counts = {'Male': 2, 'Female': 1}
		self.normalize = False
		self.ordering = ['Female', 'Male']


	def test_get_meta_counts(self):
		"""Tests that get_meta_count is sane"""
		# Checks that groups returns correctly
		known_value = {'Male':2, 'Female':1}
		test_value = get_meta_counts(mapping=self.mapping_1, 
									 category=self.category)
		self.assertEqual(known_value, test_value)

	def test_convert_counts(self):
		"""Tests that convert_counts is sane"""
		## Checks error handling
		# Tests that an error is called when counts is not a dict.
		with self.assertRaises(TypeError):
			convert_counts('self.counts')
		# Checks error calls for the handling of NA is sane. 
		with self.assertRaises(TypeError):
			convert_counts(self.counts, na_val=[self.na_val])
		with self.assertRaises(TypeError):
			convert_counts(self.counts, ignore_nas=self.na_val)
		with self.assertRaises(ValueError):
			convert_counts(self.counts, ignore_nas=True, na_val=None)
		# Checks normalize class handling is sane
		with self.assertRaises(TypeError):
			convert_counts(self.counts, normalize=[self.normalize])
		# Checks error calls for ordering handling is sane
		with self.assertRaises(TypeError):
			convert_counts(self.counts, ordering=3)
		with self.assertRaises(ValueError):
			convert_counts(self.counts, ordering='self.ordering')
		with self.assertRaises(ValueError):
			convert_counts(self.counts, ordering=['self.ordering'])

		## Tests ordering modes
		# Tests no ordering
		known_groups = ['Male', 'Female']
		known_counts = array([[2, 1]])
		[test_groups, test_counts] = convert_counts(self.counts)
		self.assertEqual(test_groups, known_groups)
		self.assertTrue((test_counts==known_counts).all())
		# Alphabetical sorting
		ordering='ALPHA'
		known_groups = ['Female', 'Male']
		known_counts = array([[1, 2]])
		[test_groups, test_counts] = convert_counts(self.counts, 
													ordering=ordering)
		self.assertEqual(test_groups, known_groups)
		self.assertTrue((test_counts==known_counts).all())

		# Reverse alphabetical sorting
		ordering='ALPHA_REV'
		known_groups = ['Male', 'Female']
		known_counts = array([[2, 1]])
		[test_groups, test_counts] = convert_counts(self.counts, 
													ordering=ordering)
		self.assertEqual(test_groups, known_groups)
		self.assertTrue((test_counts==known_counts).all())
		# Frequency sorting
		ordering = 'FREQ'
		known_groups = ['Female', 'Male']
		known_counts = array([[1, 2]])
		[test_groups, test_counts] = convert_counts(self.counts, 
													ordering=ordering)
		self.assertEqual(test_groups, known_groups)
		self.assertTrue((test_counts==known_counts).all())
		# Reverse Frequency sorting 
		ordering = 'FREQ_REV'
		known_groups = ['Male', 'Female']
		known_counts = array([[2, 1]])
		[test_groups, test_counts] = convert_counts(self.counts, 
													ordering=ordering)
		self.assertEqual(test_groups, known_groups)
		self.assertTrue((test_counts==known_counts).all())
		# Tests list ordering
		known_groups = ['Female', 'Male']
		known_counts = array([[1, 2]])
		[test_groups, test_counts] = convert_counts(self.counts, 
													ordering=self.ordering)
		self.assertEqual(test_groups, known_groups)
		self.assertTrue((test_counts==known_counts).all())

		## Checks the normalization argument
		known_groups = ['Male', 'Female']
		known_counts = array([[2/3, 1/3]])
		[test_groups, test_counts] = convert_counts(self.counts, 
													normalize=True)
		self.assertEqual(test_groups, known_groups)
		self.assertTrue((test_counts==known_counts).all())

		## Checks the no category skipping works
		counts = get_meta_counts(self.mapping_1, category='LIVING_PARENTS')
		known_groups = ['Arthur and Molly', 'Helen and Unknown']
		known_counts = array([[1, 1]])
		[test_groups, test_counts] = convert_counts(counts, 
													na_val='NA', 
													ignore_nas=True)
		self.assertEqual(test_groups, known_groups)
		self.assertTrue((test_counts==known_counts).all())

	def test_compare_counts(self):
		"""Checks that compare counts renders sane outputs"""
		## Preforms error checking
		# Checks that maps is handled sanely
		with self.assertRaises(TypeError):
			compare_counts(self.mapping_1, self.category)
		# Checks that NA handling is sane
		with self.assertRaises(TypeError):
			compare_counts(self.maps, self.category, na_val=[self.na_val])
		with self.assertRaises(TypeError):
			compare_counts(self.maps, self.category, ignore_nas=self.na_val)
		with self.assertRaises(ValueError):
			compare_counts(self.maps, self.category, ignore_nas=True)
		# Checks normal class handling is sane
		with self.assertRaises(TypeError):
			compare_counts(self.maps, self.category, normalize=[self.normalize])
		# Checks error calls for order handling is sane
		with self.assertRaises(TypeError):
			compare_counts(self.maps, self.category, ordering=3)
		with self.assertRaises(ValueError):
			compare_counts(self.maps, self.category, ordering='self.ordering')

		## Checks default case handling
		known_groups = ['Male', 'Female']
		known_counts = array([[2, 2],
							  [1, 1]])
		[test_groups, test_counts] = compare_counts(self.maps, self.category)
		self.assertEqual(test_groups, known_groups)
		self.assertTrue((test_counts ==known_counts).all())

		## Checks ordering modes
		# Alphabetical ordering
		ordering='ALPHA'
		known_groups = ['Female', 'Male']
		known_counts = array([[1, 1],
							  [2, 2]])
		[test_groups, test_counts] = compare_counts(self.maps, self.category, 
													ordering=ordering)
		self.assertEqual(test_groups, known_groups)
		self.assertTrue((test_counts ==known_counts).all())
		# Reverse alphabetical sorting
		ordering = ordering='ALPHA_REV'
		known_groups = ['Male', 'Female']
		known_counts = array([[2, 2],
							  [1, 1]])
		[test_groups, test_counts] = compare_counts(self.maps, self.category, 
													ordering=ordering)
		self.assertEqual(test_groups, known_groups)
		self.assertTrue((test_counts ==known_counts).all())
		# Frequency based sorting
		ordering='FREQ'
		known_groups = ['Female', 'Male']
		known_counts = array([[1, 1],
							  [2, 2]])
		[test_groups, test_counts] = compare_counts(self.maps, self.category, 
													ordering=ordering)
		self.assertEqual(test_groups, known_groups)
		self.assertTrue((test_counts ==known_counts).all())
		# Reverse alphabetical sorting
		ordering = ordering='FREQ_REV'
		known_groups = ['Male', 'Female']
		known_counts = array([[2, 2],
							  [1, 1]])
		[test_groups, test_counts] = compare_counts(self.maps, self.category, 
													ordering=ordering)
		self.assertEqual(test_groups, known_groups)
		self.assertTrue((test_counts ==known_counts).all())
		# Tests list ordering
		known_groups = ['Female', 'Male']
		known_counts = array([[1, 1],
							  [2, 2]])
		[test_groups, test_counts] = compare_counts(self.maps, self.category, 
													ordering=self.ordering)
		self.assertEqual(test_groups, known_groups)
		self.assertTrue((test_counts ==known_counts).all())

		## Checks the normalization argument
		known_groups = ['Male', 'Female']
		known_counts = array([[2/3, 2/3],
							  [1/3, 1/3]])
		[test_groups, test_counts] = compare_counts(self.maps, self.category,
													normalize=True)
		print test_counts
		self.assertEqual(test_groups, known_groups)
		self.assertTrue((test_counts ==known_counts).all())

		## Checks the category skipping argument works
		category = 'AGE'
		known_groups = ['19', '18', '40', '17']
		known_counts = array([[0, 1], [1, 0], [0, 1], [2, 0]])
		[test_groups, test_counts] = compare_counts(self.maps, category,
													na_val=self.na_val,
													ignore_nas=True)
		self.assertEqual(test_groups, known_groups)
		self.assertTrue((test_counts ==known_counts).all())

		




if __name__ == '__main__':
	main()

