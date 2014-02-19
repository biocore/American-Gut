#!/usr/bin/env python 
# test_build_cat_table.py

from unittest import TestCase, main
from numpy import array, flipud
from build_cat_table import (check_category,
                             identify_groups,
                             get_category_position,
                             sort_alphabetically,
                             fuzzy_match,
                             identify_sample_group,
                             build_sub_table,
                             sort_samples,
                             sort_categories)

__author__ = "Justine Debelius"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Justine Debelius"]
__license__ = "BSD" 
__version__ = "unversioned"
__maintainer__ = "Justine Debelius"
__email__ = "Justine.Debelius@colorado.edu"

class BuildCatTableTest(TestCase):

    def setUp(self):
        # Sets up the values to test
        self.common_cats = [(u'k__Bacteria', u' p__Firmicutes'),
                            (u'k__Bacteria', u' p__Bacteroidetes'),
                            (u'k__Bacteria', u' p__Proteobacteria'),
                            (u'k__Bacteria', u' p__Actinobacteria'),
                            (u'k__Bacteria', u' p__Verrucomicrobia'),
                            (u'k__Bacteria', u' p__Tenericutes'),
                            (u'k__Bacteria', u' p__Cyanobacteria'),
                            (u'k__Bacteria', u' p__Fusobacteria')]

        self.sample_ids = ['00009223', '00009814', '00001459', '00008071', 
                           '00007117', '00009112', '00005455', '00008247']

        self.data = array([[0.4738, 0.5646, 0.6382, 0.6170, 0.5180, 0.5609, 
                       0.6557, 0.5105],
                      [0.3755, 0.3670, 0.2232, 0.3114, 0.3991, 0.3434, 
                       0.2122, 0.2694],
                      [0.0801, 0.0099, 0.0135, 0.0330, 0.0821, 0.0135, 
                       0.0675, 0.0872],
                      [0.0511, 0.0448, 0.0239, 0.0000, 0.0000, 0.0365,
                       0.0344, 0.0376],
                      [0.0159, 0.0000, 0.0249, 0.0000, 0.0000, 0.0085, 
                       0.0025, 0.0000],
                      [0.0036, 0.0137, 0.0000, 0.0200, 0.0000, 0.0065, 
                       0.0041, 0.0072],
                      [0.0000, 0.0000, 0.0089, 0.0081, 0.0008, 0.0036, 
                       0.0055, 0.0038],
                      [0.0000, 0.0000, 0.0676, 0.0105, 0.0000, 0.0270, 
                       0.0181, 0.0842]])

        self.mapping = {'00009223':{'#SAMPLEID':'00009223', 'SEX':'male', 
                               'VERSE':'Marvel', 'AGE':'40', 'POWER':'No',
                               'HOME':'New York', 'SERIES':'Avengers'}, 
                        '00009814':{'#SAMPLEID':'00009814', 'SEX':'female',
                               'VERSE':'Marvel', 'AGE':'35', 'POWER':'Yes',
                               'HOME':'New York', 'SERIES':'Avengers'}, 
                        '00001459':{'#SAMPLEID':'00001459', 'SEX':'female', 
                               'VERSE':'Wheedon', 'AGE':'38', 'POWER':'No',
                               'HOME':'Serenity', 'SERIES':'Firefly'}, 
                        '00008071':{'#SAMPLEID':'00008071', 'SEX':'male', 
                               'VERSE':'DC', 'AGE':'19', 'POWER':'NA', 
                               'HOME':'Central City', 'SERIES':'Arrow'}, 
                        '00007117':{'#SAMPLEID':'00007117', 'SEX':'female', 
                               'VERSE':'DC', 'AGE':'22', 'POWER':'No', 
                               'HOME':'Starling City', 'SERIES':'Arrow'}, 
                        '00009112':{'#SAMPLEID':'00009112', 'SEX':'male', 
                               'VERSE':'Wheedon', 'AGE':'30', 'POWER':'No',
                               'HOME':'Serenity', 'SERIES':'Firefly'}, 
                        '00005455':{'#SAMPLEID':'00005455', 'SEX':'male', 
                               'VERSE':'Marvel', 'AGE':'60', 'POWER':'Yes',
                               'HOME':'Graymalkin Lane', 'SERIES':'X-Men'}, 
                        '00008247':{'#SAMPLEID':'00008247', 'SEX':'male', 
                               'VERSE':'Marvel', 'AGE':'15', 'POWER':'Yes',
                               'HOME':'Graymalkin Lane', 'SERIES':'X-Men'}}

    def test_check_category(self):
        """Tests that check_category throws sane errors"""
        # Sets up test veriables
        category_exists = 'VERSE'
        category_absent = 'ALIGNMENT'
        oneDdict = self.mapping['00009223']
        missing_cat = self.mapping
        missing_cat['00009223'] = {'#SAMPLEID':'00009223', 'SEX':'male'}
        # Checks that an error is raised when mapping is not a 2D dictionary
        with self.assertRaises(TypeError):
            check_category(mapping='self.mapping', category=category_exists)
        with self.assertRaises(TypeError):
            check_category(mapping=oneDdict, category=category_exists)
        
        # Checks an error is thrown when category is not a string
        with self.assertRaises(TypeError):
            check_category(mapping=self.mapping, category=[category_exists])
        
        # Checks that an error gets thrown when a mapping category is wrong.
        with self.assertRaises(ValueError):
            check_category(mapping=missing_cat, category=category_exists)
        
        # Checks an error is thrown when category is not in mapping
        with self.assertRaises(ValueError):
            check_category(mapping=self.mapping, category=category_absent)

    def test_identify_groups(self):
        # Sets up the known value
        known_groups = set(['Marvel', 'Wheedon', 'DC'])
        category = 'VERSE'

        # Checks that an error is raised on a problematic input
        with self.assertRaises(ValueError):
            identify_groups(mapping = self.mapping, 
                            category = 'NAME')

        # Identifies the groups from the supplied category
        test_groups = identify_groups(mapping = self.mapping, 
                                      category = category)

        self.assertEqual(test_groups, known_groups)

    def test_get_category_position(self):
        # Sets up test inputs
        test_proteo_cat = 'Proteo'
        test_list_all_cats = ['Castor', 'Pollux', 'Regina', 'Garfield']
        test_list_cat = 'Poll'
        # Sets up known values
        known_pos_proteo = 2
        known_pos_poll = 1

        # Tests error assertions
        with self.assertRaises(ValueError):
            get_category_position(self.common_cats, 'Bacteria')

        # Calculates the test values
        test_pos_proteo = get_category_position(all_cats = self.common_cats, 
                                                cat_descr = test_proteo_cat)
        test_pos_poll = get_category_position(all_cats = test_list_all_cats,
                                              cat_descr = test_list_cat)

        # Test the calculated values
        self.assertEqual(test_pos_proteo, known_pos_proteo)
        self.assertEqual(test_pos_poll, known_pos_poll)

    def test_sort_alphabetically(self):
        # Defines the test input
        test_unord = ['Rose', 'Jack', 'Martha', 'Donna', 'River', 'Amy', 'Rory', 'Clara']

        # Defines the known values
        known_alpha = ['Amy', 'Clara', 'Donna', 'Jack', 'Martha', 'River', 'Rory', 'Rose']
        known_order = [5, 7, 3, 1, 2, 4, 6, 0]

        # Caclulates the test values
        (test_alpha, test_order) = sort_alphabetically(unord = test_unord)

        # Checks the known and test values are equal
        self.assertEqual(known_alpha, test_alpha)
        self.assertEqual(known_order, test_order)
       
    def test_fuzzy_match(self):
        # Defines values to test
        target_tuple = 'Fuso'
        target_string = 'Don'
        possible_strings = ['Rose', 'Jack', 'Martha', 'Donna', 'River', 'Amy', 
                           'Rory', 'Clara']
        possible_integers = [0, 3, 2, 4]

        # Defines the known value
        known_pos_tup = 7
        known_match_tup = (u'k__Bacteria', u' p__Fusobacteria')
        known_pos_str = 3
        known_match_str = 'Donna'

        # Checks that errors are called appropriately
        with self.assertRaises(ValueError):
            fuzzy_match(target = 'Wilf', 
                        possible_matches = possible_strings)

        with self.assertRaises(ValueError):
            fuzzy_match(target = 'Bacteria', 
                        possible_matches = self.common_cats)

        with self.assertRaises(TypeError):
            fuzzy_match(target = 3,
                        possible_matches = possible_integers)

        # Calculates the test values
        [pos_tup, match_tup] = fuzzy_match(target = target_tuple, 
                                           possible_matches = self.common_cats)

        [pos_str, match_str] = fuzzy_match(target = target_string, 
                                           possible_matches = possible_strings)

        # Checks the values are sane
        self.assertEqual(pos_tup, known_pos_tup)
        self.assertEqual(pos_str, known_pos_str)
        self.assertEqual(match_tup, known_match_tup)
        self.assertEqual(match_str, known_match_str)

    def test_identify_sample_group(self):
        # Sets up values for testing
        category = 'VERSE'
        wrong_sample_ids = ['00100', '00009112']

        # Sets up the known value
        known_groups = {'Marvel':  ['00009223', '00009814', '00005455', 
                                    '00008247'],
                        'DC':      ['00008071', '00007117'],
                        'Wheedon': ['00001459', '00009112']}

        # Test that errors are appropriately raised
        with self.assertRaises(ValueError):
            identify_sample_group(sample_ids = wrong_sample_ids, 
                                  mapping = self.mapping, 
                                  category = category)

        # Calculates the test value
        test_groups = identify_sample_group(sample_ids = self.sample_ids, 
                                            mapping = self.mapping, 
                                            category = category)

        # Verifies the outputs are the sample
        self.assertEqual(test_groups, known_groups)

    def test_build_sub_table(self):
        # Sets up the values for testing
        test_ids = ['00008071', '00007117']
        # Sets up the known value
        known_data_out = self.data[:,[3,4]]

        # Runs a sanity check when incorrect values are passed
        with self.assertRaises(ValueError):
            build_sub_table(data = self.data, 
                            sample_ids = test_ids, 
                            target_ids = test_ids)

        with self.assertRaises(TypeError):
            build_sub_table(data=self.data,
                            sample_ids= test_ids,
                            target_ids = 3)

        with self.assertRaises(ValueError):
            build_sub_table(data = self.data, 
                            sample_ids = self.sample_ids,
                            target_ids = ['Cat', 'Dog'])

        # Calculates the test value
        test_data_out = build_sub_table(data = self.data,
                                        sample_ids = self.sample_ids,
                                        target_ids = test_ids)

        # Checks the output and the known are equal
        self.assertTrue((known_data_out == test_data_out).all())

    def test_sort_samples(self):
        # Sets up expected output values
        known_ids_alpha = ['00001459', '00005455', '00007117', '00008071',
                           '00008247', '00009112', '00009223', '00009814']
        
        known_data_alpha = array([[0.6382, 0.6557, 0.5180, 0.6170, 0.5105, 
                                   0.5609, 0.4738, 0.5646],
                                  [0.2232, 0.2122, 0.3991, 0.3114, 0.2694, 
                                   0.3434, 0.3755, 0.3670],
                                  [0.0135, 0.0675, 0.0821, 0.0330, 0.0872, 
                                   0.0135, 0.0801, 0.0099],
                                  [0.0239, 0.0344, 0.0000, 0.0000, 0.0376, 
                                   0.0365, 0.0511, 0.0448],
                                  [0.0249, 0.0025, 0.0000, 0.0000, 0.0000, 
                                   0.0085, 0.0159, 0.0000],
                                  [0.0000, 0.0041, 0.0000, 0.0200, 0.0072, 
                                   0.0065, 0.0036, 0.0137],
                                  [0.0089, 0.0055, 0.0008, 0.0081, 0.0038, 
                                   0.0036, 0.0000, 0.0000],
                                  [0.0676, 0.0181, 0.0000, 0.0105, 0.0842, 
                                   0.0270, 0.0000, 0.0000]])
        
        known_ids_row_3 = ['00009814', '00001459', '00009112', '00008071',
                           '00005455', '00009223', '00007117', '00008247']

        known_data_row_3 = array([[0.5646, 0.6382, 0.5609, 0.6170, 0.6557,
                                   0.4738, 0.5180, 0.5105],
                                  [0.3670, 0.2232, 0.3434, 0.3114, 0.2122,
                                   0.3755, 0.3991, 0.2694],
                                  [0.0099, 0.0135, 0.0135, 0.0330, 0.0675,
                                   0.0801, 0.0821, 0.0872],
                                  [0.0448, 0.0239, 0.0365, 0.0000, 0.0344,
                                   0.0511, 0.0000, 0.0376],
                                  [0.0000, 0.0249, 0.0085, 0.0000, 0.0025,
                                   0.0159, 0.0000, 0.0000],
                                  [0.0137, 0.0000, 0.0065, 0.0200, 0.0041,
                                   0.0036, 0.0000, 0.0072],
                                  [0.0000, 0.0089, 0.0036, 0.0081, 0.0055,
                                   0.0000, 0.0008, 0.0038],
                                  [0.0000, 0.0676, 0.0270, 0.0105, 0.0181,
                                   0.0000, 0.0000, 0.0842]])



        # Checks errors are raised when data is not appropriate passed
        with self.assertRaises(ValueError):
            sort_samples(data = self.data, 
                         sample_ids = self.sample_ids[:3])

        with self.assertRaises(ValueError):
            sort_samples(data = self.data, 
                         sample_ids = self.sample_ids,
                         sort_key = 12)

        # Generates outputs using test data
        [t_alpha_ids, t_alpha_data] = sort_samples(data = self.data,
                                                   sample_ids = self.sample_ids,
                                                   sort_key = None)

        [t_rev_ids, t_rev_data] = sort_samples(data = self.data,
                                               sample_ids = self.sample_ids,
                                               reverse = True)

        [t_row_ids, t_row_data] = sort_samples(data = self.data,
                                               sample_ids = self.sample_ids,
                                               sort_key = 2)

        # Checks the outputs are correct        
        self.assertEqual(known_ids_alpha, t_alpha_ids)
        self.assertTrue((known_data_alpha==t_alpha_data).all())

        self.assertEqual(known_ids_alpha[::-1], t_rev_ids)
        self.assertTrue((known_data_alpha[:,::-1] == t_rev_data).all())
        
        self.assertEqual(known_ids_row_3, t_row_ids)
        self.assertTrue((known_data_row_3 == t_row_data).all())

    def test_sort_categories(self):
        # Defines the values to test
        first_cat = 'Proteobacteria'
        sort_method_2 = 'ALPHA'
        sort_method_3 = 'RETAIN'
        sort_method_4 = 'CUSTOM'
        sort_order = ['Proteo', 'Firmicutes', 'detes', 'Tenericutes', 
                      'Actino', 'Cyano', 'Prot']

        # Defines the known values
        def_cats = [(u'k__Bacteria', u' p__Proteobacteria'), 
                     (u'k__Bacteria', u' p__Tenericutes'), 
                     (u'k__Bacteria', u' p__Verrucomicrobia'), 
                     (u'k__Bacteria', u' p__Actinobacteria'), 
                     (u'k__Bacteria', u' p__Bacteroidetes'), 
                     (u'k__Bacteria', u' p__Cyanobacteria'), 
                     (u'k__Bacteria', u' p__Firmicutes'),
                     (u'k__Bacteria', u' p__Fusobacteria')]

        m_2_cats = [(u'k__Bacteria', u' p__Proteobacteria'), 
                     (u'k__Bacteria', u' p__Actinobacteria'), 
                     (u'k__Bacteria', u' p__Bacteroidetes'), 
                     (u'k__Bacteria', u' p__Cyanobacteria'), 
                     (u'k__Bacteria', u' p__Firmicutes'),
                     (u'k__Bacteria', u' p__Fusobacteria'),
                     (u'k__Bacteria', u' p__Tenericutes'), 
                     (u'k__Bacteria', u' p__Verrucomicrobia')]

        m_3_cats = [(u'k__Bacteria', u' p__Proteobacteria'),
                    (u'k__Bacteria', u' p__Firmicutes'),
                    (u'k__Bacteria', u' p__Bacteroidetes'),                   
                    (u'k__Bacteria', u' p__Actinobacteria'),
                    (u'k__Bacteria', u' p__Verrucomicrobia'),
                    (u'k__Bacteria', u' p__Tenericutes'),
                    (u'k__Bacteria', u' p__Cyanobacteria'),
                    (u'k__Bacteria', u' p__Fusobacteria')]

        m_4_cats = [(u'k__Bacteria', u' p__Proteobacteria'),
                    (u'k__Bacteria', u' p__Firmicutes'),
                    (u'k__Bacteria', u' p__Bacteroidetes'),
                    (u'k__Bacteria', u' p__Tenericutes'),
                    (u'k__Bacteria', u' p__Actinobacteria'),
                    (u'k__Bacteria', u' p__Cyanobacteria'),
                    (u'k__Bacteria', u' p__Proteobacteria')]

        def_data = array([[0.0801, 0.0099, 0.0135, 0.0330, 0.0821, 0.0135, 
                           0.0675, 0.0872],
                          [0.0036, 0.0137, 0.0000, 0.0200, 0.0000, 0.0065, 
                           0.0041, 0.0072],
                          [0.0159, 0.0000, 0.0249, 0.0000, 0.0000, 0.0085, 
                           0.0025, 0.0000],
                          [0.0511, 0.0448, 0.0239, 0.0000, 0.0000, 0.0365, 
                           0.0344, 0.0376],
                          [0.3755, 0.3670, 0.2232, 0.3114, 0.3991, 0.3434, 
                           0.2122, 0.2694],
                          [0.0000, 0.0000, 0.0089, 0.0081, 0.0008, 0.0036, 
                           0.0055, 0.0038],
                          [0.4738, 0.5646, 0.6382, 0.6170, 0.5180, 0.5609, 
                           0.6557, 0.5105],
                          [0.0000, 0.0000, 0.0676, 0.0105, 0.0000, 0.0270, 
                           0.0181, 0.0842]])
                                
        rev_data = flipud(def_data)

        m_2_data = array([[0.0801, 0.0099, 0.0135, 0.0330, 0.0821, 0.0135, 
                           0.0675, 0.0872],
                          [0.0511, 0.0448, 0.0239, 0.0000, 0.0000, 0.0365, 
                           0.0344, 0.0376],
                          [0.3755, 0.3670, 0.2232, 0.3114, 0.3991, 0.3434, 
                           0.2122, 0.2694],
                          [0.0000, 0.0000, 0.0089, 0.0081, 0.0008, 0.0036, 
                           0.0055, 0.0038],
                          [0.4738, 0.5646, 0.6382, 0.6170, 0.5180, 0.5609, 
                           0.6557, 0.5105],
                          [0.0000, 0.0000, 0.0676, 0.0105, 0.0000, 0.0270, 
                           0.0181, 0.0842],
                          [0.0036, 0.0137, 0.0000, 0.0200, 0.0000, 0.0065, 
                           0.0041, 0.0072],
                          [0.0159, 0.0000, 0.0249, 0.0000, 0.0000, 0.0085, 
                           0.0025, 0.0000]])

        
        m_3_data = array([[0.0801, 0.0099, 0.0135, 0.0330, 0.0821, 0.0135, 
                           0.0675, 0.0872],
                          [0.4738, 0.5646, 0.6382, 0.6170, 0.5180, 0.5609, 
                           0.6557, 0.5105],
                          [0.3755, 0.3670, 0.2232, 0.3114, 0.3991, 0.3434, 
                           0.2122, 0.2694],
                          [0.0511, 0.0448, 0.0239, 0.0000, 0.0000, 0.0365, 
                           0.0344, 0.0376],
                          [0.0159, 0.0000, 0.0249, 0.0000, 0.0000, 0.0085, 
                           0.0025, 0.0000],
                          [0.0036, 0.0137, 0.0000, 0.0200, 0.0000, 0.0065, 
                           0.0041, 0.0072],
                          [0.0000, 0.0000, 0.0089, 0.0081, 0.0008, 0.0036, 
                           0.0055, 0.0038],
                          [0.0000, 0.0000, 0.0676, 0.0105, 0.0000, 0.0270, 
                           0.0181, 0.0842]])

        m_4_data = array([[0.0801, 0.0099, 0.0135, 0.0330, 0.0821, 0.0135, 
                           0.0675, 0.0872],
                          [0.4738, 0.5646, 0.6382, 0.6170, 0.5180, 0.5609, 
                           0.6557, 0.5105],
                          [0.3755, 0.3670, 0.2232, 0.3114, 0.3991, 0.3434, 
                           0.2122, 0.2694],
                          [0.0036, 0.0137, 0.0000, 0.0200, 0.0000, 0.0065, 
                           0.0041, 0.0072],
                          [0.0511, 0.0448, 0.0239, 0.0000, 0.0000, 0.0365, 
                           0.0344, 0.0376],
                          [0.0000, 0.0000, 0.0089, 0.0081, 0.0008, 0.0036, 
                           0.0055, 0.0038],
                          [0.0801, 0.0099, 0.0135, 0.0330, 0.0821, 0.0135, 
                           0.0675, 0.0872]])

        # Gets the test value
        [test_def_cats, test_def_data] = sort_categories(self.data, 
            self.common_cats, first_cat)

        [test_rev_cats, test_rev_data] = sort_categories(self.data, 
            self.common_cats, first_cat, reverse = True)       

        [test_m_2_cats, test_m_2_data] = sort_categories(self.data, 
            self.common_cats, first_cat, sort_method = sort_method_2)

        [test_m_3_cats, test_m_3_data] = sort_categories(self.data, 
            self.common_cats, first_cat, sort_method = sort_method_3)

        [test_m_4_cats, test_m_4_data] = sort_categories(self.data, 
            self.common_cats, first_cat, sort_method = sort_method_4, 
            category_order = sort_order)
        

        # Checks the categories are coming out correctly ordered
        self.assertEqual(test_def_cats, def_cats)
        self.assertEqual(test_rev_cats, def_cats[::-1])
        self.assertEqual(test_m_2_cats, m_2_cats)
        self.assertEqual(test_m_3_cats, m_3_cats)
        self.assertEqual(test_m_4_cats, m_4_cats)

        self.assertTrue((test_def_data == def_data).all())
        self.assertTrue((test_rev_data == rev_data).all())
        self.assertTrue((test_m_2_data == m_2_data).all())
        self.assertTrue((test_m_3_data == m_3_data).all())
        self.assertTrue((test_m_4_data == m_4_data).all())

if __name__ == '__main__':
    main()