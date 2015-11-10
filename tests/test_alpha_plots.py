#!/usr/bin/env python

from unittest import TestCase, main

import numpy as np
import pandas as pd

import numpy.testing as npt

import pandas.util.testing as pdt

from americangut.alpha_plots import (horizontal_to_longtable,
                                     modify_alpha_shannon_pd,
                                     update_map,
                                     ag_summary_violin)

__author__ = "Justine Debelius"
__copyright__ = "Copyright 2015"
__credits__ = ["Justine Debelius"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Justine Debelius"
__email__ = "jdebelius@ucsd.edu"


class AlphaPlotTest(TestCase):

    def setUp(self):
        # Sets up the original dataframe with the data
        data = np.array([['Gryffendor', 'unknown', 'M',
                          'UBERON:feces', '20s', 'Occasionally', 'Never',
                          'Overweight', 'Rarely', 'Right'],
                         ['Gryffendor', 'yes', 'M', 'UBERON:feces',
                          '20s', 'Occasionally', 'Never', 'Normal',
                          'Regularly', 'Right'],
                         ['Gryffendor', 'no', 'M', 'UBERON:nose',
                          '20s', 'Rarely', 'Never', 'Overweight', 'Regularly',
                          'Right'],
                         ['Gryffendor', 'yes', 'M', 'UBERON:feces',
                          'teens', 'Occasionally', 'Never', 'Normal', 'Daily',
                          'Left'],
                         ['Gryffendor', 'yes', 'M', 'UBERON:feces',
                          'teens', 'Occasionally', 'Never', 'Normal', 'Rarely',
                          'Right'],
                         ['Gryffendor', 'yes', 'M', 'UBERON:feces',
                          'teens', 'Never', 'Never', 'Normal', 'Never',
                          'Left'],
                         ['Gryffendor', 'yes', 'F', 'UBERON:feces',
                          'teens', 'Never', 'Never', 'Normal', 'Regularly',
                          'Right'],
                         ['Gryffendor', 'yes', 'M', 'UBERON:skin',
                          'teens', 'Rarely', 'Last year', 'Underweight',
                          'Never', 'Right'],
                         ['Gryffendor', 'no', 'F',
                          'UBERON:oral cavity', 'teens', 'Never', 'Last year',
                          'Unknown', 'Daily', 'Right']
                         ])
        columns = ['House', 'Quiddich', 'Sex', 'BODY_HABITAT',
                   'AGE_CAT', 'ALCOHOL_FREQUENCY', 'ANTIBIOTIC_HISTORY',
                   'BMI_CAT', 'FLOSSING_FREQUENCY', 'DOMINANT_HAND']
        samples = ['Bill', 'Charlie', 'Percy', 'Fred', 'George', 'Ron',
                   'Ginny', 'Harry', 'Hermione']
        self.wide = pd.DataFrame(data, columns=columns, index=samples)
        for column in self.wide.columns:
            self.wide[column] = self.wide[column].astype(str)

        # Identifies the sample, categories, and metrics to test
        self.sample = 'Ron'
        self.categories = ['Quiddich', 'Sex', 'House']
        self.metrics = ['shannon_1000', 'pd_whole_tree_1000']

        # Sets up the alpha diversity results
        a_columns = ['Unnamed: 0', 'sequences per sample', 'iteration',
                     'Bill', 'Charlie', 'Percy', 'Fred', 'George', 'Ron',
                     'Ginny', 'Hermione', 'Harry']
        pd_whole_tree = np.array([
            ['alpha_rarefaction_1000_0.txt', 1000, 0, 15.069, 14.664, 13.741,
             23.810, 21.033, 20.297, 29.491, 30.203, 10.239],
            ['alpha_rarefaction_1000_1.txt', 1000, 1, 14.268, 14.978, 13.786,
             25.398, 21.023, 19.611, 26.715, 31.302,  8.203],
            ['alpha_rarefaction_1000_2.txt', 1000, 2, 13.126, 15.757, 14.824,
             26.489, 19.735, 18.164, 30.083, 32.203,  9.203],
            ['alpha_rarefaction_1000_3.txt', 1000, 3, 15.465, 14.825, 15.007,
             27.416, 20.634, 18.953, 29.311, 31.032,  9.293],
            ['alpha_rarefaction_1000_4.txt', 1000, 4, 15.415, 15.539, 14.763,
             27.604, 19.717, 17.380, 28.958, 29.358,  8.203],
            ['alpha_rarefaction_1000_5.txt', 1000, 5, 15.244, 14.917, 14.403,
             24.731, 21.128, 18.099, 27.847, 29.028,  10.203],
            ['alpha_rarefaction_1000_6.txt', 1000, 6, 15.771, 16.578, 13.443,
             26.693, 20.936, 19.648, 27.476, 31.305,  8.302],
            ['alpha_rarefaction_1000_7.txt', 1000, 7, 13.845, 15.243, 15.243,
             23.744, 20.473, 20.187, 28.257, 32.305,  10.393],
            ['alpha_rarefaction_1000_8.txt', 1000, 8, 14.128, 15.365, 14.587,
             25.076, 21.026, 18.719, 27.933, 27.302,   8.294],
            ['alpha_rarefaction_1000_9.txt', 1000, 9, 15.078, 14.312, 14.607,
             25.415, 21.100, 19.108, 28.811, 28.302,   7.201]
            ], dtype=object)
        pd_whole_tree = pd.DataFrame(pd_whole_tree, columns=a_columns)
        shannon = np.array([
            ['alpha_rarefaction_1000_0.txt', 1000, 0, 5.045, 4.958, 4.582,
             6.206, 5.039, 5.343, 6.987, 5.324, 8.302],
            ['alpha_rarefaction_1000_1.txt', 1000, 1, 5.152, 5.086, 4.56,
             6.270, 5.157, 5.165, 6.751, 5.235, 8.320],
            ['alpha_rarefaction_1000_2.txt', 1000, 2, 5.093, 5.054, 4.633,
             6.334, 5.022, 5.182, 6.922, 3.181, 8.205],
            ['alpha_rarefaction_1000_3.txt', 1000, 3, 5.126, 4.971, 4.496,
             6.290, 5.02, 5.173, 6.897,  4.173, 8.201],
            ['alpha_rarefaction_1000_4.txt', 1000, 4, 5.028, 5.233, 4.640,
             6.368, 5.05, 5.145, 6.842,  6.273, 8.105],
            ['alpha_rarefaction_1000_5.txt', 1000, 5, 5.244, 5.097, 4.558,
             6.344, 5.139, 5.059, 6.884, 6.582, 8.105],
            ['alpha_rarefaction_1000_6.txt', 1000, 6, 5.306, 5.045, 4.616,
             6.360, 4.969, 5.283, 6.766, 5.239, 7.924],
            ['alpha_rarefaction_1000_7.txt', 1000, 7, 5.100, 5.001, 4.696,
             6.054, 5.023, 5.394, 6.811, 5.530, 7.305],
            ['alpha_rarefaction_1000_8.txt', 1000, 8, 5.229, 5.053, 4.738,
             6.425, 5.071, 5.152, 6.818, 6.203, 8.205],
            ['alpha_rarefaction_1000_9.txt', 1000, 9, 5.192, 5.003, 4.529,
             6.299, 5.059, 5.203, 6.903, 3.205, 8.402]
            ], dtype=object)
        shannon = pd.DataFrame(shannon, columns=a_columns)
        self.alpha = {'pd_whole_tree': pd_whole_tree,
                      'shannon': shannon}
        alpha_data = np.array([['Gryffendor', 'unknown', 'M',
                                'UBERON:feces', '20s', 'Occasionally', 'Never',
                                'Overweight', 'Rarely', 'Right', 14.7409,
                                5.1515],
                               ['Gryffendor', 'yes', 'M', 'UBERON:feces',
                                '20s', 'Occasionally', 'Never', 'Normal',
                                'Regularly', 'Right', 15.2178, 5.0501],
                               ['Gryffendor', 'no', 'M', 'UBERON:nose',
                                '20s', 'Rarely', 'Never', 'Overweight',
                                'Regularly', 'Right', 14.4404, 4.6048],
                               ['Gryffendor', 'yes', 'M', 'UBERON:feces',
                                'teens', 'Occasionally', 'Never', 'Normal',
                                'Daily', 'Left', 25.6376, 6.295],
                               ['Gryffendor', 'yes', 'M', 'UBERON:feces',
                                'teens', 'Occasionally', 'Never', 'Normal',
                                'Rarely', 'Right', 20.6805, 5.0549],
                               ['Gryffendor', 'yes', 'M', 'UBERON:feces',
                                'teens', 'Never', 'Never', 'Normal', 'Never',
                                'Left', 19.0166, 5.2099],
                               ['Gryffendor', 'yes', 'F', 'UBERON:feces',
                                'teens', 'Never', 'Never', 'Normal',
                                'Regularly', 'Right', 28.4882, 6.8581],
                               ['Gryffendor', 'yes', 'M', 'UBERON:skin',
                                'teens', 'Rarely', 'Last year', 'Underweight',
                                'Never', 'Right', 8.9534,  8.1074],
                               ['Gryffendor', 'no', 'F',
                                'UBERON:oral cavity', 'teens', 'Never',
                                'Last year', 'Unknown', 'Daily', 'Right',
                                30.2340,  5.0945]
                               ])
        alpha_columns = ['House', 'Quiddich', 'Sex', 'BODY_HABITAT',
                         'AGE_CAT', 'ALCOHOL_FREQUENCY', 'ANTIBIOTIC_HISTORY',
                         'BMI_CAT', 'FLOSSING_FREQUENCY', 'DOMINANT_HAND',
                         'pd_whole_tree_1000', 'shannon_1000']
        self.alpha_map = pd.DataFrame(alpha_data,
                                      columns=alpha_columns,
                                      index=samples)
        for column in ['House', 'Quiddich', 'Sex']:
            self.alpha_map[column] = self.alpha_map[column].astype(str)
        for column in ['pd_whole_tree_1000', 'shannon_1000']:
            self.alpha_map[column] = self.alpha_map[column].astype(float)

        # Sets up the longform table
        long_ = np.array([
            ['Charlie',  'Quiddich',       'yes',     5.0501,  'shannon'],
            ['Fred',     'Quiddich',       'yes',     6.2950,  'shannon'],
            ['George',   'Quiddich',       'yes',     5.0549,  'shannon'],
            ['Ron',      'Quiddich',       'yes',     5.2099,  'shannon'],
            ['Ginny',    'Quiddich',       'yes',     6.8581,  'shannon'],
            ['Harry',    'Quiddich',       'yes',     8.1074,  'shannon'],
            ['Charlie',  'Quiddich',       'yes',    15.2178,       'pd'],
            ['Fred',     'Quiddich',       'yes',    25.6376,       'pd'],
            ['George',   'Quiddich',       'yes',    20.6805,       'pd'],
            ['Ron',      'Quiddich',       'yes',    19.0166,       'pd'],
            ['Ginny',    'Quiddich',       'yes',    28.4882,       'pd'],
            ['Harry',    'Quiddich',       'yes',     8.9534,       'pd'],
            ['Bill',          'Sex',         'M',     5.1515,  'shannon'],
            ['Charlie',       'Sex',         'M',     5.0501,  'shannon'],
            ['Percy',         'Sex',         'M',     4.6048,  'shannon'],
            ['Fred',          'Sex',         'M',     6.2950,  'shannon'],
            ['George',        'Sex',         'M',     5.0549,  'shannon'],
            ['Ron',           'Sex',         'M',     5.2099,  'shannon'],
            ['Harry',         'Sex',         'M',     8.1074,  'shannon'],
            ['Bill',          'Sex',         'M',    14.7409,       'pd'],
            ['Charlie',       'Sex',         'M',    15.2178,       'pd'],
            ['Percy',         'Sex',         'M',    14.4404,       'pd'],
            ['Fred',          'Sex',         'M',    25.6376,       'pd'],
            ['George',        'Sex',         'M',    20.6805,       'pd'],
            ['Ron',           'Sex',         'M',    19.0166,       'pd'],
            ['Harry',         'Sex',         'M',     8.9534,       'pd'],
            ['Bill',        'House',  'Gryffendor',   5.1515,  'shannon'],
            ['Charlie',     'House',  'Gryffendor',   5.0501,  'shannon'],
            ['Percy',       'House',  'Gryffendor',   4.6048,  'shannon'],
            ['Fred',        'House',  'Gryffendor',   6.2950,  'shannon'],
            ['George',      'House',  'Gryffendor',   5.0549,  'shannon'],
            ['Ron',         'House',  'Gryffendor',   5.2099,  'shannon'],
            ['Ginny',       'House',  'Gryffendor',   6.8581,  'shannon'],
            ['Harry',       'House',  'Gryffendor',   8.1074,  'shannon'],
            ['Hermione',    'House',  'Gryffendor',   5.0945,  'shannon'],
            ['Bill',        'House',  'Gryffendor',  14.7409,       'pd'],
            ['Charlie',     'House',  'Gryffendor',  15.2178,       'pd'],
            ['Percy',       'House',  'Gryffendor',  14.4404,       'pd'],
            ['Fred',        'House',  'Gryffendor',  25.6376,       'pd'],
            ['George',      'House',  'Gryffendor',  20.6805,       'pd'],
            ['Ron',         'House',  'Gryffendor',  19.0166,       'pd'],
            ['Ginny',       'House',  'Gryffendor',  28.4882,       'pd'],
            ['Harry',       'House',  'Gryffendor',   8.9534,       'pd'],
            ['Hermione',    'House',  'Gryffendor',  30.2340,       'pd'],
            ], dtype=object)
        long_cols = ['#SampleID', 'category', 'group', 'alpha', 'metric']
        self.long_ = pd.DataFrame(data=long_,
                                  columns=long_cols)
        self.long_.set_index('#SampleID', inplace=True)
        for column in ['category', 'group', 'alpha', 'metric']:
            self.long_[column] = self.long_[column].astype(str)
        self.long_['alpha'] = self.long_['alpha'].astype(float)
        self.long_.index.name = None

        # Sets up the modified longform table
        self.mod = 3
        self.long_mod = self.long_.copy()
        self.long_mod.loc[self.long_mod['metric'] == 'pd', 'alpha_mod'] = \
            (self.long_mod.loc[self.long_mod['metric'] == 'pd', 'alpha'] /
             self.mod)
        self.long_mod.loc[self.long_mod['metric'] == 'shannon', 'alpha_mod'] = \
            self.long_mod.loc[self.long_mod['metric'] == 'shannon', 'alpha']
        self.long_mod.loc[self.long_mod.metric == 'pd', 'metric'] = \
            'PD Whole Tree Diversity'
        self.long_mod.loc[self.long_mod.metric == 'shannon', 'metric'] = \
            'Shannon Diversity'

    def test_update_map(self):
        test = update_map(self.wide, self.alpha)
        pdt.assert_frame_equal(test, self.alpha_map)

    def test_horizontal_to_longtable(self):
        test = horizontal_to_longtable(sample=self.sample,
                                       map_=self.alpha_map,
                                       categories=self.categories,
                                       metrics=self.metrics)
        pdt.assert_frame_equal(test, self.long_)

    def test_horizontal_to_longtable_sample_error(self):
        with self.assertRaises(ValueError):
            horizontal_to_longtable(sample='Draco',
                                    map_=self.alpha_map,
                                    categories=self.categories,
                                    metrics=self.metrics)

    def test_horizontal_to_longtable_metric_error(self):
        with self.assertRaises(ValueError):
            horizontal_to_longtable(sample=self.sample,
                                    map_=self.alpha_map,
                                    categories=self.categories,
                                    metrics=['chao1_1000'])

    def test_horizontal_to_longtable_category_error(self):
        with self.assertRaises(ValueError):
            horizontal_to_longtable(sample=self.sample,
                                    map_=self.alpha_map,
                                    categories=['Mischief'],
                                    metrics=self.metrics)

    def test_modify_alpha_shannon_pd(self):
        test = modify_alpha_shannon_pd(self.long_, self.mod)
        pdt.assert_frame_equal(test, self.long_mod)

    def test_ag_summary_violin_fecal(self):
        # Defines the known values
        known_categories = ['BODY_HABITAT', 'AGE_CAT', 'ALCOHOL_FREQUENCY',
                            'ANTIBIOTIC_HISTORY', 'BMI_CAT',
                            'FLOSSING_FREQUENCY']
        known_metrics = ['shannon_1000', 'pd_whole_tree_1000']
        known_labels = ['All Fecal Samples', 'Simillar Age', 'Alcohol Use',
                        'Last Antibiotic Use', 'Simillar BMI',
                        'Similar Flossing \nFrequency']
        known_mod = 3
        known_ylim = [0, 12]
        known_shannon_yticks = np.arange(0, 13, 2)
        known_pd_yticks = np.arange(0, 39, 6)
        known_pdline = 19.0166/known_mod
        known_shannon = 5.2099
        known_pdvalue = 19.0166
        # Performs the test
        (categories, labels, metrics, mod, ylim, shannon_yticks, pd_yticks,
         pdline, shannon, pdvalue) = ag_summary_violin(
            'Ron', self.wide, self.alpha, debug=True
            )
        self.assertEqual(known_categories, categories)
        self.assertEqual(known_labels, labels)
        self.assertEqual(known_metrics, metrics)
        self.assertEqual(known_mod, mod)
        self.assertEqual(known_ylim, ylim)
        npt.assert_array_equal(known_shannon_yticks, shannon_yticks)
        npt.assert_array_equal(known_pd_yticks, pd_yticks)
        npt.assert_almost_equal(known_pdvalue, pdvalue, 3)
        npt.assert_almost_equal(known_pdline, pdline, 3)
        npt.assert_almost_equal(known_shannon, shannon, 3)

    def test_ag_summary_violin_skin(self):
        # Defines the known values
        known_categories = ['BODY_HABITAT', 'AGE_CAT', 'ALCOHOL_FREQUENCY',
                            'ANTIBIOTIC_HISTORY', 'BMI_CAT', 'DOMINANT_HAND']
        known_metrics = ['shannon_1000', 'pd_whole_tree_1000']
        known_labels = ['All Skin Samples', 'Simillar Age', 'Alcohol Use',
                        'Last Antibiotic Use', 'Simillar BMI',
                        'Same Dominant hand']
        known_mod = 4
        known_ylim = [0, 12.5]
        known_shannon_yticks = np.arange(0, 13, 2.5)
        known_pd_yticks = np.arange(0, 51, 10)
        known_pdline = 8.9534/known_mod
        known_shannon = 8.1074
        known_pdvalue = 8.9534
        # Performs the test
        (categories, labels, metrics, mod, ylim, shannon_yticks, pd_yticks,
         pdline, shannon, pdvalue) = ag_summary_violin(
            'Harry', self.wide, self.alpha, debug=True
            )
        self.assertEqual(known_categories, categories)
        self.assertEqual(known_labels, labels)
        self.assertEqual(known_metrics, metrics)
        self.assertEqual(known_mod, mod)
        self.assertEqual(known_ylim, ylim)
        npt.assert_array_equal(known_shannon_yticks, shannon_yticks)
        npt.assert_array_equal(known_pd_yticks, pd_yticks)
        npt.assert_almost_equal(known_pdvalue, pdvalue, 3)
        npt.assert_almost_equal(known_pdline, pdline, 3)
        npt.assert_almost_equal(known_shannon, shannon, 3)

    def test_ag_summary_violin_oral(self):
        # Defines the known values
        known_categories = ['BODY_HABITAT', 'AGE_CAT', 'ALCOHOL_FREQUENCY',
                            'ANTIBIOTIC_HISTORY', 'BMI_CAT',
                            'FLOSSING_FREQUENCY']
        known_metrics = ['shannon_1000', 'pd_whole_tree_1000']
        known_labels = ['All Oral Samples', 'Simillar Age', 'Alcohol Use',
                        'Last Antibiotic Use', 'Simillar BMI',
                        'Similar Flossing \nFrequency']
        known_mod = 2.5
        known_ylim = [0, 10]
        known_shannon_yticks = np.arange(0, 11, 2)
        known_pd_yticks = np.arange(0, 26, 5)
        known_pdline = 30.2340/known_mod
        known_shannon = 5.0945
        known_pdvalue = 30.2340
        # Performs the test
        (categories, labels, metrics, mod, ylim, shannon_yticks, pd_yticks,
         pdline, shannon, pdvalue) = ag_summary_violin(
            'Hermione', self.wide, self.alpha, debug=True
            )
        self.assertEqual(known_categories, categories)
        self.assertEqual(known_labels, labels)
        self.assertEqual(known_metrics, metrics)
        self.assertEqual(known_mod, mod)
        self.assertEqual(known_ylim, ylim)
        npt.assert_array_equal(known_shannon_yticks, shannon_yticks)
        npt.assert_array_equal(known_pd_yticks, pd_yticks)
        npt.assert_almost_equal(known_pdvalue, pdvalue, 3)
        npt.assert_almost_equal(known_pdline, pdline, 3)
        npt.assert_almost_equal(known_shannon, shannon, 3)

    def test_ag_summary_violin_nose(self):
        with self.assertRaises(ValueError):
            ag_summary_violin('Percy', self.wide, self.alpha, debug=True)

if __name__ == '__main__':
    main()
