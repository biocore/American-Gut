#!/usr/bin/env python

from __future__ import division

__author__ = "Sam Way"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Sam Way"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Sam Way"
__email__ = "samuel.way@colorado.edu"

from os.path import realpath, dirname

from numpy import array, array_equal
from unittest import TestCase, main

from americangut.parse import parse_mapping_file_to_dict, \
    get_filtered_taxa_summary

TEST_DIR = dirname(realpath(__file__))
TEST_MAPPING_FILE = TEST_DIR + '/files/test_mapping.txt'
TEST_TAXA_FILE = TEST_DIR + '/files/test_taxa.txt'


class TestMappingFileParse(TestCase):
    def setUp(self):
        self.categories = ['TEST_CATEGORY', 'AWESOME_CATEGORY']
        self.sample_ids = ['sample_a', 'sample_b']
        self.metadata_dict = {'sample_a':
                             {'TEST_CATEGORY': '1',
                              'AWESOME_CATEGORY': 'super'},
                              'sample_b':
                             {'TEST_CATEGORY': '2',
                              'AWESOME_CATEGORY': 'totally'}}

    def test_mapping_file(self):
        mapping_dict, comments = \
            parse_mapping_file_to_dict(open(TEST_MAPPING_FILE, 'U'))
        for sample_id, sample_dict in mapping_dict.iteritems():
            # Does the sample dictionary contain
            # all of the metadata categories?
            self.assertEqual(set(sample_dict), set(self.categories))

            # Are all metadata values correct?
            for category in sample_dict.keys():
                self.assertEqual(sample_dict[category],
                                 self.metadata_dict[sample_id][category])


class TestTaxaSummaryFileParse(TestCase):
    def setUp(self):
        self.sample_ids = ['sample_a', 'sample_b']
        self.table = array([[0.11, 0.15],
                            [0.12, 0.14],
                            [0.13, 0.13],
                            [0.14, 0.12],
                            [0.15, 0.11],
                            [0.16, 0.10],
                            [0.09, 0.09],
                            [0.10, 0.16]])
        self.taxa_ids = ['Firmicutes', 'Bacteroidetes',
                         'Proteobacteria', 'Verrucomicrobia',
                         'Actinobacteria', 'Tenericutes',
                         'Cyanobacteria']
        self.taxa_labels = self.taxa_ids + ['Other']
        self.metadata_category = 'TEST_CATEGORY'
        self.metadata_value = '1'

    def test_taxa_file(self):
        filtered_sample_ids, taxa_labels, collapsed_taxa_table = \
            get_filtered_taxa_summary(TEST_MAPPING_FILE, TEST_TAXA_FILE,
                                      self.metadata_category,
                                      self.metadata_value,
                                      select_taxa=self.taxa_ids)

        # Make sure we only get out the matching sample
        self.assertEqual(len(filtered_sample_ids), 1)
        self.assertEqual(filtered_sample_ids[0], 'sample_a')

        # Should get back our desired labels + "Others"
        self.assertEqual(set(self.taxa_labels), set(taxa_labels))

        # Should get the slide of the table corresponding
        # to the matching sample
        self.assertTrue(array_equal(collapsed_taxa_table,
                        self.table[:, 0, None]))

if __name__ == '__main__':
    main()
