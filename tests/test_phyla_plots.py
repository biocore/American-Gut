#!/usr/bin/env python

from numpy import array
from unittest import TestCase, main
from biom.table import table_factory
from americangut.phyla_plots import (table_by_category_value,
        category_value_lookup, drop_sample, sort_by_taxa,
        average_per_observation, make_collapse_f, most_common_taxa)


__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Daniel McDonald", "Justine Debelius"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"


class PhylaPlotsTests(TestCase):
    def setUp(self):
        self.table = table
        self.table.addSampleMetadata(metadata)

    def test_table_by_category_value_normalized(self):
        """Get tables by category value"""
        test_t = self.table.copy()
        exp = [test_t.filterSamples(lambda v, i, md: i in ['s1', 's3']),
               test_t.filterSamples(lambda v, i, md: i in ['s4']),
               test_t.filterSamples(lambda v, i, md: i in ['s2'])]
        obs, order = table_by_category_value(self.table, 'cat1')
        for o, e in zip(obs, exp):
            self.assertEqual(o, e)
        self.assertEqual(order, ['10', '15', '20'])

    def test_category_value_lookup(self):
        """Construct a lookup for all tables of interest"""
        def make_filter_f(ids):
            def filter_f(v, i, md):
                return i in ids
            return filter_f
        norm = lambda x: x

        t = self.table.copy()
        obs = category_value_lookup(self.table, ['cat1', 'cat2', 'cat3'])
        exp = {'cat1': {
                   '10': norm(t.filterSamples(make_filter_f(['s1', 's3']))),
                   '15': norm(t.filterSamples(make_filter_f(['s4']))),
                   '20': norm(t.filterSamples(make_filter_f(['s2'])))},
               'cat2': {
                    'a': norm(t.filterSamples(make_filter_f(['s1', 's2']))),
                    'b': norm(t.filterSamples(make_filter_f(['s3']))),
                    'c': norm(t.filterSamples(make_filter_f(['s4'])))},
               'cat3': {
                  'xyz': norm(t.filterSamples(make_filter_f(['s1']))),
              'no_data': norm(t.filterSamples(make_filter_f(['s2']))),
                  '123': norm(t.filterSamples(make_filter_f(['s3']))),
                   None: norm(t.filterSamples(make_filter_f(['s4'])))}
               }
        self.assertEqual(obs, exp)

    def test_drop_sample(self):
        """Remove a sample from a table"""
        f = lambda v, i, md: i in ['s1', 's3', 's4']
        exp = self.table.copy().filterSamples(f)
        obs = drop_sample(self.table, 's2')
        self.assertEqual(obs, exp)

    def test_sort_by_taxa(self):
        """Sort a vector by taxa"""
        desired_order = array(['x', 'y', 'z', '1', '2', '3'])
        vector = array([0.1, 0.9, 0.2, 0.5, 0.3, 0.3])
        vector_order = array(['1', 'x', '2', 'y', '3', 'z'])
        exp = array([0.9, 0.5, 0.3, 0.1, 0.2, 0.3])
        obs = sort_by_taxa(vector, vector_order, desired_order)
        self.assertTrue((obs == exp).all())

    def test_average_per_observation(self):
        """Get the average per observation in a table"""
        exp = array([0.5 + 0.0 + (1/3.) + (1/3.),
                     0.5 + 0.0 + 0.0 + (1/3.),
                     0.0 + 1.0 + (1/3.) + 0.0,
                     0.0 + 0.0 + (1/3.) + (1/3.)]) / 4.0
        obs = average_per_observation(table.normObservationBySample())
        self.assertTrue((obs == exp).all())
        self.assertAlmostEqual(obs.sum(), 1.0)

    def test_make_collapse_f(self):
        """Make a collapsing function"""
        f = make_collapse_f(['a', 'b', 'c'])
        obs = [f({'id': c}) for c in ['a', 'z', 'c', 'd', 'b', 'b', 'q']]
        exp = ['a', 'other', 'c', 'other', 'b', 'b', 'other']
        self.assertEqual(obs, exp)

    def test_most_common_taxa(self):
        """Determine the most common taxa"""
        obs = most_common_taxa(self.table.normObservationBySample(), 2)
        exp = ['k__Bacteria;p__what;c__cool', 'k__Bacteria;p__foo;c__bar']
        self.assertEqual(obs, exp)


table = table_factory(array([[1, 0, 1, 1],
                             [1, 0, 0, 1],
                             [0, 1, 1, 0],
                             [0, 0, 1, 1]]),
                      ['s1', 's2', 's3', 's4'],
                      ['k__Bacteria;p__foo;c__bar',
                       'k__Bacteria;p__foo;c__other',
                       'k__Bacteria;p__what;c__cool',
                       'k__Bacteria;p__what;c__boom'])


metadata = {'s1': {'cat1': '10', 'cat2': 'a', 'cat3': 'xyz'},
            's2': {'cat1': '20', 'cat2': 'a', 'cat3': 'no_data'},
            's3': {'cat1': '10', 'cat2': 'b', 'cat3': '123'},
            's4': {'cat1': '15', 'cat2': 'c', 'cat3': None}}


if __name__ == '__main__':
    main()
