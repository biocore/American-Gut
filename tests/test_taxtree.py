#!/usr/bin/env python

from unittest import TestCase, main

from numpy import array
from copy import deepcopy
from biom import Table

from americangut.taxtree import create_node, add_node, get_node, update_tree, \
        get_rare_unique, traverse, sample_rare_unique, \
        build_persample_tree_from_taxontable, set_relative_freqs, \
        update_per_sample_tree

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Daniel McDonald"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"

class TaxTreeTests(TestCase):
    def setUp(self):
        ta = {'name': 'root',
              'count': 36,
              'freq': 1.0,
              'children': [{'name': 'k__1',
                            'count': 36,
                            'freq': 1.0,
                            'children': [{'name': 'p__x',
                                          'count': 12,
                                          'freq': 12.0 / 36.0,
                                          'children': [{'name': 'c__1',
                                                        'count': 1,
                                                        'freq': 1.0 / 36,
                                                        'children': []},
                                                       {'name': 'c__2',
                                                        'count': 4,
                                                        'freq': 4.0 / 36,
                                                        'children': []},
                                                       {'name': 'c__',
                                                        'count': 7,
                                                        'freq': 7.0 / 36,
                                                        'children': []}]},
                                         {'name': 'p__y',
                                          'count': 24,
                                          'freq': 24.0 / 36.0,
                                          'children': [{'name': 'c__3',
                                                        'count': 10,
                                                        'freq': 10.0 / 36.0,
                                                        'children': []},
                                                       {'name': 'c__',
                                                        'count': 14,
                                                        'freq': 14.0 / 36.0,
                                                        'children': []}]}]}]}
        self.persample = ta

    def test_build_persample_tree_from_taxontable(self):
        exp_a = ('a', self.persample)
        obs_a = build_persample_tree_from_taxontable(table).next()
        self.assertEqual(obs_a, exp_a)

    def test_set_relative_freqs(self):
        persample_copy = deepcopy(self.persample)
        self.persample['freq'] = None
        for node in traverse(self.persample):
            node['freq'] = None
        set_relative_freqs(self.persample)
        self.assertEqual(persample_copy, self.persample)

    def test_update_per_sample_tree(self):
        taxons = (('k__a', 'p__x'),
                  ('k__a', 'p__x'),
                  ('k__a', 'p__y'),  # indep of p__x; c__1
                  ('k__b', 'p__z'))
        counts = (5, 7, 9, 11)
        exp = {'name': 'root',
               'count': 32,
               'freq': None,
               'children': [{'name': 'k__a',
                             'count': 21,
                             'freq': None,
                             'children': [{'name': 'p__x',
                                           'count': 12,
                                           'freq': None,
                                           'children': []},
                                          {'name': 'p__y',
                                           'count': 9,
                                           'freq': None,
                                           'children': []}]},
                            {'name': 'k__b',
                             'count': 11,
                             'freq': None,
                             'children': [{'name': 'p__z',
                                           'count': 11,
                                           'freq': None,
                                           'children': []}]}]}

        tree = create_node('root', count=0, freq=None)
        update_per_sample_tree(tree, taxons[0], counts[0])
        update_per_sample_tree(tree, taxons[1], counts[1])
        update_per_sample_tree(tree, taxons[2], counts[2])
        update_per_sample_tree(tree, taxons[3], counts[3])
        self.assertEqual(tree, exp)

    def test_sample_rare_unique(self):
        t = update_tree(None, tax_strings_by_sample)
        tax_by_sample = {'a':tax_strings_by_sample[0],
                         'b':tax_strings_by_sample[1],
                         'c':tax_strings_by_sample[2]}
        exp = [('a', None, [['k__1','p__x','c__'],['k__1','p__y','c__3']],
                           [['k__1','p__x','c__1'],['k__1','p__x','c__2']]),
               ('b', None, [['k__1','p__x','c__'],['k__1','p__y','c__3']], []),
               ('c', None, [], [])]
        obs = sample_rare_unique(t, None, tax_by_sample, 0.7)
        self.assertEqual(sorted(obs), exp)

        table_a = Table(array([[14,15,16]]),
                               ['k__1; p__y; c__'],
                               ['a','b','c'])
        table_b = Table(array([[1,2,3],
                               [4,5,6],
                               [14,15,16]]),
                        ['k__1; p__x; c__1',
                         'k__1; p__x; c__2',
                         'k__1; p__y; c__'],
                        ['a','b','c'], )
        table_c = Table(array([[1,2,3],
                               [4,5,6],
                               [7,8,9],
                               [10,11,12],
                               [14,15,16]]),
                        ['k__1; p__x; c__1',
                         'k__1; p__x; c__2',
                         'k__1; p__x; c__',
                         'k__1; p__y; c__3',
                         'k__1; p__y; c__'],
                        ['a','b','c'])

        exp = [('a', table_a, [['k__1','p__x','c__'],['k__1','p__y','c__3']],
                           [['k__1','p__x','c__1'],['k__1','p__x','c__2']]),
               ('b', table_b, [['k__1','p__x','c__'],['k__1','p__y','c__3']], []),
               ('c', table_c, [], [])]

        obs = sample_rare_unique(t, table, tax_by_sample, 0.7)
        for o,e in zip(sorted(obs), exp):
            self.assertEqual(o[0], e[0])
            self.assertEqual(o[1], e[1])
            self.assertEqual(o[2], e[2])
            self.assertEqual(o[3], e[3])

    def test_create_node(self):
        exp = {'name':'foo','children':[]}
        obs = create_node('foo')
        self.assertEqual(obs, exp)

    def test_traverse(self):
        t = update_tree(None, tax_strings_by_sample)
        exp = ['c__1','c__2','p__x','c__3','p__y','k__1','root']
        obs = [n['name'] for n in traverse(t)]
        self.assertEqual(obs, exp)

    def test_add_node(self):
        exp = {'name':'foo','popcount':1,'children':[
                        {'name':'bar', 'popcount':2, 'children':[]}
                     ]
              }
        foo = create_node('foo', popcount=1)
        bar = create_node('bar', popcount=2)
        add_node(foo,bar)
        self.assertEqual(foo, exp)

    def test_get_node(self):
        tree = {'name':'foo','popcount':0,'children':[
                        {'name':'bar', 'popcount':0, 'children':[]}
                    ]}
        self.assertEqual(get_node(tree, 'does not exist'), None)
        self.assertEqual(get_node(tree, 'bar'), tree['children'][0])

    def test_create_tree(self):
        """take tax strings, create a tree"""
        exp = {'name':'root',
                'popcount':3,
                'children':[
                    {'name':'k__1',
                     'popcount':3,
                     'children':[
                        {'name':'p__x',
                         'popcount':2,
                         'children':[
                            {'name':'c__1',
                             'popcount':1,
                             'children':[]},
                            {'name':'c__2',
                             'popcount':1,
                             'children':[]}
                            ]
                        },
                        {'name':'p__y',
                         'popcount':3,
                         'children':[
                             {'name':'c__3',
                              'popcount':2,
                              'children':[]}
                             ]
                         }
                    ]}
                ]}
        tax_strings = deepcopy(tax_strings_by_sample)
        tax_strings[0].append(['k__[1]','p__x','c__'])
        obs = update_tree(None, tax_strings)
        self.assertEqual(obs, exp)

    def test_create_tree_ignore_contested(self):
        """take tax strings, create a tree"""
        exp = {'name':'root',
                'popcount':3,
                'children':[
                    {'name':'k__1',
                     'popcount':3,
                     'children':[
                        {'name':'p__x',
                         'popcount':2,
                         'children':[
                            {'name':'c__1',
                             'popcount':1,
                             'children':[]},
                            {'name':'c__2',
                             'popcount':1,
                             'children':[]}
                            ]
                        },
                        {'name':'p__y',
                         'popcount':3,
                         'children':[
                             {'name':'c__3',
                              'popcount':2,
                              'children':[]}
                             ]
                         }
                    ]}
                ]}
        obs = update_tree(None, tax_strings_by_sample)
        self.assertEqual(obs, exp)
    def test_get_rare_unique(self):
        t = update_tree(None, tax_strings_by_sample)
        exp_unique = [['k__1','p__x','c__2']]
        exp_rare = [['k__1','p__x','c__'],['k__1','p__y','c__3']]
        obs_rare, obs_unique = get_rare_unique(t, sample_taxa, 0.7)

        self.assertEqual(obs_rare, exp_rare)
        self.assertEqual(obs_unique, exp_unique)

        exp_unique = [['k__1','p__x','c__2']]
        exp_rare = []
        obs_rare, obs_unique = get_rare_unique(t, sample_taxa, 0.4)

        self.assertEqual(obs_rare, exp_rare)
        self.assertEqual(obs_unique, exp_unique)

tax_strings_by_sample = [
        [
            ['k__1','p__x','c__1'],
            ['k__1','p__x','c__2'],
            ['k__1','p__x','c__'],
            ['k__1','p__y','c__3']
        ],
        [
            ['k__1','p__x','c__'],
            ['k__1','p__y','c__3']
        ],
        [
            ['k__1','p__y','c__']
        ]
    ]

table = Table(array([[1,2,3],
                     [4,5,6],
                     [7,8,9],
                     [10,11,12],
                     [14,15,16]]),
              ['k__1; p__x; c__1',
               'k__1; p__x; c__2',
               'k__1; p__x; c__',
               'k__1; p__y; c__3',
               'k__1; p__y; c__'],
              ['a','b','c'])

sample_taxa = [
        ['k__1','p__x','c__2'],
        ['k__1','p__x','c__'],
        ['k__1','p__y','c__3']
    ]


if __name__ == '__main__':
    main()
