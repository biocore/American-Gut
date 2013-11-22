#!/usr/bin/env python

from unittest import TestCase, main
from americangut.taxtree import create_node, add_node, get_node, update_tree, \
        get_rare_unique, traverse, sample_rare_unique
from biom.table import table_factory
from numpy import array

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Daniel McDonald"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"

class TaxTreeTests(TestCase):
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

        table_a = table_factory(array([[14,15,16]]), ['a','b','c'], 
                                ['k__1; p__y; c__'])
        table_b = table_factory(array([[1,2,3],
                                       [4,5,6],
                                       [14,15,16]]), ['a','b','c'], 
                                    ['k__1; p__x; c__1',
                                     'k__1; p__x; c__2',
                                     'k__1; p__y; c__'])
        table_c = table_factory(array([[1,2,3],
                                       [4,5,6],
                                       [7,8,9],
                                       [10,11,12],
                                       [14,15,16]]), ['a','b','c'], 
                                    ['k__1; p__x; c__1',
                                     'k__1; p__x; c__2',
                                     'k__1; p__x; c__',
                                     'k__1; p__y; c__3',
                                     'k__1; p__y; c__'])

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
        exp = {'name':'foo','popcount':0,'children':[]}
        obs = create_node('foo')
        self.assertEqual(obs, exp)

    def test_traverse(self):
        t = update_tree(None, tax_strings_by_sample)
        exp = ['c__1','c__2','p__x','c__3','p__y','k__1','root']
        obs = [n['name'] for n in traverse(t)]
        self.assertEqual(obs, exp)

    def test_add_node(self):
        exp = {'name':'foo','popcount':0,'children':[
                        {'name':'bar', 'popcount':0, 'children':[]}
                     ]
              }
        foo = create_node('foo')
        bar = create_node('bar')
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

table = table_factory(array([[1,2,3],
           [4,5,6],
           [7,8,9],
           [10,11,12],
           [14,15,16]]), ['a','b','c'],['k__1; p__x; c__1', 
                                        'k__1; p__x; c__2', 
                                        'k__1; p__x; c__',
                                        'k__1; p__y; c__3',
                                        'k__1; p__y; c__'])

sample_taxa = [
        ['k__1','p__x','c__2'],
        ['k__1','p__x','c__'],
        ['k__1','p__y','c__3']
    ]
if __name__ == '__main__':
    main()
