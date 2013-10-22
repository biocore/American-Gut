#!/usr/bin/env python

from unittest import TestCase, main
from americangut.util import pick_rarifaction_level, slice_mapping_file, \
        parse_mapping_file, verify_subset
from biom.table import table_factory
from numpy import array
from StringIO import StringIO

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Daniel McDonald"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"

class UtilTests(TestCase):
    def test_pick_rarifaction_level(self):
        ids_10k = {'a':'a.1', '000001000':'000001000.123'}
        ids_1k = {'a':'a.1', '000001000':'000001000.123', 'b':123}
        
        exp_a = '10k'
        exp_b = '1k'
        exp_c = None

        obs_a = pick_rarifaction_level('a', [('10k',ids_10k), ('1k',ids_1k)])
        obs_b = pick_rarifaction_level('b', [('10k',ids_10k), ('1k',ids_1k)])
        obs_c = pick_rarifaction_level('c', [('10k',ids_10k), ('1k',ids_1k)])
        
        self.assertEqual(obs_a, exp_a)
        self.assertEqual(obs_b, exp_b)
        self.assertEqual(obs_c, exp_c)

    def test_verify_subset(self):
        metadata = [('a','other stuff\tfoo'), ('b', 'asdasdasd'), 
                    ('c','123123123')]
        table = table_factory(array([[1,2,3],[4,5,6]]), ['a','b','c'], ['x','y'])
        self.assertTrue(verify_subset(table, metadata))
        table = table_factory(array([[1,2],[3,4]]), ['a','b'], ['x','y'])
        self.assertTrue(verify_subset(table, metadata))
        table = table_factory(array([[1,2,3],[4,5,6]]), ['a','b','x'], ['x','y'])
        self.assertFalse(verify_subset(table, metadata))

    def test_slice_mapping_file(self):
        header, metadata = parse_mapping_file(StringIO(test_mapping))
        table = table_factory(array([[1,2],[4,5]]), ['a','c'], ['x','y'])
        exp = ["a\t1\t123123", "c\tpoop\tdoesn't matter"]
        obs = slice_mapping_file(table, metadata)
        self.assertEqual(obs,exp)

    def test_parse_mapping_file(self):
        exp = ("#SampleIDs\tfoo\tbar", [['a','1\t123123'],
                                        ['b','yy\txxx'],
                                        ['c',"poop\tdoesn't matter"]])
        obs = parse_mapping_file(StringIO(test_mapping))
        self.assertEqual(obs, exp)

test_mapping = """#SampleIDs\tfoo\tbar
a\t1\t123123
b\tyy\txxx
c\tpoop\tdoesn't matter
"""

if __name__ == '__main__':
    main()
