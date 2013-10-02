#!/usr/bin/env python

from unittest import TestCase, main
from americangut.util import pick_rarifaction_level

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Daniel McDonald"]
__license__ = "GPLv2"
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

if __name__ == '__main__':
    main()
