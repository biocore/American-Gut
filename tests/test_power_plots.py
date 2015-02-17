from __future__ import division
from unittest import TestCase, main

import numpy as np

from americangut.power_plots import collate_effect_size, plot_effects

__author__ = "Justine Debelius"
__copyright__ = "Copyright 2014"
__credits__ = ["Justine Debelius"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Justine Debelius"
__email__ = "Justine.Debelius@colorado.edu"


class PowerPlotTest(TestCase):
    def setUp(self):
        self.alpha = 0.05
        self.counts = []
        self.powers = []

    def test_collate_effect_size(self):
        pass

    def test_collate_effect_size_unequal_counts(self):
        
        with self.assertRaises(ValueError):
            collate_effect_size(self.counts=[np.array([1, 2, 3])])

    def test_collate_effect_size_count_shape(self):
        pass

    def test_collate_effect_size_count_for_power_1d(self):
        pass

    def test_collate_effect_size_count_for_power_2d(self):
        pass

    def test_collate_effect_size_1_power(self):
        pass

    def test_collate_effect_size_2_power(self):
        pass


if __name__ == '__main__':
    main()
