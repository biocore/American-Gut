from __future__ import division
from unittest import TestCase, main

import numpy as np
import numpy.testing as npt

from americangut.power_plots import collate_effect_size, summarize_effect

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
        self.counts = [np.array([5, 10, 15, 20]),
                       np.array([5, 10, 15, 20, 25, 30, 35, 40, 45])]
        self.power1 = [np.array([0.202, 0.434, 0.684, 0.954]),
                       np.array([0.146, 0.300, 0.442, 0.658, 0.824, 0.922,
                                 0.972, 1.000, 1.000])]
        self.power2 = [np.array([[0.202, 0.434, 0.684, 0.954],
                                 [0.200, 0.460, 0.688, 0.950],
                                 [0.214, 0.440, 0.684, 0.956]]),
                       np.array([[0.146, 0.300, 0.442, 0.658, 0.824, 0.922,
                                  0.972, 1.000, 1.000],
                                 [0.154, 0.310, 0.490, 0.628, 0.770, 0.918,
                                  0.976, 1.000, 1.000],
                                 [0.154, 0.300, 0.466, 0.630, 0.788, 0.916,
                                  0.984, 0.996, 1.000]])]
        self.effects = np.array([0.72802675, 0.5891959])
        self.bounds = np.array([0.05568197,  0.03144488])
        self.labels = np.array(['25_samples', '50_samples'])
        self.order = np.array([0, 1])
        tab = ['<table style="border-style:hidden;\n'
               '              border-collapse:collapse\n'
               '              line-height:120%\n'
               '              ">\n'
               '\t<tr>\n'
               '\t\t<th style="text-align:center;\n'
               '\t\t           background-color:black;\n'
               '\t\t           color:white\n'
               '\t\t           ">\n'
               '\t\t\tCategory\t\t</th>\n'
               '\t\t<th style="text-align:center;\n'
               '\t\t           background-color:black;\n'
               '\t\t           color:white";\n'
               '\t\t    colspan=3>\n'
               '\t\t\tAlpha\n'
               '\t\t</th>\n'
               '\t\t<td style="border-hidden;\n'
               '\t\t           background-color:black;\n'
               '\t\t           padding:20px">\n'
               '\t\t<th style="text-align:center;\n'
               '\t\t           background-color:black;\n'
               '\t\t         color:white";\n'
               '\t\t  colspan=3>\n'
               '\t\t\tBeta\n'
               '\t\t</th>\n'
               '\t</tr>\n'
               '\t<tr>\n'
               '\t\t<td style="border-top:hidden;\n'
               '\t\t           border-bottom:hidden;\n'
               '\t\t           border-left: hidden;\n'
               '\t\t           border-bottom: hidden;\n'
               '\t\t           padding:10px;\n'
               '\t\t           text-align:left\n'
               '\t\t          ">\n'
               '\t\t\t2\n'
               '\t\t</td>\n'
               '\t\t<td style="border-top:hidden;\n'
               '\t\t           border-bottom:hidden;\n'
               '\t\t           border-left: hidden;\n'
               '\t\t           border-bottom: hidden;\n'
               '\t\t           text-align:right\n'
               '\t\t          ">\n'
               '\t\t\t20\n'
               '\t\t</td>\n'
               '\t\t<td style="border-top:hidden;\n'
               '\t\t           border-bottom:hidden;\n'
               '\t\t           border-left: hidden;\n'
               '\t\t           border-bottom: hidden;\n'
               '\t\t           text-align:center\n'
               '\t\t          ">\n'
               '\t\t\t&plusmn;\n'
               '\t\t</td>\n'
               '\t\t<td style="border-top:hidden;\n'
               '\t\t           border-bottom:hidden;\n'
               '\t\t           border-left: hidden;\n'
               '\t\t           border-bottom: hidden;\n'
               '\t\t           text-align:right\n'
               '\t\t          ">\n'
               '\t\t\t5\n'
               '\t\t</td>\n'
               '\t\t<td style="border-top:hidden;\n'
               '\t\t           border-bottom:hidden;\n'
               '\t\t           border-left: hidden;\n'
               '\t\t           border-bottom: hidden;\n'
               '\t\t           padding:20px">\n'
               '\t\t</td>\n'
               '\t\t<td style="border-top:hidden;\n'
               '\t\t           border-bottom:hidden;\n'
               '\t\t           border-left: hidden;\n'
               '\t\t           border-bottom: hidden;\n'
               '\t\t           text-align:right\n'
               '\t\t          ">\n'
               '\t\t\t20\n'
               '\t\t</td>\n'
               '\t\t<td style="border-top:hidden;\n'
               '\t\t           border-bottom:hidden;\n'
               '\t\t           border-left: hidden;\n'
               '\t\t           border-bottom: hidden;\n'
               '\t\t           text-align:center\n'
               '\t\t          ">\n'
               '\t\t\t&plusmn;\n'
               '\t\t</td>\n'
               '\t\t<td style="border-top:hidden;\n'
               '\t\t           border-bottom:hidden;\n'
               '\t\t           border-left: hidden;\n'
               '\t\t           border-bottom: hidden;\n'
               '\t\t           text-align:right\n'
               '\t\t          ">\n'
               '\t\t\t5\n'
               '\t\t</td>\n'
               '\t<tr>\n'
               '\t\t<td style="border-top:hidden;\n'
               '\t\t           border-bottom:hidden;\n'
               '\t\t           border-left: hidden;\n'
               '\t\t           border-bottom: hidden;\n'
               '\t\t           padding:10px;\n'
               '\t\t           text-align:left\n'
               '\t\t          ">\n'
               '\t\t\t5\n'
               '\t\t</td>\n'
               '\t\t<td style="border-top:hidden;\n'
               '\t\t           border-bottom:hidden;\n'
               '\t\t           border-left: hidden;\n'
               '\t\t           border-bottom: hidden;\n'
               '\t\t           text-align:right\n'
               '\t\t          ">\n'
               '\t\t\t25\n'
               '\t\t</td>\n'
               '\t\t<td style="border-top:hidden;\n'
               '\t\t           border-bottom:hidden;\n'
               '\t\t           border-left: hidden;\n'
               '\t\t           border-bottom: hidden;\n'
               '\t\t           text-align:center\n'
               '\t\t          ">\n'
               '\t\t\t&plusmn;\n'
               '\t\t</td>\n'
               '\t\t<td style="border-top:hidden;\n'
               '\t\t           border-bottom:hidden;\n'
               '\t\t           border-left: hidden;\n'
               '\t\t           border-bottom: hidden;\n'
               '\t\t           text-align:right\n'
               '\t\t          ">\n'
               '\t\t\t5\n'
               '\t\t</td>\n'
               '\t\t<td style="border-top:hidden;\n'
               '\t\t           border-bottom:hidden;\n'
               '\t\t           border-left: hidden;\n'
               '\t\t           border-bottom: hidden;\n'
               '\t\t           padding:20px">\n'
               '\t\t</td>\n'
               '\t\t<td style="border-top:hidden;\n'
               '\t\t           border-bottom:hidden;\n'
               '\t\t           border-left: hidden;\n'
               '\t\t           border-bottom: hidden;\n'
               '\t\t           text-align:right\n'
               '\t\t          ">\n'
               '\t\t\t25\n'
               '\t\t</td>\n'
               '\t\t<td style="border-top:hidden;\n'
               '\t\t           border-bottom:hidden;\n'
               '\t\t           border-left: hidden;\n'
               '\t\t           border-bottom: hidden;\n'
               '\t\t           text-align:center\n'
               '\t\t          ">\n'
               '\t\t\t&plusmn;\n'
               '\t\t</td>\n'
               '\t\t<td style="border-top:hidden;\n'
               '\t\t           border-bottom:hidden;\n'
               '\t\t           border-left: hidden;\n'
               '\t\t           border-bottom: hidden;\n'
               '\t\t           text-align:right\n'
               '\t\t          ">\n'
               '\t\t\t5\n'
               '\t\t</td>\n'
               '</table>'
               ]
        self.table = tab[0]

    def test_collate_effect_size_unequal_counts_error(self):
        with self.assertRaises(ValueError):
            collate_effect_size(counts=[np.array([1, 2, 3])],
                                powers=self.power2,
                                alpha=self.alpha)

    def test_collate_effect_size_count_shape_error(self):
        with self.assertRaises(TypeError):
            collate_effect_size(counts=self.power2,
                                powers=self.power2,
                                alpha=self.alpha)

    def test_collate_effect_size_count_for_power_1d(self):
        with self.assertRaises(ValueError):
            collate_effect_size(counts=[self.counts[0]],
                                powers=self.power1[1],
                                alpha=0.05)

    def test_collate_effect_size_count_for_power_2d(self):
        with self.assertRaises(ValueError):
            collate_effect_size(counts=[self.counts[0]],
                                powers=self.power2[1],
                                alpha=0.05)

    def test_collate_effect_size_1_power(self):
        known_means = np.array([0.72374193, 0.5795861])
        known_bounds = np.array([0.17654622, 0.06264818])
        test_means, test_bounds = collate_effect_size(counts=self.counts,
                                                      powers=self.power1,
                                                      alpha=self.alpha)
        npt.assert_almost_equal(known_means, test_means, 4)
        npt.assert_almost_equal(known_bounds, test_bounds, 4)

    def test_collate_effect_size_2_power(self):
        # The known effect is used as a reference elsewhere
        test_means, test_bounds = collate_effect_size(counts=self.counts,
                                                      powers=self.power2,
                                                      alpha=self.alpha)
        npt.assert_almost_equal(self.effects, test_means, 4)
        npt.assert_almost_equal(self.bounds, test_bounds, 4)

    def test_collate_effect_size_nan(self):
        power = [count*0 for count in self.counts]

        test_mean, test_bounds = collate_effect_size(counts=self.counts,
                                                     powers=power,
                                                     alpha=self.alpha)

        self.assertTrue(np.isnan(test_mean[0]))
        self.assertTrue(np.isnan(test_bounds[0]))

    def test_summarize_effect(self):
        test_table = summarize_effect(order=self.order,
                                      fecal_cats=self.labels,
                                      a_eff_means=self.effects,
                                      a_eff_bounds=self.bounds,
                                      b_eff_means=self.effects,
                                      b_eff_bounds=self.bounds)
        self.assertEqual(self.table, test_table)

if __name__ == '__main__':
    main()
