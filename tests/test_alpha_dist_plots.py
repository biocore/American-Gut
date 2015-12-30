from unittest import TestCase, main

import numpy as np
import pandas as pd
import matplotlib

import numpy.testing as npt

from americangut.alpha_dist_plots import plot_alpha


class AlphaPlotTest(TestCase):
    def setUp(self):
        self.map_ = pd.DataFrame(
            data=np.array([
                ['skin', '1990', 'female', 'Verity', 'US', 12.5],
                ['fecal', '1990', 'female', 'Verity', 'US', 8.6],
                ['fecal', '1987', 'male', 'Alex', 'US', 7.9],
                ['fecal', '1993', 'female', 'Annie', 'US', 7.5],
                ['skin', '1989', 'male', 'Dominic', 'UK', 14.0],
                ['fecal', '1986', 'female', 'Sarah', 'US', 15.0],
                ['oral', '1988', 'female', 'Shelby', 'AUS', 4.2],
                ]),
            index=['VeP0', 'VeP1', 'AxP0', 'AnP0', 'DoD0', 'SaZ0', 'ShT0'],
            columns=['SIMPLE_BODY_SITE', 'BIRTH_YEAR', 'SEX',
                     'HOST_SUBJECT_ID', 'NATIONALITY', 'alpha'],
            )

    def test_plot_alpha(self):
        sample = 'VeP0'
        kgroup = 'skin'
        kgalpha = pd.Series([12.5, 14.0],
                            index=['VeP0', 'DoD0'],
                            name='alpha',
                            dtype=object)
        ksalpha = 12.5
        kxlabel = 'alphadiversity'

        (tgroup, tgalpha, tsalpha, txlabel) = \
            plot_alpha(sample=sample,
                       alpha_map=self.map_,
                       alpha_field='alpha',
                       debug=True)

        self.assertEqual(tgroup, kgroup)
        npt.assert_array_equal(tgalpha.values, kgalpha.values)
        self.assertEqual(tsalpha, ksalpha)
        self.assertEqual(txlabel, kxlabel)

    def test_plot_alpha_no_sample(self):
        sample = 'AlP0'
        kvalue = 'alpha does not have an alpha diversity value for AlP0.'
        tvalue = plot_alpha(sample=sample,
                            alpha_map=self.map_,
                            alpha_field='alpha')
        self.assertEqual(tvalue, kvalue)

    def test_plot_alpha_alpha_field_error(self):
        with self.assertRaises(ValueError):
            plot_alpha(sample='VeP0',
                       alpha_map=self.map_,
                       alpha_field='InCryptid')

    def test_plot_alpha_group_field_error(self):
        with self.assertRaises(ValueError):
            plot_alpha(sample='VeP0',
                       alpha_map=self.map_,
                       alpha_field='alpha',
                       group_field='BODY_HABITAT')

if __name__ == '__main__':
    main()
