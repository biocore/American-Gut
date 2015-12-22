from unittest import TestCase, main

import numpy as np
import matplotlib
import pandas as pd

import numpy.testing as npt
import pandas.util.testing as pdt

from americangut.ag_alpha_plots import plot_alpha_distribution


class AlphaPlotTest(TestCase):
    def setUp(self):
        self.map_ = pd.DataFrame(
            data=np.array([
                ['UBERON:skin', '1990', 'female', 'Verity', 'US', 12.5],
                ['UBERON:feces', '1990', 'female', 'Verity', 'US', 8.6],
                ['UBERON:feces', '1987', 'male', 'Alex', 'US', 7.9],
                ['UBERON:feces', '1993', 'female', 'Annie', 'US', 7.5],
                ['UBERON:skin', '1989', 'male', 'Dominic', 'UK', 14.0],
                ['UBERON:feces', '1986', 'female', 'Sarah', 'US', 15.0],
                ['UBERON:oral cavity', '1988', 'female', 'Shelby', 'AUS', 4.2],
                ]),
            index=['VeP0', 'VeP1', 'AxP0', 'AnP0', 'DoD0', 'SaZ0', 'ShT0'],
            columns=['BODY_HABITAT', 'BIRTH_YEAR', 'SEX', 'HOST_SUBJECT_ID',
                     'NATIONALITY', 'alpha'],
            )

    def test_plot_alpha_distribution(self):
        kgroup = 'UBERON:skin'
        sample = 'VeP0'
        kgalpha = pd.Series([12.5, 14.0],
                            index=['VeP0', 'DoD0'],
                            name='alpha',
                            dtype=object)
        ksalpha = 12.5
        kxlabel = 'alpha'

        (tgroup, tgalpha, tsalpha, tax, txlabel) = \
            plot_alpha_distribution(sample,
                                    self.map_,
                                    debug=True)
            

        self.assertEqual(tgroup, kgroup)
        npt.assert_array_equal(tgalpha.values, kgalpha.values)
        self.assertEqual(tsalpha, ksalpha)
        self.assertTrue(isinstance(tax, matplotlib.axes.Axes))
        self.assertEqual(txlabel, kxlabel)

if __name__ == '__main__':
    main()
