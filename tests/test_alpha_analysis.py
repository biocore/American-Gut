#!/usr/bin/env python

from __future__ import division
from unittest import TestCase, main
from os import rmdir
from os.path import realpath, dirname, join as pjoin, exists
from numpy import nan
from pandas import Series, DataFrame, Index
from pandas.util.testing import assert_index_equal
from matplotlib import use
use('agg', warn=False)
from americangut.alpha_analysis import pad_index, check_dir


__author__ = "Justine Debelius"
__copyright__ = "Copyright 2014"
__credits__ = ["Justine Debelius"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Justine Debelius"
__email__ = "Justine.Debelius@colorado.edu"

# Determines the location fo the reference files
TEST_DIR = dirname(realpath(__file__))


class AlphaAnalysisTest(TestCase):

    def setUp(self):
        """Initializes data for each test instance"""
        # Sets up lists for the data frame
        self.ids = ['000011224.1210932', '000005034.1076350',
                    '000011351.1210669', '000002013.1131550',
                    '000010163.1213006', '000007749.1130433',
                    '000001221.1075847', '000004758.1076374',
                    '000001111.1076323', '000012418.1210549',
                    '000013510.1257164', '000001571.1212966',
                    '000012994.1257089', '000004086.1213010',
                    '000013116.1257202', '000001793.1076039',
                    '000007044.1130122', '000002290.1210781',
                    '000002772.1130083', '000009628.1130042',
                    '000010900.1212969', '000003369.1076434',
                    '000007748.1130523', '000006080.1131591',
                    '000013982.1257320', '000003969.1075720',
                    '000011225.1210593', '000011246.1210942',
                    '000003430.1075947', '000003730.1130454',
                    '000009445.1130058', '000004659.1210528',
                    '000002030.1053854', '000001546.1131928',
                    '000009419.1130338', '000008951.1076310',
                    '000012265.1210552', '000001150.1076054',
                    '000005576.1076167', '000011119.1210964',
                    '000001711.1075845', '000004145.1076463',
                    '000011335.1210545', '000001613.1213143',
                    '000001549.1212975', '000012157.1210674',
                    '000011977.1257387', '000007037.1131746',
                    '000002276.1076156', '000010940.1210797',
                    '000003895.1076443', '000007668.1130386',
                    '000009670.1213024', '000004705.1130352',
                    '000010827.1210606', '000007005.1131779',
                    '000007151.1130047', '000005927.1131701',
                    '000009074.1213041', '000005033.1076435']
        self.raw_ids = ['11224.1210932', '5034.1076350', '11351.1210669',
                        '2013.1131550', '10163.1213006', '7749.1130433',
                        '1221.1075847', '4758.1076374', '1111.1076323',
                        '12418.1210549', '13510.1257164', '1571.1212966',
                        '12994.1257089', '4086.1213010', '13116.1257202',
                        '1793.1076039', '7044.1130122', '2290.1210781',
                        '2772.1130083', '9628.1130042', '10900.1212969',
                        '3369.1076434', '7748.1130523', '6080.1131591',
                        '13982.1257320', '3969.1075720', '11225.1210593',
                        '11246.1210942', '3430.1075947', '3730.1130454',
                        '9445.1130058', '4659.1210528', '2030.1053854',
                        '1546.1131928', '9419.1130338', '8951.1076310',
                        '12265.1210552', '1150.1076054', '5576.1076167',
                        '11119.1210964', '1711.1075845', '4145.1076463',
                        '11335.1210545', '1613.1213143', '1549.1212975',
                        '12157.1210674', '11977.1257387', '7037.1131746',
                        '2276.1076156', '10940.1210797', '3895.1076443',
                        '7668.1130386', '9670.1213024', '4705.1130352',
                        '10827.1210606', '7005.1131779', '7151.1130047',
                        '5927.1131701', '9074.1213041', '5033.1076435']
        self.age = ['50s', '30s', '50s', '30s', '60+', '20s', '50s', '40s',
                    '50s', '40s', '60+', '40s', '60+', '40s', '50s', '40s',
                    '50s', '40s', '50s', '50s', '30s', '60+', '30s', '50s',
                    '30s', '40s', '50s', '50s', '30s', '40s', '30s', '60+',
                    '30s', '20s', '30s', '20s', '60+', '60+', '30s', '50s',
                    '60+', '30s', '60+', '20s', '60+', '60+', '50s', '20s',
                    '40s', '50s', '40s', '40s', '60+', '30s', '50s', '20s',
                    '30s', '20s', '20s', '30s']
        self.gender = ['female', 'female', 'male', 'male', 'male', 'female',
                       'female', 'female', 'male', 'male', 'male', nan,
                       'male', 'female', nan, 'male', 'female', 'male',
                       'male', 'female', 'male', 'male', 'female', 'male',
                       nan, 'male', 'male', 'female', 'female', 'female',
                       'female', 'female', 'female', 'female', 'female',
                       'male', 'female', 'male', 'male', 'female', 'male',
                       'male', 'male', 'male', 'male', 'female', 'female',
                       'male', 'female', 'male', 'female', 'female', 'female',
                       'female', 'male', 'female', nan, 'male', 'male',
                       'female']
        self.season = ['Winter', 'Spring', 'Winter', 'Summer', 'Fall',
                       'Fall', 'Spring', 'Spring', 'Spring', 'Winter',
                       'Winter', 'Fall', 'Winter', 'Fall', 'Winter',
                       'Spring', 'Summer', 'Winter', 'Summer', 'Fall',
                       'Fall', 'Summer', 'Fall', 'Summer', 'Spring',
                       'Spring', 'Winter', 'Winter', 'Spring',  'Summer',
                       'Fall', 'Fall', 'Summer', 'Summer', 'Fall',
                       'Summer', 'Winter', 'Spring', 'Spring', 'Winter',
                       'Spring', 'Spring', 'Winter', 'Fall', 'Spring',
                       'Winter', 'Winter', 'Summer', 'Spring', 'Fall',
                       'Spring', 'Summer', 'Fall', 'Summer', 'Winter',
                       'Summer', 'Summer', 'Summer', 'Fall', 'Summer']
        self.diversity = [35.859310, 24.812545, 27.124025, 37.737965,
                          21.372255, 28.351530, 33.039095, 31.414730,
                          32.477700, 28.428650, 36.156570, 43.669460,
                          35.108995, 40.236090, 36.765910, 34.283980,
                          38.353880, 39.716300, 38.772970, 37.266020,
                          36.671575, 39.574585, 28.102095, 42.354890,
                          36.896660, 20.573455, 35.285450, 38.543620,
                          28.068870, 34.339165, 33.854125, 38.554955,
                          32.109585, 33.154860, 34.184950, 21.879765,
                          25.228975, 37.243695, 35.061765, 40.311775,
                          39.252515, 22.082885, 36.343690, 34.477320,
                          30.286625, 26.227675, 42.740985, 27.927150,
                          33.649625, 38.621705, 33.673875, 36.506995,
                          26.550305, 37.921820, 32.541755, 19.373825,
                          26.243955, 39.665195, 31.505100, 28.442290]

        # Creates a data frame object
        self.df = DataFrame({'AGE': Series(self.age, index=self.ids),
                             'GENDER': Series(self.gender, index=self.ids),
                             'SEASON': Series(self.season, index=self.ids),
                             'DIVERSITY': Series(self.diversity,
                                                 index=self.ids)})

    def test_pad_index_default(self):
        """Tests that a set of sample ids can be update sanely for defaults"""
        # Creates a data frame with raw ids and no sample column
        df = DataFrame({'#SampleID': Series(self.raw_ids),
                        'AGE': Series(self.age),
                        'GENDER': Series(self.gender)})
        # Pads the raw text
        df = pad_index(df)
        assert_index_equal(self.df.index, df.index)

    def test_pad_index_custom_index(self):
        """Tests index column can be set with pad index"""
        # Creates a data frame with raw ids and no sample column
        df = DataFrame({'RawID': Series(self.raw_ids),
                        'AGE': Series(self.age),
                        'GENDER': Series(self.gender)})
        # Pads the raw text
        df = pad_index(df, index_col='RawID')
        assert_index_equal(self.df.index, df.index)

    def test_pad_index_number(self):
        """Tests index column can be padded with """
        # Creates a data frame with raw ids and no sample column
        df = DataFrame({'#SampleID': Series(self.raw_ids),
                        'AGE': Series(self.age),
                        'GENDER': Series(self.gender)})
        # Pads the raw text
        df = pad_index(df, nzeros=4)
        assert_index_equal(Index(self.raw_ids), df.index)

    def test_check_dir(self):
        """Checks a directory is created if appropriate"""
        # Sets up a dummy directory that does not exist
        does_not_exist = pjoin(TEST_DIR, 'this_dir_does_not_exist')
        # Checks the directory does not currently exist
        self.assertFalse(exists(does_not_exist))
        # checks the directory
        check_dir(does_not_exist)
        # Checks the directory exists now
        self.assertTrue(exists(does_not_exist))
        # Removes the directory
        rmdir(does_not_exist)

if __name__ == '__main__':
    main()
