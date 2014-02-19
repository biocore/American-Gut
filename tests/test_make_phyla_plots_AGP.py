#!/usr/bin/env python

from unittest import TestCase, main
from StringIO import StringIO
from numpy import array
from numpy.testing import assert_almost_equal
from biom.table import table_factory, SparseOTUTable
from americangut.make_phyla_plots import (map_to_2D_dict,
                                          identify_most_common_categories,
                                          summarize_common_categories,
                                          calculate_dimensions_rectangle,
                                          calculate_dimensions_bar,
                                          translate_colors)
from matplotlib.transforms import Bbox


__author__ = "Justine Debelius"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Justine Debelius"]
__license__ = "BSD" 
__version__ = "unversioned"
__maintainer__ = "Justine Debelius"
__email__ = "Justine.Debelius@colorado.edu"



class MakePhylaPlotsAGPTest(TestCase):
    
    def setUp(self):
        # Creates an otu table for testing which corresponds with the 
        sample_ids = ['00010', '00100', '00200', '00111', '00112', '00211']

        observation_ids = ['1001', '2001', '2002', '2003', '3001', '3003', 
                           '4001', '5001', '6001', '7001', '8001', '9001', 
                           '9002', '9003']

        observation_md = [{'taxonomy': (u'k__Bacteria', u' p__Bacteroidetes',
                                        u'c__Bacteroidia')},
                          {'taxonomy': (u'k__Bacteria', u' p__Firmicutes',
                                        u'c__Clostridia')},
                          {'taxonomy': (u'k__Bacteria', u' p__Firmicutes',
                                        u'c__Erysipelotrichi')},
                          {'taxonomy': (u'k__Bacteria', u' p__Firmicutes',
                                        u'c__Bacilli')},
                          {'taxonomy': (u'k__Bacteria', u' p__Proteobacteria',
                                        u'c__Alphaproteobacteria')},
                          {'taxonomy': (u'k__Bacteria', u' p__Proteobacteria', 
                                        u'c__Gammaproteobacteria')},
                          {'taxonomy': (u'k__Bacteria', u' p__Tenericutes', 
                                        u'c__Mollicutes')},
                          {'taxonomy': (u'k__Bacteria', u' p__Actinobacteria', 
                                        u'c__Coriobacteriia')},
                          {'taxonomy': (u'k__Bacteria', u' p__Verrucomicrobia', 
                                        u'c__Verrucomicrobiae')},
                          {'taxonomy': (u'k__Bacteria', u' p__Cyanobacteria', 
                                        u'c__4C0d-2')},
                          {'taxonomy': (u'k__Bacteria', u' p__Fusobacteria', 
                                        u'c__Fusobacteriia')},
                          {'taxonomy': (u'k__Bacteria', u' p__TM7', 
                                        u'c__TM7-2')},
                          {'taxonomy': (u'k__Bacteria', u' p__Acidobacteria', 
                                        u'c__Chloracidobacteria')},
                          {'taxonomy': (u'k__Bacteria', u' p__', u'c__')}]

        data = array([[ 1691,  3004, 18606,  6914,  1314, 22843],
                      [ 2019,  1091,  8163,  1112,   738,  2362],
                      [   67,     4,  2835,   310,    85,   161],
                      [  731,   407, 18240,  1924,   492,   522],
                      [    8,     1,     0,    53,     8,   275],
                      [  105,   179,     0,   504,    79,  2771],
                      [  451,     0,     0,    33,     0,     0],
                      [  282,    60,   106,    11,     0,     0],
                      [  113,   481,  2658,    22,   146,     0],
                      [    6,   120,     0,     0,     0,     0],
                      [   45,     0,   106,     0,     0,  1523],
                      [   39,   341,  1761,   139,    18,     0],
                      [   21,   268,   153,     8,    15,     0],
                      [   59,    51,   531,   120,    25,     0]])                     

        self.otu_table = table_factory(data = data, 
                                       sample_ids = sample_ids,
                                       observation_ids = observation_ids,
                                       observation_metadata = observation_md,
                                       constructor = SparseOTUTable)

        self.common_cats = [(u'k__Bacteria', u' p__Firmicutes'),
                            (u'k__Bacteria', u' p__Bacteroidetes'),
                            (u'k__Bacteria', u' p__Proteobacteria'),
                            (u'k__Bacteria', u' p__Actinobacteria'),
                            (u'k__Bacteria', u' p__Verrucomicrobia'),
                            (u'k__Bacteria', u' p__Tenericutes'),
                            (u'k__Bacteria', u' p__Cyanobacteria'),
                            (u'k__Bacteria', u' p__Fusobacteria')]

    def test_map_to_2D_dict(self):

        # Creates a pseudo-opening function
        test_map = StringIO(\
            '#SampleID\tBIRTH_YEAR\tDEATH_YEAR\tSEX\tPROFESSION\tHOME_STATE\n'\
            '00010\t1954\t2006\tmale\tMechanic\tKansas\n'\
            '00100\t1954\t1983\tfemale\tHunter\tKansas\n'\
            '00200\tNA\t2009\tfemale\tNurse\tMinnesota\n'\
            '00111\t1979\t2007\tmale\tHunter\tImpala\n'\
            '00112\t1983\t2006\tmale\tHunter\tImpala\n'\
            '00211\t1990\t2009\tmale\tStudent\tMinnesota\n')

        # Sets up the known dictionary
        known_dict = {'00010': {'#SampleID': '00010', 'BIRTH_YEAR':'1954', \
                                'DEATH_YEAR':'2006', 'SEX':'male', \
                                'PROFESSION':'Mechanic', 'HOME_STATE':'Kansas'},
                      '00100': {'#SampleID': '00100', 'BIRTH_YEAR':'1954', \
                                 'DEATH_YEAR':'1983', 'SEX':'female', \
                                 'PROFESSION':'Hunter', 'HOME_STATE':'Kansas'},
                      '00200': {'#SampleID': '00200', 'BIRTH_YEAR':'NA', \
                                 'DEATH_YEAR':'2009', 'SEX':'female', \
                                 'PROFESSION':'Nurse', \
                                 'HOME_STATE':'Minnesota'},
                      '00111': {'#SampleID': '00111', 'BIRTH_YEAR':'1979', \
                                'DEATH_YEAR':'2007', 'SEX':'male',\
                                'PROFESSION':'Hunter', 'HOME_STATE':'Impala'},
                      '00112': {'#SampleID': '00112', 'BIRTH_YEAR':'1983', \
                                'DEATH_YEAR':'2006', 'SEX':'male', \
                                'PROFESSION':'Hunter', 'HOME_STATE':'Impala'},
                      '00211': {'#SampleID': '00211', 'BIRTH_YEAR':'1990', \
                                'DEATH_YEAR':'2009', 'SEX':'male', \
                                'PROFESSION':'Student', 'HOME_STATE':\
                                'Minnesota'}}

        # Checks the test dictionary is loaded properly and equals the known
        test_dict = map_to_2D_dict(test_map)
        self.assertEqual(test_dict,known_dict)

    def test_identify_most_common_categories(self):
        # Sets up known values
        known_cats_comp = [(u'k__Bacteria', u' p__Bacteroidetes'),
                           (u'k__Bacteria', u' p__Firmicutes'),
                           (u'k__Bacteria', u' p__Proteobacteria'),
                           (u'k__Bacteria', u' p__Verrucomicrobia'),
                           (u'k__Bacteria', u' p__TM7'),
                           (u'k__Bacteria', u' p__Acidobacteria'),
                           (u'k__Bacteria', u' p__Actinobacteria'),
                           (u'k__Bacteria', u' p__'),
                           (u'k__Bacteria', u' p__Fusobacteria')]

        known_cats_aver = [(u'k__Bacteria', u' p__Bacteroidetes'), 
                           (u'k__Bacteria', u' p__Firmicutes')]

        known_cat_count = [(u'k__Bacteria', u' p__Bacteroidetes'), 
                           (u'k__Bacteria', u' p__Firmicutes'), 
                           (u'k__Bacteria', u' p__Proteobacteria'), 
                           (u'k__Bacteria', u' p__Verrucomicrobia'), 
                           (u'k__Bacteria', u' p__TM7'), 
                           (u'k__Bacteria', u' p__Acidobacteria'), 
                           (u'k__Bacteria', u' p__'), 
                           (u'k__Bacteria', u' p__Actinobacteria')]

        known_cats_none = [(u'k__Bacteria', u' p__'),
                           (u'k__Bacteria', u' p__Acidobacteria'),
                           (u'k__Bacteria', u' p__Actinobacteria'),
                           (u'k__Bacteria', u' p__Bacteroidetes'),
                           (u'k__Bacteria', u' p__Cyanobacteria'),
                           (u'k__Bacteria', u' p__Firmicutes'),
                           (u'k__Bacteria', u' p__Fusobacteria'),
                           (u'k__Bacteria', u' p__Proteobacteria'),
                           (u'k__Bacteria', u' p__TM7'),
                           (u'k__Bacteria', u' p__Tenericutes'),
                           (u'k__Bacteria', u' p__Verrucomicrobia')]
        
        known_scores_comp = [[(u'k__Bacteria', u' p__Bacteroidetes'),   0.4950, 
                              1.0000, 4950.00], 
                             [(u'k__Bacteria', u' p__Firmicutes'),      0.3584, 
                              1.0000, 3584.00], 
                             [(u'k__Bacteria', u' p__Proteobacteria'),  0.0383, 
                              0.8333,  319.15], 
                             [(u'k__Bacteria', u' p__Verrucomicrobia'), 0.0337, 
                              0.8333,  280.82], 
                             [(u'k__Bacteria', u' p__TM7'),             0.0192, 
                              0.8333,  159.99], 
                             [(u'k__Bacteria', u' p__Acidobacteria'),   0.0095, 
                              0.8333,   79.16], 
                             [(u'k__Bacteria', u' p__Actinobacteria'),  0.0105, 
                              0.6667,   70.00], 
                             [(u'k__Bacteria', u' p__'),                0.0080, 
                              0.8333,   66.66], 
                             [(u'k__Bacteria', u' p__Fusobacteria'),    0.0100, 
                              0.5000,  50.00],
                             [(u'k__Bacteria', u' p__Tenericutes'),     0.0138, 
                              0.3333,   46.00], 
                             [(u'k__Bacteria', u' p__Cyanobacteria'),   0.0035, 
                              0.3333,   11.67]]    

        known_scores_aver = [[(u'k__Bacteria', u' p__Bacteroidetes'),   0.4950, 
                               1.0000, 4950.00], 
                             [(u'k__Bacteria', u' p__Firmicutes'),      0.3584, 
                               1.0000, 3584.00], 
                             [(u'k__Bacteria', u' p__Proteobacteria'),  0.0383, 
                               0.8333,  319.15], 
                             [(u'k__Bacteria', u' p__Verrucomicrobia'), 0.0337, 
                               0.8333,  280.82], 
                             [(u'k__Bacteria', u' p__TM7'),             0.0192, 
                               0.8333,  159.99], 
                             [(u'k__Bacteria', u' p__Tenericutes'),     0.0138, 
                               0.3333,   46.00], 
                             [(u'k__Bacteria', u' p__Actinobacteria'),  0.0105, 
                               0.6667,   70.00], 
                             [(u'k__Bacteria', u' p__Fusobacteria'),    0.0100, 
                               0.5000,   50.00], 
                             [(u'k__Bacteria', u' p__Acidobacteria'),   0.0095, 
                               0.8333,   79.16], 
                             [(u'k__Bacteria', u' p__'),                0.0080, 
                               0.8333,   66.66], 
                             [(u'k__Bacteria', u' p__Cyanobacteria'),   0.0035, 
                               0.3333,   11.67]]

        known_score_count = [[(u'k__Bacteria', u' p__Bacteroidetes'),   0.4950, 
                               1.0000, 4950.00], 
                             [(u'k__Bacteria', u' p__Firmicutes'),      0.3584, 
                               1.0000, 3584.00], 
                             [(u'k__Bacteria', u' p__Proteobacteria'),  0.0383, 
                               0.8333,  319.15], 
                             [(u'k__Bacteria', u' p__Verrucomicrobia'), 0.0337, 
                               0.8333,  280.82], 
                             [(u'k__Bacteria', u' p__TM7'),             0.0192, 
                               0.8333,  159.99], 
                             [(u'k__Bacteria', u' p__Acidobacteria'),   0.0095, 
                               0.8333,   79.16], 
                             [(u'k__Bacteria', u' p__'),                0.0080, 
                               0.8333,   66.66], 
                             [(u'k__Bacteria', u' p__Actinobacteria'),  0.0105, 
                               0.6667,   70.00], 
                             [(u'k__Bacteria', u' p__Fusobacteria'),    0.0100, 
                               0.5000,   50.00], 
                             [(u'k__Bacteria', u' p__Tenericutes'),     0.0138, 
                               0.3333,   46.00], 
                             [(u'k__Bacteria', u' p__Cyanobacteria'),   0.0035, 
                               0.3333,   11.67]]

        known_scores_none = [[(u'k__Bacteria', u' p__'),                0.0080, 
                                0.8333,   66.66],
                             [(u'k__Bacteria', u' p__Acidobacteria'),   0.0095, 
                                0.8333,   79.16],
                             [(u'k__Bacteria', u' p__Actinobacteria'),  0.0105, 
                                0.6667,   70.00],
                             [(u'k__Bacteria', u' p__Bacteroidetes'),   0.4950, 
                                1.0000, 4950.00],
                             [(u'k__Bacteria', u' p__Cyanobacteria'),   0.0035, 
                                0.3333,   11.67],
                             [(u'k__Bacteria', u' p__Firmicutes'),      0.3584, 
                                1.0000, 3584.00],
                             [(u'k__Bacteria', u' p__Fusobacteria'),    0.0100,
                                0.5000,   50.00],
                             [(u'k__Bacteria', u' p__Proteobacteria'),  0.0383, 
                                0.8333,  319.15],
                             [(u'k__Bacteria', u' p__TM7'),             0.0192, 
                                0.8333,  159.99],
                             [(u'k__Bacteria', u' p__Tenericutes'),     0.0138, 
                                0.3333,   46.00],
                             [(u'k__Bacteria', u' p__Verrucomicrobia'), 0.0337, 
                                0.8333,  280.82]]

        # Tests code
        [test_cats_none, test_scores_none] = \
            identify_most_common_categories(biom_table = self.otu_table, 
                                            level = 2,
                                            metadata_category = 'taxonomy', 
                                            limit_mode = 'NONE')

        [test_cats_comp, test_scores_comp] = \
            identify_most_common_categories(biom_table = self.otu_table, 
                                            level = 2,
                                            metadata_category = 'taxonomy', 
                                            limit_mode = 'COMPOSITE',
                                            limit = 49)

        [test_cats_aver, test_scores_aver] = \
            identify_most_common_categories(biom_table = self.otu_table, 
                                            level = 2,
                                            metadata_category = 'taxonomy', 
                                            limit_mode = 'AVERAGE',
                                            limit = 0.1)        

        [test_cat_count, test_score_count] = \
            identify_most_common_categories(biom_table = self.otu_table, 
                                            level = 2,
                                            metadata_category = 'taxonomy', 
                                            limit_mode = 'COUNTS',
                                            limit = 0.5)

        # Checks that appropriate errors are called
        with self.assertRaises(ValueError):
            identify_most_common_categories(biom_table = self.otu_table, 
                                        level = 2,
                                        limit_mode = 'This is a test')

        with self.assertRaises(ValueError):
            identify_most_common_categories(biom_table = self.otu_table,
                                        level = 2,
                                        limit = 100000)

        # Checks that output values are correct
        self.assertEqual(test_cats_none, known_cats_none)
        self.assertEqual(test_scores_none, known_scores_none)
       
        self.assertEqual(test_cats_comp, known_cats_comp)
        self.assertEqual(test_scores_comp, known_scores_comp)

        self.assertEqual(test_cats_aver, known_cats_aver)
        self.assertEqual(test_scores_aver, known_scores_aver)
      
        self.assertEqual(test_cat_count, known_cat_count)
        self.assertEqual(test_score_count, known_score_count)
    
    def test_summarize_common_categories(self):
        # Defines the known values
        known_ids = ('00010', '00100', '00200', '00111', '00112', '00211')

        table_known = array([[ 0.49973390, 0.25004162, 0.55001035, 
                               0.30008969, 0.45034247, 0.09997702],
                             [ 0.29998226, 0.50008324, 0.35000658, 
                               0.62008969, 0.45000000, 0.75000821],
                             [ 0.02004612, 0.02996504, 0.00000000, 
                               0.04995516, 0.02979452, 0.10000985],
                             [ 0.05002661, 0.00998835, 0.00199402, 
                               0.00098655, 0.00000000, 0.00000000],
                             [ 0.02004612, 0.08007325, 0.05000094, 
                               0.00197309, 0.05000000, 0.00000000],
                             [ 0.08000710, 0.00000000, 0.00000000, 
                               0.00295964, 0.00000000, 0.00000000],
                             [ 0.00106440, 0.01997669, 0.00000000, 
                               0.00000000, 0.00000000, 0.00000000],
                             [ 0.00798297, 0.00000000, 0.00199402, 
                               0.00000000, 0.00000000, 0.05000492],
                             [ 0.02111052, 0.10987182, 0.04599409, 
                               0.02394619, 0.01986301, 0.00000000]])

        known_common_cats = [(u'k__Bacteria', u' p__Firmicutes'),
                             (u'k__Bacteria', u' p__Bacteroidetes'),
                             (u'k__Bacteria', u' p__Proteobacteria'),
                             (u'k__Bacteria', u' p__Actinobacteria'),
                             (u'k__Bacteria', u' p__Verrucomicrobia'),
                             (u'k__Bacteria', u' p__Tenericutes'),
                             (u'k__Bacteria', u' p__Cyanobacteria'),
                             (u'k__Bacteria', u' p__Fusobacteria'),
                             (u'k__Bacteria', u' p__Other')]

        # Checks that appropriate errors are raised when the wrong type of 
        # argument is passed.
        with self.assertRaises(ValueError):
            summarize_common_categories(biom_table = self.otu_table,
                                        level = 2,
                                        common_categories = self.common_cats,
                                        metadata_category = 'Billy_Joel_Song')


        # Calculates the test values
        [test_ids, test_table, test_common_cats] = summarize_common_categories(
                                        biom_table = self.otu_table, 
                                        level = 2, 
                                        common_categories = self.common_cats)        
        
        # Checks that all the outputs are correct
        self.assertEqual(test_ids, known_ids)
        assert_almost_equal(test_table, table_known, decimal = 4)
        self.assertEqual(test_common_cats, known_common_cats)

    def test_calculate_dimensions_rectangle(self):
        # Sets up known values
        known_figure_dimensions_def = (5.2, 4.45)
        known_axis_dimensions_def = Bbox(array([[0.01923077, 0.02247191],
                                                [0.78846154, 0.92134831]]))

        known_figure_dims_in = (4.2, 3.45)
        known_axis_dims_in = Bbox(array([[0.26190476, 0.31884058],
                                         [0.73809524, 0.89855072]]))

        known_figure_dims_cm = (1.6535433, 1.3582677)
        known_axis_dims_cm = Bbox(array([[0.26190476, 0.31884058],
                                         [0.73809524, 0.89855072]]))
        # Sets up test values
        test_axis_side = 2
        test_border = 0.1
        test_xlab = 1
        test_ylab = 1

        # Tests that an error is raised if the units are not sane
        with self.assertRaises(ValueError):
            calculate_dimensions_rectangle(unit = 'Demons')

        # Calculates the test values
        (test_axis_df, test_fig_df) = calculate_dimensions_rectangle()
        
        (test_axis_in, test_fig_in) = calculate_dimensions_rectangle(
                                                axis_width = test_axis_side,
                                                axis_height = test_axis_side,
                                                border = test_border,
                                                xlab = test_xlab,
                                                ylab = test_ylab)

        (test_axis_cm, test_fig_cm) = calculate_dimensions_rectangle(
                                                axis_width = test_axis_side,
                                                axis_height = test_axis_side,
                                                border = test_border,
                                                xlab = test_xlab,
                                                ylab = test_ylab,
                                                unit = 'cm')

        assert_almost_equal(test_fig_df, known_figure_dimensions_def, 
                            decimal = 5)
        assert_almost_equal(test_axis_df, known_axis_dimensions_def, 
                            decimal = 5)

        assert_almost_equal(test_fig_in, known_figure_dims_in, decimal = 5)
        assert_almost_equal(test_axis_in, known_axis_dims_in, decimal = 5)

        assert_almost_equal(test_fig_cm, known_figure_dims_cm, decimal = 5)
        assert_almost_equal(test_axis_cm, known_axis_dims_cm, decimal = 5)

    def test_calculate_dimensions_bar(self):
        # Sets up known values
        known_axis_df = array([[0.01388889, 0.02380952],
                               [0.70833333, 0.73809524]])

        known_axis_in = array([[0.06179775, 0.21153846],
                               [0.43258427, 0.98076923]])

        known_axis_cm = array([[0.06179775, 0.21153846],
                               [0.43258427, 0.98076923]])

        known_fig_df = (7.2, 4.2)
        known_fig_in = (8.9, 2.6)
        known_fig_cm = (8.9, 2.6)
        
        # Sets up test values
        test_num_bars = 10
        test_bar_width = 0.33
        test_axis_height = 2
        test_border = 0.05
        test_title = 0
        test_legend = 5
        test_xlab = 0.5
        test_ylab = 0.5


        # Checks that errors are raised when improper arguments are passed. (A 
        # number of bars that is not an integer less than one, or a unit other 
        # than 'in' or 'cm'.)
        with self.assertRaises(ValueError):
            calculate_dimensions_bar(num_bars = -3)

        with self.assertRaises(ValueError):
            calculate_dimensions_bar(num_bars = 3.1435)

        with self.assertRaises(ValueError):
            calculate_dimensions_bar(num_bars = 3, unit = 'Angel')


        # Calculates the test values
        (test_axis_df, test_fig_df) = calculate_dimensions_bar(test_num_bars)

        (test_axis_in, test_fig_in) = calculate_dimensions_bar(
                                                num_bars = test_num_bars,
                                                bar_width = test_bar_width,
                                                border = test_border, 
                                                axis_height = test_axis_height,
                                                xlab = test_xlab,
                                                ylab = test_ylab,
                                                title = test_title,                                       
                                                legend = test_legend)

        (test_axis_cm, test_fig_cm) = calculate_dimensions_bar(
                                                num_bars = test_num_bars,
                                                bar_width = test_bar_width, 
                                                border = test_border,
                                                axis_height = test_axis_height,
                                                title = test_title,
                                                xlab = test_xlab,
                                                ylab = test_ylab,
                                                legend = test_legend,
                                                unit = 'cm')

        # Checks that values are appropriate
        self.assertEqual(known_fig_df, test_fig_df)
        self.assertEqual(known_fig_in, test_fig_in)
        self.assertEqual(known_fig_cm, test_fig_cm)

        assert_almost_equal(known_axis_df, test_axis_df, decimal = 5)
        assert_almost_equal(known_axis_in, test_axis_in, decimal = 5)
        assert_almost_equal(known_axis_cm, test_axis_cm, decimal = 5)

    def test_translate_colors(self):
        # Sets up knowns values
        known_def_8 = array([[0.83529412, 0.24313725, 0.30980392],
                             [0.95686275, 0.42745098, 0.26274510],
                             [0.99215686, 0.68235294, 0.38039216],
                             [0.99607843, 0.87843137, 0.54509804],
                             [0.90196078, 0.96078431, 0.59607843],
                             [0.67058824, 0.86666667, 0.64313725],
                             [0.40000000, 0.76078431, 0.64705882],
                             [0.19607843, 0.53333333, 0.74117647]])

        known_PuRd_9 = array([[0.96862745, 0.95686275, 0.97647059],
                              [0.90588235, 0.88235294, 0.93725490],
                              [0.83137255, 0.72549020, 0.85490196],
                              [0.78823529, 0.58039216, 0.78039216],
                              [0.87450980, 0.39607843, 0.69019608],
                              [0.90588235, 0.16078431, 0.54117647],
                              [0.80784314, 0.07058824, 0.33725490],
                              [0.59607843, 0.00000000, 0.26274510],
                              [0.40392157, 0.00000000, 0.12156863]])
        
        # Test the calls
        with self.assertRaises(ValueError):
            translate_colors(5, 'Winchester')

        with self.assertRaises(ValueError):
            translate_colors(13)

        def_map = translate_colors(8)
        PuRd_map = translate_colors(9, 'PuRd')

        # Checks the outputs are sane
        assert_almost_equal(known_def_8, def_map)
        assert_almost_equal(known_PuRd_9, PuRd_map)

if __name__ == '__main__':
    main()