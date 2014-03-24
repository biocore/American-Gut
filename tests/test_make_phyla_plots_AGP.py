# #!/usr/bin/env python

from unittest import TestCase, main
from StringIO import StringIO
from numpy import array
from biom.table import table_factory, SparseOTUTable
from americangut.make_phyla_plots import (map_to_2D_dict,
                                          most_common_taxa_gg_13_5,
                                          summarize_human_taxa)


__author__ = "Justine Debelius"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Justine Debelius"]
__license__ = "BSD" 
__version__ = "unversioned"
__maintainer__ = "Justine Debelius"
__email__ = "Justine.Debelius@colorado.edu"



class MakePhylaPlotsAGPTest(TestCase):
    
    def setUP(self):
        pass     

    def tearDown(self):
        #
        pass
        

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
 
    def test_most_common_taxa_gg_13_5(self):
        # Sets up the known value for phylum level taxa. If taxa change the 
        # values here will change.        
        phyla_known = [(u'k__Bacteria', u' p__Firmicutes'),
                       (u'k__Bacteria', u' p__Bacteroidetes'),
                       (u'k__Bacteria', u' p__Proteobacteria'),
                       (u'k__Bacteria', u' p__Actinobacteria'),
                       (u'k__Bacteria', u' p__Verrucomicrobia'),
                       (u'k__Bacteria', u' p__Tenericutes'),
                       (u'k__Bacteria', u' p__Cyanobacteria'),
                       (u'k__Bacteria', u' p__Fusobacteria'),
                       (u'k__Bacteria', u' p__Other')]

        # Checks that a level of 2 gives the phyla
        test_phyla = most_common_taxa_gg_13_5()
        self.assertEqual(test_phyla, phyla_known)

    def test_summarize_human_taxa(self):
       # Defines the known values
        known_ids = ['00010', '00100', '00200', '00111', '00112', '00211']

        known_taxa = [(u'k__Bacteria', u' p__Firmicutes'),
                      (u'k__Bacteria', u' p__Bacteroidetes'),
                      (u'k__Bacteria', u' p__Proteobacteria'),
                      (u'k__Bacteria', u' p__Actinobacteria'),
                      (u'k__Bacteria', u' p__Verrucomicrobia'),
                      (u'k__Bacteria', u' p__Tenericutes'),
                      (u'k__Bacteria', u' p__Cyanobacteria'),
                      (u'k__Bacteria', u' p__Fusobacteria'),
                      (u'k__Bacteria', u' p__Other')]
        
        table_known = array([[0.5000, 0.2500, 0.5500, 0.3000, 0.45000, 0.1000],
                             [0.3000, 0.5000, 0.3500, 0.6200, 0.45000, 0.7500],
                             [0.0000, 0.0000, 0.0000, 0.0000, 0.00000, 0.0000],
                             [0.0800, 0.0000, 0.0000, 0.0030, 0.00000, 0.0000],
                             [0.0500, 0.0100, 0.0020, 0.0010, 0.00000, 0.0000],
                             [0.0200, 0.0800, 0.0500, 0.0020, 0.05000, 0.0000],
                             [0.0010, 0.0200, 0.0000, 0.0000, 0.00010, 0.0000],
                             [0.0080, 0.0000, 0.0020, 0.0000, 0.00000, 0.0500],
                             [0.0210, 0.1100, 0.0460, 0.0240, 0.01990, 0.0000]])

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

        otu_table = table_factory(data = data, sample_ids = sample_ids, \
                                  observation_ids = observation_ids, \
                                  observation_metadata = observation_md, \
                                  constructor = SparseOTUTable)

        # Tests that the table is summarized correctly
        (test_taxa, test_ids, test_table) = summarize_human_taxa(otu_table, 2)

        # Checks that all the outputs are correct
        self.assertEqual(test_ids, known_ids)
        self.assertEqual(test_taxa, known_taxa)
        self.assertEqual(test_table.all(), table_known.all())


    def test_format_date(self):
        """Test formating the date information for a metadata dictionary"""
        # Sets the locations of the time information in the metadata
        date_field = 'Sample_Date'
        time_field = 'Sample_Time'

        # Checks that errors are called appropriately
        with self.assertRaises(ValueError):
            format_date(self.meta)

        with self.assertRaises(ValueError):
            format_date(self.meta, date_field='Tardis', d_form_in='%Y')

        with self.assertRaises(ValueError):
            format_date(self.meta, date_field=date_field)

        with self.assertRaises(ValueError):
            format_date(self.meta, time_field='Smith', t_form_in='%H:%M')

        with self.assertRaises(ValueError):
            format_date(self.meta, time_field=time_field)

        # Generates and checks the test values for a single date
        known = '17 Feb 1963'
        format_in = '%m/%d/%Y'
        format_out = '%d %b %Y'
        test = format_date(self.meta, date_field=date_field, 
                           d_form_in=format_in, format_out=format_out)
        self.assertEqual(test, known)

        # Generates and checks the test value for a time only
        known = '03:27'
        format_in = '%I:%M %p'
        format_out = '%H:%M'
        test = format_date(self.meta, time_field=time_field, 
                           t_form_in=format_in, format_out=format_out)
        self.assertEqual(test, known)

        # Checks the combination of reading a date and time format
        known = '17 Feb 1963 03:27 AM'
        d_format_in = '%m/%d/%Y'
        t_format_in = '%I:%M %p'
        format_out = '%d %b %Y %I:%M %p'
        test = format_date(self.meta, date_field=date_field, 
                           time_field=time_field, d_form_in=d_format_in, 
                           t_form_in=t_format_in, format_out=format_out)

        self.assertEqual(test, known)






if __name__ == '__main__':
    main()