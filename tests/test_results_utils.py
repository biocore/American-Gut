#!/usr/bin/env python

import os
from StringIO import StringIO
from unittest import TestCase, main
from americangut.results_utils import filter_mapping_file, massage_mapping


class ResultsUtilsTests(TestCase):
    def test_filter_mapping_file(self):
        output = StringIO()
        # filter to just fecal samples, keep the age and title_acronym columns
        criteria = {'SIMPLE_BODY_SITE': lambda x: x == 'FECAL',
                    'AGE': lambda x: float(x) > 20,
                    'TITLE_ACRONYM': None}
        filter_mapping_file(filter_mapping_testdata, output, criteria)
        output.seek(0)

        # parse output
        test_mapping = [l.strip().split('\t') for l in output]

        # fish header, verify sanity of it
        test_header = test_mapping[0]
        self.assertEqual(len(test_header), 4)
        self.assertEqual(test_header[0], '#SampleID')
        self.assertEqual(sorted(test_header), sorted(['#SampleID',
                                                      'SIMPLE_BODY_SITE',
                                                      'AGE', 'TITLE_ACRONYM']))

        # check each record
        test_sbs = test_header.index('SIMPLE_BODY_SITE')
        test_age = test_header.index('AGE')
        for l in test_mapping[1:]:
            self.assertEqual(len(l), 4)
            self.assertEqual(l[test_sbs], 'FECAL')
            self.assertTrue(float(l[test_age]) > 20)

    def test_massage_mapping(self):
        """Exercise the massage mapping code, verify expected results"""
        out = StringIO()
        massage_mapping(massage_mapping_testdata, out, 'body_site', 'test')
        out.seek(0)

        # verify the resulting header structure
        test_mapping = [l.strip().split('\t') for l in out]
        test_header = test_mapping[0]
        self.assertEqual(test_header[-6:], ['SIMPLE_BODY_SITE', 'TITLE_ACRONYM',
                                            'TITLE_BODY_SITE', 'HMP_SITE',
                                            'AGE_CATEGORY', 'BMI_CATEGORY'])

        self.assertEqual(test_mapping[1][:], ['A', 'w00t', '43.0',
                                              'UBERON_mucosa_of_tongue', '5',
                                              'ORAL', 'test', 'test-ORAL',
                                              'ORAL', '40s', 'Underweight'])
        self.assertEqual(test_mapping[2][:], ['B', 'left', '51.0',
                                              'UBERON:FECES', '10',
                                              'FECAL', 'test', 'test-FECAL',
                                              'FECAL', '50s', 'Underweight'])
        self.assertEqual(test_mapping[3][:], ['C', 'right', '12.0',
                                              'UBERON_FECES', '15',
                                              'FECAL', 'test', 'test-FECAL',
                                              'FECAL', 'Child', 'Underweight'])
        self.assertEqual(test_mapping[4][:], ['E', 'stuff', '56.0',
                                              'UBERON:SKIN', '37',
                                              'SKIN', 'test', 'test-SKIN',
                                              'SKIN', '50s', 'Severely obese'])

filter_mapping_testdata = StringIO(
"""#SampleID	COUNTRY	TITLE_ACRONYM	AGE	SIMPLE_BODY_SITE
A	United States of America	AGP	43.0	ORAL
B	United States of America	foo	51.0	FECAL
C	United States of America	bar	12.0	FECAL
D	United States of America	AGP	32.0	SKIN
E	United States of America	AGP	56.0	FECAL
""")
massage_mapping_testdata = StringIO(
"""#SampleID	COUNTRY	AGE	BODY_SITE	BMI
A	GAZ:w00t	43.0	UBERON_mucosa_of_tongue	5
B	GAZ:left	51.0	UBERON:FECES	10
C	GAZ:right	12.0	UBERON_FECES	15
D	GAZ:stuff	32.0	unknown	26
E	GAZ:stuff	56.0	UBERON:SKIN	37
""")
if __name__ == '__main__':
    main()
