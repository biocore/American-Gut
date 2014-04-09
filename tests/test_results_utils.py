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

    def test_massage_mapping(ag_100nt):
        """Exercise the massage mapping code, verify expected results"""
        # output test file
        ag_100nt_m_TEST_fp = 'AG_100nt_TEST.txt'

        massage_mapping(ag_100nt, ag_100nt_m_TEST_fp, 'body_site', 'AGP')

        # verify the resulting header structure
        test_mapping = [l.strip().split('\t') for l in open(ag_100nt_m_TEST_fp)]
        test_header = test_mapping[0]
        test_header_length = len(test_header)
        assert test_header[-6:] == ['SIMPLE_BODY_SITE', 'TITLE_ACRONYM',
                                    'TITLE_BODY_SITE', 'HMP_SITE',
                                    'AGE_CATEGORY', 'BMI_CATEGORY']

        # verify each line in the test file
        for l in test_mapping[1:]:
            assert l[-6] in ['FECAL', 'SKIN', 'ORAL']
            assert l[-5] == 'AGP'
            acro, site = l[-4].split('-')
            assert acro == 'AGP'
            assert site in ['FECAL', 'SKIN', 'ORAL']
            assert l[-3] in ['FECAL', 'SKIN', 'ORAL']
            assert len(l) == test_header_length

filter_mapping_testdata = StringIO(
"""#SampleID	COUNTRY	TITLE_ACRONYM	AGE	SIMPLE_BODY_SITE
A	United States of America	AGP	43.0	ORAL
B	United States of America	foo	51.0	FECAL
C	United States of America	bar	12.0	FECAL
D	United States of America	AGP	32.0	SKIN
E	United States of America	AGP	56.0	FECAL
""")
if __name__ == '__main__':
    main()
