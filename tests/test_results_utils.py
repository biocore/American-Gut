#!/usr/bin/env python

from StringIO import StringIO
from unittest import TestCase, main
from collections import defaultdict

import numpy as np
import numpy.testing as npt
import pandas as pd

from americangut.results_utils import (
    filter_mapping_file, count_unique_sequences_per_otu,
    write_bloom_fasta, plot_alpha
)


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

    def test_count_unique_sequences_per_otu(self):
        input_fasta = StringIO(test_fasta)
        otu_map = StringIO(test_otu_map)

        otu_ids = set(['otu1', 'otu2'])

        result = count_unique_sequences_per_otu(otu_ids, otu_map, input_fasta)

        expected = {x:defaultdict(int) for x in otu_ids}
        expected['otu1']['ATCG'] = 3
        expected['otu2']['AT'] = 2
        expected['otu2']['A'] = 1

        self.assertEqual(expected, result)

    def test_write_bloom_fasta(self):
        otu_ids = set(['otu1', 'otu2'])
        unique_counts = {x:defaultdict(int) for x in otu_ids}
        unique_counts['otu1']['ATCG'] = 3
        unique_counts['otu2']['AT'] = 2
        unique_counts['otu2']['A'] = 1

        result = StringIO()
        write_bloom_fasta(unique_counts, result, 0.67)

        result.seek(0)
        self.assertEqual(result.read(), '>otu1_1\nATCG\n')

    def test_plot_alpha(self):
        map_ = pd.DataFrame(
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
                       alpha_map=map_,
                       alpha_field='alpha',
                       debug=True)

        self.assertEqual(tgroup, kgroup)
        npt.assert_array_equal(tgalpha.values, kgalpha.values)
        self.assertEqual(tsalpha, ksalpha)
        self.assertEqual(txlabel, kxlabel)


filter_mapping_testdata = StringIO(
"""#SampleID	COUNTRY	TITLE_ACRONYM	AGE	SIMPLE_BODY_SITE
A	United States of America	AGP	43.0	ORAL
B	United States of America	foo	51.0	FECAL
C	United States of America	bar	12.0	FECAL
D	United States of America	AGP	32.0	SKIN
E	United States of America	AGP	56.0	FECAL
""")

# Inputs for count_unique_seqs_per_otu
test_fasta = """>sample1_1 yea
ATCG
>sample1_2 awyea
ATCG
>sample2_1 dumb
ATCG
>sample2_2 dummy
AT
>sample2_3 wow
AT
>sample2_4 wowagain
A
>sample9_1
ATGC
>sample9_2
A
"""

test_otu_map = """otu1	sample1_1	sample1_2	sample2_1
otu2	sample2_2	sample2_3	sample2_4
otu3	sample9_1	smaple9_2
"""

if __name__ == '__main__':
    main()
