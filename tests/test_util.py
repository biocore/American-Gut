#!/usr/bin/env python

import os

import pandas as pd

from StringIO import StringIO
from unittest import TestCase, main

from numpy import array, nan, arange
from numpy.testing import assert_almost_equal
from biom import Table
from pandas.util.testing import assert_frame_equal

from americangut.util import (
    slice_mapping_file, parse_mapping_file,
    verify_subset, concatenate_files, trim_fasta, count_samples,
    count_seqs, count_unique_participants, clean_and_reformat_mapping,
    add_alpha_diversity, get_single_id_lists, collapse_taxonomy, collapse_full
)

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Daniel McDonald", "Adam Robbins-Pianka"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"


class UtilTests(TestCase):
    def test_count_samples(self):
        test_mapping = ["#SampleID\tfoo\tbar",
                        "A\t1\t2",
                        "B\t1\t3",
                        "C\t2\t4",
                        "D\t3\t5",
                        "E\t2\t6"]
        obs = count_samples(iter(test_mapping))
        exp = 5
        self.assertEqual(obs, exp)

        obs = count_samples(iter(test_mapping), criteria={'foo': '2'})
        exp = 2

    def test_count_seqs(self):
        test_seqs = [">a b",
                     "aattggcc",
                     ">b.xyz stuff",
                     "asdasd",
                     ">c",
                     "qweasd",
                     ">d.foo",
                     "qweasdasd"]

        obs = count_seqs(iter(test_seqs))
        exp = 4
        self.assertEqual(obs, exp)
        obs = count_seqs(iter(test_seqs), subset=['b', 'c', 'foo'])
        exp = 2
        self.assertEqual(obs, exp)

    def test_count_unique_participants(self):
        test_mapping = ["#SampleID\tfoo\tbar\tHOST_SUBJECT_ID",
                        "A\t1\t2\tx",
                        "B\t1\t3\tx",
                        "C\t2\t4\ty",
                        "D\t3\t5\tz",
                        "E\t2\t6\tw"]
        obs = count_unique_participants(iter(test_mapping))
        exp = 4
        self.assertEqual(obs, exp)

        obs = count_unique_participants(iter(test_mapping),
                                        criteria={'foo': '1'})
        exp = 1
        self.assertEqual(obs, exp)

        obs = count_unique_participants(iter(test_mapping),
                                        criteria={'foo': '2'})
        exp = 2
        self.assertEqual(obs, exp)

    def test_verify_subset(self):
        metadata = [('a','other stuff\tfoo'), ('b', 'asdasdasd'),
                    ('c','123123123')]
        table = Table(array([[1,2,3],[4,5,6]]),
                      ['x', 'y'],
                      ['a', 'b', 'c'])
        self.assertTrue(verify_subset(table, metadata))

        table = Table(array([[1,2],[3,4]]),
                      ['x','y'],
                      ['a','b'])
        self.assertTrue(verify_subset(table, metadata))

        table = Table(array([[1,2,3],[4,5,6]]),
                      ['x','y'],
                      ['a','b','x'])
        self.assertFalse(verify_subset(table, metadata))

    def test_slice_mapping_file(self):
        header, metadata = parse_mapping_file(StringIO(test_mapping))
        table = Table(array([[1,2],[4,5]]),
                      ['x','y'],
                      ['a','c'])

        exp = ["a\t1\t123123", "c\tpoop\tdoesn't matter"]
        obs = slice_mapping_file(table, metadata)
        self.assertEqual(obs,exp)

    def test_parse_mapping_file(self):
        exp = ("#SampleIDs\tfoo\tbar", [['a','1\t123123'],
                                        ['b','yy\txxx'],
                                        ['c',"poop\tdoesn't matter"]])
        obs = parse_mapping_file(StringIO(test_mapping))
        self.assertEqual(obs, exp)

    def test_concatenate_files(self):
        expected_output = concat_test_input + concat_test_input

        input_files = [StringIO(concat_test_input),
                       StringIO(concat_test_input)]
        output_file = StringIO()
        concatenate_files(input_files, output_file)
        output_file.seek(0)
        self.assertEqual(expected_output, output_file.read())

        # try again with a tiny chunk size
        input_files = [StringIO(concat_test_input),
                       StringIO(concat_test_input)]
        output_file = StringIO()
        concatenate_files(input_files, output_file, 2)
        output_file.seek(0)
        self.assertEqual(expected_output, output_file.read())

    def test_trim_fasta(self):
        infasta = StringIO(test_fasta)

        # Trim length 10
        expected = (">seq1\n"
                    "0123456789\n"
                    ">seq2\n"
                    "0123456789\n"
                    ">seq3\n"
                    "012345\n")
        outfasta = StringIO()
        trim_fasta(infasta, outfasta, 10)
        outfasta.seek(0)
        self.assertEqual(expected, outfasta.read())

    def test_clean_and_reformat_mapping(self):
        """Exercise the reformat mapping code, verify expected results"""
        out = StringIO()
        reformat_mapping_testdata.seek(0)
        clean_and_reformat_mapping(reformat_mapping_testdata, out, 'body_site',
                                   'test')
        out.seek(0)

        # verify the resulting header structure
        test_mapping = [l.strip().split('\t') for l in out]
        test_header = test_mapping[0]
        self.assertEqual(test_header[-4:], ['SIMPLE_BODY_SITE',
                                            'TITLE_ACRONYM', 'TITLE_BODY_SITE',
                                            'HMP_SITE'])

        self.assertEqual(test_mapping[1][:], ['A', 'w00t', '43.0',
                                              'UBERON_mucosa_of_tongue', '5',
                                              'ORAL', 'test', 'test-ORAL',
                                              'ORAL'])
        self.assertEqual(test_mapping[2][:], ['B', 'left', '51.0',
                                              'UBERON:FECES', '10',
                                              'FECAL', 'test', 'test-FECAL',
                                              'FECAL'])
        self.assertEqual(test_mapping[3][:], ['C', 'right', '12.0',
                                              'UBERON_FECES', '15',
                                              'FECAL', 'test', 'test-FECAL',
                                              'FECAL'])
        self.assertEqual(test_mapping[4][:], ['E', 'stuff', '56.0',
                                              'UBERON:SKIN', '37',
                                              'SKIN', 'test', 'test-SKIN',
                                              'SKIN'])

    def test_clean_and_reformat_mapping_nopgp(self):
        """Exercise the reformat mapping code, verify expected results"""
        out = StringIO()
        reformat_mapping_testdata.seek(0)
        clean_and_reformat_mapping(reformat_mapping_testdata, out, 'body_site',
                                   'test')
        out.seek(0)

        # verify the resulting header structure
        test_mapping = [l.strip().split('\t') for l in out]
        test_header = test_mapping[0]
        self.assertEqual(test_header[-4:], ['SIMPLE_BODY_SITE',
                                            'TITLE_ACRONYM', 'TITLE_BODY_SITE',
                                            'HMP_SITE'])

        self.assertEqual(test_mapping[1][:], ['A', 'w00t', '43.0',
                                              'UBERON_mucosa_of_tongue', '5',
                                              'ORAL', 'test', 'test-ORAL',
                                              'ORAL'])
        self.assertEqual(test_mapping[2][:], ['B', 'left', '51.0',
                                              'UBERON:FECES', '10',
                                              'FECAL', 'test', 'test-FECAL',
                                              'FECAL'])
        self.assertEqual(test_mapping[3][:], ['C', 'right', '12.0',
                                              'UBERON_FECES', '15',
                                              'FECAL', 'test', 'test-FECAL',
                                              'FECAL'])
        self.assertEqual(test_mapping[4][:], ['E', 'stuff', '56.0',
                                              'UBERON:SKIN', '37',
                                              'SKIN', 'test', 'test-SKIN',
                                              'SKIN'])

    def test_clean_and_reformat_mapping_allpgp(self):
        """Exercise the reformat mapping code, verify expected results"""
        out = StringIO()
        reformat_mapping_testdata.seek(0)
        clean_and_reformat_mapping(reformat_mapping_testdata, out, 'body_site',
                                   'test')
        out.seek(0)

        # verify the resulting header structure
        test_mapping = [l.strip().split('\t') for l in out]
        test_header = test_mapping[0]
        self.assertEqual(test_header[-4:], ['SIMPLE_BODY_SITE',
                                            'TITLE_ACRONYM', 'TITLE_BODY_SITE',
                                            'HMP_SITE'])

        self.assertEqual(test_mapping[1][:], ['A', 'w00t', '43.0',
                                              'UBERON_mucosa_of_tongue', '5',
                                              'ORAL', 'test', 'test-ORAL',
                                              'ORAL'])
        self.assertEqual(test_mapping[2][:], ['B', 'left', '51.0',
                                              'UBERON:FECES', '10',
                                              'FECAL', 'test', 'test-FECAL',
                                              'FECAL'])
        self.assertEqual(test_mapping[3][:], ['C', 'right', '12.0',
                                              'UBERON_FECES', '15',
                                              'FECAL', 'test', 'test-FECAL',
                                              'FECAL'])
        self.assertEqual(test_mapping[4][:], ['E', 'stuff', '56.0',
                                              'UBERON:SKIN', '37',
                                              'SKIN', 'test', 'test-SKIN',
                                              'SKIN'])

    def test_add_alpha_diversity(self):
        map_ = pd.DataFrame(
            array([
                ['GAZ:w00t', '43.0', 'UBERON_mucosa_of_tongue',  '5'],
                ['GAZ:left', '51.0', 'UBERON:FECES', '10'],
                ['GAZ:right', '12.0', 'UBERON_FECES', '15'],
                ['GAZ:stuff', '32.0', 'unknown', '26'],
                ['GAZ:stuff', '56.0', 'UBERON:SKIN', '37'],
            ]),
            columns=['COUNTRY', 'AGE', 'BODY_SITE', 'BMI'],
            index=['A', 'B', 'C', 'D', 'E']
        )
        alpha = {
            'alpha_1': pd.DataFrame(
                array([
                    ['0', '1', '2', '3', '4'],
                    ['100', '100', '100', '100', '100'],
                    [nan, nan, nan, nan, nan],
                    ['14.5', '14.0', '15.1', '14.7', '14.4'],
                    ['12.1', '15.2', '13.1', '14.1', '12.8'],
                    ['16.2', '16.5', '16.9', '15.9', '16.2'],
                    ['10.1',  '9.8', '10.5', '10.0', '10.2'],
                    ]),
                columns=[
                    'alpha_rarefaction_100_0.txt',
                    'alpha_rarefaction_100_1.txt',
                    'alpha_rarefaction_100_2.txt',
                    'alpha_rarefaction_100_3.txt',
                    'alpha_rarefaction_100_4.txt',
                ],
                index=['sequences per sample', 'iteration',
                       'A', 'B', 'C', 'D', 'E']
                )
        }
        expected = pd.DataFrame(
            array([
                ['GAZ:left', '51.0', 'UBERON:FECES', '10', 14.54],
                ['GAZ:right', '12.0', 'UBERON_FECES', '15', 13.46],
                ['GAZ:stuff', '32.0', 'unknown', '26', 16.34],
                ['GAZ:stuff', '56.0', 'UBERON:SKIN', '37', 10.12]
                ]),
            index=['B', 'C', 'D', 'E'],
            columns=['COUNTRY', 'AGE', 'BODY_SITE', 'BMI', 'alpha_1']
            )
        expected['alpha_1'] = expected['alpha_1'].astype(float)

        test = add_alpha_diversity(map_, alpha)
        assert_frame_equal(expected, test)

    def test_get_single_id_list(self):
        map_ = pd.DataFrame(
            array([
                ['GAZ:w00t', '43.0', 'UBERON_mucosa_of_tongue',  '5', 'A',
                 '12'],
                ['GAZ:left', '51.0', 'UBERON:FECES', '10', 'B', '1500'],
                ['GAZ:right', '12.0', 'UBERON_FECES', '15', 'C', '121'],
                ['GAZ:stuff', '32.0', 'unknown', '26', 'D', '150'],
                ['GAZ:stuff', '56.0', 'UBERON:SKIN', '37', 'E', '201'],
            ]),
            columns=['COUNTRY', 'AGE', 'BODY_SITE', 'BMI', 'HOST_SUBJECT_ID',
                     'depth'],
            index=['A', 'B', 'C', 'D', 'E']
        )
        depths = [100]
        test = get_single_id_lists(map_, depths)
        known = {100: ['B', 'C', 'D', 'E'],
                 'unrare': ['A', 'B', 'C', 'D', 'E']}
        self.assertEqual(test, known)

    def test_collapse_taxonomy(self):
        obs = collapse_taxonomy(table)
        exp = Table(array([[100.0,  105.0,  110.0,  115.0],
                           [44.0,   46.0,   48.0,   50.0],
                           [36.0,   39.0,   42.0,   45.0]]),
                    ['Bacteria; Firmicutes', 'Bacteria; Bacteroidetes',
                     'Bacteria; Proteobacteria'], sample_ids,
                    sample_metadata=sample_metadata, observation_metadata=[
                    {'collapsed_ids': ['O0', 'O1', 'O7', 'O8', 'O9']},
                    {'collapsed_ids': ['O5', 'O6']},
                    {'collapsed_ids': ['O2', 'O3', 'O4']}])
        self.assertEqual(obs, exp)

    def test_collapse_full(self):
        obs = collapse_full(table)
        exp = Table(array([[0.00769230769231], [0.0282051282051],
                           [0.0487179487179], [0.0692307692308],
                           [0.0897435897436], [0.110256410256],
                           [0.130769230769], [0.151282051282],
                           [0.171794871795], [0.192307692308]]),
                    observ_ids, ['average'],
                    observation_metadata=observ_metadata)
        for r in range(10):
            assert_almost_equal(obs[r, 0],  exp[r, 0])
        self.assertEqual(obs.ids(), exp.ids())
        self.assertItemsEqual(obs.ids('observation'), exp.ids('observation'))

        obs_meta = []
        for _, _, m in obs.iter(axis='observation'):
            obs_meta.append(m)
        self.assertItemsEqual(obs_meta, observ_metadata)

test_mapping = """#SampleIDs\tfoo\tbar
a\t1\t123123
b\tyy\txxx
c\tpoop\tdoesn't matter
"""

concat_test_input="""This is
a
test file that is used

in the concatenation test. The file will be concatenated to itself."""

test_fasta = """>seq1
0123456789
>seq2
0123456789AB
>seq3
012345"""

reformat_mapping_testdata = StringIO(
"""#SampleID	COUNTRY	AGE	BODY_SITE	BMI
A	GAZ:w00t	43.0	UBERON_mucosa_of_tongue	5
B	GAZ:left	51.0	UBERON:FECES	10
C	GAZ:right	12.0	UBERON_FECES	15
D	GAZ:stuff	32.0	unknown	26
E	GAZ:stuff	56.0	UBERON:SKIN	37
""")

data = arange(40).reshape(10, 4)
sample_ids = ['S%d' % i for i in range(4)]
observ_ids = ['O%d' % i for i in range(10)]
sample_metadata = [{'environment': 'A'}, {'environment': 'B'},
                   {'environment': 'A'}, {'environment': 'B'}]
observ_metadata = [{'taxonomy': ['Bacteria', 'Firmicutes']},
                   {'taxonomy': ['Bacteria', 'Firmicutes']},
                   {'taxonomy': ['Bacteria', 'Proteobacteria']},
                   {'taxonomy': ['Bacteria', 'Proteobacteria']},
                   {'taxonomy': ['Bacteria', 'Proteobacteria']},
                   {'taxonomy': ['Bacteria', 'Bacteroidetes']},
                   {'taxonomy': ['Bacteria', 'Bacteroidetes']},
                   {'taxonomy': ['Bacteria', 'Firmicutes']},
                   {'taxonomy': ['Bacteria', 'Firmicutes']},
                   {'taxonomy': ['Bacteria', 'Firmicutes']}]
table = Table(data, observ_ids, sample_ids, observ_metadata,
              sample_metadata, table_id='Example Table')


if __name__ == '__main__':
    main()
