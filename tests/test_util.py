#!/usr/bin/env python

import os
from StringIO import StringIO
from unittest import TestCase, main

from numpy import array
from biom import Table

from americangut.util import (
    slice_mapping_file,parse_mapping_file,
    verify_subset, concatenate_files, trim_fasta, count_samples,
    count_seqs, count_unique_participants, clean_and_reformat_mapping
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

if __name__ == '__main__':
    main()
