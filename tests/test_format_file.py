#!/usr/bin/env python

from unittest import TestCase, main
from StringIO import StringIO
from os.path import join

from americangut.format_file import do_replacements

__author__ = "Adam Robbins-Pianka"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Adam Robbins-Pianka"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Adam Robbins-Pianka"
__email__ = "adam.robbinspianka@colorado.edu"

class FormatFileTests(TestCase):
    def setUp(self):
        self.support_files_dir = join('support_files', 'format_file')
        self.insertion_file_path = join(self.support_files_dir,
            'test_insertion_file.txt')
        self.test_target_string = TARGET_FILE_TEXT

    def testDoReplacements(self):
        """Test the functionality of find-and-replace and find-and-insert
        functionality of format_file.do_replacements.
        """
        kfr = ['#VALUE#']
        vfr = ['WORD']

        kfi = ['#HERE#']
        vfi = [self.insertion_file_path]
        result = do_replacements(self.test_target_string, kfr, vfr, kfi, vfi)
        self.assertEqual(result, EXPECTED_TEXT)

TARGET_FILE_TEXT = '''There should be a #VALUE# (WORD) placed here after
do_replacements.
<
#HERE#
>
Above, between <angle brackets>, is the contents of test_insertion_file.txt,
including the substitution that is in that file.'''

EXPECTED_TEXT = '''There should be a WORD (WORD) placed here after
do_replacements.
<
this is

a test

insertion file that has a place for a WORD (WORD) in it

>
Above, between <angle brackets>, is the contents of test_insertion_file.txt,
including the substitution that is in that file.'''

if __name__ == '__main__':
    main()
