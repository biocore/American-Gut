#!/usr/bin/env python

from __future__ import division

from americangut.format_file import ArgumentError, do_replacements

from argparse import ArgumentParser

__author__ = "Adam Robbins-Pianka"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Adam Robbins-Pianka"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Adam Robbins-Pianka"
__email__ = "adam.robbinspianka@colorado.edu"

parser = ArgumentParser()
parser.add_argument('-i', '--input_file', help=('The file to be formatted'),
    required=True),
parser.add_argument('-k', '--keys_for_find_and_replace', help=('A comma '
    'separated list of string substitution keys.'))
parser.add_argument('-v', '--values_for_find_and_replace', help=('A comma '
    'separated list of values that will replace the keys specified with -k')),
parser.add_argument('-K', '--keys_for_find_and_insert', help=('A comma '
    'separated list of string substitution keys.')),
parser.add_argument('-V', '--values_for_find_and_insert', help=('A comma '
    'separated list of file paths. The contents of the files will replace the '
    'keys specified with -K.')),
parser.add_argument('-o', '--output_file', help=('Specify the location '
    'of the output formatted file.'), required=True)

if __name__ == '__main__':
    args = parser.parse_args()

    output_fp = args.output_file
    input_fp = args.input_file
    kfr = args.keys_for_find_and_replace
    vfr = args.values_for_find_and_replace
    kfi = args.keys_for_find_and_insert
    vfi = args.values_for_find_and_insert

    # Validate arguments
    if (bool(kfr) != bool(vfr)):
        raise ArgumentError, ('If you specify -k, you must also specify -v.')
    if (bool(kfi) != bool(vfi)):
        raise ArgumentError, ('If you specify -K, you must also specify -V.')

    if kfr:
        kfr = kfr.split(',')
        vfr = vfr.split(',')
        if len(kfr) != len(vfr):
            raise ArgumentError, ('Number of keys (-k) does not match number '
                'of values (-v)')
    else:
        kfr = []
        vfr = []

    if kfi:
        kfi = kfi.split(',')
        vfi = vfi.split(',')
        if len(kfi) != len(vfi):
            raise ArgumentError, ('Number of Keys (-K) does not match number '
                'of Values (-V)')
    else:
        kfi = []
        vfi = []

    output_string = do_replacements(open(input_fp).read(), kfr, vfr, kfi, vfi)

    output_fh = open(output_fp, 'w')
    output_fh.write(output_string)
    output_fh.close()
