#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser()
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

class ArgumentError(Exception):
    """Generic error that can be raised if there is a problem detected with
    one or more of the arguments passed by the user.
    """
    pass

def do_replacements(input_fp, kfr, vfr, kfi, vfi):
    """Does find-and-replace and find-and-insert operations on an input file

    input_fp is a filepath to an input file

    kfr (keys for find-and-replace) is a list of keys
    vfr (values for find-and-replae) is a list of values that will replace
        the keys specified in kfr

    kfi (keys for find-and-insert) is a list of keys
    vfi (values for find-and-insert) is a list of file paths. The contents
        of those files will replace the keys specified in kfi

    The insertions are done first, so that the files specified in vfi can have
    keys specified in kfr, if desiered.

    Returns a string that is the contents of the input file with all
    replacements performed.
    """
    input_string = open(input_fp, 'U').read()

    for k, v in zip(kfi, vfi):
        contents = open(v, 'U').read()
        input_string = input_string.replace(k, contents, -1)

    for k, v in zip(kfr, vfr):
        input_string = input_string.replace(k, v, -1)

    return input_string

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

    output_string = do_replacements(input_fp, kfr, vfr, kfi, vfi)

    output_fh = open(output_fp, 'w')
    output_fh.write(output_string)
    output_fh.close()
