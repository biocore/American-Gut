#!/usr/bin/env python

from __future__ import division

from argparse import ArgumentParser
from os.path import join, isdir
from os import mkdir

from americangut.util import parse_mapping_file_dict
from americangut.create_titles import MappingFileError, COUNTRY_ADJECTIVES

__author__ = "Adam Robbins-Pianka"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Adam Robbins-Pianka"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Adam Robbins-Pianka"
__email__ = "adam.robbinspianka@colorado.edu"

parser = ArgumentParser()
parser.add_argument('-m', '--mapping_file', help='The mapping file. Must '
    'have columns TITLE and COUNTRY.', required=True)
parser.add_argument('-f', '--force', help="Force creation of the sample's "
    "output directory if it does not exist", action="store_true")

if __name__ == '__main__':
    args = parser.parse_args()

    input_mapping_file = args.mapping_file
    force = args.force

    # parse the mapping file
    m = parse_mapping_file_dict(open(input_mapping_file))

    # check if it has the necessary columns
    _, test_data = m.iteritems().next()
    if 'TITLE' not in test_data:
        raise MappingFileError, ("Mapping file must have TITLE column!")
    if 'COUNTRY' not in test_data:
        raise MappingFileError, ("Mapping file must have COUNTRY column!")

    # go through the mapping file, find those that are American Gut samples,
    # and generate a 'title.txt' file that contains the proper text in that
    #sample's directory
    for sample_id, data in m.iteritems():
        if data['TITLE'] != 'American Gut': continue

        if not isdir(sample_id):
            if force:
                mkdir(sample_id)
            else:
                raise IOError, ("Output directory not found for sample %s. "
                    "Pass -f or --force to automatically create sample output "
                    "directories.") % sample_id

        title_file = open(join(sample_id, 'title.txt'), 'w')
        if data['COUNTRY'] in COUNTRY_ADJECTIVES:
            title = ' %s ' % COUNTRY_ADJECTIVES[data['COUNTRY']]
        else:
            title = ' '

        title_file.write(title)
        title_file.close()
