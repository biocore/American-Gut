#!/usr/bin/env python

from os.path import join
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--mapping_file', help=('The mapping file. Must '
    'have columns TITLE and COUNTRY.'), required=True)

class MappingFileError(Exception):
    """Generic error type that can be raised if an error is found in a
    mapping file.
    """
    pass

def parse_mapping_file(mapping_file):
    """mapping_file is an open file handle to the mapping file

    Returns a dict of dicts. Outer keys are sample IDs, inner keys are column
    headers.
    """
    old_pos = mapping_file.tell()
    mapping_file.seek(0)
    headers = mapping_file.readline().strip().split('\t')[1:]
    m = {}
    for line in mapping_file:
        elements = line.strip().split('\t')
        m[elements[0]] = dict(zip(headers, elements[1:]))

    mapping_file.seek(old_pos)
    return m

# These are the only countries that are in the American Gut samples
# for this release of module 2
COUNTRY_ADJECTIVES = {
    'GAZ:United States of America': 'American',
    'GAZ:Canada': 'Canadian',
    'GAZ:Spain': 'Spanish',
    'GAZ:Norway': 'Norwegian'
}

if __name__ == '__main__':
    args = parser.parse_args()
    input_mapping_file = args.mapping_file

    # parse the mapping file
    m = parse_mapping_file(open(input_mapping_file))

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

        title_file = open(join(sample_id, 'title.txt'), 'w')
        if data['COUNTRY'] in COUNTRY_ADJECTIVES:
            title = ' %s ' % COUNTRY_ADJECTIVES[data['COUNTRY']]
        else:
            title = ' '

        title_file.write(title)
        title_file.close()
