#!/usr/bin/env python

"""Subset a mapping file based on Sample Ids in a table

1) verify the sample ids in the table are a subset of the samples in the
metadata
2) slice the metadata as to make a 100% overlap between the samples in the
table and the samples in the metadata
"""

from americangut.util import parse_mapping_file, verify_subset, \
    slice_mapping_file
from biom import load_table
from sys import argv, exit


__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Daniel McDonald"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"


if __name__ == '__main__':
    if len(argv) != 4:
        print "usage: python %s mappingfile table output" % argv[0]
        exit(1)

    header, mapping_file = parse_mapping_file(open(argv[1]))

    table = load_table(argv[2])

    if not verify_subset(table, mapping_file):
        print "****"
        print set([i[0] for i in mapping_file]) - set(table.ids())
        print set(table.ids()) - set([i[0] for i in mapping_file])
        raise ValueError("The table is not a subset of the mapping file!")

    sliced = slice_mapping_file(table, mapping_file)

    f = open(argv[3], 'w')
    f.write(header)
    f.write('\n')
    f.write('\n'.join(sliced))
    f.write('\n')
    f.close()
