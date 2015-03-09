#!/usr/bin/env python

from sys import argv
from americangut.util import parse_mapping_file
from biom import load_table


__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2015, The American Gut Project"
__credits__ = ["Daniel McDonald"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"


if __name__ == '__main__':
    if len(argv) != 4:
        print "usage: python %s mappingfile table threshold" % argv[0]
        exit(1)

    header, mapping_file = parse_mapping_file(open(argv[1]))

    table = load_table(argv[2])
    threshold = float(argv[3])

    map_ids = {i[0].split('.')[0] for i in mapping_file}
    table_ids = {i.split('.')[0] for i in table.ids()}

    missing = {i for i in (map_ids - table_ids) if 'blank' not in i.lower()}

    counts = table.sum(axis='sample')

    below_threshold = set()
    for id_, count in zip([i.split('.')[0] for i in table.ids()], counts):
        if 'blank' in id_.lower():
            continue

        if count < threshold:
            below_threshold.add(id_)

    if missing:
        print "The following IDs do not have any reads:"
        for id_ in sorted(missing):
            print id_
        print

    if below_threshold:
        print "The following have fewer than %d reads" % threshold
        for id_ in sorted(below_threshold):
            print id_
