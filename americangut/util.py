#!/usr/bin/env python

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Daniel McDonald"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"

def pick_rarifaction_level(id_, lookups):
    """Determine which lookup has the appropriate key
    
    id_ is a barcode, e.g., '000001000'
    lookups is a list of tuples, e.g., [('10k',{'000001000':'000001000.123'})]

    The order of the lookups matters. The first lookup found with the key will
    be returned.

    None is returned if the key is not found
    """
    for name, lookup in lookups:
        if id_ in lookup:
            return name
    return None

def parse_mapping_file(open_file):
    """return (header, [(sample_id, all_other_fields)])
    
    """
    header = open_file.readline().strip()
    res = []

    for l in open_file:
        res.append(l.strip().split('\t',1))
    
    return (header, res)

def verify_subset(table, mapping):
    """Returns True/False if the table is a subset"""
    ids = set([i[0] for i in mapping])
    t_ids = set(table.SampleIds)

    return t_ids.issubset(ids)

def slice_mapping_file(table, mapping):
    """Returns a new mapping corresponding to just the ids in the table"""
    t_ids = set(table.SampleIds)
    res = []
    
    for id_, l in mapping:
        if id_ in t_ids:
            res.append('\t'.join([id_, l]))

    return res
