#!/usr/bin/env python

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Daniel McDonald"]
__license__ = "GPLv2"
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
