#!/usr/bin/env python

from __future__ import division

__author__ = "Adam Robbins-Pianka"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Adam Robbins-Pianka"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Adam Robbins-Pianka"
__email__ = "adam.robbinspianka@colorado.edu"

class MappingFileError(Exception):
    """Generic error type that can be raised if an error is found in a
    mapping file.
    """
    pass

# These are the only countries that are in the American Gut samples
# for this release of module 2
COUNTRY_ADJECTIVES = {
    'United States of America': 'American',
    'Canada': 'Canadian',
    'Spain': 'Spanish',
    'Norway': 'Norwegian'
}
