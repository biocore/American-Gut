#!/usr/bin/env python

from __future__ import division

__author__ = "Adam Robbins-Pianka"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Adam Robbins-Pianka"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Adam Robbins-Pianka"
__email__ = "adam.robbinspianka@colorado.edu"

class ArgumentError(Exception):
    """Generic error that can be raised if there is a problem detected with
    one or more of the arguments passed by the user.
    """
    pass

def do_replacements(input_string, kfr, vfr, kfi, vfi):
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
    for k, v in zip(kfi, vfi):
        contents = open(v, 'U').read()
        input_string = input_string.replace(k, contents, -1)

    for k, v in zip(kfr, vfr):
        input_string = input_string.replace(k, v, -1)

    return input_string

