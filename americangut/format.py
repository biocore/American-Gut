#!/usr/bin/env python

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Daniel McDonald"]
__license__ = "GPLv2"
__version__ = "unversioned"
__maintainer__ = "Yoshiki Vazquez Baeza"
__email__ = "yoshiki.vazquezbaeza@colorado.edu"

from re import compile, findall, escape

def format_print_for_magnified_sample(sample_id, per_sample_file_string,
        global_file_string):
    """ """

    escaped_sample_id = escape(sample_id)

    re = compile('\\<path\ id\=\"'+escaped_sample_id+'\"\ d\=\\"M\\ ([A-Za-z0-9\\.\\ \\-\\,]+)\\" style\\=\\"([A-Za-z0-9\\.\\ \\-\\,\\(\\)\\:\\;]+)\\"\\>\\<\\/path\\>')

    big_sphere_contents = findall(re, per_sample_file_string)
    small_sphere_contents = findall(re, global_file_string)

    # this indicates an internal inconsistency so let the user know
    if big_sphere_contents == None or small_sphere_contents == None:
        return None

    temp = global_file_string

    matchable = '<path id="%s"' % sample_id
    matchable = matchable + ' d="M %s" style="%s"></path>'

    # iter over each of the matches
    for index in range(0, min([len(small_sphere_contents), len(big_sphere_contents)])):
        temp = temp.replace(matchable % small_sphere_contents[index], matchable % big_sphere_contents[index])

    return temp
