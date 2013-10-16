#!/usr/bin/env python
# File created on 04 Oct 2013
from __future__ import division

__author__ = "Yoshiki Vazquez Baeza"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Yoshiki Vazquez Baeza"]
__license__ = "GPLv2"
__version__ = "unversioned"
__maintainer__ = "Yoshiki Vazquez Baeza"
__email__ = "yoshiki.vazquezbaeza@colorado.edu"
__status__ = "Development"


from shutil import copy
from os import listdir, makedirs
from os.path import join, exists
from optparse import OptionParser, OptionGroup
from americangut.format import format_print_for_magnified_sample

def main():
    parser = OptionParser()
    options = OptionGroup(parser, "options")
    options.add_option("-i", "--input_dir", type="string", help="input "
        "directory containing the SVG formatted files")
    options.add_option("--prefix", type="string", help="prefix filename for the"
        " SVG files in the folder i. e. figure_1")
    options.add_option("-o", "--output_dir", type="string", help="output "
        "directory where you want to store the PDF formatted files")
    parser.add_option_group(options)
    opts, args = parser.parse_args()

    input_directory = opts.input_dir
    file_prefix = opts.prefix
    output_directory = opts.output_dir

    if exists(input_directory) == False:
        parser.error("The input directory doesn't exist")

    all_the_files = [element for element in listdir(input_directory) if element.startswith(file_prefix)]
    global_file_path = file_prefix + '.GLOBAL'

    if not len(all_the_files):
        parser.error("no files were found")

    if global_file_path not in all_the_files:
        parser.error("cannot continue ... no global file found")

    try:
        makedirs(output_directory)
    except:
        pass

    fd = open(join(input_directory, global_file_path), 'U')
    global_file_contents = ''.join(fd.readlines())
    fd.close()

    for element in all_the_files:
        # avoid processing the global files
        if element == global_file_path:
            continue

        # avoid processing files for label files
        if element == file_prefix+'.labels':
            copy(join(input_directory, element), join(output_directory, element+'.svg'))
            continue

        # build the sample id from the file name
        sample_id = element.replace(file_prefix+'.', '').replace('_huge', '')

        # extract the contents of this unique sample SVG file
        fd = open(join(input_directory, element), 'U')
        per_sample_file = ''.join(fd.readlines())
        fd.close()

        try:
            temp = format_print_for_magnified_sample(sample_id, per_sample_file,
            global_file_contents)
        except RuntimeError:
            print 'Problem found with %s, skipping it for now' % sample_id

        fd_out = open(join(output_directory, element+'.svg'), 'w')
        fd_out.write(temp)
        fd_out.close()



if __name__ == "__main__":
    main()
 