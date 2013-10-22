#!/usr/bin/env python
# File created on 04 Oct 2013
from __future__ import division

__author__ = "Yoshiki Vazquez Baeza"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Yoshiki Vazquez Baeza", "Adam Robbins-Pianka"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Yoshiki Vazquez Baeza"
__email__ = "yoshiki.vazquezbaeza@colorado.edu"
__status__ = "Development"


from shutil import copy
from os.path import join, exists
from subprocess import Popen, PIPE
from os import listdir, makedirs, remove
from optparse import OptionParser, OptionGroup
from americangut.format import format_print_for_magnified_sample

def main():
    usage = "usage: %prog -i svg_files --prefix figure_1 --sample_id 000000"+\
        "001.314159 -o results/000000001.314159/"
    parser = OptionParser(usage=usage)
    options = OptionGroup(parser, "options")
    options.add_option("-i", "--input_dir", type="string", help="input "
        "directory containing the SVG formatted files")
    options.add_option("--prefix", type="string", help="prefix filename for the"
        " SVG files in the folder i. e. figure_1")
    options.add_option("--sample_id", type="string", help="sample name to be "
        "processed and converted to PDF")
    options.add_option("-o", "--output_dir", type="string", help="output "
        "directory where you want to store the PDF formatted files")
    parser.add_option_group(options)
    opts, args = parser.parse_args()

    input_directory = opts.input_dir
    file_prefix = opts.prefix
    output_directory = opts.output_dir
    sample_id = opts.sample_id

    if exists(input_directory) == False:
        parser.error("The input directory doesn't exist")

    all_the_files = [element for element in listdir(input_directory)
        if element.startswith(file_prefix)]
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

    sample_name = file_prefix+'.'+sample_id+'_huge'

    # extract the contents of this unique sample SVG file
    fd = open(join(input_directory, sample_name), 'U')
    per_sample_file = ''.join(fd.readlines())
    fd.close()

    try:
        temp = format_print_for_magnified_sample(sample_id, per_sample_file,
            global_file_contents)
    except RuntimeError:
        parser.error('Problem found with sample %s ' % sample_id)

    transformed_svg_file = join(output_directory, sample_name+'.svg')

    fd_out = open(transformed_svg_file, 'w')
    fd_out.write(temp)
    fd_out.close()

    # in case something went wrong just make sure that the file exists
    assert exists(transformed_svg_file), "Something went wrong the file does "+\
        "exist (%s)." % transformed_svg_file

    # -f takes any input file and -A will convert it into a PDF
    inkscape_command = 'inkscape -f %s -A %s' % (transformed_svg_file,
        transformed_svg_file[:-3]+'pdf')

    # based on pyqi/util.pyqi_system_call
    process = Popen(inkscape_command, shell=True, stdout=PIPE, stderr=PIPE)
    _, _ = process.communicate()
    if process.returncode != 0:
        parser.error('Could not convert the file from SVG to PDF, exit status '
            'code is %d', process.returncode)

    # remove the svg file as it is just transient
    try:
        remove(transformed_svg_file)
    except OSError:
        parser.error('Could not delete the SVG file')


if __name__ == "__main__":
    main()
 