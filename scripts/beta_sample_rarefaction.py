#!/usr/bin/env python
# File created on 01 Mar 2012
from __future__ import division

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Development"

from matplotlib import use
use('Agg')
from qiime.util import parse_command_line_parameters, make_option
from qiime.parse import parse_distmat
from random import sample
from numpy import mean, std, inf
from matplotlib.pyplot import (figure, subplot, grid, title, axis, savefig,
                               ylabel, xlabel)

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("", "", "")]
script_info['output_description'] = ""
script_info['required_options'] = [
    make_option('-i', '--input_path', type="existing_filepaths",
                help='the input distance matrix file(s)'),
    make_option('-l', '--labels', type=str,
                help='legend labels for the input files'),
    make_option('-t', '--title', type=str,
                help='plot title'),
    make_option('-y', '--ylabel', type=str,
                help='ylabel'),
    make_option('-o', '--output_path', type="new_filepath",
                help='the output file [default: %default]',
                default='plot.pdf'),
    make_option('-n', '--iterations', type=int,
                help="Number of iterations: %default", default=100),
]
script_info['optional_options'] = [
    make_option('--y_max', type='float',
                help='max y value [default: %default]',
                default=None),
]
script_info['version'] = __version__


def main():
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)

    input_path = opts.input_path
    output_path = opts.output_path
    iterations = opts.iterations
    verbose = opts.verbose
    y_max = opts.y_max
    labels = opts.labels.split(',')

    results = {}
    for input_file in input_path:
        if verbose:
            print input_file

        # Reading OTU/biom table
        samples, distmat = parse_distmat(open(input_file, 'U'))
        possible_samples = range(len(distmat[0]))

        result_iteration = []
        for iteration in range(iterations):
            iter_vals = []
            for n in possible_samples:
                if n < 1:
                    continue
                curr_samples = sample(possible_samples, n+1)

                curr_vals = []
                for curr_i, i in enumerate(curr_samples):
                    for j in curr_samples[curr_i+1:]:
                        curr_vals.append(distmat[i][j])
                iter_vals.append(min(curr_vals))
            result_iteration.append(iter_vals)

        results[input_file] = [mean(result_iteration, axis=0),
                               std(result_iteration, axis=0)]

        if verbose:
            f = open(output_path + '.txt', 'a')
            f.write('\t'.join(map(str,results[input_file][0])))
            f.write('\n')
            f.write('\t'.join(map(str,results[input_file][1])))
            f.write('\n')
            f.close()

    # generating plot, some parts taken from
    # http://stackoverflow.com/questions/4700614
    fig = figure()
    ax = subplot(111)

    max_x, max_y = -inf, -inf
    for i, (label, input_file) in enumerate(zip(labels, input_path)):
        len_x = len(results[input_file][0])
        len_y = max(results[input_file][0])
        if max_x < len_x:
            max_x = len_x
        if max_y < len_y:
            max_y = len_y
        if i % 2 == 0:
            coloring = (215/255.0, 48/255.0, 39/255.0)
        else:
            coloring = (69/255.0, 177/255.0, 180/255.0)

        ax.errorbar(range(1, len_x+1), results[input_file][0],
                    yerr=results[input_file][1], fmt='o', color=coloring,
                    label=label)

    if y_max:
        axis([0, max_x, 0, max_y])
    else:
        axis([0, max_x, 0, y_max])

    # Shink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    title(opts.title)
    xlabel('Samples')
    ylabel(opts.ylabel)
    grid(True)
    savefig(output_path)


if __name__ == "__main__":
    main()
