#!/usr/bin/env python

from __future__ import division

__author__ = "Sam Way"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Sam Way"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Sam Way"
__email__ = "samuel.way@colorado.edu"

import argparse

import brewer2mpl

from americangut.parse import get_filtered_taxa_summary
from americangut.plots import make_stack_plot, \
    make_pie_chart, make_legend


def interface():
    args = argparse.ArgumentParser()
    args.add_argument('-m', '--mapping-file',
                      help='Mapping file',
                      required=True)
    args.add_argument('-t', '--taxa-file',
                      help='Taxa summary file',
                      required=True)
    args.add_argument('-c', '--metadata-category',
                      help='Metadata category',
                      required=True)
    args.add_argument('-v', '--metadata-value',
                      help='Specific metadata value',
                      default=None)
    args.add_argument('-k', '--key-taxa-file',
                      help='List of taxa to examine')
    args.add_argument('-o', '--output-prefix',
                      help='Output file prefix',
                      default='./out')
    args.add_argument('-f', '--output-type',
                      help='Output file type',
                      default='pdf')
    args.add_argument('-s', '--samples-file',
                      help='Sample ids to be labeled')
    args.add_argument('-l', '--ylabel',
                      help='Y-axis label',
                      default='Phylum')
    args = args.parse_args()
    return args


def get_sample_ids_to_label(samples_file):
    """ Get sample id, label tuples to be highlighted
        Each non-empty, non-comment line of the samples_file should
        contain a tab-separated sample id/label pair

        Inputs:
        samples_file - file containing sample id/label lines

        Outputs:
        sample_label_tuples - a list containing tuples extracted from
                              samples_file.  Tuples are of the form:
                              (sample_id, sample_label)

    """
    sample_label_tuples = []
    lines = open(samples_file, 'U')
    for line in lines:
        if line.startswith('#'):
            continue
        line_pieces = [x.strip() for x in line.split('\t')]
        if len(line_pieces) == 2:
            sample_label_tuples.append(tuple(line_pieces[0:2]))
    lines.close()

    return sample_label_tuples


def get_key_taxa(key_taxa_file):
    """ Load taxa to be plotted from file.
        Each non-comment line of the file is interpretted as a taxonomy label.

        Inputs:
        key_taxa_file - file containing a list of taxonomy labels,
                        one per line.

        Outputs:
        key_taxa - list of taxonomy labels extracted from key_taxa_file

    """
    key_taxa = []
    lines = open(key_taxa_file, 'U')
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        key_taxa.append(line)
    lines.close()
    return key_taxa


if __name__ == "__main__":
    args = interface()

    # Sample ticks to be labeled in the stack plot
    if args.samples_file:
        special_labels = get_sample_ids_to_label(args.samples_file)
    else:
        special_labels = []

    # Specify taxa, otherwise the top N most abundant
    if args.key_taxa_file:
        select_taxa = get_key_taxa(args.key_taxa_file)
    else:
        select_taxa = None

    filtered_sample_ids, taxa_labels, collapsed_taxa_table = \
        get_filtered_taxa_summary(args.mapping_file, args.taxa_file,
                                  args.metadata_category, args.metadata_value,
                                  select_taxa=select_taxa)

    colors = brewer2mpl.get_map('Spectral', 'Diverging',
                                len(taxa_labels)).mpl_colors

    # Create stack plot
    output = args.output_prefix + 'stack.' + args.output_type
    make_stack_plot(collapsed_taxa_table, filtered_sample_ids, args.ylabel,
                    colors, output, sample_ticks=special_labels)

    # Create pie chart
    output = args.output_prefix + 'pie.' + args.output_type
    make_pie_chart(collapsed_taxa_table, colors, output)

    # Create figure legend
    output = args.output_prefix + 'legend.' + args.output_type
    make_legend(taxa_labels, colors, output)
