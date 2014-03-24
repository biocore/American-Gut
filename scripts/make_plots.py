#!/usr/bin/env python

__author__ = "Sam Way"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Sam Way"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Sam Way"
__email__ = "samuel.way@colorado.edu"

import argparse

import brewer2mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from numpy import array, cumsum, arange

from americangut.parse import parse_mapping_file_to_dict, \
    get_filtered_taxa_summary


def interface():
    args = argparse.ArgumentParser()
    args.add_argument('-m', '--mapping-file',
                      help='Mapping file',
                      required=True)
    args.add_argument('-t', '--taxa-file',
                      help='Taxa summary file',
                      required=True)
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
    args.add_argument('-c', '--metadata-category',
                      help='Metadata category',
                      default='SIMPLE_MATTER')
    args.add_argument('-v', '--metadata-value',
                      help='Specific metadata value',
                      default=None)
    args.add_argument('-l', '--ylabel',
                      help='Y-axis label',
                      default='Phylum')
    args = args.parse_args()
    return args


def make_stack_plot(taxa_table, sample_ids, ylabel,
                    colors, output_file, sample_ticks=None):
    """ Makes a stack plot. Sample ids in sample_ticks
        are labeled in the figure.

        Inputs:
        taxa_table - table of taxa counts by sample ids
        sample_ids - sample ids for the taxonomy table
        ylabel - label for the y axis
        colors - colors to be used for each taxonomy group
        output_file - filename for the output figure
        sample_ticks - list of tuples to be labeled
                       where each tuple is
                       (sample_id, sample_label).
                       NOTE: if sample_id is found in the
                       list of sample_ids, it is
                       labeled in the figure.

        Outputs:
        (None - output figure is written to output_file)

    """

    #N=taxa groups, M=num samples
    N, M = collapsed_taxa_table.shape
    x = arange(M)
    cumulative = cumsum(collapsed_taxa_table, axis=0)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    # Get xticks, if any
    xticks = []
    xtick_labels = []
    if sample_ticks is not None:
        for sample_id, sample_label in sample_ticks:
            if sample_id in sample_ids:
                sample_index = sample_ids.index(sample_id)
                xticks.append(sample_index)
                xtick_labels.append(sample_label)

    ax1.fill_between(x, 1, 1-cumulative[0, :], color=colors[0])
    for i in xrange(1, N):
        ax1.fill_between(x, 1-cumulative[i-1, :],
                         1-cumulative[i, :], color=colors[i])
    yticks = arange(0.0, 1.25, .25)
    ytick_labels = [str(int(ytick*100)) for ytick in yticks]
    ax1.set_xticks(xticks)
    ax1.set_xticklabels(xtick_labels, rotation=40, ha='right')
    ax1.xaxis.set_tick_params(width=1, length=10, pad=7,
                              direction='out', top='off')
    ax1.set_yticks(yticks)
    ax1.set_yticklabels(ytick_labels)
    ax1.yaxis.set_tick_params(width=1, length=5, direction='out', right='off')
    plt.ylabel('%s Frequency (%%)' % (ylabel), fontsize=16)
    plt.ylim([0, 1])
    plt.xlim([1, M])
    plt.subplots_adjust(bottom=0.25)
    plt.savefig(output_file)


def make_pie_chart(table, colors, output_file):
    """ Create a simple pie chart.

        Inputs:
        table - table of counts number of categories by number of samples
        colors - list of matplotlib-compatible colors, equal in length
                 to the number of rows/categories in table
        output_file - filename for the resulting output figure

        Outputs:
        (None - output figure is written to output_file)

     """
    if table.shape[0] != len(colors):
        raise ValueError("Number of colors must equal number of rows!")
    fractions = [100*x for x in table.mean(axis=1)]
    fig = plt.figure()
    wedges, texts = plt.pie(fractions, colors=colors)

    for w in wedges:
        w.set_linewidth(0)

    plt.axis('equal')
    plt.savefig(output_file, bbox_inches='tight')


def make_legend(labels, colors, output_file):
    """ Hack to generate a legend as a separate image.  A dummy plot
        is created in one figure, and its legend gets saved to an
        output file.

        Inputs:
        labels - text labels for the legend
        colors - list of colors for the legend (same length as labels)
        output_file - filename for the resulting output figure

        Outputs:
        (None - output figure is written to output_file)

    """
    if len(labels) != len(colors):
        raise ValueError("Lists of labels and colors "
                         " should have the same length")
    fig = plt.figure()
    font_prop = FontProperties()
    font_prop.set_size('xx-large')
    font_prop.set_family('sans-serif')
    seperate_legend = plt.figure(figsize=(5, 3.5))
    ax = fig.add_subplot(111)
    N = len(labels)
    # Make a dummy pie chart so we can steal its legend
    wedges, texts = ax.pie([100/N]*N, colors=colors)
    seperate_legend.legend(wedges, labels, 'center',
                           prop=font_prop, frameon=False)
    seperate_legend.savefig(output_file)


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
    for line in open(samples_file, 'U'):
        if line.startswith('#'):
            continue
        line_pieces = [x.strip() for x in line.split('\t')]
        if len(line_pieces) == 2:
            sample_label_tuples.append(tuple(line_pieces[0:2]))

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
    for line in open(key_taxa_file, 'U'):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        key_taxa.append(line)
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
