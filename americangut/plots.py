#!/usr/bin/env python

from __future__ import division

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from numpy import cumsum, arange

__author__ = "Sam Way"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Sam Way"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Sam Way"
__email__ = "samuel.way@colorado.edu"


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

    # N=taxa groups, M=num samples
    N, M = taxa_table.shape
    x = arange(M)
    cumulative = cumsum(taxa_table, axis=0)
    fig = plt.figure()
    ax1 = fig.gca()

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
    plt.ylabel('%s Frequency (%%)' % ylabel, fontsize=16)
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
    plt.figure()
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
