#!/usr/bin/env python

from __future__ import division
from os import mkdir
from os.path import isfile, exists, join as pjoin
from numpy import array, vstack
from argparse import ArgumentParser
from biom.parse import parse_biom_table
from matplotlib.font_manager import FontProperties
from make_phyla_plots import (translate_colors,
                              calculate_dimensions_rectangle,
                              render_single_pie)
from generate_otu_signifigance_tables import(calculate_abundance,
                                             clean_greengenes_string)
from taxtree import (build_tree_from_taxontable,
                     sample_rare_unique)

__author__ = "Justine Debelius"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Justine Debelius", 'Daniel McDonald']
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Justine Debelius"
__email__ = "Justine.Debelius@colorado.edu"


def main(tax_table, output_dir, samples_to_analyze=None):
    """Generates pie chart of the most abundant twelve taxa in the sample
    INPUTS:
        otu_table -- a biom formatted taxonomy table at the desired level of
        resolution

        output_dir -- the location of the directory where output files should
                    be stored.

        samples_to_analyze -- a list of sample ids to plot. If no value is
                    passed, then all samples in the biom table are analyzed.

    OUTPUTS:
        A pdf of the piechart summarizing the most abundant taxa will be
        generated and saved to the output directory. These will follow the
        naming convention PIECHART_<SAMPLEID>.pdf.
    """

    # Creates the text around hte file name
    FILENAME_BEFORE = 'piechart_'
    FILENAME_AFTER = '.pdf'

    # Handles string cleaning
    RENDER = 'LATEX'
    UNCLASSIFIED = False

    # Sets up the rare threshhold for
    RARE_THRESH = 0.0
    SUM_MIN = 1

    # Sets up axis parameters
    AXIS_LENGTH = 7.25
    AXIS_BORDER = 0.01
    AXIS_TITLE = 0
    AXIS_LEGEND = 6.253
    # Modifies the axis limits
    AX_LIMS = [-1.05, 1.05]
    # Sets up constants for getting the colormap and plotting
    MAP_NAME = 'BrBG'
    NUM_SHOW = 12
    OTHER_COLOR = array([[85/255, 85/255, 85/255]])
    # Sets up plotting parameters
    FIG_LEGEND = True
    FIG_COLOR_EDGE = False
    FIG_LEG_FRAME = False
    FIG_LEG_OFFSET = [1.9, 0.5]
    # Sets up the the legend font
    LEG_FONT = FontProperties()
    LEG_FONT.set_size(30)
    LEG_FONT.set_family('sans-serif')
    # Sets the general font properties
    use_latex = True
    rc_font_family = 'sans-serif'
    rc_font = ['Helvetica', 'Arial']

    # Sets up the colormap
    colormap = translate_colors((NUM_SHOW-1), MAP_NAME)
    colormap = vstack((colormap, OTHER_COLOR))

     # Sets up plotting constants
    (axis_dims, fig_dims) = calculate_dimensions_rectangle(
        axis_width=AXIS_LENGTH, axis_height=AXIS_LENGTH, border=AXIS_BORDER,
        title=AXIS_TITLE, legend=AXIS_LEGEND)

    # Walks over a taxa tree and prioritizes based on taxonomy
    (tree, all_taxa) = build_tree_from_taxontable(tax_table)

    # Sets up samples for which tables are being generated
    if samples_to_analyze is not None:
        samples_to_test = samples_to_analyze
    else:
        samples_to_test = all_taxa.keys()

    # Checks the samples exist
    if samples_to_test:
        samples_to_test = set(samples_to_test)
        tmp = {k: v for k, v in all_taxa.items() if k in samples_to_test}
        all_taxa = tmp
        if not samples_to_test:
            raise ValueError("No samples!")

    # Walks over the table
    filt_fun = lambda v, i, md: v.sum() > 0
    for samp, filtered_table, rare, unique in sample_rare_unique(tree,
                                                                 tax_table,
                                                                 all_taxa,
                                                                 RARE_THRESH):
        # abund_fun = lambda v, i, md: i in all_taxa[samp]
        filtered_table = tax_table.filterObservations(filt_fun)
        sample_data = filtered_table.sampleData(samp)
        taxa = filtered_table.ObservationIds

        # Calculates abundance and limits to the top n samples.
        abund_rank = calculate_abundance(sample=sample_data,
                                         taxa=taxa,
                                         sum_min=SUM_MIN)

        abund_rank = abund_rank[:(NUM_SHOW-1)]

        # Cleans the greengenes strings and adds an "Other" Category for
        # missing taxa
        [sample_tax, sample_freq] = [list(a) for a in zip(*abund_rank)]
        clean_tax = [clean_greengenes_string(tax, RENDER,
                                             unclassified=UNCLASSIFIED)
                     for tax in sample_tax]
        clean_tax.append('Other')
        print samp
        print clean_tax
        sample_freq.append(1-sum(sample_freq))

        # Sets up the sample filename
        filename = pjoin(output_dir, '%s%s%s' % (FILENAME_BEFORE, samp,
                                                 FILENAME_AFTER))


        # Creates the pie chart
        render_single_pie(data_vec=sample_freq,
                          group_names=clean_tax,
                          axis_dims=axis_dims,
                          fig_dims=fig_dims,
                          file_out=filename,
                          legend=FIG_LEGEND,
                          colors=colormap,
                          show_edge=FIG_COLOR_EDGE,
                          legend_frame=FIG_LEG_FRAME,
                          rc_font=rc_font,
                          legend_offset=FIG_LEG_OFFSET,
                          rc_fam=rc_font_family,
                          legend_font=LEG_FONT,
                          use_latex=use_latex,
                          x_lims=AX_LIMS,
                          y_lims=AX_LIMS)


# Sets up command line parsing
parser = ArgumentParser(description='Creates a pie chart of the most abundant'
                        ' taxa in a supplied table.')

parser.add_argument('-i', '--input',
                    help='Path to taxonomy biom table. [REQUIRED]')
parser.add_argument('-o', '--output',
                    help='Path to output directory. [REQUIRED]')
parser.add_argument('-s', '--samples',
                    default=None,
                    help='Sample IDs to be analyzed. If no value is '
                    'supplied, all samples in the taxonomy file will be'
                    ' analyzed.')

if __name__ == '__main__':
    args = parser.parse_args()

    # Checks the tax table file is sane and loads it.
    if not args.input:
        parser.error('An input taxonomy table is required.')
    elif not isfile(args.input):
        parser.error('The supplied taxonomy file does not exist in the path.')
    else:
        tax_table = parse_biom_table(open(args.input, 'U'))

    # Checks the output directory is sane
    if not args.output:
        parser.error('An output directory must be supplied.')
    elif not exists(args.output):
        mkdir(args.output)

    output_dir = args.output

    # Parses the sample IDs
    if args.samples:
        samples_to_analyze = args.samples.split(',')
    else:
        samples_to_analyze = None

    main(tax_table=tax_table,
         output_dir=output_dir,
         samples_to_analyze=samples_to_analyze)
