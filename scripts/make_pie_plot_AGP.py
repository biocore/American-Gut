#!/usr/bin/env python

from __future__ import division
from os import mkdir
from os.path import isfile, exists, join as pjoin
from numpy import array, vstack
from argparse import ArgumentParser
from biom.parse import parse_biom_table
from americangut.make_phyla_plots import (translate_colorbrewer,
                                          calculate_dimensions_rectangle,
                                          render_single_pie)
from americangut.generate_otu_signifigance_tables import(calculate_abundance,
                                                        clean_greengenes_string)
from americangut.taxtree import (build_tree_from_taxontable, 
                                 sample_rare_unique)

__author__ = "Justine Debelius"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Justine Debelius", 'Daniel McDonald']
__license__ = "BSD" 
__version__ = "unversioned"
__maintainer__ = "Justine Debelius"
__email__ = "Justine.Debelius@colorado.edu"

def main(tax_table, output_dir, samples_to_plot = None, continue_sparse = False):
    """
    INPUTS:
        otu_table -- a biom formatted taxonomy_table
    """

    # Sets up constants for getting the colormap and plotting
    MAP_NAME = 'Spectral'
    NUM_SHOW = 12
    
    # Creates the text around hte file name
    FILENAME_BEFORE = 'figure_5_'
    FILENAME_AFTER = '.pdf'

    # Handles string cleaning
    RENDER = 'RAW'
    UNCLASSIFIED = False

    # Sets up the rare threshhold for 
    RARE_THRESH = 0.001

    # Sets the "other" color to dark grey
    OTHER_COLOR = array([[85/255, 85/255, 85/255]])

    # Walks over a taxa tree and prioritizes based on taxonomy
    (tree, all_taxa) = build_tree_from_taxontable(tax_table)

    # Sets up samples for which tables are being generated    
    if not samples_to_analyze == None:
        samples_to_test = samples_to_analyze
    else:
        samples_to_test = all_taxa.keys()

    # Checks the samples exist
    if samples_to_test:
        samples_to_test = set(samples_to_test)
        tmp = {k:v for k,v in all_taxa.items() if k in samples_to_test}
        all_taxa = tmp
        if not samples_to_test:
            raise ValueError, "No samples!"

    # Walks over the table
    for samp, filtered_table, rare, unique in sample_rare_unique(tree, 
                                                                 tax_table, 
                                                                 all_taxa, 
                                                                 RARE_THRESH):
        filtered_table = filtered_table.filterObservations(lambda v,i,md:\
            v.sum() > 0)
        sample_data = filtered_table.sampleData(samp)
        taxa = filtered_table.ObservationIds

        # Calculates abundance and limits to the top n samples.
        abund_rank = calculate_abundance(sample = sample_data,
                                         taxa = taxa,
                                         abundance_threshhold = 1)
        abund_rank = abund_rank[:(NUM_SHOW-1)]

        # Cleans the greengenes strings and adds an "Other" Category for missing 
        # taxa
        [sample_tax, sample_freq] = [list(a) for a in zip(*abund_rank)]
        clean_tax = [clean_greengenes_string(tax, RENDER, 
                            unclassified = UNCLASSIFIED)for tax in sample_tax]
        clean_tax.append('Other')
        sample_freq.append(1-sum(sample_freq))

        # Sets up plotting constants
        (axis_dims, fig_dims) = calculate_dimensions_rectangle(axis_width = 6,
                                                           axis_height = 6, 
                                                           border = 0.1, 
                                                           title = 0,
                                                           legend = 6)

        colormap = translate_colorbrewer((NUM_SHOW-1), MAP_NAME)
        colormap = vstack((colormap, OTHER_COLOR))

        # Sets up the sample filename
        filename = pjoin(output_dir, '%s%s%s' % (FILENAME_BEFORE, samp, 
            FILENAME_AFTER))

        # Creates the pie chart
        render_single_pie(cats_vec = sample_freq,
                          cat_names = clean_tax,
                          axis_dims = axis_dims,
                          fig_dims = fig_dims,
                          file_out = filename,
                          legend = True,
                          colors = colormap,
                          show_edge = False,
                          fontsize = 20,
                          legend_offset_x = 1.875)








# Sets up command line parsing
parser = ArgumentParser(description = 'Creates a pie chart of the most abundant taxa in a supplied table.')

parser.add_argument('-i', '--input',
                    help = 'Path to taxonomy biom table. [REQUIRED]')
parser.add_argument('-o', '--output',
                    help = 'Path to output directory. [REQUIRED]')
parser.add_argument('-s', '--samples',
                    default = None,
                    help = 'Sample IDs to be analyzed. If no value is '\
                    'supplied, all samples in the taxonomy file will be'\
                    ' analyzed.')

if __name__ == '__main__':
    args = parser.parse_args()

    # Checks the tax table file is sane and loads it.
    if not args.input:
        parser.error('An input taxonomy table is required.')
    elif not isfile(args.input):
        parser.error('The supplied taxonomy file does not exist in the path.')
    else:
        tax_table = parse_biom_table(open(args.input))

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

    main(tax_table = tax_table, 
         output_dir = output_dir, 
         samples_to_plot = samples_to_analyze)



