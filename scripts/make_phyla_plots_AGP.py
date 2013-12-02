#!/usr/bin/env python

from __future__ import division
from os.path import isfile, exists
from os.path import join as pjoin
from os import mkdir
from numpy import zeros, mean, array
from argparse import ArgumentParser
from biom.parse import parse_biom_table
from make_phyla_plots import (map_to_2D_dict,
                              render_barchart,
                              summarize_common_categories,
                              load_category_files,
                              parse_category_files)

__author__ = "Justine Debelius"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Justine Debelius", "Daniel McDonald"]
__license__ = "BSD" 
__version__ = "unversioned"
__maintainer__ = "Justine Debelius"
__email__ = "Justine.Debelius@colorado.edu"


def main(otu_table, mapping_data, cat_tables, output_dir, \
    samples_to_plot = None, legend = False, xaxis = True):
    """Creates stacked bar plots for an otu table
    INPUTS:
        otu_table -- an open OTU table

        mapping_data -- a tab delimited string containing the mapping data 
                    passed from the mapping file.

        categories -- a dictionary keying a mapping category to the 
                    corresponding biom table

        output_dir -- the location of the directory where output files should be
                    saved. If this directory does not exist, it will be created.

        samples_to_plot -- a list of sample ids to plot. If no value is passed, 
                    then all samples in the biom table are analyzed.

    OUTPUTS:
        A pdf of stacked taxonomy will be generated for each sample and saved 
        in the output directory. These will follow the file name format 
        Figure_4_<SAMPLEID>.pdf
    """
    # Sets constants for analyzing the data
    LEVEL = 2
    CATEGORY = 'taxonomy'
    NUM_TAXA = 9
    NUM_CATS_TO_PLOT = 7
    MICHAEL_POLLAN = '000007108.1075657'

    # Sets up file name constants
    FILEPREFIX = 'Figure_4_'
    FILE_END = '.pdf'

    # Common taxa are designated before processing to remain constant.
    COMMON_TAXA = [(u'k__Bacteria', u' p__Firmicutes'),
                   (u'k__Bacteria', u' p__Bacteroidetes'),
                   (u'k__Bacteria', u' p__Proteobacteria'),
                   (u'k__Bacteria', u' p__Actinobacteria'),
                   (u'k__Bacteria', u' p__Verrucomicrobia'),
                   (u'k__Bacteria', u' p__Tenericutes'),
                   (u'k__Bacteria', u' p__Cyanobacteria'),
                   (u'k__Bacteria', u' p__Fusobacteria')]
    
    # Names categories being plotted
    CAT_LIST = ['You', 'Average', 'Similar Diet', ' Similar BMI', 
                        'Same Gender', 'Similar Age', 
                        'Michael Pollan']


    # Sets up plotting constants
    COLORMAP = array([[0.8353, 0.2421, 0.3098],
                      [0.9569, 0.4275, 0.2627],
                      [0.9922, 0.6824, 0.3804],
                      [0.9961, 0.8784, 0.5351],
                      [0.9020, 0.9608, 0.5961],
                      [0.6706, 0.8667, 0.6431],
                      [0.4000, 0.7608, 0.6471],
                      [0.1961, 0.5333, 0.7412],
                      [0.3333, 0.3333, 0.3333]])

    FIG_DIMS = (4.44444, 3.33333)
    AXIS_DIMS = array([[0.05, 0.05],
                       [0.95, 0.95]])
    
    
    # Creates the mapping file
    map_dict = map_to_2D_dict(mapping_data)

    # Generates the category file dictionary
    categories = parse_category_files(raw_tables = cat_tables, 
                                      common_groups = COMMON_TAXA, 
                                      level = LEVEL, 
                                      metadata = CATEGORY)

    (whole_sample_ids, whole_summary, new_common_taxa) = \
    summarize_common_categories(biom_table = otu_table, 
                                level = LEVEL, 
                                common_categories = COMMON_TAXA[:8], 
                                metadata_category = CATEGORY)

    # Converts final taxa to a clean list
    common_phyla = []
    for taxon in new_common_taxa: 
        common_phyla.append(taxon[1].strip(' p__').strip('[').strip(']'))
    new_common_taxa = common_phyla
   
    # Checks that the correct sample ids are plotted
    if samples_to_plot is None:
        sample_ids = whole_sample_ids
    else:
        sample_ids = samples_to_plot

    # Identifies Michael Pollan's pre-ABX sample
    mp_sample_pos = whole_sample_ids.index(MICHAEL_POLLAN)
    mp_sample_taxa = whole_summary[:,mp_sample_pos]

    # Generates a figure for each sample
    for idx, sample_id in enumerate(whole_sample_ids):
        if sample_id in sample_ids:

            # Preallocates a numpy array for the plotting data
            tax_array = zeros((NUM_TAXA, NUM_CATS_TO_PLOT))
            meta_data = map_dict[sample_id]            

            tax_array[:,0] = whole_summary[:,idx]
            tax_array[:,1] = mean(whole_summary, 1)
        
            cat_watch = 2
            # Identifies the appropriate metadata categories
            for cat in categories:
                # Pulls metadata for the sample and category
                mapping_key = meta_data[cat]
                # Pulls taxonomic summary and group descriptions for the category
                tax_summary = categories[cat]['Summary']
                group_descriptions = categories[cat]['Groups']          
                # Appends plotting tables
                try:
                    mapping_col = group_descriptions.index(mapping_key)
                except:
                    raise ValueError, 'The %s cannot be found in %s.' \
                    % (mapping_key, cat)

                tax_array[:,cat_watch] = tax_summary[:,mapping_col]

                cat_watch = cat_watch + 1

            tax_array[:,-1] = mp_sample_taxa

            # Plots the data
            filename = pjoin(output_dir, '%s%s%s' \
                % (FILEPREFIX, sample_id, FILE_END))
            
            render_barchart(data_table = tax_array, x_axis = False,
                group_names = new_common_taxa, sample_names = CAT_LIST,
                axis_dims = AXIS_DIMS, fig_dims = FIG_DIMS, file_out = filename,
                colors = COLORMAP, show_edge = False, legend = False,  
                y_axis = False)

# Sets up the command line interface

# Creates the parser object
parser = ArgumentParser(description = 'Creates stacked bar plots for an OTU'\
                        ' table.')
parser.add_argument('-i', '--input', required = True, \
                    help = 'OTU table file path [REQUIRED]')
parser.add_argument('-m', '--mapping', required = True, \
                    help = 'Mapping file path [REQUIRED]')
parser.add_argument('-o', '--output', required = True, \
                    help = 'Path to the output directory [REQUIRED]'),
parser.add_argument('-c', '--categories', \
                    help = 'Category associations with a collapsed OTU file '\
                    'path. The string should be associated with a colon, for'\
                    ' example, "SEX:sex.biom,DIET_TYPE:diet.biom"'),
parser.add_argument('-s', '--samples_to_plot', default = None, \
                    help = 'Sample IDs you wish to plot. If no value is '\
                    'specified, all samples are plotted.')

if __name__ == '__main__':
    args = parser.parse_args()

    # Checks the biom table is sane
    if not args.input:
        parser.error("An input BIOM table is required.")
    elif not isfile(args.input):
        parser.error('The supplied biom table does not exist in the path.')
    else:
        otu_table = parse_biom_table(open(args.input, 'U'))     

    # Checks the mapping file is sane
    if not args.mapping:
        parser.error("An input mapping file is required.")
    elif not isfile(args.mapping):
        parser.error('he supplied file does not exist in the path')        
    else:
        mapping = open(args.mapping, 'U')
        
    # Checks the output directory is sane   
    if not args.output:
        parser.error("An output directory is required.")
    elif not exists(args.output):
        mkdir(args.output)
    output_dir = args.output

    # Parses the category argument
    if not args.categories:
        categories = {}     
    else:
        category_fp = dict([c.strip().split(':') for c in \
            args.categories.split(',')])
        categories = load_category_files(category_files = category_fp)
    
    # Deals with the sample list
    if args.samples_to_plot:
        samples = args.samples_to_plot
        samples = samples.split(',')
    else:
        samples = None

    main(otu_table, mapping, output_dir = output_dir, \
        cat_tables = categories, samples_to_plot = samples)


### Commentary on the selection of common taxa:
# Common taxa can be calculated using the function, 
# identify_most_common_categories. When this was run on rounds 1, 2, and 3 of 
# the American Gut for fecal and all sites equally weighted, and for the HMP for
# fecal only and equal weights on the fecal, skin and oral sites. 


