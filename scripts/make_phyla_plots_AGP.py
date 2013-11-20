#!/usr/bin/env python

from __future__ import division
from os.path import isfile, exists
from os.path import join as pjoin
from os import mkdir
from argparse import ArgumentParser
from americangut.make_phyla_plots import (map_to_2D_dict,
                              summarize_human_taxa,
                              plot_american_gut,
                              load_category_files)

__author__ = "Justine Debelius"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Justine Debelius"]
__license__ = "BSD" 
__version__ = "unversioned"
__maintainer__ = "Justine Debelius"
__email__ = "Justine.Debelius@colorado.edu"


def main(otu_table, mapping_data, categories, output_dir, \
    samples_to_plot = None, legend = False, xaxis = True):
    """Creates stacked bar plots for an otu table
    INPUTS:
        otu_table -- an open OTU table

        mapping_data -- a tab delimited string containing the mapping data 
                    passed from the mapping file.

        categories -- a dictionary keying a mapping category to the 
                    corresponding sample IDs and taxonomy for a collapsed 
                    biom table

        output_dir -- the location of the directory where output files should be
                    saved. If this directory does not exist, it will be created.

        samples_to_plot -- a list of sample ids to plot. If no value is passed, 
                    then all samples in the biom table are analyzed.

    OUTPUTS:
        A pdf of stacked taxonomy will be generated for each sample and saved 
        in the output directory. These will follow the file name format 
        Figure_4_<SAMPLEID>.pdf
    """
    # Sets constants
    LEVEL = 2
    FILEPREFIX = 'Figure_4_'
    MICHAEL_POLLAN = '000007108.1075657'
    NUM_TAXA = 9
    NUM_CATS_TO_PLOT = 7
    
    # Loads the mapping file
    map_dict = map_to_2D_dict(mapping_data)
    
    (common_taxa, whole_sample_ids, whole_summary) = \
        summarize_human_taxa(otu_table, LEVEL)

    # Converts final taxa to a clean list
    common_phyla = []
    for taxon in common_taxa: 
        common_phyla.append(taxon[1].strip(' p__').strip('[').strip(']'))
    common_taxa = common_phyla
   
    # Checks that the correct sample ids are plotted
    if samples_to_plot == None:
        sample_ids = whole_sample_ids
    else:
        sample_ids = samples_to_plot

    # Identifies Michael Pollan's pre-ABX sample
    mp_sample_pos = whole_sample_ids.index(MICHAEL_POLLAN)
    mp_sample_taxa = whole_summary[:,mp_sample_pos]

    # Loads the category dictionary
    categories = load_category_files(category_fp, LEVEL)

    # Generates a figure for each sample
    for idx, sample_id in enumerate(whole_sample_ids):
        if sample_id in sample_ids:
            # Preallocates a numpy array for the plotting data
            tax_array = zeros((NUM_TAXA, NUM_CATS_TO_PLOT))        
            meta_data = map_dict[sample_id] 
            cat_list = ['You', 'Average', 'Similar Diet', ' Similar BMI', 
                        'Same Gender', 'Similar Age', 
                        'Michael Pollan', '']

            #cat_list.append('Your Fecal Sample')
            #cat_list.append('Average Fecal Samples')
        
            tax_array[:,0] = whole_summary[:,idx]
            tax_array[:,1] = mean(whole_summary, 1)
        
            cat_watch = 2
            # Identifies the appropriate metadata categories
            for cat in categories:                      
                # Pulls metadata for the sample and category
                mapping_key = meta_data[cat]
                # Pulls taxonomic summary and group descriptions for the category
                tax_summary = categories[cat]['Taxa Summary']
                group_descriptions = categories[cat]['Groups']               
                # Amends plotting tables
                try:
                    mapping_col = group_descriptions.index(mapping_key)
                except:
                    raise ValueError, 'The %s cannot be found in %s.' \
                    % (mapping_key, cat)
                tax_array[:,cat_watch] = tax_summary[:,mapping_col]

                cat_watch = cat_watch + 1

            tax_array[:,-1] = mp_sample_taxa
            # Plots the data
            filename = pjoin(output_dir, '%s%s.pdf' \
                % (FILEPREFIX, sample_id))
            plot_american_gut(tax_array, filename)

# Sets up the command line interface
## This uses argparse instead of optparse since optparse is being phased out 
## with new versions of Python.

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
    # Sets the plotting level at phylum
    LEVEL = 2

    args = parser.parse_args()

    # Checks the biom table is sane
    if not args.input:
        parser.error("An input BIOM table is required.")
    elif not isfile(args.input):
        parse.error('The supplied biom table does not exist in the path.')
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
        category_fp = dict([c.strip().split(':') \
            for c in args.categories.split(',')])
        categories = load_category_files(category_fp, LEVEL)
    
    # Deals with the sample list
    if args.samples_to_plot:
        samples = args.samples_to_plot
        samples = samples.split(',')
    else:
        samples = None

    main(otu_table, mapping, output_dir = output_dir, \
        categories = categories, samples_to_plot = samples)
