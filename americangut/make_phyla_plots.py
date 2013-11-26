#!/usr/bin/env python

from __future__ import division
from matplotlib import use
use('agg')
from os.path import isfile, exists
from os.path import join as pjoin
from os import mkdir
from biom.parse import parse_biom_table, table_factory
from numpy import array, zeros, arange, shape, ones, around
import matplotlib.pyplot as plt
from matplotlib import font_manager, rc 
from matplotlib.transforms import Bbox
from argparse import ArgumentParser

__author__ = "Justine Debelius"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Justine Debelius"]
__license__ = "BSD" 
__version__ = "unversioned"
__maintainer__ = "Justine Debelius"
__email__ = "Justine.Debelius@colorado.edu"

def map_to_2D_dict(mapping_data):
    """ Converts mapping file to 2D dictionary

    INPUT:
        mapping_data -- a tab delimited string from the opened mapping file

    OUTPUT:
        D2 -- a two dimensional dictionary where each sample ID is keyed to a 
                    dictionary of meta data containing headers and observations.
    """
    lines = [l.strip().split('\t') for l in mapping_data]
    header = lines[0]
    D2 = {}
    for l in lines[1:]:
        inner = {k:v for k,v in zip(header, l)}
        sample_id = inner['#SampleID']
        D2[sample_id] = inner

    return D2

def most_common_taxa_gg_13_5():
    """Identifies the most common taxa at a given phylogenetic level

    INPUT:
        level -- a number (2), specifying the level of taxa

    OUTPUTS:
        common_taxa -- a list of common taxa

    To identify the most common phyla in human samples (hard coded version), the
    average frequency for the phyla (+/- standard deviation), the fraction of 
    samples containing the phyla and the product of the two values (multiplied 
    by a constant) were calculated for fecal samples and the equally weighted 
    average of fecal, oral and skin samples. For the HMP, fecal samples were 
    identified using the HMPbodysupersite metadata field. Airway and Urovaginal 
    samples in the HMP were not considered (Tables 1-4).

    Samples were ranked on the composite score. The composite scores by phylum 
    and study are summarized in table 5. Eight phyla were identified using this 
    table which appear in every sample, have a composite score of atleast 0.50 
    in all studies and combined represent on average 99% of the otus in a given 
    sample.

    Table 1. HMP Fecal Samples
    ---------------------------------------------------------------------------
        Phylum          Average Frequency     % Present  Composite score
    ---------------------------------------------------------------------------
    Bacteroidetes       0.6680 +/- 0.1856         100       6679.94
    Firmicutes          0.2944 +/- 0.1812         100       2943.61
    Proteobacteria      0.0262 +/- 0.0395          93        242.29
    Tenericutes         0.0060 +/- 0.0160          32         19.37
    Actinobacteria      0.0018 +/- 0.0050          44          7.95
    Verrucomicrobia     0.0022 +/- 0.0099          21          4.77 
    Cyanobacteria       0.0006 +/- 0.0030          31          0.50
    Fusobacteria        0.0007 +/- 0.0049           5          0.33
    Lentisphaerae       0.0001 +/- 0.0009           5          0.06
    ---------------------------------------------------------------------------

    Table 2. HMP Equally Weighted Average
    ---------------------------------------------------------------------------
       Phylum          Average Frequency     % Present  Composite score
    ---------------------------------------------------------------------------
    Firmicutes          0.3312 +/- 0.3472           100       3311.63
    Bacteroidetes       0.2876 +/- 0.2195            94       2711.00
    Actinobacteria      0.2531 +/- 0.2783            80       2025.61
    Proteobacteria      0.0968 +/- 0.1793            95        924.57
    Fusobacteria        0.0228 +/- 0.0736            47        106.16
    Cyanobacteria       0.0028 +/- 0.0503            17          4.61
    Tenericutes         0.0020 +/- 0.0163            18          3.98
    Spirochaetes        0.0019 +/- 0.0191            15          2.95
    Verrucomicrobia     0.0008 +/- 0.0100             8          0.58
    Synergistetes       0.0002 +/- 0.0039             4          0.09
    TM7                 0.0001 +/- 0.0016             6          0.07
    GN02                0.0001 +/- 0.0017             4          0.05
    Lentisphaerae       0.0000 +/- 0.0009             2          0.01
    Chloroflexi         0.0000 +/- 0.0010             1          0.01
    Thermi              0.0000 +/- 0.0009             2          0.01
    ---------------------------------------------------------------------------

    Table 3. AGP Feces
    ---------------------------------------------------------------------------
         Phylum        Average Frequency     % Present   Composite score       
    ---------------------------------------------------------------------------
    Firmicutes          0.4364 +/- 0.2083           100         4363.90
    Bacteroidetes       0.4218 +/- 0.2123           100         4217.64
    Proteobacteria      0.0838 +/- 0.1361            98          822.30
    Actinobacteria      0.0217 +/- 0.0565            87          189.11
    Verrucomicrobia     0.0283 +/- 0.0717            58          165.31
    Tenericutes         0.0060 +/- 0.0180            44           26.18
    Cyanobacteria       0.0009 +/- 0.0058            23            2.10
    Euryarchaeota       0.0005 +/- 0.0036            10            0.47
    Fusobacteria        0.0005 +/- 0.0040             7            0.34
    Lentisphaerae       0.0001 +/- 0.0007             6            0.06
    Synergistetes       0.0001 +/- 0.0005             4            0.03
    ---------------------------------------------------------------------------

    Table 4. AGP Equally Weighted
    ---------------------------------------------------------------------------
         Phylum        Average Frequency     % Present   Composite score       
    ---------------------------------------------------------------------------
    Firmicutes          0.4267 +/- 0.3281           100         4267.45
    Bacteroidetes       0.1899 +/- 0.2593           100         1893.06
    Proteobacteria      0.1696 +/- 0.2813            99         1685.63
    Actinobacteria      0.1544 +/- 0.2100            96         1478.62
    Cyanobacteria       0.0325 +/- 0.1241            52          170.22
    Fusobacteria        0.0120 +/- 0.0418            60           72.24
    Verrucomicrobia     0.0098 +/- 0.0717            31           30.59
    Tenericutes         0.0022 +/- 0.0180            22            4.77
    SR1                 0.0005 +/- 0.0045            13            0.69
    Acidobacteria       0.0005 +/- 0.0037            12            0.55
    Euryarchaeota       0.0007 +/- 0.0077             7            0.46
    TM7                 0.0002 +/- 0.0012            12            0.26
    Spirochaetes        0.0002 +/- 0.0014             9            0.16
    Planctomycetes      0.0002 +/- 0.0015             8            0.14
    Thermi              0.0001 +/- 0.0009             6            0.07
    Chloroflexi         0.0001 +/- 0.0007             7            0.07
    Gemmatimonadetes    0.0001 +/- 0.0008             5            0.04
    Synergistetes       0.0001 +/- 0.0011             4            0.04
    Lentisphaerae       0.0000 +/- 0.0007             2            0.01
    ---------------------------------------------------------------------------

    Table 5. Composite Scores for fecal and equally weighted samples by project
    ---------------------------------------------------------------
    Phylum          AGP fecal   AGP Equal   HMP fecal   HMP Equal
    ---------------------------------------------------------------
    Bacteroidetes   4217.64     1893.06     6679.94     2711.00
    Firmicutes      4363.90     4267.45     2943.61     3311.63
    Actinobacteria   189.11     1478.62        7.95     2025.61
    Proteobacteria   822.30     1685.63      242.29      924.57
    Verrucomicrobia  165.31       30.59        4.77        0.58
    Fusobacteria       0.34       72.24        0.33      106.16
    Cyanobacteria      2.10      170.22        0.50        4.61
    Tenericutes       26.18        4.77       19.37        3.98
    Spirochaetes       0           0.16        0           2.95
    Euryarchaeota      0.47        0.46        0           0
    TM7                0           0.26        0           0.07
    Synergistetes      0.03        0.04        0           0.09
    Lentisphaerae      0.06        0.01        0.06        0.01
    Thermi             0           0.07        0           0.07
    Planctomycetes     0           0.14        0           0
    Chloroflexi        0           0.07        0           0.01
    GN02               0           0           0.05        0
    Gemmatimonadetes   0           0.04        0           0
    --------------------------------------------------------------
    * Composite scores of less than 0.01 are listed as 0 here.

    
    """
    ### Current plan is to implement a method for identifying the most common 
    ### taxa in a data set based on the composite score. This will probably be 
    ### thresholded either based on the composite score (default 1.00), or at an
    ### average frequency (default likely at 95%). However, this still needs to 
    ### written. 


    common_taxa = [(u'k__Bacteria', u' p__Firmicutes'),
                   (u'k__Bacteria', u' p__Bacteroidetes'),
                   (u'k__Bacteria', u' p__Proteobacteria'),
                   (u'k__Bacteria', u' p__Actinobacteria'),
                   (u'k__Bacteria', u' p__Verrucomicrobia'),
                   (u'k__Bacteria', u' p__Tenericutes'),
                   (u'k__Bacteria', u' p__Cyanobacteria'),
                   (u'k__Bacteria', u' p__Fusobacteria'),
                   (u'k__Bacteria', u' p__Other')]

    return common_taxa

def summarize_human_taxa(otu_table, level, metadata_category = 'taxonomy'):
    """Determines the frequency of major human taxa in an OTU at a preset level

    INPUTS:
        otu_table -- a sparse biom table to be summarized

        level -- an integer corresponding to the taxonomic level (or other meta 
                    data category level) at which data should be summarized.
                    The argument will only take a single integer.

        metadata_category -- the feature of the table which should be used for 
                    summarization.
        

    OUTPUTS:
        common_taxa -- a list of common taxonomy at the specified level.

        sample_ids -- a list of the sample ids found in the OTU table. These 

        tax_summary -- a numpy array 
    """
    common_taxa = most_common_taxa_gg_13_5()

    num_taxa = len(common_taxa)

    # Gets the sample ids from the otu table
    sample_ids = list(otu_table.SampleIds)
    
    # Determines the total number of counts in the table
    table_total = otu_table.sum('sample')

    # Collapses OTUs into taxonomic summaries from using the correct levels
    chunked = [(bin, table) for bin, table in \
        otu_table.binObservationsByMetadata(lambda x: x[metadata_category][:level])]
   
    tax_summary = zeros([num_taxa, len(otu_table.SampleIds)])
    tax_other = zeros((len(chunked), len(otu_table.SampleIds)))
    for idx, (bin, table) in enumerate(chunked):
        if bin in common_taxa:
            tax_summary[common_taxa.index(bin)] = table.sum('sample')
        else:
            tax_other[idx] = table.sum('sample')

    tax_summary[num_taxa - 1] = (tax_other.sum(0))

    tax_summary = around(tax_summary/table_total,4)

    return common_taxa, sample_ids, tax_summary

def plot_stacked_phyla(taxonomy_table, taxonomy_headers, sample_labels, \
  file_out, sample_ids=None, legend = True, x_axis=True):
    """Creates a stacked taxonomy plot at the phylum level

    INPUTS:
        taxonomy_table -- a numpy array with sample information in the columns 
                        and phylum frequency information in the rows

        taxonomy_headers -- a  of the phyla which to the rows in the 
                        taxonomy_table

        sample_labels -- a list of the sample labels which correspond to the 
                        rows in the taxonomy_table

        file_out -- a string describing the filename for the output filename

    OUTPUT:
        The rendered figure is saved as a pdf in the at the file_out location.
    """
    # Colorbrewer python package. Can be loaded

    # Colormap is taken from the colorbrewer    
    COLORMAP = array([[0.8353, 0.2421, 0.3098],
                      [0.9569, 0.4275, 0.2627],
                      [0.9922, 0.6824, 0.3804],
                      [0.9961, 0.8784, 0.5351],
                      [0.9020, 0.9608, 0.5961],
                      [0.6706, 0.8667, 0.6431],
                      [0.4000, 0.7608, 0.6471],
                      [0.1961, 0.5333, 0.7412],
                      [0.3333, 0.3333, 0.3333]])

    X_TICK_OFFSET = 0.6

    BAR_WIDTH = 0.8

    # Paramatizable the constants and the options for shape/size. User Passed size parameters

    if legend and x_axis:
        figure_dimensions = (8, 5)
        axis_dimensions = Bbox(array([[0.2000, 0.3000],
                                      [0.7000, 0.9000]]))
    elif legend:
        figure_dimensions = (8, 3.75)
        axis_dimensions = Bbox(array([[0.2000, 0.1000],
                                      [0.7000, 0.9000]]))
    elif x_axis:
        figure_dimensions = (6.2222, 5)
        axis_dimensions = Bbox(array([[0.2000, 0.3000],
                                      [0.9000, 0.9000]]))
        #bah
    else:
        figure_dimensions = (6.22222, 3.75)
        axis_dimensions = Bbox(array([[0.2000, 0.1000],
                                      [0.9000, 0.9000]]))


    X_MIN = -0.5
    X_TICK_INTERVAL = 1.0

    Y_MIN = 0
    Y_MAX = 1.0
    Y_TICK_INTERVAL = 0.2

    TICK_FONT_SIZE = 15
    LABEL_FONT_SIZE = 20

    [no_phyla, no_samples] = taxonomy_table.shape

    x_tick = arange(0,no_samples)
    x_max = X_MIN+no_samples

    sample_figure = plt.figure(1, figure_dimensions)


    patches_watch = []

    # Plots the data
    for plot_count, phyla in enumerate(taxonomy_table):
        already_added_index = arange(plot_count)
        bottom_bar = sum(taxonomy_table[already_added_index,:])
        faces = plt.bar(x_tick-BAR_WIDTH/2, phyla, BAR_WIDTH, \
                        bottom = bottom_bar, color = COLORMAP[plot_count,:], \
                        edgecolor=COLORMAP[plot_count,:])
        patches_watch.append(faces[1])

    # Sets up axis dimensions and limits
    ax1 = plt.gca()
    ax1.set_position(axis_dimensions)
    
    # The y-direction is reversed so the labels are in the same order as the 
    # colors in the legend
    plt.axis([X_MIN, x_max, Y_MAX, Y_MIN])

    # Sets y axis labels
    y_tick_labels = (arange(Y_MAX + Y_TICK_INTERVAL, Y_MIN, \
        -Y_TICK_INTERVAL) - Y_TICK_INTERVAL)*100
    y_tick_labels[-1] = 0

    # Converts y label to text
    y_text_labels = [str(e) for e in y_tick_labels]

    y_tick_labels = ax1.set_yticklabels(y_text_labels, size = TICK_FONT_SIZE)
    y_axis_label = ax1.set_ylabel('Frequency', size = LABEL_FONT_SIZE)

    # Sets up the x-axis labels
    x_text_labels = sample_labels[:] # copy
    x_text_labels.insert(0, '') # insert at the head of the list
    
    x_tick_labels = ax1.set_xticklabels(x_text_labels, size = TICK_FONT_SIZE,\
        rotation = 45, horizontalalignment = 'right')

    # Adds the legend
    if legend:
        plt.figlegend(patches_watch, taxonomy_headers, 'right')


    plt.savefig(file_out, format = 'pdf')

def plot_american_gut(taxonomy_table, file_out):
    """Makes American Gut specific plots
    """
    # Colormap is taken from the colorbrewer    
    COLORMAP = array([[0.8353, 0.2421, 0.3098],
                      [0.9569, 0.4275, 0.2627],
                      [0.9922, 0.6824, 0.3804],
                      [0.9961, 0.8784, 0.5351],
                      [0.9020, 0.9608, 0.5961],
                      [0.6706, 0.8667, 0.6431],
                      [0.4000, 0.7608, 0.6471],
                      [0.1961, 0.5333, 0.7412],
                      [0.3333, 0.3333, 0.3333]])

    BAR_WIDTH = 0.8

    X_MIN = -0.5

    Y_MIN = 0
    Y_MAX = 1.0

    figure_dimensions = (4.44444, 3.33333)
    axis_dimensions = Bbox(array([[0.05, 0.05],[0.95,0.95]]))

    [no_phyla, no_samples] = taxonomy_table.shape

    x_tick = arange(0,no_samples)
    x_max = X_MIN+no_samples

    sample_figure = plt.figure(1, figure_dimensions)

    patches_watch = []

    # Plots the data
    for plot_count, phyla in enumerate(taxonomy_table):
        already_added_index = arange(plot_count)
        bottom_bar = sum(taxonomy_table[already_added_index,:])
        faces = plt.bar(x_tick-BAR_WIDTH/2, phyla, BAR_WIDTH, \
                        bottom = bottom_bar, color = COLORMAP[plot_count,:], \
                        edgecolor=COLORMAP[plot_count,:])
        patches_watch.append(faces[1])

     # Sets up axis dimensions and limits
    ax1 = plt.gca()
    ax1.set_position(axis_dimensions)

    # The y-direction is reversed so the labels are in the same order as the 
    # colors in the legend
    plt.axis([X_MIN, x_max, Y_MAX, Y_MIN])

    # Removes any axis labels
    ax1.set_yticklabels('')
    ax1.set_xticklabels('')


    plt.savefig(file_out, format = 'pdf')

def load_category_files(category_files,level):
    """
    INPUTS:
         category_files -- a dictionary that associates the mapping category 
                    (key) with the file path to the otu_table summarizing that 
                    category.
    OUTPUTS:
        category_tables -- a dictionary that associates the mapping category 
                    with the summarized OTU tables for that category.
    """

    category_tables = {}
    watch_count = 0
    watch_list = []

    for category, category_file in category_files.iteritems():
        
        if isfile(category_file) == False:
            watch_list.append('The summarized OTU table file cannot be found '
                              'for %s. \n%s is not in the file path.' 
                              % (category, category_file))            
            watch_count = watch_count + 1
        else:
            cat_table = parse_biom_table(open(category_file, 'U'))
            (common_taxa, cat_ids, cat_summary)  = \
              summarize_human_taxa(cat_table,level)
            category_tables[category] = {'Groups': cat_ids, \
                                         'Taxa Summary': cat_summary}

    if watch_count > 0:        
        print 'The following category files could not be found: \n%s' \
        % '\n'.join(watch_list)
    if watch_count == len(category_files):       
        raise ValueError, 'No files could be found for any of the supplied '\
            'categories. \n%s' % '\n'.join(watch_list)

    return category_tables

