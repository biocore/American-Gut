#!/usr/bin/env python

from __future__ import division
from matplotlib import use
use('agg')
from os.path import isfile, exists
from os.path import join as pjoin
from os import mkdir
from biom.parse import parse_biom_table, table_factory
from numpy import (array, zeros, mean, arange, shape, ones, around, 
                   vstack, amax, sum as n_sum)
import matplotlib.pyplot as plt
from matplotlib import font_manager, rc 
from matplotlib.transforms import Bbox
from argparse import ArgumentParser
from operator import itemgetter

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

def idenitfy_most_common_categories(biom_table, level, 
    metadata_category = 'taxonomy', limit_mode = 'COMPOSITE', limit = 1.0):
    """Identifies the most common taxa in a population using variable limits

    This method uses a composite score to account for both the average frequency
    and the percentage of the population in which the sample is present. This 
    hopes to be more robust than either the average or fraction present alone.

    INPUTS:
        biom_table -- a sparse biom table to be summarized

        level -- an integer corresponding to the taxonomic level (or other meta 
                    data category level) at which data should be summarized.
                    The argument will only take a single integer.

        metadata_category -- a description of the metadata category over which 
                    the data will be summarized.

        limit_mode -- a string describing how the scoring limit which should be 
                    used. Options are 'COMPOSITE', 'AVERAGE', 'COUNTS' or 
                    'NUMBER', 'NONE'. 
                    COMPOSITE keeps samples whose composite scores greater than 
                        the limit. (Default limit is 1)
                    AVERAGE keeps samples whose average are greater than the 
                        specified limit. (Default limit is 0.01)
                    COUNTS keep groups which are represented in at least the 
                        specified fraction of sample (default limit is 0.01)
                    NONE sorts groups alphabetically and keeps all.                 

        limit -- the numeric lower limit for common taxa

        logging -- a binary value specifying if a table of the frequency, 
                    composite value, and metadata category should be printed and
                    retained.

        log_file -- the filename where the log file should be saved.

    OUTPUTS
        sample_ids -- a list of the sample ids contained in the matrix 
                    (column header)
        
        common_categories -- a list of the common categories in the (row header)

        category_summary -- a numpy array of the category values.

        scores -- a sorted list of all the taxa and their scores

    """ 
    # Sets up positions
    COMPOSITE_CONSTANT = 10000

    if limit_mode == 'COMPOSITE':
        score_position = 3
        score_reverse = True
        second_score = 1
    elif limit_mode == 'AVERAGE':
        score_position = 1
        score_reverse = True
        second_score = 3
    elif limit_mode == 'COUNTS':
        score_position = 2  
        score_reverse = True  
        second_score = 3
    elif limit_mode == 'NONE':
        score_position = 0
        score_reverse = False
        second_score = 3
    else:
        raise ValueError ("limit_mode is not a supported option. \n"\
            "Options for limit_mode are 'COMPOSITE', 'AVERAGE', 'COUNTS', or"\
            " 'NONE'.")

    # Gets the Sample IDs
    sample_ids = list(biom_table.SampleIds)
    num_samples = len(sample_ids)

    # Prealocates output objects
    scoring_all = []
    common_categories = []
    table_row = []
    other = []


    # Normalizes the data by the number of observations so relative frequencies 
    # are used.
    biom_table_norm = biom_table.normObservationBySample()
    
    # Collapses the OTUs into category summaries using the correct levels
    for (bin, table) in biom_table_norm.binObservationsByMetadata(lambda x: 
            x[metadata_category][:level]):      
        # Pulls out the sample data for the group
        group_value = array(table.sum('sample'))
        group_binary = group_value > 0

        # Calculates presence scores
        average_freq = round(mean(group_value),4)
        fraction_pres = round(float(sum(group_binary))/float(num_samples),4)
        composite = round(average_freq*fraction_pres*COMPOSITE_CONSTANT,2)
        score_row = [bin, average_freq, fraction_pres, composite]

        # Adds the scores to the watch matrix
        if fraction_pres > 0:
            scoring_all.append(score_row)

    # Sorts based on scoring method
    scores_all = sorted(sorted(scoring_all, 
                               key = itemgetter(second_score), 
                               reverse = True),
                        key = itemgetter(score_position),
                        reverse = score_reverse)

    # Identifies rows which meet the scoring criteria
    scores = []
    watch = True
    for idx, score in enumerate(scores_all):
        scores.append([score[0], round(score[1],4), round(score[2], 4), 
            round(score[3], 2)])

        if score[score_position] > limit:
            common_categories.append(score[0])                         

    # Raises an error if necessary
    if len(common_categories) == 0:
        raise ValueError, ('Limit too high! No common categories could be'\
            ' identified.')

    common_categories.append(u'Other')   

    # Returns the values
    return common_categories, scores

def summarize_common_categories(biom_table, level, common_categories, 
    metadata_category = 'taxonomy'):
    """Determines the frequency of common categories present in a biom table

    INPUTS:
        biom_table -- a sparse biom table to be evaluated

        level -- an integer corresponding to the taxonomic level (or other meta 
                    data category level) at which data should be summarized.
                    The argument will only take a single integer.

        common_categories -- a list of values which define the common categories.

        metadata_category -- a description of the metadata category over which 
                    the data will be summarized.

    OUTPUTS:
        sample_ids -- a list of the sample ids found in the OTU table.

        cat_summary -- a numpy array describing the frequency of the common 
                    categories with an additional frequencies collapsed into 
                    the "other" category.
    """
    num_cats = len(common_categories)

    sample_ids = list(biom_table.SampleIds)
    num_samples = len(sample_ids)

    # Normalizes the biom table
    biom_norm = biom_table.normObservationBySample()

    # Prealocates numpy objects (because that makes life fun!). tax_other is 
    # set up as a row array because this is ultimately summed
    cat_summary = zeros([num_cats, num_samples])
    cat_other = zeros([1, num_samples])

    # Collapses the OTU table using the category at the correct level
    chunked = [(bin, table) for bin, table in 
               biom_norm.binObservationsByMetadata(\
               lambda x: x[metadata_category][:level])]

    for idx, (bin, table) in enumerate(chunked):
        if bin in common_categories:
            cat_summary[common_categories.index(bin)] = table.sum('sample') 
        else:
            cat_other = vstack((cat_other, table.sum('sample')))

    cat_summary = vstack((cat_summary, sum(cat_other)))

    return sample_ids, cat_summary

def calculate_dimensions_square(axis_size = 4, border = 0.1, title = 1, 
    legend = 2, xlab = 0, ylab = 0, unit = 'inches'):
    """Determines the appriate axis and figure dimensions for square axis.

    INPUTS:
        axis_size -- a number specifying the side length the axis. DEFAULT: 4

        border -- the fraction of the axis size to be added as white space 
                    around the edge of the axis.

        title -- the height to add to the top of the figure for a title. This 
                    is separate from the border, which is added by default. 
                    DEFAULT: 1

        legend -- the width of the legend to the be added to the figure. 
                    DEFAULT: 2

        xlab -- the height to be added for labels along the x axis. DEFAULT: 0

        ylab -- the width to be added for labels along the y axis. DEFAULT: 0

        unit -- a string ('inches' or 'cm'), specifying the unit to be used 
                    in image generation. DEFAULT: inches

    OUTPUTS:
        axis_dimensions -- a Bbox class describing the axis position in the 
            figure

        figure_dimensions -- a 2 element tuple giving the width and height of 
            the figure in inches
    """

    # Specifies a value for converting between units
    if unit == 'cm':
        conversion = float(1)/float(2.54)
    elif unit == 'inches':
        conversion = 1.0
    else:
        raise ValueError, 'unit must be "inches" or "cm".'

    # Determines the figure dimensions
    fig_width = float(axis_size*(1+border*2)+legend+ylab)
    fig_height = float(axis_size*(1+border*2)+title+xlab)

    figure_dimensions = (fig_width*conversion, fig_height*conversion)

    # Determines the axis bounds
    axis_left = float(border*axis_size+ylab)/fig_width
    axis_right = float((1+border)*axis_size+ylab)/fig_width
    axis_bottom = float(border*axis_size+xlab)/fig_height
    axis_top = float((1+border)*axis_size+xlab)/fig_height

    axis_dimensions = Bbox(array([[axis_left, axis_bottom],
                                  [axis_right, axis_top]]))

    return axis_dimensions, figure_dimensions

def calculate_dimensions_bar(num_bars, bar_width = 0.5, axis_height = 3, 
    border = 0.1, title = 1, legend = 2, xlab = 0, ylab = 0, unit = 'inches'):
    """Determines the axis and figure dimensions for a bar chart.

    INPUTS:
        num_bars 

        border -- the size of white space to be added around the axis. 
                    DEFAULT: 0.1"

        title -- the height to add to the top of the figure for a title. This 
                    is separate from the border, which is added by default. 
                    DEFAULT: 1

        legend -- the width of the legend to the be added to the figure. 
                    DEFAULT: 2

        xlab -- the height to be added for labels along the x axis. DEFAULT: 0

        ylab -- the width to be added for labels along the y axis. DEFAULT: 0

        unit -- a string ('inches' or 'cm'), specifying the unit to be used 
                    in image generation. DEFAULT: inches

    OUTPUTS:
        axis_dimensions -- a Bbox class describing the axis position in the 
            figure

        figure_dimensions -- a 2 element tuple giving the width and height of 
            the figure in inches
    """

    # Preforms some sanity checks. 
    if num_bars < 1:
        raise ValueError ('There must be at least one group to plot.\nnum_bars'\
            ' must be a whole number greater than 1.')
    elif not round(num_bars) == num_bars:
        raise ValueError ('There cannot be partial categories. \nnum_bars must'\
            ' be a whole number greater than 1.')
 

   # Specifies a value for converting between units
    if unit == 'cm':
        conversion = float(1)/float(2.54)
    elif unit == 'inches':
        conversion = 1.0
    else:
        raise ValueError, 'unit must be "inches" or "cm".'

    # Determines the figure width
    axis_width = float(num_bars*bar_width)
    figure_width = float(border*2+ylab+axis_width+legend)
    figure_height = float(border*2+xlab+axis_height+title)

    figure_dimensions = (figure_width, figure_height)

    axis_left = float(ylab+border)/figure_width
    axis_right = float(ylab+border+axis_width)/figure_width
    axis_bottom = float(xlab+border)/figure_height
    axis_top = float(xlab+border+axis_height)/figure_height

    axis_dimensions = Bbox(array([[axis_left, axis_bottom],
                                  [axis_right, axis_top]]))

    return axis_dimensions, figure_dimensions

def render_single_pie(cats_vec, cat_names, axis_dims, fig_dims, 
    file_out = 'piechart', filetype = 'PDF', colormap = [[1, 1, 1]], 
    legend = True, title = None, labels = None, label_distance = 1.1, 
    start_angle = 90, radius = 1, fontsize = 15):
    """Creates a pie chart summarizing the category data

    INPUTS:
        cats_vec -- a vector which sums to 1 describing the fraction of the chart
                    represented by each group in cat_names.

        cat_names -- a list of the groups in cat_vec. (i.e. Firmictues, 
                    Bacteriodetes, Proteobacteria for a category of Taxonomy).

        file_out -- a string giving the file path where the pie plot should be 
                    saved.

        colormap -- an n x 3 numpy array giving the desired colormap Default is 
                    to color each wedge white.

        axis_dims -- a 2 x 2 numpy array giving the fraction of the figure 
                    which should bound the axis. (row 1: [left, bottom], row 2: 
                    [right, top] as a percentage of the figure space)

        legend -- a binary value indicating whether or not a legend should be 
                    shown. This must be accounted for in the figure dimensions. 
                    Default is to show the legend.

        filetype -- a string specifying the type of file to be generated. 
                    Options include 'PNG', 'PDF', etc. Default is PDF. 

        title -- a string giving a title to append to the figure. This must be 
                    accounted for in the figure dimensions.
                    Default is no title.

        labels -- a list of values to be used as labels, if they are desired. 
                    Default is no labels (None)

        label_distance -- the radial distance from the chart where the label 
                    should be placed. Default is 1.1

        start_angle -- the angle from the x axis which the plot should start 
                    plotting. Default is 90 (start at the y-axis and plot 
                    clockwise)

    OUTPUTS:
        The rendered figure is saved in the at the file_out location.
    """

    # Creates the figure
    sample_figure = plt.figure(1, fig_dims)

    # Plots the data clockwise
    pie_plot = pie(x = cat_vec[-1:], 
                   labels = labels, 
                   labeldistance = label_distance, 
                   shadow = False,
                   startangle = start_angle,
                   color = colormap)

    # Sets up the axis dimensions
    ax1 = plt.gca()
    ax1.set_position(axis_dims)

    # Adds the legend if necessary
    if legend:
        plt.figlegend(pie_plot, cat_names, 'right')

    # Adds the title if desired
    if title.__class__ == str:
        plt.title(title)

    # Saves the output figure
    plt.savefig(file_out, format = filetype)

def render_bar_char(cat_table, cat_names, sample_names, axis_dims, 
    fig_dims, file_out = 'barchart', filetype = 'PDF', colormap = [[1, 1, 1]],
    legend = True, x_axis = True, x_tick_offset = 0.6, 
    bar_width = 0.8, x_min = -0.5, x_tick_interval = 1.0, y_lims = [0, 1], 
    y_tick_interval = 0.2, tick_font_size = 15, label_font_size = 20, 
    font_angle = 45, font_alignment = 'right'):
    """Creates a stacked barchart using the data in the category table.

    INPUTS:

    OUTPUT:
        The rendered figure is saved in the at the file_out location.       
    """

    # Preforms a sanity checks that the provided data is good
    (table_height, table_width) = cat_table.shape
    num_cats = len(cat_names)
    num_samples = len(sample_names)

    if not table_height == num_cats:
        raise ValueError ('The number of provided categories differ.')
    elif not table_width == num_samples:
        raise ValueError ('The number of samples differ.')

    # Sets up the x ticks 
    x_tick = arrange(0, num_samples)
    x_max = x_min + no_samples
    bar_left = x_tick - bar_width*0.5

    # Creates the figure
    fig = plt.figure(1, fig_dims)

    
    # Plots the data
    patches_watch = []
    for plot_count, category in enumerate(cat_table):
        bottom_bar = sum(category_table[:plot_count,:])
        faces = plt.bar(left = bar_left, 
                        height = category, 
                        width = bar_width, 
                        bottom = bottom_bar,
                        color = colormap[plot_count,:],
                        edgecolor = colormap[plot_count,:])
        patches_watch.append(faces[0])

    # Sets up the axis dimesnions and limits.
    ax1 = plt.gca()
    x1.set_position(axis_dims)

     # The y-direction is reversed so the labels are in the same order as the 
    # colors in the legend
    plt.axis([x_min, x_max, y_lims[1], y_lims[0]])

     # Sets y axis labels
    y_tick_labels = (arange(y_max + y_tick_interval, y_min, \
        -y_tick_interval) - y_tick_interval)*100
    y_tick_labels[-1] = 0

    # Converts y label to text
    y_text_labels = [str(e) for e in y_tick_labels]

    y_tick_labels = ax1.set_yticklabels(y_text_labels, 
                                        size = tick_font_size)
    y_axis_label = ax1.set_ylabel('Frequency', 
                                  size = label_font_size)

    # Sets up the x-axis labels
    x_text_labels = sample_labels[:] # copy
    x_text_labels.insert(0, '') # insert at the head of the list
    
    x_tick_labels = ax1.set_xticklabels(x_text_labels, 
                                        size = tick_font_size,
                                        rotation = font_angle, 
                                        horizontalalignment = font_alignment)
    if legend:
        plt.figlegend(patches_watch, category_names, 'right')

    plt.savefig(file_out, format = file_format)

   

# def plot_stacked_phyla(taxonomy_table, taxonomy_headers, sample_labels, \
#   file_out, sample_ids=None, legend = True, x_axis=True):
#     """Creates a stacked taxonomy plot at the phylum level

#     INPUTS:
#         taxonomy_table -- a numpy array with sample information in the columns 
#                         and phylum frequency information in the rows

#         taxonomy_headers -- a  of the phyla which to the rows in the 
#                         taxonomy_table

#         sample_labels -- a list of the sample labels which correspond to the 
#                         rows in the taxonomy_table

#         file_out -- a string describing the filename for the output filename

#     OUTPUT:
#         The rendered figure is saved as a pdf in the at the file_out location.
#     """
#     # Colorbrewer python package. Can be loaded

#     # Colormap is taken from the colorbrewer    
#     COLORMAP = array([[0.8353, 0.2421, 0.3098],
#                       [0.9569, 0.4275, 0.2627],
#                       [0.9922, 0.6824, 0.3804],
#                       [0.9961, 0.8784, 0.5351],
#                       [0.9020, 0.9608, 0.5961],
#                       [0.6706, 0.8667, 0.6431],
#                       [0.4000, 0.7608, 0.6471],
#                       [0.1961, 0.5333, 0.7412],
#                       [0.3333, 0.3333, 0.3333]])

#     X_TICK_OFFSET = 0.6

#     BAR_WIDTH = 0.8

#     # Paramatizable the constants and the options for shape/size. User Passed size parameters

#     if legend and x_axis:
#         figure_dimensions = (8, 5)
#         axis_dimensions = Bbox(array([[0.2000, 0.3000],
#                                       [0.7000, 0.9000]]))
#     elif legend:
#         figure_dimensions = (8, 3.75)
#         axis_dimensions = Bbox(array([[0.2000, 0.1000],
#                                       [0.7000, 0.9000]]))
#     elif x_axis:
#         figure_dimensions = (6.2222, 5)
#         axis_dimensions = Bbox(array([[0.2000, 0.3000],
#                                       [0.9000, 0.9000]]))
#         #bah
#     else:
#         figure_dimensions = (6.22222, 3.75)
#         axis_dimensions = Bbox(array([[0.2000, 0.1000],
#                                       [0.9000, 0.9000]]))


#     X_MIN = -0.5
#     X_TICK_INTERVAL = 1.0

#     Y_MIN = 0
#     Y_MAX = 1.0
#     Y_TICK_INTERVAL = 0.2

#     TICK_FONT_SIZE = 15
#     LABEL_FONT_SIZE = 20

#     [no_phyla, no_samples] = taxonomy_table.shape

#     x_tick = arange(0,no_samples)
#     x_max = X_MIN+no_samples

#     sample_figure = plt.figure(1, figure_dimensions)


#     patches_watch = []

#     # Plots the data
#     for plot_count, phyla in enumerate(taxonomy_table):
#         already_added_index = arange(plot_count)
#         bottom_bar = sum(taxonomy_table[already_added_index,:])
#         faces = plt.bar(x_tick-BAR_WIDTH/2, phyla, BAR_WIDTH, \
#                         bottom = bottom_bar, color = COLORMAP[plot_count,:], \
#                         edgecolor=COLORMAP[plot_count,:])
#         patches_watch.append(faces[1])

#     # Sets up axis dimensions and limits
#     ax1 = plt.gca()
#     ax1.set_position(axis_dimensions)
    
#     # The y-direction is reversed so the labels are in the same order as the 
#     # colors in the legend
#     plt.axis([X_MIN, x_max, Y_MAX, Y_MIN])

#     # Sets y axis labels
#     y_tick_labels = (arange(Y_MAX + Y_TICK_INTERVAL, Y_MIN, \
#         -Y_TICK_INTERVAL) - Y_TICK_INTERVAL)*100
#     y_tick_labels[-1] = 0

#     # Converts y label to text
#     y_text_labels = [str(e) for e in y_tick_labels]

#     y_tick_labels = ax1.set_yticklabels(y_text_labels, size = TICK_FONT_SIZE)
#     y_axis_label = ax1.set_ylabel('Frequency', size = LABEL_FONT_SIZE)

#     # Sets up the x-axis labels
#     x_text_labels = sample_labels[:] # copy
#     x_text_labels.insert(0, '') # insert at the head of the list
    
#     x_tick_labels = ax1.set_xticklabels(x_text_labels, size = TICK_FONT_SIZE,\
#         rotation = 45, horizontalalignment = 'right')

#     # Adds the legend
#     if legend:
#         plt.figlegend(patches_watch, taxonomy_headers, 'right')


#     plt.savefig(file_out, format = 'pdf')

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

