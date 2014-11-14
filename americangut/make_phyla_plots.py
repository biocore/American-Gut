#!/usr/bin/env python

from __future__ import division
from matplotlib import use
use('agg')
from os.path import isfile
from biom.parse import parse_biom_table
from biom.util import biom_open
from numpy import (array, zeros, mean, ones, vstack, arange, ndarray)
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from matplotlib.font_manager import FontProperties
from matplotlib import rc
from operator import itemgetter
import colorbrewer

# Colors from www.ColorBrewer.org by Cynthia A. Brewer, Geography,
#    Pennsylvania State University.
# Copyright (c) 2002 Cynthia Brewer, Mark Harrower, and The Pennsylvania State
# University.

__author__ = "Justine Debelius"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Justine Debelius", "Daniel McDonald"]
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
                    dictionary of meta data containing headers and observations
    """

    lines = [l.strip().split('\t') for l in mapping_data]
    header = lines[0]
    D2 = {}
    for l in lines[1:]:
        inner = {k: v for k, v in zip(header, l)}
        sample_id = inner['#SampleID']
        D2[sample_id] = inner

    return D2


def load_category_files(category_files):
    """Loads the category tables as biom files

    INPUTS:
        category_files -- a dictionary that associates the mapping category
                    (key) with the file path to the otu_table summarizing that

    OUTPUTS:
        category_tables -- a dictionary that associates the mapping category
                    with the summarized otu table for the category.
    """

    category_tables = {}
    watch_count = 0
    watch_list = []

    for (category, category_file) in category_files.iteritems():
        if isfile(category_file):
            with biom_open(category_file, 'U') as fp:
                cat_table = parse_biom_table(fp)
            category_tables[category] = cat_table
        else:
            watch_list.append('The summarized OTU table file cannot be found '
                              'for %s. \n%s is not in the file path.'
                              % (category, category_file))
            watch_count = watch_count + 1

    if watch_count > 0:
        print 'The following category files could not be found: \n%s' \
            % '\n'.join(watch_list)
    if watch_count == len(category_files):
        raise ValueError('No files could be found for any of the supplied '
                         'categories. \n%s' % '\n'.join(watch_list))

    return category_tables


def parse_category_files(raw_tables, common_groups, level=2,
                         metadata='taxonomy'):
    """ Collapses categeory tables using the most common OUTPUTS

    INPUTS:
    category_tables -- a dictionary keying the category name in the mapping
                file to the biom otu table of the collapsed data.

    common_group -- the reference groups in the metadata category which
                    should be used to summarize the data.

    level -- the level at which data should be summarized

    metadata -- the metadata category which should be used to summarize the
                data


    OUTPUTS:
        category_data -- a dictionary that associates the mapping category
                    with the data summarized using common_categories."""

    num_g = len(common_groups)
    category_data = {}
    for (cat, cat_table) in raw_tables.items():
        [ids, data, cats] = \
            summarize_common_categories(biom_table=cat_table,
                                        level=level,
                                        common_categories=common_groups[:num_g],
                                        metadata_category=metadata)
        category_data.update({cat: {'Groups': ids,
                                    'Summary': data}})
    return category_data


def identify_most_common_categories(biom_table, level, limit_mode='COMPOSITE',
                                    metadata_category='taxonomy', limit=1.0):
    """Identifies the most common taxa in a population using variable limits

    This method uses a composite score to account for both the average
    frequency and the percentage of the population in which the sample is
    present. This hopes to be more robust than either the average or fraction
    present alone.

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
                    composite value, and metadata category should be printed
                    and retained.

    OUTPUTS
        common_categories -- a list of the common categories in the (row
                    header)

        scores -- a sorted list of all the taxa and their scores"""

    # Sets up positions
    COMPOSITE_CONSTANT = 10000

    LIMIT_MODES = {'COMPOSITE': [3, True, 1],
                   'AVERAGE': [1, True, 3],
                   'COUNTS': [2, True, 3],
                   'NONE': [0, False, 3]}
    score_position, score_reverse, second_score = \
        LIMIT_MODES.get(limit_mode, [None, None, None])

    if score_position is None:
        raise ValueError("limit_mode is not a supported option. \n"
                         "Options for limit_mode are 'COMPOSITE', 'AVERAGE', "
                         "'COUNTS', or 'NONE'.")

    # Gets the Sample IDs
    sample_ids = list(biom_table.ids())
    num_samples = len(sample_ids)

    # Prealocates output objects
    scoring_all = []
    common_categories = []

    # Normalizes the data by the number of observations so relative frequencies
    # are used.
    biom_table_norm = biom_table.norm(inplace=False) #biom_table.normObservationBySample()

    # Collapses the OTUs into category summaries using the correct levels
    bin_fun = lambda y, x: x[metadata_category][:level]
    for (bin, table) in biom_table_norm.partition(bin_fun, axis='observation'): #binObservationsByMetadata(bin_fun):
        # Pulls out the sample data for the group
        group_value = array(table.sum('sample'))
        group_binary = group_value > 0

        # Calculates presence scores
        average_freq = round(mean(group_value), 4)
        fraction_pres = round(sum(group_binary)/(num_samples), 4)
        composite = round(average_freq*fraction_pres*COMPOSITE_CONSTANT, 2)
        score_row = [bin, average_freq, fraction_pres, composite]

        # Adds the scores to the watch matrix
        if fraction_pres > 0:
            scoring_all.append(score_row)

    # Sorts based on scoring method
    scores_all = sorted(sorted(scoring_all,
                               key=itemgetter(second_score),
                               reverse=True),
                        key=itemgetter(score_position),
                        reverse=score_reverse)

    # Identifies rows which meet the scoring criteria
    scores = []
    for score in scores_all:
        scores.append(score)
        if score[score_position] > limit:
            common_categories.append(score[0])

    # Raises an error if necessary
    if len(common_categories) == 0:
        raise ValueError('Limit too high! No common categories.')

    # Returns the values
    return common_categories, scores


def summarize_common_categories(biom_table, level, common_categories,
                                metadata_category='taxonomy'):
    """Determines the frequency of common categories present in a biom table

    INPUTS:
        biom_table -- a sparse biom table to be evaluated

        level -- an integer corresponding to the taxonomic level (or other meta
                    data category level) at which data should be summarized.
                    The argument will only take a single integer.

        common_categories -- a list of values which define the common
                    categories.

        metadata_category -- a description of the metadata category over which
                    the data will be summarized.

    OUTPUTS:
        sample_ids -- a list of the sample ids found in the OTU table.

        cat_summary -- a numpy array describing the frequency of the common
                    categories with an additional frequencies collapsed into
                    the "other" category.

        common_cats -- a summary of the common categories with an "other"
                    category appended."""

    # Checks that input biom table can be processed
    all_cats = biom_table.metadata(axis='observation')
    cats_zip = zip(*all_cats)

    # Gets the set of all categories
    group_all = ()
    for cat in cats_zip:
        group_all = group_all+(cat)

    categories = set(group_all)

    if not metadata_category in categories:
        raise ValueError('The biom table cannot be summarized; supplied '
                         'category does not exist.')

    # Identifies the common categories and removes extraneous characters
    num_cats = len(common_categories)

    for idx, cat in enumerate(common_categories):
        temp_cat = []
        for i in cat:
            temp_cat.append(i.strip())
        common_categories[idx] = tuple(temp_cat)

    sample_ids = biom_table.ids()
    num_samples = len(sample_ids)

    # Sets up the "other category name"
    summary_name = all_cats[0][metadata_category]
    other_name = [summary_name[0]]
    if len(summary_name) > 2:
        for cat_des in summary_name[1:(level)]:
            other_name.append('%s__%s' % (cat_des.split('__')[0], 'Other'))
    other_name = [tuple(other_name)]

    # Normalizes the biom table
    biom_norm = biom_table.norm(inplace=False) #normObservationBySample()

    # Prealocates numpy objects (because that makes life fun!). tax_other is
    # set up as a row array because this is ultimately summed
    cat_summary = zeros([num_cats, num_samples])

    # Collapses the OTU table using the category at the correct level
    bin_fun = lambda y, x: x[metadata_category][:level]
    for (bin_, table) in biom_norm.partition(bin_fun, axis='observation'):#.binObservationsByMetadata(bin_fun):
        new_bin = []
        for item in bin_:
            new_bin.append(item.strip())
        new_bin = tuple(new_bin)

        if new_bin in common_categories:
            current = cat_summary[common_categories.index(new_bin)]
            cat_summary[common_categories.index(new_bin)] \
                = table.sum('sample') + current

    cat_summary = vstack((cat_summary, 1 - sum(cat_summary)))
    common_cats = common_categories

    common_cats.extend(other_name)

    return sample_ids, cat_summary, common_cats


def translate_colors(num_colors, map_name='Spectral'):
    """Gets a colorbrewer colormap and sets it up for plotting in matplotlib

    OUTPUTS:
        colormap -- a numpy array with the colorbrewer map formatted for use in
                    matplotlib.
    """
    try:
        raw_map = getattr(colorbrewer, map_name)
    except:
        raise ValueError('%s is not a valid colorbrewer map name. '
                         '\nSee http://colorbrewer2.org for valid map names.'
                         % map_name)

    if num_colors not in raw_map:
        raise ValueError('Too many colors. \n%i wanted %i possible. \n'
                         'Pick fewer colors.'
                         % (num_colors, max(raw_map.keys())))

    map_ar = array(raw_map[num_colors])
    # Corrects for colorbrewer's 0 to 255 scaling and matplotlib's use of 0 to
    # 1 color sclaing.
    colormap = map_ar.astype(float)/255

    return colormap


def calculate_dimensions_rectangle(axis_width=4, axis_height=4, border=0.1,
                                   title=0.25, legend=1, xlab=0, ylab=0,
                                   unit='in'):
    """Determines the appriate axis and figure dimensions for square axis.

    INPUTS:
        axis_size -- a number specifying the side length the axis. DEFAULT: 4

        border -- the width of the border around the figure

        title -- the height to add to the top of the figure for a title. This
                    is separate from the border, which is added by default.
                    DEFAULT: 1

        legend -- the width of the legend to the be added to the figure.
                    DEFAULT: 2

        xlab -- the height to be added for labels along the x axis. DEFAULT: 0

        ylab -- the width to be added for labels along the y axis. DEFAULT: 0

        unit -- a string ('inches' or 'cm'), specifying the unit to be used
                    in image generation. DEFAULT: in

    OUTPUTS:
        axis_dimensions -- a Bbox class describing the axis position in the
            figure

        figure_dimensions -- a 2 element tuple giving the width and height of
            the figure in inches
    """

    # Specifies a value for converting between units
    if unit == 'cm':
        conversion = 1/2.54
    elif unit == 'in':
        conversion = 1.0
    else:
        raise ValueError('unit must be "in" or "cm".')

    # Determines the figure dimensions
    fig_width = axis_width+(border*2)+legend+ylab
    fig_height = axis_height+(border*2)+title+xlab

    figure_dimensions = (fig_width*conversion, fig_height*conversion)

    # Determines the axis bounds
    axis_left = (border+ylab)/fig_width
    axis_right = (border+axis_width+ylab)/fig_width
    axis_bottom = (border+xlab)/fig_height
    axis_top = (border+axis_height+xlab)/fig_height

    axis_dimensions = array([[axis_left,  axis_bottom],
                             [axis_right, axis_top   ]])

    return axis_dimensions, figure_dimensions


def calculate_dimensions_bar(num_bars, bar_width=0.5, axis_height=3,
                             border=0.1, title=1, legend=2, xlab=0, ylab=0,
                             unit='in'):
    """Determines the axis and figure dimensions for a bar chart.

    INPUTS:
        num_bars -- the number of bars in the bar chart being created.

        bar_width -- the width of the plotted bars in units DEFAULT: 0.5.

        axis_heigth -- the height of the axes DEFAULT: 3

        border -- the size of white space to be added around the axis.
                    DEFAULT: 0.1

        title -- the height to add to the top of the figure for a title. This
                    is separate from the border, which is added by default.
                    DEFAULT: 1

        legend -- the width of the legend to the be added to the figure.
                    DEFAULT: 2

        xlab -- the height to be added for labels along the x axis. DEFAULT: 0

        ylab -- the width to be added for labels along the y axis. DEFAULT: 0

        unit -- a string ('in' or 'cm'), specifying the unit to be used
                    in image generation. DEFAULT: in

    OUTPUTS:
        axis_dimensions -- a Bbox class describing the axis position in the
            figure

        figure_dimensions -- a 2 element tuple giving the width and height of
            the figure"""

    # Preforms some sanity checks.
    if num_bars < 1:
        raise ValueError('There must be at least one group to plot.\nnum_bars'
                         ' must be a whole number greater than 1.')
    elif not round(num_bars) == num_bars:
        raise ValueError('There cannot be partial categories. \nnum_bars must'
                         ' be a whole number greater than 1.')

    # Specifies a value for converting between units
    if unit == 'cm':
        conversion = 1/2.54
    elif unit == 'in':
        conversion = 1.0
    else:
        raise ValueError('unit must be "in" or "cm".')

    # Determines the figure width
    axis_width = (num_bars*bar_width)
    figure_width = (border*2+ylab+axis_width+legend)*conversion
    figure_height = (border*2+xlab+axis_height+title)*conversion

    figure_dimensions = (figure_width, figure_height)

    axis_left = (ylab+border)/figure_width*conversion
    axis_right = (ylab+border+axis_width)/figure_width*conversion
    axis_bottom = (xlab+border)/figure_height*conversion
    axis_top = (xlab+border+axis_height)/figure_height*conversion

    axis_dimensions = array([[axis_left,  axis_bottom],
                             [axis_right, axis_top]])

    return axis_dimensions, figure_dimensions


def render_single_pie(data_vec, group_names, axis_dims, fig_dims,
                      file_out='piechart', filetype='PDF', colors=None,
                      show_edge=True, axis_on=False, plot_ccw=False,
                      start_angle=90, x_lims=[-1.1, 1.1], y_lims=[-1.1, 1.1],
                      legend=True, legend_offset=None, legend_font=None,
                      legend_frame=False, title=None, title_font=None,
                      labels=None, label_distance=1.1, label_font=None,
                      use_latex=False, rc_fam='sans-serif',
                      rc_font=['Helvetica']):
    """Creates a pie chart summarizing the category data

    INPUTS:
        data_vec -- a vector which sums to 1 describing the fraction of the
                    chart represented by each group in group_names.

        group_names -- a list of the groups in cat_vec. (i.e. Firmictues,
                    Bacteriodetes, Proteobacteria for a category of Taxonomy).

        file_out -- a string giving the file path where the pie plot should be
                    saved.

        filetype -- a string describing the file format to save the output
                    file. Possible values include 'PNG', 'PDF', 'EPS', and
                    'SVG'.
                    DEFAULT: PDF

        colormap -- an n x 3 or n x 4 numpy array giving the desired colormap.
                    Default is to fill all wedges in white.
                    DEFAULT: array([[1,1,1]])

        show_edge -- a binary value dictating whether or the edge should be
                    outlined in black.
                    DEFAULT: True

        axis_on -- a binary value indicating whether or not the axes should be
                    displayed.
                    DEFAULT: False

        axis_frame -- a binary value indicating whether the frame should be
                    displayed on the axes
                    DEFAULT: False

        plot_ccw -- a binary value indicating whether whether the data should
                    be plotted clockwise (False) or counter-clockwise (True).
                    DEFAULT: False

        start_angle -- the angle from the x-axis (horizontal = 0) in degrees
                    at which to start plotting.
                    DEFAULT: 90

        x_lims -- the limits on the x-axis in the form of (min, max)
                    DEFAULT: [-1.1, 1.1]

        y_lims -- the limit sont he y axis in the form of (min, max)
                    DEFAULT: [-1.1, 1.1]

        legend -- a binary value indicating whether or not a legend should be
                    shown. This must be accounted for in the figure dimensions.
                    Default is to show the legend.

        legend_offset -- a two-element list giving the offset for the axes.
                    If this has a value of None, the legend will remain in the
                    original position.
                    DEFAULT: None

        legend_font -- a FontProperties object (dictionary) describing the
                    properties of the font. If None is supplied, the default is
                    to use the 15 pt normal sans-serif.
                    DEFAULT: None

        legend_frame -- a binary value indicating whether or not a the legend
                    should display a frame.

        title -- a string giving a title to append to the figure. This must be
                    accounted for in the figure dimensions.
                    DEFAULT: None

        title_font -- a FontProperties object (dictionary) describing the
                    properties of the font. If None is supplied, the default is
                    to use the 36 pt normal sans-serif.
                    DEFAULT: None

        labels -- a list of labels to be added to the figure. A value of None
                    will turn of labels.
                    DEFAULT: None

        label_distance -- the distance from the origin at which the labels
                    should be displayed.
                    DEFAULT: 1.1

        label_font -- a FontProperties object (dictionary) describing the font.
                    If None is supplied, the default is to use size 15 normal
                    sans-serif.
                    DEFAULT: None

        use_latex -- a binary value indicating if matplotlib's ability to
                    render using LaTeX should be considered or not. If this is
                    the case, the supplied rc_fam and rc_font will be used.
                    DEFAULT: False

        rc_family -- the font family which should be used for LaTeX rendering.
                    Options are 'sans-serif', 'serif', 'cursive'.
                    DEFAULT: 'sans-serif'

        rc_font -- a list of the font(s) which should be used with latex
                    rendering.
                    DEFAULT: ['Helvetica']

    OUTPUTS:
        The rendered figure is saved in the at the file_out location.
    """

    # Sets up the colormap
    num_wedges = len(data_vec)
    num_colors = len(colors[:, 0])

    if num_colors is None:
        colormap = ones((num_wedges, 3))
    elif not isinstance(colors, ndarray):
        raise TypeError('The colormap must be a numpy array.')
    elif num_colors == 1:
        colormap = colors*ones((num_wedges, 1))
    elif num_colors >= num_wedges:
        colormap = colors
    else:
        raise ValueError('The color map cannot be determined. \nColors must '
                         'be a a list of n x 3 lists where n is the number of '
                         'patches being supplied or a single color to be used'
                         ' for all patches.')

    # Sets up the font properties for each of the label objects
    if label_font is None:
        label_font = FontProperties()
        label_font.set_size(20)
        label_font.set_family('sans-serif')

    if legend_font is None:
        legend_font = FontProperties()
        legend_font.set_size(15)
        legend_font.set_family('sans-serif')

    if title_font is None:
        title_font = FontProperties()
        title_font.set_size(30)
        title_font.set_family('sans-serif')

    # Sets up LateX rendering if desired
    if use_latex:
        rc('text', usetex=True)
        rc('font', **{'family': rc_fam, rc_fam: rc_font})

    # Creates the figure
    fig = plt.gcf()
    fig.set_size_inches(fig_dims)
    ax1 = plt.axes(Bbox(axis_dims))
    ax1.set_position(Bbox(axis_dims))
    if axis_on:
        ax1.set_axis_on()

    # Plots the data clockwise
    [pie_patches, pie_text] = ax1.pie(x=data_vec,
                                      labels=labels,
                                      labeldistance=label_distance,
                                      shadow=False,
                                      startangle=start_angle)

    # Colors the data so its pretty!
    for idx, patch in enumerate(pie_patches):
        # Sets the face color
        patch.set_facecolor(colormap[idx, :])
        if not show_edge:
            patch.set_edgecolor(colormap[idx, :])

    # Sets the label properties
    if not labels is None:
        for lab in pie_text:
            lab.set_set_font_properties(label_font)

    # Sets the axis and figure dimensions
    plt.draw()

    # Reverses the axis dimensions for a counter-clockwise plot
    if not plot_ccw:
        plt.axis([x_lims[1], x_lims[0], y_lims[0], y_lims[1]])

    plt.draw()

    # Adds the legend if necessary
    if legend:
        leg = plt.legend(pie_patches, group_names,
                         loc='center right',
                         prop=legend_font,
                         frameon=legend_frame)
        plt.draw()

        if not legend_offset is None and len(legend_offset) == 2:
            leg.set_bbox_to_anchor((legend_offset[0], legend_offset[1]))
        elif not legend_offset is None:
            leg.set_bbox_to_anchor(tuple(legend_offset))

    plt.draw()

    # Adds the title if desired
    if isinstance(title, str):
        plt.title(title, prop=title_font)

    plt.draw()

    # Saves the output figure
    plt.savefig(file_out, format=filetype)
    plt.clf()


def render_barchart(data_table, group_names, sample_names, axis_dims,
    fig_dims, file_out='barchart', filetype='PDF', colors=None,
    show_edge=True, legend=True, title=None, match_legend=True,
    frame=True, bar_width=0.8, x_axis=True, x_label=None,
    x_min=-0.5, x_tick_interval=1.0, y_axis=True,
    y_lims=[0, 1], y_tick_interval=0.2, y_tick_labels=None,
    y_label=None, legend_frame=False,
    legend_offset=None, font_angle=45, font_alignment='right',
    tick_font=None, label_font=None, legend_font=None,
    title_font=None, use_latex=False, rc_fam='sans-serif',
    rc_font=['Helvetica']):
    """Creates a stacked bar chart using the data in the category table.

    A single value bar chart can be created using a vector for data_table
    instead of a table.

    INPUTS:
        data_table -- a numpy array of the category information to be plotted
                    where the rows are the groups in the category and the
                    columns are the samples.
        group_names -- a list of the groups, corresponding to rows in the
                    data table

        sample_names -- a list of the sample names, corresponding to the
                    columns in the data table. The sample names must be string
                    objects.

        axis_dims -- a 2 x 2 numpy array giving the fraction of the figure
                    which should bound the axis. (row 1: [left, bottom], row 2:
                    [right, top] as a percentage of the figure space)

        fig_dims -- a 2 element tuple giving the width and height in inches of
                    the output figure.

        file_out -- a string giving the file path where the pie plot should be
                    saved.

        filetype -- a string describing the file format to save the output
                    file.  Possible values include 'PNG', 'PDF', 'EPS', and
                    'SVG'.
                    DEFAULT: PDF

        colors -- an n x 3 or n x 4 numpy array giving the desired color map.
                    Default is to color each wedge white. (array([[1, 1, 1]]))

        show_edge -- a binary value dictating whether or the edge should be
                    outlined in black.
                    DEFAULT: True

        legend -- a binary value indicating whether or not a legend should be
                    shown. This must be accounted for in the figure dimensions.
                    DEFAULT: True

        title -- a string giving a title to append to the figure. This must be
                    accounted for in the figure dimensions.
                    DEFAULT: None

        match_legend -- a binary value indicating whether the order of colors
                    in the plot should match the order of colors in the legend.
                    DEFAULT: True

        frame -- a binary value indicating whether or not the a frame should be
                    displayed around the axis
                    DEFAULT: True

        bar_width -- the fraction of the bar width to be occupied by the data.
                    DEFAULT: 0.8

        x_axis -- a binary value indicating whether or not the x-axis should be
                    labeled.
                    DEFAULT: True

        x_label -- a string describing the the data in the plot's x-axis. A
                    value of None leaves the axis label off.
                    DEFAULT: None

        x_min -- the minimum value for the x-axis. The maximum is determined by
                    the number of bars, where each bar is one unit away from
                    the next.
                    DEFAULT: 0.5

        x_tick_interval -- the spacing between the the plotted bars.
                    DEFAULT: 1.0

        y_axis -- a binary value indicating whether or not tick labels should
                    be shown on the y axis.
                    DEFAULT: True

        y_lims -- a 2 element list giving the minimum and maximum values for
                    the y axis.
                    DEFAULT: [0, 1]

        y_tick_interval -- the spacing between ticks on the y-axis.
                    DEFAULT: 0.2

        y_tick_labels -- a string with the labels for the y-ticks. If no value
                    is supplied, the labels are set up using the y-tick values.
                    DEFAULT: None

        y_label -- a string describing the data plotted on the y-axis. If None,
                    no string will be present.
                    DEFAULT: None

        legend_frame -- a binary value indicating whether a box will be
                    displayed around the legend.
                    DEFAULT: False

        legend_offset -- a two-element list giving the offset for the axes.
                    If this has a value of None, the legend will remain in the
                    original position.
                    DEFAULT: None

        font_angle -- the angle in degrees at which the x-axis text should be
                    displayed.
                    DEFAULT: 45

        font_alignment -- the horizontal alignment of the x axis labels. Values
                    may be 'left', 'right' or 'center'.
                    DEFAULT: 'right'

        tick_font -- a FontProperties object (dictionary) describing the
                    properties of the font. If None is supplied, the default is
                    to use the 15 pt normal sans-serif.
                    DEFAULT: None

        label_font -- a FontProperties object (dictionary) describing the
                    properties of the font. If None is supplied, the default is
                    to use the 20 pt italic sans-serif.
                    DEFAULT: None

        legend_font -- a FontProperties object (dictionary) describing the
                    properties of the font. If None is supplied, the default is
                    to use the 15 pt normal sans-serif.
                    DEFAULT: None

        title_font -- a FontProperties object (dictionary) describing the
                    properties of the font. If None is supplied, the default is
                    to use the 36 pt normal sans-serif.
                    DEFAULT: None

         use_latex -- a binary value indicating if matplotlib's ability to
                    render using LaTeX should be considered or not. If this is
                    the case, the supplied rc_fam and rc_font will be used.
                    DEFAULT: False

        rc_family -- the font family which should be used for LaTeX rendering.
                    Options are 'sans-serif', 'serif', 'cursive'.
                    DEFAULT: 'sans-serif'

        rc_font -- a list of the font(s) which should be used with latex
                    rendering.
                    DEFAULT: ['Helvetica']

    OUTPUT:
        The rendered figure is saved in the file_out location.
    """

    # Preforms a sanity checks that the provided data is good
    (table_height, table_width) = data_table.shape
    num_cats = len(group_names)
    num_samples = len(sample_names)

    if not table_height == num_cats:
        raise ValueError('The number of provided categories differ.')
    elif not table_width == num_samples:
        raise ValueError('The number of samples differ.')

     # Sets up the colormap
    num_colors = len(colors[:, 0])
    if num_colors is None:
        colormap = ones((table_height, 3))
    elif not isinstance(colors, ndarray):
        raise TypeError('The colormap must be a numpy array.')
    elif num_colors == 1:
        colormap = colors*ones((table_height, 1))
    elif num_colors >= table_height:
        colormap = colors
    else:
        raise ValueError('The color map cannot be determined. \nColors must '
                         'be a a list of n x 3 lists where n is the number of '
                         'patches being supplied or a single color to be used'
                         ' for all patches.')

    # Sets up the edge colormap
    if show_edge:
        edgecolor = zeros((num_cats, 3))
    else:
        edgecolor = colormap

    # Sets up the font properties for each of the label objects
    if label_font is None:
        label_font = FontProperties()
        label_font.set_size(20)
        label_font.set_family('sans-serif')
        label_font.set_style('italic')

    if legend_font is None:
        legend_font = FontProperties()
        legend_font.set_size(15)
        legend_font.set_family('sans-serif')

    if tick_font is None:
        tick_font = FontProperties()
        tick_font.set_size(15)
        tick_font.set_family('sans-serif')

    if title_font is None:
        title_font = FontProperties()
        title_font.set_size(30)
        title_font.set_family('sans-serif')

     # Sets up LateX rendering if desired
    if use_latex:
        rc('text', usetex=True)
        rc('font', **{'family': rc_fam, rc_fam: rc_font})

    # Sets up the x ticks.
    # Bar width is divided by two because the tick is assumed to be at the
    # center of the bar.
    x_tick = arange(num_samples*x_tick_interval)
    x_max = x_min + num_samples*x_tick_interval
    bar_left = x_tick - bar_width/2

    # Creates the x tick labels.
    if x_axis:
        x_text_labels = map(str, sample_names)
    else:
        x_text_labels = ['']*num_samples

    # Creates the y tick labels
    if y_tick_labels is None:
        y_tick_labels = arange(y_lims[1] + y_tick_interval, y_lims[0],
                               -y_tick_interval)
        y_tick_labels = y_tick_labels - y_tick_interval
        y_tick_labels[-1] = y_lims[0]

    num_y_ticks = len(y_tick_labels)

    if y_axis:
        y_text_labels = map(str, y_tick_labels)
    else:
        y_text_labels = ['']*num_y_ticks

    # Plots the data
    fig = plt.figure()
    fig.set_size_inches(fig_dims)
    ax1 = plt.axes(Bbox(axis_dims))

    patches_watch = []
    for plot_count, category in enumerate(data_table):
        bottom_bar = sum(data_table[0:plot_count, :], 0)
        faces = ax1.bar(bar_left, category, bar_width, bottom_bar,
                        color=colormap[plot_count, :],
                        edgecolor=edgecolor[plot_count, :])
        patches_watch.append(faces[0])

    # The y-axis is reversed so the labels are in the same order as the
    # colors in the legend
    if match_legend:
        plt.axis([x_min, x_max, y_lims[1], y_lims[0]])

    # Sets up y labels if they are desired.

    y_tick_labels = ax1.set_yticklabels(y_text_labels,
                                        fontproperties=tick_font)
    if not y_label is None:
        ax1.set_ylabel(y_label, fontproperties=label_font)

    # Set the x-axis labels
    ax1.set_xticks(x_tick)
    ax1.set_xticklabels(x_text_labels,
                        rotation=font_angle,
                        horizontalalignment='right',
                        fontproperties=label_font)

    if x_label is not None:
            ax1.set_xlabel(x_label, fontproperties=label_font)

    if legend:
        leg = plt.legend(patches_watch, group_names, prop=legend_font)
        if not legend_offset is None:
            leg.set_bbox_to_anchor((legend_offset[0], legend_offset[1]))

    if isinstance(title, str):
        plt.title(title, fontproperties=title_font)

    plt.savefig(file_out, format=filetype)
