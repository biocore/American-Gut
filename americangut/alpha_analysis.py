from __future__ import division
from os import mkdir
from os.path import exists
from numpy import nan, arange
from scipy.stats import kruskal
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pandas import DataFrame, Series

# Sets up plotting parameters so that the default setting is use to Helvetica
# in plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Helvetica', 'Arial']
rcParams['text.usetex'] = True

__author__ = "Justine Debelius"
__copyright__ = "Copyright 2014, The American Gut Project"
__credits__ = ["Justine Debelius"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Justine Debelius"
__email__ = "Justine.Debelius@colorado.edu"


# Creates a quick inline function
def check_dir(dir_):
    """Creates the specified directory if it does not exist

    Parameters
    ----------
    dir : str
        the directory to be checked
    """

    if not exists(dir_):
        mkdir(dir_)


def pad_index(df, index_col='#SampleID', nzeros=9):
    """Adds zeros to the sample ID strings

    Parameters
    ----------
    df : dataframe

    index_col : {#SampleID, str}
        the name of the column containing the index data

    n_zeros : {9, int}
        the number of zeros to add before the string

    Returns
    -------
    df : dataframe
        the dataframe with an appropriate index column
    """

    # Gets the sample IDs
    samples = df[index_col].values
    new_samples = []

    # Pads the zeros on the id
    for samp in samples:
        if not isinstance(samp, str):
            samp = str(samp)
        splits = samp.split('.')
        first_clean = [splits[0].zfill(nzeros)]
        first_clean.extend(splits[1:])
        new_samples.append('.'.join(first_clean))

    # Sets the column as the index
    df.index = new_samples
    del df[index_col]

    # Returns the dataframe
    return df


def pretty_boxplot(grouped, order, cat, **kwargs):
    """
    Creates a more attractive boxplot than pandas

    Parameters
    ----------
    grouped : pandas grouped object
        the dataframe grouped by the category of interest

    order : {list, None}
        if category order matters, this argument sets the order. Otherwise,
        default sorting is used.

    cat : str
        the dataframe column header being plotted. Ideally, this should be
        continous data.

    Returns
    -------
    fig : matplotlib figure
        the figure object displaying the graph

    features : dictionary
        a dictionary keyed to the axis, boxplot features, kruskal-wallis test
        statitics, numerical test objects, and p-value text.

    Other Parameters
    ----------------
    axis_dims : {tuple}
        the portion of the figure covered by the axis. The four element tuple
        provides the positions of the (Left, Bottom, Width, Height) for the
        axis.

    fig_dims : {None, tuple}
        the size of the figure, in inches

    interval : {0.5, float}
        Spacing between each boxplot

    notch : {True, False}
        Have the boxplot display the 95% confidence interval notch

    show_n : {True, False}
        display group sizes on the figure

    n_xs : {None, list}
        the x-positions to display counts. None will place the counts at the
        x-ticks.

    n_y : {None, float}
        the y-position to display the counts. None will place the counts at the
        bottom of the axis.

    n_size : {12, int}
        the text size for the displayed counts

    show_p : {True, False}
        display the p value on the axis

    p_x : {None, float}
        the x-position for the p-value string. None will place it to the far
        right

    p_y : {None, float}
        the y-position for the p-value string. None will place the p text at
        the top of the figure.

    p-size : {12, int}
        the font size for the displayed p-value

    title : {'', str}
        title for the figure

    title_size : {15, int}
        the font size for the displayed title

    xtick_size : {12, int}
        the text size for xtick labels

    xfont_align : {'left', 'center', 'right'}
        the alignment for the xtick labels. Right alignment is recommend with
        rotated text

    xfont_angle : {0, int}
        the angle for xtick labels.

    xlabel : {'', str}
        x-axis label

    xlabel_size : {15, int}
        the text size for the xlabel

    show_xgrid : {False, True}
        display vertical grid lines on the axis

    ylims : {list, tuple, array}
        Upper and lower limits for the y-axis

    ytick_size : {12, int}
        the text size for ytick labels

    ylabel : {'', str}
        y-axis label

    ylabel_size : {15, int}
        the font size for the y axis label

    show_ygrid : {True, False}
        display horizontal grid lines on the axis

    Also See
    --------
    boxplot : maplotlib.pyplot.boxplot

    kruskal-wallis : scipy.stats.kruskal
    """

    # Handles keyword arguemnts
    keywords = {'fig_dims': None,
                'axis_dims': (0.1, 0.1, 0.8, 0.8),
                'notch': True,
                'interval': 0.5,
                'show_p': True,
                'show_n': True,
                'title': '',
                'tick_names': None,
                'xlabel': '',
                'ylabel': '',
                'ylims': None,
                'p_size': 12,
                'p_x': None,
                'p_y': None,
                'n_xs': None,
                'n_y': None,
                'n_size': 12,
                'title_size': 15,
                'xtick_size': 12,
                'xlabel_size': 15,
                'xfont_align': 'center',
                'xfont_angle': 0,
                'ytick_size': 12,
                'ylabel_size': 15,
                'show_xgrid': False,
                'show_ygrid': True}

    for key, val in kwargs.iteritems():
        if key in keywords:
            keywords[key] = val
        else:
            raise ValueError('%s is not an input for plot_with_dist.' % key)

    # Prealocates a dictionary of plotting features
    features = {}

    # Calculates plotting features
    if order is None:
        order = grouped.groups.keys()

    num_cats = len(order)
    xlim = [-keywords['interval']/2,
            keywords['interval']*(num_cats-1)+keywords['interval']/2]

    # Creates the figure
    fig = plt.figure()
    if keywords['fig_dims'] is not None:
        fig.set_size_inches(keywords['fig_dims'])
    ax = fig.add_axes(keywords['axis_dims'])

    # Sets up plotting constants
    ticks = arange(0, keywords['interval']*num_cats, keywords['interval'])
    names = []
    values = []
    counts = []

    # Appends the information to the features
    features['fig'] = fig
    features['axis'] = ax

    for idx, group in enumerate(order):
        tick = [ticks[idx]]
        value = grouped.get_group(group)[cat].values
        values.append(value)
        counts.append(len(value))
        bp = ax.boxplot(value,
                        positions=tick,
                        notch=keywords['notch'])

        # Handles group naming
        if keywords['tick_names'] is not None:
            names.append(keywords['tick_names'][group])
        else:
            names.append(group)

    # Adds the boxplot to the features dictionary
    features['boxplot'] = bp

    # Calculates the test statitics
    (features['h'], features['p']) = kruskal(*values)

    # Sets up the axis properties
    ax.set_xlim(xlim)
    if keywords['ylims'] is not None:
        ax.set_ylim(keywords['ylims'])
    ax.set_xticks(ticks)
    ax.set_xticklabels(names,
                       rotation=keywords['xfont_angle'],
                       horizontalalignment=keywords['xfont_align'],
                       size=keywords['xtick_size'])
    ax.set_yticklabels(ax.get_yticks(), size=keywords['ytick_size'])
    ax.set_xlabel(keywords['xlabel'], size=keywords['xlabel_size'])
    ax.set_ylabel(keywords['ylabel'], size=keywords['ylabel_size'])

    # Adds grid, if appropriate
    if keywords['show_xgrid']:
        ax.xaxis.grid()
    if keywords['show_ygrid']:
        ax.yaxis.grid()

    # Adds the title, if appropriate
    ax.set_title(keywords['title'], size=keywords['title_size'])

    text = []
    # Adds the group counts if desired
    if keywords['show_n']:
        ylim = ax.get_ylim()
        if keywords['n_xs'] is None:
            keywords['n_xs'] = ticks
        if keywords['n_y'] is None:
            keywords['n_y'] = ylim[0]+(ylim[1]-ylim[0])*0.025

        for idx, count in enumerate(counts):
            text.append(ax.text(x=keywords['n_xs'][idx],
                                y=keywords['n_y'],
                                s='(%i)' % count,
                                horizontalalignment='center',
                                size=keywords['n_size']))

    # Adds a p-value to the plot if desired
    if keywords['show_p']:
        # Gets the positions if necessary
        if keywords['p_x'] is None:
            xticks = ax.get_xticks()
            x = ax.get_xlim()[1] - (xticks[-1] - xticks[-2])/20
        else:
            x = keywords['p_x']

        if keywords['p_y'] is None:
            yticks = ax.get_yticks()
            y = yticks[-2]+(yticks[-1]-yticks[-2])/2
        else:
            y = keywords['p_y']

        # Adds the text
        if features['p'] >= 0.005:
            p_text = ax.text(s='p = %1.2f' % features['p'],
                             x=x, y=y,
                             size=keywords['p_size'],
                             horizontalalignment='right')
        else:
            p_text = ax.text(s='p = %1.1e' % features['p'],
                             x=x, y=y,
                             size=keywords['p_size'],
                             horizontalalignment='right')
        features['p_text'] = p_text

    fig = plt.gcf()

    return fig, features
