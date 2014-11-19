from __future__ import division

from os import mkdir
from os.path import exists

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import skbio

from matplotlib import rcParams
from scipy.stats import kruskal
from skbio.stats.power import _check_strs
from statsmodels.sandbox.stats.multicomp import multipletests

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
        splits = samp.split('.', 1)
        first_clean = [splits[0].zfill(nzeros)]
        first_clean.extend(splits[1:])
        new_samples.append('.'.join(first_clean))

    # Sets the column as the index
    df.index = new_samples
    del df[index_col]

    # Returns the dataframe
    return df


def boxplot(vecs, ax=None, **kwargs):
    """Makes a more attractive boxplot

    Parameters
    ----------
    vecs : list
        The list of arrays to plot as boxplots. The list format allows the
        arrays to be of uneven length.
    ax : matplotlib axis, optional
        The axis where data should be plotted. If none, a new axis instance
        will be created.
    interval : float, optional
        The spacing between the boxplot instances on the axes
    notch : bool, optional
        Displays the parametric 95% confidence interval around the mean.
    show_n : bool, optional
        Shows the size of the groups below each plot on the x-axis
    p_value : float, optional
        Default is None. When supplied, the signfigance value will be displayed
        on the plot in the upper right hand corner by default.
    show_xgrid: bool, optional
        Default is False. Adds vertical lines at each major x-tick.
    show_ygrid: bool, optional
        Default is True. Adds horizonal lines at each major y-tick.
    title: str, optional
        The title to be placed on the graph.
    ylims : list
        The limits for the y-axis.
    ylabel : str
        The label text for the y-axis. Every time you leave off appropriate
        labels and units, a science grad student grading lab reports cries
        another bitter tear into their bottle of craft beer.

    Returns
    -------
    ax : axes
        A matplotlib axes containing the plotted data
    feats : dict
        A dictionary with features of the plot, including the boxplot handles
        and text objects.

    Other Parameters
    ----------------
    hide_x_ticks : bool, optional
        Display x-tick symbols on the plot
    hide_y_ticks : bool, optional
        Display y-tick symbols on the plot
    p_x : float
        The x position of the critical value text
    p_y : float
        The y position of the critical value text
    p_size : int
        The font size for hte critical value text
    title_size: int
        The font size for the title
    xticklabels : list
        The strings to label each point on the x-axis.
    xfont_angle : float
        The angle in degrees for the x tick label text.
    xfont_align : {'left', 'right', 'center'}
        The horizonal alignment of the x tick label text. For rotated text,
        an alignment or 'right' is recommended.
    xlabel_size : int
        The font size of the x-axis label.
    xtick_size : int
        The font size for the xtick labels
    yticks : array_like
        The positions where ticks should appear on the y-axis.
    yticklabels : list
        The text to be displayed at each y tick.
    ylabel_size : int
        The font size of the y-axis label.
    ytick_size : int
        The font size for the ytick labels
    """

    # Sets up an axes instance if necessary
    if ax is None:
        ax = plt.axes()

    # Handles keyword arguments
    kwds = {'boxplot_props': {},
            'hide_xticks': False,
            'hide_yticks': False,
            'interval': 0.5,
            'notch': True,
            'n_xs': None,
            'n_y': None,
            'n_size': 11,
            'p_value': None,
            'p_x': None,
            'p_y': None,
            'p_size': 12,
            'show_n': True,
            'show_xgrid': False,
            'show_ygrid': True,
            'title': '',
            'title_size': 18,
            'xlabel': '',
            'xticklabels': None,
            'xfont_align': 'center',
            'xfont_angle': 0,
            'xtick_size': 12,
            'xlabel_size': 15,
            'ylims': None,
            'yticks': None,
            'yticklabels': None,
            'ylabel': '',
            'ylabel_size': 15,
            'ytick_size': 12}

    for key, val in kwargs.iteritems():
        if key in kwds:
            kwds[key] = val
        else:
            raise ValueError('%s is not an input for boxplot.' % key)

    features = {}

    # Determines the plotting locations
    num_cats = len(vecs)
    xlim = [-kwds['interval']/2,
            kwds['interval']*(num_cats-1)+kwds['interval']/2]

    # Sets up the plotting constants
    ticks = np.arange(0, kwds['interval']*num_cats, kwds['interval'])
    counts = []

    boxes = []
    # Loops through the data
    for idx, vec in enumerate(vecs):
        # Gets vector characteristics
        tick = [ticks[idx]]
        counts.append(len(vec))
        # Plots the data
        bp = ax.boxplot(vec,
                        positions=tick,
                        notch=kwds['notch'],
                        **kwds['boxplot_props'])
        boxes.append(bp)

    features['boxplot'] = boxes

    # Sets axis limits
    ax.set_xlim(xlim)
    if kwds['ylims'] is not None:
        ax.set_ylim(kwds['ylims'])
    else:
        ylim = ax.get_ylim()
        ax.set_ylim([ylim[0], ylim[1]*1.2])

    # Adds axis ticks
    ax.set_xticks(ticks)
    if kwds['yticks'] is not None:
        ax.set_yticks(kwds['yticks'])

    # Handles x-axis tick names
    if kwds['xticklabels'] is None:
        xtl = ''
    else:
        xtl = kwds['xticklabels']
    ax.set_xticklabels(xtl,
                       ha=kwds['xfont_align'],
                       rotation=kwds['xfont_angle'],
                       size=kwds['xtick_size'])
    if kwds['yticklabels'] is None:
        ax.set_yticklabels(ax.get_yticks(), size=kwds['ytick_size'])
    else:
        ax.set_yticklabels(kwds['yticklabels'], size=kwds['ytick_size'])

    # Sets the axis labels
    ax.set_xlabel(kwds['xlabel'], size=kwds['xlabel_size'])
    ax.set_ylabel(kwds['ylabel'], size=kwds['ylabel_size'])

    # Adds a grid, if appropriate
    if kwds['show_xgrid']:
        ax.xaxis.grid()
    if kwds['show_ygrid']:
        ax.yaxis.grid()

    # Hides axis ticks, if desired
    if kwds['hide_xticks']:
        for tic in ax.xaxis.get_major_ticks():
            tic.tick1On = False
            tic.tick2On = False

    if kwds['hide_yticks']:
        for tic in ax.yaxis.get_major_ticks():
            tic.tick1On = False
            tic.tick2On = False

    # Adds a title, if appropriate
    ax.set_title(kwds['title'], size=kwds['title_size'])

    # Adds group counts, if desired
    txt = []
    if kwds['show_n']:
        # Gets the positions for the count labels
        ylim = ax.get_ylim()
        if kwds['n_xs'] is None:
            kwds['n_xs'] = ticks
        if kwds['n_y'] is None:
            kwds['n_y'] = ylim[0]+(ylim[1]-ylim[0])*0.025

        # Adds the count labels
        for idx, count in enumerate(counts):
            txt.append(ax.text(kwds['n_xs'][idx],
                               kwds['n_y'],
                               '(%i)' % count,
                               horizontalalignment='center',
                               size=kwds['n_size']))
        features['n_text'] = txt
    else:
        features['n_text'] = None

    # Adds the p-value, if desired
    if kwds['p_value'] is not None:
        p = kwds['p_value']
        # Auto-calculates the position, if necessary
        if kwds['p_x'] is None:
            x = ax.get_xlim()[1] - (ax.get_xlim()[1] - ax.get_xlim()[0])/100
        else:
            x = kwds['p_x']
        if kwds['p_y'] is None:
            yticks = ax.get_yticks()
            y = ax.get_ylim()[1] - (yticks[-1] - yticks[-2])/2.5
        else:
            y = kwds['p_y']

        # Adds the text
        if p >= 0.005:
            p_str = 'p = %1.2f' % p
        else:
            p_str = 'p = %1.1e' % p
        p_text = ax.text(x, y, p_str,
                         size=kwds['p_size'],
                         horizontalalignment='right')
        features['p_text'] = p_text
    else:
        features['p_text'] = None

    return ax, features


def pretty_pandas_boxplot(meta, group, cat, order=None, ax=None,
    **boxplot_props):
    """Creates a more attractive poxplot than pandas

    Parameters
    ----------
    meta : pandas dataframe
        The metadata for the variable containing a column with a continous
        varaible, designated in `cat`, and a categorical variable, `group`
        with categories given by `order`.
    group : str
        The name of a column in meta which is a categorical predictor variable
    cat : str
        A column in meta which contains a continous response variable
    order : list, optional
        The order of categories in `group`. This can be used to limit
        the categories plotted. For instance, if there are three categories in
        `group`: A, B and C, and you only wish to compare A and C, you can
        list order as ['A', 'C'] to limit the categories.
    interval : float, optional
        The spacing between the boxplot instances on the axes
    notch : bool, optional
        Displays the parametric 95% confidence interval around the mean.
    show_n : bool, optional
        Shows the size of the groups below each plot on the x-axis
    p_value : float, optional
        Default is None. When supplied, the signfigance value will be displayed
        on the plot in the upper right hand corner by default.
    show_xgrid: bool, optional
        Default is False. Adds vertical lines at each major x-tick.
    show_ygrid: bool, optional
        Default is True. Adds horizonal lines at each major y-tick.
    title: str, optional
        The title to be placed on the graph.
    ylims : list
        The limits for the y-axis.
    ylabel : str
        The label text for the y-axis. Every time you leave off appropriate
        labels and units, a science grad student grading lab reports cries
        another bitter tear into their bottle of craft beer.

    Returns
    -------
    ax : axes
        A matplotlib axes containing the plotted data
    feats : dict
        A dictionary with features of the plot, including the boxplot handles
        and text objects.
    """

    grouped = meta.groupby(group)

    # Sets up the plotting order
    if order is None:
        order = grouped.groups.keys()

    # Gets the data vectors
    vecs = [grouped.get_group(g)[cat].values for g in order]

    # Checks the group names
    if 'xticklabels' not in boxplot_props:
        boxplot_props['xticklabels'] = order

    # Calculates the p value
    h, p = kruskal(*vecs)

    # Sets the boxplot properties
    ax, feats = boxplot(vecs=vecs, ax=ax, p_value=p, **boxplot_props)

    return ax, feats


def post_hoc_pandas(meta, group, cat, order=None, correct=None):
    """Preforms an post-hoc comparison between two groups

    Parameters
    ----------
    meta : pandas DataFrame
        the metadata object for the data
    group : str
        the metadata category being interograted
    cat : str
        the name of the column with the result
    order : None, list
        Default is None. The order of groups in the category.
    correct : None, str
        Method for multiple hypothesis correction using
        `statsmodels.sandbox.stats.multicomp.multipletests`. Methods you're
        likely to use are `bonferroni` and `fdr_bh`.

    Returns
    -------
    post_hoc : dataframe
    """

    # Groups the data
    grouped = meta.groupby(group)

    # Gets the order
    if order is None:
        order = grouped.groups.keys()

    # Prealocates an output frame
    stats = pd.DataFrame({'Counts': grouped[cat].count(),
                          'Mean':  grouped[cat].mean(),
                          'Stdv': grouped[cat].std(),
                          'Median': grouped[cat].median()})

    # Preforms ad-hoc comparisons
    comparison = {}

    for g1_name in order[:-1]:
        # Determines the position of the group
        pos = order.index(g1_name)
        g1_data = grouped.get_group(g1_name)[cat]
        compare = []
        index = []
        for id2, g2_name in enumerate(order):
            if id2 <= pos:
                index.append(g2_name)
                compare.append(np.nan)
            else:
                g2_data = grouped.get_group(g2_name)[cat]
                index.append(g2_name)
                compare.append(kruskal(g1_data, g2_data)[1])
        add_series = pd.Series(compare, index=index)
        comparison[g1_name] = add_series

    # Converts the data to a dataframe
    compare = pd.DataFrame(comparison)
    post_hoc = stats.join(compare[order[:-1]])
    post_hoc = post_hoc.reindex(order)

    # Performs the multiple hypothesis correction
    if correct is not None:
        post_hoc = multiple_correct_post_hoc(post_hoc, order, correct)

    return post_hoc


def multiple_correct_post_hoc(raw_ph, order, alphafwer=0.05,
    method='bonferroni'):
    """Performs multiple hypothesis correction on post hoc test matrices"""
    # Gets the positon matrix
    num_rows = len(order)
    num_cols = num_rows - 1
    pos = np.vstack([np.arange(0, num_cols) + i*num_cols for i in
                     range(num_rows)])

    # Draws the results and reshapes
    raw_ps = raw_ph.loc[order, order[:-1]].values
    ps_all = raw_ps.reshape(pos.max() + 1)
    pos_all = pos.reshape(pos.max() + 1)

    # Sorts the data
    ps_ord = np.argsort(ps_all)
    ps_sort = ps_all[ps_ord]
    pos_sort = pos_all[ps_ord]
    # Identifies the position of missing values in the sorted matrix
    ps_nan = np.isnan(ps_sort) == False

    # Corrects for multiple hypotheses
    reject, p_corr, asidak, abonf = multipletests(ps_sort[ps_nan],
                                                  alpha=0.05,
                                                  method=method)

    # Sorts the data back into its orginal order
    ps_sort[ps_nan] = p_corr
    pos_reord = np.argsort(pos_sort)
    ps_corr = ps_sort[pos_reord]

    raw_ph.loc[order, order[:-1]] = ps_corr.reshape(raw_ps.shape)

    return raw_ph


def barchart(height, interval=0.5, width=0.4, ax=None, errors=None, **kwargs):
    """Renders a barchart

    Parameters
    ----------
    height : array_like
        The height of each of the bars
    interval : float, optional
        The spacing between the bars
    width : float, optional
        The width of each bars. Should be less than or equal to the interal.
    ax : matplotlib axis, optional
        The axis where data should be plotted. If none, a new axis instance
        will be created.
    errors : array_like, optional
        The error bars assoicated with the groups being plotted.
    colormap : array-like, optional
        An n x 3 or n x 4 array of colors. If none is supplied, the bar
        facecolors will be white.
    match_colors: bool, optional
        If the error bars and plot edges should be the same color as the bar
        (True) or black (False).
    p_value : float, optional
        Default is None. When supplied, the signfigance value will be displayed
        on the plot in the upper right hand corner by default.
    show_xgrid: bool, optional
        Default is False. Adds vertical lines at each major x-tick.
    show_ygrid: bool, optional
        Default is True. Adds horizonal lines at each major y-tick.
    title: str, optional
        The title to be placed on the graph.
    ylims : list
        The limits for the y-axis.
    ylabel : str
        The label text for the y-axis. Every time you leave off appropriate
        labels and units, a science grad student grading lab reports cries
        another bitter tear into their bottle of craft beer.

    Returns
    -------
    ax : axes
        A matplotlib axes containing the plotted data
    feats : dict
        A dictionary with features of the plot, including the boxplot handles
        and text objects.
    xpos : array
        the center position for each bar

    Other Parameters
    ----------------
    elinewidth : int, optional
        The thickness of the errorbar in pixels
    e_capthickness : int, optional
        The tickness of the error bar cap in pixels
    hide_x_ticks : bool, optional
        Display x-tick symbols on the plot
    hide_y_ticks : bool, optional
        Display y-tick symbols on the plot
    p_x : float
        The x position of the critical value text
    p_y : float
        The y position of the critical value text
    p_size : int
        The font size for hte critical value text
    title_size: int
        The font size for the title
    xticklabels : list
        The strings to label each point on the x-axis.
    xfont_angle : float
        The angle in degrees for the x tick label text.
    xfont_align : {'left', 'right', 'center'}
        The horizonal alignment of the x tick label text. For rotated text,
        an alignment or 'right' is recommended.
    xlabel_size : int
        The font size of the x-axis label.
    xtick_size : int
        The font size for the xtick labels
    yticks : array_like
        The positions where ticks should appear on the y-axis.
    yticklabels : list
        The text to be displayed at each y tick.
    ylabel_size : int
        The font size of the y-axis label.
    ytick_size : int
        The font size for the ytick labels
    """
    # Sets up an axes instance if necessary
    if ax is None:
        ax = plt.axes()

    # Handles keyword arguments
    kwds = {'colormap': None,
            'edgecolors': None,
            'elinewidth': 2,
            'e_capthickness': 2,
            'hide_xticks': True,
            'hide_yticks': False,
            'ebar_ticks': None,
            'match_colors': True,
            'p_value': None,
            'p_x': None,
            'p_y': None,
            'p_size': 12,
            'show_xgrid': False,
            'show_ygrid': False,
            'title': '',
            'title_size': 18,
            'text_ticks': None,
            'ebar_ticks': None,
            'xlabel': '',
            'xticklabels': None,
            'xfont_align': 'center',
            'xfont_angle': 0,
            'xtick_size': 12,
            'xlabel_size': 15,
            'ylims': None,
            'yticks': None,
            'yticklabels': None,
            'ylabel': '',
            'ylabel_size': 15,
            'ytick_size': 12}

    for key, val in kwargs.iteritems():
        if key in kwds:
            kwds[key] = val
        else:
            raise ValueError('%s is not an input for boxplot.' % key)

    # Sets the colormap and egdecolor, if necessary
    if kwds['colormap'] is None:
        kwds['colormap'] = np.array([[0.5, 0.5, 0.5]]*len(height))

    if kwds['match_colors']:
        kwds['edgecolors'] = kwds['colormap']
    else:
        kwds['edgecolors'] = np.array([[0, 0, 0]]*len(height))

    # Gets an axis instance
    if ax is None:
        ax = plt.axes()

    feats = {}

    # Gets the xposition, for errorbars, if unspecified.
    if kwds['ebar_ticks'] is None:
        kwds['ebar_ticks'] = np.arange(0, len(height))*interval + interval/2
    if kwds['text_ticks'] is None:
        kwds['text_ticks'] = kwds['ebar_ticks']
    xpos = kwds['ebar_ticks']
    xleft = np.arange(0, len(height)) * interval + (interval - width)/2
    xlims = [0, len(height)*interval]

    # Plots the errorbars
    if errors is not None:
        e_bars = []
        for idx, err in enumerate(errors):
            eb = ax.errorbar(x=xpos[idx],
                             y=height[idx],
                             yerr=err,
                             fmt='none',
                             ecolor=kwds['edgecolors'][idx, :],
                             elinewidth=kwds['elinewidth'],
                             capthick=kwds['e_capthickness'])
        e_bars.append(eb)

    # Plots the bars
    bars = ax.bar(xleft, height, width=width)
    for idx, bar in enumerate(bars):
        bar.set_facecolor(kwds['colormap'][idx, :])
        bar.set_edgecolor(kwds['edgecolors'][idx, :])

    feats['errorbars'] = e_bars
    feats['bars'] = bars

    # Sets x-axis properties
    ax.set_xlim(xlims)
    ax.set_xticks(kwds['text_ticks'])
    if kwds['xticklabels'] is None:
        xtl = ''
    else:
        xtl = kwds['xticklabels']
    ax.set_xticklabels(xtl,
                       ha=kwds['xfont_align'],
                       rotation=kwds['xfont_angle'],
                       size=kwds['xtick_size'])
    ax.set_xlabel(kwds['xlabel'], size=kwds['xlabel_size'])
    if kwds['show_xgrid']:
        ax.xaxis.grid()
    if kwds['hide_xticks']:
        for tic in ax.xaxis.get_major_ticks():
            tic.tick1On = False
            tic.tick2On = False

    # Sets the y-axis parameters
    ax.set_ylim(kwds['ylims'])
    if kwds['yticks'] is not None:
        ax.set_yticks(kwds['yticks'])
    if kwds['yticklabels'] is None:
        ax.set_yticklabels(ax.get_yticks(), size=kwds['ytick_size'])
    else:
        ax.set_yticklabels(kwds['yticklabels'], size=kwds['ytick_size'])
    ax.set_ylabel(kwds['ylabel'], size=kwds['ylabel_size'])
    if kwds['show_ygrid']:
        ax.yaxis.grid()
    if kwds['hide_yticks']:
        for tic in ax.yaxis.get_major_ticks():
            tic.tick1On = False
            tic.tick2On = False

    # Adds a title, if appropriate
    ax.set_title(kwds['title'], size=kwds['title_size'])

    # Adds p-value text
    if kwds['p_value'] is not None:
        p = kwds['p_value']
        # Auto-calculates the position, if necessary
        if kwds['p_x'] is None:
            x = ax.get_xlim()[1] - (ax.get_xlim()[1] - ax.get_xlim()[0])/100
        else:
            x = kwds['p_x']
        if kwds['p_y'] is None:
            yticks = ax.get_yticks()
            y = ax.get_ylim()[1] - (yticks[-1] - yticks[-2])/2.5
        else:
            y = kwds['p_y']

        # Adds the text
        if p >= 0.005:
            p_str = 'p = %1.2f' % p
        else:
            p_str = 'p = %1.1e' % p
        p_text = ax.text(x, y, p_str,
                         size=kwds['p_size'],
                         horizontalalignment='right')
        feats['p_text'] = p_text
    else:
        feats['p_text'] = None

    return ax, feats, xpos


def add_comparison_bars(centers, tops, p_values, ax, space=None,
    interval=None):
    """Adds p_value bars

    The assumes that comparison bars are being introduced for a
    group of error bars, assuming a reference and comparison.

    Parameters
    ----------
    centers : array-like
        The center of the barplot.
    tops : array_like
        The maximum height of the groups to be compared. There must be a
        top for each center supplied.
    p_values : array-like
        The critical values between the groups in the bars. There should be
        one less p value than there are bars.
    ax : matplotlib axis, optional
        The axis where data should be plotted. If none, a new axis instance
        will be created.
    space : float
        The y-distance between each the comparison line.
    interval : float
        The space between the bars

    Returns
    -------
    lines : list
        A list of the lines and text objects which have been plotted

    """

    # Checks the shapes of the inputs
    if not centers.shape == tops.shape:
        raise ValueError('centers and tops must be the same length')
    if not centers.shape[0] == (p_values.shape[0] + 1):
        raise ValueError('there must be a p-value for each center')

    # Deterines the bar locations
    max_hi = tops.max()
    if max_hi < 1:
        fudge = np.power(10, -np.floor(np.log(max_hi)))
    else:
        fudge = np.power(10, -np.ceil(np.log(max_hi)))
    first = np.ceil(max_hi * fudge)/fudge

    correct = np.power(10, -np.log10(fudge))
    # Sets the spacing
    if space is None:
        space = correct/10.
    if interval is None:
        interval = correct/3.

    # Sets the tops of hte bars
    bar_levels = np.arange(0., len(p_values))*interval + first
    # Identifies the center positions for the p text
    p_cents = [(centers[1] - centers[0])/2. + centers[0]]
    if len(centers) > 2:
        for idx, center in enumerate(centers[2:]):
            p_cents.append((center - centers[0])/2. + centers[0])
    # Determines the p text
    p_text = []
    for p in p_values:
        if p < 0.01:
            p_text.append('%1.e' % p)
        else:
            p_text.append('%1.2f' % p)
    lines = []

    # Plots the first comparison
    l1 = []
    l1.append(ax.plot([centers[0]]*2, [tops[0] + space, bar_levels[-1]], 'k-'))
    l1.append(ax.plot([centers[1]]*2, [tops[1] + space, bar_levels[0]], 'k-'))
    l1.append(ax.plot([centers[0], centers[1]], [bar_levels[0]]*2, 'k-'))
    l1.append(ax.text(x=p_cents[0], y=bar_levels[0] + 0.2*space, s=p_text[0],
                      size=11, ha='center'))
    lines.append(l1)

    if len(centers) > 2:
        for idx, p_cent in enumerate(p_cents[:-1]):
            ln = []
            ln.append(ax.plot([centers[idx + 2]]*2,
                              [tops[idx + 2] + space, bar_levels[idx + 1]],
                              'k-'))
            ln.append(ax.plot([centers[0], centers[idx + 2]],
                              [bar_levels[idx + 1]]*2, 'k-'))
            ln.append(ax.text(x=p_cents[idx + 1],
                              y=bar_levels[idx + 1] + 0.2*space,
                              s=p_text[idx+1], size=11, ha='center'))
            lines.append(ln)

    return lines


def get_distance_vectors(dm, df, group, order=None):
    """Extracts the distance information for all samples in a group

    Parameters
    ----------
    dm : skbio DistanceMatrix
        A distance matrix object with the samples corresponding to those in
        the data frame with the metadata.
    df : pandas DataFrame
        A dataframe containing the metadata associated withed with the object
    group : str
        the metadata category being interograted, used to group the data
    order : None, list
        Default is None. The order of groups in group.

    Returns
    -------
    within : list
        A list of the group names for the within group distance
    w_dist : list
        A list of arrays with the within group distances
    between : list
        A list of the group names for between group distance
    b_dist : list
        A list of arrays for between group distance

    """
    # Checks the column is supported
    if group not in df:
        raise ValueError('%s is not a metadata category.' % group)

    # Gets the sample ids associated with the group
    if order is None:
        order = df.groupby(group).groups
    else:
        order = {o: df.groupby(group).groups[o] for o in order}

    dist_ids = dm.ids
    dist_data = dm.data

    within = []
    w_dist = []
    between = []
    b_dist = []
    for id1, o1 in enumerate(order.keys()):

        # Gets the intragroup distance
        within.append(o1)
        w_dist.append(dm.filter(order[o1]).condensed_form())

        # Gets the intergroup distances
        for o2 in order.keys()[(id1+1):]:
            loc1 = [dist_ids.index(i) for i in order[o1]]
            loc2 = [dist_ids.index(i) for i in order[o2]]
            data1 = dist_data[loc1, :]
            # Adds the intragroup distance
            between.append({o1, o2})
            b_dist.append((data1[:, loc2]).flatten())

    return within, w_dist, between, b_dist


def beta_diversity_bars(dm, meta, group, order=None, ref_group=None,
    num_iter=999, p_crit=0.01, p_table=None,
    p_tab_col='Parametric p-value (Bonferroni-corrected)',
    tails="two", ax=None, interval=0.1, width=0.1, show_seperation=True,
    hide_axes=True, bar_props={}):
    """
    Parameters
    ----------
    dm : skbio DistanceMatrix
        A distance matrix object with the samples corresponding to those in
        the data frame with the metadata.
    meta : pandas dataframe
        The metadata associated with the samples in the DataFrame. The samples
        in the DataFrame may be a superset of the samples in the distance
        matrix, `dm`.
    group : str
        the metadata category being interograted, used to group the data
    order : list, optional
        The order of categories in group. If no `order` is specified, all
        categories in group will be used.
    ref_group : str, optional
        The group within `order` to which all other groups should be compared.
        If group is specified, the first group in `order` will be used.
    num_iter : int
        The number of times to run the permanova to test signifigance between
        the groups
    p_crit : float
        The permnova p-value must be less than p-crit to generate a plot.
    p_table : pandas DataFrame
        A comparison listing the groups and their critical values. This should
        be generated by the Qiime Script, `make_distance_boxplots.py`
    p_tab_col : str
        The column containing the p value to be used. Default is
        'Parametric p-value (Bonferroni-corrected)'.
    tails : {'two', 'left', 'right'}
        Describes the way p values should be calculated with reguard to the
        critical value. Critical values are calculated as two-tailed, but
        will be converted to one tailed, so that if `tail` is "left",
        the mean for the ref_group must be smaller than the mean for the
        group being compared.
    ax : matplotlib axis, optional
        The axis where data should be plotted. If none, a new axis instance
        will be created.
    interval : float, optional
        The spacing between the bars
    width : float, optional
        The width of each bars. Should be less than or equal to the interal.
    boxplot_props : dict, optional
        Features of the barchart.

    Returns
    -------
    ax : axes
        A matplotlib axes containing the plotted data
    feats : dict
        A dictionary with features of the plot, including the boxplot handles
        and text objects.

    """

    # Removes any undefined groups
    map_ = meta.groupby(meta[group].apply(_check_strs)
                        ).get_group(True)
    dm = dm.filter(map_.index)

    # Orders the objects
    if order is None:
        order = map_.groupby(group).groups.keys()
    if ref_group is None:
        ref_group = order[0]

    # Determines hte position of the reference group
    ref_loc = order.index(ref_group)

    # Checks for an overall signifigant differences
    perma_res = skbio.stats.distance.permanova(dm, map_, group, num_iter)
    all_p = perma_res['p-value']
    if all_p > p_crit:
        raise ValueError('There is not a signfigant difference at p < %s.'
                         '\np = %1.2f' % (p_crit, all_p))
    if 'p_value' not in bar_props:
        bar_props['p_value'] = all_p

    # Gets the distance vectors
    w_groups, w_dist, b_groups, b_dist = \
        get_distance_vectors(dm, meta, group, order)

    # Adds the distance vector to the means
    dist_bar = np.zeros((len(order)))
    dist_std = np.zeros((len(order)))

    dist_bar[ref_loc] = w_dist[w_groups.index(ref_group)].mean()
    dist_std[ref_loc] = w_dist[w_groups.index(ref_group)].std()

    # Sets up the pvalues
    if p_table is None:
        p_values = p_table
    else:
        p_values = np.array([])
        sub_p = pd.concat([p_table.loc[p_table['Group 1'] == '%s vs. %s'
                                       % (ref_group, ref_group)],
                           p_table.loc[p_table['Group 2'] == '%s vs. %s'
                                       % (ref_group, ref_group)]])

    for idx, group in enumerate(order):
        if group == ref_group:
            continue
        dist_bar[idx] = b_dist[b_groups.index({ref_group, group})].mean()
        dist_std[idx] = b_dist[b_groups.index({ref_group, group})].std()

        if p_table is not None:
            if tails == 'left' and dist_bar[idx] < dist_bar[0]:
                p_values = np.hstack((p_values, np.array([1])))
            elif tails == 'right' and dist_bar[idx] > dist_bar[0]:
                p_values = np.hstack((p_values, np.array([1])))
            elif '%s vs. %s' % (ref_group, group) in sub_p['Group 2'].values:
                p_values = np.hstack((p_values,
                                      sub_p.loc[sub_p['Group 2'] ==
                                                '%s vs. %s'
                                                % (ref_group, group),
                                                p_tab_col]))
            elif '%s vs. %s' % (ref_group, group) in sub_p['Group 1'].values:
                p_values = np.hstack((p_values,
                                      sub_p.loc[sub_p['Group 1'] ==
                                                '%s vs. %s'
                                                % (ref_group, group),
                                                p_tab_col]))
            elif '%s vs. %s' % (group, ref_group) in sub_p['Group 2'].values:
                p_values = np.hstack((p_values,
                                      sub_p.loc[sub_p['Group 2'] ==
                                                '%s vs. %s'
                                                % (group, ref_group),
                                                p_tab_col]))
            elif '%s vs. %s' % (group, ref_group) in sub_p['Group 1'].values:
                p_values = np.hstack((p_values,
                                      sub_p.loc[sub_p['Group 1'] ==
                                                '%s vs. %s'
                                                % (group, ref_group),
                                                p_tab_col]))
            else:
                raise ValueError('%s vs. %s is not a defined group'
                                 % (ref_group, group))

    if 'ytick_size' not in bar_props:
        bar_props['ytick_size'] = 12

    dist_bar = np.array(dist_bar)
    dist_std = np.array(dist_std)

    ax, feats, xpos = barchart(dist_bar,
                               ax=ax,
                               errors=dist_std,
                               interval=interval,
                               width=width,
                               xticklabels=order,
                               **bar_props)
    if p_values is not None and (ref_loc == 0 or len(dist_bar) == 2):
        bars = add_comparison_bars(xpos, dist_bar+dist_std, p_values, ax)
    elif p_values is not None and ref_loc == (len(dist_bar) - 1):
        print 'last'
        bars = add_comparison_bars(xpos[::-1], (dist_bar+dist_std)[::-1],
                                   p_values[::-1], ax)
    elif p_values is not None:
        bars = add_comparison_bars(xpos[ref_loc:],
                                   dist_bar[ref_loc:] + dist_std[ref_loc:],
                                   p_values[(ref_loc):],
                                   ax)
        bars2 = add_comparison_bars(xpos[:(ref_loc+1)][::-1],
                                    dist_bar[:(ref_loc+1)][::-1] +
                                    dist_std[:(ref_loc+1)][::-1],
                                    p_values[:(ref_loc)][::-1],
                                    ax)
        bars.extend(bars2)
    else:
        bars = None

    feats['error_bars'] = bars
    feats['mean_height'] = dist_bar
    feats['height_std'] = dist_std

    if show_seperation:
        ax.plot(np.arange(-0.25, 0.75, 0.01)*(1 + np.floor(len(xpos))/10),
                np.array([0.4725, 0.4775]*50), 'k-', linewidth=7)
        ax.plot(np.arange(-0.25, 0.75, 0.01)*(1 + np.floor(len(xpos))/10),
                np.array([0.4725, 0.4775]*50), 'w-', linewidth=3)
        yticklabels = ax.get_yticks()
        yticklabels[0] = 0
        ax.set_yticklabels(yticklabels, size=bar_props['ytick_size'])

    if hide_axes:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')

    return ax, feats


def split_taxa(taxa, level=7):
    """Splits a greengenes taxonomy string into a dataframe

    Parameters
    ----------
    taxa : list
        a list of greengenes strings
    level: {2, 3, 4, 5, 6, 7}
        the taxonomic level for summary

    Returns
    -------
    splits : array
        the taxonomy data summarized to the specified level. If data is missing
        for a particular level, the string will specify the last known string.
        (i.e. k__kingdom;p__phylum;c__pclass;o__porder;f__family;g__;s__
        outputs with a genus and species value "f. family").

    levels: list
        the name of the taxonomic levels corresponding to columns in the table.
    """

    if level not in {2, 3, 4, 5, 6, 7}:
        raise ValueError('%r is not a supported taxonomic level.' % level)

    # Prealocates a holding object
    splits = np.zeros((1, level))

    # Sets up the names of phylogenetic levels
    levels = ['kingdom', 'phylum', 'p_class', 'p_order', 'family', 'genus',
              'species']

    # Loops through the list of taxa
    for taxon in taxa:
        if taxon == 'Unassigned':
            clean = ['Unassigned']*level
        else:
            # Splits the greengenes string
            rough = taxon.split(';')
            while len(rough) < level:
                rough.append(rough[-1])

            # Watches the last space in the data
            last = rough[0].replace('__', '. ')
            clean = []

            # Cleans the taxa
            for t in rough[:level]:
                val = t.split('__')[1].replace('_', ' ').replace('-', ' ')
                if val == '':
                    clean.append(last)
                elif '[' in val:
                    v_corr = 'cont. %s' % val.replace('[', '').replace(']', '')
                    clean.append(v_corr)
                    last = '%s. %s' % (t.split('__')[0], v_corr)
                else:
                    clean.append(val)
                    last = t.replace('__', '. ')
        # Adds the information to the taxa watch
        splits = np.vstack((splits, np.array(clean)))

    return splits[1:], levels[:level]


def get_ratio_heatmap(data, ref_pos=None, log=None):
    """Calculates a ratio array

    Parameters
    ----------
    data : array
        The data to be plotted
    ref_pos : int, optional
        The column to be used as the denomenator. If nothing is supplied,
        the ratio will be calculated compared to the arithmatic mean of
        all values
    log : float
        The value to use to caluclate the log of the data. If no value is
        specified, no log will be taken.

    Returns
    -------
    ratio : array
        A numpy array with the the ratio of the data values to the reference
        group.

    """
    if ref_pos is None:
        ref = data.mean(1)
    else:
        ref = data[:, ref_pos]

    ref = ((ref*np.ones((1, data.shape[0]))) *
           np.ones((data.shape[1], 1))).transpose()

    if log is None:
        ratio = data / ref
    else:
        ratio = np.log(data / ref) / np.log(log)

    return ratio


def heatmap(data, ax=None, xticklabels=None, yticklabels=None,
    cmap='RdBu_r', clims=None, xfont_angle=0, xfont_align='center',
    yfont_angle=0, yfont_align='right', xfont_size=12, yfont_size=12,
    cbar_size=11):
    """Wraps matplotlib's heatmap and formats the colorbar

    Parameters
    ----------
    data : array
        The data to be plotted in the heat map
   ax : matplotlib axis, optional
        The axis where data should be plotted. If none, a new axis instance
        will be created.
    xticklabels : list, optional
        The labels to be used for the x-axis, to describe the groups of
        data.
    yticklabels : list, optional
        The labels for each group on the y-axis.
    cmap : str, optional
        The colormap name to use for plotting.
    clims : list, optional
        The limits on the colormap
    xfont_angle : float, optional
        The angle in degrees for the x tick label text.
    xfont_align : {'left', 'right', 'center'}, optional
        The horizonal alignment of the x tick label text. For rotated text,
        an alignment or 'right' is recommended.
    yfont_angle : float, optional
        The angle in degrees for the x tick label text.
    yfont_align : {'left', 'right', 'center'}, optional
        The horizonal alignment of the x tick label text. For rotated text,
        an alignment or 'right' is recommended.
    xfont_size : int, optional
        Default is 12. The size of the x tick labels
    yfont_size : int, optional
        Default is 12. The size of the y tick labels
    cbar_size : int, optional
        Default is 11. The size of the labels on the colorbar.

    Returns
    -------
    ax: matplotlib axis
        The axis with the data plotted.
    cbar : matplotlib colorbar
        The colorbar instance on the plot

    """
    # Checks the shape of the data is sane
    mat_shape = data.shape
    if xticklabels is not None and not len(xticklabels) == mat_shape[1]:
        raise ValueError('There must be a label for each column in data.')
    if yticklabels is not None and not len(yticklabels) == mat_shape[0]:
        raise ValueError('There must be a label for each row in data.')

    # Gets the axis
    if ax is None:
        ax = plt.axes()
    # Gets the associated figure
    fig = ax.get_figure()
    # Plots the data
    im = ax.pcolor(data, cmap=cmap)
    if clims is not None:
        im.set_clim(clims)
    cbar = fig.colorbar(im, ax=ax)

    # Formats the x-axis
    xlim = [0, mat_shape[1]]
    ax.set_xlim(xlim)
    ax.set_xticks(np.arange(0.5, mat_shape[1], 1))
    if xticklabels is None:
        xticklabels = ['']*len(ax.get_xticks())
    ax.set_xticklabels(xticklabels, rotation=xfont_angle, ha=xfont_align,
                       size=xfont_size)
    for tic in ax.xaxis.get_major_ticks():
        tic.tick1On = False
        tic.tick2On = False

    # Formats the y axis
    ylim = [mat_shape[0], 0]
    ax.set_ylim(ylim)
    ax.set_yticks(np.arange(0.5, mat_shape[0], 1))
    if yticklabels is None:
        yticklabels = ['']*len(ax.get_xticks())
    ax.set_yticklabels(yticklabels, rotation=yfont_angle, ha=yfont_align,
                       size=yfont_size)
    for tic in ax.yaxis.get_major_ticks():
        tic.tick1On = False
        tic.tick2On = False

    # Formats the colorbar ylabels
    new_labels = ['%s' % t.get_text() for t in cbar.ax.get_yticklabels()]
    new_labels = [t.replace('$', '') for t in new_labels]
    cbar.ax.set_yticklabels(new_labels, size=cbar_size)
    cbar.outline.set_color('none')

    return ax, cbar


def make_duel_heatmaps(gs, order=None, axes=None, **kwargs):
    """Creates side by side abundance and log ratio heatmaps"""

    # Handles the keyword arguments
    kwds = {'p_column': 'Bonferroni_P',
            'p_crit': 0.05,
            'mode': 'RAW',
            'ref': None,
            '_ref_loc': None,
            'ratio_base': np.e,
            'cmap1': 'Greens',
            'cmap2': 'RdBu_r',
            'clims1': None,
            'clims2': [-2, 2],
            'label': 'INDEX',
            'sort_by_taxa': True,
            'xfont_angle': 0,
            'xfont_align': 'center',
            'yfont_angle': 0,
            'yfont_align': 'right',
            'xfont_size': 9,
            'yfont_size': 9,
            'cbar_size': 11}
    for key, value in kwargs.iteritems():
        if key not in kwds:
            raise ValueError('%s is not a supported keyword argument.' % key)
        else:
            kwds[key] = value

    # Creates axes if not specified
    if axes is None:
        ax1 = plt.subplot(1, 2, 1)
        ax2 = plt.subplot(1, 2, 2)
    else:
        ax1 = axes[0]
        ax2 = axes[1]

    # Checks the order is sane
    if order is None:
        order = list(gs.columns[4:-1].values)
    if kwds['ref'] is not None:
        kwds['_ref_loc'] = order.index[kwds['ref']]

    # Sets up the group labels
    g_labels = [o.replace('_mean', '').replace('_', ' ')
                for o in order]

    # Gets the signifigant data frame
    crit = gs[kwds['p_column']] < kwds['p_crit']
    if not crit.any():
        raise ValueError('There are no signifigant taxa')
    sig = gs.loc[crit]

    # Adds taxa to the signifigant dataframe
    taxa, levels = split_taxa(sig.taxonomy.values, 7)
    sig = sig.join(pd.DataFrame(taxa, columns=levels, index=sig.index))

    # Creates an index column
    if kwds['label'] == 'INDEX':
        sig['label'] = sig['genus'] + '(' + \
            sig['OTU'].apply(lambda x: str((x))) + ')'
    else:
        sig['label'] = sig[kwds['label']]

    # If approprate orders the tables taxonomically
    if kwds['sort_by_taxa']:
        sig = sig.sort(['phylum', 'p_class', 'p_order', 'family',
                        'genus', 'species'])

    # Gets the raw data values
    raw = sig[order].values
    ratio = get_ratio_heatmap(sig[order].values, kwds['_ref_loc'],
                              kwds['ratio_base'])
    ratio[np.isinf(ratio)] = np.nan

    # Plots the data on a heatmap
    ax1, cbar1 = heatmap(data=raw,
                         ax=ax1,
                         xticklabels=g_labels,
                         yticklabels=sig['label'].values,
                         cmap=kwds['cmap1'],
                         clims=kwds['clims1'],
                         xfont_angle=kwds['xfont_angle'],
                         xfont_align=kwds['xfont_align'],
                         yfont_angle=kwds['yfont_angle'],
                         yfont_align=kwds['yfont_align'],
                         xfont_size=kwds['xfont_size'],
                         yfont_size=kwds['yfont_size'],
                         cbar_size=kwds['cbar_size'])
    ax2, cbar2 = heatmap(data=ratio,
                         ax=ax2,
                         xticklabels=g_labels,
                         yticklabels=['']*len(sig['label'].values),
                         cmap=kwds['cmap2'],
                         clims=kwds['clims2'],
                         xfont_angle=kwds['xfont_angle'],
                         xfont_align=kwds['xfont_align'],
                         yfont_angle=kwds['yfont_angle'],
                         yfont_align=kwds['yfont_align'],
                         xfont_size=kwds['xfont_size'],
                         yfont_size=kwds['yfont_size'],
                         cbar_size=kwds['cbar_size'])
    ax1.set_xlabel('Raw Data', size=15)
    ax2.set_xlabel('Ratio Data', size=15)

    feats = {'raw_data': raw,
             'ratio_data': ratio,
             'row_labels': sig['label'].values,
             'col_labels': g_labels,
             'cbar1': cbar1,
             'cbar2': cbar2}

    return [ax1, ax2], feats

# def pretty_boxplot(grouped, order, cat, **kwargs):
#     """
#     Creates a more attractive boxplot than pandas

#     Parameters
#     ----------
#     grouped : pandas grouped object
#         the dataframe grouped by the category of interest

#     order : {list, None}
#         if category order matters, this argument sets the order. Otherwise,
#         default sorting is used.

#     cat : str
#         the dataframe column header being plotted. Ideally, this should be
#         continous data.

#     Returns
#     -------
#     fig : matplotlib figure
#         the figure object displaying the graph

#     features : dictionary
#         a dictionary keyed to the axis, boxplot features, kruskal-wallis test
#         statitics, numerical test objects, and p-value text.

#     Other Parameters
#     ----------------
#     axis_dims : {tuple}
#         the portion of the figure covered by the axis. The four element tuple
#         provides the positions of the (Left, Bottom, Width, Height) for the
#         axis.

#     fig_dims : {None, tuple}
#         the size of the figure, in inches

#     interval : {0.5, float}
#         Spacing between each boxplot

#     notch : {True, False}
#         Have the boxplot display the 95% confidence interval notch

#     show_n : {True, False}
#         display group sizes on the figure

#     n_xs : {None, list}
#         the x-positions to display counts. None will place the counts at the
#         x-ticks.

#     n_y : {None, float}
#         the y-position to display the counts. None will place the counts at the
#         bottom of the axis.

#     n_size : {12, int}
#         the text size for the displayed counts

#     show_p : {True, False}
#         display the p value on the axis

#     p_x : {None, float}
#         the x-position for the p-value string. None will place it to the far
#         right

#     p_y : {None, float}
#         the y-position for the p-value string. None will place the p text at
#         the top of the figure.

#     p-size : {12, int}
#         the font size for the displayed p-value

#     title : {'', str}
#         title for the figure

#     title_size : {15, int}
#         the font size for the displayed title

#     xtick_size : {12, int}
#         the text size for xtick labels

#     xfont_align : {'left', 'center', 'right'}
#         the alignment for the xtick labels. Right alignment is recommend with
#         rotated text

#     xfont_angle : {0, int}
#         the angle for xtick labels.

#     xlabel : {'', str}
#         x-axis label

#     xlabel_size : {15, int}
#         the text size for the xlabel

#     show_xgrid : {False, True}
#         display vertical grid lines on the axis

#     ylims : {list, tuple, array}
#         Upper and lower limits for the y-axis

#     ytick_size : {12, int}
#         the text size for ytick labels

#     ylabel : {'', str}
#         y-axis label

#     ylabel_size : {15, int}
#         the font size for the y axis label

#     show_ygrid : {True, False}
#         display horizontal grid lines on the axis

#     Also See
#     --------
#     boxplot : maplotlib.pyplot.boxplot

#     kruskal-wallis : scipy.stats.kruskal
#     """

#     # Handles keyword arguemnts
#     keywords = {'fig_dims': None,
#                 'axis_dims': (0.1, 0.1, 0.8, 0.8),
#                 'notch': True,
#                 'interval': 0.5,
#                 'show_p': True,
#                 'show_n': True,
#                 'title': '',
#                 'tick_names': None,
#                 'xlabel': '',
#                 'ylabel': '',
#                 'ylims': None,
#                 'p_size': 12,
#                 'p_x': None,
#                 'p_y': None,
#                 'n_xs': None,
#                 'n_y': None,
#                 'n_size': 12,
#                 'title_size': 15,
#                 'xtick_size': 12,
#                 'xlabel_size': 15,
#                 'xfont_align': 'center',
#                 'xfont_angle': 0,
#                 'ytick_size': 12,
#                 'ylabel_size': 15,
#                 'show_xgrid': False,
#                 'show_ygrid': True}

#     for key, val in kwargs.iteritems():
#         if key in keywords:
#             keywords[key] = val
#         else:
#             raise ValueError('%s is not an input for plot_with_dist.' % key)

#     # Prealocates a dictionary of plotting features
#     features = {}

#     # Calculates plotting features
#     if order is None:
#         order = grouped.groups.keys()

#     num_cats = len(order)
#     xlim = [-keywords['interval']/2,
#             keywords['interval']*(num_cats-1)+keywords['interval']/2]

#     # Creates the figure
#     fig = plt.figure()
#     if keywords['fig_dims'] is not None:
#         fig.set_size_inches(keywords['fig_dims'])
#     ax = fig.add_axes(keywords['axis_dims'])

#     # Sets up plotting constants
#     ticks = arange(0, keywords['interval']*num_cats, keywords['interval'])
#     names = []
#     values = []
#     counts = []

#     # Appends the information to the features
#     features['fig'] = fig
#     features['axis'] = ax

#     for idx, group in enumerate(order):
#         tick = [ticks[idx]]
#         value = grouped.get_group(group)[cat].values
#         values.append(value)
#         counts.append(len(value))
#         bp = ax.boxplot(value,
#                         positions=tick,
#                         notch=keywords['notch'])

#         # Handles group naming
#         if keywords['tick_names'] is not None:
#             names.append(keywords['tick_names'][group])
#         else:
#             names.append(group)

#     # Adds the boxplot to the features dictionary
#     features['boxplot'] = bp

#     # Calculates the test statitics
#     (features['h'], features['p']) = kruskal(*values)

#     # Sets up the axis properties
#     ax.set_xlim(xlim)
#     if keywords['ylims'] is not None:
#         ax.set_ylim(keywords['ylims'])
#     ax.set_xticks(ticks)
#     ax.set_xticklabels(names,
#                        rotation=keywords['xfont_angle'],
#                        horizontalalignment=keywords['xfont_align'],
#                        size=keywords['xtick_size'])
#     ax.set_yticklabels(ax.get_yticks(), size=keywords['ytick_size'])
#     ax.set_xlabel(keywords['xlabel'], size=keywords['xlabel_size'])
#     ax.set_ylabel(keywords['ylabel'], size=keywords['ylabel_size'])

#     # Adds grid, if appropriate
#     if keywords['show_xgrid']:
#         ax.xaxis.grid()
#     if keywords['show_ygrid']:
#         ax.yaxis.grid()

#     # Adds the title, if appropriate
#     ax.set_title(keywords['title'], size=keywords['title_size'])

#     text = []
#     # Adds the group counts if desired
#     if keywords['show_n']:
#         ylim = ax.get_ylim()
#         if keywords['n_xs'] is None:
#             keywords['n_xs'] = ticks
#         if keywords['n_y'] is None:
#             keywords['n_y'] = ylim[0]+(ylim[1]-ylim[0])*0.025

#         for idx, count in enumerate(counts):
#             text.append(ax.text(x=keywords['n_xs'][idx],
#                                 y=keywords['n_y'],
#                                 s='(%i)' % count,
#                                 horizontalalignment='center',
#                                 size=keywords['n_size']))

#     # Adds a p-value to the plot if desired
#     if keywords['show_p']:
#         # Gets the positions if necessary
#         if keywords['p_x'] is None:
#             xticks = ax.get_xticks()
#             x = ax.get_xlim()[1] - (xticks[-1] - xticks[-2])/20
#         else:
#             x = keywords['p_x']

#         if keywords['p_y'] is None:
#             yticks = ax.get_yticks()
#             y = yticks[-2]+(yticks[-1]-yticks[-2])/2
#         else:
#             y = keywords['p_y']

#         # Adds the text
#         if features['p'] >= 0.005:
#             p_text = ax.text(s='p = %1.2f' % features['p'],
#                              x=x, y=y,
#                              size=keywords['p_size'],
#                              horizontalalignment='right')
#         else:
#             p_text = ax.text(s='p = %1.1e' % features['p'],
#                              x=x, y=y,
#                              size=keywords['p_size'],
#                              horizontalalignment='right')
#         features['p_text'] = p_text

#     fig = plt.gcf()

#     return fig, features
