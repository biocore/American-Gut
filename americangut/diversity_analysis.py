from __future__ import division
from os import mkdir
from os.path import exists

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import skbio

from scipy.stats import kruskal
from skbio.stats.power import _check_strs
from statsmodels.sandbox.stats.multicomp import multipletests

__author__ = "Justine Debelius"
__copyright__ = "Copyright 2015, The American Gut Project"
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
        the data frame without an index column
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


def boxplot(vecs, ax=None, notch=True, interval=0.5, boxplot_props={},
    show_counts=True, **kwargs):
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
    show_counts : bool, optional
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
        The font size for the critical value text
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

    # Determines the plotting locations
    num_cats = len(vecs)
    xlim = [-interval/2,
            interval*(num_cats-1)+interval/2]

    # Sets up the plotting constants
    ticks = np.arange(0, interval*num_cats, interval)
    counts = []

    # Loops through the data
    for idx, vec in enumerate(vecs):
        # Gets vector characteristics
        tick = [ticks[idx]]
        counts.append(len(vec))
        # Plots the data
        ax.boxplot(vec,
                   positions=tick,
                   notch=notch,
                   **boxplot_props)

    # Sets up axis formatting
    if show_counts:
        kwargs['counts'] = counts
    if 'xlim' not in kwargs:
        kwargs['xlim'] = xlim
    if 'xticks' not in kwargs:
        kwargs['xticks'] = ticks

    _format_axis(ax, **kwargs)

    return ax


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

    """

    grouped = meta.groupby(group)

    # Sets up the plotting order
    if order is None:
        order = grouped.groups.keys()

    # Gets the data vectors
    vecs = [grouped.get_group(g)[cat].values for g in order]

    # Formats the axis, if not already done
    if 'xticklabels' not in boxplot_props:
        boxplot_props['xticklabels'] = [g.split('(')[0] for g in order]
    if 'show_xticks' not in boxplot_props:
        boxplot_props['show_xticks'] = False
    if 'show_ygrid' not in boxplot_props:
        boxplot_props['show_ygrid'] = True

    # Calculates the p value
    h, p = kruskal(*vecs)

    # Sets the boxplot properties
    ax = boxplot(vecs=vecs, ax=ax, p_value=p, **boxplot_props)

    return ax


def post_hoc_pandas(meta, group, cat, order=None, correct=None,
    show_stats=True):
    """Preforms an post-hoc comparison between two groups

    Parameters
    ----------
    meta : pandas DataFrame
        the metadata object for the data
    group : str
        the metadata category being interograted
    cat : str
        the name of the column with the result
    order : None, list, optional
        Default is None. The order of groups in the category.
    correct : None, str, optional
        Method for multiple hypothesis correction using
        `statsmodels.sandbox.stats.multicomp.multipletests`. Methods you're
        likely to use are `bonferroni` and `fdr_bh`.
    show_stats : bool, optional
        When `show_stats` is True, a summary of each group will be displayed
        along with the p values. 

    Returns
    -------
    post_hoc : dataframe
        `post_hoc` summarizes the results of the post-hoc test. It includes
        statitics about each distribution, as well as the comparison matrix
        of p-values.

    """

    # Groups the data
    grouped = meta.groupby(group)

    # Gets the order
    if order is None:
        order = grouped.groups.keys()

    # Prealocates an output frame
    if show_stats:
        stats = pd.DataFrame({'Counts': grouped[cat].count(),
                              'Mean':  grouped[cat].mean(),
                              'Stdv': grouped[cat].std(),
                              'Median': grouped[cat].median()})

    # Preforms ad-hoc comparisons
    comparison = {}

    for pos, g1_name in enumerate(order[:-1]):
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
    if show_stats:
        post_hoc = stats.join(compare[order[:-1]])
    else:
        post_hoc = compare[order[:-1]]
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


def barchart(height, interval=0.5, width=0.4, ax=None, errors=None,
    colormap=None, match_colors=True, elinewidth=2, ecapwidth=2, offset=0,
    **kwargs):
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
    elinewidth : int
        The weight in points of the errorbar on top of the plot.
    ecapwidth : int
        The weight in points of the line over the errorbar.
    offset : float
        The space on the axis for which the data should be offset. This allows
        plotting multiple barcharts on the same axis with spacing, or plotting
        the bar chart away from the origin.
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
        The font size for the critical value text
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

    # Sets the colormap and egdecolor, if necessary
    if colormap is None:
        colormap = np.array([[0.5, 0.5, 0.5]]*len(height))
    elif isinstance(colormap, str):
        colormap = segment_colormap(colormap, len(height))

    if match_colors:
        edgecolors = colormap
    else:
        edgecolors = np.array([[0, 0, 0]]*len(height))

    # Gets an axis instance
    if ax is None:
        ax = plt.axes()

    # Gets the xposition, for errorbars
    xpos = np.arange(0, len(height))*interval + interval/2 + offset

    xleft = np.arange(0, len(height)) * interval + (interval - width)/2 + offset
    xlims = [0, len(height)*interval]

    if 'xlim' not in kwargs or kwargs['xlim'] is None:
        kwargs['xlim'] = xlims

    # Plots the errorbars
    if errors is not None:
        e_bars = []
        for idx, err in enumerate(errors):
            eb = ax.errorbar(x=xpos[idx],
                             y=height[idx],
                             yerr=err,
                             fmt='none',
                             ecolor=edgecolors[idx, :],
                             elinewidth=elinewidth,
                             capthick=ecapwidth)
        e_bars.append(eb)

    # Plots the bars
    bars = ax.bar(xleft, height, width=width)
    for idx, bar in enumerate(bars):
        bar.set_facecolor(colormap[idx, :])
        bar.set_edgecolor(edgecolors[idx, :])

    _format_axis(ax)

    return ax, xpos


def add_comparison_bars(centers, tops, p_values, ax, space=None,
    interval=None, lowest=None, factor=5, label_size=10,
    show_value=True):
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
    lowest : float, optional
        If mutliple sets of comparison bars are being added, this allows
        the user to set the position of the lowest bar in the group
    factor : unsigned int, optional
        The ones-place value to which values will be rounded. For example,
        for a `factor` of 5, a value of 0.12 will be rounded to 0.15 and a
        value of 2.7 will be rounded to 3.0.
    label_size : unsigned int, optional
        The font size for displaying signifigance labels
    show_value : bool, optional
        When True, the signigance bars will display the actual p value.
        Otherwise, p values will be coded as (p <= 0.1: '+', p <= 0.05: '*',
        p <= 0.01: '**', p <= 0.001: '***', p <= 0.0001).

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
    if lowest is None:
        lowest, fudge = _get_bar_height(tops, factor)
    else:
        __, fudge = _get_bar_height(tops, factor)

    correct = np.power(10, -np.log10(fudge))

    # Sets the spacing
    if space is None:
        space = correct/15.
    if interval is None:
        interval = correct/3.

    # Sets the tops of hte bars
    bar_levels = np.arange(0., len(p_values))*interval + lowest

    # Identifies the center positions for the p text
    p_cents = [(centers[1] - centers[0])/2. + centers[0]]

    if len(centers) > 2:
        for idx, center in enumerate(centers[2:]):
            p_cents.append((center - centers[0])/2. + centers[0])
    # Determines the p text
    p_text = []
    for p in p_values:
        if p < 0.01 and show_value:
            p_text.append('%1.e' % p)
        elif show_value:
            p_text.append('%1.2f' % p)
        elif p < 0.0001:
            p_text.append('****')
        elif p < 0.001:
            p_text.append('***')
        elif p < 0.01:
            p_text.append('**')
        elif p < 0.05:
            p_text.append('*')
        elif p < 0.1:
            p_text.append('+')
        else:
            p_text.append('')

    lines = []

    # Plots the first comparison
    l1 = []
    l1.append(ax.plot([centers[0]]*2, [tops[0] + space, bar_levels[-1]], 'k-'))
    l1.append(ax.plot([centers[1]]*2, [tops[1] + space, bar_levels[0]], 'k-'))
    l1.append(ax.plot([centers[0], centers[1]], [bar_levels[0]]*2, 'k-'))
    l1.append(ax.text(x=p_cents[0], y=bar_levels[0] + 0.2*space, s=p_text[0],
                      size=label_size, ha='center'))
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
                              s=p_text[idx+1], size=label_size, ha='center'))
            lines.append(ln)

    return lines


def segment_colormap(cm_name, n_colors, n_pad=None, start=None):
    """Segments a matplotlib colormap into descrete colors

    Parameters
    ----------
    cm_name : str
        The name of the continous matplotlib colormap which will be segmented.
    n_colors : int
        The number of colors needed in the colormap
    n_pad : int, optional
        The number of total colors to have in the colormap. By default, this
        will be one more than the number of colors provided.
    start : int, optional
        The number of colors to skip over before display starts. By default,
        this is `n_pad` - `n_colors`.

    Returns
    -------
    new_map : array
        A segmented array containing the colormap

    """

    # Sets parameters if necessary
    if n_pad is None:
        n_pad = n_colors + 1
    if start is None:
        start = n_pad - n_colors

    # Gets the colormap
    cm = mpl.cm.get_cmap(cm_name)
    # Sets up the new map
    new_map = np.array([list(cm(1. * (i + start) / n_pad)) for i in
                        xrange(n_colors)])

    return new_map


def _get_bar_height(tops, factor=5):
    """Calculates the lowest bar height"""
    max_hi = tops.max()
    # Gets the correct order of magnitdue
    if max_hi < 1:
        fudge = np.power(10, -np.floor(np.log(max_hi)))
    else:
        fudge = np.power(10, -np.ceil(np.log(max_hi)))
    # Gets the correct rounding of the factor
    if int(max_hi*fudge) < 5:
        lowest = np.ceil(max_hi*fudge/factor)/(fudge)*factor
    else:
        lowest = np.ceil(max_hi*fudge*10/factor)/(fudge*10)*factor

    return lowest, fudge


def _get_p_value(sub_p, ref_group, group, p_tab_col):
    """Determines the comparison value within a post-hoc table"""
    if '%s vs. %s' % (ref_group, group) in sub_p['Group 2'].values:
        p_value = sub_p.loc[sub_p['Group 2'] == '%s vs. %s'
                            % (ref_group, group), p_tab_col].values[0]
    elif '%s vs. %s' % (ref_group, group) in sub_p['Group 1'].values:
        p_value = sub_p.loc[sub_p['Group 1'] == '%s vs. %s'
                            % (ref_group, group), p_tab_col].values[0]
    elif '%s vs. %s' % (group, ref_group) in sub_p['Group 2'].values:
        p_value = sub_p.loc[sub_p['Group 2'] == '%s vs. %s'
                            % (group, ref_group), p_tab_col].values[0]
    elif '%s vs. %s' % (group, ref_group) in sub_p['Group 1'].values:
        p_value = sub_p.loc[sub_p['Group 1'] == '%s vs. %s'
                            % (group, ref_group), p_tab_col].values[0]
    else:
        raise ValueError('%s vs. %s is not a defined group'
                         % (ref_group, group))

    return p_value


def _correct_p_value(tail, p_value, ref_val, current_val):
    """Determines if p-values should be tail corrected"""
    if tail and ref_val > current_val:
        return 1
    else:
        return p_value


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
        order = list(df.groupby(group).groups)
    ordered_ids = {o: df.groupby(group).groups[o] for o in order}

    # Gets the data
    dist_ids = np.array(list(dm.ids))
    dist_data = dm.data

    # Prealocates objects for return
    within = []
    w_dist = []
    between = []
    b_dist = []

    # Loops through the groups
    for id1, o1 in enumerate(order):
        # Gets the intragroup distance
        within.append(o1)
        w_dist.append(dm.filter(ordered_ids[o1]).condensed_form())
        for o2 in order[(id1+1):]:
            loc1 = np.array([dist_ids == i for i in ordered_ids[o1]]).any(0)
            loc2 = np.array([dist_ids == i for i in ordered_ids[o2]]).any(0)
            data1 = dist_data[loc1, :]
            # Adds the intragroup distance
            between.append({o1, o2})
            b_dist.append((data1[:, loc2]).flatten())

    return within, w_dist, between, b_dist


def beta_diversity_bars(dm, meta, group, order=None, ref_groups=None,
    num_iter=999, p_crit=0.01, p_table=None,
    p_tab_col='Parametric p-value (Bonferroni-corrected)',
    ref_less=True, ax=None, interval=0.1, width=0.1,
    show_seperation=True, colormap=None, match_colors=True,
    elinewidth=2, ecapwidth=2, show_p=False, lowest=None,
    sep_size=0.035, label_size=12, show_p_value=False, **kwargs):
    """Creates a barchart of the beta diversity distances

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
    ref_less : bool, optional
        Indicates the reference group should be the smallest bar among the
        others, and p_values are caluclated accordingly. That is, even if a
        difference is signifigant, if the reference group is greater than the
        group being examined, the returned p-value will be 1.
    ax : matplotlib axis, optional
        The axis where data should be plotted. If none, a new axis instance
        will be created.
    interval : float, optional
        The spacing between the bars
    width : float, optional
        The width of each bars. Should be less than or equal to the interal.
    show_seperation : bool, optional
        Inserts a set of wavy lines to suggest the bottom of the distance
        axes is at 0, and not whatever convenient x-lim is used.
    colormap : array, optional
        The colors to be used for the barchart. This must be an nx3 or nx4
        array.
    match_color: bool, optional
        When True, the edges and errorbar lines are the same color as the
        bar fills. When False, these are outlined in black.
    elinewidth : int, optional
        The linethickness for the errorbar
    ecapwidth: int, optional
        The width of the line on top of the errorbar
    show_p: bool, optional
        When True, displays the overall p value on the plot
    lowest : float, optional
        The position for the lowest comparison bar, if showing comparison bars
        on the plot. If None, this will be calculated.
    sep_size : float, optional
        When `show_seperation` is true, this is used to set the size of the
        seperation lines.
    label_size : int, optional
        Sets the size of the signfigance labels
    show_p_value : bool, optional
        When True, the signigance bars will display the actual p value.
        Otherwise, p values will be coded as (p <= 0.1: '+', p <= 0.05: '*',
        p <= 0.01: '**', p <= 0.001: '***', p <= 0.0001).
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
    """

    # Removes any undefined groups
    map_ = meta.groupby(meta[group].apply(_check_strs)
                        ).get_group(True)
    dm = dm.filter(map_.index)

    # Orders the objects
    if order is None:
        order = map_.groupby(group).groups.keys().sort()

    # Gets the reference groups
    if isinstance(ref_groups, str):
        ref_groups = [ref_groups]
    if ref_groups is None:
        ref_groups = order

    # Checks for an overall signifigant differences
    if show_p:
        perma_res = skbio.stats.distance.permanova(dm, map_, group, num_iter)
        kwargs['p_value'] = perma_res['p-value']

    # Gets the distance vectors
    w_groups, w_dist, b_groups, b_dist = \
        get_distance_vectors(dm, meta, group, order)

    bar_counter = 0

    # Loops through the data to make boxplots compared to the reference group
    for id1, ref_group in enumerate(ref_groups):
        # Determines the offset
        offset = (id1 * len(order) + bar_counter)*width
        # Determines the position of the reference group
        ref_loc = order.index(ref_group)

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

        for id2, group in enumerate(order):
            if group == ref_group:
                continue
            # Gets the distance vector
            dist_bar[id2] = b_dist[b_groups.index({ref_group, group})].mean()
            dist_std[id2] = b_dist[b_groups.index({ref_group, group})].std()
            p_value = _get_p_value(sub_p, ref_group, group, p_tab_col)
            p_values = np.hstack((p_values,
                                  _correct_p_value(ref_less, p_value,
                                                   dist_bar[ref_loc],
                                                   dist_bar[id2])))

        dist_bar = np.array(dist_bar)
        dist_std = np.array(dist_std)

        # Creates the boxplot
        ax, xpos = barchart(dist_bar,
                            ax=ax,
                            errors=dist_std,
                            interval=interval,
                            width=width,
                            offset=offset,
                            xticklabels=order,
                            colormap=colormap,
                            match_colors=match_colors,
                            elinewidth=elinewidth,
                            ecapwidth=ecapwidth)
        # Watches the xlimits
        if id1 == 0:
            xlim = [xpos[0] - interval*0.75]

        # Gets the critical lines to display on the figure
        if lowest is None:
            lowest, __ = _get_bar_height(dist_bar+dist_std)
        if p_values is not None and (ref_loc == 0 or len(dist_bar) == 2):
            bars = add_comparison_bars(xpos,
                                       dist_bar+dist_std,
                                       p_values,
                                       ax,
                                       lowest=lowest,
                                       label_size=label_size,
                                       show_value=show_p_value)
        elif p_values is not None and ref_loc == (len(dist_bar) - 1):
            bars = add_comparison_bars(xpos[::-1],
                                       (dist_bar+dist_std)[::-1],
                                       p_values[::-1],
                                       ax,
                                       lowest=lowest,
                                       label_size=label_size,
                                       show_value=show_p_value)
        elif p_values is not None:
            bars = add_comparison_bars(xpos[ref_loc:],
                                       dist_bar[ref_loc:] + dist_std[ref_loc:],
                                       p_values[(ref_loc):],
                                       ax, lowest=lowest,
                                       label_size=label_size,
                                       show_value=show_p_value)
            bars2 = add_comparison_bars(xpos[:(ref_loc+1)][::-1],
                                        dist_bar[:(ref_loc+1)][::-1] +
                                        dist_std[:(ref_loc+1)][::-1],
                                        p_values[:(ref_loc)][::-1],
                                        ax, lowest=lowest,
                                        label_size=label_size,
                                        show_value=show_p_value)
            bars.extend(bars2)

        # Advances the bar count
        bar_counter = bar_counter + 1
    xlim.append(xpos[-1] + interval*0.75)

    # Sets up axis formatting defaults, if not supplied
    if 'show_frame' not in kwargs:
        kwargs['show_frame'] = False
    if 'ytick_size' not in kwargs:
        kwargs['ytick_size'] = 12
    if 'show_xticks' not in kwargs:
        kwargs['show_xticks'] = False
    if 'xticks' not in kwargs:
        # Calculates the width of a group of bars
        width = (xpos[-1] - xpos[0])+interval*2
        left = xlim[0] + width/2 - interval*0.25
        kwargs['xticks'] = np.arange(left, left+width*len(ref_groups), width)

    if 'xticklabels' not in kwargs:
        kwargs['xticklabels'] = [r.split('(')[0].replace('_', ' ') for r
                                 in ref_groups]
    if 'xtick_size' not in kwargs:
        kwargs['xtick_size'] = 12
    if 'xlabel' not in kwargs:
        kwargs['xlabel'] = 'Reference Group'
    if 'xlim' not in kwargs:
        kwargs['xlim'] = xlim

    # Formats the axis, to make it pretty
    _format_axis(ax, **kwargs)

    # Shows the seperation line if desired
    if show_seperation:
        xlim = ax.get_xlim()
        lower_y, upper_y = ax.get_ylim()
        tick_dist = (upper_y - lower_y)/5.
        ax.plot(np.arange(-0.25, 0.75, 0.01) *
                (1 + np.floor(len(kwargs['xticks']))),
                np.array([-tick_dist*sep_size, tick_dist*sep_size]*50) +
                lower_y + tick_dist*0.4,
                'k-', linewidth=7)
        ax.plot(np.arange(-0.25, 0.75, 0.01) *
                (1 + np.floor(len(kwargs['xticks']))),
                np.array([-tick_dist*sep_size, tick_dist*sep_size]*50) +
                lower_y + tick_dist*0.4,
                'w-', linewidth=3)
    else:
        xlim = ax.get_xlim()

    # Checks the format, just to be sure... Formats the axis, to make it pretty
    ax.set_xlim(xlim)

    # Formats the yaxis to show seperation
    if show_seperation:
        yticklabels = ax.get_yticks()
        yticklabels[0] = 0
        ax.set_yticklabels(yticklabels, size=kwargs['ytick_size'])

    # Adds a legend to the top to distinguish groups
    # Determines the horizontal ordering for the groups
    num_groups = len(order)
    num_leg_cols = int(np.ceil(float(num_groups)/2))
    leg_order = np.hstack([np.array([i, i+num_leg_cols])
                           for i in xrange(num_leg_cols)])[:num_groups]
    # Gets the patches and labels
    patches = np.array(ax.patches[:num_groups])[leg_order]
    labels = np.array([l.split('(')[0].replace('_', ' ')
                       for l in order])[leg_order]
    # Adds the legend
    ax.legend(patches, labels,
              loc="upper center",
              ncol=num_leg_cols,
              fontsize=12,
              frameon=False)

    return ax


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
    splits = np.zeros((1, (level)))

    # Sets up the names of phylogenetic levels
    levels = ['kingdom', 'phylum', 'p_class', 'p_order', 'family', 'genus',
              'species']

    # Loops through the list of taxa
    for taxon in taxa:
        if taxon == 'Unassigned':
            clean = ['Unassigned']*level
        else:
            # Splits the greengenes string
            rough = [v.strip() for v in taxon.split(';')]
            # Watches the last space in the data
            try:
                last = rough[level - 1].replace('__', '. ')
            except:
                raise ValueError('There are %i levels in your taxa string and'
                                 ' you would like to look at %i levels.'
                                 % (len(rough), level))
            clean = []

            # Cleans the taxa
            for t in rough[:(level)]:
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


def heatmap(data, ax=None,  cmap='RdBu_r', clims=None, cbar_size=11, **kwargs):
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
    if ('xticklabels' in kwargs and kwargs['xticklabels'] is not None and not
        len(kwargs['xticklabels']) == mat_shape[1]):
        raise ValueError('There must be a label for each column in '
                         'data.')

    if ('yticklabels' in kwargs and kwargs['yticklabels'] is not None and not
        len(kwargs['yticklabels']) == mat_shape[0]):
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

    # Gets the axis limits
    if 'xlim' not in kwargs:
        kwargs['xlim'] = [0, mat_shape[1]]
    if 'ylim' not in kwargs:
        kwargs['ylim'] = [mat_shape[0], 0]
    if 'xticks' not in kwargs:
        kwargs['xticks'] = np.arange(0, mat_shape[1]) + 0.5
    if 'yticks' not in kwargs:
        kwargs['yticks'] = np.arange(0, mat_shape[0]) + 0.5

    # Formats the axis
    _format_axis(ax, **kwargs)

    # Formats the colorbar ylabels
    new_labels = ['%s' % t.get_text() for t in cbar.ax.get_yticklabels()]
    new_labels = [t.replace('$', '') for t in new_labels]
    cbar.ax.set_yticklabels(new_labels, size=cbar_size)
    cbar.outline.set_color('none')

    return ax, cbar


def make_dual_heatmaps(gs, order=None, axes=None, p_column='Bonferroni_P',
    p_crit=0.05, ref=None, ratio_base=np.e, cmap1='Greens',
    cmap2='RdBu_r', clims1=None, clims2=[-2, 2],
    label='INDEX', sort_by_taxa=True, cbar_size=12,
    consistant_size=True, width=11, height=8.5, **kwargs):
    """Creates side by side abundance and log ratio heatmaps

    Parameters
    ----------
    gs : pandas DataFrame
        The results of Qiime's `group_signifigance.py` read into pandas.
    order : list, optional
        The order in which columns from `gs` should be plotted in the heatmap.
        Groups in order should contain the suffix, `'_mean'`, which is
        standard in the group_signifigance output.
    p_column : {'Bonferroni_P', 'FDR_P', 'P'}, optional
        The column name from which signifigance values should be drawn. It is
        recommended that 'Bonferroni_P' or 'FDR_P' be used.
    p_crit : float, optional
        The critical p value. Comparison p values must be less than this
        for the comparisons to be displayed on the heatmap.
    ref : str, optional
        The name of the column which should serve as the refence for the ratio
        heatmap. If none is supplied, the arethmatic mean will be used.
    ratio_base : float, optional
        The logarthemic base for the ratio heatmap. If `None`, then no
        logarthem will be taken.
    cmap1 : str, optional
        A name of a matplotlib colormap instance, used for the display of the
        raw data heatmap.
    cmap2 : str, optional
        A name of a matplotlib colormap instance, used for the display of the
        ratio data heatmap.
    clims1 : list, optional
        The limits on the colormap for the raw heatmap
    clims2 : list, optional
        The limits on the colormap for the ratio data heatmap. If using a two
        color, diversing colormap (i.e `RdBu`), it is suggested that these be
        positive and negative values of the same magnitude (i.e. `[-2, 2]`).
    label : str, optional
        The name of the column to be used as the label. The string, `INDEX`
        indicates that the OTU ID should be used, in conjunction with a
        genus-level taxonomic description.
    sort_by_taxa : bool, optional
        Orders the groups using some taxonomic information, so OTUs are grouped
        alphabetically by phylum, class, order, family, genus, and species.
    cbar_size : int, optional
        Font size for the colorbar text.
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

    Returns
    -------
    ax1, ax2 : matplotlib axes
        The matplotlib axis instances with the raw and ratio data, respectively
    cbar1, cbar2 : matplotlib colorbars
        The colobar instances for the raw and ratio data, respectively.

    """

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
    if ref is not None:
        ref_loc = order.index[ref]
    else:
        ref_loc = ref

    # Sets up the group labels
    g_labels = [o.replace('_mean', '').replace('_', ' ').split('(')[0]
                for o in order]

    # Gets the signifigant data frame
    crit = gs[p_column] < p_crit
    if not crit.any():
        raise ValueError('There are no signifigant taxa')
    sig = gs.loc[crit]

    # Adds taxa to the signifigant dataframe
    taxa, levels = split_taxa(sig.taxonomy.values, 7)
    sig = sig.join(pd.DataFrame(taxa, columns=levels, index=sig.index))

    # Creates an index column
    if label == 'INDEX':
        sig['label'] = sig['genus'] + ' (' + \
            sig['OTU'].apply(lambda x: str((x))) + ')'
    else:
        sig['label'] = sig[label]

    # If approprate orders the tables taxonomically
    if sort_by_taxa:
        sig = sig.sort(['phylum', 'p_class', 'p_order', 'family',
                        'genus', 'species'])

    # Gets the raw data values
    raw = sig[order].values
    ratio = get_ratio_heatmap(sig[order].values, ref_loc, ratio_base)
    ratio[np.isinf(ratio)] = np.nan

    # Plots the data on a heatmap
    ax1, cbar1 = heatmap(data=raw,
                         ax=ax1,
                         cmap=cmap1,
                         clims=clims1,
                         cbar_size=cbar_size,
                         xticklabels=g_labels,
                         yticklabels=sig['label'].values,
                         **kwargs)
    ax2, cbar2 = heatmap(data=ratio,
                         ax=ax2,
                         cmap=cmap2,
                         clims=clims2,
                         cbar_size=cbar_size,
                         xticklabels=g_labels,
                         yticklabels=['']*len(sig['label'].values),
                         **kwargs)
    ax1.set_xlabel('Raw Data', size=15)
    ax2.set_xlabel('Ratio Data', size=15)

    if consistant_size:
        # Checks lower padding needed for the figure. If text is not horizonal,
        # an extra inch of padding is provided.
        xangled = ('xfont_angle' in kwargs and
                   ((isinstance(kwargs['xfont_angle'], int) and not
                    kwargs['xfont_angle'] == 0) or kwargs['xfont_angle'] ==
                    'vertical'))
        if xangled:
            bottom_pad = 2.0
        else:
            bottom_pad = 1.0
        # The left side is padded with 2.5 inches
        left_pad = 2.5
        if len(order) > 6:
            ax_width = 0.25*len(order)
        else:
            ax_width = 0.50*len(order)
        cbar_pad = 0.1
        cbar_width = 0.25
        total_pad = 0.85 + ax_width + left_pad
        ax_height = 0.12*sig.shape[0]

        if ax_height < 2:
            cbar_height = 2
        else:
            cbar_height = ax_height

        if (ax_height + bottom_pad + 0.5) > height:
            height = (ax_height + bottom_pad + 0.5)

        a1p = (left_pad/width, bottom_pad/height,
               ax_width/width, ax_height/height)
        a2p = (total_pad/width, bottom_pad/height,
               ax_width/width, ax_height/height)
        c1p = ((left_pad + ax_width + cbar_pad)/width, bottom_pad/height,
               cbar_width/width, cbar_height/height)
        c2p = ((total_pad + ax_width + cbar_pad)/width, bottom_pad/height,
               cbar_width/width, cbar_height/height)

        fig = ax1.get_figure()
        fig.set_size_inches((width, height))
        ax1.set_position(a1p)
        ax2.set_position(a2p)
        cbar1.ax.set_position(c1p)
        cbar2.ax.set_position(c2p)

    return [ax1, ax2], [cbar1, cbar2]


def _format_axis(ax, **kwargs):
    """Beautifies the axis

    Parameters
    ----------
    ax : matplotlib axis
        The axis where data should be plotted. If none, a new axis instance
        will be created.
    counts: array_like, optional
        A list of the number of samples in each plotting location, for use
        with a barchart or boxplt
    p_value : float, optional
        Default is None. When supplied, the signfigance value will be displayed
        on the plot in the upper right hand corner by default.
    show_frame: bool, optional
        When true, the frame around the axis is displayed. When false, only the
        lower and left axes will be dispalyed.
    show_xticks : bool, optional
        Default is True. Display a small black tick symbol at the top and
        bottom of the graph at each major tick mark. This does not affect
        whether text will be displayed at the same location.
    show_yticks : bool, optional
        Default is True. Display a small black tick symbol at the top and
        bottom of the graph at each major tick mark. This does not affect
        whether text will be displayed at the same location.
    show_xgrid: bool, optional
        Default is False. Adds vertical grid lines at each major x-tick.
    show_ygrid: bool, optional
        Default is True. Adds horizonal grid lines at each major y-tick.
    title: str, optional
        The title to be placed on the graph.
    xlim : list, optional
        The limits for the x-axis
    ylim : list, optional
        The limits for the y-axis.
    xlabel : str
        The label text describing the contents of the x-axis.
    ylabel : str
        The label text for the y-axis. Every time you leave off appropriate
        labels and units, a science grad student grading lab reports cries
        another bitter tear into their bottle of craft beer.
    xticklabels : None, "text", list, optional
        The labels to be displayed at the x-tick positions. If None, no
        labels will be shown. If "text", the numeric value of each tick
        will be displayed according to the mapping given in `xtick_format`.
        Otherwise, the listed values will be used.
        Text labels will be formatted using the rotation supplied by
        `xfont_angle` and, the alignment from `xfont_align` and the font size
        given by `xfont_size`. It is recommended that rotated labels are
        right aligned.
    yticklabels : None, "text", list, optional
        The labels to be displayed at the y-tick positions. If None, no
        labels will be shown. If "text", the numeric value of each tick
        will be displayed according to the mapping given in `ytick_format`.
        Otherwise, the listed values will be used.
        Text labels will be formatted using the rotation supplied by
        `yfont_angle` and, the alignment from `yfont_align` and the font size
        given by `yfont_size`. It is recommended that rotated labels are
        right aligned.
    n_xs : array_like, optional
        The position for the count values. If none is supplied, a connt will
        be plotted at each x tick position.
    n_y : float, optional
        The y-position for the count values. If none is supplied, counts will
        be displayed just above the lower axis.
    n_size : int, optional
        Default is 12. The font size for the counts text in points.
    p_x : float, optional
        The x position of the critical value text. Default is to place the
        text in the upper-right hand corner.
    p_y : float, optional
        The y position of the critical value text. By default, the text will
        be placed 90% of the way up the plot.
    p_size : int, optional
        Default is 12. The font size for the critical value text in points.
    title_size: int, optional
        Default is 18. The font size for the title text.
    xtick_size : int, optional
        Default is 12. The font size for the x-tick labels.
    xlabel_size : int, optional
        Default is 15. The fontsize for the x-axis label.
    ytick_size : int, optional
        Default is 12. The font size for the y-tick labels.
    ylabel_size : int, optional
        Default is 15. The fontsize for the y-axis label.
    xfont_angle : float
        The angle in degrees for the x-tick label text.
    xfont_align : {'left', 'right', 'center'}
        The horizonal alignment of the x-tick label text. For rotated text,
        an alignment or 'right' is recommended.
    yfont_angle : float
        The angle in degrees for the x-tick label text.
    yfont_align : {'left', 'right', 'center'}
        The horizonal alignment of the x-tick label text. For rotated text,
        an alignment or 'right' is recommended.
    xtick_format : class, optional
        Default is str. The format for xtick label text if using the
        default tick values as the tick labels.
    ytick_format : class, optional
        Default is str. The format for xtick label text if using the
        default tick values as the tick labels.
    """
    kwds = {'counts': None,
            'n_xs': None,
            'n_y': None,
            'n_size': 11,
            'p_value': None,
            'p_x': None,
            'p_y': None,
            'p_size': 12,
            'show_frame': True,
            'show_xticks': True,
            'show_yticks': True,
            'show_xgrid': False,
            'show_ygrid': False,
            'title': '',
            'title_size': 18,
            'xlim': None,
            'xlabel': '',
            'xticks': None,
            'xticklabels': None,
            'xtick_format': '%s',
            'xfont_align': 'center',
            'xfont_angle': 0,
            'xtick_size': 12,
            'xlabel_size': 15,
            'ylim': None,
            'yticks': None,
            'yticklabels': "text",
            'ytick_format': '%s',
            'ylabel': '',
            'yfont_angle': 0,
            'yfont_align': 'right',
            'ylabel_size': 15,
            'ytick_size': 12}

    for key, val in kwargs.iteritems():
        if key in kwds:
            kwds[key] = val
        else:
            raise ValueError('%s is not an input for axis formatting.' % key)

    # Sets axis ticks
    if kwds['xticks'] is not None:
        ax.set_xticks(kwds['xticks'])

    if kwds['yticks'] is not None:
        ax.set_yticks(kwds['yticks'])

    # Checks the max and min values values
    if kwds['xlim'] is not None and not len(kwds['xlim']) == 2:
        raise ValueError('xlim must specify a xmin and xmax value.')

    if kwds['ylim'] is not None and not len(kwds['ylim']) == 2:
        raise ValueError('xlim must specify a ymin and ymax value.')

    # Sets axis limits
    if kwds['xlim'] is not None:
        ax.set_xlim(kwds['xlim'])
    if kwds['ylim'] is not None:
        ax.set_ylim(kwds['ylim'])

    # Formats the tick labels
    if kwds['xticklabels'] is None:
        xtls = ['']*len(ax.get_xticks())
    elif (isinstance(kwds['xticklabels'], str) and
          kwds['xticklabels'].lower() == 'text'):
        xtls = [kwds['xtick_format'] % t for t in ax.get_xticks()]
    elif not isinstance(kwds['xticklabels'], (list, tuple, np.ndarray)):
        raise TypeError('xticklabels must be None, "text" or an iterable.')
    elif not len(kwds['xticklabels']) == len(ax.get_xticks()):
        raise ValueError('There must be a label for each xtick.')
    else:
        xtls = kwds['xticklabels']

    ax.set_xticklabels(xtls,
                       ha=kwds['xfont_align'],
                       rotation=kwds['xfont_angle'],
                       size=kwds['xtick_size'])

    if kwds['yticklabels'] is None:
        ytls = ['']*len(ax.get_yticks())
    elif (isinstance(kwds['yticklabels'], str) and
          kwds['yticklabels'].lower() == 'text'):
        ytls = [kwds['ytick_format'] % t for t in ax.get_yticks()]
    elif not isinstance(kwds['yticklabels'], (list, tuple, np.ndarray)):
        raise TypeError('yticklabels must be None, "text" or an iterable.')
    elif not len(kwds['yticklabels']) == len(ax.get_yticks()):
        raise ValueError('There must be a label for each ytick.')
    else:
        ytls = kwds['yticklabels']

    ax.set_yticklabels(ytls,
                       ha=kwds['yfont_align'],
                       rotation=kwds['yfont_angle'],
                       size=kwds['ytick_size'])

    # Adds axis labels, if appropriate
    ax.set_xlabel(kwds['xlabel'], size=kwds['xlabel_size'])
    ax.set_ylabel(kwds['ylabel'], size=kwds['ylabel_size'])

    # Adds a grid, if appropriate
    if kwds['show_xgrid']:
        ax.xaxis.grid()
    if kwds['show_ygrid']:
        ax.yaxis.grid()

    # Hides the axis tick marks, if desirable
    if not kwds['show_xticks']:
        for tic in ax.xaxis.get_major_ticks():
            tic.tick1On = False
            tic.tick2On = False
    if not kwds['show_yticks']:
        for tic in ax.yaxis.get_major_ticks():
            tic.tick1On = False
            tic.tick2On = False

    # Removes the axis frame, if desired.
    if not kwds['show_frame']:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        for tic in ax.xaxis.get_major_ticks():
            tic.tick2On = False
        for tic in ax.yaxis.get_major_ticks():
            tic.tick2On = False

    # Adds a title, if appropriate
    ax.set_title(kwds['title'], size=kwds['title_size'])

    # Adds group-size text
    if kwds['counts'] is not None:
        # Gets the positions for the count labels
        ylim = ax.get_ylim()
        if kwds['n_xs'] is None:
            kwds['n_xs'] = ax.get_xticks()
        if kwds['n_y'] is None:
            kwds['n_y'] = ylim[0]+(ylim[1]-ylim[0])*0.03

        # Checks there is a position for each count
        if not len(kwds['counts']) == len(kwds['n_xs']):
            raise ValueError('There must be a position for each count '
                             'being displayed.')
        elif not len(kwds['n_xs']) == len(ax.get_xticks()):
            raise ValueError('There must be a count position for each xtick.')

        # Adds the count labels
        for idx, count in enumerate(kwds['counts']):
            ax.text(kwds['n_xs'][idx],
                    kwds['n_y'],
                    '(%i)' % count,
                    horizontalalignment='center',
                    size=kwds['n_size'])

    # Adds the p-value, if desired
    if kwds['p_value'] is not None:
        p = kwds['p_value']
        # Auto-calculates the position, if necessary
        if kwds['p_x'] is None:
            x = ax.get_xlim()[1] - (ax.get_xlim()[1] - ax.get_xlim()[0])/100
        else:
            x = kwds['p_x']
        if kwds['p_y'] is None:
            y = ax.get_ylim()[1] - (ax.get_ylim()[1] - ax.get_ylim()[0])*0.075
        else:
            y = kwds['p_y']

        # Adds the text
        if p >= 0.005:
            p_str = 'p = %1.2f' % p
        else:
            p_str = 'p = %1.1e' % p
        ax.text(x, y, p_str,
                size=kwds['p_size'],
                horizontalalignment='right')
