import numpy as np

from future.utils import viewitems

from statsmodels.stats.power import FTestAnovaPower
from skbio.stats.power import confidence_bound

from matplotlib import rcParams

ft = FTestAnovaPower()

__author__ = "Justine Debelius"
__copyright__ = "Copyright 2015, The American Gut Project"
__credits__ = ["Justine Debelius"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Justine Debelius"
__email__ = "Justine.Debelius@colorado.edu"

# Gives plots pretty helvetica text!
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Helvetica', 'Arial']
rcParams['text.usetex'] = True


def collate_effect_size(counts, powers, alpha):
    """Calculates the effects for power values

    Parameters
    ----------
    counts : array
        the number of samples used to calculate the power
    powers : {list, ndarray}
        list of arrays of power values. If there are multiple power arrays,
        each power array should have the same dimensions. If samples are
        missing, these should be denoted by a nan.
    alpha : float
        the critical value used to calculate the power.

    Returns
    -------
    effect_means : array
        A 1-dimensional array of the mean effect sizes for the counts array.
    effect_bounds : array
        A 1-dimensional array with the confidence bound for the effect size
        array. The lower bound is `effect_mean - effect_bounds` and the upper
        bound is `effect_means - effect_bounds1`.

    Raises
    ------
    TypeError
        if counts is not a one-dimensional array
    ValueError
        if the arrays in powers have different shapes
    ValueError
        if the length of the power arrays and the length of the count arrays
        are different
    """

    # Checks the power and counts iteratibility
    if isinstance(powers, np.ndarray):
        powers = [powers]

    num_powers = len(powers)
    if isinstance(counts, np.ndarray):
        counts = [counts]*num_powers

    # Checks there is a count for each power
    if not len(counts) == len(powers):
        raise ValueError('There must be a counts array for each power array.')

    # Checks the shape array
    for idx in xrange(num_powers):
        count_shape = counts[idx].shape
        power_shape = powers[idx].shape
        # Checks the count array is 1d
        if not len(count_shape) == 1:
            raise TypeError('Each count array must be a 1d array.')

        if len(power_shape) == 1:
            if not count_shape[0] == power_shape[0]:
                raise ValueError('There must be a sample count for each '
                                 'power.')
        elif not count_shape[0] == power_shape[1]:
            raise ValueError('There must be a sample count for each power.')

    # Prealocates the output arrays
    effect_means = np.zeros((num_powers))
    effect_bounds = np.zeros((num_powers))

    # Iterates through the powers and calculates the effect sizes
    for idp, pwr in enumerate(powers):
        count = counts[idp]
        pwr_shape = pwr.shape
        # Calculates the effect size for the power array
        eff = np.zeros(pwr_shape)*np.nan
        for id2, cnt in enumerate(count):
            if len(pwr_shape) == 1:
                try:
                    eff[id2] = ft.solve_power(effect_size=None,
                                              nobs=cnt,
                                              alpha=alpha,
                                              power=pwr[id2])
                except:
                    pass
            else:
                for id1 in xrange(pwr_shape[0]):
                    try:
                        eff[id1, id2] = ft.solve_power(effect_size=None,
                                                       nobs=cnt,
                                                       alpha=alpha,
                                                       power=pwr[id1, id2])
                    except:
                        pass
        # Caluclates the mean and bound
        if np.isnan(eff).all():
            effect_means[idp] = np.nan
            effect_bounds[idp] = np.nan
        else:
            effect_means[idp] = np.nanmean(eff, None)
            effect_bounds[idp] = confidence_bound(eff, alpha, None)

    return effect_means, effect_bounds


def summarize_effect(order, fecal_cats, a_eff_means, b_eff_means, b_eff_bounds,
                     a_eff_bounds):
    """Creates a pretty html formatted table of the results

    The output table can be displayed in an IPython notebook as an HTML object,
    rather than building a pandas table. Hopefully, this will allow for a
    beautiful display. The numeric cells are split into three columns to allow
    for numerical alignment.

    The effect sizes are rounded to 5 samples, and the standard devation are
    rounded to 1 sample.

    Parameters
    ----------
    order : array
        The order in which categories should be displayed.
    fecal_cats : list or iterables
        A list giving the fecal category and order for each category to be
        displayed.
    a_eff_means, b_eff_means : list of arrays
        An array with the mean effects. The output of `collate_effect_size`.
    a_eff_bounds, b_eff_bounds : list of arrays
         An array with the mean bound. The confidence interval for the mean
         is given by [eff_mean - eff_bound, eff_mean + eff_bound].
         This is an output of `collate_effect_size`.

    Returns
    -------
    table : str
        An html-formmated table which can be rendered in IPython.


    """
    # Sets up the html for the header
    table = ['<table style="border-style:hidden;',
             '              border-collapse:collapse',
             '              line-height:120%',
             '              ">',
             '\t<tr>',
             '\t\t<th style="text-align:center;',
             '\t\t           background-color:black;',
             '\t\t           color:white',
             '\t\t           ">',
             '\t\t\tCategory'
             '\t\t</th>',
             '\t\t<th style="text-align:center;',
             '\t\t           background-color:black;',
             '\t\t           color:white";',
             '\t\t    colspan=3>',
             '\t\t\tAlpha',
             '\t\t</th>',
             '\t\t<td style="border-hidden;',
             '\t\t           background-color:black;',
             '\t\t           padding:20px">',
             '\t\t<th style="text-align:center;',
             '\t\t           background-color:black;',
             '\t\t         color:white";',
             '\t\t  colspan=3>',
             '\t\t\tBeta',
             '\t\t</th>',
             '\t</tr>']
    # Loops through each row
    for idx in order:
        # Gets the effect sizes
        cat = fecal_cats[idx][0]
        a_eff_m = a_eff_means[idx]
        a_eff_b = a_eff_bounds[idx]
        b_eff_m = b_eff_means[idx]
        b_eff_b = b_eff_bounds[idx]
        # Gets the fitted values
        a_eff_fit = np.ceil(ft.solve_power(a_eff_m, nobs=None, power=0.8,
                                           alpha=0.05)/5.)*5
        a_eff_rnd = np.ceil(np.abs(ft.solve_power(a_eff_m, nobs=None,
                                                  power=0.8, alpha=0.05) -
                                   ft.solve_power(a_eff_m - a_eff_b,
                                                  nobs=None, power=0.8,
                                                  alpha=0.05))/5.)*5
        b_eff_fit = np.ceil(ft.solve_power(b_eff_m, nobs=None, power=0.8,
                                           alpha=0.05)/5.)*5
        b_eff_rnd = np.ceil(np.abs(ft.solve_power(b_eff_m, nobs=None,
                                                  power=0.8, alpha=0.05) -
                                   ft.solve_power(b_eff_m - b_eff_b, nobs=None,
                                                  power=0.8, alpha=0.05))/5.)*5
        # Fills in the html text
        row = ['\t<tr>',
               '\t\t<td style="border-top:hidden;',
               '\t\t           border-bottom:hidden;',
               '\t\t           border-left: hidden;',
               '\t\t           border-bottom: hidden;',
               '\t\t           padding:10px;',
               '\t\t           text-align:left',
               '\t\t          ">',
               '\t\t\t%s' % cat.replace('_', ' ').title(),
               '\t\t</td>',
               '\t\t<td style="border-top:hidden;',
               '\t\t           border-bottom:hidden;',
               '\t\t           border-left: hidden;',
               '\t\t           border-bottom: hidden;',
               '\t\t           text-align:right',
               '\t\t          ">',
               '\t\t\t%i' % a_eff_fit,
               '\t\t</td>',
               '\t\t<td style="border-top:hidden;',
               '\t\t           border-bottom:hidden;',
               '\t\t           border-left: hidden;',
               '\t\t           border-bottom: hidden;',
               '\t\t           text-align:center',
               '\t\t          ">',
               '\t\t\t&plusmn;',
               '\t\t</td>',
               '\t\t<td style="border-top:hidden;',
               '\t\t           border-bottom:hidden;',
               '\t\t           border-left: hidden;',
               '\t\t           border-bottom: hidden;',
               '\t\t           text-align:right',
               '\t\t          ">',
               '\t\t\t%i' % a_eff_rnd,
               '\t\t</td>',
               '\t\t<td style="border-top:hidden;',
               '\t\t           border-bottom:hidden;',
               '\t\t           border-left: hidden;',
               '\t\t           border-bottom: hidden;',
               '\t\t           padding:20px">',
               '\t\t</td>',
               '\t\t<td style="border-top:hidden;',
               '\t\t           border-bottom:hidden;',
               '\t\t           border-left: hidden;',
               '\t\t           border-bottom: hidden;',
               '\t\t           text-align:right',
               '\t\t          ">',
               '\t\t\t%i' % b_eff_fit,
               '\t\t</td>',
               '\t\t<td style="border-top:hidden;',
               '\t\t           border-bottom:hidden;',
               '\t\t           border-left: hidden;',
               '\t\t           border-bottom: hidden;',
               '\t\t           text-align:center',
               '\t\t          ">',
               '\t\t\t&plusmn;',
               '\t\t</td>',
               '\t\t<td style="border-top:hidden;',
               '\t\t           border-bottom:hidden;',
               '\t\t           border-left: hidden;',
               '\t\t           border-bottom: hidden;',
               '\t\t           text-align:right',
               '\t\t          ">',
               '\t\t\t%i' % b_eff_rnd,
               '\t\t</td>',
               ]
        table.append('\n'.join(row))
    table.append('</table>')

    return '\n'.join(table)


def plot_effects(effect_means, effect_bounds, labels, sample_counts, **kwargs):

    """Makes a power curve plot

    Parameters
    ----------
    effect_means: 1d array
        the mean effect sizes to plots.
    effect_bounds : {None, 1d array}
        the range used for the confidence interval. If there is no effect to
        show, this should be None.
    labels : 1d array
        a list of formatted strings describing the effects, to be used in the
        legend.
    sample_counts : 1d array
        the counts where power should be calculated.
    alpha : int, optional
        Default is 0.05. The critical value for the power curves.
    colormap : {None, array}, optional
        Default is None. A colormap to use for the lines. Each color
        designation must appear in a new row. If no colormap is supplied, the
        defualt colormap will be used.
    grid : bool, optional
        Default is True. Show grid.
    show_bound : bool
        Default is True. Shows the confidence bounds on the effect size. If
        `effect_bounds` is None, no bounds will be shown.

    Returns
    -------
    fig : figure
        a figure with the power curves plotted.

    Other parameters
    ----------------
    leg_offset : tuple
        Changes the legend position.
    tick_size : usigned int
        sets the font size for tick labels
    label_size : unsigned int
        sets the font size for the axis labels
    title_size : unsigned int
        sets the font size for the title
    legend_size : unsigned int
        sets the font size for enteries in the legend

    """

    # Sets the keyword properties
    kwds = {'alpha': 0.05,
            'colormap': None,
            'grid': True,
            'title': '',
            'show_bound': True,
            'leg_offset': None,
            'tick_size': 12,
            'label_size': 15,
            'title_size': 18,
            'legend_size': 11}
    for key, value in viewitems(kwargs):
        if key in kwds:
            kwds[key] = value
        else:
            raise ValueError('%s is not a property of plot_effects.' % key)

    # Checks the effect, bound, and mean argument is sane
    mean_shape = effect_means.shape
    if effect_bounds is None:
        kwds['show_bound'] = False
        effect_bounds = np.zeros(mean_shape)
    bound_shape = effect_bounds.shape
    label_shape = labels.shape

    if not len(mean_shape) == 1:
        raise ValueError('Effect Mean must be a 1d numpy array')
    elif mean_shape != bound_shape or mean_shape != label_shape:
        raise ValueError('There must be a label and bound for each effect.')

    # Plots the the lower bound data
    fig = ft.plot_power(dep_var='nobs',
                        nobs=sample_counts,
                        effect_size=effect_means - effect_bounds,
                        alpha=kwds['alpha'])
    # Gets the axis of the first plot and its position
    lax = fig.axes[0]
    # Makes the lower bound lines dashed and thin, and changes the color if
    # desired
    for idx, l in enumerate(lax.get_lines()):
        l.set_linestyle(':')
        l.set_linewidth(1.5)
        if kwds['colormap'] is not None and len(kwds['colormap'].shape) == 1:
            l.set_color(kwds['colormap'][idx])
        elif kwds['colormap'] is not None:
            l.set_color(kwds['colormap'][idx, :])
    # Hides the x ticks and labels
    lax.set_title('')
    lax.set_xticklabels('')
    lax.set_yticklabels('')
    lax.set_xlabel('')
    # Hides the legend
    lax.get_legend().set_visible(False)

    # Plots the upper bound data
    uax = fig.add_axes(lax.get_position())
    fig = ft.plot_power('nobs', sample_counts, effect_means + effect_bounds,
                        alpha=kwds['alpha'], ax=uax)
    # Makes the lower bound axes visable, if desired
    if kwds['show_bound']:
        uax.set_axis_bgcolor('none')
    # Makes the lower bound lines dashed and thin, and changes the color if
    # desired
    for idx, l in enumerate(uax.get_lines()):
        l.set_linestyle(':')
        l.set_linewidth(1.5)
        if kwds['colormap'] is not None and len(kwds['colormap'].shape) == 1:
            l.set_color(kwds['colormap'][idx])
        elif kwds['colormap'] is not None:
            l.set_color(kwds['colormap'][idx, :])
    # Hides the x ticks and labels
    uax.set_title('')
    uax.set_xticklabels('')
    uax.set_yticklabels('')
    uax.set_xlabel('')
    # Hides the legend
    uax.get_legend().set_visible(False)

    # Plots the mean data
    axm = fig.add_axes(lax.get_position())
    fig = ft.plot_power('nobs', sample_counts, effect_means, ax=axm,
                        alpha=kwds['alpha'])

    # Shows the confidence bounds, if desired
    if kwds['show_bound']:
        axm.set_axis_bgcolor('none')

    # Recolors the lines, if desired
    if kwds['colormap'] is not None and len(kwds['colormap'].shape) == 1:
        for idx, l in enumerate(axm.get_lines()):
            l.set_color(kwds['colormap'][idx])
    elif kwds['colormap'] is not None:
        for idx, l in enumerate(axm.get_lines()):
            l.set_color(kwds['colormap'][idx, :])

    # Sets up the labels
    axm.set_xticklabels(map(int, axm.get_xticks()), size=kwds['tick_size'])
    axm.set_yticklabels(axm.get_yticks(), size=kwds['tick_size'])
    axm.set_xlabel('Number of Observations', size=kwds['label_size'])
    axm.set_ylabel('Power of the Test', size=kwds['label_size'])
    axm.set_title(kwds['title'], size=kwds['title_size'])

    # Adds the grid, if desired
    if kwds['grid']:
        axm.grid()

    leg = axm.get_legend()
    # Sets the legend position
    if kwds['leg_offset'] is not None:
        leg.set_bbox_to_anchor(kwds['leg_offset'])
    # Sets up the legend text
    for idx, txt in enumerate(leg.get_texts()):
        txt.set_text(labels[idx])
        txt.set_size(kwds['legend_size'])

    # Returns the figure
    return fig


def trace_bounds(power, counts):
    """Gets the mean, upper and lower confidence interval for counts

    Parameters
    ----------
    power : ndarray
        An array with multiple power returns where each row corresponds to a
        permutation run and each column corresponds to a number of samples
        drawn.
    counts : ndarray
        The number of samples drawn per group for each power calculation.

    Returns
    -------
        pwr_mean, pwr_lower, pwr_upper : ndarray
            The mean, upper and lower bounds for the power array.

    """

    pwr_mean = np.nanmean(power, 0)
    pwr_bound = confidence_bound(power, axis=0)

    # Prevents the confidence interval of being less than 0 or more than 1
    pwr_lower = np.vstack((pwr_mean - pwr_bound,
                           np.zeros(pwr_mean.shape))).max(0)
    pwr_upper = np.vstack((pwr_mean + pwr_bound,
                           np.ones(pwr_mean.shape))).min(0)

    return pwr_mean, pwr_lower, pwr_upper


def add_average_trace(fig, power, counts, labels, **kwargs):
    """Adds a trace of averaged points

    This can be used with power curves which are too steep to use to estimate
    the effect size.

    Parameters
    ----------
    fig : figure
        A matplitlib figure instance generated by `plot_effects` which has
        three overlayed axes.
    power : array
        An array of the power estimation results from
        `skbio.stats.power.subsample_power` or
        `skbio.stats.power.subsample_paired_power`.
    counts : array
        A 1 dimensional numpy array which corresponds to the number of counts
        used to estimate the power value at the corresponding position in
        `power`.
    labels : array
        The labels for all the curves show in `fig`, including the curve being
        added in the appropriate order.
    figure_size : tuple
        The figure size in inches.
    axis_size : tuple
        The axis size in inches, specifying the (left, bottom, width, height)
        for the final axes.
    legend_pad : tuple
        The axis size to space a legend power. The tuple species the
        (left, bottom, width, height) for the dummy axes.
    legend_font_size : int
        The size in pixels of the legend text to be displayed.
    legend_position : tuple
        The position of the legend relative to the axes.
    color : list or array
        The color of the line being added.

    """

    # Sets up the keyword arguements
    kwds = {'figure_size': (7., 3.5),
            'axis_size': (0.6, 0.6, 4., 2.5),
            'legend_pad': (5., 0.6, 2., 2.5),
            'legend_font_size': 11,
            'legend_position': (1.05, 0.95),
            'color': [0.5, 0.5, 0.5]}

    # Updates the keyword arguments
    for kwd, val in viewitems(kwargs):
        if kwd not in kwds:
            raise ValueError('%s is not a supported add_bodysite key word.')
        else:
            kwds[kwd] = val

    # Sets up the figure size and dimensions
    fig_dims = (kwds['axis_size'][0] / kwds['figure_size'][0],
                kwds['axis_size'][1] / kwds['figure_size'][1],
                kwds['axis_size'][2] / kwds['figure_size'][0],
                kwds['axis_size'][3] / kwds['figure_size'][1])

    pad_dims = (kwds['legend_pad'][0] / kwds['figure_size'][0],
                kwds['legend_pad'][1] / kwds['figure_size'][1],
                kwds['legend_pad'][2] / kwds['figure_size'][0],
                kwds['legend_pad'][3] / kwds['figure_size'][1])

    # Calculates the average and bound for the bodysite
    site_mean, site_lower, site_upper = trace_bounds(power, counts)

    # Adds the body site line to the appropriate axes
    axes = fig.axes
    axes[0].plot(counts, site_lower, ':', color=kwds['color'])
    axes[1].plot(counts, site_upper, ':', color=kwds['color'])
    axes[2].plot(counts, site_mean, linewidth=2, color=kwds['color'],
                 label=labels)

    # Adjusts the figure size to allow for the new legend placement
    fig.set_size_inches(kwds['figure_size'])

    # Adjusts the axis size for the sub axes
    axes[0].set_position(fig_dims)
    axes[1].set_position(fig_dims)
    axes[2].set_position(fig_dims)

    # Adds a new axis to pad for the legend display
    pad_ax = fig.add_axes(pad_dims)
    pad_ax.set_visible(False)

    # Hides the current axis legend so a new legend showing body site can be
    # added.
    axes[2].get_legend().set_visible(False)

    # Adjusts the position of the filte vector so the bodysite line will appear
    # first in the legend
    handles = list(axes[2].get_legend_handles_labels())
    handles[0].insert(0, handles[0].pop(-1))

    # Adds a new legend to the figure
    leg2 = fig.legend(handles=handles[0],
                      labels=labels,
                      fontsize=kwds['legend_font_size'])

    # Adjusts the legend position
    leg2.set_bbox_to_anchor(kwds['legend_position'])
