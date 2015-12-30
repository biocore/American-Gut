import os

import numpy as np
import pandas as pd
import seaborn as sn

import matplotlib.pyplot as plt

from matplotlib import rcParams

# Formats the axes using seabron so they will be white, and have ticks
# on the bottom of the axes.
sn.set_style('ticks', {'axes.facecolor': 'none'})

# Sets up plotting parameters so that the default setting is use to Helvetica
# in plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Helvetica', 'Arial']
rcParams['text.usetex'] = True


__author__ = "Justine Debelius"
__copyright__ = "Copyright 2015"
__credits__ = ["Justine Debelius"]
__license__ = "BSD"
__version__ = "0.0.1"
__maintainer__ = "Justine Debelius"
__email__ = "Justine.Debelius@colorado.edu"


def plot_alpha(sample, alpha_map, alpha_field, group_field='SIMPLE_BODY_SITE',
               output_dir=None, xlabel=None, debug=False):
    """Generates a distrbution plot for the data

    Parameters
    ----------
    sample : str
        The sample ID to be plotted
    alpha_map_fp : str, dataframe
        The filepath of a comma-seperated file or the pandas dataframe where
        the sample ID are given in the `'#SampleID'` column, a column with
        the name given by `alpha_field` contains alpha diversity values,
        and the `group_field` column specifying the groups which should be
        used to seperate the data for making the distribution plot.
    alpha_field : str
        The name of the column in `alpha_map` which includes the alpha
        diversity values.
    group_field : str
        Default is 'SIMPLE_BODY_SITE'. The name of the column in `alpha_map`
        which provides the grouping for generating distribution plots.
    output_dir : str
        The location where the alpha diversity figures should be saved.
    xlabel : str
        Text describing the quantity on the x-axis.

    Returns
    -------
    If the sample is not included in the the mapping file, a string is returned
    stating this fact.

    If the sample is present, a matplotlib figure with the alpha diversity
    distribution and a line indicating the sample value is returned.

    If debug is passed, the following parameters are returned:
        group : str
            The value of the `group_field` for the sample
        group_alpha : ndarray
            The alpha diversity values assoicated with the group
        sample_alpha : float
            The alpha diversity for the sample
        xlabel : str
            The label used for the x-axis of the plot.

    Raises
    ------
    ValueError
        If the alpha_field is not in alpha_map
    ValueError
        If the group_field is not in alpha_map

    """

    # Checks the alpha_field is in the mapping file
    if alpha_field not in alpha_map.columns:
        raise ValueError('%s is not a valid alpha diversity field name.'
                         % alpha_field)
    # Checks the group_field is in the mapping file
    if group_field not in alpha_map.columns:
        raise ValueError('%s is not a valid field name.' % group_field)
    # Checks the same is in the mapping file
    if sample not in alpha_map.index:
        return ('%s does not have an alpha diversity value for %s.'
                % (alpha_field, sample))

    # Explicitly casts the alpha diversity to a float
    alpha_map[alpha_field] = alpha_map[alpha_field].astype(float)

    # Draws the observations and group
    group = alpha_map.loc[sample, group_field]
    group_alpha = alpha_map.loc[alpha_map[group_field] == group, alpha_field]
    sample_alpha = alpha_map.loc[sample, alpha_field]

    if xlabel is None:
        xlabel = '%sdiversity' % alpha_field.split('1')[0].replace('_', ' ')

    if debug:
        return group, group_alpha, sample_alpha, xlabel

    # Defines the group color. This is currently hardcoded, although the
    # longer term plan is to substitute in function which will define the color
    # based on the relationship between the sample and a yet to be written
    # predicted value.
    group_color = '#1f78b4'
    sample_color = '#525252'

    # Sets up the axis for plotting
    ax = plt.axes()

    # Plots the distribution
    sn.kdeplot(group_alpha,
               ax=ax,
               legend=False,
               color=group_color)
    ylim = ax.get_ylim()

    # Plots the individual line
    ax.plot([sample_alpha, sample_alpha], [-1, 1], color=sample_color)

    # Returns the y-limits to the original value and removes the ticks
    ax.set_ylim(ylim)
    ax.set_yticks([])
    # Removes the spines on the axis and offsets the plot
    sn.despine(offset=5, trim=True, top=True, left=True, right=True)
    # Updates the xticks to match the correct font
    ax.set_xticklabels(map(int, ax.get_xticks()), size=11)
    ax.set_xlabel(xlabel, size=13)

    # Adds text describing the sample
    ax.text(x=ax.get_xticks().max(),
            y=ax.get_ylim()[1]*0.85,
            s='Your Sample:\t%1.1f\nAverage:\t%1.1f'
            % (sample_alpha, group_alpha.mean()),
            ha='right',
            size=11,
            )

    # Sets the figure size
    fig = ax.figure
    fig.set_size_inches((5, 2.5))
    ax.set_position((0.125, 0.375, 0.75, 0.5))

    return fig
