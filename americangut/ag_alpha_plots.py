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


def plot_alpha_distribution(sample, alpha_map, alpha_name='alpha',
                            ax=None, xlabel=None, group_field='BODY_HABITAT',
                            debug=False):
    """Generates a distrbution plot for the data"""

    # Explicitly casts the alpha diversity to a float
    alpha_map[alpha_name] = alpha_map[alpha_name].astype(float)

    # Draws the observations and group
    group = alpha_map.loc[sample, group_field]
    group_alpha = alpha_map.loc[alpha_map[group_field] == group, alpha_name]
    sample_alpha = alpha_map.loc[sample, alpha_name]

    # Sets up the defaults
    if ax is None:
        ax = plt.axes()

    if xlabel is None:
        xlabel = alpha_name

    # Returns values if debug mode is true (basically to allow for testing)
    if debug:
        return group, group_alpha, sample_alpha, ax, xlabel

    # Defines the group color. This is currently hardcoded, although the
    # longer term plan is to substitute in function which will define the color
    # based on the relationship between the sample and a yet to be written
    # predicted value.
    group_color = '#1f78b4'
    sample_color = '#525252'

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
            % (sample_alpha, sample_alpha.mean()),
            ha='right',
            size=11,
            )

    # Sets the figure size
    fig = ax.figure
    fig.set_size_inches((5, 2.5))
    ax.set_position((0.125, 0.375, 0.75, 0.5))

    return fig