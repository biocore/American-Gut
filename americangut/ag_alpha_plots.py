import os

import click

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


@click.group()
def alpha_plots():
    pass


@alpha_plots.command()
@click.option('--sample', required=True, type=str, help='sample to be plotted')
@click.option('--alpha_map',
              required=True,
              help="The mapping file with alpha diversity. The sample id "
              "column should be named `'#SampleID'`, and the alpha diversity"
              " metric should be specified by `alpha_field`.",
              type=click.Path(resolve_path=True, readable=True, exists=True,))
@click.option('--output_dir', required=True,
              help='The location where the alpha diversity figures should be'
              ' saved.',
              type=click.Path(exists=True))
@click.option('--alpha_field', required=True,
              type=click.Choice(['PD_whole_tree_1k', 'PD_whole_tree_10k',
                                 'chao1_1k', 'chao1_10k',
                                 'observed_otus_1k', 'observed_otus_10k',
                                 'shannon_1k', 'shannon_10k', 'alpha']),
              help='The name of the column in `alpha_map` which gives the '
              'alpha diversity values')
@click.option('--group_field', type=str, default='SIMPLE_BODY_SITE',
              help='The name of the column in `alpha_map` which provides '
              'the grouping for generating distribution plots. For example,'
              ' if a mapping file encompasses multiple bodysites, it may be '
              'useful to split by bodysite.')
@click.option('--xlabel', type=str, help='The label for the x axis.')
def plot_alpha(sample, alpha_map, alpha_field='alpha',
               group_field='SIMPLE_BODY_SITE', output_dir=None, xlabel=None):
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
        diversity values. Options include 'PD_whole_tree_1k',
        'PD_whole_tree_10k', 'chao1_1k', 'chao1_10k','observed_otus_1k',
        'observed_otus_10k', 'shannon_1k', 'shannon_10k', or 'alpha'
    group_field : str
        Default is 'SIMPLE_BODY_SITE'. The name of the column in `alpha_map`
        which provides the grouping for generating distribution plots.
    output_dir : str
        The location where the alpha diversity figures should be saved.
    xlabel : str
        Text describing the quantity on the x-axis.

    Returns
    -------
    A plot with the alpha diversity distribution and the sample value
    highlighted.

    Raises
    ------
    TypeError
        If alpha_map is not a filepath string or pandas dataframe.
    ValueError
        If the sample is not in the mapping file.
    ValueError
        If the alpha_field is not in alpha_map
    ValueError
        If the group_field is not in alpha_map

    """
    

    alpha_map = pd.read_csv(alpha_map, sep='\t', dtype=str)
    alpha_map.set_index('#SampleID', inplace=True)

    click.echo(alpha_map.columns)

    # Checks the same is in the mapping file
    if sample not in alpha_map.index:
        raise ValueError('%s is not a valid sample.' % sample)
    # Checks the alpha_field is in the mapping file
    if alpha_field not in alpha_map.columns:
        raise ValueError('%s is not a valid alpha diversity field name.'
                         % alpha_field)
    # Checks the group_field is in the mapping file
    if group_field not in alpha_map.columns:
        raise ValueError('%s is not a valid field name.' % group_field)

    # Explicitly casts the alpha diversity to a float
    alpha_map[alpha_field] = alpha_map[alpha_field].astype(float)

    # Draws the observations and group
    group = alpha_map.loc[sample, group_field]
    group_alpha = alpha_map.loc[alpha_map[group_field] == group, alpha_field]
    sample_alpha = alpha_map.loc[sample, alpha_field]

    if xlabel is None:
        xlabel = '%sdiversity' % alpha_field.split('1')[0].replace('_', ' ')

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

    if output_dir is not None:
        out_fp = os.path.join(output_dir, '%s_%s.pdf' % (alpha_field, sample))
        fig.savefig(out_fp, dpi=300)
    else:
        return fig

if __name__ == '__main__':
    alpha_plots()
