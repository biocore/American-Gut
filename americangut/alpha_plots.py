#!/usr/bin/env python

import numpy as np
import pandas as pd
import seaborn as sb

import matplotlib.pyplot as plt

from matplotlib.collections import PolyCollection

# Sets up all plots generated to use the "whitegrid" style in seaborn.
# See (https://stanford.edu/~mwaskom/software/seaborn/tutorial/aesthetics.html)
# for more information about seaborn styles.
sb.set_style("whitegrid")


def update_map(map_, alpha):
    """Adds alpha diversity to the mapping file

    Parameters
    ----------
    map_ : DataFrame
        The metadata for samples, where the index is the sample ID represented
        as a string.
    alpha : dict
        A dictionary keying the name of the alpha diversity metric a DataFrame
        read from QIIME's `collate_alpha.py`.
    Returns
    -------
    map_ : DataFrame
        The updated metadata including the alpha diversity results.
    """

    # Iterates over the metrics in the dataframe
    for metric, df in alpha.iteritems():
        # Iterates over the depths in the metric
        depths = df.groupby('sequences per sample').groups
        for depth, ids in depths.iteritems():
            # Draws out the alpha diversity values
            data = df.loc[ids, df.columns[3:]].mean(0)
            # Adds the updated data to the mapping file
            data.name = '%s_%i' % (metric.lower(), depth)
            map_ = map_.join(data)

    return map_


def horizontal_to_longtable(sample, map_, categories, metrics, compare=False):
    """Converts the map into a long form table for plotting

    Parameters
    ----------
    sample : str
        The sample barcode of interest.
    map_ : dataframe
        A horizonal (wideform) dataframe containing the information from
        the mapping file. The index should be sample barcodes. If alpha
        diversity is included in `map_`, `metric` should identify these
        columns and `alpha` should be None.
    categories : list
        A list of the metadata categories (columns in `map_`) which should
        be examined.
    metrics : list
        The alpha diversity metrics to be plotted. These should either be
        columns in `map_`.
    compare : bool, optional
        If a category is undefined, should the data be compared to all values.

    Returns
    -------
    long_ : dataframe
        A longform table keyed to a random id, with columns specifying the
        metadata category, metadata group, alpha diversity metric and
        alpha diversity value.

    Raises
    ------
    ValueError
        If `metrics` are not columns in `map_`
    ValueError
        If `categories` are not columns in `map_`
    ValueError
        If `sample` not in the index of `map_`.
    """
    # Checks the sample is present
    if sample not in map_.index:
        raise ValueError('%r is not a valid sample.' % sample)

    # Checks the metrics are present
    map_cats = set(map_.columns)
    met_check = np.array([metric in map_cats for metric in metrics])
    cat_check = np.array([cat in map_cats for cat in categories])
    if not met_check.all():
        raise ValueError('Values for %s cannot be found.'
                         % (', ').join(np.array(metrics)[~met_check]))
    if not cat_check.all():
        raise ValueError('Values for %s cannot be found.'
                         % (', ').join(np.array(categories)[~cat_check]))

    # Sets up the longform table
    long_ = pd.DataFrame(columns=['category', 'group', 'alpha', 'metric'])

    # Creates the metric work
    for category in categories:
        for metric in metrics:
            df = map_.loc[map_[category] == map_.loc[sample, category],
                          [category, metric]]
            df['group'] = map_.loc[sample, category]
            df['metric'] = metric.split('_')[0]
            df['category'] = category.replace('_', ' ').title()
            df.rename(columns={metric: 'alpha'}, inplace=True)
            long_ = pd.concat((long_, df))

    long_ = long_[['category', 'group', 'alpha', 'metric']]

    return long_


def modify_alpha_shannon_pd(long_, mod):
    """Modifies PD 1k to plot next to shannon diversity

    To be able to plot alpha diversity side-by-side, the PD diversity must
    be rescaled. This performs the rescaling for 1k rarefaction. The code
    uses 3.5 as a magic number; this is based on rounds 1-21 of the American
    Gut data. The code is intended as a one-off for making alpha diversity
    violin plots until a more permenant solution is introduced.

    Parameters
    ----------
    long_ : dataframe
        The longform metadata with metrics `shannon` and `PD`.

    Returns
    -------
    long_ : dataframe
        The longform medadata frame with the adjusted metdata.

    Raises
    ------
    ValueError
        When `shannon` and `PD` are not groups in the `metric` column in long_.

    """

    if 'metric' not in long_.columns:
        raise ValueError('metric must be a column in long_')
    metrics = set(long_['metric'].values)
    if 'shannon' not in metrics and 'pd' not in metrics:
        raise ValueError('shannon and pd must be metrics in long_')

    shannon_index = long_.metric == 'shannon'
    pd_index = long_.metric == 'pd'

    long_.loc[shannon_index, 'alpha_mod'] = long_.loc[shannon_index, 'alpha']
    long_.loc[pd_index, 'alpha_mod'] = (long_.loc[pd_index, 'alpha'])/mod
    long_.loc[shannon_index, 'metric'] = 'Shannon Diversity'
    long_.loc[pd_index, 'metric'] = 'PD Whole Tree Diversity'

    return long_


def ag_summary_violin(sample, map_, alpha, filepath=None, ax=None,
                      categories=None, metrics=None, debug=False):
    """Creates a single sample Violin Plot for American Gut participants"""

    if sample not in map_.index:
        raise ValueError('%s is not in map' % sample)

    # Builds the alpha map
    alpha_map = update_map(map_, alpha)

    # Gets the alpha diversity map for the bodysite
    body_habitat = alpha_map.loc[sample, 'BODY_HABITAT']
    alpha_map = update_map(map_.loc[map_['BODY_HABITAT'] == body_habitat],
                           alpha)

    if body_habitat == 'UBERON:feces':
        categories = ['BODY_HABITAT', 'AGE_CAT', 'ALCOHOL_FREQUENCY',
                      'ANTIBIOTIC_HISTORY', 'BMI_CAT', 'FLOSSING_FREQUENCY']
        labels = ['All Fecal Samples', 'Simillar Age', 'Alcohol Use',
                  'Last Antibiotic Use', 'Simillar BMI',
                  'Similar Flossing \nFrequency']
        mod = 3
        ylim = [0, 12]
        shannon_yticks = np.arange(0, 13, 2)
        pd_yticks = np.arange(0, 39, 6)

    elif body_habitat == 'UBERON:oral cavity':
        categories = ['BODY_HABITAT', 'AGE_CAT', 'ALCOHOL_FREQUENCY',
                      'ANTIBIOTIC_HISTORY', 'BMI_CAT', 'FLOSSING_FREQUENCY']
        labels = ['All Oral Samples', 'Simillar Age', 'Alcohol Use',
                  'Last Antibiotic Use', 'Simillar BMI',
                  'Similar Flossing \nFrequency']
        mod = 2.5
        ylim = [0, 10]
        shannon_yticks = np.arange(0, 11, 2)
        pd_yticks = np.arange(0, 26, 5)

    elif body_habitat == 'UBERON:skin':
        categories = ['BODY_HABITAT', 'AGE_CAT', 'ALCOHOL_FREQUENCY',
                      'ANTIBIOTIC_HISTORY', 'BMI_CAT', 'DOMINANT_HAND']
        labels = ['All Skin Samples', 'Simillar Age', 'Alcohol Use',
                  'Last Antibiotic Use', 'Simillar BMI',
                  'Same Dominant hand']
        mod = 4
        ylim = [0, 12.5]
        shannon_yticks = np.arange(0, 13, 2.5)
        pd_yticks = np.arange(0, 51, 10)

    else:
        raise ValueError('%s is not a supported body habitat.' % body_habitat)

    # Supplies default values
    if metrics is None:
        metrics = ['shannon_1000', 'pd_whole_tree_1000']

    # Gets the long form table
    long_ = horizontal_to_longtable(sample, alpha_map, categories, metrics)
    modify_alpha_shannon_pd(long_, mod)

    # Identifes the alpha diversity and modified diversity value for the sample
    (pdline, shannon) = long_.loc[sample].groupby('metric').mean()['alpha_mod']
    pdvalue = long_.loc[sample].groupby('metric').mean(
        ).loc['PD Whole Tree Diversity', 'alpha']

    # If the code is in debugging mode, the values are returned
    if debug:
        return (categories, labels, metrics, mod, ylim, shannon_yticks,
                pd_yticks, pdline, shannon, pdvalue)

    # Sets up axis limits.
    # The axis limits here are selected as magic numbers. Yes, I, the author
    # acknowelge magic numbers are generally a bad idea. However, in this case,
    # the values were selected based on 1k rarefaction AG 1-21. Given that
    # this is one-off plotting code for American Gut violin plots, and its
    # typically difficult to test plotting code, a set value and a magic number
    # rather than a parameter seems the lesser of two evils.

    # Creates an axis instance for plotting
    if ax is None:
        ax = plt.axes()
    ax.set_ylim(ylim)

    # Creates the initial violin plot
    sb.violinplot(x=long_['category'],
                  y=long_['alpha_mod'],
                  hue=long_['metric'],
                  split=True,
                  inner='quartile',
                  ax=ax
                  )
    xlim = ax.get_xlim()

    # Gets the face colors and uses alpha to make the fill lighter
    two_face_color = []
    for idx, child in enumerate(ax.get_children()):
        if isinstance(child, PolyCollection):
            if len(two_face_color) == 0:
                two_face_color.append(child.get_facecolor()[0])
                two_face_color.append(
                    ax.get_children()[idx+1].get_facecolor()[0]
                    )
            child.set_alpha(0.4)

    leg = ax.get_legend()

    # Idenfies the metrics
    [m0, m1] = [t.get_text() for t in leg.get_texts()]
    for t in leg.get_texts():
        t.set_size(10)
    leg.set_title('')

    # Formats the inner axis
    ax.set_ylim(ylim)
    ax.set_yticks(shannon_yticks)
    ax.set_yticklabels(shannon_yticks,
                       color=two_face_color[0],
                       size=13,
                       weight='bold')
    ax.set_ylabel(m0,
                  size=15,
                  color=two_face_color[0])

    # Updates the xaxis values to contain the label text and be vertical
    ax.set_xticklabels(labels, size=11, ha='right', rotation=90)
    ax.set_xlabel('')

    # Formats the second axis
    ax2 = ax.twinx()
    ax2.set_ylim(ylim)
    ax2.set_yticks(shannon_yticks)
    ax2.set_yticklabels(pd_yticks, color=two_face_color[1], size=13,
                        weight='bold')
    ax2.set_ylabel(m1.replace('_1000', '').replace('_', ' '),
                   size=15,
                   color=two_face_color[1])
    ax2.grid(False)
    ax2.set_frame_on(False)

    # Plots the alpha diversity value on the second axis
    ax2.plot([-5, 15], [shannon]*2,
             color=two_face_color[0],
             linestyle='-'
             )
    ax2.plot([-5, 15], [pdline]*2,
             color=two_face_color[1],
             linestyle='--')
    ax2.set_xlim(xlim)

    # Adds a legend with the alpha diversity value
    ax2.legend(loc=2,
               labels=['Shannon Diveristy (%1.2f)' % shannon,
                       'PD Whole Tree Diversity (%1.1f)' % pdvalue],
               fontsize=10)

    fig = ax.figure
    if filepath is not None:
        fig.savefig(filepath)
    else:
        return fig
