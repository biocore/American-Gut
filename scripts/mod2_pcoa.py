#!/usr/bin/env python

import os
import click
from matplotlib import use
use('Agg')  # noqa
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from skbio import read, DistanceMatrix
from skbio.stats import isubsample
from skbio.stats.ordination import OrdinationResults
from collections import defaultdict
from collections import OrderedDict


ALPHA = 1.0
LINE_WIDTH = 0.3
LINE_WIDTH_WHITE = 2.0
LINE_WIDTH_BLACK = 1.0


@click.group()
def mod2_pcoa():
    pass


@mod2_pcoa.command()
@click.option('--coords', required=True, type=click.Path(
              resolve_path=True, readable=True, exists=True),
              help='Coordinates file')
@click.option('--mapping_file', required=True, type=click.Path(
              resolve_path=True, readable=True, exists=True),
              help='Mapping file')
@click.option('--output', required=True, type=click.Path(exists=True,
              writable=True, resolve_path=True), help='Output directory')
@click.option('--prefix', required=True, type=str, help='Output file prefix')
@click.option('--samples', required=False, type=str,
              help='Comma separated list of samples to print')
def body_site(coords, mapping_file, output, prefix, samples):
    """Generates as many figures as samples in the coordinates file"""
    o = read(coords, into=OrdinationResults)

    # coordinates
    c_df = pd.DataFrame(o.site, o.site_ids)

    # mapping file
    mf = pd.read_csv(mapping_file, '\t', converters=defaultdict(str),
                     index_col='#SampleID')
    mf = mf.loc[o.site_ids]

    if samples is None:
        samples = mf.index
    else:
        samples = set(samples.split(',')).intersection(set(o.site_ids))
        samples = mf.loc[samples].index

    color_hmp_fecal = sns.color_palette('Paired', 12)[10]  # light brown
    color_agp_fecal = sns.color_palette('Paired', 12)[11]  # dark brown
    color_hmp_oral = sns.color_palette('Paired', 12)[0]    # light blue
    color_agp_oral = sns.color_palette('Paired', 12)[1]    # dark blue
    color_hmp_skin = sns.color_palette('Paired', 12)[2]    # light green
    color_agp_skin = sns.color_palette('Paired', 12)[3]    # dark green

    grp_colors = {'AGP-FECAL': color_agp_fecal,
                  'AGP-ORAL':  color_agp_oral,
                  'AGP-SKIN':  color_agp_skin,
                  'HMP-FECAL': color_hmp_fecal,
                  'GG-FECAL':  color_hmp_fecal,
                  'PGP-FECAL': color_hmp_fecal,
                  'HMP-ORAL':  color_hmp_oral,
                  'PGP-ORAL':  color_hmp_oral,
                  'HMP-SKIN':  color_hmp_skin,
                  'PGP-SKIN':  color_hmp_skin}

    for sample in samples:

        # plot categories as 50 slices with random zorder
        for grp, color in grp_colors.iteritems():
            sub_coords = c_df[mf.TITLE_BODY_SITE == grp].values
            for i in np.array_split(sub_coords, 50):
                plt.scatter(i[:, 0], i[:, 1], color=color,
                            edgecolor=np.asarray(color)*0.6, lw=LINE_WIDTH,
                            alpha=ALPHA, zorder=np.random.rand())

        # plot participant's dot
        plt.scatter(c_df.loc[sample][0], c_df.loc[sample][1],
                    color=grp_colors[mf.loc[sample]['TITLE_BODY_SITE']],
                    s=270, edgecolor='w', zorder=1, lw=LINE_WIDTH_WHITE)
        plt.scatter(c_df.loc[sample][0], c_df.loc[sample][1],
                    color=grp_colors[mf.loc[sample]['TITLE_BODY_SITE']],
                    s=250, edgecolor=np.asarray(
                    grp_colors[mf.loc[sample]['TITLE_BODY_SITE']])*0.6,
                    zorder=2, lw=LINE_WIDTH_BLACK)

        plt.axis('off')
        my_dpi = 72
        figsize = (1000 / my_dpi, 1000 / my_dpi)
        out_file = os.path.join(output, '.'.join([prefix, sample, 'pdf']))
        plt.savefig(out_file, figsize=figsize, dpi=my_dpi)
        plt.close()


@mod2_pcoa.command()
@click.option('--distmat', required=True, type=click.Path(resolve_path=True,
                                                          readable=True,
                                                          exists=True),
              help='Input distance matrix to subsample nearest sample')
@click.option('--mapping_file', required=True, type=click.Path(
              resolve_path=True, readable=True, exists=True),
              help='Mapping file')
@click.option('--max', required=True, type=int,
              help='Max number of samples per category value')
@click.option('--category', required=True, type=str,
              help='The category to subsample in (likely COUNTRY)')
@click.option('--output', required=True, type=click.Path(exists=False,
              writable=True, resolve_path=True), help='Output file')
def subsample_dm(distmat, mapping_file, max, category, output):
    """Subsample the distmat to max samples per category value"""
    mf = pd.read_csv(mapping_file, '\t', converters=defaultdict(str),
                     index_col='#SampleID')
    id_to_cat = dict(mf[category])

    def bin_f(x):
        return id_to_cat[x]

    dm = read(distmat, into=DistanceMatrix)
    dm = dm.filter([id for _, id in isubsample(dm.ids, max, bin_f=bin_f)])
    dm.to_file(output)


@mod2_pcoa.command()
@click.option('--coords', required=True, type=click.Path(resolve_path=True,
              readable=True, exists=True), help='Coordinates file')
@click.option('--mapping_file', required=True, type=click.Path(
              resolve_path=True, readable=True, exists=True),
              help='Mapping file')
@click.option('--output', required=True, type=click.Path(exists=True,
              writable=True, resolve_path=True), help='Output directory')
@click.option('--prefix', required=True, type=str, help='Output file prefix')
@click.option('--samples', required=False, type=str,
              help='Comma separated list of samples to print')
@click.option('--distmat', required=True, type=click.Path(resolve_path=True,
                                                          readable=True,
                                                          exists=True),
              help=('Input distance matrix to find nearest sample (if not '
                    'present in the coordinates'))
def country(coords, mapping_file, output, prefix, samples, distmat):
    """Generates as many figures as samples in the coordinates file"""
    o = read(coords, into=OrdinationResults)
    o_id_lookup = set(o.site_ids)

    dm = read(distmat, into=DistanceMatrix)
    dm_id_lookup = {i: idx for idx, i in enumerate(dm.ids)}
    coord_samples_in_dm = {idx for idx, i in enumerate(dm.ids)
                           if i in o_id_lookup}

    # we'll be computing min values, so we need to avoid catching the  diagonal
    np.fill_diagonal(dm._data, np.inf)

    x, y = o.site[:, 0], o.site[:, 1]

    # coordinates
    c_df = pd.DataFrame(o.site, o.site_ids)

    # mapping file
    mf = pd.read_csv(mapping_file, '\t', converters=defaultdict(str),
                     index_col='#SampleID')
    # mf = mf.loc[o.site_ids]

    if samples is None:
        samples = dm.ids[:]
    else:
        samples = set(samples.split(',')).intersection(set(dm.ids))
        samples = mf.loc[samples].index

    color_Venezuela = sns.color_palette('Paired', 12)[10]
    color_Malawi = sns.color_palette('Paired', 12)[1]
    color_Western = sns.color_palette('Paired', 12)[4]
    color_Highlight = sns.color_palette('Paired', 12)[5]
    color_no_data = (0.5, 0.5, 0.5)

    grp_colors = OrderedDict()
    grp_colors['no_data'] = color_no_data
    grp_colors['Australia'] = color_Western
    grp_colors['Belgium'] = color_Western
    grp_colors['Canada'] = color_Western
    grp_colors['China'] = color_Western
    grp_colors['Finland'] = color_Western
    grp_colors['France'] = color_Western
    grp_colors['Germany'] = color_Western
    grp_colors['Great Britain'] = color_Western
    grp_colors['Ireland'] = color_Western
    grp_colors['Japan'] = color_Western
    grp_colors['Netherlands'] = color_Western
    grp_colors['New Zealand'] = color_Western
    grp_colors['Norway'] = color_Western
    grp_colors['Scotland'] = color_Western
    grp_colors['Spain'] = color_Western
    grp_colors['Switzerland'] = color_Western
    grp_colors['Thailand'] = color_Western
    grp_colors['United Arab Emirates'] = color_Western
    grp_colors['United Kingdom'] = color_Western
    grp_colors['United States of America'] = color_Western
    grp_colors['Malawi'] = color_Malawi
    grp_colors['Venezuela'] = color_Venezuela

    for sample_to_plot in samples:
        if sample_to_plot in o_id_lookup:
            sample = sample_to_plot
        else:
            # find the closest sample in the distance matrix that is in the
            # coordinates data
            sample = None
            for i in dm[dm_id_lookup[sample_to_plot]].argsort():
                if i in coord_samples_in_dm:
                    sample = dm.ids[i]
                    break

            # this should not ever happen
            if sample is None:
                raise ValueError("Unable to find a similar sample?")

        # countour plot superimposed
        sns.kdeplot(x, y, cmap='bone')
        sns.set_context(rc={"lines.linewidth": 0.75})

        # change particapant's country's color to color_Highlight unless
        # country is Venezuela or Malawi
        if (mf.loc[sample_to_plot]['COUNTRY'] != 'Malawi') & (
                mf.loc[sample_to_plot]['COUNTRY'] != 'Venezuela'):
            grp_colors[mf.loc[sample_to_plot]['COUNTRY']] = color_Highlight

        # plot each country except participant's according to colors above
        for grp, color in grp_colors.iteritems():
            if grp == mf.loc[sample_to_plot]['COUNTRY']:
                continue
            sub_coords = c_df[mf.COUNTRY == grp]
            plt.scatter(sub_coords[0], sub_coords[1], color=color,
                        edgecolor=np.asarray(color)*0.6, lw=LINE_WIDTH,
                        alpha=ALPHA)

        # now plot participant's country
        grp = mf.loc[sample_to_plot]['COUNTRY']
        color = grp_colors[grp]
        sub_coords = c_df[mf.COUNTRY == grp]
        plt.scatter(sub_coords[0], sub_coords[1], color=color,
                    edgecolor=np.asarray(color)*0.6, lw=LINE_WIDTH,
                    alpha=ALPHA)

        # plot participant's dot
        plt.scatter(c_df.loc[sample][0], c_df.loc[sample][1],
                    color=grp_colors[mf.loc[sample_to_plot]['COUNTRY']],
                    s=270, edgecolor='w', zorder=1, lw=LINE_WIDTH_WHITE)
        plt.scatter(c_df.loc[sample][0], c_df.loc[sample][1],
                    color=grp_colors[mf.loc[sample_to_plot]['COUNTRY']],
                    s=250, edgecolor=np.asarray(grp_colors[
                        mf.loc[sample_to_plot]['COUNTRY']])*0.6,
                    zorder=2, lw=LINE_WIDTH_BLACK)

        # reset particapant's country's color to color_Western unless country
        # is Venezuela or Malawi
        if (mf.loc[sample_to_plot]['COUNTRY'] != 'Malawi') & (
                mf.loc[sample_to_plot]['COUNTRY'] != 'Venezuela'):
            grp_colors[mf.loc[sample_to_plot]['COUNTRY']] = color_Western

        plt.axis('off')
        my_dpi = 72
        figsize = (1000 / my_dpi, 1000 / my_dpi)
        out_file = os.path.join(output, '.'.join([prefix, sample, 'pdf']))
        plt.savefig(out_file, figsize=figsize, dpi=my_dpi)
        plt.close()


@mod2_pcoa.command()
@click.option('--coords', required=True, type=click.Path(resolve_path=True,
              readable=True, exists=True), help='Coordinates file')
@click.option('--mapping_file', required=True, type=click.Path(
              resolve_path=True, readable=True, exists=True),
              help='Mapping file')
@click.option('--color', required=True, type=str,
              help='Metadata category to set color by')
@click.option('--output', required=True, type=click.Path(exists=True,
              writable=True, resolve_path=True), help='Output directory')
@click.option('--prefix', required=True, type=str, help='Output file prefix')
@click.option('--samples', required=False, type=str,
              help='Comma separated list of samples to print')
def gradient(coords, mapping_file, color, output, prefix, samples):
    """Generates as many figures as samples in the coordinates file"""
    o = read(coords, into=OrdinationResults)

    # coordinates
    c_df = pd.DataFrame(o.site, o.site_ids)

    # mapping file
    mf = pd.read_csv(mapping_file, '\t', converters=defaultdict(str),
                     index_col='#SampleID')
    mf = mf.loc[o.site_ids]
    mf[color] = mf[color].convert_objects(convert_numeric=True)

    if samples is None:
        samples = mf.index
    else:
        samples = set(samples.split(',')).intersection(set(o.site_ids))
        samples = mf.loc[samples].index

    numeric = mf[~pd.isnull(mf[color])]
    non_numeric = mf[pd.isnull(mf[color])]

    color_array = plt.cm.RdBu(numeric[color]/max(numeric[color]))

    for sample in samples:

        # plot numeric metadata as colored gradient
        ids = numeric.index
        x, y = c_df.loc[ids][0], c_df.loc[ids][1]
        plt.scatter(x, y, c=numeric[color], cmap=plt.get_cmap('RdBu'),
                    alpha=ALPHA, lw=LINE_WIDTH, edgecolor=color_array*0.6)

        # plt.colorbar()

        # plot non-numeric metadata as gray
        ids = non_numeric.index
        x, y = c_df.loc[ids][0], c_df.loc[ids][1]
        plt.scatter(x, y, c='0.5', alpha=ALPHA, lw=LINE_WIDTH, edgecolor='0.3')

        # plot individual's dot
        try:
            color_index = numeric.index.tolist().index(sample)
        except ValueError:
            color_index = None

        if color_index is None:
            _color = (0.5, 0.5, 0.5)
        else:
            _color = color_array[color_index]

        plt.scatter(c_df.loc[sample][0], c_df.loc[sample][1],
                    color=_color, s=270, edgecolor='w', lw=LINE_WIDTH_WHITE)
        plt.scatter(c_df.loc[sample][0], c_df.loc[sample][1],
                    color=_color, s=250, edgecolor=np.asarray(_color)*0.6,
                    lw=LINE_WIDTH_BLACK)

        plt.axis('off')
        my_dpi = 72
        figsize = (1000 / my_dpi, 1000 / my_dpi)
        out_file = os.path.join(output, '.'.join([prefix, sample, 'pdf']))
        plt.savefig(out_file, figsize=figsize, dpi=my_dpi)
        plt.close()

if __name__ == '__main__':
    mod2_pcoa()
