#!/usr/bin/env python

import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from skbio.stats.ordination import OrdinationResults
from collections import defaultdict

ALPHA = 1.0
LINE_WIDTH = 0.3

@click.group()
def mod1_pcoa():
    pass

@mod1_pcoa.command()
@click.option('--coords', required=True, type=click.Path(resolve_path=True, readable=True, exists=True), help='Coordinates file')
@click.option('--mapping_file', required=True, type=click.Path(resolve_path=True, readable=True, exists=True), help='Mapping file')
def body_site(coords,mapping_file):
    """Generates as many figures as samples in the coordinates file (Figure 1)"""
    o = OrdinationResults.from_file(coords)
    x, y = o.site[:,0], o.site[:,1]

    # coordinates
    c_df = pd.DataFrame(o.site, o.site_ids)

    # mapping file
    mf = pd.read_csv(mapping_file, '\t', converters=defaultdict(str), index_col='#SampleID')
    mf = mf.loc[o.site_ids]

    color_hmp_fecal = sns.color_palette('Paired', 12)[10] # light brown
    color_agp_fecal = sns.color_palette('Paired', 12)[11] # dark brown
    color_hmp_oral  = sns.color_palette('Paired', 12)[0]  # light blue
    color_agp_oral  = sns.color_palette('Paired', 12)[1]  # dark blue
    color_hmp_skin  = sns.color_palette('Paired', 12)[2]  # light green
    color_agp_skin  = sns.color_palette('Paired', 12)[3]  # dark green

    cat_colors = {'HMP-FECAL': color_hmp_fecal,
                  'GG-FECAL':  color_hmp_fecal,
                  'PGP-FECAL': color_hmp_fecal, 
                  'AGP-FECAL': color_agp_fecal, 
                  'HMP-ORAL':  color_hmp_oral,
                  'PGP-ORAL':  color_hmp_oral,
                  'AGP-ORAL':  color_agp_oral, 
                  'HMP-SKIN':  color_hmp_skin,
                  'PGP-SKIN':  color_hmp_skin,
                  'AGP-SKIN':  color_agp_skin}

    for sample in mf.index:
        for cat, color in cat_colors.iteritems():
            sub_coords = c_df[mf.TITLE_BODY_SITE == cat]

            plt.scatter(sub_coords[0], sub_coords[1], color=color, edgecolor='k', lw=LINE_WIDTH, alpha=ALPHA)

        plt.scatter(c_df.loc[sample][0], c_df.loc[sample][1],
                    color=cat_colors[mf.loc[sample]['TITLE_BODY_SITE']],
                    s=270, edgecolor='w')
        plt.scatter(c_df.loc[sample][0], c_df.loc[sample][1],
                    color=cat_colors[mf.loc[sample]['TITLE_BODY_SITE']],
                    s=250, edgecolor='k')
        plt.axis('off')
        my_dpi = 72
        plt.savefig(sample+'.pdf', figsize=(1000/my_dpi, 1000/my_dpi), dpi=my_dpi)
        plt.close()

@mod1_pcoa.command()
@click.option('--coords', required=True, type=click.Path(resolve_path=True, readable=True, exists=True), help='Coordinates file')
@click.option('--mapping_file', required=True, type=click.Path(resolve_path=True, readable=True, exists=True), help='Mapping file')
def country(coords,mapping_file):
    """Generates as many figures as samples in the coordinates file (Figure 1)"""
    o = OrdinationResults.from_file(coords)
    x, y = o.site[:,0], o.site[:,1]

    # coordinates
    c_df = pd.DataFrame(o.site, o.site_ids)

    # mapping file
    mf = pd.read_csv(mapping_file, '\t', converters=defaultdict(str), index_col='#SampleID')
    mf = mf.loc[o.site_ids]

    color_Venezuela = sns.color_palette('Paired', 12)[10]
    color_Malawi    = sns.color_palette('Paired', 12)[1]
    color_Western   = sns.color_palette('Paired', 12)[4]
    color_Highlight = sns.color_palette('Paired', 12)[5]
    color_no_data   = '0.5'

    cat_colors = {'Australia':                color_Western,
                  'Belgium':                  color_Western,
                  'Canada':                   color_Western,
                  'China':                    color_Western,
                  'Finland':                  color_Western,
                  'France':                   color_Western,
                  'Germany':                  color_Western,
                  'Great Britain':            color_Western,
                  'Ireland':                  color_Western,
                  'Japan':                    color_Western,
                  'Malawi':                   color_Malawi,
                  'Netherlands':              color_Western,
                  'New Zealand':              color_Western,
                  'Norway':                   color_Western,
                  'Scotland':                 color_Western,
                  'Spain':                    color_Western,
                  'Switzerland':              color_Western,
                  'Thailand':                 color_Western,
                  'United Arab Emirates':     color_Western,
                  'United Kingdom':           color_Western,
                  'United States of America': color_Western,
                  'Venezuela':                color_Venezuela,
                  'no_data':                  color_Western}

    for sample in mf.index:

        # countour plot superimposed
        sns.kdeplot(x, y, cmap='bone')
        sns.set_context(rc={"lines.linewidth": 0.75})

        # change particapant's country's color to color_Highlight unless country is Venezuela or Malawi
        if (mf.loc[sample]['COUNTRY'] != 'Malawi') & (mf.loc[sample]['COUNTRY'] != 'Venezuela'):
            cat_colors[mf.loc[sample]['COUNTRY']] = color_Highlight

        # plot each country according to colors above
        for cat, color in cat_colors.iteritems():
            sub_coords = c_df[mf.COUNTRY == cat]
            plt.scatter(sub_coords[0], sub_coords[1], color=color, edgecolor='k', lw=LINE_WIDTH, alpha=ALPHA)

        # plot participant's dot
        plt.scatter(c_df.loc[sample][0], c_df.loc[sample][1],
                    color=cat_colors[mf.loc[sample]['COUNTRY']],
                    s=270, edgecolor='w', zorder=1)
        plt.scatter(c_df.loc[sample][0], c_df.loc[sample][1],
                    color=cat_colors[mf.loc[sample]['COUNTRY']],
                    s=250, edgecolor='k', zorder=2)      
        
        # reset particapant's country's color to color_Western unless country is Venezuela or Malawi
        if (mf.loc[sample]['COUNTRY'] != 'Malawi') & (mf.loc[sample]['COUNTRY'] != 'Venezuela'):
            cat_colors[mf.loc[sample]['COUNTRY']] = color_Western

        plt.axis('off')
        my_dpi = 72
        plt.savefig(sample+'.pdf', figsize=(1000/my_dpi, 1000/my_dpi), dpi=my_dpi)
        plt.close()

@mod1_pcoa.command()
@click.option('--coords', required=True, type=click.Path(resolve_path=True, readable=True, exists=True), help='Coordinates file')
@click.option('--mapping_file', required=True, type=click.Path(resolve_path=True, readable=True, exists=True), help='Mapping file')
@click.option('--color', required=True, type=str, help='Metadata category to set color by')
def gradient(coords,mapping_file, color):
    """Generates as many figures as samples in the coordinates file (Figures 2 & 3)"""
    o = OrdinationResults.from_file(coords)
    c_df = pd.DataFrame(o.site, o.site_ids)

    # mapping file
    mf = pd.read_csv(mapping_file, '\t', converters=defaultdict(str), index_col='#SampleID')
    mf = mf.loc[o.site_ids]
    mf[color] = mf[color].convert_objects(convert_numeric=True)

    numeric = mf[~pd.isnull(mf[color])]
    non_numeric = mf[pd.isnull(mf[color])]

    color_array = plt.cm.RdBu(numeric[color]/max(numeric[color]))

    for sample in mf.index:

        # NUMERIC
        ids = numeric.index
        x, y = c_df.loc[ids][0], c_df.loc[ids][1]
        plt.scatter(x, y, c=numeric[color], cmap=plt.get_cmap('RdBu'),
                    alpha=ALPHA, lw=LINE_WIDTH)

        #plt.colorbar()

        # NON-NUMERIC
        ids = non_numeric.index
        x, y = c_df.loc[ids][0], c_df.loc[ids][1]
        plt.scatter(x, y, c='0.5', alpha=ALPHA, lw=LINE_WIDTH)

        # INDIVIDUAL BIG DOT
        try:
            color_index = numeric.index.tolist().index(sample)
        except ValueError:
            color_index = None

        if color_index is None:
            _color = '0.5'
        else:
            _color = color_array[color_index]

        plt.scatter(c_df.loc[sample][0], c_df.loc[sample][1],
            color=_color,
            s=270, edgecolor='w')
        plt.scatter(c_df.loc[sample][0], c_df.loc[sample][1],
            color=_color,
            s=250, edgecolor='k')

        plt.axis('off')
        my_dpi = 72
        plt.savefig(sample+'.pdf', figsize=(1000/my_dpi, 1000/my_dpi), dpi=my_dpi)
        plt.close()

@mod1_pcoa.command()
@click.option('--coords', required=True, type=click.Path(resolve_path=True, readable=True, exists=True), help='Coordinates file')
@click.option('--mapping_file', required=True, type=click.Path(resolve_path=True, readable=True, exists=True), help='Mapping file')
@click.option('--color', required=True, type=str, help='Metadata category to set color by')
def double_gradient(coords,mapping_file, color):
    """Generates as many figures as samples in the coordinates file (Figures 2 & 3)"""
    o = OrdinationResults.from_file(coords)
    c_df = pd.DataFrame(o.site, o.site_ids)

    # mapping file
    mf = pd.read_csv(mapping_file, '\t', converters=defaultdict(str), index_col='#SampleID')
    mf = mf.loc[o.site_ids]
    mf[color] = mf[color].convert_objects(convert_numeric=True)

    numeric = mf[~pd.isnull(mf[color])]
    non_numeric = mf[pd.isnull(mf[color])]

    numeric_UK = numeric[numeric.COUNTRY == 'United Kingdom']
    numeric_Other = numeric[numeric.COUNTRY != 'United Kingdom']

    color_array_UK = plt.cm.spring(numeric_UK[color]/max(numeric_UK[color]))
    color_array_Other = plt.cm.winter(numeric_Other[color]/max(numeric_Other[color]))

    for sample in mf.index:

        # NUMERIC OTHER
        ids = numeric_Other.index
        x, y = c_df.loc[ids][0], c_df.loc[ids][1]
        plt.scatter(x, y, c=numeric_Other[color], cmap=plt.get_cmap('winter'),
                    alpha=ALPHA, lw=LINE_WIDTH)

        # NUMERIC UK
        ids = numeric_UK.index
        x, y = c_df.loc[ids][0], c_df.loc[ids][1]
        plt.scatter(x, y, c=numeric_UK[color], cmap=plt.get_cmap('spring'),
                    alpha=ALPHA, lw=LINE_WIDTH)

        #plt.colorbar()

        # NON-NUMERIC
        ids = non_numeric.index
        x, y = c_df.loc[ids][0], c_df.loc[ids][1]
        plt.scatter(x, y, c='0.5', alpha=ALPHA, lw=LINE_WIDTH)

        # INDIVIDUAL BIG DOT
        try:
            if mf.loc[sample].COUNTRY == 'United Kingdom':
                color_index = numeric_UK.index.tolist().index(sample)
            else:
                color_index = numeric_Other.index.tolist().index(sample)
        except ValueError:
            color_index = None

        if color_index is None:
            _color = '0.5'
        else:
            if mf.loc[sample].COUNTRY == 'United Kingdom':
                _color = color_array_UK[color_index]
            else:
                _color = color_array_Other[color_index]

        plt.scatter(c_df.loc[sample][0], c_df.loc[sample][1],
            color=_color,
            s=270, edgecolor='w')
        plt.scatter(c_df.loc[sample][0], c_df.loc[sample][1],
            color=_color,
            s=250, edgecolor='k')

        plt.axis('off')
        my_dpi = 72
        plt.savefig(sample+'.pdf', figsize=(1000/my_dpi, 1000/my_dpi), dpi=my_dpi)
        plt.close()

if __name__ == '__main__':
    mod1_pcoa()
