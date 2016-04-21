from os.path import join
from copy import copy

from biom import load_table
import pandas as pd

import americangut.notebook_environment as agenv
from .util import collapse_full, collapse_taxonomy, get_existing_path
from .results_utils import plot_alpha


def cat_taxa_summaries():
    """Creates taxa summary files for each available summary category per site
    """
    paths = copy(agenv.paths['collapsed']['notrim']['1k'])
    out_dir = get_existing_path(
        agenv.path['populated-templates']['result-taxa'])
    del paths['ag-biom']
    for name, path in paths.items():
        # consistent naming as stool for all participant items
        name = name.replace('-biom', '').replace('fecal', 'stool')
        table = load_table(get_existing_path(path))
        if len(name.split('-')) == 2:
            # Have entire cohort of samples for site, so need to get averages
            table = collapse_full(table)
        table = collapse_taxonomy(table)
        ids = table.ids(axis='observation')
        for col in table.ids():
            if col == 'Unknown':
                continue
            cleaned_col = col.split('(')[0].strip().replace(' ', '_')
            filename = '-'.join([name, cleaned_col]) + '.txt'

            with open(join(out_dir, filename), 'w') as f:
                for otu, val in zip(ids, table.data(col)):
                    if val == 0.0:
                        continue
                    f.write('%s\t%s\n' % (otu, str(val)))


def cat_alpha_plots():
    out_dir = get_existing_path(
        agenv.paths['populated-templates']['alpha-files'])
    alpha_map = pd.read_csv(
        get_existing_path(agenv.paths['collapsed']['100nt']['alpha-map']),
        sep='\t', dtype=str)
    alpha_map.set_index('#SampleID', inplace=True)

    # Checks the group_field is in the mapping file, and groups sites
    if 'SIMPLE_BODY_SITE' not in alpha_map.columns:
        raise ValueError('SIMPLE_BODY_SITE is not a valid field name.')
    sites = set(alpha_map['SIMPLE_BODY_SITE'])

    alpha_metrics = ['shannon_1k', 'PD_whole_tree_1k']
    alpha_map[alpha_metrics] = alpha_map[alpha_metrics].astype(float)

    categories = ['AGE_CAT', 'BMI_CAT', 'SEX', 'DIET_TYPE', 'COSMETICS',
                  'DOMINANT_HAND']
    sample_color = '#1f78b4'

    for cat in categories:
        cat_groups = alpha_map.groupby(cat).groups
        for site in sites:
            # consistent naming as stool for all participant items
            site_name = site.lower().replace('fecal', 'stool')
            for group in list(cat_groups.keys()):
                sample_name = '-'.join([site_name, cat, group])
                pd_samps = alpha_map.loc[cat_groups[group], 'PD_whole_tree_1k']
                pd_mean = pd_samps.mean()
                pd_stdev = pd_samps.std()
                shannon_samps = alpha_map.loc[cat_groups[group], 'shannon_1k']
                sh_mean = shannon_samps.mean()
                sh_stdev = shannon_samps.std()

                # Add the new group row row to the metadata
                new_df = pd.DataFrame(
                    [site, sh_mean, pd_mean], index=[sample_name],
                    columns=['SIMPLE_BODY_SITE', 'shannon_1k',
                             'PD_whole_tree_1k'])
                new_map = alpha_map.concat(new_df)

                # Generates the shannon diversity figure
                shannon_path = join(out_dir, 'shannon_%s-%s-%s.png' % (
                    site_name, cat, group))
                plot_alpha(sample_name, new_map, 'shannon_1k',
                           xlabel='Shannon Diversity',
                           fp=shannon_path,
                           sample_color=sample_color,
                           highlight_range=[sh_mean-sh_stdev,
                                            sh_mean+sh_stdev],
                           categorical=True)

                # Generates the pd whole tree diversity figure
                pd_path = join(out_dir, 'pd_%s-%s-%s.png' %
                               (site_name, cat, group))
                plot_alpha(sample_name, new_map, 'PD_whole_tree_1k',
                           xlabel='PD Whole Tree Diversity',
                           fp=pd_path,
                           sample_color=sample_color,
                           highlight_range=[pd_mean-pd_stdev,
                                            pd_mean+pd_stdev],
                           categorical=True)
