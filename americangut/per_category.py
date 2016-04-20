from os.path import join
from copy import copy

from biom import load_table

import americangut.notebook_environment as agenv
from .util import collapse_full, collapse_taxonomy, get_existing_path


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
