#!/usr/bin/env python

import numpy as np
from collections import defaultdict
from argparse import ArgumentParser

from biom.parse import parse_biom_table

from americangut.parse import parse_mapping_file_to_dict


__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Daniel McDonald", "Justine Debelius"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"


def table_by_category_value(table, cat):
    """Get tables by category value

    Tables are filtered to the respective category value, dropping all samples
    that do not hold the category value. The tables and the order of the
    category values is returned.

    Parameters
    ----------

    table : Table, the table to operate on
    cat : str, the category to operate on

    Returns
    -------
    tuple : (tables, values) where tables is a list of Table and values is
        a list of str
    """
    values = sorted({md[cat] for md in table.SampleMetadata})
    tables = []
    for val in values:
        tables.append(table.filterSamples(lambda v, i, md: md[cat] == val))
    return tables, values


def category_value_lookup(table, categories):
    """Construct a lookup for each category and value of interest

    Parameters
    ----------

    table : Table, the table to operate on
    categories : list of str, the categories to operate on

    Returns
    -------

    dict of dict : the outer dict is keyed by the category and valued by a
        dict. the inner dict is keyed by the category value and valued by the
        corresponding table
    """
    lookup = defaultdict(dict)
    for cat in categories:
        tables, order = table_by_category_value(table, cat)
        for t, o in zip(tables, order):
            lookup[cat][o] = t
    return lookup


def drop_sample(table, sample):
    """Drop a sample from a table

    Parameters
    ----------

    table : Table, the table to operate on
    sample : str, a sample ID

    Returns
    -------

    Table : the table without the sample
    """
    return table.filterSamples(lambda v, i, md: i != sample)


def sort_by_taxa(vector, vector_order, desired_order):
    """Sort the vector by the desired order

    Parameters
    ----------

    vector : np.array, a vector of abundances
    vector_order : the order of the labels represented
    desired_order : the desired order of the vector

    Returns
    -------

    np.array : the abundances in the desired order
    """
    lookup = {v: i for i, v in enumerate(vector_order)}
    order = np.array([lookup[v] for v in desired_order])
    return vector[order]


def average_per_observation(table):
    """Returns the average value of an observation in the table

    Parameters
    ----------

    table : Table, the table to operate on

    Returns
    -------

    np.array : the average abundance per observation
    """
    normalized_sums = table.sum('observation')
    mean_values = normalized_sums / len(table.SampleIds)

    if not np.isclose(mean_values.sum(), 1.0):
        raise ValueError("Averages do not appear to sum to 1!")

    return mean_values


def most_common_taxa(table, n):
    """Determine the most common taxa by relative abundance in the table

    Parameters
    ----------

    table : Table, the table to operate on
    n : int, the number of taxa to return

    Returns
    -------
    list : the n most common taxa sorted by abundance
    """
    obs_ids = np.array(table.ObservationIds)
    avg = np.array([v.mean() for v in table.iterObservationData()])
    sort_order = np.argsort(avg)[::-1]

    return obs_ids[sort_order][:n].tolist()


def make_collapse_f(taxa):
    """Group taxa

    This will create a collapsing function such that the taxa of interest will
    be represented by a single observation, and all other taxa will be lumped
    into a single "other" observation

    Parameters
    ----------

    taxa : iterable of str, the tax to care about

    Returns
    -------

    function : a function that can be applied to
        Table.collapseObservationsByMetadata
    """
    taxa = set(taxa)
    def f(md):
        if md['id'] in taxa:
            return md['id']
        else:
            return 'other'
    return f


def _make_plots_debug_printer(name, vector):
    out = [name]
    out.extend(map(str, vector))
    print '\t'.join(out)


def make_plot(sample, sample_vector, cat_vectors, categories, taxa_order,
              identified_vector, identified_person, debug=True):
    """Construct a stacked bar chart the vectors of interest

    Parameters
    ----------

    sample : str, the sample ID
    sample_vector : np.array, the taxa abundances in the sample
    cat_vectors : list of np.array, the category abundances
    categories : list of str, the categories
    taxa_order : list of str, the order of the taxa in the vectors
    identified_vector : np.array or None, the abundances for the identified
    identified_person : str or None, the name of the identified
    debug : bool, if true, dump out a text representation of the results
    """
    if debug:
        _make_plots_debug_printer('taxa_order', taxa_order)
        _make_plots_debug_printer(sample, sample_vector)

        for cat, v in zip(categories, cat_vectors):
            _make_plots_debug_printer(cat, v)

        if identified_person is not None:
            _make_plots_debug_printer(identified_person, identified_vector)

        _make_plots_debug_printer("*****", ["****"] * len(taxa_order))
    else:
        raise ValueError("Graphics not supported yet")


COMMON_TAXA = []


def main(args):
    table = parse_biom_table(open(args.table))
    metadata, commends = parse_mapping_file_to_dict(open(args.metadata))
    identified_person = args.identified_person
    identified_sample_id = args.identified_sample_id
    categories = args.categories
    samples = args.samples

    # sanity check: verify our categories exist in the metadata
    for cat in categories:
        if cat not in metadata[metadata.keys()[0]]:
            raise ValueError("%s does not appear in the metadata!" % cat)

    # sanity check: verify the samples to operate on exist in the table
    for s in samples:
        if not table.sampleExists(s):
            raise ValueError("%s does not exist in the table!" % s)

    # determine what taxa we're interested in
    if args.use_existing_common_taxa:
        common_taxa = COMMON_TAXA
    else:
        common_taxa = most_common_taxa(table, args.n_most_common_taxa)
        common_taxa.append('other')

    if 'other' not in common_taxa:
        raise ValueError("other is not in our list of taxa!")

    # pack in the sample metadata. we're adding the observation ID as metadata
    # as the BIOM collapse method operates on metadata, and we are assuming
    # that the provided table to this script has already been run through
    # QIIME's summarize_taxa.py script.
    table.addSampleMetadata(metadata)
    table.addObservationMetadata({i: {'id': i} for i in table.ObservationIds})

    # collapse to our taxa of interested, with everything else in 'other'
    collapse_f = make_collapse_f(common_taxa)
    table = table.collapseObservationsByMetadata(collapse_f, norm=False,
                                                 min_group_size=1)

    # snag a view on our identified sample
    if identified_sample_id is not None:
        identified_vector = table.sampleData(identified_sample_id)
    else:
        identified_vector = None

    # construct our per-category, per-value lookup table
    cat_tables = category_value_lookup(table, categories)

    # for each sample to operate on
    for s in samples:
        # fetch the corresponding sample metadata
        sample_metadata = metadata[s]

        # determine the abundances for the full sample
        vector = table.sampleData(s)
        vector_order = table.ObservationIds
        full_vector = sort_by_taxa(vector, vector_order, common_taxa)

        # determine the averages of each category group
        cat_vectors = []
        for cat in categories:
            # get the category value and table (e.g., BMI_CATEGORY: overweight)
            cat_value = sample_metadata[cat]
            cat_table = cat_tables[cat][cat_value]

            # remove the sample from the category table
            cat_table = drop_sample(cat_table, s)

            # compute the averages per observation, and the observation order
            vector = average_per_observation(cat_table)
            vector_order = cat_table.ObservationIds

            # store the vector in sorted order
            cat_vectors.append(sort_by_taxa(vector, vector_order, common_taxa))

        # make a pretty picture
        make_plot(s, full_vector, cat_vectors, categories, common_taxa,
                  identified_vector, identified_person)


if __name__ == '__main__':
    parser = ArgumentParser(description='Stacked bar charts per individual '
                                        'and corresponding categories')
    parser.add_argument('-i', '--table', required=True,
                        help='A summarized taxa table')
    parser.add_argument('-m', '--metadata', required=True,
                        help='Sample metadata')
    parser.add_argument('--identified-person', required=False, default=None,
                        help='The name of an identified person')
    parser.add_argument('--identified-sample-id', required=False, default=None,
                        help='The sample ID of an identified person')
    parser.add_argument('--samples', required=True,
                        help='The samples to produce plots for')
    parser.add_argument('--categories', required=True,
                        help='The mapping file categories to use')
    parser.add_argument('--use-existing-common-taxa', required=False,
                        default=False, action='store_true',
                        help='Use the existing list of common taxa')
    parser.add_argument('--n-most-common-taxa', required=False,
                        default=None, help='Use the N most common taxa',
                        type=int)

    args = parser.parse_args()

    if args.use_existing_common_taxa and args.n_most_common_taxa is not None:
        raise parser.error("--use-existing-common-taxa is not compatible with "
                           "--n-most-common-taxa")

    args.samples = args.samples.split(',')
    args.categories = args.categories.split(',')

    main(args)
