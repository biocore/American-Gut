import os

from matplotlib import use, rcParams
use('Agg')  # noqa

import biom
import matplotlib.pyplot as plt
import pandas as pd
from qiime.util import qiime_system_call
import seaborn as sn

import americangut.util as agu
import americangut.notebook_environment as agenv
import americangut.results_utils as agru

# Sets up plotting parameters so that the default setting is use to Helvetica
# in plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']


def create_opts(sample_type, chp_path, gradient_color_by, barchart_categories):
    """Create a dict of options for processing functions

    Parameters
    ----------
    sample_type : str, {fecal, oral, skin}
        The sample type.

    chp_path : str
        The path to the current processing chapter.

    gradient_color_by : str
        The taxon to color by for the gradient PCoA.

    barchart_categories : Iterable of str
        The relevant barchart categories to use.

    Returns
    -------
    dict
        The paths and desired options.
    """
    def attach_dict(d, items):
        return [(d, i) for i in items]

    opts = {}
    queue = attach_dict(opts, agenv.paths.items())

    while queue:
        d, (key, value) = queue.pop()
        if isinstance(value, dict):
            # get the relevant dict or a new dict
            nested_d = d.get(key, {})
            d[key] = nested_d

            queue.extend(attach_dict(nested_d, value.items()))
        else:
            d[key] = agu.get_path(value)

    opts['sample_type'] = sample_type
    opts['gradient_color_by'] = gradient_color_by
    opts['chp-path'] = chp_path
    opts['rarefaction-depth'] = int(agenv.get_rarefaction_depth()[0])

    barchart_map = {'diet': 'DIET_TYPE',
                    'sex':  'SEX',
                    'age':  'AGE_CAT',
                    'bmi':  'BMI_CAT',
                    'flossing': 'FLOSSING_FREQUENCY',
                    'cosmetics': 'COSMETICS_FREQUENCY',
                    'hand': 'DOMINANT_HAND'}

    fmt = []
    for cat in barchart_categories:
        if cat not in barchart_map:
            raise KeyError("%s is not known" % cat)

        table_key = 'ag-%s-%s-biom' % (sample_type, cat)
        cat_path = opts['collapsed']['100nt']['1k'][table_key]
        fmt.append("%s:%s" % (barchart_map[cat], cat_path))

    opts['barchart_categories'] = '"%s"' % ', '.join(fmt)
    return opts


def sample_type_processor(functions, opts, ids):
    """Execute a series of functions per sample

    Parameters
    ----------
    functions : list of function
        A list of functions where each function is expected to have the
        following signature:

        dict <- f(dict, list)

        The passed in list is expected to contain a set of sample IDs to
        operate on. The passed in dict contains path information and processing
        specific details (e.g., taxon to use for gradient coloring).  The
        return from a function is expected to be a dict in which every input
        sample ID from the passed in list is represented as a key, and the
        value is a list of str where each str indicates an error occured -- an
        empty list would indicate no observed errors for the sample.

    opts : dict
        Paths and relevant processing options.

    ids : Iterable of str
        The list of sample IDs to examine

    Returns
    -------
    dict
        A dict keyed by sample ID and valued by a list. The list contains
        all errors observed for the sample, or the empty list if no errors
        were observed.
    """
    return merge_error_reports(*[f(opts, ids) for f in functions])


def _result_path(opts, id_):
    """Form a ID specific result path"""
    return os.path.join(opts['per-sample']['results'], id_)


def _base_barcode(id_):
    """Get the actual barcode from a sample ID"""
    # NOTE: old ebi accessions and processing do not follow this format.
    # expectation is <study_id>.<barcode>

    return id_.split('.')[1]


def partition_samples_by_bodysite(md, site_to_function):
    """Yield the processing function and associated sample IDs

    Parameters
    ----------
    md : pd.DataFrame
        The metadata to examine. It is expected that this DataFrame contain a
        column named SIMPLE_BODY_SITE.
    site_to_function : list of tuple
        A mapping of a body site category value to a processing function.

    Returns
    -------
    generator
        function
            The processing functione.
        Iterable of str
            The sample IDs to process.
    """
    for site, func in site_to_function:
        df_subset = md[md.SIMPLE_BODY_SITE == site]
        yield func, df_subset.index


def merge_error_reports(*reports):
    """Merge error reports

    Parameters
    ----------
    *reports : list of dict
        The error reports

    Returns
    -------
    dict
        Keyed by sample ID, valued by the list of observed errors. An empty
        list is associted with a sample ID if no errors were observed.
        {str: [str]}
    """
    result = {}

    for report in reports:
        for id_, value in report.items():
            if id_ not in result:
                result[id_] = []

            if value is not None:
                result[id_].append(value)

    return result


def _iter_ids_over_system_call(cmd_fmt, sample_ids, opts):
    """Iteratively execute a system call over sample IDs

    Parameters
    ----------
    cmd_fmt : str
        The format of the command to execute. It is expected that there
        is a single string format to be done and it should take a sample
        ID
    sample_ids : Iterable of str
        A list of sample IDs of interest

    Returns
    -------
    dict
        A dict containing each sample ID and any errors observed or None if
        no error was observed for the sample. {str: str or None}
    """
    results = {}

    for id_ in sample_ids:
        cmd = cmd_fmt % {'result_path': _result_path(opts, id_),
                         'id': id_}
        stdout, stderr, return_value = qiime_system_call(cmd)

        if return_value != 0:
            msg = stderr.splitlines()
            results[id_] = 'FAILED (%s): %s' % (msg[-1] if msg else '', cmd)
        else:
            results[id_] = None

    return results


def taxa_summaries(opts, sample_ids):
    """Produce digestable taxonomy summaries per sample

    Parameters
    ----------
    opts : dict
        A dict of relevant opts.

    sample_ids : Iterable of str
        A list of sample IDs of interest

    Returns
    -------
    dict
        A dict containing each sample ID and any errors observed or None if
        no error was observed for the sample. {str: str or None}
    """
    results = {}
    path = opts['taxa']['notrim']['L6']['ag-%s-biom' % opts['sample_type']]
    site_table = biom.load_table(path)
    table_taxon_ids = site_table.ids(axis='observation')

    for id_ in sample_ids:
        if not site_table.exists(id_):
            results[id_] = 'ID not found'
        else:
            results[id_] = None
            taxa_path = os.path.join(_result_path(opts, id_),
                                     '%s.txt' % id_)

            with open(taxa_path, 'w') as fp:
                fp.write("#taxon\trelative_abundance\n")

                v = site_table.data(id_, dense=True)
                for sorted_v, taxa in sorted(zip(v, table_taxon_ids))[::-1]:
                    if sorted_v:
                        fp.write("%s\t%f\n" % (taxa, sorted_v))
    return results


def alpha_plot(opts, sample_ids):
    """Produces digestable alpha diversity distribution plots per sample

    Parameters
    ----------
    opts : dict
        A dict of relevant opts.
    sample_ids : Iterable of str
        A list of sample IDs of interest

    Returns
    -------
    dict
        A dict containing each sample ID and any errors observed or None if
        no error was observed for the sample. {str: str or None}
    """

    results = {}
    alpha_map = pd.read_csv(
        agu.get_existing_path(opts['collapsed']['100nt']['alpha-map']),
        sep='\t',
        dtype=str,
        )

    alpha_metrics = ['shannon_1k', 'PD_whole_tree_1k']

    # Checks the alpha_field is in the mapping file
    for metric in alpha_metrics:
        if metric not in alpha_map.columns:
            raise ValueError('%s is not a valid alpha diversity field name.'
                             % metric)
    # Checks the group_field is in the mapping file
    if 'SIMPLE_BODY_SITE' not in alpha_map.columns:
        raise ValueError('SIMPLE_BODY_SITE is not a valid field name.')

    alpha_map[alpha_metrics] = alpha_map[alpha_metrics].astype(float)
    alpha_map.set_index('#SampleID', inplace=True)

    results = {}
    for id_ in sample_ids:
        if id_ not in alpha_map.index:
            results[id_] = 'ID not found'
        else:
            results[id_] = None
            shannon_path = os.path.join(_result_path(opts, id_),
                                        'shannon_%s.png' % id_)
            _plot_alpha(id_, alpha_map, 'shannon_1k',
                        xlabel='Shannon Diversity',
                        fp=shannon_path)

            # Generates the pd whole tree diversity figure
            pd_path = os.path.join(_result_path(opts, id_),
                                   'pd_%s.png' % id_)
            _plot_alpha(id_, alpha_map, 'PD_whole_tree_1k',
                        xlabel='PD Whole Tree Diversity',
                        fp=pd_path)

    return results


def sufficient_sequence_counts(opts, sample_ids):
    """Errors if the sequence counts post filtering are < 1000

    Parameters
    ----------
    opts : dict
        A dict of relevant opts.

    sample_ids : Iterable of str
        A list of sample IDs of interest

    Returns
    -------
    dict
        A dict containing each sample ID and any errors observed or None if
        no error was observed for the sample. {str: str or None}
    """
    results = {}
    table = biom.load_table(opts['otus']['100nt']['ag-biom'])

    minimum_depth = opts['rarefaction-depth']

    for id_ in sample_ids:
        results[id_] = None

        if table.exists(id_):
            counts = table.data(id_).sum()
            if counts < minimum_depth:
                results[id_] = '%d seqs after filtering for blooms' % counts
        else:
            results[id_] = '0 seqs after filtering for blooms'

    return results


def taxon_significance(opts, sample_ids):
    """Produce OTU significance results

    Parameters
    ----------
    opts : dict
        A dict of relevant opts.
    sample_ids : Iterable of str
        A list of sample IDs of interest

    Returns
    -------
    dict
        A dict containing each sample ID and any errors observed or None if
        no error was observed for the sample. {str: str or None}
    """
    path = opts['taxa']['notrim']['L6']['ag-%s-biom' % opts['sample_type']]
    map_path = opts['meta']['ag-cleaned-md']

    cmd_fmt = 'generate_otu_signifigance_tables_AGP.py -i %s ' % path
    cmd_fmt += '-o %(result_path)s -s %(id)s '
    cmd_fmt += '-m %s' % map_path

    return _iter_ids_over_system_call(cmd_fmt, sample_ids, opts)


def body_site_pcoa(opts, sample_ids):
    """Produce the per-sample all body site PCoA

    Parameters
    ----------
    opts : dict
        A dict of relevant opts.
    sample_ids : Iterable of str
        A list of sample IDs of interest

    Returns
    -------
    dict
        A dict containing each sample ID and any errors observed or None if
        no error was observed for the sample. {str: str or None}
    """
    coords = opts['beta']['100nt']['1k']['ag-pgp-hmp-gg-unifrac-pc']
    mapping = opts['meta']['ag-pgp-hmp-gg-cleaned-md']
    cmd_fmt = ' '.join(["mod2_pcoa.py body_site",
                        "--coords %s" % coords,
                        "--mapping_file %s" % mapping,
                        "--output %(result_path)s",
                        "--filename figure1.pdf",
                        "--sample %(id)s"])

    return _iter_ids_over_system_call(cmd_fmt, sample_ids, opts)


def country_pcoa(opts, sample_ids):
    """Produce the per-sample subsampled country PCoA

    Parameters
    ----------
    opts : dict
        A dict of relevant opts.
    sample_ids : Iterable of str
        A list of sample IDs of interest

    Returns
    -------
    dict
        A dict containing each sample ID and any errors observed or None if
        no error was observed for the sample. {str: str or None}
    """
    beta1k = opts['beta']['100nt']['1k']
    coords = beta1k['ag-gg-subsampled-unifrac-pc']
    cmd_fmt = ' '.join(["mod2_pcoa.py country",
                        "--distmat %s" % beta1k['ag-gg-unifrac'],
                        "--coords %s" % coords,
                        "--mapping_file %s" % opts['meta']['ag-gg-cleaned-md'],
                        "--output %(result_path)s",
                        "--filename figure2.pdf",
                        "--sample %(id)s"])

    return _iter_ids_over_system_call(cmd_fmt, sample_ids, opts)


def gradient_pcoa(opts, sample_ids):
    """Produce a gradient PCoA

    Parameters
    ----------
    opts : dict
        A dict of relevant opts.
    sample_ids : Iterable of str
        A list of sample IDs of interest

    Returns
    -------
    dict
        A dict containing each sample ID and any errors observed or None if
        no error was observed for the sample. {str: str or None}
    """
    mapping = opts['taxa']['notrim']['L2']['ag-md']
    coords_key = 'ag-%s-unifrac-pc' % opts['sample_type']
    coords = opts['beta']['100nt']['1k'][coords_key]
    cmd_fmt = ' '.join(["mod2_pcoa.py gradient",
                        "--coords %s" % coords,
                        "--mapping_file %s" % mapping,
                        "--output %(result_path)s",
                        "--filename figure3.pdf",
                        "--color %s" % opts['gradient_color_by'],
                        "--sample %(id)s"])

    return _iter_ids_over_system_call(cmd_fmt, sample_ids, opts)


def pie_plot(opts, sample_ids):
    """Produce a pie chart

    Parameters
    ----------
    opts : dict
        A dict of relevant opts.
    sample_ids : Iterable of str
        A list of sample IDs of interest

    Returns
    -------
    dict
        A dict containing each sample ID and any errors observed or None if
        no error was observed for the sample. {str: str or None}
    """
    cmd_fmt = ' '.join(['make_pie_plot_AGP.py',
                        '-i %s' % opts['taxa']['notrim']['L3']['ag-tsv'],
                        '-o %(result_path)s',
                        '-s %(id)s'])
    return _iter_ids_over_system_call(cmd_fmt, sample_ids, opts)


def bar_chart(opts, sample_ids):
    """Produce a bar chart

    Parameters
    ----------
    opts : dict
        A dict of relevant opts.
    sample_ids : Iterable of str
        A list of sample IDs of interest

    Returns
    -------
    dict
        A dict containing each sample ID and any errors observed or None if
        no error was observed for the sample. {str: str or None}
    """
    sample_type = opts['sample_type']
    path = opts['collapsed']['notrim']['1k']['ag-%s-biom' % sample_type]
    cmd_fmt = ' '.join(['make_phyla_plots_AGP.py',
                        '-i %s' % path,
                        '-m %s' % opts['meta']['ag-cleaned-md'],
                        '-o %(result_path)s',
                        '-c %s' % opts['barchart_categories'],
                        '-t %s' % opts['sample_type'],
                        '-s %(id)s'])
    return _iter_ids_over_system_call(cmd_fmt, sample_ids, opts)


def per_sample_directory(opts, sample_ids):
    """Create a per sample directory

    Parameters
    ----------
    opts : dict
        A dict of relevant opts.
    sample_ids : Iterable of str
        A list of sample IDs of interest

    Returns
    -------
    dict
        A dict containing each sample ID and any errors observed or None if
        no error was observed for the sample. {str: str or None}
    """
    result = {}

    for id_ in sample_ids:
        result[id_] = None

        path = _result_path(opts, id_)

        if not os.path.exists(path):
            try:
                os.mkdir(path)
            except Exception as e:
                result[id_] = e.args[0]

    return result


def stage_per_sample_specific_statics(opts, sample_ids):
    """Items like the Latex templates

    Parameters
    ----------
    opts : dict
        A dict of relevant opts.
    sample_ids : Iterable of str
        A list of sample IDs of interest

    Returns
    -------
    dict
        A dict containing each sample ID and any errors observed or None if
        no error was observed for the sample. {str: str or None}
    """
    result = {}
    sample_type = opts['sample_type'].lower()
    statics_src = opts['per-sample']['statics-%s' % sample_type]
    for id_ in sample_ids:
        result[id_] = None
        path = _result_path(opts, id_)
        template_path = os.path.join(path, id_ + '.tex')
        statics_path = os.path.join(path, 'statics')

        try:
            agru.stage_static_latex(opts['sample_type'], template_path)
        except:
            result[id_] = "Cannot stage template."
            continue

        try:
            os.symlink(statics_src, statics_path)
        except:
            result[id_] = "Cannot symlink for statics."

    return result


def _plot_alpha(sample, alpha_map, alpha_field, group_field='SIMPLE_BODY_SITE',
                output_dir=None, xlabel=None, fp=None, debug=False):
    """Generates a distrbution plot for the data

    Parameters
    ----------
    sample : str
        The sample ID to be plotted
    alpha_map_fp : pandas DataFrame
        A pandas dataframe containing the sample metadata. The sample ID
        should be given in the `'#SampleID'` column, a column with
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
    If the sample is present, a matplotlib figure with the alpha diversity
    distribution and a line indicating the sample value is returned. If a
    file path is specified, the figure will be saved at the filepath instead
    of returning.

    If debug is passed, the following parameters are returned:
        group : str
            The value of the `group_field` for the sample
        group_alpha : ndarray
            The alpha diversity values associated with the group
        sample_alpha : float
            The alpha diversity for the sample
        xlabel : str
            The label used for the x-axis of the plot.

    """

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

    with sn.axes_style('ticks', {'axes.facecolor': 'none'}):
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
        # Removes the spine
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

    if fp is None:
        return fig
    else:
        fig.savefig(fp, dpi=300)
        fig.clear()
