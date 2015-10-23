import os
import shutil

import biom
from qiime.util import qiime_system_call

import americangut.util as agu
import americangut.notebook_environment as agenv


def create_opts(sample_type, chp_path, gradient_color_by, barchart_categories):
    """Create a dict of options for processing functions

    Paramters
    ---------
    sample_type : str, {fecal, oral, skin}
        The sample type.

    chp_path : str
        The path to the current processing chapter.

    gradient_color_by : str
        The taxon to color by for the gradient PCoA.

    barchart_categories : iterable
        The relevant barchart categories to use.

    Returns
    -------
    dict
        The paths and desired options.
    """
    opts = {k: agu.get_path(v) for k, v in agenv.paths.items()}
    opts['sample_type'] = sample_type
    opts['gradient_color_by'] = gradient_color_by
    opts['chp-path'] = chp_path

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

        cat_path = opts['ag-100nt-1k-%s-%s-biom' % (sample_type, cat)]
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

    ids : iterable
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
    return os.path.join(opts['per-sample-results'], id_)


def _base_barcode(id_):
    """Get the actual barcode from a sample ID"""
    # NOTE: old ebi accessions and processing do not follow this format.
    # expectation is <study_id>.<barcode>

    return id_.split('.')[1]


def partition_samples_by_bodysite(md, site_to_function):
    """Yield the processing function and associated sample IDs

    Parameters
    ----------
    md : DataFrame
        The metadata to examine. It is expected that this DataFrame contain a
        column named SIMPLE_BODY_SITE.
    site_to_function : list of tuple
        A mapping of a body site category value to a processing function.

    Returns
    -------
    generator
        (fucntion, iterable)
        The processing function and the sample IDs to examine.
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
    sample_ids : iterable
        A list of sample IDs of interest

    Returns
    -------
    dict
        A dict containing each sample ID and any errors observed or None if
        no error was observed for the sample.
    """
    results = {}

    for id_ in sample_ids:
        cmd = cmd_fmt % {'result_path': _result_path(opts, id_),
                         'id': id_}
        stdout, stderr, return_value = qiime_system_call(cmd)

        if return_value != 0:
            msg = stderr.splitlines()[-1]
            results[id_] = 'FAILED (%s): %s' % (msg if msg else '', cmd)
        else:
            results[id_] = None

    return results


def taxa_summaries(opts, sample_ids):
    """Produce digestable taxonomy summaries per sample

    Paramters
    ---------
    opts : dict
        A dict of relevant opts.

    sample_ids : iterable
        A list of sample IDs of interest

    Returns
    -------
    dict
        A dict containing each sample ID and any errors observed or None if
        no error was observed for the sample.
    """
    results = {}
    table_path = opts['ag-L6-taxa-%s-biom' % opts['sample_type']]
    site_table = biom.load_table(table_path)
    table_taxon_ids = site_table.ids(axis='observation')

    for id_ in sample_ids:
        if not site_table.exists(id_):
            results[id_] = 'ID not found'
        else:
            results[id_] = None
            taxa_path = os.path.join(_result_path(opts, id_), '%s.txt')

            with open(taxa_path % _base_barcode(id_), 'w') as fp:
                fp.write("#taxon\trelative_abundance\n")

                v = site_table.data(id_, dense=True)
                for sorted_v, taxa in sorted(zip(v, table_taxon_ids))[::-1]:
                    if sorted_v:
                        fp.write("%s\t%f\n" % (taxa, sorted_v))
    return results


def taxon_significance(opts, sample_ids):
    """Produce OTU significance results

    Paramters
    ---------
    opts : dict
        A dict of relevant opts.
    sample_ids : iterable
        A list of sample IDs of interest

    Returns
    -------
    dict
        A dict containing each sample ID and any errors observed or None if
        no error was observed for the sample.
    """
    table_path = opts['ag-L6-taxa-%s-biom' % opts['sample_type']]

    cmd_fmt = 'generate_otu_signifigance_tables_AGP.py -i %s ' % table_path
    cmd_fmt += '-o %(result_path)s -s %(id)s'

    return _iter_ids_over_system_call(cmd_fmt, sample_ids, opts)


def body_site_pcoa(opts, sample_ids):
    """Produce the per-sample all body site PCoA

    Paramters
    ---------
    opts : dict
        A dict of relevant opts.
    sample_ids : iterable
        A list of sample IDs of interest

    Returns
    -------
    dict
        A dict containing each sample ID and any errors observed or None if
        no error was observed for the sample.
    """
    coords = opts['ag-pgp-hmp-gg-100nt-1k-unifrac-pc']
    cmd_fmt = ' '.join(["mod2_pcoa.py body_site",
                        "--coords %s" % coords,
                        "--mapping_file %s" % opts['ag-cleaned-md'],
                        "--output %(result_path)s",
                        "--prefix Figure_1",
                        "--samples %(id)s"])

    return _iter_ids_over_system_call(cmd_fmt, sample_ids, opts)


def country_pcoa(opts, sample_ids):
    """Produce the per-sample subsampled country PCoA

    Paramters
    ---------
    opts : dict
        A dict of relevant opts.
    sample_ids : iterable
        A list of sample IDs of interest

    Returns
    -------
    dict
        A dict containing each sample ID and any errors observed or None if
        no error was observed for the sample.
    """
    coords = opts['ag-gg-100nt-1k-subsampled-unifrac-pc']
    cmd_fmt = ' '.join(["mod2_pcoa.py country",
                        "--distmat %s" % opts['ag-gg-100nt-1k-bdiv-unifrac'],
                        "--coords %s" % coords,
                        "--mapping_file %s" % opts['ag-gg-cleaned-md'],
                        "--output %(result_path)s",
                        "--prefix Figure_2",
                        "--samples %(id)s"])

    return _iter_ids_over_system_call(cmd_fmt, sample_ids, opts)


def gradient_pcoa(opts, sample_ids):
    """Produce a gradient PCoA

    Paramters
    ---------
    opts : dict
        A dict of relevant opts.
    sample_ids : iterable
        A list of sample IDs of interest

    Returns
    -------
    dict
        A dict containing each sample ID and any errors observed or None if
        no error was observed for the sample.
    """
    coords = opts['ag-100nt-%s-1k-unifrac-pc' % opts['sample_type']]
    cmd_fmt = ' '.join(["mod2_pcoa.py gradient",
                        "--coords %s" % coords,
                        "--mapping_file %s" % opts['ag-L2-taxa-md'],
                        "--output %(result_path)s",
                        "--prefix Figure_3",
                        "--color %s" % opts['gradient-color-by'],
                        "--samples %(id)s"])

    return _iter_ids_over_system_call(cmd_fmt, sample_ids, opts)


def pie_plot(opts, sample_ids):
    """Produce a pie chart

    Paramters
    ---------
    opts : dict
        A dict of relevant opts.
    sample_ids : iterable
        A list of sample IDs of interest

    Returns
    -------
    dict
        A dict containing each sample ID and any errors observed or None if
        no error was observed for the sample.
    """
    cmd_fmt = ' '.join(['make_pie_plot_AGP.py',
                        '-i %s' % opts['ag-L3-taxa-tsv'],
                        '-o %(result_path)s',
                        '-s %(id)s'])
    return _iter_ids_over_system_call(cmd_fmt, sample_ids, opts)


def bar_chart(opts, sample_ids):
    """Produce a bar chart

    Paramters
    ---------
    opts : dict
        A dict of relevant opts.
    sample_ids : iterable
        A list of sample IDs of interest

    Returns
    -------
    dict
        A dict containing each sample ID and any errors observed or None if
        no error was observed for the sample.
    """
    table_path = opts['ag-100nt-1k-%s-biom' % opts['sample_type']]
    cmd_fmt = ' '.join(['make_phyla_plots_AGP.py',
                        '-i %s' % table_path,
                        '-m %s' % opts['ag-cleaned-md'],
                        '-o %(result_path)s',
                        '-c %s' % opts['barchart_categories'],
                        '-t %s' % opts['sample_type'],
                        '-s %(id)s'])
    return _iter_ids_over_system_call(cmd_fmt, sample_ids, opts)


def per_sample_directory(opts, sample_ids):
    """Create a per sample directory

    Paramters
    ---------
    opts : dict
        A dict of relevant opts.
    sample_ids : iterable
        A list of sample IDs of interest
    """
    for id_ in sample_ids:
        path = _result_path(opts, id_)
        if not os.path.exists(path):
            os.mkdir(path)
    return {}


def stage_per_sample_specific_statics(opts, sample_ids):
    """Items like the Latex templates

    Paramters
    ---------
    opts : dict
        A dict of relevant opts.
    sample_ids : iterable
        A list of sample IDs of interest

    Returns
    -------
    dict
        A dict containing each sample ID and any errors observed or None if
        no error was observed for the sample.
    """
    chp_path = opts['chp-path']
    for id_ in sample_ids:
        path = _result_path(opts, id_)

        if opts['sample_type'] == 'fecal':
            shutil.copy(os.path.join(chp_path, 'template_gut.tex'), path)
        elif opts['sample_type'] in ('oral', 'skin'):
            shutil.copy(os.path.join(chp_path, 'template_oralskin.tex'), path)
        else:
            raise ValueError('Unknown sample type: %s' % opts['sample_type'])
    return {}
