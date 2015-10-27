# establish any minimals for the notebook environment
import os
import shutil
from distutils.spawn import find_executable

import qiime
import qiime_default_reference as qdr

import americangut as ag
from americangut.results_utils import get_repository_dir
from americangut.util import get_existing_path


_EBI_ACCESSIONS = ['ERP012803']
_TEST_ACCESSIONS = ['ag_testing']


# essential paths relative to the working_dir to be used between notebooks

paths = {
    # raw files
    'raw': {
        'sequences': '01-raw/sequences.fna',
        'metadata': '01-raw/metadata.txt',
        },

    # sequences filtered for blooms
    'filtered': {
        'sequences': '02-filtered/sequences.fna',
        'sequences-100nt': '02-filtered/sequences-100nt.fna',

        # only fecal sequences (for filtering for blooms)
        'fecal-sequences': '02-filtered/fecal-sequences.fna',

        # observed bloom sequences in samples
        'observed-blooms': '02-filtered/observed-blooms',
        'observed-blooms-biom': '02-filtered/observed-blooms/otu_table.biom',
        'observed-blooms-otu-map':
            '02-filtered/observed-blooms/sortmerna_picked_otus/fecal-sequences_otus.txt',
        },

    # resulting OTU data
    'otus': {
        'ag': '03-otus/gg-13_8-97-percent',
        'ag-100nt': '03-otus/gg-13_8-97-percent-100nt',
        'ag-biom': '03-otus/gg-13_8-97-percent/otu_table.biom',
        'ag-100nt-biom': '03-otus/gg-13_8-97-percent-100nt/otu_table.biom',
        },

    # merged files for diversity analyses
    'meta': {
        'ag-gg-100nt-biom': '04-meta/ag-gg-100nt.biom',
        'pgp-hmp-100nt-biom': '04-meta/pgp-hmp-100nt.biom',
        'ag-pgp-hmp-gg-100nt-biom': '04-meta/ag-pgp-hmp-gg-100nt.biom',
        'ag-cleaned-md': '04-meta/ag-cleaned.txt',
        'gg-cleaned-md': '04-meta/gg-cleaned.txt',
        'pgp-cleaned-md': '04-meta/pgp-cleaned.txt',
        'hmp-cleaned-md': '04-meta/hmp-cleaned.txt',
        'ag-gg-cleaned-md': '04-meta/ag-gg-cleaned.txt',
        'pgp-hmp-cleaned-md': '04-meta/pgp-hmp-cleaned.txt',
        'ag-pgp-hmp-gg-cleaned-md': '04-meta/ag-pgp-hmp-gg-cleaned.txt',
        },

    # alpha diversity analysis files
    'alpha': {
        '1k': {
            'ag-pgp-hmp-gg-100nt-multiple': '05-alpha/1k/ag-pgp-hmp-gg-100nt-multiple',
            'ag-pgp-hmp-gg-100nt': '05-alpha/1k/ag-pgp-hmp-gg-100nt',
            'ag-pgp-hmp-gg-100nt-pd':
                '05-alpha/1k/ag-pgp-hmp-gg-100nt/PD_whole_tree.txt',
            'ag-pgp-hmp-gg-100nt-chao1':
                '05-alpha/1k/ag-pgp-hmp-gg-100nt/chao1.txt',
            'ag-pgp-hmp-gg-100nt-observedotus':
                '05-alpha/1k/ag-pgp-hmp-gg-100nt/observed_otus.txt',
            },
        '10k': {
            'ag-pgp-hmp-gg-100nt-multiple': '05-alpha/10k/ag-pgp-hmp-gg-100nt-multiple',
            'ag-pgp-hmp-gg-100nt': '05-alpha/10k/ag-pgp-hmp-gg-100nt',
            'ag-pgp-hmp-gg-100nt-pd':
                '05-alpha/10k/ag-pgp-hmp-gg-100nt/PD_whole_tree.txt',
            'ag-pgp-hmp-gg-100nt-chao1':
                '05-alpha/10k/ag-pgp-hmp-gg-100nt/chao1.txt',
            'ag-pgp-hmp-gg-100nt-observedotus':
                '05-alpha/10k/ag-pgp-hmp-gg-100nt/observed_otus.txt',
            }
        },

    # beta diversity analysis files
    'beta': {
        '1k': {
            'ag-pgp-hmp-gg-100nt-biom': '06-beta/1k/ag-pgp-hmp-gg-100nt.biom',
            'ag-pgp-hmp-gg-100nt': '06-beta/1k/ag-pgp-hmp-gg-100nt',

            'ag-100nt-unifrac':
                ('06-beta/1k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-100nt.txt'),
            'ag-100nt-unifrac-pc':
                ('06-beta/1k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-100nt-pc.txt'),

            'ag-100nt-oral-unifrac':
                ('06-beta/1k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-100nt-oral.txt'),
            'ag-100nt-oral-unifrac-pc':
                ('06-beta/1k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-100nt-oral-pc.txt'),

            'ag-100nt-fecal-unifrac':
                ('06-beta/1k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-100nt-fecal.txt'),
            'ag-100nt-fecal-unifrac-pc':
                ('06-beta/1k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-100nt-fecal-pc.txt'),

            'ag-100nt-skin-unifrac':
                ('06-beta/1k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-100nt-skin.txt'),
            'ag-100nt-skin-unifrac-pc':
                ('06-beta/1k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-100nt-skin-pc.txt'),

            'ag-pgp-hmp-gg-100nt-unifrac':
                ('06-beta/1k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-pgp-hmp-gg-100nt.txt'),
            'ag-pgp-hmp-gg-100nt-unifrac-pc':
                ('06-beta/1k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-pgp-hmp-gg-100nt-pc.txt'),

            'ag-pgp-hmp-gg-100nt-wunifrac':
                ('06-beta/1k/ag-pgp-hmp-gg-100nt/'
                 'weighted_unifrac_ag-pgp-hmp-gg-100nt.txt'),
            'ag-pgp-hmp-gg-100nt-wunifrac-pc':
                ('06-beta/1k/ag-pgp-hmp-gg-100nt/'
                 'weighted_unifrac_ag-pgp-hmp-gg-100nt-pc.txt'),

            'ag-gg-100nt-unifrac':
                ('06-beta/1k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-gg-100nt.txt'),
            'ag-gg-100nt-unifrac-pc':
                ('06-beta/1k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-gg-100nt-pc.txt'),

            'ag-gg-100nt-subsampled-unifrac':
                ('06-beta/1k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-gg-100nt-subsampled.txt'),
            'ag-gg-100nt-subsampled-unifrac-pc':
                ('06-beta/1k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-gg-100nt-subsampled-pc.txt'),

            'ag-gg-100nt-wunifrac':
                ('06-beta/1k/ag-pgp-hmp-gg-100nt/'
                 'weighted_unifrac_ag-gg-100nt.txt'),
            'ag-gg-100nt-wunifrac-pc':
                ('06-beta/1k/ag-pgp-hmp-gg-100nt/'
                 'weighted_unifrac_ag-gg-100nt-pc.txt'),
            },
        '10k': {
            'ag-pgp-hmp-gg-100nt-biom': '06-beta/10k/ag-pgp-hmp-gg-100nt.biom',
            'ag-pgp-hmp-gg-100nt': '06-beta/10k/ag-pgp-hmp-gg-100nt',

            'ag-100nt-unifrac':
                ('06-beta/10k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-100nt.txt'),
            'ag-100nt-unifrac-pc':
                ('06-beta/10k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-100nt-pc.txt'),

            'ag-100nt-oral-unifrac':
                ('06-beta/10k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-100nt-oral.txt'),
            'ag-100nt-oral-unifrac-pc':
                ('06-beta/10k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-100nt-oral-pc.txt'),

            'ag-100nt-fecal-unifrac':
                ('06-beta/10k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-100nt-fecal.txt'),
            'ag-100nt-fecal-unifrac-pc':
                ('06-beta/10k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-100nt-fecal-pc.txt'),

            'ag-100nt-skin-unifrac':
                ('06-beta/10k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-100nt-skin.txt'),
            'ag-100nt-skin-unifrac-pc':
                ('06-beta/10k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-100nt-skin-pc.txt'),

            'ag-pgp-hmp-gg-100nt-unifrac':
                ('06-beta/10k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-pgp-hmp-gg-100nt.txt'),
            'ag-pgp-hmp-gg-100nt-unifrac-pc':
                ('06-beta/10k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-pgp-hmp-gg-100nt-pc.txt'),

            'ag-pgp-hmp-gg-100nt-wunifrac':
                ('06-beta/10k/ag-pgp-hmp-gg-100nt/'
                 'weighted_unifrac_ag-pgp-hmp-gg-100nt.txt'),
            'ag-pgp-hmp-gg-100nt-wunifrac-pc':
                ('06-beta/10k/ag-pgp-hmp-gg-100nt/'
                 'weighted_unifrac_ag-pgp-hmp-gg-100nt-pc.txt'),

            'ag-gg-100nt-unifrac':
                ('06-beta/10k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-gg-100nt.txt'),
            'ag-gg-100nt-unifrac-pc':
                ('06-beta/10k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-gg-100nt-pc.txt'),

            'ag-gg-100nt-subsampled-unifrac':
                ('06-beta/10k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-gg-100nt-subsampled.txt'),
            'ag-gg-100nt-subsampled-unifrac-pc':
                ('06-beta/10k/ag-pgp-hmp-gg-100nt/'
                 'unweighted_unifrac_ag-gg-100nt-subsampled-pc.txt'),

            'ag-gg-100nt-wunifrac':
                ('06-beta/10k/ag-pgp-hmp-gg-100nt/'
                 'weighted_unifrac_ag-gg-100nt.txt'),
            'ag-gg-100nt-wunifrac-pc':
                ('06-beta/10k/ag-pgp-hmp-gg-100nt/'
                 'weighted_unifrac_ag-gg-100nt-pc.txt'),
            }
        },

    # taxonomy summaries
    'taxa': {
        'L2': {
            'ag-tsv': '07-taxa/otu_table_L2.txt',
            'ag-biom': '07-taxa/otu_table_L2.biom',
            'ag-md': '07-taxa/ag-cleaned_L2.txt',
            'ag-skin-biom': '07-taxa/otu_table_skin_L2.biom',
            'ag-oral-biom': '07-taxa/otu_table_oral_L2.biom',
            'ag-fecal-biom': '07-taxa/otu_table_fecal_L2.biom',
            },
        'L3': {
            'ag-tsv': '07-taxa/otu_table_L3.txt',
            'ag-biom': '07-taxa/otu_table_L3.biom',
            'ag-md': '07-taxa/ag-cleaned_L3.txt',
            'ag-skin-biom': '07-taxa/otu_table_skin_L3.biom',
            'ag-oral-biom': '07-taxa/otu_table_oral_L3.biom',
            'ag-fecal-biom': '07-taxa/otu_table_fecal_L3.biom',
            },
        'L6': {
            'ag-tsv': '07-taxa/otu_table_L6.txt',
            'ag-biom': '07-taxa/otu_table_L6.biom',
            'ag-md': '07-taxa/ag-cleaned_L6.txt',
            'ag-skin-biom': '07-taxa/otu_table_skin_L6.biom',
            'ag-oral-biom': '07-taxa/otu_table_oral_L6.biom',
            'ag-fecal-biom': '07-taxa/otu_table_fecal_L6.biom',
            }
        },

    # collapsed samples
    'collapsed': {
        'ag-100nt-1k-biom': '08-collapsed/ag-100nt-1k.biom',
        'ag-100nt-1k-fecal-biom': '08-collapsed/ag-100nt-1k-fecal.biom',
        'ag-100nt-1k-skin-biom': '08-collapsed/ag-100nt-1k-oral.biom',
        'ag-100nt-1k-oral-biom': '08-collapsed/ag-100nt-1k-skin.biom',

        'ag-100nt-1k-fecal-sex-biom':  '08-collapsed/ag-100nt-1k-fecal-sex.biom',
        'ag-100nt-1k-fecal-diet-biom': '08-collapsed/ag-100nt-1k-fecal-diet.biom',
        'ag-100nt-1k-fecal-age-biom':  '08-collapsed/ag-100nt-1k-fecal-age.biom',
        'ag-100nt-1k-fecal-bmi-biom':  '08-collapsed/ag-100nt-1k-fecal-bmi.biom',

        'ag-100nt-1k-oral-sex-biom':  '08-collapsed/ag-100nt-1k-oral-sex.biom',
        'ag-100nt-1k-oral-diet-biom': '08-collapsed/ag-100nt-1k-oral-diet.biom',
        'ag-100nt-1k-oral-age-biom':  '08-collapsed/ag-100nt-1k-oral-age.biom',
        'ag-100nt-1k-oral-flossing-biom':  '08-collapsed/ag-100nt-1k-oral-flossing.biom',

        'ag-100nt-1k-skin-sex-biom':  '08-collapsed/ag-100nt-1k-skin-sex.biom',
        'ag-100nt-1k-skin-cosmetics-biom': '08-collapsed/ag-100nt-1k-skin-cosmetics.biom',
        'ag-100nt-1k-skin-age-biom':  '08-collapsed/ag-100nt-1k-skin-age.biom',
        'ag-100nt-1k-skin-hand-biom':  '08-collapsed/ag-100nt-1k-skin-hand.biom',
    },

    # per-sample results
    'per-sample': {
        'successful-ids': '09-per-sample/successful_ids.txt',
        'unsuccessful-ids': '09-per-sample/unsuccessful_ids.txt',
        'per-sample-results': '09-per-sample/results',
        'statics-fecal': '09/statics-fecal',
        'statics-oral': '09/statics-oral',
        'statics-skin': '09/statics-skin',
    },

    'populated-results': {
        'result-pdfs': '10/result_pdfs/',
        'result-taxa': '10/result_taxa/',
        'successful-pdfs': '10/successful_pdfs.txt',
        'unsuccessful-pdfs': '10/unsuccessful_pdfs.txt',
    },

    'demux': {},
    'otu-map': {},
}


def _assert_environment():
    if qiime.__version__ != '1.9.1':
        obs_version = qiime.__version__
        raise ImportError("QIIME 1.9.1 is not in the environment, found "
                          "%s." % obs_version)

    if find_executable('print_qiime_config.py') is None:
        raise EnvironmentError("The QIIME executables are not in $PATH.")

    if find_executable('mod2_pcoa.py') is None:
        raise EnvironmentError("The AG scripts are not in $PATH.")
_assert_environment()


def activate(chp):
    """Activate a chapter

    Parameters
    ----------
    chp : str
        The chapter.

    Returns
    -------
    str
        The path to the activated directory

    Notes
    -----
    Activation creates the chapter processing directory if it does
    not already exist.
    """
    path = os.path.join(ag.WORKING_DIR, chp)
    if not os.path.exists(path):
        os.mkdir(path)
    return path


def get_sortmerna_index():
    """Return the absolute path of a SortMeRNA index if available"""
    return os.environ.get('AG_SMR_INDEX')


def get_rarefaction_depth():
    """Return the rarefaction depth to use"""
<<<<<<< HEAD
    if _TEST_ENV:
        return ("100", "200")
=======
    if ag.is_test_env():
        return "100"
>>>>>>> master
    else:
        return ("1000", "10000")


def get_reference_set():
    """Get the reference set to use for OTU picking

    Returns
    -------
    str
        The file path to the reference sequences.
    str
        The file path to the reference taxonomy.
    """
    if ag.is_test_env():
        repo = get_repository_dir()
        ref_seqs = os.path.join(repo, 'tests/data/otus.fna')
        ref_tax = os.path.join(repo, 'tests/data/otus.txt')
        return ref_seqs, ref_tax
    else:
        return qdr.get_reference_sequences(), qdr.get_reference_taxonomy()


def get_cpu_count():
    """Get the CPU count to use

    Returns
    -------
    int
        The number of CPUs available for use.

    Notes
    -----
    This method provides a layer of indirection in the event that a) a user
    does not want to allow all cores to be used or b) for testing purposes
    as Travis CI appears to be unreliable under multiprocessing and QIIME.
    """
    if os.environ.get('AG_CPU_COUNT') is not None:
        return int(os.environ.get('AG_CPU_COUNT'))
    else:
        # import is scoped as far its use is contained here.
        import multiprocessing
        return multiprocessing.cpu_count()


def get_hmp():
    """Get the HMP 100nt table and mapping"""
    return _get_data('HMP', 'HMPv35_100nt')


def get_pgp():
    """Get the PGP 100nt table and mapping"""
    return _get_data('PGP', 'PGP_100nt')


def get_global_gut():
    """Get the Global Gut table and mapping"""
    return _get_data('GG', 'GG_100nt')


def _get_data(data_dir, tag):
    """Get a non-AG table and mapping file

    Parameters
    ----------
    data_dir : str
        The base data path
    tag : str
        The filetag (e.g., HMPv35_100nt)

    Returns
    -------
    (str, str)
        The filepath to the table, and the filepath to the mapping file.

    Notes
    -----
    If $AG_TESTING == 'True', then the data returned will correspond to the
    test dataset.

    Raises
    ------
    IOError
        If the filepaths are not accessible
    """
    repo = get_repository_dir()
    data = 'tests/data' if ag.is_test_env() else 'data'
    base = os.path.join(repo, data)

    table = os.path.join(base, data_dir, '%s.biom' % tag)
    mapping = os.path.join(base, data_dir, '%s.txt' % tag)

    if not os.path.exists(table):
        raise IOError("Unable to access: %s" % table)
    if not os.path.exists(mapping):
        raise IOError("Unable to access: %s" % table)

    return table, mapping


def get_study_accessions():
    """Get the accessions to use, or redirect to test data

    Returns
    -------
    list of str
        The accessions, which are expected to be basenames for the actual data.
        For instance, the accession "foo" would have sequences as "foo.fna" and
        metadata as "foo.txt".

    Notes
    -----
    If $AG_TESTING == 'True', then the accessions returned will
    correspond to the test dataset.
    """
    if ag.is_test_env():
        _stage_test_accessions()
        return _TEST_ACCESSIONS[:]
    else:
        return _EBI_ACCESSIONS[:]


def get_files(rootdir, suffix):
    """Get the filepaths with a given suffix

    Parameters
    ----------
    rootdir : str
        The root directory to look under
    suffix : str
        The file suffix of the files to keep

    Returns
    -------
    fps : list, str
        List of file paths for all of the
        sample fasta files

    Note
    ----
    This only looks at the directory names under the
    root directory.  This assumes that the sample names
    correspond to the folders within the base folder
    """
    fps = []
    for root, dirs, files in os.walk(rootdir):
        for _file in files:
            if _file.endswith(".%s" % suffix):
                fps.append(os.path.join(root, _file))
    return fps


def get_bloom_sequences():
    """Get the filepath to the bloom sequences

    Returns
    -------
    str
        The filepath to the bloom sequences

    Raises
    ------
    IOError
        If the path does not exist
    """
    repo = get_repository_dir()
    return get_existing_path(os.path.join(repo, 'data/AG/BLOOM.fasta'))


def _stage_test_accessions():
    """Stage test data

    Notes
    -----
    Staging copies the test dataset into the working directory. This "tricks"
    the fetch_study mechanism as it'll appear that the data have already been
    sourced from EBI.
    """
    repo = get_repository_dir()
    for acc in _TEST_ACCESSIONS:
        src = os.path.join(repo, 'tests/data/%s' % acc)
        dst = os.path.join(ag.WORKING_DIR, '01/%s' % acc)

        if not os.path.exists(os.path.join(ag.WORKING_DIR, '01')):
            os.mkdir('01')

<<<<<<< HEAD
        shutil.copy(src_fna, os.path.join(ag.WORKING_DIR, '1-raw'))
        shutil.copy(src_map, os.path.join(ag.WORKING_DIR, '1-raw'))
=======
        shutil.copytree(src, dst)
>>>>>>> master
