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
    'raw-sequences': '1/raw-sequences.fna',
    'raw-metadata': '1/raw-metadata.txt',

    # sequences filtered for blooms
    'filtered-sequences': '2/filtered-sequences.fna',
    'filtered-sequences-100nt': '2/filtered-sequences-100nt.fna',

    # only fecal sequences (for filtering for blooms)
    'fecal-sequences': '2/fecal-sequences.fna',

    # observed bloom sequences in samples
    'observed-blooms': '2/observed-blooms',
    'observed-blooms-biom': '2/observed-blooms/otu_table.biom',
    'observed-blooms-otu-map':
        '2/observed-blooms/sortmerna_picked_otus/fecal-sequences_otus.txt',

    # resulting OTU directories
    'ag-otus': '3/otus/gg-13_8-97-per-otus',
    'ag-otus-100nt': '3/otus/gg-13_8-97-per-otus-with-100nt',
    'ag-biom': '3/otus/gg-13_8-97-per-otus/otu_table.biom',
    'ag-100nt-biom': '3/otus/gg-13_8-97-per-otus-with-100nt/otu_table.biom',

    # merged files for diversity analyses
    'ag-gg-100nt-biom': '4/ag-gg-100nt.biom',
    'pgp-hmp-100nt-biom': '4/pgp-hmp-100nt.biom',
    'ag-pgp-hmp-gg-100nt-biom': '4/ag-pgp-hmp-gg-100nt.biom',
    'ag-cleaned-md': '4/ag-cleaned.txt',
    'gg-cleaned-md': '4/gg-cleaned.txt',
    'pgp-cleaned-md': '4/pgp-cleaned.txt',
    'hmp-cleaned-md': '4/hmp-cleaned.txt',
    'ag-gg-cleaned-md': '4/ag-gg-cleaned.txt',
    'pgp-hmp-cleaned-md': '4/pgp-hmp-cleaned.txt',
    'ag-pgp-hmp-gg-cleaned-md': '4/ag-pgp-hmp-gg-cleaned.txt',

    # alpha diversity analysis files
    'ag-pgp-hmp-gg-100nt-1k-multiple': '5/ag-pgp-hmp-gg-100nt-1k-multiple',


    'ag-pgp-hmp-gg-100nt-1k-adiv': '5/ag-pgp-hmp-gg-100nt-1k-adiv',
    'ag-pgp-hmp-gg-100nt-1k-adiv-pd':
        '5/ag-pgp-hmp-gg-100nt-1k-adiv/PD_whole_tree.txt',
    'ag-pgp-hmp-gg-100nt-1k-adiv-chao1':
        '5/ag-pgp-hmp-gg-100nt-1k-adiv/chao1.txt',
    'ag-pgp-hmp-gg-100nt-1k-adiv-observedotus':
        '5/ag-pgp-hmp-gg-100nt-1k-adiv/observed_otus.txt',

    # beta diversity analysis files
    'ag-pgp-hmp-gg-100nt-1k-biom': '6/ag-pgp-hmp-gg-100nt-1k.biom',
    'ag-pgp-hmp-gg-100nt-1k-bdiv': '6/ag-pgp-hmp-gg-100nt-1k-bdiv',

    'ag-100nt-1k-bdiv-unifrac':
        ('6/ag-pgp-hmp-gg-100nt-1k-bdiv/'
         'unweighted_unifrac_ag-100nt-1k.txt'),
    'ag-100nt-1k-unifrac-pc':
        ('6/ag-pgp-hmp-gg-100nt-1k-bdiv/'
         'unweighted_unifrac_ag-100nt-1k-pc.txt'),

    'ag-100nt-oral-1k-bdiv-unifrac':
        ('6/ag-pgp-hmp-gg-100nt-1k-bdiv/'
         'unweighted_unifrac_ag-100nt-oral-1k.txt'),
    'ag-100nt-oral-1k-unifrac-pc':
        ('6/ag-pgp-hmp-gg-100nt-1k-bdiv/'
         'unweighted_unifrac_ag-100nt-oral-1k-pc.txt'),

    'ag-100nt-fecal-1k-bdiv-unifrac':
        ('6/ag-pgp-hmp-gg-100nt-1k-bdiv/'
         'unweighted_unifrac_ag-100nt-fecal-1k.txt'),
    'ag-100nt-fecal-1k-unifrac-pc':
        ('6/ag-pgp-hmp-gg-100nt-1k-bdiv/'
         'unweighted_unifrac_ag-100nt-fecal-1k-pc.txt'),

    'ag-100nt-skin-1k-bdiv-unifrac':
        ('6/ag-pgp-hmp-gg-100nt-1k-bdiv/'
         'unweighted_unifrac_ag-100nt-skin-1k.txt'),
    'ag-100nt-skin-1k-unifrac-pc':
        ('6/ag-pgp-hmp-gg-100nt-1k-bdiv/'
         'unweighted_unifrac_ag-100nt-skin-1k-pc.txt'),

    'ag-pgp-hmp-gg-100nt-1k-bdiv-unifrac':
        ('6/ag-pgp-hmp-gg-100nt-1k-bdiv/'
         'unweighted_unifrac_ag-pgp-hmp-gg-100nt-1k.txt'),
    'ag-pgp-hmp-gg-100nt-1k-unifrac-pc':
        ('6/ag-pgp-hmp-gg-100nt-1k-bdiv/'
         'unweighted_unifrac_ag-pgp-hmp-gg-100nt-1k-pc.txt'),

    'ag-pgp-hmp-gg-100nt-1k-bdiv-wunifrac':
        ('6/ag-pgp-hmp-gg-100nt-1k-bdiv/'
         'weighted_unifrac_ag-pgp-hmp-gg-100nt-1k.txt'),
    'ag-pgp-hmp-gg-100nt-1k-wunifrac-pc':
        ('6/ag-pgp-hmp-gg-100nt-1k-bdiv/'
         'weighted_unifrac_ag-pgp-hmp-gg-100nt-1k-pc.txt'),

    'ag-gg-100nt-1k-bdiv-unifrac':
        ('6/ag-pgp-hmp-gg-100nt-1k-bdiv/'
         'unweighted_unifrac_ag-gg-100nt-1k.txt'),
    'ag-gg-100nt-1k-unifrac-pc':
        ('6/ag-pgp-hmp-gg-100nt-1k-bdiv/'
         'unweighted_unifrac_ag-gg-100nt-1k-pc.txt'),

    'ag-gg-100nt-1k-bdiv-subsampled-unifrac':
        ('6/ag-pgp-hmp-gg-100nt-1k-bdiv/'
         'unweighted_unifrac_ag-gg-100nt-1k-subsampled.txt'),
    'ag-gg-100nt-1k-subsampled-unifrac-pc':
        ('6/ag-pgp-hmp-gg-100nt-1k-bdiv/'
         'unweighted_unifrac_ag-gg-100nt-1k-subsampled-pc.txt'),

    'ag-gg-100nt-1k-bdiv-wunifrac':
        ('6/ag-pgp-hmp-gg-100nt-1k-bdiv/'
         'weighted_unifrac_ag-gg-100nt-1k.txt'),
    'ag-gg-100nt-1k-wunifrac-pc':
        ('6/ag-pgp-hmp-gg-100nt-1k-bdiv/'
         'weighted_unifrac_ag-gg-100nt-1k-pc.txt'),

    # taxonomy summaries
    'ag-taxa': '7/taxa',
    'ag-L2-taxa-tsv': '7/taxa/otu_table_L2.txt',
    'ag-L2-taxa-biom': '7/taxa/otu_table_L2.biom',
    'ag-L2-taxa-md': '7/taxa/ag-cleaned_L2.txt',
    'ag-L3-taxa-tsv': '7/taxa/otu_table_L3.txt',
    'ag-L3-taxa-biom': '7/taxa/otu_table_L3.biom',
    'ag-L6-taxa-tsv': '7/taxa/otu_table_L6.txt',
    'ag-L6-taxa-biom': '7/taxa/otu_table_L6.biom',
    'ag-L2-taxa-skin-biom': '7/taxa/otu_table_skin_L2.biom',
    'ag-L2-taxa-oral-biom': '7/taxa/otu_table_oral_L2.biom',
    'ag-L2-taxa-fecal-biom': '7/taxa/otu_table_fecal_L2.biom',
    'ag-L6-taxa-skin-biom': '7/taxa/otu_table_skin_L6.biom',
    'ag-L6-taxa-oral-biom': '7/taxa/otu_table_oral_L6.biom',
    'ag-L6-taxa-fecal-biom': '7/taxa/otu_table_fecal_L6.biom',

    # collapsed samples
    'ag-100nt-1k-biom': '8/ag-100nt-1k.biom',
    'ag-100nt-1k-fecal-biom': '8/ag-100nt-1k-fecal.biom',
    'ag-100nt-1k-skin-biom': '8/ag-100nt-1k-oral.biom',
    'ag-100nt-1k-oral-biom': '8/ag-100nt-1k-skin.biom',

    'ag-100nt-1k-fecal-sex-biom':  '8/ag-100nt-1k-fecal-sex.biom',
    'ag-100nt-1k-fecal-diet-biom': '8/ag-100nt-1k-fecal-diet.biom',
    'ag-100nt-1k-fecal-age-biom':  '8/ag-100nt-1k-fecal-age.biom',
    'ag-100nt-1k-fecal-bmi-biom':  '8/ag-100nt-1k-fecal-bmi.biom',

    'ag-100nt-1k-oral-sex-biom':  '8/ag-100nt-1k-oral-sex.biom',
    'ag-100nt-1k-oral-diet-biom': '8/ag-100nt-1k-oral-diet.biom',
    'ag-100nt-1k-oral-age-biom':  '8/ag-100nt-1k-oral-age.biom',
    'ag-100nt-1k-oral-flossing-biom':  '8/ag-100nt-1k-oral-flossing.biom',

    'ag-100nt-1k-skin-sex-biom':  '8/ag-100nt-1k-skin-sex.biom',
    'ag-100nt-1k-skin-cosmetics-biom': '8/ag-100nt-1k-skin-cosmetics.biom',
    'ag-100nt-1k-skin-age-biom':  '8/ag-100nt-1k-skin-age.biom',
    'ag-100nt-1k-skin-hand-biom':  '8/ag-100nt-1k-skin-hand.biom',

    # per-sample results
    'successful-ids': '9/successful_ids.txt',
    'unsuccessful-ids': '9/unsuccessful_ids.txt',
    'per-sample-results': '9/per-sample-results',
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
    if ag.is_test_env():
        return "100"
    else:
        return "1000"


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
        dst = os.path.join(ag.WORKING_DIR, '1/%s' % acc)
        shutil.copytree(src, dst)
