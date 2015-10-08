# establish any minimals for the notebook environment
import os
import shutil
from distutils.spawn import find_executable

import qiime
import qiime_default_reference as qdr

import americangut as ag
from americangut.results_utils import get_repository_dir
from americangut.util import get_existing_path


_EBI_ACCESSIONS = ['ERP012511']
_TEST_ACCESSIONS = ['ag_testing']


# essential paths relative to the working_dir to be used between notebooks

paths = {
    # raw files
    'raw-sequences': 'raw-sequences.fna',
    'raw-metadata': 'raw_metadata.txt',

    # sequences filtered for blooms
    'filtered-sequences': 'filtered-sequences.fna',
    'filtered-sequences-100nt': 'filtered-sequences-100nt.fna',

    # only fecal sequences (for filtering for blooms)
    'fecal-sequences': 'fecal-sequences.fna',

    # observed bloom sequences in samples
    'observed-blooms': 'observed-blooms',
    'observed-blooms-biom': 'observed-blooms/otu_table.biom',
    'observed-blooms-otu-map':
        'observed-blooms/sortmerna_picked_otus/fecal-sequences_otus.txt',

    # resulting OTU directories
    'gg-otus': 'otus/gg-13_8-97-percent-otus',
    'gg-otus-100nt': 'otus/gg-13_8-97-percent-otus-with-100nt',
    'gg-otus-biom': 'otus/gg-13_8-97-percent-otus/otu_table.biom',
    'gg-otus-100nt-biom':
        'otus/gg-13_8-97-percent-otus-with-100nt/otu_table.biom',
}


def _assert_environment():
    if qiime.__version__ != '1.9.1':
        raise ImportError("QIIME 1.9.1 is not in the environment.")

    if find_executable('print_qiime_config.py') is None:
        raise EnvironmentError("The QIIME executables are not in $PATH.")

    if find_executable('mod2_pcoa.py') is None:
        raise EnvironmentError("The AG scripts are not in $PATH.")
_assert_environment()


def get_sortmerna_index():
    """Return the absolute path a SortMeRNA index if available"""
    return os.environ.get('AG_SMR_INDEX')


def get_reference_set():
    """Get the reference set to use for OTU picking

    Returns
    -------
    (str, str)
        The file paths to the reference sequences and the reference taxonomy.
    """
    if os.environ.get('AG_TESTING') == 'True':
        repo = get_repository_dir()
        ref_seqs = os.path.join(repo, 'tests/data/otus.fna')
        ref_tax = os.path.join(repo, 'tests/data/otus.txt')
        return ref_seqs, ref_tax
    else:
        return qdr.get_reference_sequences(), qdr.get_reference_taxonomy()


def get_accessions():
    """Get the accessions to use, or redirect to test data

    Notes
    -----
    If os.environ['AG_TESTING'] == 'True', then the accessions returned will
    correspond to the test dataset.

    Returns
    -------
    list of str
        The accessions, which are expected to be basenames for the actual data.
        For instance, the accession "foo" would have sequences as "foo.fna" and
        metadata as "foo.txt".
    """
    if os.environ.get('AG_TESTING') == 'True':
        _stage_test_accessions()
        return _TEST_ACCESSIONS[:]
    else:
        return _EBI_ACCESSIONS[:]


def get_bloom_sequences():
    """Get the filepath to the bloom sequences

    Raises
    ------
    IOError
        If the path does not exist

    Returns
    -------
    str
        The filepath to the bloom sequences
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
        src_fna = os.path.join(repo, 'tests/data/%s.fna' % acc)
        src_map = os.path.join(repo, 'tests/data/%s.txt' % acc)

        shutil.copy(src_fna, ag.working_dir)
        shutil.copy(src_map, ag.working_dir)
