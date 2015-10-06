# establish any minimals for the notebook environment
import os
import shutil
from distutils.spawn import find_executable

import qiime

import americangut as ag
from americangut.results_utils import get_repository_dir
from americangut.util import get_path, get_new_path, get_existing_path


_EBI_ACCESSIONS = ['ERP012511']
_TEST_ACCESSIONS = ['ag_testing']


# essential filenames that may be used between notebooks
filenames = {
    'raw_sequences': 'agp_sequences.fna',
    'raw_metadata' : 'agp_metadata.txt'
}


def _assert_environment():
    if qiime.__version__ != '1.9.1':
        raise ImportError("QIIME 1.9.1 is not in the environment.")

    if find_executable('print_qiime_config.py') is None:
        raise EnvironmentError("The QIIME executables are not in $PATH.")

    if find_executable('mod2_pcoa.py') is None:
        raise EnvironmentError("The AG scripts are not in $PATH.")
_assert_environment()


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
