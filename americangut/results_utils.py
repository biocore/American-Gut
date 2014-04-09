#!/usr/bin/env python

import os
import shutil
import zipfile

# These are the data files in the American-Gut repository that are used for
# results processing
_data_files = [
        ('AG', 'AG_100nt.biom.gz'),
        ('AG', 'AG_100nt.txt'),
        ('PGP', 'PGP_100nt.biom.gz'),
        ('PGP', 'PGP_100nt.txt'),
        ('HMP', 'HMPv35_100nt.biom.gz'),
        ('HMP', 'HMPv35_100nt.txt'),
        ('GG', 'GG_100nt.biom.gz'),
        ('GG', 'GG_100nt.txt')
        ]


# These are the Latex templates for the different results types
_templates = {
        'fecal': ('template_gut.tex', 'macros_gut.tex'),
        'oralskin': ('template_oralskin.tex', 'macros_oralskin.tex')
        }


def check_file(f):
    """Verify a file (or directory) exists"""
    if not os.path.exists(f):
        raise IOError("Cannot continue! The file %s does not exist!" % f)


def get_path(d, f):
    """Check and get a path, or throw IOError"""
    path = os.path.join(d, f)
    check_file(path)
    return path


def get_repository_dir():
    """Get the root of the American-Gut repository"""
    expected = os.path.abspath(__file__).rsplit('/', 2)[0]

    # get_path verifies the existance of these directories
    _ = get_path(expected, 'data')
    _ = get_path(expected, 'latex')

    return expected


def get_repository_data():
    """Get the path to the data"""
    return get_path(get_repository_dir(), 'data')


def get_repository_latex():
    """Get the path to the latex directory"""
    return get_path(get_repository_dir(), 'latex')


def get_repository_latex_pdfs(sample_type):
    """Get the Latex static PDFs directory"""
    latex_dir = get_repository_latex()

    if sample_type == 'oralskin':
        pdfs_dir = get_path(latex_dir, 'pdfs-oralskin')
    elif sample_type == 'fecal':
        pdfs_dir = get_path(latex_dir, 'pdfs-gut')
    else:
        raise ValueError("Unknown sample type: %s" % sample_type)

    check_file(pdfs_dir)

    return pdfs_dir


def _stage_static_latex(sample_type, working_dir):
    latex_dir = get_repository_latex()

    for item in _templates[sample_type]:
        src = get_path(latex_dir, item)
        shutil.copy(src, working_dir)


def _stage_static_pdfs(sample_type, working_dir):
    pdfs_dir = get_repository_latex_pdfs(sample_type)

    for f in os.listdir(pdfs_dir):
        src = get_path(pdfs_dir, f)
        shutil.copy(src, working_dir)


def _stage_static_data(working_dir):
    data_dir = get_repository_data()

    for d, f in _data_files:
        src = get_path(get_path(data_dir, d), f)
        shutil.copy(src, working_dir)


def stage_static_files(sample_type, working_dir):
    """Stage static files in the current working directory"""
    _stage_static_data(working_dir)
    _stage_static_latex(sample_type, working_dir)
    _stage_static_pdfs(sample_type, working_dir)


# use participant names only if the data are available.
# NOTE: these data are _not_ part of the github repository for
#       privacy reasons.
def parse_identifying_data(path, passwd, embedded_file='participants.txt'):
    """Process identifying data if available

    The expected format of the file is a passworded zipfile that contains
    an embedded, tab delimited file. The format of the tab delimited file
    is expected to be barcode TAB participant name
    """
    if path is not None:
        zf = zipfile.ZipFile(path)
        zf.setpassword(passwd)

        participants = {}
        for l in zf.read(embedded_file).splitlines():
            if l.startswith('#'):
                continue

            bc, name = l.strip().split('\t')[:2]
            participants[bc] = name.replace(",", "")

        print "Using identified data!"
    else:
        participants = None

    return participants
