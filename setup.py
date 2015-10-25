import re
import ast
from glob import glob

from setuptools import find_packages, setup

# https://github.com/mitsuhiko/flask/blob/master/setup.py
_version_re = re.compile(r'__version__\s+=\s+(.*)')

with open('americangut/__init__.py', 'rb') as f:
    hit = _version_re.search(f.read().decode('utf-8')).group(1)
    version = str(ast.literal_eval(hit))

classes = """
    Development Status :: 1 - Planning
    Topic :: Software Development :: Libraries
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python
    Programming Language :: Python :: 3
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classes.split('\n') if s]

setup(name='americangut',
      version=version,
      description="American Gut analyses",
      long_description="American Gut open-access data "
                       "and IPython Notebooks",
      author="American Gut development team",
      author_email="info@americangut.org",
      maintainer="American Gut development team",
      maintainer_email="info@americangut.org",
      url='https://github.com/biocore/American-Gut/',
      test_suite='nose.collector',
      packages=find_packages(),
      scripts=glob('scripts/*py'),
      setup_requires=['numpy >= 1.9.2'],
      install_requires=[
          'IPython < 4.0.0',
          'matplotlib >= 1.4.3',
          'numpy >= 1.9.2',
          'pandas >= 0.15',
          'scipy >= 0.15.1',
          'lxml',
          'h5py>=2.3.1',
          'notebook',
          'scikit-bio==0.2.3',
          'biom-format',
          'colorbrewer',
          'seaborn',
          'click',
          'qiime',
          'runipy',
          'ipymd'
      ],
      extras_require={'test': ["nose", "pep8", "flake8"]},
      classifiers=classifiers,
      package_data={'americangut': ['tests',
                                    'tests/data',
                                    'latex',
                                    'latex/pdfs-gut',
                                    'latex/pdfs-oralskin']}
      )
