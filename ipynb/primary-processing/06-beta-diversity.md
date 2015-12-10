In this chapter, we'll compute unweighted and weighted UniFrac distances over the samples and additionally generate principal coordinates from the resulting distance matrices.

```python
>>> import os
...
>>> import americangut.notebook_environment as agenv
>>> import americangut.util as agu
...
>>> import qiime_default_reference as qdr
...
>>> chp_path = agenv.activate('06-beta')
>>> greengenes_tree = qdr.get_reference_tree()
>>> cpu_count = agenv.get_cpu_count()
```

Let's first setup the files that we need for processing.

```python
>>> ag_100nt_biom            = agu.get_existing_path(agenv.paths['otus']['100nt']['ag-biom'])
>>> ag_notrim_biom           = agu.get_existing_path(agenv.paths['otus']['notrim']['ag-biom'])
>>> ag_gg_100nt_biom         = agu.get_existing_path(agenv.paths['meta']['ag-gg-100nt-biom'])
>>> ag_pgp_hmp_gg_100nt_biom = agu.get_existing_path(agenv.paths['meta']['ag-pgp-hmp-gg-100nt-biom'])
>>> ag_cleaned_md            = agu.get_existing_path(agenv.paths['meta']['ag-cleaned-md'])
>>> ag_gg_cleaned_md         = agu.get_existing_path(agenv.paths['meta']['ag-gg-cleaned-md'])
>>> ag_pgp_hmp_gg_cleaned_md = agu.get_existing_path(agenv.paths['meta']['ag-pgp-hmp-gg-cleaned-md'])
```

Like alpha diversity, beta diversity is sensitive to sequencing effort. And while rarefaction is a random process, a single rarefaction gives a rough idea of how samples relate to each other. It is possible to perform beta diversity calculations over multiple rarefactions but that is outside of the scope of this notebook.

Before we start computation, we're going to get the depth to sample at. The depth is abstracted out for integration testing purposes.

```python
>>> low_depth, high_depth = agenv.get_rarefaction_depth()
```

```python
>>> def parameter_iterator():
...     """A helper iterator for getting the different parameters"""
...     # for the depth and the rarefaction label
...     for trim, trim_table in zip(['100nt', 'notrim'], [ag_pgp_hmp_gg_100nt_biom, ag_notrim_biom]):
...         if not os.path.exists(os.path.join(chp_path, trim)):
...             os.mkdir(os.path.join(chp_path, trim))
...
...         keybases = ['ag']
...         if trim == '100nt':
...             keybases += ['ag-pgp-hmp-gg']
...
...         for depth, rarefaction in zip([low_depth, high_depth], ['1k', '10k']):
...             if not os.path.exists(os.path.join(chp_path, trim, rarefaction)):
...                 os.mkdir(os.path.join(chp_path, trim, rarefaction))
...
...             # for each table and the base key label
...             for keybase in keybases:
...                 rarefied_table  = agu.get_path(agenv.paths['beta'][trim][rarefaction][keybase + '-biom'])
...                 bdiv_directory  = agu.get_path(agenv.paths['beta'][trim][rarefaction][keybase])
...
...                 yield (trim_table, depth, rarefaction, rarefied_table, bdiv_directory, keybase)
```

```python
>>> for table, depth, rarefaction, rarefied_table, bdiv_directory, keybase in parameter_iterator():
...     !single_rarefaction.py -i $table \
...                            -o $rarefied_table \
...                            -d $depth
```

Next, we'll compute the distances between samples. By default, this method will produce unweighted and weighted UniFrac distance matrices.

```python
>>> for table, depth, rarefaction, rarefied_table, bdiv_directory, keybase in parameter_iterator():
...     # on test data, single threaded is 15s cputime vs. 15m with 7 threads via parallel_beta_diversity. :(
...     !beta_diversity.py -i $rarefied_table \
...                        -o $bdiv_directory \
...                        -t $greengenes_tree
```

These set of next cells are going to filter down the distance matrix to interesting subsets. Notably, the AG and GG subset, just the AG samples, the oral subset within AG , the skin subset within AG and the fecal subset within AG.

```python
>>> # Create a AG + GG distance matrices, only doing this for 100nt as GG is 100nt
... for rarefaction in ['1k', '10k']:
...     for metric in ['unifrac', 'wunifrac']:
...         input_dm  = agu.get_existing_path(agenv.paths['beta']['100nt'][rarefaction]['ag-pgp-hmp-gg-%s' % metric])
...         output_dm = agu.get_new_path(agenv.paths['beta']['100nt'][rarefaction]['ag-gg-%s' % metric])
...
...         !filter_distance_matrix.py -i $input_dm \
...                                    -o $output_dm \
...                                    --sample_id_fp $ag_gg_cleaned_md
```

```python
>>> # Iterate over the distance matrices and generate parameters for per-body site filtering
... def dm_parameter_iterator():
...     """A helper iterator for getting the different parameters"""
...     # for the depth and the rarefaction label
...     for trim, keybases in zip(['100nt', 'notrim'], [['ag-pgp-hmp-gg', 'ag-gg', 'ag'], ['ag']]):
...         for rarefaction in ['1k', '10k']:
...             paths = agenv.paths['beta'][trim][rarefaction]
...             for keybase in keybases:
...                 for site in ['ORAL', 'SKIN', 'FECAL']:
...                     for metric in ['unifrac', 'wunifrac']:
...                         input_dm  = agu.get_existing_path(paths[keybase + '-%s' % metric])
...                         output_dm = agu.get_new_path(paths[keybase + '-%s-%s' % (site.lower(), metric)])
...
...                         yield (input_dm, output_dm, site)
```

```python
>>> # filter each distance matrix of interest to each body site
... for input_dm, output_dm, site in dm_parameter_iterator():
...     site_filter = "SIMPLE_BODY_SITE:%s" % site
...     !filter_distance_matrix.py -i $input_dm \
...                                -o $output_dm \
...                                -m $ag_cleaned_md \
...                                -s $site_filter
```

For one of the figures, the volume of US samples dominates those from other countries. To help mitigate the sample size effect, we're going to subsample the distance matrix prior to producing principal coordinates.

```python
>>> ag_gg_100nt_1k_bdiv_u     = agu.get_existing_path(agenv.paths['beta']['100nt']['1k']['ag-gg-unifrac'])
>>> ag_gg_100nt_10k_bdiv_u    = agu.get_existing_path(agenv.paths['beta']['100nt']['10k']['ag-gg-unifrac'])
>>> ag_gg_100nt_1k_ss_bdiv_u  = agu.get_new_path(agenv.paths['beta']['100nt']['1k']['ag-gg-subsampled-unifrac'])
>>> ag_gg_100nt_10k_ss_bdiv_u = agu.get_new_path(agenv.paths['beta']['100nt']['10k']['ag-gg-subsampled-unifrac'])
...
>>> !mod2_pcoa.py subsample_dm --distmat $ag_gg_100nt_1k_bdiv_u \
...                            --max 500 \
...                            --category COUNTRY \
...                            --mapping_file $ag_gg_cleaned_md \
...                            --output $ag_gg_100nt_1k_ss_bdiv_u
...
>>> !mod2_pcoa.py subsample_dm --distmat $ag_gg_100nt_10k_bdiv_u \
...                            --max 500 \
...                            --category COUNTRY \
...                            --mapping_file $ag_gg_cleaned_md \
...                            --output $ag_gg_100nt_10k_ss_bdiv_u
Traceback (most recent call last):
  File "/Users/jwdebelius/anaconda/envs/americangut/bin/mod2_pcoa.py", line 4, in <module>
    __import__('pkg_resources').require('americangut==0.0.1')
  File "/Users/jwdebelius/anaconda/envs/americangut/lib/python2.7/site-packages/setuptools-18.5-py2.7.egg/pkg_resources/__init__.py", line 3095, in <module>
  File "/Users/jwdebelius/anaconda/envs/americangut/lib/python2.7/site-packages/setuptools-18.5-py2.7.egg/pkg_resources/__init__.py", line 3081, in _call_aside
  File "/Users/jwdebelius/anaconda/envs/americangut/lib/python2.7/site-packages/setuptools-18.5-py2.7.egg/pkg_resources/__init__.py", line 3108, in _initialize_master_working_set
  File "/Users/jwdebelius/anaconda/envs/americangut/lib/python2.7/site-packages/setuptools-18.5-py2.7.egg/pkg_resources/__init__.py", line 660, in _build_master
  File "/Users/jwdebelius/anaconda/envs/americangut/lib/python2.7/site-packages/setuptools-18.5-py2.7.egg/pkg_resources/__init__.py", line 673, in _build_from_requirements
  File "/Users/jwdebelius/anaconda/envs/americangut/lib/python2.7/site-packages/setuptools-18.5-py2.7.egg/pkg_resources/__init__.py", line 846, in resolve
pkg_resources.DistributionNotFound: The 'IPython<4.0.0' distribution was not found and is required by americangut
Traceback (most recent call last):
  File "/Users/jwdebelius/anaconda/envs/americangut/bin/mod2_pcoa.py", line 4, in <module>
    __import__('pkg_resources').require('americangut==0.0.1')
  File "/Users/jwdebelius/anaconda/envs/americangut/lib/python2.7/site-packages/setuptools-18.5-py2.7.egg/pkg_resources/__init__.py", line 3095, in <module>
  File "/Users/jwdebelius/anaconda/envs/americangut/lib/python2.7/site-packages/setuptools-18.5-py2.7.egg/pkg_resources/__init__.py", line 3081, in _call_aside
  File "/Users/jwdebelius/anaconda/envs/americangut/lib/python2.7/site-packages/setuptools-18.5-py2.7.egg/pkg_resources/__init__.py", line 3108, in _initialize_master_working_set
  File "/Users/jwdebelius/anaconda/envs/americangut/lib/python2.7/site-packages/setuptools-18.5-py2.7.egg/pkg_resources/__init__.py", line 660, in _build_master
  File "/Users/jwdebelius/anaconda/envs/americangut/lib/python2.7/site-packages/setuptools-18.5-py2.7.egg/pkg_resources/__init__.py", line 673, in _build_from_requirements
  File "/Users/jwdebelius/anaconda/envs/americangut/lib/python2.7/site-packages/setuptools-18.5-py2.7.egg/pkg_resources/__init__.py", line 846, in resolve
pkg_resources.DistributionNotFound: The 'IPython<4.0.0' distribution was not found and is required by americangut
```

And finally, we'll produce principal coordinates from the different distance matrices produced above.

```python
>>> for directory, subdirectories, files in os.walk(chp_path):
...     for filename in files:
...         if filename.endswith('.txt'):
...             dm_basename, ext = os.path.splitext(filename)
...
...             dm_filename = os.path.join(directory, filename)
...             pc_filename = os.path.join(directory, dm_basename + '-pc.txt')
...
...             !principal_coordinates.py -i $dm_filename -o $pc_filename
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -8.36228927449e-06 and the largest is 0.201037662373.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.00585800148341 and the largest is 0.150213923426.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.0570613541406 and the largest is 2.44188989496.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.000563013884505 and the largest is 0.0501361977778.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.0634852183095 and the largest is 4.41030436665.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.0125769761048 and the largest is 0.176897688215.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.0657716582459 and the largest is 2.05188229841.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.0125769761048 and the largest is 0.176897688215.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.0673975590609 and the largest is 2.41280504429.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -7.92008482334e-05 and the largest is 0.0377396239742.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.0197936957677 and the largest is 2.69925100996.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -7.92008482334e-05 and the largest is 0.0377396239742.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.0730382493733 and the largest is 4.2540617994.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -7.25715612149e-06 and the largest is 0.286896774631.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.0545446154772 and the largest is 2.40748198756.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.00103158503795 and the largest is 0.0526739537408.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.0733557672067 and the largest is 4.26067098905.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -7.19916657783e-06 and the largest is 0.214818378104.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.0720774466297 and the largest is 1.89223776186.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -7.19916657783e-06 and the largest is 0.214818378104.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.0790847402418 and the largest is 2.60265317828.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.000121838575273 and the largest is 0.0247020018299.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.00158319774517 and the largest is 0.0926805661722.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.0407323272701 and the largest is 2.84345262978.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.000121838575273 and the largest is 0.0247020018299.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.00158319774517 and the largest is 0.0926805661722.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.0751636902393 and the largest is 4.47786788005.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.0575209209403 and the largest is 1.25141269592.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.000510788317219 and the largest is 0.0389040597018.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.0220203589103 and the largest is 2.38103151051.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.0474971610626 and the largest is 1.13288747823.
  RuntimeWarning
/Users/jwdebelius/.virtualenvs/this_is_fucking_ridiculous/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:107: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.0168165244083 and the largest is 2.23870030046.
  RuntimeWarning
```

To end, as always, let's sanity check and make sure the files of specific interest exist and contain data.

```python
>>> error = False
>>> for key, value in (agenv.paths['beta']['notrim']['1k'].items() +
...                    agenv.paths['beta']['notrim']['10k'].items() +
...                    agenv.paths['beta']['100nt']['1k'].items() +
...                    agenv.paths['beta']['100nt']['10k'].items()):
...     if key.endswith('-pc'):
...         path = agu.get_path(value)
...         if not os.path.exists(path):
...             print "Could not find: %s" % path
...             error = True
...         else:
...             if os.stat(path).st_size == 0:
...                 print "File appears empty: %s" % path
...                 error = True
>>> assert not error
Could not find: /Users/jwdebelius/Repositories/American-Gut/ipynb/primary-processing/agp_processing/06-beta/100nt/1k/ag-pgp-hmp-gg/unweighted_unifrac_ag-gg-subsampled-pc.txt
Could not find: /Users/jwdebelius/Repositories/American-Gut/ipynb/primary-processing/agp_processing/06-beta/100nt/10k/ag-pgp-hmp-gg/unweighted_unifrac_ag-gg-subsampled-pc.txt
```

```python

```

```python

```
