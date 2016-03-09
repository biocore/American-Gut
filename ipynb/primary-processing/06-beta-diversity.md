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
```

To end, as always, let's sanity check and make sure the files of specific interest exist and contain data.

```python
>>> error = False
>>> message = []
>>> for key, value in (agenv.paths['beta']['notrim']['1k'].items() +
...                    agenv.paths['beta']['notrim']['10k'].items() +
...                    agenv.paths['beta']['100nt']['1k'].items() +
...                    agenv.paths['beta']['100nt']['10k'].items()):
...     if key.endswith('-pc'):
...         path = agu.get_path(value)
...         if not os.path.exists(path):
...             message.append("Could not find: %s" % path)
...             error = True
...         else:
...             if os.stat(path).st_size == 0:
...                 message.append("File appears empty: %s" % path)
...                 error = True
>>> if error:
...     raise RuntimeError('\n'.join(message))
```
