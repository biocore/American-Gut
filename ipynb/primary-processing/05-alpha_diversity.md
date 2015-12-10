In this next chapter, we're going to calculate alpha diversity over all of the samples in the meta-analysis. As always, let's first get our environment setup.

```python
>>> import os
...
>>> import americangut.notebook_environment as agenv
>>> import americangut.util as agu
...
>>> import qiime_default_reference as qdr
...
>>> chp_path = agenv.activate('05-alpha')
>>> greengenes_tree = qdr.get_reference_tree()
>>> cpu_count = agenv.get_cpu_count()
```

And then we'll pull in the files we need for processing.

```python
>>> ag_notrim_biom           = agu.get_existing_path(agenv.paths['otus']['notrim']['ag-biom'])
>>> ag_md                    = agu.get_existing_path(agenv.paths['raw']['metadata'])
>>> ag_pgp_hmp_gg_100nt_biom = agu.get_existing_path(agenv.paths['meta']['ag-pgp-hmp-gg-100nt-biom'])
>>> ag_gg_100nt_biom         = agu.get_existing_path(agenv.paths['meta']['ag-gg-100nt-biom'])
>>> ag_pgp_hmp_gg_cleaned_md = agu.get_existing_path(agenv.paths['meta']['ag-pgp-hmp-gg-cleaned-md'])
```

Now for the fun bits. We're going to make a call to `parallel_multiple_rarefactions.py` which will produce multiple tables normalized for sequencing effort at 1000 sequences per sample. This process is called rarefaction, and we perform it as the number of sequences produced per sample can vary quite a bit. This variation can greatly influence the resulting diversity scores, so in order to reduce bias driven by sequencing effort, we make an attempt to normalize. The process involves random subsampling, so we're going to do it a few times to get rough bounds on the diversity of each sample.

Before we start computation, we're going to get two depths to sample at. The depth is abstracted out for integration testing purposes. The `low_depth` is useful for comparisons against the Human Microbiome Project where the sequencing depth was relatively low. The `high_depth` is useful for more accurate perspectives on diversity within the American Gut, Global Gut and Personal Genome Project, all of which have on average much more sequence per sample than the Human Microbiome Project.

```python
>>> low_depth, high_depth = agenv.get_rarefaction_depth()
```

Next, we're going to define a helper function that will give us the various parameters we'd like to vary, such as the rarefaction level and the input BIOM table.

```python
>>> def parameter_iterator():
...     """A helper iterator for getting the different parameters"""
...     # for the depth and the rarefaction label
...     for depth, rarefaction in zip([low_depth, high_depth], ['1k', '10k']):
...         # for each table and the base key label
...         for table, keybase in zip([ag_pgp_hmp_gg_100nt_biom, ag_notrim_biom],
...                                   ['ag-pgp-hmp-gg-100nt', 'ag-notrim']):
...             multi_directory = agu.get_path(agenv.paths['alpha'][rarefaction][keybase + '-multiple'])
...             adiv_directory  = agu.get_path(agenv.paths['alpha'][rarefaction][keybase])
...
...             yield (table, depth, rarefaction, multi_directory, adiv_directory, keybase)
```

```python
>>> for table, depth, rarefaction, multi_directory, adiv_directory, keybase in parameter_iterator():
...     !multiple_rarefactions_even_depth.py -i $table \
...                                          -o $multi_directory \
...                                          -d $depth
```

Once we have the rarefactions, we can then compute the diversity of every sample within our dataset.

```python
>>> for table, depth, rarefaction, multi_directory, adiv_directory, keybase in parameter_iterator():
...     !parallel_alpha_diversity.py -i $multi_directory \
...                                  -o $adiv_directory \
...                                  -t $greengenes_tree \
...                                  -m PD_whole_tree,chao1,observed_otus,shannon \
...                                  -O $cpu_count
```

We're going to aggregate the diversity calculations from the multiple rarefactions.

```python
>>> for table, depth, rarefaction, multi_directory, adiv_directory, keybase in parameter_iterator():
...     # sometimes QIIME leaves this directory, unsure why.
...     if os.path.exists(os.path.join(adiv_directory, 'ALDIV_*')):
...         !rmdir $adiv_directory/ALDIV_*
...     !collate_alpha.py -i $adiv_directory \
...                       -o $adiv_directory
Traceback (most recent call last):
  File "/Users/jwdebelius/anaconda/envs/americangut/bin/collate_alpha.py", line 130, in <module>
    main()
  File "/Users/jwdebelius/anaconda/envs/americangut/bin/collate_alpha.py", line 96, in main
    file_name_table = map(parse_rarefaction_fname, file_names)
  File "/Users/jwdebelius/anaconda/envs/americangut/lib/python2.7/site-packages/qiime/parse.py", line 379, in parse_rarefaction_fname
    iters = int(root_list.pop())
ValueError: invalid literal for int() with base 10: 'chao1'
Traceback (most recent call last):
  File "/Users/jwdebelius/anaconda/envs/americangut/bin/collate_alpha.py", line 130, in <module>
    main()
  File "/Users/jwdebelius/anaconda/envs/americangut/bin/collate_alpha.py", line 96, in main
    file_name_table = map(parse_rarefaction_fname, file_names)
  File "/Users/jwdebelius/anaconda/envs/americangut/lib/python2.7/site-packages/qiime/parse.py", line 379, in parse_rarefaction_fname
    iters = int(root_list.pop())
ValueError: invalid literal for int() with base 10: 'chao1'
Traceback (most recent call last):
  File "/Users/jwdebelius/anaconda/envs/americangut/bin/collate_alpha.py", line 130, in <module>
    main()
  File "/Users/jwdebelius/anaconda/envs/americangut/bin/collate_alpha.py", line 96, in main
    file_name_table = map(parse_rarefaction_fname, file_names)
  File "/Users/jwdebelius/anaconda/envs/americangut/lib/python2.7/site-packages/qiime/parse.py", line 379, in parse_rarefaction_fname
    iters = int(root_list.pop())
ValueError: invalid literal for int() with base 10: 'chao1'
Traceback (most recent call last):
  File "/Users/jwdebelius/anaconda/envs/americangut/bin/collate_alpha.py", line 130, in <module>
    main()
  File "/Users/jwdebelius/anaconda/envs/americangut/bin/collate_alpha.py", line 96, in main
    file_name_table = map(parse_rarefaction_fname, file_names)
  File "/Users/jwdebelius/anaconda/envs/americangut/lib/python2.7/site-packages/qiime/parse.py", line 379, in parse_rarefaction_fname
    iters = int(root_list.pop())
ValueError: invalid literal for int() with base 10: 'chao1'
```

Finally, we're going to

To end, we're going to verify that the files of interest exist and contain data.

```python
>>> for table, depth, rarefaction, multi_directory, adiv_directory, keybase in parameter_iterator():
...     for metric in ['pd', 'chao1', 'observedotus', 'shannon']:
...         metric_file = agu.get_existing_path(agenv.paths['alpha'][rarefaction][keybase + '-%s' % metric])
...         assert os.stat(metric_file).st_size > 0
```

```python

```
