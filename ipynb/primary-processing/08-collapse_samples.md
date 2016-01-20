Some of the per-results figures require various slices and perspectives of the data. This notebook performs these necessary additional summarizations.

```python
>>> import os
...
>>> import pandas as pd
...
>>> import americangut.notebook_environment as agenv
>>> import americangut.util as agu
...
>>> import qiime_default_reference as qdr
...
>>> chp_path = agenv.activate('08-collapsed')
```

Let's make sure we have the paths we need.

```python
>>> ag_cleaned_md  = agu.get_existing_path(agenv.paths['meta']['ag-cleaned-md'])
```

We're going to start by generating per-bodysite mapping files with alpha diversity attached.

```python
>>> low_depth, high_depth = agenv.get_rarefaction_depth()
```

```python
>>> for trim in ['ag-pgp-hmp-gg-100nt', 'ag-notrim']:
...     trimpath = os.path.join(chp_path, trim.split('-')[-1])
...     alpha = {}
...     if not os.path.exists(trimpath):
...         os.mkdir(trimpath)
...     for depth, rarefaction in zip([low_depth, high_depth], ['1k', '10k']):
...         rarepath = os.path.join(trimpath, rarefaction)
...         if not os.path.exists(rarepath):
...             os.mkdir(rarepath)
...         for met in agenv.alpha_metrics.keys():
...             alpha_fp = agu.get_existing_path(
...                 agenv.paths['alpha'][rarefaction]['%s-%s' % (trim, met)]
...             )
...             alpha_df = pd.read_csv(alpha_fp, sep='\t', dtype=str).set_index('Unnamed: 0').transpose()
...             alpha['%s_%s' % (agenv.alpha_metrics[met], rarefaction)] = alpha_df
...     map_ = pd.read_csv(ag_cleaned_md, sep='\t', dtype=str).set_index('#SampleID')
...     alpha_map = agu.add_alpha_diversity(map_, alpha)
...     alpha_map.to_csv(
...         agu.get_new_path(agenv.paths['collapsed'][trim.split('-')[-1]]['alpha-map']),
...         sep='\t',
...         index_label='#SampleID'
...     )
```

Next, we're going to operate on rarefied data again.

```python
>>> def rarefaction_parameters():
...     for trim in ['100nt', 'notrim']:
...         trimpath = os.path.join(chp_path, trim)
...
...         for depth, rarefaction in zip([low_depth, high_depth], ['1k', '10k']):
...             rarepath = os.path.join(trimpath, rarefaction)
...
...             table  = agu.get_existing_path(agenv.paths['otus'][trim]['ag-biom'])
...             output = agu.get_new_path(agenv.paths['collapsed'][trim][rarefaction]['ag-biom'])
...
...             yield (table, output, depth)
```

```python
>>> for table, output, depth in rarefaction_parameters():
...     !single_rarefaction.py -i $table \
...                            -o $output \
...                            -d $depth
```

Then, we're going to partition the data into per-body site tables.

```python
>>> def filter_parameters():
...     for trim in ['100nt', 'notrim']:
...         for rarefaction in ['1k', '10k']:
...             for site in ['oral', 'fecal', 'skin']:
...                 table  = agu.get_existing_path(agenv.paths['collapsed'][trim][rarefaction]['ag-biom'])
...                 output = agu.get_new_path(agenv.paths['collapsed'][trim][rarefaction]['ag-%s-biom' % site])
...                 criteria = "SIMPLE_BODY_SITE:%s" % site.upper()
...
...                 yield (table, output, criteria)
```

```python
>>> for table, output, criteria in filter_parameters():
...     !filter_samples_from_otu_table.py -i $table \
...                                       -o $output \
...                                       -m $ag_cleaned_md \
...                                       -s $criteria
```

Finally, within each body site, we're going to collapse over categories of interest.

```python
>>> def collapse_parameters():
...     # [(site, (category, file-label))]
...     categories = [('fecal', (('SEX', 'sex'),
...                              ('DIET_TYPE', 'diet'),
...                              ('AGE_CAT', 'age'),
...                              ('BMI_CAT', 'bmi'))),
...                   ('oral',  (('SEX', 'sex'),
...                              ('DIET_TYPE', 'diet'),
...                              ('AGE_CAT', 'age'),
...                              ('FLOSSING_FREQUENCY', 'flossing'))),
...                   ('skin',  (('SEX', 'sex'),
...                              ('COSMETICS_FREQUENCY', 'cosmetics'),
...                              ('AGE_CAT', 'age'),
...                              ('DOMINANT_HAND', 'hand')))]
...
...     for trim in ['100nt', 'notrim']:
...         for rarefaction in ['1k', '10k']:
...             paths = agenv.paths['collapsed'][trim][rarefaction]
...
...             for site, cats in categories:
...                 table = agu.get_existing_path(paths['ag-%s-biom' % site])
...
...                 for cat, cat_label in cats:
...                     output = agu.get_new_path(paths['ag-%s-%s-biom' % (site, cat_label)])
...
...                     yield (table, output, cat)
```

```python
>>> ignored = 'foo'  # a necessary but unused parameter
>>> for table, output, cat in collapse_parameters():
...     !collapse_samples.py -m $ag_cleaned_md \
...                          --output_biom_fp $output \
...                          --normalize \
...                          -b $table \
...                          --collapse_fields $cat \
...                          --output_mapping_fp $ignored
```

As usual, let's make sure we have files.

```python
>>> error = False
...
>>> for trim in ['100nt', 'notrim']:
...     for rarefaction in ['1k', '10k']:
...         for path in agenv.paths['collapsed'][trim][rarefaction].values():
...             filepath = agu.get_path(path)
...             if not os.path.exists(filepath):
...                 print "Could not find: %s" % filepath
...                 error = True
...             else:
...                 if os.stat(filepath) == 0:
...                     print "File appears to be empty: %s" % filepath
...                     error = True
...
>>> assert not error
```
