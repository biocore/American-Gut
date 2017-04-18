In order to provide taxonomy information for participants, we need to summarize the taxonomies of the individual OTUs. Within AG and for participants, we're interested in providing information at 3 levels: phylum, class and genus. These summarizations are used in different areas of the participant results.

```python
>>> import os
...
>>> import americangut.notebook_environment as agenv
>>> import americangut.util as agu
...
>>> import qiime_default_reference as qdr
...
>>> chp_path = agenv.activate('07-taxa')
```

For this notebook, we only need a single input file.

```python
>>> ag_100nt_biom  = agu.get_existing_path(agenv.paths['otus']['100nt']['ag-biom'])
>>> ag_notrim_biom = agu.get_existing_path(agenv.paths['otus']['notrim']['ag-biom'])
>>> ag_cleaned_md  = agu.get_existing_path(agenv.paths['meta']['ag-cleaned-md'])
```

But like the other notebooks, we'll be generating a variety of output files. So, let's get those paths setup.

```python
>>> ag_100nt_base  = agu.get_new_path(agenv.paths['taxa']['100nt']['base'])
>>> ag_notrim_base = agu.get_new_path(agenv.paths['taxa']['notrim']['base'])
```

Now, let's perform the taxonomy summarization.

```python
>>> for table, base in zip([ag_100nt_biom, ag_notrim_biom], [ag_100nt_base, ag_notrim_base]):
...     !summarize_taxa.py -i $table \
...                        -o $base \
...                        -L 2,3,6
```

We're also going to need the relative abundances of each phyla available within a mapping file. This is used when "painting" one of the principal coordinates plots.

```python
>>> for table, base in zip([ag_100nt_biom, ag_notrim_biom], [ag_100nt_base, ag_notrim_base]):
...     !summarize_taxa.py -i $table \
...                        -o $base \
...                        -L 2,3,6 \
...                        -m $ag_cleaned_md
```

The per-sample analysis that assesses OTU significance needs to operate within a given sample type, so let's also partition the tables by sample type.

```python
>>> def parameter_iterator():
...     for trim in ['100nt', 'notrim']:
...         for level in ['L2', 'L3', 'L6']:
...             for site in ['oral', 'skin', 'fecal']:
...                 table          = agu.get_existing_path(agenv.paths['taxa'][trim][level]['ag-biom'])
...                 filtered_table = agu.get_new_path(agenv.paths['taxa'][trim][level]['ag-%s-biom' % site])
...                 criteria       = "SIMPLE_BODY_SITE:%s" % site.upper()
...
...                 yield (table, filtered_table, criteria)
```

```python
>>> for table, filtered_table, filter_criteria in parameter_iterator():
...     !filter_samples_from_otu_table.py -i $table \
...                                       -o $filtered_table \
...                                       -m $ag_cleaned_md \
...                                       -s $filter_criteria
```

And finally, we'll sanity check that the files appear to be there.

```python
>>> error = False
>>> message = []
>>> for trim in ['100nt', 'notrim']:
...     for level in ['L2', 'L3', 'L6']:
...         for filename in agenv.paths['taxa'][trim][level].values():
...             path = agu.get_path(filename)
...             if not os.path.exists(path):
...                 message.append("Could not find: %s" % path)
...                 error = True
...             else:
...                 if os.stat(path).st_size == 0:
...                     message.append("File appears empty: %s" % path)
...                     error = True
>>> if error:
...     raise RuntimeError('\n'.join(message))
```
