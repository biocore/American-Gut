Pick OTUs (for full and trimmed data) at approximately genus level resolution (97% similarity) using SortMeRNA, closed reference against Greengenes 13_8.

```python
>>> import os
>>> import multiprocessing
...
>>> import americangut.util as agu
>>> import americangut.notebook_environment as agenv
...
>>> chp_path = agenv.activate('03-otus')
```

We're going to now setup a parameters file for the OTU picking runs. It is possible to specify a precomputed SortMeRNA index by indicating it's path as the environment variable `$AG_SMR_INDEX`. The reason we're using an environment variable is that it makes it much easier to inject an index during continuous integration testing.

```python
>>> _params_file = os.path.join(chp_path, 'sortmerna_pick_params.txt')
...
>>> with open(_params_file, 'w') as f:
...     f.write("pick_otus:otu_picking_method sortmerna\n")
...     f.write("pick_otus:threads %d\n" % agenv.get_cpu_count())
...
...     if agenv.get_sortmerna_index():
...         f.write("pick_otus:sortmerna_db %s\n" % agenv.get_sortmerna_index())
```

Determine reference set (in the event of testing).

```python
>>> ref_seqs, ref_tax = agenv.get_reference_set()
```

And now we can actually pick the OTUs. This will take sometime. Note, we're issuing two separate commands as we're picking against the untrimmed and the trimmed data.

```python
>>> for trim in ['notrim', '100nt']:
...     seqs    = agu.get_existing_path(agenv.paths['filtered']['sequences-%s' % trim])
...     ag_otus = agu.get_new_path(agenv.paths['otus'][trim]['ag'])
...
...     !pick_closed_reference_otus.py -i $seqs \
...                                    -o $ag_otus \
...                                    -r $ref_seqs \
...                                    -t $ref_tax \
...                                    -p $_params_file
```

And we'll end with some sanity checking of the outputs.

```python
>>> for trim in ['notrim', '100nt']:
...     ag_biom = agu.get_existing_path(agenv.paths['otus'][trim]['ag-biom'])
...     summary = !biom summarize-table -i $ag_biom
...     print "Trim: %s" % trim
...     print '\n'.join(summary[:10])
...     print
Trim: notrim
Num samples: 15
Num observations: 11
Total count: 6065
Table density (fraction of non-zero values): 0.485

Counts/sample summary:
 Min: 2.0
 Max: 1273.0
 Median: 318.000
 Mean: 404.333

Trim: 100nt
Num samples: 15
Num observations: 11
Total count: 6065
Table density (fraction of non-zero values): 0.485

Counts/sample summary:
 Min: 2.0
 Max: 1273.0
 Median: 318.000
 Mean: 404.333
```

```python

```
