Pick OTUs (for full and trimmed data) at approximately genus level resolution (97% similarity) using SortMeRNA, closed reference against Greengenes 13_8.

```python
>>> import os
>>> import multiprocessing
...
>>> import americangut.util as agu
>>> import americangut.notebook_environment as agenv
...
>>> chp_path = agenv.activate('3')
```

Before we go too far, let's make sure the files we need are present.

```python
>>> filtered_sequences       = agu.get_existing_path(agenv.paths['filtered-sequences'])
>>> filtered_sequences_100nt = agu.get_existing_path(agenv.paths['filtered-sequences-100nt'])
```

And, let's make sure that the output files we need do not already exist.

```python
>>> ag_otus       = agu.get_new_path(agenv.paths['ag-otus'])
>>> ag_biom       = agu.get_new_path(agenv.paths['ag-biom'])
>>> ag_otus_100nt = agu.get_new_path(agenv.paths['ag-otus-100nt'])
>>> ag_100nt_biom = agu.get_new_path(agenv.paths['ag-100nt-biom'])
```

We're going to now setup a parameters file for the OTU picking runs. It is possible to specify a precomputed SortMeRNA index by indicating it's path as the environment variable `$AG_SMR_INDEX`. The reason we're using an environment variable is that it makes it much easier to inject an index during continuous integration testing.

```python
>>> _params_file = os.path.join(chp_path, 'sortmerna_pick_params.txt')
...
>>> with open(_params_file, 'w') as f:
...     f.write("pick_otus:otu_picking_method sortmerna\n")
...     f.write("pick_otus:threads %d\n" % multiprocessing.cpu_count())
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
>>> !pick_closed_reference_otus.py -i $filtered_sequences \
...                                -o $ag_otus \
...                                -r $ref_seqs \
...                                -t $ref_tax \
...                                -p $_params_file
```

```python
>>> !pick_closed_reference_otus.py -i $filtered_sequences_100nt \
...                                -o $ag_otus_100nt \
...                                -r $ref_seqs \
...                                -t $ref_tax \
...                                -p $_params_file
```

And we'll end with some sanity checking of the outputs.

```python
>>> !biom summarize-table -i $ag_biom | head -n 25
```

```python
>>> !biom summarize-table -i $ag_100nt_biom | head -n 25
```
