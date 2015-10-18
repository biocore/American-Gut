The sequences and the metadata. All of the deidentified American Gut data are deposited into the European Bioinformatics Institute sequence repository, and can be retrieved by pulling down the data associated with the appropriate accessions.

First, let's setup and sanity check our environment.

```python
>>> import os
>>> import americangut.notebook_environment as agenv
>>> import americangut.util as agu
...
>>> chp_path = agenv.activate('1')
...
>>> agp_sequences = agu.get_new_path(agenv.paths['raw-sequences'])
>>> agp_metadata  = agu.get_new_path(agenv.paths['raw-metadata'])
```

Now that we have what appears to be a sane environment, let's setup a variable that defines the American Gut accessions.

```python
>>> accessions = agenv.get_accessions()
```

Now let's actually fetch the study data. `fetch_study` will only pull down accessions that do not appear in the current working directory.

```python
>>> for accession in accessions:
...     agu.fetch_study(accession, chp_path)
```

Now that we have the sequences and sample information, let's merge all the data into a single file to ease downstream processing.

```python
>>> form_path = lambda acc, ext: agenv.get_existing_path(os.path.join('1', '%s.%s' % (acc, ext)))
...
>>> bash_compatible_sequence_files = ' '.join([form_path(acc, 'fna') for acc in accessions])
>>> !cat $bash_compatible_sequence_files > $agp_sequences
...
>>> if len(accessions) > 1:
...     qiime_compatible_mapping_files = ','.join([form_path(acc, 'txt') for acc in accessions])
...     !merge_mapping_files.py -m $qiime_compatible_mapping_files -o $agp_metadata
>>> else:
...     src = form_path(accessions[0], 'txt')
...     !cp $src $agp_metadata
```

And finally, let's verify that the files we expect were created.

```python
>>> assert os.stat(agp_sequences).st_size > 0
>>> assert os.stat(agp_metadata).st_size > 0
```
