The sequences and the metadata. All of the deidentified American Gut data are deposited into the European Bioinformatics Institute sequence repository, and can be retrieved by pulling down the data associated with the appropriate accessions.

First, let's setup and sanity check our environment.

```python
>>> import os
>>> import americangut.notebook_environment as agenv
>>> import americangut.util as agu
...
>>> chp_path = agenv.activate('01-raw')
...
>>> agp_sequences = agu.get_new_path(agenv.paths['raw']['sequences'])
>>> agp_metadata  = agu.get_new_path(agenv.paths['raw']['metadata'])
```

Now that we have what appears to be a sane environment, let's setup a variable that defines the American Gut accessions.

```python
>>> study_accessions = agenv.get_study_accessions()
```

Now let's actually fetch the study data. `fetch_study` will only pull down accessions that do not appear in the current working directory.

```python
>>> for accession in study_accessions:
...     new_samples = agu.fetch_study(accession, chp_path)
>>> print "New samples received: %d" % new_samples
```

Now that we have the sequences and sample information, let's merge all the data into a single file to ease downstream processing.

```python
>>> sample_sequence_files = agenv.get_files(chp_path, suffix='fna')
>>> for f in sample_sequence_files:
...     !cat $f >> $agp_sequences
...
>>> mapping_files = agenv.get_files(chp_path, suffix='txt')
>>> agu.from_xmls_to_mapping_file(mapping_files, agp_metadata)
```

And finally, let's verify that the files we expect were created.

```python
>>> if not os.stat(agp_sequences).st_size > 0:
...     raise RuntimeError('Sequences were not downloaded')
...
>>> if not os.stat(agp_metadata).st_size > 0:
...     raise RuntimeError('Metadata was not downloaded')
```
