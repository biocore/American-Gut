Early in the American Gut Project, it was observed that some organisms bloomed likely as a result of increased shipping time and delay between when samples were collected and when they were put on ice (more detail can be found [here](http://americangut.org/?page_id=277)). The purpose of this notebook is to apply the filter developed in order to bioinformatically subtract these observed bloom sequences from fecal samples. It is important to apply this filter when combining data with the American Gut as to remove a potential study-effect bias as all fecal data in the American Gut has had this filter applied. In addition, we're going to produce a read-length trimmed set of sequences as to reduce a potential study-effect bias in the subsequent analyses we'll be doing. The specific steps covered in this notebook are:

* Filter demultiplexed sequence data to only fecal samples
* Determine what sequences in the fecal samples recruit to the observed bloom sequences
* Remove the recruited bloom sequences from the demultiplexed sequence data
* Trim the filtered demultiplexed data (to reduce study bias when combining with short reads, such as the data in [Yatsunenko et al 2012](http://www.nature.com/nature/journal/v486/n7402/abs/nature11053.html))

The filtering is only intended to be applied to fecal data. As such, this notebook allows you to describe what metadata column and value to use so that only fecal samples are used. This will be done in a few cells.

```python
>>> import os
>>> import multiprocessing
>>> import americangut.notebook_environment as agenv
>>> import americangut.util as agu
...
>>> chp_path = agenv.activate('02-filtered')
```

Next, we'll establish the paths we will be creating.

```python
>>> fecal_sequences          = agu.get_new_path(agenv.paths['filtered']['fecal-sequences'])
>>> filtered_sequences       = agu.get_new_path(agenv.paths['filtered']['sequences-notrim'])
>>> filtered_sequences_100nt = agu.get_new_path(agenv.paths['filtered']['sequences-100nt'])
>>> observed_blooms          = agu.get_new_path(agenv.paths['filtered']['observed-blooms'])
>>> observed_blooms_biom     = agu.get_new_path(agenv.paths['filtered']['observed-blooms-biom'])
>>> observed_blooms_otu_map  = agu.get_new_path(agenv.paths['filtered']['observed-blooms-otu-map'])
```

This next call will setup and verify the path to the bloom sequences used for filtering.

```python
>>> bloom_sequences = agenv.get_bloom_sequences()
```

Now let's setup the paths to the sequences to filter. We need the metadata as well in order to reduce the data to just the fecal samples. Please replace these variables with your own paths if you wish to filter your data for blooms (a necessary precursor if you wish to combine data with the American Gut).

```python
>>> # If you are filtering your own data, please update these filepath variables as necessary
... sequences = agenv.get_existing_path(agenv.paths['raw']['sequences'])
>>> metadata  = agenv.get_existing_path(agenv.paths['raw']['metadata'])
```

We also need to specify what specific metadata category and value indicate which samples are fecal. It is possible that these values are study specific, so please modify these as needed if you filtering other datasets.

```python
>>> # If you are filtering your own data, please update these variables to reflect your mapping file
... metadata_category = 'BODY_SITE'
>>> metadata_value    = 'UBERON:feces'
```

Now that we know what sequences to focus on, we can filter the input data down to just those that need to be considered for filtering.

```python
>>> _fecal_states = ':'.join([metadata_category, metadata_value])
...
>>> !filter_fasta.py -f $sequences \
...                  -o $fecal_sequences \
...                  --mapping_fp $metadata \
...                  --valid_states $_fecal_states
```

The next thing we need to do is setup the parameters for SortMeRNA, which is the method we'll use to compare all the input data to our reference of bloom sequences.

```python
>>> _params_file = agu.get_path('sortmerna_pick_params.txt')
>>> with open(_params_file, 'w') as f:
...     f.write("pick_otus:otu_picking_method sortmerna\n")
...     f.write("pick_otus:threads %d\n" % agenv.get_cpu_count())
```

```python
>>> !pick_closed_reference_otus.py -i $fecal_sequences \
...                                -o $observed_blooms \
...                                -r $bloom_sequences \
...                                -p $_params_file \
...                                -f
```

And now, we can remove the blooms from the input sequences.

```python
>>> !filter_fasta.py -f $sequences \
...                  -m $observed_blooms_otu_map \
...                  -n \
...                  -o $filtered_sequences
```

As the data have now been filtered for blooms, we can now trim the reads back to 100nt to minimize a potential study effect when combining with the Global Gut.

```python
>>> with open(filtered_sequences) as in_, open(filtered_sequences_100nt, 'w') as out:
...     agu.trim_fasta(in_, out, 100)
```

Finally, let's do a quick sanity check that we have sequence data, and that the number of sequences is the same between both trimmed and untrimmed. We'll also dump out summary information about how many reads per sample recruited to the blooms.

```python
>>> if not (os.stat(filtered_sequences).st_size > 0 and os.stat(filtered_sequences_100nt).st_size > 0):
...     raise RuntimeError('There are no sequences follow filtering')
```
