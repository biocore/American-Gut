Now that we have our American Gut OTU table, we can begin to combine it with other studies in preparation for a meta-analysis. This will let us explore American Gut samples in the context of projects like the Human Microbiome Project.

As before, let's first sanity check our environment.

```python
>>> import americangut.notebook_environment as agenv
>>> import americangut.util as agu
...
>>> chp_path = agenv.activate('04-meta')
...
>>> ag_md          = agu.get_existing_path(agenv.paths['raw']['metadata'])
>>> ag_100nt_biom  = agu.get_existing_path(agenv.paths['otus']['100nt']['ag-biom'])
>>> ag_notrim_biom = agu.get_existing_path(agenv.paths['otus']['notrim']['ag-biom'])
```

But, we also want to include other interesting datasets, specifically the Human Microbiome Project, microbiome samples from the Personal Genome Project, and the Global Gut samples.

```python
>>> hmp_100nt_biom, hmp_md = agenv.get_hmp()
>>> pgp_100nt_biom, pgp_md = agenv.get_pgp()
>>> gg_100nt_biom, gg_md   = agenv.get_global_gut()
```

We're also going to generate some new files, so let's get them setup.

```python
>>> # the merged BIOM tables
... ag_gg_100nt_biom         = agu.get_new_path(agenv.paths['meta']['ag-gg-100nt-biom'])
>>> pgp_hmp_100nt_biom       = agu.get_new_path(agenv.paths['meta']['pgp-hmp-100nt-biom'])
>>> ag_pgp_hmp_gg_100nt_biom = agu.get_new_path(agenv.paths['meta']['ag-pgp-hmp-gg-100nt-biom'])
...
>>> # the cleaned and simplified metadata
... ag_cleaned_md  = agu.get_new_path(agenv.paths['meta']['ag-cleaned-md'])
>>> gg_cleaned_md  = agu.get_new_path(agenv.paths['meta']['gg-cleaned-md'])
>>> pgp_cleaned_md = agu.get_new_path(agenv.paths['meta']['pgp-cleaned-md'])
>>> hmp_cleaned_md = agu.get_new_path(agenv.paths['meta']['hmp-cleaned-md'])
...
>>> # the merged simplified metadata
... ag_gg_cleaned_md         = agu.get_new_path(agenv.paths['meta']['ag-gg-cleaned-md'])
>>> pgp_hmp_cleaned_md       = agu.get_new_path(agenv.paths['meta']['pgp-hmp-cleaned-md'])
>>> ag_pgp_hmp_gg_cleaned_md = agu.get_new_path(agenv.paths['meta']['ag-pgp-hmp-gg-cleaned-md'])
```

The first step in the process is to merge the the individual tables into larger ones, and then to merge the larger tables into a final one. We're not using QIIME's `parallel_merge_otu_tables.py` here as we also need one of the intermediate tables for subsequent processing.

The first merge we'll do is between the Global Gut and the American Gut.

```python
>>> to_merge = ','.join([ag_100nt_biom, gg_100nt_biom])
>>> !merge_otu_tables.py -i $to_merge -o $ag_gg_100nt_biom
```

The second merge we'll do is between the Personal Genome Project microbiome samples and the Human Microbiome Project v35 16S samples.

```python
>>> to_merge = ','.join([pgp_100nt_biom, hmp_100nt_biom])
>>> !merge_otu_tables.py -i $to_merge -o $pgp_hmp_100nt_biom
```

And the last merge will bring them all into a single table.

```python
>>> to_merge = ','.join([ag_gg_100nt_biom, pgp_hmp_100nt_biom])
>>> !merge_otu_tables.py -i $to_merge -o $ag_pgp_hmp_gg_100nt_biom
```

Before we proceed, let's dump out a little information about the two merged tables we'll be using for subsequent analysis.

```python
>>> summary = !biom summarize-table -i $ag_pgp_hmp_gg_100nt_biom
>>> print '\n'.join(summary[:10])
```

```python
>>> summary = !biom summarize-table -i $ag_gg_100nt_biom
>>> print '\n'.join(summary[:10])
```

We also need to make sure the metadata (the information about the samples) are also merged and consistent. Prior to merge, we're going to add in some additional detail about every sample, such as a column in the mapping file that is the combination of the study title and the body site. We're also going to "generalize" body sites to the type of site they're from (e.g., the back of the hand is just "skin"). This process will also clean the metadata to remove blanks and unknown sample types.

```python
>>> original_clean_category_acronym = \
...     [(ag_md, ag_cleaned_md, 'AGP', 'body_site'),
...      (gg_md, gg_cleaned_md, 'GG', 'body_site'),
...      (pgp_md, pgp_cleaned_md, 'PGP', 'body_site'),
...      (hmp_md, hmp_cleaned_md, 'HMP', 'bodysite')]
...
>>> for original, cleaned, acronym, category in original_clean_category_acronym:
...     with open(original, 'U') as ori_fp, open(cleaned, 'w') as clean_fp:
...         agu.clean_and_reformat_mapping(ori_fp, clean_fp, category, acronym)
```

Now let's merge the mapping files so we can move on to the diversity analyses!

```python
>>> to_merge = ','.join([ag_cleaned_md, gg_cleaned_md])
>>> !merge_mapping_files.py -m $to_merge -o $ag_gg_cleaned_md
```

```python
>>> to_merge = ','.join([pgp_cleaned_md, hmp_cleaned_md])
>>> !merge_mapping_files.py -m $to_merge -o $pgp_hmp_cleaned_md
```

```python
>>> to_merge = ','.join([pgp_hmp_cleaned_md, ag_gg_cleaned_md])
>>> !merge_mapping_files.py -m $to_merge -o $ag_pgp_hmp_gg_cleaned_md
```

And, we'll finish up with a few sanity checks on the resulting metadata.

```python
>>> !print_metadata_stats.py -m $ag_pgp_hmp_gg_cleaned_md -c TITLE_BODY_SITE
```

```python
>>> !print_metadata_stats.py -m $ag_gg_cleaned_md -c TITLE_BODY_SITE
```
