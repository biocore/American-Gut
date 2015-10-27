In order to provide taxonomy information for participants, we need to summarize the taxonomies of the individual OTUs. Within AG and for participants, we're interested in providing information at 3 levels: phylum, class and genus. These summarizations are used in different areas of the participant results.

```python
>>> import os
...
>>> import americangut.notebook_environment as agenv
>>> import americangut.util as agu
...
>>> import qiime_default_reference as qdr
...
>>> chp_path = agenv.activate('07')
```

For this notebook, we only need a single input file.

```python
>>> ag_100nt_biom = agu.get_existing_path(agenv.paths['ag-100nt-biom'])
>>> ag_cleaned_md = agu.get_existing_path(agenv.paths['ag-cleaned-md'])
```

But like the other notebooks, we'll be generating a variety of output files. So, let's get those paths setup.

```python
>>> ag_taxa               = agu.get_new_path(agenv.paths['ag-taxa'])
>>> ag_l2_taxa_tsv        = agu.get_new_path(agenv.paths['ag-L2-taxa-tsv'])
>>> ag_l2_taxa_biom       = agu.get_new_path(agenv.paths['ag-L2-taxa-biom'])
>>> ag_l2_taxa_md         = agu.get_new_path(agenv.paths['ag-L2-taxa-md'])
>>> ag_l3_taxa_tsv        = agu.get_new_path(agenv.paths['ag-L3-taxa-tsv'])
>>> ag_l3_taxa_biom       = agu.get_new_path(agenv.paths['ag-L3-taxa-biom'])
>>> ag_l6_taxa_tsv        = agu.get_new_path(agenv.paths['ag-L6-taxa-tsv'])
>>> ag_l6_taxa_biom       = agu.get_new_path(agenv.paths['ag-L6-taxa-biom'])
>>> ag_l2_taxa_oral_biom  = agu.get_new_path(agenv.paths['ag-L2-taxa-oral-biom'])
>>> ag_l2_taxa_skin_biom  = agu.get_new_path(agenv.paths['ag-L2-taxa-skin-biom'])
>>> ag_l2_taxa_fecal_biom = agu.get_new_path(agenv.paths['ag-L2-taxa-fecal-biom'])
>>> ag_l6_taxa_oral_biom  = agu.get_new_path(agenv.paths['ag-L6-taxa-oral-biom'])
>>> ag_l6_taxa_skin_biom  = agu.get_new_path(agenv.paths['ag-L6-taxa-skin-biom'])
>>> ag_l6_taxa_fecal_biom = agu.get_new_path(agenv.paths['ag-L6-taxa-fecal-biom'])
```

Now, let's perform the taxonomy summarization.

```python
>>> !summarize_taxa.py -i $ag_100nt_biom \
...                    -o $ag_taxa \
...                    -L 2,3,6
```

We're also going to need the relative abundances of each phyla available within a mapping file.

```python
>>> !summarize_taxa.py -i $ag_100nt_biom \
...                    -o $ag_taxa \
...                    -L 2 \
...                    -m $ag_cleaned_md
```

The per-sample analysis that assesses OTU significance needs to operate within a given sample type, so let's also partition the relevant table by sample type.

```python
>>> !filter_samples_from_otu_table.py -i $ag_l6_taxa_biom \
...                                   -o $ag_l6_taxa_oral_biom \
...                                   -m $ag_cleaned_md \
...                                   -s "SIMPLE_BODY_SITE:ORAL"
```

```python
>>> !filter_samples_from_otu_table.py -i $ag_l2_taxa_biom \
...                                   -o $ag_l2_taxa_oral_biom \
...                                   -m $ag_cleaned_md \
...                                   -s "SIMPLE_BODY_SITE:ORAL"
```

```python
>>> !filter_samples_from_otu_table.py -i $ag_l6_taxa_biom \
...                                   -o $ag_l6_taxa_skin_biom \
...                                   -m $ag_cleaned_md \
...                                   -s "SIMPLE_BODY_SITE:SKIN"
```

```python
>>> !filter_samples_from_otu_table.py -i $ag_l2_taxa_biom \
...                                   -o $ag_l2_taxa_skin_biom \
...                                   -m $ag_cleaned_md \
...                                   -s "SIMPLE_BODY_SITE:SKIN"
```

```python
>>> !filter_samples_from_otu_table.py -i $ag_l6_taxa_biom \
...                                   -o $ag_l6_taxa_fecal_biom \
...                                   -m $ag_cleaned_md \
...                                   -s "SIMPLE_BODY_SITE:FECAL"
```

```python
>>> !filter_samples_from_otu_table.py -i $ag_l2_taxa_biom \
...                                   -o $ag_l2_taxa_fecal_biom \
...                                   -m $ag_cleaned_md \
...                                   -s "SIMPLE_BODY_SITE:FECAL"
```

And finally, we'll sanity check that the files appear to be there.

```python
>>> assert os.stat(ag_l2_taxa_tsv).st_size > 0
>>> assert os.stat(ag_l2_taxa_biom).st_size > 0
>>> assert os.stat(ag_l2_taxa_md).st_size > 0
>>> assert os.stat(ag_l3_taxa_tsv).st_size > 0
>>> assert os.stat(ag_l3_taxa_biom).st_size > 0
>>> assert os.stat(ag_l6_taxa_tsv).st_size > 0
>>> assert os.stat(ag_l6_taxa_biom).st_size > 0
>>> assert os.stat(ag_l2_taxa_oral_biom).st_size > 0
>>> assert os.stat(ag_l2_taxa_skin_biom).st_size > 0
>>> assert os.stat(ag_l2_taxa_fecal_biom).st_size > 0
>>> assert os.stat(ag_l6_taxa_oral_biom).st_size > 0
>>> assert os.stat(ag_l6_taxa_skin_biom).st_size > 0
>>> assert os.stat(ag_l6_taxa_fecal_biom).st_size > 0
```
