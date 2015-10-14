In order to provide taxonomy information for participants, we need to summarize the taxonomies of the individual OTUs. Within AG and for participants, were interested in providing information at 3 levels: phylum, class and genus. These summarizations are used in different areas of the participant results.

```python
>>> import os
...
>>> import americangut.notebook_environment as agenv
>>> import americangut.util as agu
...
>>> import qiime_default_reference as qdr
...
>>> chp_path = agenv.activate('7')
```

For this notebook, we only need a single input file.

```python
>>> ag_100nt_biom = agu.get_existing_path(agenv.paths['ag-100nt-biom'])
```

But like the other notebooks, we'll be generating a variety of output files. So, let's get those paths setup.

```python
>>> ag_taxa         = agu.get_new_path(agenv.paths['ag-taxa'])
>>> ag_l2_taxa_tsv  = agu.get_new_path(agenv.paths['ag-L2-taxa-tsv'])
>>> ag_l2_taxa_biom = agu.get_new_path(agenv.paths['ag-L2-taxa-biom'])
>>> ag_l3_taxa_tsv  = agu.get_new_path(agenv.paths['ag-L3-taxa-tsv'])
>>> ag_l3_taxa_biom = agu.get_new_path(agenv.paths['ag-L3-taxa-biom'])
>>> ag_l6_taxa_tsv  = agu.get_new_path(agenv.paths['ag-L6-taxa-tsv'])
>>> ag_l6_taxa_biom = agu.get_new_path(agenv.paths['ag-L6-taxa-biom'])
```

Now, let's perform the taxonomy summarization.

```python
>>> !summarize_taxa.py -i $ag_100nt_biom \
...                    -o $ag_taxa \
...                    -L 2,3,6
```

And finally, we'll sanity check that the files appear to be there.

```python
>>> assert os.stat(ag_l2_taxa_tsv).st_size > 0
>>> assert os.stat(ag_l2_taxa_biom).st_size > 0
>>> assert os.stat(ag_l3_taxa_tsv).st_size > 0
>>> assert os.stat(ag_l3_taxa_biom).st_size > 0
>>> assert os.stat(ag_l6_taxa_tsv).st_size > 0
>>> assert os.stat(ag_l6_taxa_biom).st_size > 0
```
