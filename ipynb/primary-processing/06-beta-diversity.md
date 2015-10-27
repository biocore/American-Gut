In this chapter, we'll compute unweighted and weighted UniFrac distances over the samples and additionally generate principal coordinates from the resulting distance matrices.

```python
>>> import os
...
>>> import americangut.notebook_environment as agenv
>>> import americangut.util as agu
...
>>> import qiime_default_reference as qdr
...
>>> chp_path = agenv.activate('06')
>>> greengenes_tree = qdr.get_reference_tree()
>>> cpu_count = agenv.get_cpu_count()
```

Let's first setup the files that we need for processing.

```python
>>> ag_100nt_biom            = agu.get_existing_path(agenv.paths['ag-100nt-biom'])
>>> ag_gg_100nt_biom         = agu.get_existing_path(agenv.paths['ag-gg-100nt-biom'])
>>> ag_pgp_hmp_gg_100nt_biom = agu.get_existing_path(agenv.paths['ag-pgp-hmp-gg-100nt-biom'])
>>> ag_cleaned_md            = agu.get_existing_path(agenv.paths['ag-cleaned-md'])
>>> ag_gg_cleaned_md         = agu.get_existing_path(agenv.paths['ag-gg-cleaned-md'])
```

Next, we'll setup the paths that we're going to create and need.

```python
>>> ag_pgp_hmp_gg_100nt_1k_biom    = agu.get_new_path(agenv.paths['ag-pgp-hmp-gg-100nt-1k-biom'])
...
>>> # beta diversity distance matrices
... ag_100nt_1k_bdiv_u             = agu.get_new_path(agenv.paths['ag-100nt-1k-bdiv-unifrac'])
>>> ag_100nt_fecal_1k_bdiv_u       = agu.get_new_path(agenv.paths['ag-100nt-fecal-1k-bdiv-unifrac'])
>>> ag_100nt_oral_1k_bdiv_u        = agu.get_new_path(agenv.paths['ag-100nt-oral-1k-bdiv-unifrac'])
>>> ag_100nt_skin_1k_bdiv_u        = agu.get_new_path(agenv.paths['ag-100nt-skin-1k-bdiv-unifrac'])
>>> ag_gg_100nt_1k_bdiv_u          = agu.get_new_path(agenv.paths['ag-gg-100nt-1k-bdiv-unifrac'])
>>> ag_gg_100nt_1k_bdiv_u_ss       = agu.get_new_path(agenv.paths['ag-gg-100nt-1k-bdiv-subsampled-unifrac'])
>>> ag_gg_100nt_1k_bdiv_wu         = agu.get_new_path(agenv.paths['ag-gg-100nt-1k-bdiv-wunifrac'])
>>> ag_pgp_hmp_gg_100nt_1k_bdiv    = agu.get_new_path(agenv.paths['ag-pgp-hmp-gg-100nt-1k-bdiv'])
>>> ag_pgp_hmp_gg_100nt_1k_bdiv_u  = agu.get_new_path(agenv.paths['ag-pgp-hmp-gg-100nt-1k-bdiv-unifrac'])
>>> ag_pgp_hmp_gg_100nt_1k_bdiv_wu = agu.get_new_path(agenv.paths['ag-pgp-hmp-gg-100nt-1k-bdiv-wunifrac'])
...
>>> # principal coordinates
... ag_100nt_1k_u_pc               = agu.get_new_path(agenv.paths['ag-100nt-1k-unifrac-pc'])
>>> ag_100nt_fecal_1k_u_pc         = agu.get_new_path(agenv.paths['ag-100nt-fecal-1k-unifrac-pc'])
>>> ag_100nt_oral_1k_u_pc          = agu.get_new_path(agenv.paths['ag-100nt-oral-1k-unifrac-pc'])
>>> ag_100nt_skin_1k_u_pc          = agu.get_new_path(agenv.paths['ag-100nt-skin-1k-unifrac-pc'])
>>> ag_gg_100nt_1k_u_pc            = agu.get_new_path(agenv.paths['ag-gg-100nt-1k-unifrac-pc'])
>>> ag_gg_100nt_1k_wu_pc           = agu.get_new_path(agenv.paths['ag-gg-100nt-1k-wunifrac-pc'])
>>> ag_gg_100nt_1k_u_ss_pc         = agu.get_new_path(agenv.paths['ag-gg-100nt-1k-subsampled-unifrac-pc'])
>>> ag_pgp_hmp_gg_100nt_1k_u_pc    = agu.get_new_path(agenv.paths['ag-pgp-hmp-gg-100nt-1k-unifrac-pc'])
>>> ag_pgp_hmp_gg_100nt_1k_wu_pc   = agu.get_new_path(agenv.paths['ag-pgp-hmp-gg-100nt-1k-wunifrac-pc'])
```

Like alpha diversity, beta diversity is sensitive to sequencing effort. And while rarefaction is a random process, a single rarefaction gives a rough idea of how samples relate to each other. It is possible to perform beta diversity calculations over multiple rarefactions but that is outside of the scope of this notebook.

Before we start computation, we're going to get the depth to sample at. The depth is abstracted out for integration testing purposes.

```python
>>> depth = agenv.get_rarefaction_depth()
```

```python
>>> !single_rarefaction.py -i $ag_pgp_hmp_gg_100nt_biom \
...                        -o $ag_pgp_hmp_gg_100nt_1k_biom \
...                        -d $depth
```

Next, we'll compute the distances between samples. By default, this method will produce unweighted and weighted UniFrac distance matrices.

```python
>>> # on test data, single threaded is 15s cputime vs. 15m with 7 threads via parallel_beta_diversity. :(
... !beta_diversity.py -i $ag_pgp_hmp_gg_100nt_1k_biom \
...                    -o $ag_pgp_hmp_gg_100nt_1k_bdiv \
...                    -t $greengenes_tree
```

These next set of cells are going to filter down the distance matrix to interesting subsets. Notably, the AG and GG subset, just the AG samples, the oral subset within AG , the skin subset within AG and the fecal subset within AG.

```python
>>> # unweighted to ag + gg
... !filter_distance_matrix.py -i $ag_pgp_hmp_gg_100nt_1k_bdiv_u \
...                            -o $ag_gg_100nt_1k_bdiv_u \
...                            -t $ag_gg_100nt_biom
```

```python
>>> # unweighted to just ag
... !filter_distance_matrix.py -i $ag_gg_100nt_1k_bdiv_u \
...                            -o $ag_100nt_1k_bdiv_u \
...                            -t $ag_100nt_biom
```

```python
>>> # unweighted to just ag oral
... !filter_distance_matrix.py -i $ag_100nt_1k_bdiv_u \
...                            -o $ag_100nt_oral_1k_bdiv_u \
...                            -m $ag_cleaned_md \
...                            -s "SIMPLE_BODY_SITE:ORAL"
```

```python
>>> # unweighted to just ag skin
... !filter_distance_matrix.py -i $ag_100nt_1k_bdiv_u \
...                            -o $ag_100nt_skin_1k_bdiv_u \
...                            -m $ag_cleaned_md \
...                            -s "SIMPLE_BODY_SITE:SKIN"
```

```python
>>> # unweighted to just ag skin
... !filter_distance_matrix.py -i $ag_100nt_1k_bdiv_u \
...                            -o $ag_100nt_fecal_1k_bdiv_u \
...                            -m $ag_cleaned_md \
...                            -s "SIMPLE_BODY_SITE:FECAL"
```

For one of the figures, the volume of US samples dominates those from other countries. To help mitigate the sample size effect, we're going to subsample the distance matrix prior to producing principal coordinates.

```python
>>> !mod2_pcoa.py subsample_dm --distmat $ag_gg_100nt_1k_bdiv_u \
...                            --max 500 \
...                            --category COUNTRY \
...                            --mapping_file $ag_gg_cleaned_md \
...                            --output $ag_gg_100nt_1k_bdiv_u_ss
```

And finally, we'll produce principal coordinates from the different distance matrices produced above.

```python
>>> !principal_coordinates.py -i $ag_pgp_hmp_gg_100nt_1k_bdiv_u \
...                           -o $ag_pgp_hmp_gg_100nt_1k_u_pc
```

```python
>>> !principal_coordinates.py -i $ag_pgp_hmp_gg_100nt_1k_bdiv_wu \
...                           -o $ag_pgp_hmp_gg_100nt_1k_wu_pc
```

```python
>>> !principal_coordinates.py -i $ag_gg_100nt_1k_bdiv_u \
...                           -o $ag_gg_100nt_1k_u_pc
```

```python
>>> !principal_coordinates.py -i $ag_100nt_1k_bdiv_u \
...                           -o $ag_100nt_1k_u_pc
```

```python
>>> !principal_coordinates.py -i $ag_100nt_oral_1k_bdiv_u \
...                           -o $ag_100nt_oral_1k_u_pc
```

```python
>>> !principal_coordinates.py -i $ag_100nt_skin_1k_bdiv_u \
...                           -o $ag_100nt_skin_1k_u_pc
```

```python
>>> !principal_coordinates.py -i $ag_100nt_fecal_1k_bdiv_u \
...                           -o $ag_100nt_fecal_1k_u_pc
```

```python
>>> !principal_coordinates.py -i $ag_gg_100nt_1k_bdiv_u_ss \
...                           -o $ag_gg_100nt_1k_u_ss_pc
```

To end, as always, let's sanity check and make sure the files of specific interest exist and contain data.

```python
>>> assert os.stat(ag_pgp_hmp_gg_100nt_1k_u_pc).st_size > 0
>>> assert os.stat(ag_pgp_hmp_gg_100nt_1k_wu_pc).st_size > 0
>>> assert os.stat(ag_gg_100nt_1k_u_pc).st_size > 0
>>> assert os.stat(ag_100nt_1k_u_pc).st_size > 0
>>> assert os.stat(ag_100nt_fecal_1k_u_pc ).st_size > 0
>>> assert os.stat(ag_100nt_oral_1k_u_pc ).st_size > 0
>>> assert os.stat(ag_100nt_skin_1k_u_pc).st_size > 0
>>> assert os.stat(ag_gg_100nt_1k_u_ss_pc).st_size > 0
```
