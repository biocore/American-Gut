In this next chapter, we're going to calculate alpha diversity over all of the samples in the meta-analysis. As always, let's first get our environment setup.

```python
>>> import os
...
>>> import americangut.notebook_environment as agenv
>>> import americangut.util as agu
...
>>> import qiime_default_reference as qdr
...
>>> chp_path = agenv.activate('05')
>>> greengenes_tree = qdr.get_reference_tree()
>>> cpu_count = agenv.get_cpu_count()
```

And then we'll pull in the files we need for processing.

```python
>>> ag_pgp_hmp_gg_100nt_biom = agu.get_existing_path(agenv.paths['ag-pgp-hmp-gg-100nt-biom'])
>>> ag_gg_100nt_biom         = agu.get_existing_path(agenv.paths['ag-gg-100nt-biom'])
>>> ag_pgp_hmp_gg_cleaned_md = agu.get_existing_path(agenv.paths['ag-pgp-hmp-gg-cleaned-md'])
```

As well as establish the new paths that we'll be creating.

```python
>>> ag_pgp_hmp_gg_100nt_1k_multi   = agu.get_new_path(agenv.paths['ag-pgp-hmp-gg-100nt-1k-multiple'])
...
>>> ag_pgp_hmp_gg_100nt_1k_adiv    = agu.get_new_path(agenv.paths['ag-pgp-hmp-gg-100nt-1k-adiv'])
>>> ag_pgp_hmp_gg_100nt_1k_adiv_pd = agu.get_new_path(agenv.paths['ag-pgp-hmp-gg-100nt-1k-adiv-pd'])
>>> ag_pgp_hmp_gg_100nt_1k_adiv_ch = agu.get_new_path(agenv.paths['ag-pgp-hmp-gg-100nt-1k-adiv-chao1'])
>>> ag_pgp_hmp_gg_100nt_1k_adiv_oo = agu.get_new_path(agenv.paths['ag-pgp-hmp-gg-100nt-1k-adiv-observedotus'])
```

Now for the fun bits. We're going to make a call to `parallel_multiple_rarefactions.py` which will produce multiple tables normalized for sequencing effort at 1000 sequences per sample. This process is called rarefaction, and we perform it as the number of sequences produced per sample can vary quite a bit. This variation can greatly influence the resulting diversity scores, so in order to reduce bias driven by sequencing effort, we make an attempt to normalize. The process involves random subsampling, so we're going to do it a few times to get rough bounds on the diversity of each sample.

Before we start computation, we're going to get the depth to sample at. The depth is abstracted out for integration testing purposes.

```python
>>> depth = agenv.get_rarefaction_depth()
```

```python
>>> !multiple_rarefactions_even_depth.py -i $ag_pgp_hmp_gg_100nt_biom \
...                                      -o $ag_pgp_hmp_gg_100nt_1k_multi \
...                                      -d $depth \
```

Once we have the rarefactions, we can then compute the diversity of every sample within our dataset.

```python
>>> !parallel_alpha_diversity.py -i $ag_pgp_hmp_gg_100nt_1k_multi \
...                              -o $ag_pgp_hmp_gg_100nt_1k_adiv \
...                              -t $greengenes_tree \
...                              -O $cpu_count
```

And finally, we're going to aggregate the diversity calculations from the multiple rarefactions.

```python
>>> !collate_alpha.py -i $ag_pgp_hmp_gg_100nt_1k_adiv \
...                   -o $ag_pgp_hmp_gg_100nt_1k_adiv
```

To end, we're going to verify that the files of interest exist and contain data.

```python
>>> assert os.stat(ag_pgp_hmp_gg_100nt_1k_adiv_pd).st_size > 0
>>> assert os.stat(ag_pgp_hmp_gg_100nt_1k_adiv_ch).st_size > 0
>>> assert os.stat(ag_pgp_hmp_gg_100nt_1k_adiv_oo).st_size > 0
```
