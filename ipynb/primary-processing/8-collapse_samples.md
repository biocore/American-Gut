Some of the per-results figures require various slices and perspectives of the data. This notebook performs these necessary additional summarizations.

```python
>>> import os
...
>>> import americangut.notebook_environment as agenv
>>> import americangut.util as agu
...
>>> import qiime_default_reference as qdr
...
>>> chp_path = agenv.activate('8')
```

Let's make sure we have the paths we need.

```python
>>> ag_100nt_biom = agu.get_existing_path(agenv.paths['ag-100nt-biom'])
>>> ag_cleaned_md = agu.get_existing_path(agenv.paths['ag-cleaned-md'])
```

And let's setup all the paths that we're going to create.

```python
>>> ag_100nt_1k_biom       = agu.get_new_path(agenv.paths['ag-100nt-1k-biom'])
>>> ag_100nt_1k_fecal_biom = agu.get_new_path(agenv.paths['ag-100nt-1k-fecal-biom'])
>>> ag_100nt_1k_skin_biom  = agu.get_new_path(agenv.paths['ag-100nt-1k-skin-biom'])
>>> ag_100nt_1k_oral_biom  = agu.get_new_path(agenv.paths['ag-100nt-1k-oral-biom'])
...
>>> ag_100nt_1k_fecal_sex_biom  = agu.get_new_path(agenv.paths['ag-100nt-1k-fecal-sex-biom'])
>>> ag_100nt_1k_fecal_diet_biom = agu.get_new_path(agenv.paths['ag-100nt-1k-fecal-diet-biom'])
>>> ag_100nt_1k_fecal_age_biom  = agu.get_new_path(agenv.paths['ag-100nt-1k-fecal-age-biom'])
>>> ag_100nt_1k_fecal_bmi_biom  = agu.get_new_path(agenv.paths['ag-100nt-1k-fecal-bmi-biom'])
...
>>> ag_100nt_1k_oral_sex_biom       = agu.get_new_path(agenv.paths['ag-100nt-1k-oral-sex-biom'])
>>> ag_100nt_1k_oral_diet_biom      = agu.get_new_path(agenv.paths['ag-100nt-1k-oral-diet-biom'])
>>> ag_100nt_1k_oral_age_biom       = agu.get_new_path(agenv.paths['ag-100nt-1k-oral-age-biom'])
>>> ag_100nt_1k_oral_flossing_biom  = agu.get_new_path(agenv.paths['ag-100nt-1k-oral-flossing-biom'])
...
>>> ag_100nt_1k_skin_sex_biom       = agu.get_new_path(agenv.paths['ag-100nt-1k-skin-sex-biom'])
>>> ag_100nt_1k_skin_cosmetics_biom = agu.get_new_path(agenv.paths['ag-100nt-1k-skin-cosmetics-biom'])
>>> ag_100nt_1k_skin_age_biom       = agu.get_new_path(agenv.paths['ag-100nt-1k-skin-age-biom'])
>>> ag_100nt_1k_skin_hand_biom      = agu.get_new_path(agenv.paths['ag-100nt-1k-skin-hand-biom'])
...
>>> ignored = 'foo'
```

First, we're going to operate on rarefied data again.

```python
>>> depth = agenv.get_rarefaction_depth()
```

```python
>>> !single_rarefaction.py -i $ag_100nt_biom \
...                        -o $ag_100nt_1k_biom \
...                        -d $depth
```

Next, we're going to partition the data into per-body site tables.

```python
>>> !filter_samples_from_otu_table.py -i $ag_100nt_1k_biom \
...                                   -o $ag_100nt_1k_fecal_biom \
...                                   -m $ag_cleaned_md \
...                                   -s "SIMPLE_BODY_SITE:FECAL"
```

```python
>>> !filter_samples_from_otu_table.py -i $ag_100nt_1k_biom \
...                                   -o $ag_100nt_1k_skin_biom \
...                                   -m $ag_cleaned_md \
...                                   -s "SIMPLE_BODY_SITE:SKIN"
```

```python
>>> !filter_samples_from_otu_table.py -i $ag_100nt_1k_biom \
...                                   -o $ag_100nt_1k_oral_biom \
...                                   -m $ag_cleaned_md \
...                                   -s "SIMPLE_BODY_SITE:ORAL"
```

Finally, within each body site, we're going to collapse over categories of interest.

```python
>>> !collapse_samples.py -m $ag_cleaned_md \
...                      --output_biom_fp $ag_100nt_1k_fecal_sex_biom \
...                      --normalize \
...                      -b $ag_100nt_1k_fecal_biom \
...                      --collapse_fields "SEX" \
...                      --output_mapping_fp $ignored
```

```python
>>> !collapse_samples.py -m $ag_cleaned_md \
...                      --output_biom_fp $ag_100nt_1k_fecal_diet_biom \
...                      --normalize \
...                      -b $ag_100nt_1k_fecal_biom \
...                      --collapse_fields "DIET_TYPE" \
...                      --output_mapping_fp $ignored
```

```python
>>> !collapse_samples.py -m $ag_cleaned_md \
...                      --output_biom_fp $ag_100nt_1k_fecal_age_biom \
...                      --normalize \
...                      -b $ag_100nt_1k_fecal_biom \
...                      --collapse_fields "AGE_CAT" \
...                      --output_mapping_fp $ignored
```

```python
>>> !collapse_samples.py -m $ag_cleaned_md \
...                      --output_biom_fp $ag_100nt_1k_fecal_bmi_biom \
...                      --normalize \
...                      -b $ag_100nt_1k_fecal_biom \
...                      --collapse_fields "BMI_CAT" \
...                      --output_mapping_fp $ignored
```

```python
>>> !collapse_samples.py -m $ag_cleaned_md \
...                      --output_biom_fp $ag_100nt_1k_oral_sex_biom \
...                      --normalize \
...                      -b $ag_100nt_1k_oral_biom \
...                      --collapse_fields "SEX" \
...                      --output_mapping_fp $ignored
```

```python
>>> !collapse_samples.py -m $ag_cleaned_md \
...                      --output_biom_fp $ag_100nt_1k_oral_diet_biom \
...                      --normalize \
...                      -b $ag_100nt_1k_oral_biom \
...                      --collapse_fields "DIET_TYPE" \
...                      --output_mapping_fp $ignored
```

```python
>>> !collapse_samples.py -m $ag_cleaned_md \
...                      --output_biom_fp $ag_100nt_1k_oral_age_biom \
...                      --normalize \
...                      -b $ag_100nt_1k_oral_biom \
...                      --collapse_fields "AGE_CAT" \
...                      --output_mapping_fp $ignored
```

```python
>>> !collapse_samples.py -m $ag_cleaned_md \
...                      --output_biom_fp $ag_100nt_1k_oral_flossing_biom \
...                      --normalize \
...                      -b $ag_100nt_1k_oral_biom \
...                      --collapse_fields "FLOSSING_FREQUENCY" \
...                      --output_mapping_fp $ignored
```

```python
>>> !collapse_samples.py -m $ag_cleaned_md \
...                      --output_biom_fp $ag_100nt_1k_skin_sex_biom \
...                      --normalize \
...                      -b $ag_100nt_1k_skin_biom \
...                      --collapse_fields "SEX" \
...                      --output_mapping_fp $ignored
```

```python
>>> !collapse_samples.py -m $ag_cleaned_md \
...                      --output_biom_fp $ag_100nt_1k_skin_cosmetics_biom \
...                      --normalize \
...                      -b $ag_100nt_1k_skin_biom \
...                      --collapse_fields "COSMETICS_FREQUENCY" \
...                      --output_mapping_fp $ignored
```

```python
>>> !collapse_samples.py -m $ag_cleaned_md \
...                      --output_biom_fp $ag_100nt_1k_skin_age_biom \
...                      --normalize \
...                      -b $ag_100nt_1k_skin_biom \
...                      --collapse_fields "AGE_CAT" \
...                      --output_mapping_fp $ignored
```

```python
>>> !collapse_samples.py -m $ag_cleaned_md \
...                      --output_biom_fp $ag_100nt_1k_skin_hand_biom \
...                      --normalize \
...                      -b $ag_100nt_1k_skin_biom \
...                      --collapse_fields "DOMINANT_HAND" \
...                      --output_mapping_fp $ignored
```

As usual, let's make sure we have files.

```python
>>> assert os.stat(ag_100nt_1k_fecal_biom).st_size > 0
>>> assert os.stat(ag_100nt_1k_skin_biom).st_size > 0
>>> assert os.stat(ag_100nt_1k_oral_biom).st_size > 0
>>> assert os.stat(ag_100nt_1k_fecal_sex_biom).st_size > 0
>>> assert os.stat(ag_100nt_1k_fecal_diet_biom).st_size > 0
>>> assert os.stat(ag_100nt_1k_fecal_age_biom).st_size > 0
>>> assert os.stat(ag_100nt_1k_fecal_bmi_biom).st_size > 0
>>> assert os.stat(ag_100nt_1k_oral_sex_biom).st_size > 0
>>> assert os.stat(ag_100nt_1k_oral_diet_biom).st_size > 0
>>> assert os.stat(ag_100nt_1k_oral_age_biom).st_size > 0
>>> assert os.stat(ag_100nt_1k_oral_flossing_biom).st_size > 0
>>> assert os.stat(ag_100nt_1k_skin_sex_biom).st_size > 0
>>> assert os.stat(ag_100nt_1k_skin_hand_biom).st_size > 0
>>> assert os.stat(ag_100nt_1k_skin_age_biom).st_size > 0
>>> assert os.stat(ag_100nt_1k_skin_cosmetics_biom).st_size > 0
```
