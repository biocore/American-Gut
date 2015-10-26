Now that we've done all the bulk processing, let's generate the per-sample results files that will end up in the PDFs for participants. We have a large number of samples to deal with, and every single sample is independent of each other sample, so we're going to farm out the work over all the processors in the system the notebook is running on.

```python
>>> import os
>>> from functools import partial
>>> import pandas as pd
...
>>> import americangut.notebook_environment as agenv
>>> import americangut.util as agu
>>> import americangut.results_utils as agru
>>> import americangut.per_sample as agps
>>> import americangut.parallel as agpar
...
>>> chp_path = agenv.activate('9')
```

First we'll setup our existing paths that we need.

```python
>>> ag_cleaned_md = agu.get_existing_path(agenv.paths['ag-cleaned-md'])
```

Then we'll establish our new paths as well as "per-sample-results" directory where the individual figures will go.

```python
>>> successful_ids     = agu.get_new_path(agenv.paths['successful-ids'])
>>> unsuccessful_ids   = agu.get_new_path(agenv.paths['unsuccessful-ids'])
>>> per_sample_results = agu.get_new_path(agenv.paths['per-sample-results'])
...
>>> os.mkdir(per_sample_results)
```

We're also going to load up the American Gut mapping file so we can determine what samples (within the 3 major body sites at least) were processed, and what samples had errors.

```python
>>> ag_cleaned_df = pd.read_csv(ag_cleaned_md, sep='\t', dtype=str)
>>> ag_cleaned_df.index = ag_cleaned_df['#SampleID']
```

And finally, these next blocks of code support the per-sample type processing. First, for every sample type, there are common outputs to produce, such as taxonomy summaries. Second, there are some functions that are specific to a sample type. And last, there are a few sample specific options.

```python
>>> common_functions = [agps.per_sample_directory,
...                     agps.stage_per_sample_specific_statics,
...                     agps.taxa_summaries,
...                     agps.taxon_significance,
...                     agps.body_site_pcoa,
...                     agps.gradient_pcoa,
...                     agps.bar_chart]
...
>>> fecal_functions = common_functions + [agps.country_pcoa]
>>> oral_functions  = common_functions + [agps.pie_plot]
>>> skin_functions  = common_functions + [agps.pie_plot]
...
>>> fecal_opts = agps.create_opts('fecal', chp_path,
...                               gradient_color_by='k__Bacteria\;p__Firmicutes',
...                               barchart_categories=('diet', 'bmi', 'sex', 'age'))
>>> oral_opts = agps.create_opts('oral', chp_path,
...                              gradient_color_by='k__Bacteria\;p__Firmicutes',
...                              barchart_categories=('diet', 'flossing', 'sex', 'age'))
>>> skin_opts = agps.create_opts('skin', chp_path,
...                              gradient_color_by='k__Bacteria\;p__Proteobacteria',
...                              barchart_categories=('cosmetics', 'hand', 'sex', 'age'))
...
>>> process_fecal = partial(agps.sample_type_processor, fecal_functions, fecal_opts)
>>> process_oral = partial(agps.sample_type_processor, oral_functions, oral_opts)
>>> process_skin = partial(agps.sample_type_processor, skin_functions, skin_opts)
```

And before the fun starts, let's stage static aspects of the participant results. These are things like the American Gut logo, the result template, etc.

```python
>>> agru.stage_static_files('fecal', chp_path)
>>> agru.stage_static_files('oralskin', chp_path)
```

And now, let's start mass generating figures!

```python
>>> site_to_functions = [('FECAL', process_fecal),
...                      ('ORAL', process_oral),
...                      ('SKIN', process_skin)]
>>> partitions = agps.partition_samples_by_bodysite(ag_cleaned_df, site_to_functions)
...
>>> with open(successful_ids, 'w') as successful_ids_fp, open(unsuccessful_ids, 'w') as unsuccessful_ids_fp:
...     agpar.dispatcher(successful_ids_fp, unsuccessful_ids_fp, partitions)
```

And we'll end with some numbers on the number of successful and unsuccessful samples.

```python
>>> print "Number of successfully processed samples: %d" % len([l for l in open(successful_ids) if not l.startswith('#')])
>>> print "Number of unsuccessfully processed samples: %d" % len([l for l in open(unsuccessful_ids) if not l.startswith('#')])
```
