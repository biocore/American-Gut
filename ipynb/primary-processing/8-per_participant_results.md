Now that we've done all the bulk processing, let's generate the per-sample results files that will end up in the PDFs for participants. We have a large number of samples to deal with, and every single sample is independent of each other sample, so we're going to farm out the work over all the processors in the system the notebook is running on.

```python
>>> import os
>>> import multiprocessing as mp
...
>>> import biom
>>> import pandas as pd
>>> from qiime.util import qiime_system_call
...
>>> import americangut.notebook_environment as agenv
>>> import americangut.util as agu
>>> import americangut.results_utils as agru
...
>>> chp_path = agenv.activate('8')
```

First we'll setup our existing paths that we need.

```python
>>> ag_cleaned_md         = agu.get_existing_path(agenv.paths['ag-cleaned-md'])
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
>>> ag_cleaned_df = pd.read_csv(ag_cleaned_md, sep='\t', index_col='#SampleID')
```

These next functions actually perform the processing per site as well as compute support methods.

```python
>>> def _iter_ids_over_system_call(cmd_fmt, sample_ids):
...     """Iteratively execute a system call over sample IDs
...
...     Parameters
...     ----------
...     cmd_fmt : str
...         The format of the command to execute. It is expected that there
...         is a single string format to be done and it should take a sample
...         ID
...     sample_ids : iterable
...         A list of sample IDs of interest
...
...     Returns
...     -------
...     dict
...         A dict containing each sample ID and any errors observed or None if
...         no error was observed for the sample.
...     """
...     results = {}
...
...     for id_ in sample_ids:
...         cmd = cmd_fmt % id_
...         stdout, stderr, return_value = qiime_system_call(cmd)
...
...         if return_value != 0:
...             last_err = stderr.splitlines()[-1]
...             results[id_] = 'FAILED (%s): %s' % (last_err if last_err else '', cmd)
...         else:
...             results[id_] = None
...
...     return results
...
>>> def taxa_summaries(sample_ids, paths):
...     """Produce digestable taxonomy summaries per sample
...
...     Paramters
...     ---------
...     sample_ids : iterable
...         A list of sample IDs of interest
...     paths : dict
...         A dict of relevant paths.
...
...     Returns
...     -------
...     dict
...         A dict containing each sample ID and any errors observed or None if
...         no error was observed for the sample.
...     """
...     results = {}
...     site_table = biom.load_table(paths['site_table'])
...     table_taxon_ids = site_table.ids(axis='observation')
...     taxa_path = os.path.join(paths['output_path'], "Figure_6_%s.txt")
...
...     for id_ in sample_ids:
...         if not site_table.exists(id_):
...             results[id_] = 'ID not found'
...         else:
...             results[id_] = None
...
...             with open(taxa_path % id_, 'w') as fp:
...                 fp.write("#taxon\trelative_abundance\n")
...
...                 v = site_table.data(id_, dense=True)
...                 for sorted_v, taxa in sorted(zip(v, table_taxon_ids))[::-1]:
...                     if sorted_v:
...                         fp.write("%s\t%f\n" % (taxa, sorted_v))
...     return results
...
>>> def taxon_significance(sample_ids, paths):
...     """Produce OTU significance results
...
...     Paramters
...     ---------
...     sample_ids : iterable
...         A list of sample IDs of interest
...     paths : dict
...         A dict of relevant paths.
...
...     Returns
...     -------
...     dict
...         A dict containing each sample ID and any errors observed or None if
...         no error was observed for the sample.
...     """
...     site_table = paths['site_table']
...     output_path = paths['output_path']
...
...     cmd_fmt = 'generate_otu_signifigance_tables_AGP.py -i %s -o %s ' % (site_table, output_path)
...     cmd_fmt += '-s %s'
...
...     return _iter_ids_over_system_call(cmd_fmt, sample_ids)
...
>>> def body_site_pcoa(sample_ids, paths):
...     """Produce the per-sample all body site PCoA
...
...     Paramters
...     ---------
...     sample_ids : iterable
...         A list of sample IDs of interest
...     paths : dict
...         A dict of relevant paths.
...
...     Returns
...     -------
...     dict
...         A dict containing each sample ID and any errors observed or None if
...         no error was observed for the sample.
...     """
...     cmd_fmt = ' '.join(["mod2_pcoa.py body_site",
...                         "--coords %s" % paths['full_pc'],
...                         "--mapping_file %s" % paths['full_md'],
...                         "--output %s" % paths['output_path'],
...                         "--prefix Figure_1",
...                         "--samples %s"])
...
...     return _iter_ids_over_system_call(cmd_fmt, sample_ids)
...
>>> def country_pcoa(sample_ids, paths):
...     """Produce the per-sample subsampled country PCoA
...
...     Paramters
...     ---------
...     sample_ids : iterable
...         A list of sample IDs of interest
...     paths : dict
...         A dict of relevant paths.
...
...     Returns
...     -------
...     dict
...         A dict containing each sample ID and any errors observed or None if
...         no error was observed for the sample.
...     """
...     cmd_fmt = ' '.join(["mod2_pcoa.py country",
...                         "--distmat %s" % paths['ag_gg_dm'],
...                         "--coords %s" % paths['ag_gg_ss_pc'],
...                         "--mapping_file %s" % paths['ag_gg_md'],
...                         "--output %s" % paths['output_path'],
...                         "--prefix Figure_2",
...                         "--samples %s"])
...
...     return _iter_ids_over_system_call(cmd_fmt, sample_ids)
...
>>> def gradient_pcoa(sample_ids, paths):
...     """Produce a gradient PCoA
...
...     Paramters
...     ---------
...     sample_ids : iterable
...         A list of sample IDs of interest
...     paths : dict
...         A dict of relevant paths.
...
...     Returns
...     -------
...     dict
...         A dict containing each sample ID and any errors observed or None if
...         no error was observed for the sample.
...     """
...     cmd_fmt = ' '.join(["mod2_pcoa.py gradient",
...                         "--coords %s" % paths['site_pc'],
...                         "--mapping_file %s" % paths['ag-L2-taxa-md'],
...                         "--output %s" % paths['output_path'],
...                         "--prefix Figure_3",
...                         "--color %s" % paths['color-by'],
...                         "--samples %s"])
...
...     return _iter_ids_over_system_call(cmd_fmt, sample_ids)
...
>>> def pie_plot(sample_ids, paths):
...     """Produce a pie chart
...
...     Paramters
...     ---------
...     sample_ids : iterable
...         A list of sample IDs of interest
...     paths : dict
...         A dict of relevant paths.
...
...     Returns
...     -------
...     dict
...         A dict containing each sample ID and any errors observed or None if
...         no error was observed for the sample.
...     """
...     cmd_fmt = ' '.join(['make_pie_plot_AGP.py',
...                         '-i %s' % paths['ag-L3-taxa-tsv'],
...                         '-o %s' % paths['output_path'],
...                         '-s %s'])
...     return _iter_ids_over_system_call(cmd_fmt, sample_ids)
...
>>> def bar_chart(sample_ids, paths):
...     """Produce a bar chart
...
...     Paramters
...     ---------
...     sample_ids : iterable
...         A list of sample IDs of interest
...     paths : dict
...         A dict of relevant paths.
...
...     Returns
...     -------
...     dict
...         A dict containing each sample ID and any errors observed or None if
...         no error was observed for the sample.
...     """
...     cmd_fmt = ' '.join(['make_phyla_plots_AGP.py',
...                         '-i %s' % paths['site_table_1k'],
...                         '-m %s' % paths['full_md'],
...                         '-o %s' % paths['output_path'],
...                         '-c %s' % paths['barchart_categories'],
...                         '-t %s' % paths['sample_type'],
...                         '-s %s'])
...     return _iter_ids_over_system_call(cmd_fmt, sample_ids)
```

These next methods partition the samples by sample type and enable farming of work over available processors.

```python
>>> def dispatcher(success_fp, fail_fp, partitioner):
...     """Dispatch execution over a pool of processors
...
...     Parameters
...     ----------
...     success_fp : file-like object
...         A file-like object to write a list of successful sample IDs too
...     fail_fp : file-like object
...         A file-like object to write a list of unsuccessful sample IDs too,
...         and any associated error information
...     """
...     pool = mp.Pool(processes=mp.cpu_count())
...
...     success_fp.write('%s\n' % '#SampleID')
...     fail_fp.write('%s\t%s\n' % ('#SampleID', 'Error(s)'))
...
...     for func, ids in partitioner():
...         for success_details in pool.map(func, list(agru.chunk_list(ids))):
...             for id_, detail in success_details.items():
...                 if detail:
...                     fail_fp.write("%s\t%s\n" % (id_, '\t'.join(detail)))
...                 else:
...                     success_fp.write("%s\n" % id_)
...
>>> def partition_samples_by_bodysite():
...     """Yield the processing function and associated sample IDs
...
...     Returns
...     -------
...     generator
...         (fucntion, iterable)
...         The processing function and the sample IDs to examine.
...     """
...     site_and_function = [('FECAL', process_fecal),
...                          ('ORAL', process_oral),
...                          ('SKIN', process_skin)]
...
...     for site, func in site_and_function:
...         df_subset = ag_cleaned_df[ag_cleaned_df.SIMPLE_BODY_SITE == site]
...         yield func, df_subset.index
```

And finally, these next blocks of code support the per-sample type processing.

```python
>>> def merge_error_reports(*reports):
...     """Merge error reports
...
...     Parameters
...     ----------
...     *reports : list of dict
...         The error reports
...
...     Returns
...     -------
...     dict
...         Keyed by sample ID, valued by the list of observed errors. An empty
...         list is associted with a sample ID if no errors were observed.
...     """
...     result = {}
...
...     for report in reports:
...         for id_, value in report.items():
...             if id_ not in result:
...                 result[id_] = []
...
...             if value is not None:
...                 result[id_].append(value)
...
...     return result
...
>>> def get_paths_thread_local(extras):
...     """Grab paths thread local. This method likely is not necessary.
...
...     This method was put in as a precaution against the GIL. It likely
...     isn't necessary and is a refactor target.
...     """
...     import americangut.notebook_environment as agenv_local
...     import americangut.util as agu_local
...
...     paths = {'full_table': agu_local.get_existing_path(agenv_local.paths['ag-L6-taxa-biom']),
...              'full_md': agu_local.get_existing_path(agenv_local.paths['ag-cleaned-md']),
...              'full_dm': agu_local.get_existing_path(agenv_local.paths['ag-pgp-hmp-gg-100nt-1k-bdiv-unifrac']),
...              'full_pc': agu_local.get_existing_path(agenv_local.paths['ag-pgp-hmp-gg-100nt-1k-unifrac-pc']),
...              'ag_gg_dm': agu_local.get_existing_path(agenv_local.paths['ag-gg-100nt-1k-bdiv-unifrac']),
...              'ag_gg_md': agu_local.get_existing_path(agenv_local.paths['ag-gg-cleaned-md']),
...              'ag_gg_ss_pc': agu_local.get_existing_path(agenv_local.paths['ag-gg-100nt-1k-subsampled-unifrac-pc']),
...              'output_path': agu_local.get_existing_path(agenv_local.paths['per-sample-results']),
...              'ag-L2-taxa-biom': agu_local.get_existing_path(agenv_local.paths['ag-L2-taxa-biom']),
...              'ag-L2-taxa-md': agu_local.get_existing_path(agenv_local.paths['ag-L2-taxa-md']),
...              'ag-L3-taxa-tsv': agu.get_existing_path(agenv.paths['ag-L3-taxa-tsv'])
>>> }
...
...     for ksrc, kdest in extras:
...         paths[kdest] = agu_local.get_existing_path(agenv_local.paths[ksrc])
...
...     return paths
...
>>> def process_fecal(ids):
...     """Process per-sample fecal results
...
...     Parameters
...     ----------
...     ids : iterable
...         The list of sample IDs to examine
...
...     Returns
...     -------
...     dict
...         A dict keyed by sample ID and valued by a list. The list contains
...         all errors observed for the sample, or the empty list if no errors
...         were observed.
...     """
...     paths = get_paths_thread_local([('ag-L6-taxa-fecal-biom', 'site_table'),
...                                     ('ag-L2-taxa-fecal-biom', 'site_table_L2'),
...                                     ('ag-100nt-1k-fecal-biom', 'site_table_1k'),
...                                     ('ag-100nt-fecal-1k-unifrac-pc', 'site_pc'),
...                                     ('ag-100nt-1k-fecal-sex-biom', 'sex_biom'),
...                                     ('ag-100nt-1k-fecal-diet-biom', 'diet_biom'),
...                                     ('ag-100nt-1k-fecal-age-biom', 'age_biom'),
...                                     ('ag-100nt-1k-fecal-bmi-biom', 'bmi_biom'),
...                                    ])
...
...     paths['color-by'] = 'k__Bacteria\;p__Firmicutes'
...     paths['sample_type'] = 'fecal'
...     paths['barchart_categories'] = '"%s"' % ', '.join(["DIET_TYPE:%s" % paths['diet_biom'],
...                                                        "BMI_CAT:%s" % paths['bmi_biom'],
...                                                        "SEX:%s" % paths['sex_biom'],
...                                                        "AGE_CAT:%s" % paths['age_biom']])
...     functions = [taxa_summaries,
...                  taxon_significance,
...                  body_site_pcoa,
...                  country_pcoa,
...                  gradient_pcoa,
...                  bar_chart
...                 ]
...
...     return merge_error_reports(*[f(ids, paths) for f in functions])
...
>>> def process_oral(ids):
...     """Process per-sample oral results
...
...     Parameters
...     ----------
...     ids : iterable
...         The list of sample IDs to examine
...
...     Returns
...     -------
...     dict
...         A dict keyed by sample ID and valued by a list. The list contains
...         all errors observed for the sample, or the empty list if no errors
...         were observed.
...     """
...     paths = get_paths_thread_local([('ag-L6-taxa-oral-biom', 'site_table'),
...                                     ('ag-L2-taxa-oral-biom', 'site_table_L2'),
...                                     ('ag-100nt-1k-oral-biom', 'site_table_1k'),
...                                     ('ag-100nt-oral-1k-unifrac-pc', 'site_pc'),
...                                     ('ag-100nt-1k-oral-sex-biom', 'sex_biom'),
...                                     ('ag-100nt-1k-oral-diet-biom', 'diet_biom'),
...                                     ('ag-100nt-1k-oral-age-biom', 'age_biom'),
...                                     ('ag-100nt-1k-oral-flossing-biom', 'floss_biom'),
...
...                                    ])
...     paths['color-by'] = 'k__Bacteria\;p__Firmicutes'
...     paths['sample_type'] = 'oral'
...     paths['barchart_categories'] = '"%s"' % ', '.join(["DIET_TYPE:%s" % paths['diet_biom'],
...                                                        "FLOSSING_FREQUENCY:%s" % paths['floss_biom'],
...                                                        "SEX:%s" % paths['sex_biom'],
...                                                        "AGE_CAT:%s" % paths['age_biom']])
...     functions = [taxa_summaries,
...                  taxon_significance,
...                  body_site_pcoa,
...                  gradient_pcoa,
...                  pie_plot,
...                  bar_chart
...                 ]
...
...     return merge_error_reports(*[f(ids, paths) for f in functions])
...
>>> def process_skin(ids):
...     """Process per-sample skin results
...
...     Parameters
...     ----------
...     ids : iterable
...         The list of sample IDs to examine
...
...     Returns
...     -------
...     dict
...         A dict keyed by sample ID and valued by a list. The list contains
...         all errors observed for the sample, or the empty list if no errors
...         were observed.
...     """
...     paths = get_paths_thread_local([('ag-L6-taxa-skin-biom', 'site_table'),
...                                     ('ag-L2-taxa-skin-biom', 'site_table_L2'),
...                                     ('ag-100nt-1k-skin-biom', 'site_table_1k'),
...                                     ('ag-100nt-skin-1k-unifrac-pc', 'site_pc'),
...                                     ('ag-100nt-1k-skin-sex-biom', 'sex_biom'),
...                                     ('ag-100nt-1k-skin-cosmetics-biom', 'cosmetics_biom'),
...                                     ('ag-100nt-1k-skin-age-biom', 'age_biom'),
...                                     ('ag-100nt-1k-skin-hand-biom', 'hand_biom'),
...                                    ])
...     paths['color-by'] = 'k__Bacteria\;p__Proteobacteria'
...     paths['sample_type'] = 'skin'
...     paths['barchart_categories'] = '"%s"' % ', '.join(["COSMETICS_FREQUENCY:%s" % paths['cosmetics_biom'],
...                                                        "DOMINANT_HAND:%s" % paths['hand_biom'],
...                                                        "SEX:%s" % paths['sex_biom'],
...                                                        "AGE_CAT:%s" % paths['age_biom']])
...     functions = [taxa_summaries,
...                  taxon_significance,
...                  body_site_pcoa,
...                  gradient_pcoa,
...                  pie_plot,
...                  bar_chart
...                 ]
...
...     return merge_error_reports(*[f(ids, paths) for f in functions])
```

And now, let's start mass generating figures!

```python
>>> with open(successful_ids, 'w') as successful_ids_fp, open(unsuccessful_ids, 'w') as unsuccessful_ids_fp:
...     dispatcher(successful_ids_fp, unsuccessful_ids_fp, partition_samples_by_bodysite)
```

And we'll end with some numbers on the number of successful and unsuccessful samples.

```python
>>> print "Number of successfully processed samples: %d" % len([l for l in open(successful_ids) if not l.startswith('#')])
>>> print "Number of unsuccessfully processed samples: %d" % len([l for l in open(unsuccessful_ids) if not l.startswith('#')])
```
