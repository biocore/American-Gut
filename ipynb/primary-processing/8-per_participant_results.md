Now that we've done all the bulk processing, let's generate the per-sample results files that will end up in the PDFs for participants. We have a large number of samples to deal with, and every single sample is independent of each other sample, so we're going to farm out the work over all the processors in the system the notebook is running on.

```python
>>> import os
>>> import shutil
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
>>> ag_cleaned_df = pd.read_csv(ag_cleaned_md, sep='\t', index_col='#SampleID')
```

These next functions actually perform the processing per site as well as compute support methods.

```python
>>> def _result_path(paths, id_):
...     return os.path.join(paths['output_path'], id_)
...
>>> def _base_barcode(id_):
...     ### NOTE: old ebi accessions and processing do not follow this format.
...     # expectation is <study_id>.<barcode>
...
...     return id_.split('.')[1]
...
>>> def _iter_ids_over_system_call(cmd_fmt, sample_ids, paths):
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
...         cmd = cmd_fmt % {'result_path': _result_path(paths, id_),
...                          'id': id_}
...         print cmd
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
...
...     for id_ in sample_ids:
...         if not site_table.exists(id_):
...             results[id_] = 'ID not found'
...         else:
...             results[id_] = None
...             taxa_path = os.path.join(_result_path(paths, id_), '%s.txt')
...
...             with open(taxa_path % _base_barcode(id_), 'w') as fp:
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
...
...     cmd_fmt = 'generate_otu_signifigance_tables_AGP.py -i %s ' % site_table
...     cmd_fmt += '-o %(result_path)s -s %(id)s'
...
...     return _iter_ids_over_system_call(cmd_fmt, sample_ids, paths)
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
...                         "--output %(result_path)s",
...                         "--prefix Figure_1",
...                         "--samples %(id)s"])
...
...     return _iter_ids_over_system_call(cmd_fmt, sample_ids, paths)
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
...                         "--output %(result_path)s",
...                         "--prefix Figure_2",
...                         "--samples %(id)s"])
...
...     return _iter_ids_over_system_call(cmd_fmt, sample_ids, paths)
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
...                         "--output %(result_path)s",
...                         "--prefix Figure_3",
...                         "--color %s" % paths['color-by'],
...                         "--samples %(id)s"])
...
...     return _iter_ids_over_system_call(cmd_fmt, sample_ids, paths)
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
...                         '-o %(result_path)s',
...                         '-s %(id)s'])
...     return _iter_ids_over_system_call(cmd_fmt, sample_ids, paths)
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
...                         '-o %(result_path)s',
...                         '-c %s' % paths['barchart_categories'],
...                         '-t %s' % paths['sample_type'],
...                         '-s %(id)s'])
...     return _iter_ids_over_system_call(cmd_fmt, sample_ids, paths)
...
>>> def per_sample_directory(ids, paths):
...     """Create a per sample directory
...
...     Paramters
...     ---------
...     sample_ids : iterable
...         A list of sample IDs of interest
...     paths : dict
...         A dict of relevant paths.
...     """
...     for id_ in ids:
...         path = _result_path(paths, id_)
...         if not os.path.exists(path):
...             os.mkdir(path)
...     return {}
...
>>> def stage_per_sample_specific_statics(ids, paths):
...     """Items like the Latex templates
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
...     for id_ in ids:
...         path = _result_path(paths, id_)
...
...         if paths['sample_type'] == 'fecal':
...             shutil.copy(os.path.join(chp_path, 'template_gut.tex'), path)
...         elif paths['sample_type'] in ('oral', 'skin'):
...             shutil.copy(os.path.join(chp_path, 'template_oralskin.tex'), path)
...         else:
...             raise ValueError('Unknown sample type: %s' % paths['sample_type'])
...     return {}
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
...              'ag-L3-taxa-tsv': agu.get_existing_path(agenv.paths['ag-L3-taxa-tsv'])}
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
...
...     functions = [per_sample_directory,
...                  stage_per_sample_specific_statics,
...                  taxa_summaries,
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
...     functions = [per_sample_directory,
...                  stage_per_sample_specific_statics,
...                  taxa_summaries,
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
...     functions = [per_sample_directory,
...                  stage_per_sample_specific_statics,
...                  taxa_summaries,
...                  taxon_significance,
...                  body_site_pcoa,
...                  gradient_pcoa,
...                  pie_plot,
...                  bar_chart
...                 ]
...
...     return merge_error_reports(*[f(ids, paths) for f in functions])
```

And before the fun starts, let's stage static aspects of the participant results.

```python
>>> agru.stage_static_files('fecal', chp_path)
>>> agru.stage_static_files('oralskin', chp_path)
```

And now, let's start mass generating figures!

```python
>>> with open(successful_ids, 'w') as successful_ids_fp, open(unsuccessful_ids, 'w') as unsuccessful_ids_fp:
...     dispatcher(successful_ids_fp, unsuccessful_ids_fp, partition_samples_by_bodysite)
generate_otu_signifigance_tables_AGP.py -i agp_processing/7/taxa/otu_table_fecal_L6.biom -o agp_processing/8/per-sample-results/000007108.1075657 -s 000007108.1075657
generate_otu_signifigance_tables_AGP.py -i agp_processing/7/taxa/otu_table_oral_L6.biom -o agp_processing/8/per-sample-results/000013396.fixed997 -s 000013396.fixed997
generate_otu_signifigance_tables_AGP.py -i agp_processing/7/taxa/otu_table_skin_L6.biom -o agp_processing/8/per-sample-results/000007113.1075702 -s 000007113.1075702
generate_otu_signifigance_tables_AGP.py -i agp_processing/7/taxa/otu_table_fecal_L6.biom -o agp_processing/8/per-sample-results/000002035.1075856 -s 000002035.1075856generate_otu_signifigance_tables_AGP.py -i agp_processing/7/taxa/otu_table_oral_L6.biom -o agp_processing/8/per-sample-results/000002207.1210605 -s 000002207.1210605generate_otu_signifigance_tables_AGP.py -i agp_processing/7/taxa/otu_table_skin_L6.biom -o agp_processing/8/per-sample-results/000005644.1053889 -s 000005644.1053889


generate_otu_signifigance_tables_AGP.py -i agp_processing/7/taxa/otu_table_fecal_L6.biom -o agp_processing/8/per-sample-results/000002244.1076101 -s 000002244.1076101generate_otu_signifigance_tables_AGP.py -i agp_processing/7/taxa/otu_table_oral_L6.biom -o agp_processing/8/per-sample-results/000005362.1131922 -s 000005362.1131922generate_otu_signifigance_tables_AGP.py -i agp_processing/7/taxa/otu_table_skin_L6.biom -o agp_processing/8/per-sample-results/000007776.fixed893 -s 000007776.fixed893


generate_otu_signifigance_tables_AGP.py -i agp_processing/7/taxa/otu_table_fecal_L6.biom -o agp_processing/8/per-sample-results/000005844.1131797 -s 000005844.1131797generate_otu_signifigance_tables_AGP.py -i agp_processing/7/taxa/otu_table_oral_L6.biom -o agp_processing/8/per-sample-results/000007109.1075688 -s 000007109.1075688generate_otu_signifigance_tables_AGP.py -i agp_processing/7/taxa/otu_table_skin_L6.biom -o agp_processing/8/per-sample-results/000007824.fixed1092 -s 000007824.fixed1092


generate_otu_signifigance_tables_AGP.py -i agp_processing/7/taxa/otu_table_fecal_L6.biom -o agp_processing/8/per-sample-results/000009111.1131950 -s 000009111.1131950mod2_pcoa.py body_site --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-pgp-hmp-gg-100nt-1k-pc.txt --mapping_file agp_processing/4/ag-cleaned.txt --output agp_processing/8/per-sample-results/000013396.fixed997 --prefix Figure_1 --samples 000013396.fixed997generate_otu_signifigance_tables_AGP.py -i agp_processing/7/taxa/otu_table_skin_L6.biom -o agp_processing/8/per-sample-results/000009079.1130049 -s 000009079.1130049


generate_otu_signifigance_tables_AGP.py -i agp_processing/7/taxa/otu_table_fecal_L6.biom -o agp_processing/8/per-sample-results/000007118.1075682 -s 000007118.1075682mod2_pcoa.py body_site --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-pgp-hmp-gg-100nt-1k-pc.txt --mapping_file agp_processing/4/ag-cleaned.txt --output agp_processing/8/per-sample-results/000002207.1210605 --prefix Figure_1 --samples 000002207.1210605mod2_pcoa.py body_site --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-pgp-hmp-gg-100nt-1k-pc.txt --mapping_file agp_processing/4/ag-cleaned.txt --output agp_processing/8/per-sample-results/000007113.1075702 --prefix Figure_1 --samples 000007113.1075702


mod2_pcoa.py body_site --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-pgp-hmp-gg-100nt-1k-pc.txt --mapping_file agp_processing/4/ag-cleaned.txt --output agp_processing/8/per-sample-results/000007108.1075657 --prefix Figure_1 --samples 000007108.1075657mod2_pcoa.py body_site --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-pgp-hmp-gg-100nt-1k-pc.txt --mapping_file agp_processing/4/ag-cleaned.txt --output agp_processing/8/per-sample-results/000005362.1131922 --prefix Figure_1 --samples 000005362.1131922mod2_pcoa.py body_site --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-pgp-hmp-gg-100nt-1k-pc.txt --mapping_file agp_processing/4/ag-cleaned.txt --output agp_processing/8/per-sample-results/000005644.1053889 --prefix Figure_1 --samples 000005644.1053889


mod2_pcoa.py body_site --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-pgp-hmp-gg-100nt-1k-pc.txt --mapping_file agp_processing/4/ag-cleaned.txt --output agp_processing/8/per-sample-results/000002035.1075856 --prefix Figure_1 --samples 000002035.1075856mod2_pcoa.py body_site --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-pgp-hmp-gg-100nt-1k-pc.txt --mapping_file agp_processing/4/ag-cleaned.txt --output agp_processing/8/per-sample-results/000007109.1075688 --prefix Figure_1 --samples 000007109.1075688mod2_pcoa.py body_site --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-pgp-hmp-gg-100nt-1k-pc.txt --mapping_file agp_processing/4/ag-cleaned.txt --output agp_processing/8/per-sample-results/000007776.fixed893 --prefix Figure_1 --samples 000007776.fixed893


mod2_pcoa.py body_site --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-pgp-hmp-gg-100nt-1k-pc.txt --mapping_file agp_processing/4/ag-cleaned.txt --output agp_processing/8/per-sample-results/000002244.1076101 --prefix Figure_1 --samples 000002244.1076101mod2_pcoa.py gradient --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-100nt-oral-1k-pc.txt --mapping_file agp_processing/7/taxa/ag-cleaned_L2.txt --output agp_processing/8/per-sample-results/000013396.fixed997 --prefix Figure_3 --color k__Bacteria\;p__Firmicutes --samples 000013396.fixed997mod2_pcoa.py body_site --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-pgp-hmp-gg-100nt-1k-pc.txt --mapping_file agp_processing/4/ag-cleaned.txt --output agp_processing/8/per-sample-results/000007824.fixed1092 --prefix Figure_1 --samples 000007824.fixed1092


mod2_pcoa.py body_site --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-pgp-hmp-gg-100nt-1k-pc.txt --mapping_file agp_processing/4/ag-cleaned.txt --output agp_processing/8/per-sample-results/000005844.1131797 --prefix Figure_1 --samples 000005844.1131797mod2_pcoa.py gradient --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-100nt-oral-1k-pc.txt --mapping_file agp_processing/7/taxa/ag-cleaned_L2.txt --output agp_processing/8/per-sample-results/000002207.1210605 --prefix Figure_3 --color k__Bacteria\;p__Firmicutes --samples 000002207.1210605mod2_pcoa.py body_site --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-pgp-hmp-gg-100nt-1k-pc.txt --mapping_file agp_processing/4/ag-cleaned.txt --output agp_processing/8/per-sample-results/000009079.1130049 --prefix Figure_1 --samples 000009079.1130049


mod2_pcoa.py body_site --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-pgp-hmp-gg-100nt-1k-pc.txt --mapping_file agp_processing/4/ag-cleaned.txt --output agp_processing/8/per-sample-results/000009111.1131950 --prefix Figure_1 --samples 000009111.1131950mod2_pcoa.py gradient --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-100nt-oral-1k-pc.txt --mapping_file agp_processing/7/taxa/ag-cleaned_L2.txt --output agp_processing/8/per-sample-results/000005362.1131922 --prefix Figure_3 --color k__Bacteria\;p__Firmicutes --samples 000005362.1131922mod2_pcoa.py gradient --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-100nt-skin-1k-pc.txt --mapping_file agp_processing/7/taxa/ag-cleaned_L2.txt --output agp_processing/8/per-sample-results/000007113.1075702 --prefix Figure_3 --color k__Bacteria\;p__Proteobacteria --samples 000007113.1075702


mod2_pcoa.py body_site --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-pgp-hmp-gg-100nt-1k-pc.txt --mapping_file agp_processing/4/ag-cleaned.txt --output agp_processing/8/per-sample-results/000007118.1075682 --prefix Figure_1 --samples 000007118.1075682mod2_pcoa.py gradient --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-100nt-oral-1k-pc.txt --mapping_file agp_processing/7/taxa/ag-cleaned_L2.txt --output agp_processing/8/per-sample-results/000007109.1075688 --prefix Figure_3 --color k__Bacteria\;p__Firmicutes --samples 000007109.1075688mod2_pcoa.py gradient --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-100nt-skin-1k-pc.txt --mapping_file agp_processing/7/taxa/ag-cleaned_L2.txt --output agp_processing/8/per-sample-results/000005644.1053889 --prefix Figure_3 --color k__Bacteria\;p__Proteobacteria --samples 000005644.1053889


mod2_pcoa.py country --distmat agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-gg-100nt-1k.txt --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-gg-100nt-1k-subsampled-pc.txt --mapping_file agp_processing/4/ag-gg-cleaned.txt --output agp_processing/8/per-sample-results/000007108.1075657 --prefix Figure_2 --samples 000007108.1075657make_pie_plot_AGP.py -i agp_processing/7/taxa/otu_table_L3.txt -o agp_processing/8/per-sample-results/000013396.fixed997 -s 000013396.fixed997mod2_pcoa.py gradient --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-100nt-skin-1k-pc.txt --mapping_file agp_processing/7/taxa/ag-cleaned_L2.txt --output agp_processing/8/per-sample-results/000007776.fixed893 --prefix Figure_3 --color k__Bacteria\;p__Proteobacteria --samples 000007776.fixed893


mod2_pcoa.py country --distmat agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-gg-100nt-1k.txt --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-gg-100nt-1k-subsampled-pc.txt --mapping_file agp_processing/4/ag-gg-cleaned.txt --output agp_processing/8/per-sample-results/000002035.1075856 --prefix Figure_2 --samples 000002035.1075856make_pie_plot_AGP.py -i agp_processing/7/taxa/otu_table_L3.txt -o agp_processing/8/per-sample-results/000002207.1210605 -s 000002207.1210605mod2_pcoa.py gradient --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-100nt-skin-1k-pc.txt --mapping_file agp_processing/7/taxa/ag-cleaned_L2.txt --output agp_processing/8/per-sample-results/000007824.fixed1092 --prefix Figure_3 --color k__Bacteria\;p__Proteobacteria --samples 000007824.fixed1092


mod2_pcoa.py country --distmat agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-gg-100nt-1k.txt --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-gg-100nt-1k-subsampled-pc.txt --mapping_file agp_processing/4/ag-gg-cleaned.txt --output agp_processing/8/per-sample-results/000002244.1076101 --prefix Figure_2 --samples 000002244.1076101make_pie_plot_AGP.py -i agp_processing/7/taxa/otu_table_L3.txt -o agp_processing/8/per-sample-results/000005362.1131922 -s 000005362.1131922mod2_pcoa.py gradient --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-100nt-skin-1k-pc.txt --mapping_file agp_processing/7/taxa/ag-cleaned_L2.txt --output agp_processing/8/per-sample-results/000009079.1130049 --prefix Figure_3 --color k__Bacteria\;p__Proteobacteria --samples 000009079.1130049


mod2_pcoa.py country --distmat agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-gg-100nt-1k.txt --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-gg-100nt-1k-subsampled-pc.txt --mapping_file agp_processing/4/ag-gg-cleaned.txt --output agp_processing/8/per-sample-results/000005844.1131797 --prefix Figure_2 --samples 000005844.1131797make_pie_plot_AGP.py -i agp_processing/7/taxa/otu_table_L3.txt -o agp_processing/8/per-sample-results/000007109.1075688 -s 000007109.1075688make_pie_plot_AGP.py -i agp_processing/7/taxa/otu_table_L3.txt -o agp_processing/8/per-sample-results/000007113.1075702 -s 000007113.1075702


mod2_pcoa.py country --distmat agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-gg-100nt-1k.txt --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-gg-100nt-1k-subsampled-pc.txt --mapping_file agp_processing/4/ag-gg-cleaned.txt --output agp_processing/8/per-sample-results/000009111.1131950 --prefix Figure_2 --samples 000009111.1131950make_phyla_plots_AGP.py -i agp_processing/7.5/ag-100nt-1k-skin.biom -m agp_processing/4/ag-cleaned.txt -o agp_processing/8/per-sample-results/000013396.fixed997 -c "DIET_TYPE:agp_processing/7.5/ag-100nt-1k-oral-diet.biom, FLOSSING_FREQUENCY:agp_processing/7.5/ag-100nt-1k-oral-flossing.biom, SEX:agp_processing/7.5/ag-100nt-1k-oral-sex.biom, AGE_CAT:agp_processing/7.5/ag-100nt-1k-oral-age.biom" -t oral -s 000013396.fixed997make_pie_plot_AGP.py -i agp_processing/7/taxa/otu_table_L3.txt -o agp_processing/8/per-sample-results/000005644.1053889 -s 000005644.1053889


mod2_pcoa.py country --distmat agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-gg-100nt-1k.txt --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-gg-100nt-1k-subsampled-pc.txt --mapping_file agp_processing/4/ag-gg-cleaned.txt --output agp_processing/8/per-sample-results/000007118.1075682 --prefix Figure_2 --samples 000007118.1075682make_phyla_plots_AGP.py -i agp_processing/7.5/ag-100nt-1k-skin.biom -m agp_processing/4/ag-cleaned.txt -o agp_processing/8/per-sample-results/000002207.1210605 -c "DIET_TYPE:agp_processing/7.5/ag-100nt-1k-oral-diet.biom, FLOSSING_FREQUENCY:agp_processing/7.5/ag-100nt-1k-oral-flossing.biom, SEX:agp_processing/7.5/ag-100nt-1k-oral-sex.biom, AGE_CAT:agp_processing/7.5/ag-100nt-1k-oral-age.biom" -t oral -s 000002207.1210605make_pie_plot_AGP.py -i agp_processing/7/taxa/otu_table_L3.txt -o agp_processing/8/per-sample-results/000007776.fixed893 -s 000007776.fixed893


mod2_pcoa.py gradient --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-100nt-fecal-1k-pc.txt --mapping_file agp_processing/7/taxa/ag-cleaned_L2.txt --output agp_processing/8/per-sample-results/000007108.1075657 --prefix Figure_3 --color k__Bacteria\;p__Firmicutes --samples 000007108.1075657make_phyla_plots_AGP.py -i agp_processing/7.5/ag-100nt-1k-skin.biom -m agp_processing/4/ag-cleaned.txt -o agp_processing/8/per-sample-results/000005362.1131922 -c "DIET_TYPE:agp_processing/7.5/ag-100nt-1k-oral-diet.biom, FLOSSING_FREQUENCY:agp_processing/7.5/ag-100nt-1k-oral-flossing.biom, SEX:agp_processing/7.5/ag-100nt-1k-oral-sex.biom, AGE_CAT:agp_processing/7.5/ag-100nt-1k-oral-age.biom" -t oral -s 000005362.1131922make_pie_plot_AGP.py -i agp_processing/7/taxa/otu_table_L3.txt -o agp_processing/8/per-sample-results/000007824.fixed1092 -s 000007824.fixed1092


mod2_pcoa.py gradient --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-100nt-fecal-1k-pc.txt --mapping_file agp_processing/7/taxa/ag-cleaned_L2.txt --output agp_processing/8/per-sample-results/000002035.1075856 --prefix Figure_3 --color k__Bacteria\;p__Firmicutes --samples 000002035.1075856make_phyla_plots_AGP.py -i agp_processing/7.5/ag-100nt-1k-skin.biom -m agp_processing/4/ag-cleaned.txt -o agp_processing/8/per-sample-results/000007109.1075688 -c "DIET_TYPE:agp_processing/7.5/ag-100nt-1k-oral-diet.biom, FLOSSING_FREQUENCY:agp_processing/7.5/ag-100nt-1k-oral-flossing.biom, SEX:agp_processing/7.5/ag-100nt-1k-oral-sex.biom, AGE_CAT:agp_processing/7.5/ag-100nt-1k-oral-age.biom" -t oral -s 000007109.1075688make_pie_plot_AGP.py -i agp_processing/7/taxa/otu_table_L3.txt -o agp_processing/8/per-sample-results/000009079.1130049 -s 000009079.1130049


mod2_pcoa.py gradient --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-100nt-fecal-1k-pc.txt --mapping_file agp_processing/7/taxa/ag-cleaned_L2.txt --output agp_processing/8/per-sample-results/000002244.1076101 --prefix Figure_3 --color k__Bacteria\;p__Firmicutes --samples 000002244.1076101make_phyla_plots_AGP.py -i agp_processing/7.5/ag-100nt-1k-oral.biom -m agp_processing/4/ag-cleaned.txt -o agp_processing/8/per-sample-results/000007113.1075702 -c "COSMETICS_FREQUENCY:agp_processing/7.5/ag-100nt-1k-skin-cosmetics.biom, DOMINANT_HAND:agp_processing/7.5/ag-100nt-1k-skin-hand.biom, SEX:agp_processing/7.5/ag-100nt-1k-skin-sex.biom, AGE_CAT:agp_processing/7.5/ag-100nt-1k-skin-age.biom" -t skin -s 000007113.1075702

mod2_pcoa.py gradient --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-100nt-fecal-1k-pc.txt --mapping_file agp_processing/7/taxa/ag-cleaned_L2.txt --output agp_processing/8/per-sample-results/000005844.1131797 --prefix Figure_3 --color k__Bacteria\;p__Firmicutes --samples 000005844.1131797make_phyla_plots_AGP.py -i agp_processing/7.5/ag-100nt-1k-oral.biom -m agp_processing/4/ag-cleaned.txt -o agp_processing/8/per-sample-results/000005644.1053889 -c "COSMETICS_FREQUENCY:agp_processing/7.5/ag-100nt-1k-skin-cosmetics.biom, DOMINANT_HAND:agp_processing/7.5/ag-100nt-1k-skin-hand.biom, SEX:agp_processing/7.5/ag-100nt-1k-skin-sex.biom, AGE_CAT:agp_processing/7.5/ag-100nt-1k-skin-age.biom" -t skin -s 000005644.1053889

mod2_pcoa.py gradient --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-100nt-fecal-1k-pc.txt --mapping_file agp_processing/7/taxa/ag-cleaned_L2.txt --output agp_processing/8/per-sample-results/000009111.1131950 --prefix Figure_3 --color k__Bacteria\;p__Firmicutes --samples 000009111.1131950make_phyla_plots_AGP.py -i agp_processing/7.5/ag-100nt-1k-oral.biom -m agp_processing/4/ag-cleaned.txt -o agp_processing/8/per-sample-results/000007776.fixed893 -c "COSMETICS_FREQUENCY:agp_processing/7.5/ag-100nt-1k-skin-cosmetics.biom, DOMINANT_HAND:agp_processing/7.5/ag-100nt-1k-skin-hand.biom, SEX:agp_processing/7.5/ag-100nt-1k-skin-sex.biom, AGE_CAT:agp_processing/7.5/ag-100nt-1k-skin-age.biom" -t skin -s 000007776.fixed893

mod2_pcoa.py gradient --coords agp_processing/6/ag-pgp-hmp-gg-100nt-1k-bdiv/unweighted_unifrac_ag-100nt-fecal-1k-pc.txt --mapping_file agp_processing/7/taxa/ag-cleaned_L2.txt --output agp_processing/8/per-sample-results/000007118.1075682 --prefix Figure_3 --color k__Bacteria\;p__Firmicutes --samples 000007118.1075682make_phyla_plots_AGP.py -i agp_processing/7.5/ag-100nt-1k-oral.biom -m agp_processing/4/ag-cleaned.txt -o agp_processing/8/per-sample-results/000007824.fixed1092 -c "COSMETICS_FREQUENCY:agp_processing/7.5/ag-100nt-1k-skin-cosmetics.biom, DOMINANT_HAND:agp_processing/7.5/ag-100nt-1k-skin-hand.biom, SEX:agp_processing/7.5/ag-100nt-1k-skin-sex.biom, AGE_CAT:agp_processing/7.5/ag-100nt-1k-skin-age.biom" -t skin -s 000007824.fixed1092

make_phyla_plots_AGP.py -i agp_processing/7.5/ag-100nt-1k-fecal.biom -m agp_processing/4/ag-cleaned.txt -o agp_processing/8/per-sample-results/000007108.1075657 -c "DIET_TYPE:agp_processing/7.5/ag-100nt-1k-fecal-diet.biom, BMI_CAT:agp_processing/7.5/ag-100nt-1k-fecal-bmi.biom, SEX:agp_processing/7.5/ag-100nt-1k-fecal-sex.biom, AGE_CAT:agp_processing/7.5/ag-100nt-1k-fecal-age.biom" -t fecal -s 000007108.1075657make_phyla_plots_AGP.py -i agp_processing/7.5/ag-100nt-1k-oral.biom -m agp_processing/4/ag-cleaned.txt -o agp_processing/8/per-sample-results/000009079.1130049 -c "COSMETICS_FREQUENCY:agp_processing/7.5/ag-100nt-1k-skin-cosmetics.biom, DOMINANT_HAND:agp_processing/7.5/ag-100nt-1k-skin-hand.biom, SEX:agp_processing/7.5/ag-100nt-1k-skin-sex.biom, AGE_CAT:agp_processing/7.5/ag-100nt-1k-skin-age.biom" -t skin -s 000009079.1130049

make_phyla_plots_AGP.py -i agp_processing/7.5/ag-100nt-1k-fecal.biom -m agp_processing/4/ag-cleaned.txt -o agp_processing/8/per-sample-results/000002035.1075856 -c "DIET_TYPE:agp_processing/7.5/ag-100nt-1k-fecal-diet.biom, BMI_CAT:agp_processing/7.5/ag-100nt-1k-fecal-bmi.biom, SEX:agp_processing/7.5/ag-100nt-1k-fecal-sex.biom, AGE_CAT:agp_processing/7.5/ag-100nt-1k-fecal-age.biom" -t fecal -s 000002035.1075856
make_phyla_plots_AGP.py -i agp_processing/7.5/ag-100nt-1k-fecal.biom -m agp_processing/4/ag-cleaned.txt -o agp_processing/8/per-sample-results/000002244.1076101 -c "DIET_TYPE:agp_processing/7.5/ag-100nt-1k-fecal-diet.biom, BMI_CAT:agp_processing/7.5/ag-100nt-1k-fecal-bmi.biom, SEX:agp_processing/7.5/ag-100nt-1k-fecal-sex.biom, AGE_CAT:agp_processing/7.5/ag-100nt-1k-fecal-age.biom" -t fecal -s 000002244.1076101
make_phyla_plots_AGP.py -i agp_processing/7.5/ag-100nt-1k-fecal.biom -m agp_processing/4/ag-cleaned.txt -o agp_processing/8/per-sample-results/000005844.1131797 -c "DIET_TYPE:agp_processing/7.5/ag-100nt-1k-fecal-diet.biom, BMI_CAT:agp_processing/7.5/ag-100nt-1k-fecal-bmi.biom, SEX:agp_processing/7.5/ag-100nt-1k-fecal-sex.biom, AGE_CAT:agp_processing/7.5/ag-100nt-1k-fecal-age.biom" -t fecal -s 000005844.1131797
make_phyla_plots_AGP.py -i agp_processing/7.5/ag-100nt-1k-fecal.biom -m agp_processing/4/ag-cleaned.txt -o agp_processing/8/per-sample-results/000009111.1131950 -c "DIET_TYPE:agp_processing/7.5/ag-100nt-1k-fecal-diet.biom, BMI_CAT:agp_processing/7.5/ag-100nt-1k-fecal-bmi.biom, SEX:agp_processing/7.5/ag-100nt-1k-fecal-sex.biom, AGE_CAT:agp_processing/7.5/ag-100nt-1k-fecal-age.biom" -t fecal -s 000009111.1131950
make_phyla_plots_AGP.py -i agp_processing/7.5/ag-100nt-1k-fecal.biom -m agp_processing/4/ag-cleaned.txt -o agp_processing/8/per-sample-results/000007118.1075682 -c "DIET_TYPE:agp_processing/7.5/ag-100nt-1k-fecal-diet.biom, BMI_CAT:agp_processing/7.5/ag-100nt-1k-fecal-bmi.biom, SEX:agp_processing/7.5/ag-100nt-1k-fecal-sex.biom, AGE_CAT:agp_processing/7.5/ag-100nt-1k-fecal-age.biom" -t fecal -s 000007118.1075682
```

And we'll end with some numbers on the number of successful and unsuccessful samples.

```python
>>> print "Number of successfully processed samples: %d" % len([l for l in open(successful_ids) if not l.startswith('#')])
>>> print "Number of unsuccessfully processed samples: %d" % len([l for l in open(unsuccessful_ids) if not l.startswith('#')])
Number of successfully processed samples: 15
Number of unsuccessfully processed samples: 0
```

```python

```
