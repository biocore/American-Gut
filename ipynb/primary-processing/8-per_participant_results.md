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
>>> ag_l6_taxa_biom       = agu.get_existing_path(agenv.paths['ag-L6-taxa-biom'])
>>> ag_l6_taxa_oral_biom  = agu.get_existing_path(agenv.paths['ag-L6-taxa-oral-biom'])
>>> ag_l6_taxa_skin_biom  = agu.get_existing_path(agenv.paths['ag-L6-taxa-skin-biom'])
>>> ag_l6_taxa_fecal_biom = agu.get_existing_path(agenv.paths['ag-L6-taxa-fecal-biom'])
>>> ag_cleaned_md         = agu.get_existing_path(agenv.paths['ag-cleaned-md'])
...
>>> full_md               = agu.get_existing_path(agenv.paths['ag-pgp-hmp-gg-cleaned-md'])
>>> full_dm               = agu.get_existing_path(agenv.paths['ag-pgp-hmp-gg-100nt-1k-bdiv-unifrac'])
>>> full_pc               = agu.get_existing_path(agenv.paths['ag-pgp-hmp-gg-100nt-1k-unifrac-pc'])
...
>>> ag_gg_dm              = agu.get_existing_path(agenv.paths['ag-gg-100nt-1k-bdiv-unifrac'])
>>> ag_gg_pc              = agu.get_existing_path(agenv.paths['ag-gg-100nt-1k-unifrac-pc'])
>>> ag_gg_ss_dm           = agu.get_existing_path(agenv.paths['ag-gg-100nt-1k-bdiv-subsampled-unifrac'])
>>> ag_gg_ss_pc           = agu.get_existing_path(agenv.paths['ag-gg-100nt-1k-subsampled-unifrac-pc'])
...
>>> ag_skin_dm            = agu.get_existing_path(agenv.paths['ag-100nt-skin-1k-bdiv-unifrac'])
>>> ag_skin_pc            = agu.get_existing_path(agenv.paths['ag-100nt-skin-1k-unifrac-pc'])
>>> ag_oral_dm            = agu.get_existing_path(agenv.paths['ag-100nt-oral-1k-bdiv-unifrac'])
>>> ag_oral_pc            = agu.get_existing_path(agenv.paths['ag-100nt-oral-1k-unifrac-pc'])
>>> ag_fecal_dm           = agu.get_existing_path(agenv.paths['ag-100nt-fecal-1k-bdiv-unifrac'])
>>> ag_fecal_pc           = agu.get_existing_path(agenv.paths['ag-100nt-fecal-1k-unifrac-pc'])
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

These next functions actually perform the processing per site. These functions were not pushed down into the library code as to make it easier for developers to iterate and modify the processing code.

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
...             results[id_] = 'FAILED (%s): %s' % (stderr.splitlines()[-1], cmd)
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
...         A dict of relevant paths. This is expected to contain the following
...         keys:
...             {'site_table', 'site_dm', 'site_pc', 'full_table', 'full_md',
...              'full_dm', 'full_pc', 'output_path'}
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
...         A dict of relevant paths. This is expected to contain the following
...         keys:
...             {'site_table', 'site_dm', 'site_pc', 'full_table', 'full_md',
...              'full_dm', 'full_pc', 'output_path'}
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
...         A dict of relevant paths. This is expected to contain the following
...         keys:
...             {'site_table', 'site_dm', 'site_pc', 'full_table', 'full_md',
...              'full_dm', 'full_pc', 'output_path'}
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
...         A dict of relevant paths. This is expected to contain the following
...         keys:
...             {'site_table', 'site_dm', 'site_pc', 'full_table', 'full_md',
...              'full_dm', 'full_pc', 'output_path'}
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
...                         "--mapping_file %s" % paths['full_md'],
...                         "--output %s" % paths['output_path'],
...                         "--prefix Figure_2",
...                         "--samples %s"])
...
...     return _iter_ids_over_system_call(cmd_fmt, sample_ids)
...
>>> def gradient_pcoa(sample_ids, paths):
...     pass
...
>>> # NEED BAR/PIE CHARTS, AND PCOAS
... #x = {
...   #  'mod2 pcoas fig3': 'mod2_pcoa.py gradient --coords %(coords)s --mapping_file %(mapping)s --output %(output)s --prefix Figure_3 --color %(color)s --samples %(samples)s',
...    # 'Make Pie Charts': 'make_pie_plot_AGP.py -i %(input)s -o %(output)s -s %(samples)s',
...     #'mod2 pcoas fig2 subsample': '
```

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
```

```python
>>> def get_paths_thread_local(extras):
...     import americangut.notebook_environment as agenv_local
...     import americangut.util as agu_local
...
...     paths = {'full_table': agu_local.get_existing_path(agenv_local.paths['ag-L6-taxa-biom']),
...              'full_md': agu_local.get_existing_path(agenv_local.paths['ag-cleaned-md']),
...              'full_dm': agu_local.get_existing_path(agenv_local.paths['ag-pgp-hmp-gg-100nt-1k-bdiv-unifrac']),
...              'full_pc': agu_local.get_existing_path(agenv_local.paths['ag-pgp-hmp-gg-100nt-1k-unifrac-pc']),
...              'ag_gg_dm': agu_local.get_existing_path(agenv_local.paths['ag-gg-100nt-1k-bdiv-unifrac']),
...              'ag_gg_ss_pc': agu_local.get_existing_path(agenv_local.paths['ag-gg-100nt-1k-subsampled-unifrac-pc']),
...              'output_path': agu_local.get_existing_path(agenv_local.paths['per-sample-results'])}
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
...     paths = get_paths_thread_local([('ag-L6-taxa-fecal-biom', 'site_table')])
...
...     functions = [taxa_summaries,
...                  taxon_significance,
...                  body_site_pcoa,
...                  country_pcoa
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
...     paths = get_paths_thread_local([('ag-L6-taxa-oral-biom', 'site_table')])
...
...     functions = [taxa_summaries,
...                  taxon_significance,
...                  body_site_pcoa,
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
...     paths = get_paths_thread_local([('ag-L6-taxa-skin-biom', 'site_table')])
...
...     functions = [taxa_summaries,
...                  taxon_significance,
...                  body_site_pcoa,
...                 ]
...
...     return merge_error_reports(*[f(ids, paths) for f in functions])
```

And now, let's start mass generating figures!

```python
>>> with open(successful_ids, 'w') as successful_ids_fp, open(unsuccessful_ids, 'w') as unsuccessful_ids_fp:
...     dispatcher(successful_ids_fp, unsuccessful_ids_fp, partition_samples_by_bodysite)
>>> print time.time() - start
 29.3168480396
```

And we'll end with some numbers on the number of successful and unsuccessful samples.

```python
>>> print "Number of successfully processed samples: %d" % len([l for l in open(successful_ids) if not l.startswith('#')])
>>> print "Number of unsuccessfully processed samples: %d" % len([l for l in open(unsuccessful_ids) if not l.startswith('#')])
Number of successfully processed samples: 12
Number of unsuccessfully processed samples: 3
```
