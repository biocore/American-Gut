```python
>>> import os
>>> from functools import partial
>>> import pandas as pd
...
>>> import americangut.notebook_environment as agenv
>>> import americangut.util as agu
>>> import americangut.results_utils as agru
>>> import americangut.per_sample as agps
>>> import americangut.per_category as agpc
>>> import americangut.parallel as agpar
...
>>> chp_path = agenv.activate('10-populated-templates')
```

```python
>>> result_pdfs          = agu.get_new_path(agenv.paths['populated-templates']['result-pdfs'])
>>> result_taxa          = agu.get_new_path(agenv.paths['populated-templates']['result-taxa'])
>>> result_alpha         = agu.get_new_path(agenv.paths['populated-templates']['result-alpha'])
>>> successful_pdfs      = agu.get_new_path(agenv.paths['populated-templates']['successful-pdfs'])
>>> unsuccessful_pdfs    = agu.get_new_path(agenv.paths['populated-templates']['unsuccessful-pdfs'])
...
>>> os.mkdir(result_pdfs)
>>> os.mkdir(result_taxa)
>>> os.mkdir(result_alpha)
```

```python
>>> per_sample_results = agu.get_existing_path(agenv.paths['per-sample']['results'])
>>> successful_ids     = agu.get_existing_path(agenv.paths['per-sample']['successful-ids'])
```

```python
>>> ids = pd.read_csv(successful_ids, sep='\t', dtype=str)['#SampleID']
```

```python
>>> def create_pdf(opts, ids):
...     cmd_fmt = "echo '\n\def\yourname{unidentified}\n' >> %(result_path)s/macros.tex;"
...     cmd_fmt += "cd %(result_path)s; lualatex %(id)s.tex"
...     return agps._iter_ids_over_system_call(cmd_fmt, ids, opts)
...
>>> def aggregate(opts, ids):
...     cmd_fmt = "echo '\n\def\yourname{unidentified}\n' >> %(result_path)s/macros.tex;"
...     cmd_fmt += "mv %(result_path)s/%(id)s.pdf " + opts['populated-templates']['result-pdfs']
...     cmd_fmt += '; '
...     cmd_fmt += "mv %(result_path)s/%(id)s.txt " + opts['populated-templates']['result-taxa']
...     cmd_fmt += '; '
...     cmd_fmt += "mv %(result_path)s/shannon_%(id)s.png " + opts['populated-templates']['result-alpha']
...     cmd_fmt += '; '
...     cmd_fmt += "mv %(result_path)s/pd_%(id)s.png " + opts['populated-templates']['result-alpha']
...     return agps._iter_ids_over_system_call(cmd_fmt, ids, opts)
...
>>> opts = agps.create_opts('sample-agnostic', chp_path, None, [])
>>> process_pdf = partial(agps.sample_type_processor, [create_pdf, aggregate], opts)
```

We also need to write out the taxa summary files for each of the categories in the collapsed data. These will live in the same folder as the participant's taxonomy files.

```python
>>> agpc.cat_taxa_summaries()
```

Finally, we need to generate alpha diversity figures for each of the categories in the collapsed data. These allow for comparisons on the AG website.

```python
>>> agpc.cat_alpha_plots()
```

```python
>>> partitions = [(process_pdf, ids)]
>>> with open(successful_pdfs, 'w') as successful_pdfs_fp, open(unsuccessful_pdfs, 'w') as unsuccessful_pdfs_fp:
...     agpar.dispatcher(successful_pdfs_fp, unsuccessful_pdfs_fp, partitions)
```
