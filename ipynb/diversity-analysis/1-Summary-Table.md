The goal of this notebook is to provide demographic summaries for participants in the American Gut and associated projects. We look at metadata, and summarize the available information.

The information generated here will be used for table 1 of the American Gut paper.

We'll start by importing the necessary libraries.

```python
>>> import matplotlib
>>> matplotlib.use('Agg')
...
>>> import numpy as np
>>> import pandas as pd
...
>>> import americangut.notebook_environment as agenv
>>> import americangut.util as agu
```

Next, we're going to locate the mapping file and raw OTU table.

```python
>>> map_fp = agu.get_existing_path(agenv.paths['raw']['metadata'])
```

Then, we'll load the mapping file.

```python
>>> map_ = pd.read_csv(map_fp,
...                     sep='\t',
...                     dtype=str,
...                     na_values=['NA', 'no_data', 'unknown', 'Unspecified', 'Unknown'])
>>> map_.set_index('#SampleID', inplace=True)
```

Now, let's start to build the first table. To do this, we'll define the counts from the Human Microbiome Project [[PMID: 22699609](http://www.ncbi.nlm.nih.gov/pubmed/22699609)]. The Human Microbiome Project looked at samples across 16 to 18 sites in a small number of healthy adults.

```python
>>> hmp_counts = pd.DataFrame([[ 365, 230],
...                            [1367, 238],
...                            [3316, 234],
...                            [ 482, 109],
...                            [ 339, 221],
...                            [   0,   0],
...                            [np.nan, np.nan]
...                            ],
...                            index=['Feces', 'Skin', 'Oral', 'Vagina', 'Nose', 'Hair', 'Blank'],
...                           columns=['HMP Samples', 'HMP Participants'])
```

We'll do the same calculation for the American Gut Project.

We'll use the field, `BODY_HABITAT` to identify where the sample was collected. We'll also infer that if no value is supplied for `BODY_HABITAT`, then we will assume it's a blank. We'll use a helper function to rename the values in `BODY_HABITAT` so the look clean.

```python
>>> def habitat_clean(x):
...     if x == 'None':
...         return 'Blank'
...     else:
...         return x.split(' ')[0].replace('UBERON:', '').title()
...
>>> map_['BODY_HABITAT'] = map_['BODY_HABITAT'].apply(habitat_clean)
```

We'll group the data by bodysite, and count the number of samples and participants assoicated iwth each bodysite.

```python
>>> ag_participants = pd.DataFrame.from_dict(
...     {site: {'AGP Participants': len(map_[map_['BODY_HABITAT'] == site].groupby('HOST_SUBJECT_ID').groups),
...             'AGP Samples': len(map_[map_['BODY_HABITAT'] == site])}
...     for site in set(map_['BODY_HABITAT'])}, orient='index')
```

Now that we've built the tables, let's merge them.

```python
>>> ag_participants.loc['Blank', 'AGP Participants'] = np.nan
>>> ag_participants.join(hmp_counts)
```
