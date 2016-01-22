
The goal of this notebook is to provide demographic summaries for participants in the American Gut and associated projects. We look at metadata, and summarize the available information.

The information generated here will be used for tables 1 of the American Gut paper.

We'll start by importing the necessary libraries.


```python
import matplotlib
matplotlib.use('Agg')

import os

import numpy as np
import pandas as pd
```

Next, we're going to locate the mapping file and raw OTU table.


```python
processing_dir = os.path.abspath('../primary-processing/agp_processing/')
map_fp = os.path.join(processing_dir, '01-raw/metadata.txt')
otu_fp = os.path.join(processing_dir, '03-otus/100nt/gg-13_8-97-percent/otu_table.biom')
```

We'll start by reading in the mapping file.


```python
md = pd.read_csv(map_fp, sep='\t', dtype=str, na_values=['NA', 'no_data', 'unknown', 'Unspecified', 'Unknown'])
md.set_index('#SampleID', inplace=True)
```

We're also going to calculate the number of sequences per sample, and add that to our mapping file.


```python
otu_summary = !biom summarize-table -i $otu_fp
seq_depth = pd.DataFrame(np.array([l.split(': ') for l in otu_summary[15:]]), columns=['#SampleID', 'counts'])
seq_depth.set_index('#SampleID', inplace=True)

md['count'] = seq_depth
```

Now, let's start to build the first table. To do this, we'll define the counts from the Human Microbiome Project [[PMID: 22699609](http://www.ncbi.nlm.nih.gov/pubmed/22699609)]. The Human Microbiome Project looked at samples across 16 to 18 sites in a small number of healthy adults. 


```python
hmp_counts = pd.DataFrame([[ 365, 230],
                           [1367, 238],
                           [3316, 234],
                           [ 482, 109],
                           [ 339, 221],
                           [   0,   0],
                           [np.nan, np.nan]
                          ],
                          index=['Feces', 'Skin', 'Oral', 'Vagina', 'Nose', 'Hair', 'Blank'],
                          columns=['HMP Samples', 'HMP Participants'])
```

Now, we'll do the same calculation for the American Gut Project. 

We'll use the field, `BODY_HABITAT` to identify where the sample was collected. We'll also infer that if no value is supplied for `BODY_HABITAT`, then we will assume it's a blank. We'll use a helper function to rename the values in `BODY_HABITAT` so the look clean.


```python
def habitat_clean(x):
    if x == 'None':
        return 'Blank’
    else:
        return x.split(' ')[0].replace('UBERON:', '').title()
    
md['BODY_HABITAT'] = md['BODY_HABITAT'].apply(habitat_clean)
```


```python
ag_samples = md.groupby('BODY_HABITAT').count().max(1)
```

Now, we're going to group each body site by the number of participants.


```python
ag_participants = pd.DataFrame.from_dict(
    {site: {'AGP Participants': len(md[md['BODY_HABITAT'] == site].groupby('HOST_SUBJECT_ID').groups),
            'AGP Samples': len(md[md['BODY_HABITAT'] == site])} 
    for site in set(md['BODY_HABITAT'])}, orient='index')
```

Now that we've built the tables, let's merge them.


```python
ag_participants.loc['Blank', 'AGP Participants'] = np.nan
ag_participants.join(hmp_counts)
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>AGP Participants</th>
      <th>AGP Samples</th>
      <th>HMP Samples</th>
      <th>HMP Participants</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Blank</th>
      <td>NaN</td>
      <td>687</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>Feces</th>
      <td>5380</td>
      <td>5952</td>
      <td>365</td>
      <td>230</td>
    </tr>
    <tr>
      <th>Hair</th>
      <td>5</td>
      <td>5</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>Nose</th>
      <td>7</td>
      <td>7</td>
      <td>339</td>
      <td>221</td>
    </tr>
    <tr>
      <th>Oral</th>
      <td>436</td>
      <td>477</td>
      <td>3316</td>
      <td>234</td>
    </tr>
    <tr>
      <th>Skin</th>
      <td>165</td>
      <td>337</td>
      <td>1367</td>
      <td>238</td>
    </tr>
    <tr>
      <th>Vagina</th>
      <td>13</td>
      <td>15</td>
      <td>482</td>
      <td>109</td>
    </tr>
  </tbody>
</table>
</div>



