

```python
import os
import shutil
import copy
import pickle

import numpy as np
import scipy
import matplotlib.pylab as plt
import pandas as pd

import americangut.diversity_analysis as div_an

# Initializes the notebook with inline display
%matplotlib inline
```

```python
overwrite = False

working_dir = os.path.join(os.path.abspath('.'), 'picrust_analysis')
div_an.check_dir(working_dir)

download_dir = os.path.join(os.path.abspath('.'), 'ag-precomputed-rounds-1-21/fecal')

site_dir = os.path.join(working_dir, 'fecal')
div_an.check_dir(site_dir)

data_dir = os.path.join(site_dir, 'notrim/all_participants/one_sample')

data_otu_fp = os.path.join(data_dir, '10k/ag_10k_fecal.biom')
data_map_fp = os.path.join(data_dir, '10k/ag_10k_fecal.txt')
```

```python
# Gets data for the single sample per participant fecal sample directory
if overwrite or not (os.path.exists(data_map_fp)):
    # Downloads the files
    !curl -OL https://www.dropbox.com/s/la3q3zntacei1c2/all_participants_all_samples.tgz
    !curl -OL ftp://ftp.microbio.me/AmericanGut/ag-precomputed-rounds-1-21.tgz
     Extracts the data
    !tar -xzf ag-precomputed-rounds-1-21.tgz
     Moves the directory
    os.remove(os.path.join('.', 'ag-precomputed-rounds-1-21.tgz'))
    shutil.move(os.path.join('.', download_dir), working_dir)
```

Normalize the 16S abundances by predicted 16S copy number.
```python
normalize_by_copy_number.py
        -i your_otu_table.biom
        -o normalized_otus.biom
```

Predict the metagenome abundances based off of 16S abundances.
```python
predict_metagenomes.py
        -i normalized_otus.biom
        -o metagenome_predictions.biom
```

Collapse the metagenome abundances into KEGG level 3 pathways.
```python
categorize_by_function.py
        -i predicted_metagenomes.biom
        -c KEGG_Pathways
        -l 3
        -o predicted_metagenomes.L3.biom
```
