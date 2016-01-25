
The goal of this notebook is to provide demographic summaries for participants in the American Gut and associated projects. We look at metadata, and summarize the available information.

The information generated here will be used for Table 3 of the American Gut paper.

We'll start by importing the necessary libraries.


```python
import matplotlib
matplotlib.use('Agg')
```


```python
import os
from functools import partial

import numpy as np
import pandas as pd
```

Next, we're going to locate the mapping file and raw OTU table.


```python
processing_dir = os.path.abspath('../primary-processing/agp_processing/')
map_fp = os.path.join(processing_dir, '01-raw/metadata.txt')
```

Let’s read in the mapping file.


```python
md = pd.read_csv(map_fp, sep='\t', dtype=str, na_values=['NA', 'no_data', 'unknown', 'Unspecified', 'Unknown', 'None'])
md.set_index('#SampleID', inplace=True)
```

We're going to start by dropping the blank samples. These are samples where the body habitat is not defined.


```python
md = md[md['BODY_HABITAT'].notnull()]
```

Next, we're going to get a single sample for each person. The `HOST_SUBEJCT_ID` is a unique code for each participant. It can be used to group samples from the same person within the project.

```python
single = []
for indv, ids in md.groupby('HOST_SUBJECT_ID').groups.iteritems():
    single.append(ids[0])

md['Participants'] = np.nan
md['Samples'] = 1
md.loc[single, 'Participants'] = 1
```

Now, let's look at the countries where participants come from.


```python
countries = md.groupby('COUNTRY').count()[['Samples', 'Participants']]
print 'There are %i countries and sovereign states represented here.' % len(countries)
countries.sort_values('Samples', ascending=False)
```

Let's now compare the samples to the demographics of American Gut Samples from participants in the United States to the demographics of the US population, according to the 2010 census. We’ll filter to only look at participants from the US.


```python
md = md.loc[md['COUNTRY'] == 'USA']
```

We're going to look at sex, race, smoking, Diabetes and Inflammatory Bowel disease diagnosis, and Body Mass index. To do this, we're going to reformat some of the responses about age.


```python
def mapper(mapping, value):
    return mapping.get(value, value)

diabetes_values_fix = {'I do not have this condition': 'I do not have diabetes',
                       'Diagnosed by a medical professional (doctor, physician assistant)': 'I have diabetes',
                       "Diagnosed by an alternative medicine practitioner": "I have diabetes",
                       'Type I': 'I have diabetes',
                       'Type II': 'I have diabetes',
                       'Self-diagnosed': 'I have diabetes'}

ibd_values_fix = {"Crohn's disease": "I have an IBD",
                  "Diagnosed by a medical professional (doctor, physician assistant)": "I have IBD",
                  "Diagnosed by an alternative medicine practitioner": "I have IBD",
                  "I do not have this condition": "I do not have IBD",
                  "I do not have IBD": "I do not have IBD",
                  "Ulcerative colitis": "I have IBD",
                  "Self-diagnosed": "I have IBD"}

smoking_values_fix = {'Daily': 'I smoke',
                      'Never': 'I do not smoke',
                      'Occasionally (1-2 times/week)': 'I smoke',
                      'Rarely (a few times/month)': 'I smoke',
                      'Rarely (few times/month)': 'I smoke',
                      'Regularly (3-5 times/week)': 'I smoke'}

diabetes_map = partial(mapper, diabetes_values_fix)
ibd_map = partial(mapper, ibd_values_fix)
smoke_map = partial(mapper, smoking_values_fix)

md['DIABETES'] = md['DIABETES'].apply(diabetes_map)
md['IBD'] = md['IBD'].apply(ibd_map)
md['SMOKING_FREQUENCY'] = md['SMOKING_FREQUENCY'].apply(smoke_map)
```

Finally, there are some BMI values which are extremely high. We're going to cap the BMI range between 10 and 100.

We'll also exclude BMI categorization for anyone under the age of 18. According to the World Health Organization (WHO), BMI for children under 18 must be calculated based on their age and gender.


```python
md[['BMI', 'AGE_CORRECTED']] = md[['BMI', 'AGE_CORRECTED']].astype(float)

md.loc[md['BMI'] > 100, 'BMI_CAT'] = np.nan
md.loc[md['BMI'] < 10, 'BMI_CAT'] = np.nan

# md.loc[md['AGE_CORRECTED'] < 20, 'BMI_CAT'] = np.nan
```

Now, let's build the demographics table


```python
res_table = {}
n_samples = float(len(md))
cats = ['SEX', 'RACE', 'SMOKING_FREQUENCY', 'DIABETES', 'IBD', 'BMI_CAT']

for cat in cats:
    # drop out any null values
    cat_tab = md[md[cat].notnull()]
    
    # determine how many unique subjects are represented
    n_subjects = float(len(cat_tab.HOST_SUBJECT_ID.unique()))

    # for each value in (e.g., for SEX: Male, Female, Other)
    for val in cat_tab.groupby(cat).HOST_SUBJECT_ID.unique().index:
        # get the number of unique subjects
        count = cat_tab.groupby(cat).HOST_SUBJECT_ID.nunique()[val]
        
        # store the count of subjects and the percentage of the subjects represented
        res_table["%s - %s" % (cat, val)] = (count, np.round((count / n_subjects) * 100, 1))

res = pd.DataFrame.from_dict(res_table, orient='index')
res.columns = ['Count', 'Within group percentage']
```

Now, we'll define the census data


```python
# Category/value : percent in US population
census_data = {
               #from http://quickfacts.census.gov/qfd/states/00000.html
               'SEX - female': 50.8,
               'SEX - male': 49.2,  # this is an over estimate as only the % of females is described in the above URL
               'SEX - other': np.nan,  # does not appear to be tracked
               
               # from http://www.census.gov/prod/cen2010/briefs/c2010br-02.pdf
               # doesn't sum to 100% as the fields don't map exactly, so there may be some overlap represented below
               'RACE - African American': 12.6,
               'RACE - Asian or Pacific Islander': 5.0,
               'RACE - Caucasian': 63.7,
               'RACE - Hispanic': 16.3,
               'RACE - Other': 6.2,

               # from http://www.census.gov/compendia/statab/2012/tables/12s0211.pdf
               # using total, non age adjusted values

###### we probably want to filter to > 18yo for these values in the metadata
               'BMI_CAT - Normal': 31.2, 
               'BMI_CAT - Obese': 33.0,
               'BMI_CAT - Overweight': 34.0,
               'BMI_CAT - Underweight': 1.8,

               # from http://www.cdc.gov/diabetes/data/statistics/2014statisticsreport.html
               'DIABETES - I do not have diabetes': 90.7,
               'DIABETES - I have diabetes': 9.3, # This uses 21 million 

               # from http://www.cdc.gov/ibd/ibd-epidemiology.htm
               # using 1.3 million people as the estimate, and US population size for 2014 from
               # http://quickfacts.census.gov/qfd/states/00000.html
               'IBD - I do not have IBD': 99.6,
               'IBD - I have IBD': 0.4,
          
               # from http://www.cdc.gov/tobacco/data_statistics/fact_sheets/adult_data/cig_smoking/
               'SMOKING_FREQUENCY - I do not smoke': 82.6,
               'SMOKING_FREQUENCY - I smoke': 17.4
}

res['US Census/CDC/NHANES data percentages'] = pd.DataFrame.from_dict(census_data, orient='index')
```

Now, let's look at the results.


```python
res.sort_index()
```