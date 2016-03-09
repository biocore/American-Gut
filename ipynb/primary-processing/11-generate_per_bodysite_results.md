This notebook will seperate the American Gut dataset by bodysite and rarefaction depth to generate datasets which are easy to use for analysis.

```python
>>> import itertools
>>> import os
...
>>> import numpy as np
>>> import pandas as pd
...
>>> import americangut.util as agu
>>> import americangut.notebook_environment as agenv
...
>>> chp_path = agenv.activate('11-packaged')
```

We're also going to set up a series of reference variables, which we'll use through this analysis.

```python
>>> alpha_metrics = {'pd': 'PD_whole_tree',
...                  'chao1': 'chao1',
...                  'shannon': 'shannon',
...                  'observedotus': 'observed_otus'
...                  }
>>> rare_depths = agenv.get_rarefaction_depth()
>>> depth_names = ['1k', '10k']
```

We'll start by splitting the raw OTU table and metadata by the `BODY_HABITAT` field, so we have an OTU table for each bodysite.

```python
>>> for trim in ['notrim', '100nt']:
...     for depth in ['raw', '1k', '10k']:
...         meta = agu.get_existing_path(agenv.paths['package']['split']['%s-map' % trim])
...         otus = agu.get_existing_path(agenv.paths['package']['split']['%s-%s-otu' % (trim, depth)])
...         odir = agu.get_new_path(agenv.paths['package']['split']['%s-%s-dir' % (trim, depth)])
...
...         !mkdir -p $odir
...         !split_otu_table.py -i $otus \
...                             -m $meta \
...                             -o $odir \
...                             -f BODY_HABITAT
```

Next, we'll move the unrarefied tables, split by bodysite, for each of the bodysites we're interested in. For this, we're going to write a quick parameter iterator, that will let us generate the samples we need.

In this block of code, and the rest of the notebook, we'll refer to files which already exist as source files. These will be modified and moved to a final location, refered to here as sink files.

```python
>>> def iterate_bodysite():
...     for bodysite in ['fecal', 'oral', 'skin']:
...         for trim in ['notrim', '100nt']:
...             dataset = agenv.paths['package']['all_participants_all_samples'][bodysite][trim]
...
...             # Gets the original filepath
...             source_otu = '"%s"' % agu.get_existing_path(dataset['source-unrare-otu'])
...             source_map = '"%s"' % agu.get_existing_path(dataset['source-unrare-map'])
...
...             # Gets the final paths
...             sink_dir =  agu.get_path(dataset['sink-unrare-dir'])
...             sink_otu = agu.get_path(dataset['sink-unrare-otu'])
...             sink_map = agu.get_path(dataset['sink-unrare-map'])
...
...             yield bodysite, trim, source_otu, source_map, sink_dir, sink_otu, sink_map
```

Now, let's move the split files to an easier to access location.

```python
>>> for bodysite, trim, source_otu, source_map, sink_dir, sink_otu, sink_map in iterate_bodysite():
...
...     # Makes the new directory
...     !mkdir -p $sink_dir
...
...     # Copies the files over
...     !cp $source_otu $sink_otu
...     !cp $source_map $sink_map
```

We'll follow by moving over the files for all participants. We're going to have a mapping file including alpha diversity, a rarefied OTU table, and distance matrices for each rarefaction depth. Once again, we'll build a helper function to find our file names, and then move them over.

```python
>>> def rarefied_parameter_iterator():
...     for bodysite in ['fecal', 'oral', 'skin']:
...         for trim in ['notrim', '100nt']:
...             for depth_id, depth in enumerate(['1k', '10k']):
...                 dataset = agenv.paths['package']['all_participants_all_samples'][bodysite][trim]
...                 source_otu = '"%s"' % agu.get_existing_path(dataset['source-%s-otu' % depth])
...                 source_map = agu.get_existing_path(dataset['source-%s-map' % depth])
...
...                 source_alpha = {
...                     '%s_%s' % (alpha_metrics[met], dep):
...                     agu.get_existing_path(dataset['source-%s-%s' % (dep, met)])
...                     for (met, dep) in itertools.product(alpha_metrics.keys(), depth_names[:(depth_id+1)])
...                 }
...
...                 source_beta = [
...                     agu.get_existing_path(dataset['source-%s-unweighted-unifrac' % depth]),
...                     agu.get_existing_path(dataset['source-%s-weighted-unifrac' % depth]),
...                 ]
...
...                 sink_dir = agu.get_path(dataset['sink-%s-dir' % depth])
...                 sink_otu = agu.get_path(dataset['sink-%s-otu' % depth])
...                 sink_map = agu.get_path(dataset['sink-%s-map' % depth])
...                 sink_beta = [
...                     agu.get_path(dataset['sink-%s-unweighted-unifrac' % depth]),
...                     agu.get_path(dataset['sink-%s-weighted-unifrac' % depth]),
...                     ]
...
...                 yield (bodysite, depth, source_otu, source_map, source_alpha, source_beta, sink_dir, sink_otu, sink_map, sink_beta)
```

Let's use our helper function to add the alpha diversity to our mapping file, and move the rest of the data we need.

```python
>>> for (bodysite, depth, source_otu, source_map, source_alpha, source_beta, sink_dir, sink_otu, sink_map, sink_beta) in rarefied_parameter_iterator():
...     # Creates the sink directory
...     !mkdir -p $sink_dir
...
...     # Loads the mapping file
...     map_ = pd.read_csv(source_map, sep='\t', dtype=str)
...     map_.set_index('#SampleID', inplace=True)
...
...     # Loads the alpha diveristy file
...     alphas = {metric : pd.read_csv(fp, sep='\t', dtype=str).set_index('Unnamed: 0').transpose()
...               for (metric, fp) in source_alpha.iteritems()
...              }
...
...     # Adds alpha diversity to the map
...     alpha_map = agu.add_alpha_diversity(map_, alphas)
...     alpha_map.to_csv(sink_map, sep='\t', index_label='#SampleID')
...
...     # Moves the OTU table
...     !cp $source_otu $sink_otu
...
...     # Moves the beta diversity tables
...     for source, sink in zip(*(source_beta, sink_beta)):
...         !cp $source $sink
```

We're going to select a single sample for each subject at each body site. The human microbiome project [PMID: [22699609](http://www.ncbi.nlm.nih.gov/pubmed/22699609)] demonstrated that differences across bodysites are larger than interpersonal differences within a single bodysite. However, other studies (ie Caporaso 2011 [[PMID: 21624126](http://www.ncbi.nlm.nih.gov/pubmed/21624126)], Wu 2011 [[PMID: 21885731](http://www.ncbi.nlm.nih.gov/pubmed/21885731)]) suggest a high degree of stability within individual communities. A single sample can be use to avoid the bias associated with multiple closely related samples from a single individual.

We select the sample using a constrained random method. Samples which are sequenced at depths greater than the highest rarefaction depth are given priority. However, if multiple samples exist greater than this depth, the sample is selected at random. The process is repeated for each rarefaction depth, until a sample can be identified. This way, the same sample can be retained across analysis performed at multiple depths.

```python
>>> for bodysite, trim, source_otu, source_map, sink_dir, sink_otu, sink_map in iterate_bodysite():
...     # We only need to do this once
...     if not trim == '100nt':
...         continue
...
...     # Reads in the map and adds in the sequencing depth
...     map_ = pd.read_csv(sink_map, sep='\t', dtype=str).set_index('#SampleID')
...     depths = !biom summarize-table -i $sink_otu
...     depths = pd.DataFrame([d.split(': ') for d in depths[15:]], columns=['#SampleID', 'depth']).set_index('#SampleID')
...     map_ = map_.join(depths)
...
...     # Picks a single sample for each individual
...     single_ids = agu.get_single_id_lists(map_, rare_depths)
...
...     # Saves the files
...     single_fp = agu.get_new_path(agenv.paths['package']['single_ids']['%s-%s' % (bodysite, 'unrare')])
...     with open(single_fp, 'w') as f_:
...         f_.write('\n'.join(single_ids['unrare']))
...
...     for depth, depth_name in zip(*(rare_depths, depth_names)):
...         single_fp = agu.get_new_path(agenv.paths['package']['single_ids']['%s-%s' % (bodysite, depth_name)])
...         with open(single_fp, 'w') as f_:
...             f_.write('\n'.join(single_ids[depth]))
```

We're going to define another helper function, to work through the single sample files. We're then going to filter the data.

```python
>>> def single_sample_iterator():
...     for bodysite in ['fecal', 'oral', 'skin']:
...         for depth_name in depth_names:
...             single_fp = agu.get_existing_path(agenv.paths['package']['single_ids']['%s-%s' % (bodysite, depth_name)])
...             for trim in ['notrim', '100nt']:
...                 dataset = agenv.paths['package']['all_participants_one_sample'][bodysite][trim]
...
...                 source_otu = agu.get_existing_path(dataset['source-%s-otu' % depth_name])
...                 source_map = agu.get_existing_path(dataset['source-%s-map' % depth_name])
...                 source_unweighted = agu.get_existing_path(dataset['source-%s-unweighted-unifrac' % depth_name])
...                 source_weighted = agu.get_existing_path(dataset['source-%s-weighted-unifrac' % depth_name])
...
...                 sink_dir = agu.get_path(dataset['sink-%s-dir' % depth_name])
...                 sink_otu = agu.get_new_path(dataset['sink-%s-otu' % depth_name])
...                 sink_map = agu.get_new_path(dataset['sink-%s-map' % depth_name])
...                 sink_unweighted = agu.get_new_path(dataset['sink-%s-unweighted-unifrac' % depth_name])
...                 sink_weighted = agu.get_new_path(dataset['sink-%s-weighted-unifrac' % depth_name])
...
...                 yield (bodysite, depth_name, trim, single_fp,
...                        source_otu, source_map, source_unweighted, source_weighted,
...                        sink_dir, sink_otu, sink_map, sink_unweighted, sink_weighted)
```

Let's use the helper function to partition the data into a single sample from each participant.

```python
>>> for (bodysite, depth_name, trim, single_fp, source_otu, source_map, source_unweighted, source_weighted,
...      sink_dir, sink_otu, sink_map, sink_unweighted, sink_weighted) in single_sample_iterator():
...
...     # Makes the destination directory
...     !mkdir -p $sink_dir
...
...     # Filters the OTU table, and distance matrices
...     !filter_samples_from_otu_table.py -i $source_otu -o $sink_otu --sample_id_fp $single_fp
...     !filter_distance_matrix.py -i $source_unweighted -o $sink_unweighted --sample_id_fp $single_fp
...     !filter_distance_matrix.py -i $source_weighted -o $sink_weighted --sample_id_fp $single_fp
...
...     # Filters the mapping file
...     with open(single_fp, 'r') as f_:
...         ids = np.array([i for i in f_.read().split('\n')])
...     map_ = pd.read_csv(source_map, sep='\t', dtype=str)
...     map_.set_index('#SampleID', inplace=True)
...     smap = map_.loc[ids]
...     smap.to_csv(sink_map, sep='\t', index_label='#SampleID')
```

Finally, we're going to go through the fecal samples, and generate a set of samples loosely described as "healthy adults". These are individuals 20-69 with a BMI between 18.5 and 30 who report no history of diabetes, inflammatory bowel disease, and have not used antibiotics in the past year.

```python
>>> for trim in ['notrim', '100nt']:
...     for sample_number in ['all_samples', 'one_sample']:
...         for depth_name in ['1k', '10k']:
...             dataset = agenv.paths['package']['sub_participants_%s' % sample_number]['fecal'][trim]
...             source_otu = agu.get_existing_path(dataset['source-%s-otu' % depth_name])
...             source_map = agu.get_existing_path(dataset['source-%s-map' % depth_name])
...             source_unweighted = agu.get_existing_path(dataset['source-%s-unweighted-unifrac' % depth_name])
...             source_weighted = agu.get_existing_path(dataset['source-%s-weighted-unifrac' % depth_name])
...
...             sink_dir = agu.get_path(dataset['sink-%s-dir' % depth_name])
...             sink_otu = agu.get_new_path(dataset['sink-%s-otu' % depth_name])
...             sink_map = agu.get_new_path(dataset['sink-%s-map' % depth_name])
...             sink_unweighted = agu.get_new_path(dataset['sink-%s-unweighted-unifrac' % depth_name])
...             sink_weighted = agu.get_new_path(dataset['sink-%s-weighted-unifrac' % depth_name])
...
...             !mkdir -p $sink_dir
...
...             # Filters the OTU tables and distance matrices
...             !filter_samples_from_otu_table.py -i $source_otu -m $source_map -o $sink_otu -s 'SUBSET_HEALTHY:true'
...             !filter_distance_matrix.py -i $source_unweighted -m $source_map -o $sink_unweighted -s 'SUBSET_HEALTHY:true'
...             !filter_distance_matrix.py -i $source_weighted -m $source_map -o $sink_weighted -s 'SUBSET_HEALTHY:true'
...
...             # Filters the mapping file
...             map_ = pd.read_csv(source_map, sep='\t', dtype=str).set_index('#SampleID')
...             smap = map_.loc[map_.SUBSET_HEALTHY == 'true']
...             smap.to_csv(sink_map, sep='\t', index_label='#SampleID')
```

We're going to add a readme file to the directory.

```python
>>> agenv.write_readme(chp_path)
```

Finally, we're going to check that we've generated the files we expect.

```python
>>> error = False
>>> message = []
...
>>> for path in agenv.paths['package']['split'].values():
...     filepath = agu.get_path(path)
...     if not os.path.exists(filepath):
...             message.append("Could not find: %s" % filepath)
...             error = True
...
>>> for path in agenv.paths['package']['single_ids'].values():
...     filepath = agu.get_path(path)
...     if not os.path.exists(filepath):
...             message.append("Could not find: %s" % filepath)
...             error = True
...
>>> for collection in ['all_participants_all_samples', 'all_participants_one_sample',
...                    'sub_participants_all_samples', 'sub_participants_one_sample']:
...     for bodysite in ['fecal', 'oral', 'skin']:
...         if collection in {'sub_participants_all_samples', 'sub_participants_one_sample'} and bodysite in {'oral', 'skin'}:
...             continue
...         for trim in ['notrim', '100nt']:
...             for key, path in agenv.paths['package'][collection][bodysite][trim].iteritems():
...                 filepath = agu.get_path(path)
...                 if not os.path.exists(filepath):
...                     message.append("Could not find: %s" % filepath)
...                     error = True
...
>>> if error:
...     raise RuntimeError('\n'.join(message))
```
