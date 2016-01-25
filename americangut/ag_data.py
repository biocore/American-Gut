import os

import biom
import numpy as np
import pandas as pd
import skbio

import americangut as ag


class AgData:
    na_values = ['NA', 'unknown', '', 'no_data', 'None', 'Unknown']
    alpha_columns = {'10k': ['PD_whole_tree_1k', 'PD_whole_tree_10k',
                             'shannon_1k', 'shannon_10k',
                             'chao1_1k', 'chao1_10k',
                             'observed_otus_1k', 'observed_otus_10k'],
                     '1k': ['PD_whole_tree_1k', 'shannon_1k', 'chao1_1k',
                            'observed_otus_1k']}
    bodysite_limits = {'fecal': 55,
                       }
    subset_columns = {'IBD': 'SUBSET_IBD',
                      'DIABETES': 'SUBSET_DIABETES',
                      'BMI_CAT': 'SUBSET_BMI',
                      'ANTIBIOTIC_HISTORY': 'SUBSET_ANTIBIOTIC_HISTORY',
                      'AGE_CAT': 'SUBSET_AGE'}

    def __init__(self, bodysite, trim, depth, sub_participants=False,
                 one_sample=False, base_dir=None):
        """A pretty object to contain datasets for analysis

        Parameters
        ----------
        bodysite : {'fecal', 'oral', 'skin'}
            The bodysite samples for the data set
        trim : {'100nt', 'notrim'}
            The sequence length for the dataset
        depth : {'1k', '10k'}
            The rarefaction depth to analyze.
        sub_participants : bool
            Whether the healthy subset of adults should be used. This is
            defined as participants 20-69 with BMI between 18.5-30, and no
            reported history of diabetes, inflammatory bowel disease, or
            antibiotic use in the past year.
        one_sample : bool
            Whether one sample per participant should be used
        base_dir : str
            The path where datasets are currently located. The path should
            terminate with the 11-packaged folder if the data was generated
            with the American Gut primary-processing notebook set. Sub
            directories should be layed out according to the directory
            structure in those notebooks.
        """

        # Gets the base directory
        if base_dir is None:
            base_dir = os.path.join(
                ag.WORKING_DIR.split('American-Gut')[0],
                'American-Gut/ipynb/primary-processing/agp_processing/'
                '11-packaged'
                )
        # Checks the inputs
        if not os.path.exists(base_dir):
            raise ValueError('The base directory does not exist!')
        if bodysite not in {'fecal', 'oral', 'skin'}:
            raise ValueError('%s is not a supported bodysite.' % bodysite)
        if trim not in {'notrim', '100nt'}:
            raise ValueError('%s is not a supported trim length.' % trim)
        if depth not in {'1k', '10k'}:
            raise ValueError('%s is not a supported rarefaction depth' % depth)

        self.base_dir = base_dir
        self.bodysite = bodysite
        self.trim = trim
        self.depth = depth
        self.sub_participants = sub_participants
        self.one_sample = one_sample

        self.reload_files()

    def _check_dataset(self):
        """Makes sure the directory is pointing toward the correct directory"""

        self.data_set = {'bodysite': self.bodysite,
                         'sequence_trim': self.trim,
                         'rarefaction_depth': self.depth,
                         }

        if self.one_sample:
            self.data_set['samples_per_participants'] = 'one_sample'
        else:
            self.data_set['samples_per_participants'] = 'all_samples'

        if self.sub_participants:
            self.data_set['participant_set'] = 'sub'
        else:
            self.data_set['participant_set'] = 'all'

        # Gets the data directory
        self.data_dir = os.path.join(self.base_dir,
                                     '%(bodysite)s/'
                                     '%(sequence_trim)s/'
                                     '%(participant_set)s_participants/'
                                     '%(samples_per_participants)s/'
                                     '%(rarefaction_depth)s'
                                     % self.data_set
                                     )
        if not os.path.exists(self.data_dir):
            raise ValueError('The identified dataset does not exist')

    def reload_files(self):
        """Loads the map, otu table, and distance matrices from file"""
        self._check_dataset()

        # Determines the filepaths
        map_fp = os.path.join(self.data_dir,
                              'ag_%(rarefaction_depth)s_%(bodysite)s.txt'
                              % self.data_set)
        otu_fp = os.path.join(self.data_dir,
                              'ag_%(rarefaction_depth)s_%(bodysite)s.biom'
                              % self.data_set)
        unweighted_fp = os.path.join(self.data_dir,
                                     'unweighted_unifrac_ag_'
                                     '%(rarefaction_depth)s_%(bodysite)s.txt'
                                     % self.data_set)
        weighted_fp = os.path.join(self.data_dir,
                                   'weighted_unifrac_ag_%(rarefaction_depth)s'
                                   '_%(bodysite)s.txt'
                                   % self.data_set)

        self.map_ = pd.read_csv(map_fp,
                                sep='\t',
                                dtype=str,
                                na_values=self.na_values)
        self.map_.set_index('#SampleID', inplace=True)
        columns = self.alpha_columns[self.data_set['rarefaction_depth']]
        self.map_[columns] = self.map_[columns].astype(float)

        self.otu_ = biom.load_table(otu_fp)

        self.beta = {'unweighted_unifrac':
                     skbio.DistanceMatrix.read(unweighted_fp),
                     'weighted_unifrac':
                     skbio.DistanceMatrix.read(weighted_fp)
                     }

    def drop_alpha_outliers(self, metric='PD_whole_tree_10k', limit=55):
        """Removes samples with alpha diversity above the specified range

        In rounds 1-21 of the AmericanGut, there is a participant who has an
        alpha diverisity value 4 standard deviations about the next highest
        sample. Therefore, we remove this individual as an outlier.

        Parameters
        ----------
        metric : str, optional
            The alpha diversity metric used to identify the outlier.
            Default is `PD_whole_tree_1k`.
        limit : float, optional
            The alpha diversity value which should serve as the cutoff.
            Default is 55.

        """

        if metric not in self.alpha_columns[
                self.data_set['rarefaction_depth']]:
            raise ValueError('%s is not a supported metric!' % metric)

        self.outliers = self.map_.loc[(self.map_[metric] >= limit)].index
        keep = self.map_.loc[(self.map_[metric] < limit)].index

        self.map_ = self.map_.loc[keep]
        self.otu_ = self.otu_.filter(keep)
        self.beta = {k: self.beta[k].filter(keep) for k in self.beta.keys()}

    def drop_bmi_outliers(self, bmi_limits=[5, 55]):
        """Corrects BMI values outside a specified range

        In rounds 1-21 of American Gut, participants entered BMI values which
        are imprbably (greater than 200). Therefore, we limit the BMI range
        to a more probable set of values.

        Parameters
        ----------
        bmi_limits : list, optional
            The range of BMI values to keep in a new column, called
            `BMI_CORRECTED` and `BMI_CAT`. Default is between 5 and 55.
        """

        self.map_['BMI'] = self.map_['BMI'].astype(float)

        keep = ((self.map_.BMI > bmi_limits[0]) &
                (self.map_.BMI < bmi_limits[1]))
        discard = ~ keep

        self.map_['BMI_CORRECTED'] = self.map_.loc[keep, 'BMI']
        self.map_.loc[discard, 'BMI_CAT'] = np.nan

    def return_dataset(self, group):
        """Returns a map, otu table, and distance matrix for the group

        Many tests curerntly avaliable do not respond well to missing data.
        So, we filter the dataset to remove any samples not in the dataset.

        Parameters
        ----------
        group : AgQuestion
            An AqQuestion data dictionary object specifying the field
            to be used.

        Returns
        -------
        map_ : DataFrame
            The metadata represented as a pandas DataFrame
        otu_ : biom
            The Biom object containing the sample x OTU counts
        beta : dict
            A dictioary keying beta diversity metric(s) to DistanceMatrix
            object(s)

        """
        defined = self.map_[~ pd.isnull(self.map_[group.name])].index

        map_ = self.map_.loc[defined]
        otu_ = self.otu_.filter(defined, axis='sample')
        beta = {k: v.filter(defined) for k, v in self.beta.iteritems()}

        return map_, otu_, beta

    def clean_group(self, group):
        """Leverages question structure to format the specified column in map_

        Parameters
        ----------
        group : AgQuestion
            An AqQuestion data dictionary object specifying the field
            to be used.

        """

        if group.type == 'Continous':
            group.drop_outliers(self.map_)

        elif group.type in {'Bool', 'Multiple'}:
            group.remap_data_type(self.map_)
            group.convert_to_word(self.map_)

        elif group.type == 'Categorical':
            group.remove_ambiguity(self.map_)
            group.remap_groups(self.map_)
            group.drop_infrequent(self.map_)

        elif group.type == 'Clinical':
            group.remap_clinical(self.map_)

        elif group.type == 'Frequency':
            group.clean(self.map_)
            group.combine_frequency(self.map_)

        if group.type == 'Multiple':
            group.correct_unspecified(self.map_)
