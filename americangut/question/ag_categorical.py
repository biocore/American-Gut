import copy

import numpy as np
import pandas as pd

from americangut.question.ag_question import AgQuestion


class AgCategorical(AgQuestion):

    def __init__(self, name, description, dtype, order, **kwargs):
        """A question object for categorical or ordinal questions

        Parameters
        ----------
        name : str
            The name of column where the group is stored
        description : str
            A brief description of the data contained in the question
        dtype : object
            The datatype in which the responses should be represented. (i.e.
            float, int, str).
        order : list
            The list of responses to the question
        clean_name : str, optional
            A nicer version of the way the column should be named.
        remap : function
            A function used to remap the data in the question.
        frequency_cutoff : float, optional
            The minimum number of observations required to keep a sample group.
        drop_ambiguous : bool, optional
            Whether ambigious values should be removed from the dataset.
        ambigigous_values : set, optional
            A list of values which are considered ambiguous and removed when
            `drop_ambiguous` is `True`.
        free_response: bool, optional
            Whether the question is a free response question or controlled
            vocabulary
        mimmarks : bool, optional
            If the question was a mimmarks standard field
        ontology : str, optional
            The type of ontology, if any, which was used in the field value.

        """
        if dtype not in {str, bool, int, float}:
            raise ValueError('%s is not a supported datatype for a '
                             'categorical variable.' % dtype)
        # Initializes the question
        AgQuestion.__init__(self, name, description, dtype, **kwargs)

        self.type = 'Categorical'

        if len(order) == 1:
            raise ValueError('There cannot be possible response.')
        self.order = order
        self.extremes = kwargs.get('extremes', None)
        if self.extremes is None:
            self.extremes = [order[0], order[-1]]

        self.earlier_order = []
        self.frequency_cutoff = kwargs.get('frequency_cutoff', None)
        self.drop_ambiguous = kwargs.get('drop_ambiguous', True)
        self.ambiguous_values = kwargs.get(
            'ambiguous_values',
            {"none of the above", "not sure", "other", 'nan',
             "i don't know, i do not have a point of reference"
             })

    def _update_order(self, remap_):
        """Updates the order and earlier order arguments"""
        order = copy.copy(self.order)
        self.earlier_order.append(order)
        self.order = []
        for o in order:
            new_o = remap_(o)
            if new_o not in self.order and not pd.isnull(new_o):
                self.order.append(new_o)
        self.extremes = [remap_(e) for e in self.extremes]

    def remove_ambiguity(self, map_):
        """Removes ambiguous groups from the mapping file

        Parameters
        ----------
        map_ : DataFrame
            A pandas object containing the data to be analyzed. The
            Question `name` should be a column in the `map_`.

        """
        if self.drop_ambiguous:
            self.check_map(map_)
            self.remap_data_type(map_, watch=False)

            def _remap(x):
                if x.lower() in self.ambiguous_values:
                    return np.nan
                else:
                    return x

            index = map_[self.name].dropna().index
            map_[self.name] = map_.loc[index, self.name].apply(_remap)
            self._update_order(_remap)

    def remap_groups(self, map_):
        """Remaps columns in the column

        Parameters
        ----------
        map_ : DataFrame
            A pandas object containing the data to be analyzed. The
            Question `name` should be a column in the `map_`.

        """
        self.check_map(map_)
        self.remap_data_type(map_, watch=False)

        if self.remap_ is not None:
            index = map_[self.name].dropna().index
            map_[self.name] = map_.loc[index, self.name].apply(self.remap_)
            self._update_order(self.remap_)

    def drop_infrequent(self, map_):
        """Removes groups below a frequency cutoff

        Parameters
        ----------
        map_ : DataFrame
            A pandas object containing the data to be analyzed. The
            Question `name` should be a column in the `map_`.

        """
        self.check_map(map_)
        self.remap_data_type(map_, watch=False)

        counts = map_.groupby(self.name).count().max(1)[self.order]
        below_locs = (counts <= self.frequency_cutoff) | pd.isnull(counts)
        below = counts.loc[below_locs].index

        def _remap(x):
            if x in below:
                return np.nan
            else:
                return x

        index = map_[self.name].dropna().index
        map_[self.name] = map_.loc[index, self.name].apply(_remap)
        self._update_order(_remap)

    def convert_to_numeric(self, map_):
        """Converts the data to integer values

        Parameters
        ----------
        map_ : DataFrame
            A pandas object containing the data to be analyzed. The
            Question `name` should be a column in the `map_`.

        """
        self.check_map(map_)
        order = {g: i for i, g in enumerate(self.order)}

        def remap_(x):
            if isinstance(x, (int, float)):
                return x
            elif x in order:
                return order[x]
            else:
                return np.nan

        map_[self.name] = map_[self.name].apply(remap_)
        self._update_order(remap_)

    def label_order(self, map_):
        """Prefixes the data with an ordinal integer

        Parameters
        ----------
        map_ : DataFrame
            A pandas object containing the data to be analyzed. The
            Question `name` should be a column in the `map_`.

        """
        self.check_map(map_)
        order = {g: '(%i) %s' % (i, g) for i, g in enumerate(self.order)}

        def remap_(x):
            if isinstance(x, (int, float)):
                return x
            elif x in order:
                return order[x]
            else:
                return np.nan

        map_[self.name] = map_[self.name].apply(remap_)
        self._update_order(remap_)


class AgClinical(AgCategorical):
    def __init__(self, name, description, strict=True, **kwargs):
        """A question object for questions describing clincial conditions

        The American Gut addresses several clincial conditions, based on the
        diagnosis method.

        Parameters
        ----------
        name : str
            The name of column where the group is stored
        description : str
            A brief description of the data contained in the question
        strict : bool
            Whether clincial conditions should consider only those diagnosed
            by a conventional medicine practioner (`True`), or include
            alternative or self diagnosis.
        clean_name : str, optional
            A nicer version of the way the column should be named.
        remap : function
            A function used to remap the data in the question.
        """

        order = ['Diagnosed by a medical professional (doctor, physician '
                 'assistant)',
                 'Diagnosed by an alternative medicine practitioner',
                 'Self-diagnosed',
                 'I do not have this condition'
                 ]
        AgCategorical.__init__(self, name, description, str, order, **kwargs)
        self.strict = strict
        self.type = 'Clinical'
        self.extremes = ['I do not have this condition',
                         'Diagnosed by a medical professional (doctor, '
                         'physician assistant)']

    def remap_clinical(self, map_):
        """Converts the clincial condtion to a boolean value

        map_ : DataFrame
            A pandas dataframe containing the column described by the question
            name.

        """
        self.check_map(map_)
        if self.strict:
            def remap_(x):
                if x in {'Diagnosed by a medical professional (doctor, '
                         'physician assistant)', 'yes'}:
                    return 'Yes'
                elif x in {'Diagnosed by an alternative medicine practitioner',
                           'Self-diagnosed'}:
                    return np.nan
                elif x in {'I do not have this condition', 'No'}:
                    return 'No'
                else:
                    return np.nan
        else:
            def remap_(x):
                    if x in {'Diagnosed by a medical professional (doctor, '
                             'physician assistant)',
                             'Diagnosed by an alternative medicine '
                             'practitioner', 'Self-diagnosed', 'Yes'}:
                        return 'Yes'
                    elif x in {'I do not have this condition',
                               'No'}:
                        return 'No'
                    else:
                        return np.nan

        map_[self.name] = map_[self.name].apply(remap_)
        self._update_order(remap_)
        self.extremes = ['No', 'Yes']


class AgFrequency(AgCategorical):
    def __init__(self, name, description, order=None, combine=None, **kwargs):
        """
        A question object for frequency questions relating to weekly occurance

        Parameters
        ----------
        name : str
            The name of column where the group is stored
        description : str
            A brief description of the data contained in the question
        combine : {None, 'rarely', 'weekly'}
            Whether the frequency values should be combined. `None` will lead
            the data uneffected. `'rarely'` will combine infrequent (`'Never'`
            and `'Rarely (a few times/month)'`). `weekly` will combine
            occurances more than once a week
            (`'Occasionally (1-2 times/week)'`, `'Regularly (3-5 times/week)'`,
            and `'Daily'`).

        clean_name : str, optional
            A nicer version of the way the column should be named.
        remap : function
            A function used to remap the data in the question.
        """
        if order is None:
            order = ['Never',
                     'Rarely (a few times/month)',
                     'Occasionally (1-2 times/week)',
                     'Regularly (3-5 times/week)',
                     'Daily']
        AgCategorical.__init__(self, name, description, dtype=str, order=order,
                               **kwargs)
        self.type = 'Frequency'
        self.combine = combine

    def clean(self, map_):
        """Converts the frequency values to shorter strings

        map_ : DataFrame
            A pandas dataframe containing the column described by the question
            name.

        """
        self.check_map(map_)

        def remap_(x):
            if isinstance(x, str) and '(' in x:
                return x.split('(')[1].replace(')', '').capitalize()
            else:
                return x

        map_[self.name] = map_[self.name].apply(remap_)
        self._update_order(remap_)

    def combine_frequency(self, map_):
        """Combines frequency responses using `self.

         map_ : DataFrame
            A pandas dataframe containing the column described by the question
            name.
        """
        if self.combine not in {None, 'rarely', 'weekly'}:
            raise ValueError('%s is not a supported value for combine.'
                             % self.combine)

        if self.combine == 'rarely':
            def remap_(x):
                if x in {'Never', 'Rarely (a few times/month)',
                         'A few times/month',  'Less than once/week'}:
                    return 'Less than once/week'
                else:
                    return x
        elif self.combine == 'weekly':
            def remap_(x):
                if x in {'Occasionally (1-2 times/week)', '1-2 times/week',
                         'Regularly (3-5 times/week)', '3-5 times/week',
                         'Daily'}:
                    return "More than once/week"
                else:
                    return x
        elif self.combine is None:
            def remap_(x):
                return x

        map_[self.name] = map_[self.name].apply(remap_)
        self._update_order(remap_)
