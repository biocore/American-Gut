import inspect

import numpy as np


class AgQuestion:
    """A base object class for handling American Gut Data dictionary entries"""
    true_values = {'yes', 'true', 1, 1.0, True}
    false_values = {'no', 'false', 0, 0.0, False}

    def __init__(self, name, description, dtype, clean_name=None, remap=None,
                 free_response=False, mimarks=False, ontology=None):
        """A base object for describing single question outputs

        Parameters
        ----------
        name : str
            The name of column where the group is stored
        description : str
            A brief description of the data contained in the question
        dtype : object
            The datatype in which the responses should be represented. (i.e.
            float, int, str).
        clean_name : str, optional
            A nicer version of the way the column should be named.
        remap : function
            A function used to remap the data in the question.
        free_response: bool, optional
            Whether the question is a free response question or controlled
            vocabulary
        mimarks : bool, optional
            If the question was a mimarks standard field
        ontology : str, optional
            The type of ontology, if any, which was used in the field value.

        """

        # Checks the arguments
        if not isinstance(name, str):
            raise TypeError('name must be a string.')
        if not isinstance(description, str):
            raise TypeError('description must be a string')
        if not inspect.isclass(dtype):
            raise TypeError('dtype must be a class')
        if not isinstance(clean_name, str) and clean_name is not None:
            raise TypeError('If supplied, clean_name must be a string')

        # Handles the main information about the data
        self.name = name
        self.description = description
        self.dtype = dtype

        self.type = 'Question'
        if clean_name is None:
            self.clean_name = name.replace('_', ' ').title()
        else:
            self.clean_name = clean_name

        # Sets up
        self.free_response = free_response
        self.mimarks = mimarks
        self.ontology = ontology
        self.remap_ = remap

    def check_map(self, map_):
        """Checks the group exists in the metadata

        Parameters
        ----------
        map_ : DataFrame
            A pandas object containing the data to be analyzed. The
            Question `name` should be a column in the `map_`.

        """
        if self.name not in map_.columns:
            raise ValueError('%s is not a column in the supplied map!'
                             % self.name)

    def remap_data_type(self, map_, watch=True):
        """Makes sure the target column in map_ has the correct datatype

        map_ : DataFrame
            A pandas dataframe containing the column described by the question
            name.
        watch : bool
            A flag to indicate the change should be logged. Should generally
            be true, unless `remap_data_type` is called before another
            function is executed.

        """
        if self.dtype == bool:
            if not (set(map_[self.name].dropna()).issubset(
                    self.true_values.union(self.false_values))):
                raise TypeError('%s cannot be cast to a bool value.'
                                % self.name)

            def remap_(x):
                if isinstance(x, str) and x.lower() in self.true_values:
                    return True
                elif isinstance(x, str) and x.lower() in self.false_values:
                    return False
                elif np.isnan(x):
                    return x
                else:
                    try:
                        return bool(x)
                    except:
                        return np.nan
        else:
            def remap_(x):
                return self.dtype(x)

        map_[self.name] = map_[self.name].apply(remap_).astype(self.dtype)
        map_.replace('nan', np.nan, inplace=True)

        if watch and self.type in {'Categorical', 'Multiple', 'Clinical',
                                   'Bool', 'Frequency'}:
            self._update_order(remap_)
