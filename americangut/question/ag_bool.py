import numpy as np

from americangut.question.ag_categorical import AgCategorical


class AgBool(AgCategorical):
    def __init__(self, name, description, **kwargs):
        """A question object for boolean question

        Parameters
        ----------
        name : str
            The name of column where the group is stored
        description : str
            A brief description of the data contained in the question
        clean_name : str, optional
            A nicer version of the way the column should be named.
        remap : function
            A function used to remap the data in the question.
        mimmarks : bool, optional
            If the question was a mimmarks standard field
        ontology : str, optional
            The type of ontology, if any, which was used in the field value.

        """
        AgCategorical.__init__(self, name, description, bool,
                               ['true', 'false'], **kwargs)
        self.type = 'Bool'

    def convert_to_word(self, map_):
        """Converts boolean values to 'yes' and 'no'

         map_ : DataFrame
            A pandas dataframe containing the column described by the question
            name.
        """
        self.check_map(map_)
        self.remap_data_type(map_, watch=False)

        def remap_(x):
            if x in self.true_values:
                return 'yes'
            elif x in self.false_values:
                return 'no'
            else:
                return np.nan

        map_[self.name] = map_[self.name].apply(remap_)
        self._update_order(remap_)


class AgMultiple(AgBool):
    def __init__(self, name, description, unspecified, **kwargs):
        """A question object for a multiple response question

        Multiple response questions have a checkbox, where each response
        can be coded as a boolean. One of the checkboxes, the unspecified
        column, is used to indicate the participant did not respond to the
        question block.

        Parameters
        ----------
        name : str
            The name of column where the group is stored
        description : str
            A brief description of the data contained in the question
        unspecified : str
            The name of a reference column in the map_ which was used to
            indicate the participant did not respond. This should be an AgBool
            type question.
        clean_name : str, optional
            A nicer version of the way the column should be named.
        remap : function
            A function used to remap the data in the question.
        mimmarks : bool, optional
            If the question was a mimmarks standard field
        ontology : str, optional
            The type of ontology, if any, which was used in the field value.

        """

        AgBool.__init__(self, name, description, **kwargs)
        self.type = 'Multiple'
        self.unspecified = unspecified

    def correct_unspecified(self, map_):
        """Corrects column for missing responses

        map_ : DataFrame
            A pandas dataframe containing the column described by the question
            name.
        """
        self.check_map(map_)
        if self.unspecified not in map_.columns:
            raise ValueError('The unspecified reference (%s) is not a column'
                             ' in the metadata!' % self.unspecified)

        def remap_(x):
            if x in self.true_values:
                return True
            elif x in self.false_values:
                return False
            else:
                return np.nan

        map_[self.unspecified] = \
            map_[self.unspecified].apply(remap_).astype(bool)

        map_.loc[map_[self.unspecified], self.name] = np.nan
