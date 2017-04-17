import numpy as np

from americangut.question.ag_question import AgQuestion


class AgContinous(AgQuestion):

    def __init__(self, name, description, unit=None, **kwargs):
        """A Question object with continous responses

        Parameters
        ----------
        name : str
            The name of column where the group is stored
        description : str
            A brief description of the data contained in the question
        unit : str, optional
            The unit of measure for the data type, if relevant. Units
            are awesome, and should be included whenever relevant.
        clean_name : str, optional
            A nicer version of the way the column should be named.
        remap : function
            A function used to remap the data in the question.
        range : two-elemant iterable
            The range of values pertinant to analysis.
        """

        AgQuestion.__init__(self, name, description, float, **kwargs)
        self.unit = unit
        lower, upper = kwargs.get('range', [None, None])
        if lower > upper:
            raise ValueError('The lower limit cannot be greater than '
                             'the upper!')
        self.lower = lower
        self.upper = upper
        self.type = 'Continous'

    def drop_outliers(self, map_):
        """Removes datapoints outside of the `range`"""
        self.check_map(map_)
        self.remap_data_type(map_, watch=False)

        if self.lower is not None:

            def remap_(x):
                if x < self.lower:
                    return np.nan
                elif x > self.upper:
                    return np.nan
                else:
                    return x

            map_[self.name] = map_[self.name].apply(remap_)
