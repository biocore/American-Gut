from unittest import TestCase, main

import numpy as np
import pandas as pd
import pandas.util.testing as pdt

from americangut.question.ag_categorical import (AgCategorical,
                                                 AgClinical,
                                                 AgFrequency,
                                                 )


amgut_sub = np.array([
    ["10317.000006668", "71.0", "male", "I do not have this condition",
     "false", "false", "false", "true", "", "", "07/09/2013 08:30", "Never"],
    ["10317.000030344", "40.0", "female", "I do not have this condition",
     "false", "true", "false", "false", "One", "I tend to have normal formed "
     "stool", "09/08/2015 09:10", "Regularly (3-5 times/week)"],
    ["10317.000031833", "76.0", "female", "I do not have this condition",
     "true", "true", "false", "false", "One", "I tend to have normal formed "
     "stool", "08/17/2015 08:33", "Daily"],
    ["10317.000013134", "45.0", "male", "", "false", "false", "false", "true",
     "", "", "01/12/2014 00:35", "Regularly (3-5 times/week)"],
    ["10317.000022634", "40.0", "male", "I do not have this condition",
     "false", "false", "false", "true", "Less than one", "I tend to be "
     "constipated (have difficulty passing stool)", "02/02/2015 09:00",
     "Never"],
    ["10317.000020683", "34.0", "male", "I do not have this condition",
     "true", "true", "false", "false", "Less than one", "I don't know, "
     "I do not have a point of reference", "11/08/2014 10:15",
     "Regularly (3-5 times/week)"],
    ["10317.000022916", "77.0", "male", "I do not have this condition", "true",
     "false", "false", "false", "One", "I don't know, I do not have a point of"
     " reference", "01/13/2015 19:14", "Rarely (a few times/month)"],
    ["10317.000031363", "45.0", "male", "I do not have this condition", "true",
     "false", "false", "false", "One", "I tend to be constipated (have "
     "difficulty passing stool)", "08/07/2015 08:45",
     "Occasionally (1-2 times/week)"],
    ["10317.000006952", "33.0", "female", "I do not have this condition",
     "false", "false", "false", "true", "", "", "07/03/2013 12:45",
     "Rarely (a few times/month)"],
    ["10317.000027363", "52.0", "male", "I do not have this condition",
     "true", "false", "false", "false", "One", "", "07/05/2015 11:35",
     "Rarely (a few times/month)"],
    ["10317.000002337", "57.0", "female", "I do not have this condition",
     "false", "false", "false", "true", "", "", "05/08/2013 09:10", "Daily"],
    ["10317.000027574", "6.0", "female", "I do not have this condition",
     "false", "false", "false", "true", "One", "I tend to have normal formed "
     "stool", "06/17/2015 19:15", "Never"],
    ["10317.000003365", "26.0", "female", "I do not have this condition",
     "false", "false", "false", "true", "", "", "04/02/2013 10:30", "Regularly"
     " (3-5 times/week)"],
    ["10317.000002347", "45.0", "female", "I do not have this condition",
     "false", "false", "false", "true", "", "", "05/04/2013 14:40",
     "Occasionally (1-2 times/week)"],
    ["10317.000011322", "", "female", "I do not have this condition", "false",
     "false", "false", "true", "", "", "01/14/2014 09:08", "Regularly "
     "(3-5 times/week)"],
    ["10317.000020053", "35.0", "male", "Self-diagnosed", "true", "true",
     "false", "false", "Three", "I tend to have normal formed stool",
     "01/03/2015 17:15", "Occasionally (1-2 times/week)"],
    ["10317.000029572", "65.0", "male", "Diagnosed by a medical professional "
     "(doctor, physician assistant)", "false", "false", "false", "false",
     "Less than one", "I tend to be constipated (have difficulty passing "
     "stool)", "08/28/2015 14:00", "Daily"],
    ["10317.000023046", "14.0", "female", "Diagnosed by a medical professional"
     " (doctor, physician assistant)", "false", "false", "false", "true",
     "Three", "I tend to have normal formed stool", "01/27/2015 10:15",
     "Never"],
    ["10317.000028154", "54.0", "female", "Diagnosed by a medical professional"
     " (doctor, physician assistant)", "true", "true", "false", "false",
     "Five or more", "I tend to have diarrhea (watery stool)",
     "06/25/2015 21:25", "Occasionally (1-2 times/week)"],
    ["10317.000020499", "33.0", "female", "Diagnosed by a medical professional"
     " (doctor, physician assistant)", "true", "false", "false", "false",
     "Five or more", "I tend to have diarrhea (watery stool)",
     "12/02/2014 09:00", "Rarely (a few times/month)"],
    ])


class AgCategoricalTest(TestCase):

    def setUp(self):
        self.map_ = pd.DataFrame(
            amgut_sub,
            columns=["#SampleID", "AGE_YEARS", "SEX", "IBD",
                     "ALCOHOL_TYPES_BEERCIDER", "ALCOHOL_TYPES_RED_WINE",
                     "ALCOHOL_TYPES_SOUR_BEERS",
                     "ALCOHOL_TYPES_UNSPECIFIED",
                     "BOWEL_MOVEMENT_FREQUENCY", "BOWEL_MOVEMENT_QUALITY",
                     "COLLECTION_TIMESTAMP", "ALCOHOL_FREQUENCY"],
            ).set_index('#SampleID')
        self.map_.replace('', np.nan, inplace=True)

        self.name = 'BOWEL_MOVEMENT_QUALITY'
        self.description = ('Whether the participant is constipated, has '
                            'diarrhea or is regular.')
        self.dtype = str
        self.order = ["I don't know, I do not have a point of reference",
                      "I tend to have normal formed stool",
                      "I tend to be constipated (have difficulty passing "
                      "stool)",
                      "I tend to have diarrhea (watery stool)",
                      ]
        self.extremes = ["I tend to have normal formed stool",
                         "I tend to have diarrhea (watery stool)"]

        def remap_poo_quality(x):
            if x in {"I don't know, I do not have a point of reference"}:
                return 'No Reference'
            elif x == 'I tend to have normal formed stool':
                return 'Well formed'
            elif x == ('I tend to be constipated (have difficulty passing'
                       ' stool)'):
                return 'Constipated'
            elif x == 'I tend to have diarrhea (watery stool)':
                return 'Diarrhea'
            else:
                return x

        self.remap_ = remap_poo_quality

        self.ag_categorical = AgCategorical(
            name=self.name,
            description=self.description,
            dtype=self.dtype,
            order=self.order,
            extremes=self.extremes,
            remap=remap_poo_quality
            )

    def test_ag_categorical_init(self):
        self.assertEqual('No Reference', self.remap_(
            "I don't know, I do not have a point of reference"
            ))
        self.assertEqual(3, self.remap_(3))

        test = AgCategorical(self.name,
                             self.description,
                             self.dtype,
                             self.order,
                             remap=self.remap_)

        self.assertEqual(self.order, test.order)
        self.assertEqual([], test.earlier_order)
        self.assertEqual(test.extremes,
                         ["I don't know, I do not have a point of reference",
                          "I tend to have diarrhea (watery stool)"
                          ])
        self.assertEqual('No Reference', test.remap_(
            "I don't know, I do not have a point of reference"
            ))
        self.assertEqual(3, test.remap_(3))
        self.assertEqual(test.type, 'Categorical')

    def test_ag_init_error(self):
        with self.assertRaises(ValueError):
            AgCategorical(self.name, self.description, list, self.order)

    def test_ag_categorical_update_groups(self):
        self.ag_categorical._update_order(self.ag_categorical.remap_)
        self.assertEqual(self.ag_categorical.order,
                         ['No Reference', 'Well formed', 'Constipated',
                          'Diarrhea'])
        self.assertEqual(self.ag_categorical.extremes,
                         ['Well formed', 'Diarrhea'])

        self.assertEqual(self.ag_categorical.earlier_order, [self.order])

    def test_ag_categorical_remove_ambigious(self):
        self.ag_categorical.drop_ambiguous = True
        self.ag_categorical.remove_ambiguity(self.map_)

        self.assertEqual(set(self.map_[self.ag_categorical.name]) - {np.nan}, {
            "I tend to be constipated (have difficulty passing stool)",
            "I tend to have diarrhea (watery stool)",
            "I tend to have normal formed stool",
            })
        self.assertEqual(self.ag_categorical.earlier_order, [self.order])

    def test_ag_categorical_remap_groups(self):
        self.ag_categorical.remap_groups(self.map_)

        self.assertEqual(set(self.map_[self.ag_categorical.name]) - {np.nan},
                         {'No Reference', 'Well formed', 'Constipated',
                          'Diarrhea'})

        self.assertEqual(self.ag_categorical.earlier_order, [self.order])

    def test_ag_categorical_drop_infrequent(self):
        order = [
            "I don't know, I do not have a point of reference",
            "I tend to have normal formed stool",
            "I tend to be constipated (have difficulty passing stool)",
            "I tend to have diarrhea (watery stool)",
            "Solid as a Rock"]
        self.ag_categorical.order = order
        self.ag_categorical.frequency_cutoff = 4
        self.ag_categorical.drop_infrequent(self.map_)

        self.assertEqual(set(self.map_[self.ag_categorical.name]) - {np.nan}, {
            "I tend to have normal formed stool"
            })
        self.assertEqual(self.ag_categorical.earlier_order, [order])
        self.assertEqual(self.ag_categorical.order,
                         ["I tend to have normal formed stool"])


class AgFrequencyTest(TestCase):

    def setUp(self):
        self.map_ = pd.DataFrame(
            amgut_sub,
            columns=["#SampleID", "AGE_YEARS", "SEX", "IBD",
                     "ALCOHOL_TYPES_BEERCIDER", "ALCOHOL_TYPES_RED_WINE",
                     "ALCOHOL_TYPES_SOUR_BEERS",
                     "ALCOHOL_TYPES_UNSPECIFIED",
                     "BOWEL_MOVEMENT_FREQUENCY", "BOWEL_MOVEMENT_QUALITY",
                     "COLLECTION_TIMESTAMP", "ALCOHOL_FREQUENCY"],
            ).set_index('#SampleID')
        self.map_.replace('', np.nan, inplace=True)
        self.name = 'ALCOHOL_FREQUENCY'
        self.description = 'How often the participant consumes alcohol.'

        self.ag_frequency = AgFrequency(
            name=self.name,
            description=self.description,
            )

    def test_ag_frequency_init(self):
        order = ['Never', 'Rarely (a few times/month)',
                 'Occasionally (1-2 times/week)', 'Regularly (3-5 times/week)',
                 'Daily']
        test = AgFrequency(
            name='ALCOHOL_FREQUENCY',
            description=('How often the participant consumes alcohol.'),
            )

        self.assertEqual(test.order, order)
        self.assertEqual(test.dtype, str)
        self.assertEqual(test.type, 'Frequency')

    def test_ag_frequency_init_order(self):
        order = ['Never', 'Rarely (less than once/week)',
                 'Occasionally (1-2 times/week)', 'Regularly (3-5 times/week)',
                 'Daily']
        test = AgFrequency(
            name='ALCOHOL_FREQUENCY',
            description=('How often the participant consumes alcohol.'),
            order=order,
            )

        self.assertEqual(test.order, order)
        self.assertEqual(test.dtype, str)
        self.assertEqual(test.type, 'Frequency')

    def test_ag_frequency_clean(self):
        self.ag_frequency.clean(self.map_)
        self.assertEqual(set(self.map_[self.ag_frequency.name]) - {np.nan}, {
            'Never', 'A few times/month', '1-2 times/week', '3-5 times/week',
            'Daily'
            })
        self.assertEqual(
            self.ag_frequency.earlier_order,
            [['Never', 'Rarely (a few times/month)',
              'Occasionally (1-2 times/week)', 'Regularly (3-5 times/week)',
              'Daily']]
            )

    def test_ag_frequency_combine_frequency_rarely(self):
        self.ag_frequency.combine = 'rarely'
        self.ag_frequency.combine_frequency(self.map_)

        self.assertEqual(
            set(self.map_[self.ag_frequency.name]) - {np.nan},
            {'Less than once/week', 'Occasionally (1-2 times/week)',
             'Regularly (3-5 times/week)', 'Daily'}
            )
        self.assertEqual(
            self.ag_frequency.earlier_order,
            [['Never', 'Rarely (a few times/month)',
              'Occasionally (1-2 times/week)', 'Regularly (3-5 times/week)',
              'Daily']]
            )

    def test_ag_frequency_combine_frequency_weekly(self):
        self.ag_frequency.combine = 'weekly'
        self.ag_frequency.combine_frequency(self.map_)

        self.assertEqual(
            set(self.map_[self.ag_frequency.name]) - {np.nan},
            {'Never', 'Rarely (a few times/month)', "More than once/week"}
            )
        self.assertEqual(self.ag_frequency.order,
                         ['Never', 'Rarely (a few times/month)',
                          "More than once/week"])
        self.assertEqual(
            self.ag_frequency.earlier_order,
            [['Never', 'Rarely (a few times/month)',
              'Occasionally (1-2 times/week)', 'Regularly (3-5 times/week)',
              'Daily']]
            )

    def test_ag_frequency_combine_frequency_error(self):
        self.ag_frequency.combine = 'foo'
        with self.assertRaises(ValueError):
            self.ag_frequency.combine_frequency(self.map_)


class AgClinicalTest(TestCase):

    def setUp(self):
        self.map_ = pd.DataFrame(
            amgut_sub,
            columns=["#SampleID", "AGE_YEARS", "SEX", "IBD",
                     "ALCOHOL_TYPES_BEERCIDER", "ALCOHOL_TYPES_RED_WINE",
                     "ALCOHOL_TYPES_SOUR_BEERS",
                     "ALCOHOL_TYPES_UNSPECIFIED",
                     "BOWEL_MOVEMENT_FREQUENCY", "BOWEL_MOVEMENT_QUALITY",
                     "COLLECTION_TIMESTAMP", "ALCOHOL_FREQUENCY"],
            ).set_index('#SampleID')
        self.map_.replace('', np.nan, inplace=True)

        self.name = 'IBD'
        self.description = ('Has the participant been diagnosed with IBD?')

        self.ag_clinical = AgClinical(
            name=self.name,
            description=self.description,
            )

    def test_ag_clinical_init(self):
        clinical_groups = ['Diagnosed by a medical professional '
                           '(doctor, physician assistant)',
                           'Diagnosed by an alternative medicine practitioner',
                           'Self-diagnosed',
                           'I do not have this condition']
        test_class = AgClinical(
            name='TEST_COLUMN',
            description='Testing. 1, 2, 3.\nAnything but that!',
            strict=True,
            )
        self.assertEqual(test_class.dtype, str)
        self.assertEqual(test_class.order, clinical_groups)
        self.assertEqual(test_class.type, 'Clinical')

    def test_ag_clinical_remap_strict(self):
        self.ag_clinical.remap_clinical(self.map_)

        count = pd.Series([14, 4],
                          index=pd.Index(['No', 'Yes'], name='IBD'),
                          dtype=np.int64)
        pdt.assert_series_equal(count, self.map_.groupby('IBD').count().max(1))

        self.assertEqual(
            self.ag_clinical.earlier_order,
            [['Diagnosed by a medical professional (doctor, physician '
              'assistant)',
              'Diagnosed by an alternative medicine practitioner',
              'Self-diagnosed',
              'I do not have this condition']]
            )

    def test_ag_clincial_remap(self):
        self.ag_clinical.strict = False
        self.ag_clinical.remap_clinical(self.map_)

        count = pd.Series([14, 5],
                          index=pd.Index(['No', 'Yes'], name='IBD'),
                          dtype=np.int64)
        pdt.assert_series_equal(count, self.map_.groupby('IBD').count().max(1))

        self.assertEqual(
            self.ag_clinical.earlier_order,
            [['Diagnosed by a medical professional (doctor, physician '
              'assistant)',
              'Diagnosed by an alternative medicine practitioner',
              'Self-diagnosed',
              'I do not have this condition']]
            )

if __name__ == '__main__':
    main()
