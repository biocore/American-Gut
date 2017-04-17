from unittest import TestCase, main

import numpy as np
import pandas as pd

from americangut.question.ag_question import AgQuestion

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


class AgQuestionTest(TestCase):
    def setUp(self):
        self.map_ = pd.DataFrame(
            amgut_sub,
            columns=["#SampleID", "AGE_YEARS", "SEX", "IBD",
                     "ALCOHOL_TYPES_BEERCIDER", "ALCOHOL_TYPES_RED_WINE",
                     "ALCOHOL_TYPES_SOUR_BEERS", "ALCOHOL_TYPES_UNSPECIFIED",
                     "BOWEL_MOVEMENT_FREQUENCY", "BOWEL_MOVEMENT_QUALITY",
                     "COLLECTION_TIMESTAMP", "ALCOHOL_FREQUENCY"],
            ).set_index('#SampleID')
        self.map_.replace('', np.nan, inplace=True)

        self.name = 'TEST_COLUMN'
        self.description = ('"Say something, anything."\n'
                            '"Testing... 1, 2, 3"\n'
                            '"Anything but that!"')
        self.dtype = str

        self.fun = lambda x: 'foo'

        self.ag_question = AgQuestion(
            name='ALCOHOL_TYPES_UNSPECIFIED',
            description=('Has the participant described their alcohol use'),
            dtype=bool
            )

    def test_init(self):
        test = AgQuestion(self.name, self.description, self.dtype)

        self.assertEqual(test.name, self.name)
        self.assertEqual(test.description, self.description)
        self.assertEqual(test.dtype, self.dtype)
        self.assertEqual(test.clean_name, 'Test Column')
        self.assertEqual(test.free_response, False)
        self.assertEqual(test.mimarks, False)
        self.assertEqual(test.ontology, None)
        self.assertEqual(test.remap_, None)

    def test_init_kwargs(self):
        test = AgQuestion(self.name, self.description, self.dtype,
                          clean_name='This is a test.',
                          free_response=True,
                          mimarks=True,
                          ontology='GAZ',
                          remap=self.fun,
                          )

        self.assertEqual(test.name, self.name)
        self.assertEqual(test.description, self.description)
        self.assertEqual(test.dtype, self.dtype)
        self.assertEqual(test.clean_name, 'This is a test.')
        self.assertEqual(test.free_response, True)
        self.assertEqual(test.mimarks, True)
        self.assertEqual(test.ontology, 'GAZ')
        self.assertEqual(test.remap_(1), 'foo')

    def test_init_name_error(self):
        with self.assertRaises(TypeError):
            AgQuestion(1, self.description, self.dtype)

    def test_init_description_error(self):
        with self.assertRaises(TypeError):
            AgQuestion(self.name, 3, self.dtype)

    def test_init_class_error(self):
        with self.assertRaises(TypeError):
            AgQuestion(self.name, self.description, 'bool')

    def test_init_clean_name_error(self):
        with self.assertRaises(TypeError):
            AgQuestion(self.name, self.description, self.dtype, ['Cats'])

    def test_remap_data_type(self):
        self.ag_question.remap_data_type(self.map_)
        self.assertEquals(set(self.map_[self.ag_question.name]) - {np.nan},
                          {True, False})
        self.ag_question.dtype = int
        self.ag_question.remap_data_type(self.map_)
        self.assertEqual(set(self.map_[self.ag_question.name]) - {np.nan},
                         {0, 1})

    def test_remap_data_type_raises(self):
        self.map_.loc["10317.000006668", "ALCOHOL_TYPES_UNSPECIFIED"] = 'foo'
        with self.assertRaises(TypeError):
            self.ag_question.remap_data_type(self.map_)

    def test_check_map(self):
        self.ag_question.name = 'TEST_CASE'
        with self.assertRaises(ValueError):
            self.ag_question.check_map(self.map_)

if __name__ == '__main__':
    main()
