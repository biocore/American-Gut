from os import remove, mkdir
from os.path import join, exists

from unittest import TestCase, main
import americangut.notebook_environment as agenv
from americangut.util import get_new_path
from americangut.per_category import cat_taxa_summaries


class PerCategoryTests(TestCase):
    def setUp(self):
        self.taxa_base = agenv.paths['populated-templates']['result-taxa']

        try:
            get_new_path(
                agenv.paths['populated-templates']['result-taxa'])
        except IOError:
            pass

        self.taxa_files = [
            'ag-stool-average.txt',
            'ag-stool-sex-male.txt',
            'ag-stool-sex-other.txt',
            'ag-stool-sex-female.txt',
            'ag-oral-diet-Omnivore.txt',
            'ag-oral-diet-Vegetarian.txt',
            'ag-oral-diet-Vegan.txt',
            'ag-oral-diet-Omnivore_but_do_not_eat_red_meat.txt',
            'ag-oral-diet-Vegetarian_but_eat_seafood.txt',
            'ag-skin-sex-male.txt',
            'ag-skin-sex-female.txt',
            'ag-stool-diet-Omnivore.txt',
            'ag-stool-diet-Vegetarian.txt',
            'ag-stool-diet-Vegan.txt',
            'ag-stool-diet-Omnivore_but_do_not_eat_red_meat.txt',
            'ag-stool-diet-Vegetarian_but_eat_seafood.txt',
            'ag-skin-hand-I_am_left_handed.txt',
            'ag-skin-hand-I_am_ambidextrous.txt',
            'ag-skin-hand-I_am_right_handed.txt',
            'ag-oral-sex-male.txt',
            'ag-oral-sex-female.txt',
            'ag-stool-bmi-Overweight.txt',
            'ag-stool-bmi-Obese.txt',
            'ag-stool-bmi-Underweight.txt',
            'ag-stool-bmi-Normal.txt',
            'ag-oral-flossing-Never.txt',
            'ag-oral-flossing-Rarely.txt',
            'ag-oral-flossing-Daily.txt',
            'ag-oral-flossing-Occasionally.txt',
            'ag-oral-flossing-Regularly.txt',
            'ag-stool-age-teen.txt',
            'ag-stool-age-20s.txt',
            'ag-stool-age-60s.txt',
            'ag-stool-age-30s.txt',
            'ag-stool-age-child.txt',
            'ag-stool-age-40s.txt',
            'ag-stool-age-70+.txt',
            'ag-stool-age-baby.txt',
            'ag-stool-age-50s.txt',
            'ag-oral-average.txt',
            'ag-skin-cosmetics-Never.txt',
            'ag-skin-cosmetics-Rarely.txt',
            'ag-skin-cosmetics-Daily.txt',
            'ag-skin-cosmetics-Occasionally.txt',
            'ag-skin-cosmetics-Regularly.txt',
            'ag-skin-age-teen.txt',
            'ag-skin-age-20s.txt',
            'ag-skin-age-60s.txt',
            'ag-skin-age-30s.txt',
            'ag-skin-age-70+.txt',
            'ag-skin-age-40s.txt',
            'ag-skin-age-child.txt',
            'ag-skin-age-baby.txt',
            'ag-skin-age-50s.txt',
            'ag-skin-average.txt',
            'ag-oral-age-teen.txt',
            'ag-oral-age-20s.txt',
            'ag-oral-age-60s.txt',
            'ag-oral-age-30s.txt',
            'ag-oral-age-child.txt',
            'ag-oral-age-40s.txt',
            'ag-oral-age-70+.txt',
            'ag-oral-age-50s.txt']

    def tearDown(self):
        for f in self.taxa_files:
            path = join(self.taxa_base, f)
            if exists(path):
                remove(path)

    def test_cat_taxa_summaries(self):
        cat_taxa_summaries()
        # Make sure all files created
        for f in self.taxa_files:
            path = join(self.taxa_base, f)
            if not exists(path):
                raise AssertionError('File %s not generated!' % f)

        # Test file for correct format
        self.maxDiff = None
        with open(join(self.taxa_base, 'ag-skin-age-teen.txt')) as f:
            obs = f.read()

        exp = ''
        self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()
