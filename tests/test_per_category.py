from os import makedirs, listdir
from os.path import join, dirname, realpath, exists, abspath
from shutil import rmtree

from unittest import TestCase, main
import americangut.notebook_environment as agenv
from americangut.util import get_existing_path
from americangut.per_category import cat_taxa_summaries, cat_alpha_plots


class PerCategoryTests(TestCase):
    def setUp(self):
        # Make expected directory structure
        currpath = dirname(realpath(__file__))
        self.dirpath = abspath(join(currpath,
                               '../agp_processing/10-populated-templates'))
        makedirs(self.dirpath)
        makedirs(join(self.dirpath, 'taxa'))
        makedirs(join(self.dirpath, 'alpha-div'))

        self.alpha_path = agenv.paths['collapsed']['100nt']['alpha-map']
        agenv.paths['collapsed']['100nt']['alpha-map'] = \
            '../tests/data/ag_testing/cat_metadata_test.txt'

        self.path = agenv.paths['collapsed']['notrim']['1k']
        agenv.paths['collapsed']['notrim']['1k'] = \
            {'ag-biom':
                '../tests/data/ag_testing/category_test.biom',
             'ag-fecal':
                '../tests/data/ag_testing/ag-fecal.biom',
             'ag-oral-flossing':
                '../tests/data/ag_testing/ag-oral-flossing.biom'}

    def tearDown(self):
        rmtree(self.dirpath, ignore_errors=True)
        agenv.paths['collapsed']['notrim']['1k'] = self.path
        agenv.paths['collapsed']['100nt']['alpha-map'] = self.alpha_path

    def test_cat_taxa_summaries(self):
        cat_taxa_summaries()
        path = get_existing_path(
            '../agp_processing/10-populated-templates/taxa')
        exp = ['ag-oral-flossing-Daily.txt', 'ag-oral-flossing-Never.txt',
               'ag-oral-flossing-Occasionally.txt',
               'ag-oral-flossing-Rarely.txt', 'ag-oral-flossing-Regularly.txt',
               'ag-stool-average.txt']
        files = listdir(path)
        self.assertItemsEqual(files, exp)

        with open(join(path, 'ag-oral-flossing-Rarely.txt')) as f:
            obs = f.read()
        exp = """k__Bacteria; p__Fusobacteria; c__Fusobacteriia; o__Fusobacteriales; f__Fusobacteriaceae; g__Fusobacterium\t0.0110430107527
k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Planococcaceae; g__\t0.00196774193548
k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__[Weeksellaceae]; g__Chryseobacterium\t0.0125913978495
k__Bacteria; p__Bacteroidetes; c__Sphingobacteriia; o__Sphingobacteriales; f__Sphingobacteriaceae; g__Sphingobacterium\t0.00608602150538
k__Bacteria; p__Firmicutes; c__Erysipelotrichi; o__Erysipelotrichales; f__Erysipelotrichaceae; g__Bulleidia\t0.00269892473118
k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Neisseriales; f__Neisseriaceae; g__Neisseria\t0.074935483871
k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Oribacterium\t0.00749462365591
k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Xanthomonadales; f__Xanthomonadaceae; g__\t0.017311827957
k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Moraxellaceae; g__Enhydrobacter\t0.000150537634409
k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__\t0.00667741935484
k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Caulobacterales; f__Caulobacteraceae; g__Brevundimonas\t0.00996774193548
k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__[Weeksellaceae]; g__\t0.00139784946237
k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacteriales; f__Enterobacteriaceae; g__Klebsiella\t0.00174193548387
k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus\t0.00756989247312
k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Prevotellaceae; g__Prevotella\t0.0509569892473
k__Bacteria; p__Fusobacteria; c__Fusobacteriia; o__Fusobacteriales; f__Leptotrichiaceae; g__Leptotrichia\t0.0143333333333
k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus\t0.138817204301
k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Pseudomonadaceae; g__\t0.0051935483871
k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__[Paraprevotellaceae]; g__[Prevotella]\t0.00132258064516
k__Bacteria; p__Proteobacteria; c__Epsilonproteobacteria; o__Campylobacterales; f__Campylobacteraceae; g__Campylobacter\t0.00133333333333
k__Bacteria; p__Firmicutes; c__Bacilli; o__Gemellales; f__Gemellaceae; g__\t0.0097311827957
k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacteriales; f__Enterobacteriaceae; g__\t0.0135268817204
k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pasteurellales; f__Pasteurellaceae; g__Haemophilus\t0.0668279569892
k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Brucellaceae; g__Ochrobactrum\t0.00575268817204
k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Veillonellaceae; g__Veillonella\t0.0709247311828
k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Porphyromonadaceae; g__Porphyromonas\t0.0165698924731
k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pasteurellales; f__Pasteurellaceae; g__Aggregatibacter\t0.00410752688172
k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Neisseriales; f__Neisseriaceae; g__\t0.00394623655914
k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Moraxellaceae; g__Acinetobacter\t0.0247634408602
k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Aerococcaceae; g__Alloiococcus\t0.000225806451613
k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Moryella\t0.000612903225806
k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Micrococcaceae; g__Rothia\t0.114602150538
k__Bacteria; p__Actinobacteria; c__Coriobacteriia; o__Coriobacteriales; f__Coriobacteriaceae; g__Atopobium\t0.00761290322581
k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Comamonadaceae; g__Alicycliphilus\t0.00530107526882
k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Actinomycetaceae; g__Actinomyces\t0.0204408602151
k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Xanthomonadales; f__Xanthomonadaceae; g__Stenotrophomonas\t0.0443655913978
k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Corynebacteriaceae; g__Corynebacterium\t0.000440860215054
k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacteriales; f__Enterobacteriaceae; g__Morganella\t0.000161290322581
k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Carnobacteriaceae; g__Granulicatella\t0.025247311828
k__Bacteria; p__SR1; c__; o__; f__; g__\t0.00222580645161
k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Pseudomonadaceae; g__Pseudomonas\t0.0665268817204
"""
        self.assertEqual(obs, exp)

    def test_cat_alpha_plots(self):
        cat_alpha_plots()

        self.assertTrue(exists(join(self.dirpath, 'alpha-div',
                        'pd_oral-age-50s.png')))
        self.assertTrue(exists(join(self.dirpath, 'alpha-div',
                        'pd_stool-sex-male.png')))
        self.assertTrue(exists(join(self.dirpath, 'alpha-div',
                        'shannon_skin-dominant_hand-I_am_left_handed.png')))


if __name__ == '__main__':
    main()
