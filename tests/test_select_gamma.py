#!/usr/bin/env python

from unittest import TestCase, main
import select_gamma
from biom.parse import parse_biom_table

__author__ = "Amnon Amir"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Amnon Amir"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Amnon Amir"
__email__ = "amnon.amir@colorado.edu"


class SelectGammaTests(TestCase):
    def setUp(self):
        """ prepare the biom table and the expected output for the tests"""
        self.biom_table=biom_table
        self.expected_output=expected_output_gamma_3
        self.bt=parse_biom_table(self.biom_table)

    def test_select_gamma_complex(self):
        """ test the complex situation (10 gamma out of ~50, with 8 above cumulative threhold 0.03)"""
        # test the default parameter run situation - gammaproteobacteria, 3% cumulative freq
        output=select_gamma.get_high_freq_otus(self.bt,2,'c__Gammaproteobacteria',0.03)
        self.assertEqual(output,self.expected_output)

    def test_select_gamma_threshold(self):
        """ test the effect of the cumulative frequency threhold parameter"""
        # test what happens if the cumulative threshold gives 1 bacteria
        output=select_gamma.get_high_freq_otus(self.bt,2,'c__Gammaproteobacteria',0.13)
        self.assertEqual(output,[u'104\t0.02\t0.144444444444'])

        # test what happens if the cumulative threshold is too high
        output=select_gamma.get_high_freq_otus(self.bt,2,'c__Gammaproteobacteria',0.30)
        self.assertEqual(output,[])

        # test that we get all bacteria with a threhold of 0
        output=select_gamma.get_high_freq_otus(self.bt,2,'c__Gammaproteobacteria',0)
        self.assertEqual(len(output),10)

    def test_select_gamma_taxonomy(self):
        """ test the effect of the taxonomy string and taxonomy class"""
        # test that we filter the taxonomy correctly
        output=select_gamma.get_high_freq_otus(self.bt,2,'c__Betaproteobacteria',0)
        self.assertEqual(len(output),3)

        # test that we filter the taxonomy correctly
        output=select_gamma.get_high_freq_otus(self.bt,2,'c__Alphaproteobacteria',0)
        self.assertEqual(len(output),10)

        # test that we filter using the taxonomy class correctly
        output=select_gamma.get_high_freq_otus(self.bt,1,'p__Proteobacteria',0)
        self.assertEqual(len(output),26)

        # test that we filter using the taxonomy class correctly
        output=select_gamma.get_high_freq_otus(self.bt,0,'k__Bacteria',0)
        self.assertEqual(len(output),52)





biom_table='{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "BIOM-Format 1.1.2","date": "2013-11-01T17:53:20.390965","matrix_type": "sparse","matrix_element_type": "float","shape": [52, 3],"data": [[0,0,10.0],[0,1,30.0],[0,2,60.0],[1,0,20.0],[1,1,30.0],[1,2,70.0],[2,0,30.0],[2,1,30.0],[2,2,80.0],[3,0,40.0],[3,1,30.0],[3,2,90.0],[4,0,50.0],[4,1,30.0],[4,2,100.0],[5,0,40.0],[5,1,30.0],[5,2,90.0],[6,0,30.0],[6,1,30.0],[6,2,80.0],[7,0,20.0],[7,1,30.0],[7,2,70.0],[8,0,10.0],[8,1,30.0],[8,2,60.0],[9,1,30.0],[9,2,50.0],[10,0,45.0],[10,2,98.0],[11,0,49.0],[11,1,42.0],[11,2,79.0],[12,0,70.0],[12,1,73.0],[12,2,43.0],[13,0,46.0],[13,1,58.0],[13,2,55.0],[14,0,12.0],[14,1,47.0],[14,2,83.0],[15,0,63.0],[15,1,69.0],[15,2,79.0],[16,0,17.0],[16,1,59.0],[16,2,45.0],[17,0,58.0],[17,1,26.0],[17,2,65.0],[18,0,76.0],[18,1,2.0],[18,2,3.0],[19,0,58.0],[19,1,59.0],[19,2,98.0],[20,0,38.0],[20,1,41.0],[20,2,46.0],[21,0,57.0],[21,1,18.0],[21,2,57.0],[22,0,76.0],[22,1,79.0],[22,2,99.0],[23,0,26.0],[23,1,47.0],[23,2,75.0],[24,0,26.0],[24,1,37.0],[24,2,88.0],[25,0,57.0],[25,1,92.0],[25,2,86.0],[26,0,40.0],[26,1,84.0],[26,2,50.0],[27,0,99.0],[27,1,5.0],[27,2,49.0],[28,0,54.0],[28,1,12.0],[28,2,73.0],[29,0,45.0],[29,1,68.0],[29,2,46.0],[30,0,95.0],[30,1,56.0],[30,2,89.0],[31,0,63.0],[31,1,12.0],[31,2,38.0],[32,0,41.0],[32,1,99.0],[32,2,99.0],[33,0,27.0],[33,1,53.0],[33,2,8.0],[34,0,21.0],[34,1,24.0],[34,2,28.0],[35,0,12.0],[35,1,37.0],[35,2,97.0],[36,0,32.0],[36,1,87.0],[36,2,29.0],[37,0,68.0],[37,1,72.0],[37,2,27.0],[38,0,85.0],[38,1,34.0],[38,2,72.0],[39,0,66.0],[39,1,88.0],[39,2,17.0],[40,0,7.0],[40,1,30.0],[40,2,59.0],[41,0,86.0],[41,1,10.0],[41,2,7.0],[42,0,8.0],[42,1,3.0],[42,2,12.0],[43,0,42.0],[43,1,17.0],[43,2,40.0],[44,0,60.0],[44,1,19.0],[44,2,11.0],[45,0,87.0],[45,1,17.0],[45,2,65.0],[46,0,89.0],[46,1,18.0],[46,2,39.0],[47,0,11.0],[47,1,90.0],[47,2,5.0],[48,0,37.0],[48,1,64.0],[48,2,15.0],[49,0,11.0],[49,1,32.0],[49,2,73.0],[50,0,9.0],[50,1,11.0],[50,2,33.0],[51,0,781.0],[51,1,909.0],[51,2,70.0]],"rows": [{"id": "100", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Pasteurellales", "f__Pasteurellaceae", "g__Basfia", "s__1"]}},{"id": "101", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Pasteurellales", "f__Pasteurellaceae", "g__Basfia", "s__2"]}},{"id": "102", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Pasteurellales", "f__Pasteurellaceae", "g__Basfia", "s__3"]}},{"id": "103", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Pasteurellales", "f__Pasteurellaceae", "g__Basfia", "s__4"]}},{"id": "104", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Pasteurellales", "f__Pasteurellaceae", "g__Basfia", "s__5"]}},{"id": "105", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Pasteurellales", "f__Pasteurellaceae", "g__Basfia", "s__6"]}},{"id": "106", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Pasteurellales", "f__Pasteurellaceae", "g__Basfia", "s__7"]}},{"id": "107", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Pasteurellales", "f__Pasteurellaceae", "g__Basfia", "s__8"]}},{"id": "108", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Pasteurellales", "f__Pasteurellaceae", "g__Basfia", "s__9"]}},{"id": "109", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Pasteurellales", "f__Pasteurellaceae", "g__Basfia", "s__10"]}},{"id": "110", "metadata": {"taxonomy": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria (class)", "o__Actinomycetales", "f__Micromonosporaceae", "g__", "s__"]}},{"id": "111", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Bacillales", "f__Bacillaceae", "g__Bacillus", "s__"]}},{"id": "112", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Bacillales", "f__Staphylococcaceae", "g__Staphylococcus", "s__"]}},{"id": "113", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Alphaproteobacteria", "o__Sphingomonadales", "f__Sphingomonadaceae", "g__Novosphingobium", "s__"]}},{"id": "114", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Betaproteobacteria", "o__Burkholderiales", "f__", "g__Rubrivivax", "s__Rubrivivax gelatinosus"]}},{"id": "115", "metadata": {"taxonomy": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria (class)", "o__Actinomycetales", "f__Corynebacteriaceae", "g__Corynebacterium", "s__"]}},{"id": "116", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Alphaproteobacteria", "o__Rhizobiales", "f__Hyphomicrobiaceae", "g__Rhodoplanes", "s__"]}},{"id": "117", "metadata": {"taxonomy": ["k__Bacteria", "p__Bacteroidetes", "c__Flavobacteria", "o__Flavobacteriales", "f__Flavobacteriaceae", "g__", "s__"]}},{"id": "118", "metadata": {"taxonomy": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria (class)", "o__Coriobacteriales", "f__Coriobacteriaceae", "g__Adlercreutzia", "s__"]}},{"id": "119", "metadata": {"taxonomy": ["k__Bacteria", "p__Bacteroidetes", "c__Sphingobacteria", "o__Sphingobacteriales", "f__", "g__", "s__"]}},{"id": "120", "metadata": {"taxonomy": ["k__Bacteria", "p__Spirochaetes", "c__Spirochaetes (class)", "o__Spirochaetales", "f__Spirochaetaceae", "g__Treponema", "s__"]}},{"id": "121", "metadata": {"taxonomy": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria (class)", "o__Actinomycetales", "f__Microbacteriaceae", "g__Yonghaparkia", "s__"]}},{"id": "122", "metadata": {"taxonomy": ["k__Bacteria", "p__Bacteroidetes", "c__Bacteroidia", "o__Bacteroidales", "f__Rikenellaceae", "g__Alistipes", "s__"]}},{"id": "123", "metadata": {"taxonomy": ["k__Bacteria", "p__Bacteroidetes", "c__Sphingobacteria", "o__Sphingobacteriales", "f__Flexibacteraceae", "g__Hymenobacter", "s__"]}},{"id": "124", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__", "g__", "s__"]}},{"id": "125", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Lachnospiraceae", "g__Blautia", "s__"]}},{"id": "126", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Alphaproteobacteria", "o__Sphingomonadales", "f__Sphingomonadaceae", "g__Sphingomonas", "s__Sphingomonas wittichii"]}},{"id": "127", "metadata": {"taxonomy": ["k__Bacteria", "p__Tenericutes", "c__Mollicutes", "o__Entomoplasmatales", "f__Spiroplasmataceae", "g__Spiroplasma", "s__"]}},{"id": "128", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Clostridiales Family XI. Incertae Sedis", "g__Anaerococcus", "s__Anaerococcus hydrogenalis"]}},{"id": "129", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Alphaproteobacteria", "o__Rhizobiales", "f__Rhizobiaceae", "g__Rhizobium", "s__"]}},{"id": "130", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Alphaproteobacteria", "o__Rhodobacterales", "f__Rhodobacteraceae", "g__Rubellimicrobium", "s__"]}},{"id": "131", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Betaproteobacteria", "o__Neisseriales", "f__Neisseriaceae", "g__Microvirgula", "s__Microvirgula aerodenitrificans"]}},{"id": "132", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Alphaproteobacteria", "o__Rhizobiales", "f__Phyllobacteriaceae", "g__", "s__"]}},{"id": "133", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Alphaproteobacteria", "o__Rhizobiales", "f__Hyphomicrobiaceae", "g__Devosia", "s__"]}},{"id": "134", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Alphaproteobacteria", "o__Sphingomonadales", "f__", "g__", "s__"]}},{"id": "135", "metadata": {"taxonomy": ["k__Bacteria", "p__Cyanobacteria", "c__", "o__", "f__", "g__", "s__"]}},{"id": "136", "metadata": {"taxonomy": ["k__Bacteria", "p__Bacteroidetes", "c__Bacteroidia", "o__Bacteroidales", "f__", "g__", "s__"]}},{"id": "137", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Alphaproteobacteria", "o__Rhizobiales", "f__Hyphomicrobiaceae", "g__Rhodoplanes", "s__"]}},{"id": "138", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Alphaproteobacteria", "o__Sphingomonadales", "f__Erythrobacteraceae", "g__Erythromicrobium", "s__Erythromicrobium ramosum"]}},{"id": "139", "metadata": {"taxonomy": ["k__Bacteria", "p__Bacteroidetes", "c__Sphingobacteria", "o__Sphingobacteriales", "f__", "g__", "s__"]}},{"id": "140", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Deltaproteobacteria", "o__Myxococcales", "f__Myxococcaceae", "g__Corallococcus", "s__Corallococcus exiguus"]}},{"id": "141", "metadata": {"taxonomy": ["k__Bacteria", "p__Bacteroidetes", "c__Flavobacteria", "o__Flavobacteriales", "f__Flavobacteriaceae", "g__", "s__"]}},{"id": "142", "metadata": {"taxonomy": ["k__Bacteria", "p__Acidobacteria", "c__Acidobacteria (class)", "o__Acidobacteriales", "f__", "g__", "s__"]}},{"id": "143", "metadata": {"taxonomy": ["k__Bacteria", "p__Lentisphaerae", "c__Lentisphaerae (class)", "o__Victivallales", "f__Victivallaceae", "g__", "s__"]}},{"id": "144", "metadata": {"taxonomy": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria (class)", "o__Actinomycetales", "f__Micrococcaceae", "g__Arthrobacter", "s__"]}},{"id": "145", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Deltaproteobacteria", "o__NB1-j", "f__JTB38", "g__", "s__"]}},{"id": "146", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Lachnospiraceae", "g__Clostridium", "s__"]}},{"id": "147", "metadata": {"taxonomy": ["k__Bacteria", "p__Acidobacteria", "c__Acidobacteria (class)", "o__Acidobacteriales", "f__", "g__", "s__"]}},{"id": "148", "metadata": {"taxonomy": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria (class)", "o__Actinomycetales", "f__Mycobacteriaceae", "g__Mycobacterium", "s__Mycobacterium leprae"]}},{"id": "149", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Epsilonproteobacteria", "o__Campylobacterales", "f__Helicobacteraceae", "g__Helicobacter", "s__Helicobacter hepaticus"]}},{"id": "150", "metadata": {"taxonomy": ["k__Bacteria", "p__Bacteroidetes", "c__Bacteroidia", "o__Bacteroidales", "f__Porphyromonadaceae", "g__Parabacteroides", "s__Parabacteroides gordonii"]}},{"id": "151", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Betaproteobacteria", "o__", "f__", "g__", "s__"]}}],"columns": [{"id": "A", "metadata": null},{"id": "B", "metadata": null},{"id": "C", "metadata": null}]}'
expected_output_gamma_3=[u'104\t0.02\t0.144444444444',
                       '105\t0.0177777777778\t0.124444444444',
                       '103\t0.0177777777778\t0.106666666667',
                       '106\t0.0155555555556\t0.0888888888889',
                       '102\t0.0155555555556\t0.0733333333333',
                       '107\t0.0133333333333\t0.0577777777778',
                       '101\t0.0133333333333\t0.0444444444444',
                       '108\t0.0111111111111\t0.0311111111111']

if __name__ == '__main__':
    main()
