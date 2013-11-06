# make_phyla_plots_AGP_test.py

from unittest import TestCase, main
from StringIO import StringIO
from numpy import array
from generate_otu_signifigance_tables_AGP import (taxa_importer)

class GenerateOTUSignifiganceTablesTest(TestCase):
    
    def setUP(self):
        pass     

    def tearDown(self):
        #
        pass

    def test_taxa_importer(self):

        # Sets up the known data 
        known_ids = array(['Ninth', 'Tenth', 'Eleventh'])
        known_taxa = array(['Dalek', 'Human', 'TimeLord'])
        known_matrix = array([ [2, 100, 150], 
                              [10,  25,  30], 
                              [ 0,   1,   0]])
        # Creates the test file
        test_file = StringIO('\tNinth\tTenth\tEleventh\nDalek\t2\t100\t150\n'\
            'Human\t10\t25\t30\nTimeLord\t0\t1\t0\n')

        # Runs the script
        [test_taxa, test_matrix, test_ids] = taxa_importer(test_file)   
        
        # Tests that values are correct
        for [id_, sample] in enumerate(known_ids):
            self.assertEqual(known_ids[id_], test_ids[id_])
        for [id_, taxa] in enumerate(known_taxa):
            self.assertEqual(known_taxa[id_], test_taxa[id_])        
        self.assertEqual(known_matrix.all(), test_matrix.all())
    
    def test_caculate_taxa_rank_1 (self):
        # Sets up the known data
        sample = [0.1000, 0.1000, ]
        test = [[0.0253, 0.0223, 0.0120, 0.0179, 0.0071, 0.0239, 0.0035, 0.0043, 0.0069, 0, 0.0294, 0.0183, 0.0075, 0.0287, 0],
                [0.5949, 0.5879, 0.6160, 0.6156, 0.5568, 0.5985, 0.5918, 0.6314, 0.6547, 0.6555, 0.5568, 0.6039, 0.5393, 0.5443, 0.5997]]
        taxa = [(u'k__Bacteria', u' p__Proteobacteria', u'c__Gammaproteobacteria'),
                (u'k__Bacteria', u' p__Bacteroidetes', u'c__Bacteroidia'),]

        known_high = [[u'class Bacteroidia', 0.2000, 0.0138, 0.1000/0.0138, 0.0254],
                      [u'', 0.7000],
                      [u'', 0     ],
                      [u'', 0.0001],
                      [u'', 0     ],
                      [u'', 0.1000],
                      [u'', 0.0001],
                      [u'', ]






      
        

if __name__ == '__main__':
    main()