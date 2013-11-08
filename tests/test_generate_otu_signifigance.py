# make_phyla_plots_AGP_test.py

from unittest import TestCase, main
from StringIO import StringIO
from numpy import array
from generate_otu_signifigance_tables_AGP import (taxa_importer)

class GenerateOTUSignifiganceTablesTest(TestCase):
    
    def setUp(self):
        self.sample = array([0.0200, 0.7000, 0.0000, 0.00001, 0.1000, 0.0001, 
                       0.03000, 0.0000, 0.0001, 0.0050, 0.0020, 0.1427])
        self.population = []
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




      
        

if __name__ == '__main__':
    main()