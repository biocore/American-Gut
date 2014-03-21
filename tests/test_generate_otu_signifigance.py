#!/usr/bin/env python 
# make_phyla_plots_AGP_test.py

from unittest import TestCase, main
from numpy import array
from americangut.generate_otu_signifigance_tables import (calculate_abundance,
                                                          calculate_tax_rank_1,
                                                          convert_taxa,
                                                          clean_otu_string,
                                                          render_latex_header,
                                                          render_raw_header,
                                                          convert_taxa_to_list,
                                                          generate_latex_macro,
                                                          convert_taxa_to_table)

class GenerateOTUSignifiganceTablesTest(TestCase):
    
    def setUp(self):
        self.sample = array([0.0200, 0.7000, 0.0050, 0.0001, 0.0020, 
                             0.1000, 0.0002, 0.0300, 0.1427, 0.0000])

        self.taxa = ['k__Bacteria', 
                     'k__Bacteria; p__[Proteobacteria]', 
                     'k__Bacteria; p__Proteobacteria; '\
                     'c__Gammaproteobacteria', 
                     'k__Bacteria; p__Proteobacteria; '\
                     'c__Gammaproteobacteria; o__Enterobacteriales', 
                     'k__Bacteria; p__Proteobacteria; '\
                     'c__Gammaproteobacteria; o__Enterobacteriales; '\
                     'f__Enterbacteriaceae', 
                     'k__Bacteria; p__Proteobacteria; '\
                     'c__Gammaproteobacteria; o__Enterobacteriales; '\
                     'f__Enterbacteriaceae; g__Escherichia', 
                     'k__Bacteria; p__Proteobacteria; '\
                     'c__Gammaproteobacteria; o__Enterobacteriales; '\
                     'f__Enterbacteriaceae; g__Escherichia; s__coli',
                     'k__Archaea; p__Crenarchaeota; c__Thaumarchaeota; '\
                     'o__Cenarchaeales; f__Cenarchaeaceae; g__Nitrosopumilus',
                     'k__Bacteria; p__Actinobacteria; c__Coriobacteriia; '\
                     'o__Coriobacteriales; f__Coriobacteriaceae; g__',
                     'k__Bacteria; p__Actinobacteria; c__Actinobacteria; '\
                     'o__Actinomycetales; f__Dietziaceae; g__']

        self.pop = array([[ 0.0295,  0.2715,  1.0457,  0.5374,  2.0631,  
                            0.1329,  1.2567, -0.7502, -0.8432, -1.6957, 
                           -0.0425,  0.5318,  0.2917,  1.0237,  0.9305, 
                           -0.6820,  2.5324,  1.8886,  1.4483,  0.1470,  
                            2.3343,  1.4135,  1.9938, -0.2980,  2.0538],
                          [ 0.1841, -0.1294, -0.3062, -0.1337, -0.6649, 
                           -0.0647, -0.2759,  0.8383,  1.2143, -0.0176,  
                            0.2173, -0.2600,  0.9223, -0.5197,  1.1234,  
                            0.6878, -0.8453, -0.0296,  0.0553,  0.6151, 
                            1.0436,  0.2590, -0.0437, -1.2562, -0.0673],
                          [ 0.0162, -0.0024,  0.0109,  0.0088, -0.0104,  
                            0.0026,  0.0134,  0.0098,  0.0114,  0.0055, 
                           -0.0043,  0.0097,  0.0069,  0.0022, -0.0104,  
                            0.0052, -0.0036,  0.0212,  0.0096,  0.0231,  
                            0.0061,  0.0022, -0.0009, -0.0008,  0.0010],
                          [ 0.0380, -0.0131,  0.0244,  0.0078,  0.0399, 
                           -0.0173,  0.0162, -0.0105, -0.0006,  0.0094,  
                            0.0001,  0.0173,  0.0103,  0.0165, -0.0191, 
                            0.0029, -0.0005,  0.0093,  0.0087,  0.0262,  
                            0.0124,  0.0343,  0.0331, -0.0183, -0.0132],
                          [-0.0001,  0.0001, -0.0001,  0.0001,  0.0002, 
                            0.0001,  0.0000,  0.0004,  0.0000,  0.0002, 
                            0.0000,  0.0003, -0.0003,  0.0000,  0.0001, 
                           -0.0003,  0.0007,  0.0001,  0.0005,  0.0002,  
                            0.0001, -0.0002,  0.0000,  0.0001, -0.0003],
                          [-0.0069,  0.1084,  0.1284,  0.0377, -0.0716, 
                           -0.0835,  0.0551,  0.0654,  0.2287,  0.0954, 
                            0.0538,  0.0607,  0.0247, -0.0129,  0.1440, 
                            0.0320,  0.0756,  0.0834, -0.0377, -0.0518, 
                           -0.0047, -0.1675,  0.1091,  0.1293,  0.1834],
                          [-0.0004,  0.0007, -0.0004, -0.0001,  0.0001, 
                           -0.0005,  0.0004,  0.0002,  0.0002, -0.0004,  
                            0.0004,  0.0000,  0.0008,  0.0008,  0.0001,  
                            0.0003,  0.0007,  0.0007,  0.0002,  0.0011,  
                            0.0005,  0.0007,  0.0000, -0.0004,  0.0006],
                          [ 0.4108,  0.4489, -0.1190, -0.0411,  0.5914, 
                            0.0331,  0.7549, -0.1050,  0.9035, -0.4558, 
                           -0.7154,  0.2620,  0.8578,  0.5640,  0.5965,  
                            0.1027,  0.2566,  0.2922, -0.2714,  0.7830,  
                            0.0444,  0.3566,  0.0561,  0.8699,  0.7491],
                          [ 0.6829,  0.0597, -0.4357,  0.4488,  0.9149, 
                            0.2044, -0.5362,  1.4112, -0.4213,  1.1725, 
                            0.7079,  0.4245,  0.3442,  0.5083,  0.0462,  
                            0.5746,  0.5843,  0.7929, -0.3203,  0.4132,  
                            0.2063,  1.3558, -0.7080,  0.6477, -0.0417],
                          [-0.0332,  0.1087,  0.0870,  0.0829,  0.1622,
                            0.0708, -0.0424,  0.0287, -0.0314,  0.0499,  
                            0.0245,  0.0633,  0.0772,  0.1273, -0.0030,  
                            0.0703, -0.0372,  0.0931,  0.0883, -0.0154,  
                            0.0302,  0.1436, -0.0564,  0.0490, -0.0296]])

        self.test_list = ['k__Bacteria', 
                          'k__Bacteria; p__[Proteobacteria]', 
                          'k__Bacteria; p__Proteobacteria; '\
                          'c__Gammaproteobacteria', 
                          'k__Bacteria; p__Proteobacteria; '\
                          'c__Gammaproteobacteria; o__Enterobacteriales', 
                          'k__Bacteria; p__Proteobacteria; '\
                          'c__Gammaproteobacteria; o__Enterobacteriales; '\
                          'f__Enterbacteriaceae', 
                          'k__Bacteria; p__Proteobacteria; '\
                          'c__Gammaproteobacteria; o__Enterobacteriales; '\
                          'f__Enterbacteriaceae; g__Escherichia', 
                          'k__Bacteria; p__Proteobacteria; '\
                          'c__Gammaproteobacteria; o__Enterobacteriales; '\
                          'f__Enterbacteriaceae; g__Escherichia; s__coli']
        
        self.test_table = [['k__Bacteria; p__Proteobacteria; '\
                            'c__Gammaproteobacteria; o__Enterobacteriales; '\
                            'f__Enterbacteriaceae; g__Escherichia', '10.00%',\
                            '1.00%', '10'],
                           ['k__Bacteria; p__Proteobacteria; '\
                            'c__Gammaproteobacteria; o__Enterobacteriales; '\
                            'f__Enterbacteriaceae; g__Escherichia; s__coli',\
                            '0.10%', '1.00%', '0.10']] 

        self.header = ['Taxonomy', 'Doctor', 'Humans', 'Fold']  

    def test_calculate_abundance(self):
        """Checks that abundance is calculated sanely"""
        # Checks errors are thrown when the sample and taxa are different lengths
        with self.assertRaises(ValueError):
            calculate_abundance(self.sample[0:2], self.taxa)

        # Sets up known value
        known_abundance_95 = [['k__Bacteria; p__[Proteobacteria]', 0.7],
                              ['k__Bacteria; p__Actinobacteria; '\
                               'c__Coriobacteriia; o__Coriobacteriales; '\
                               'f__Coriobacteriaceae; g__', 0.1427],
                              ['k__Bacteria; p__Proteobacteria; '\
                               'c__Gammaproteobacteria; o__Enterobacteriales; '\
                               'f__Enterbacteriaceae; g__Escherichia', 0.1],
                               ['k__Archaea; p__Crenarchaeota; '\
                               'c__Thaumarchaeota; o__Cenarchaeales; '\
                               'f__Cenarchaeaceae; g__Nitrosopumilus', 0.03]]

        known_abundance_99 = [['k__Bacteria; p__[Proteobacteria]', 0.7],
                              ['k__Bacteria; p__Actinobacteria; '\
                               'c__Coriobacteriia; o__Coriobacteriales; '\
                               'f__Coriobacteriaceae; g__', 0.1427],
                              ['k__Bacteria; p__Proteobacteria; '\
                               'c__Gammaproteobacteria; o__Enterobacteriales; '\
                               'f__Enterbacteriaceae; g__Escherichia', 0.1],
                              ['k__Archaea; p__Crenarchaeota; '\
                               'c__Thaumarchaeota; o__Cenarchaeales; '\
                               'f__Cenarchaeaceae; g__Nitrosopumilus', 0.03],
                               ['k__Bacteria', 0.02]]

        known_abundance_1 = [['k__Bacteria; p__[Proteobacteria]', 0.7], 
                             ['k__Bacteria; p__Actinobacteria; '\
                              'c__Coriobacteriia; o__Coriobacteriales; '\
                              'f__Coriobacteriaceae; g__', 0.1427], 
                             ['k__Bacteria; p__Proteobacteria; '\
                              'c__Gammaproteobacteria; o__Enterobacteriales; '\
                              'f__Enterbacteriaceae; g__Escherichia', 0.1],
                             ['k__Archaea; p__Crenarchaeota; '\
                              'c__Thaumarchaeota; o__Cenarchaeales; '\
                              'f__Cenarchaeaceae; g__Nitrosopumilus', 0.03],
                             ['k__Bacteria', 0.02],
                             ['k__Bacteria; p__Proteobacteria; '\
                              'c__Gammaproteobacteria', 0.005],
                             ['k__Bacteria; p__Proteobacteria; '\
                              'c__Gammaproteobacteria; o__Enterobacteriales; '\
                              'f__Enterbacteriaceae', 0.002],
                             ['k__Bacteria; p__Proteobacteria; '\
                              'c__Gammaproteobacteria; o__Enterobacteriales; '\
                              'f__Enterbacteriaceae; g__Escherichia; s__coli', 
                               0.0002],
                             ['k__Bacteria; p__Proteobacteria; '\
                              'c__Gammaproteobacteria; o__Enterobacteriales', 
                               0.0001],
                             ['k__Bacteria; p__Actinobacteria; '\
                              'c__Actinobacteria; o__Actinomycetales; '\
                              'f__Dietziaceae; g__', 0]]

        # Generates the value for testing
        test_abundance_def = calculate_abundance(sample = self.sample, 
                                                 taxa = self.taxa)
        self.assertEqual(test_abundance_def, known_abundance_95)


        test_abundance_99 = calculate_abundance(sample = self.sample, 
                                                taxa = self.taxa, 
                                                sum_min = 0.99)
        self.assertEqual(test_abundance_99, known_abundance_99)

        test_abundance_1 = calculate_abundance(sample = self.sample, 
                                               taxa = self.taxa, 
                                               sum_min = 1.000)
        self.assertEqual(test_abundance_1, known_abundance_1)

    def test_calculate_tax_rank_1(self):
        # Sets up known values
        known_high_10 = [['k__Bacteria; p__Proteobacteria; '\
                          'c__Gammaproteobacteria; o__Enterobacteriales; '\
                          'f__Enterbacteriaceae', 0.002, 7.6e-05, 26.0, 
                          1.4507298345686689e-22], 
                         ['k__Bacteria; p__[Proteobacteria]', 0.7, 0.101852, 
                          7.0, 0.00073396297302392652], 
                         ['k__Bacteria; p__Proteobacteria; '\
                          'c__Gammaproteobacteria; o__Enterobacteriales; '\
                          'f__Enterbacteriaceae; g__Escherichia', 0.1, 0.04714, 
                          2.0, 0.068040408640427208]]
        known_high_05 = [['k__Bacteria; p__Proteobacteria; '\
                          'c__Gammaproteobacteria; o__Enterobacteriales; '\
                          'f__Enterbacteriaceae', 0.002, 7.6e-05, 26.0, 
                          1.4507298345686689e-22], 
                         ['k__Bacteria; p__[Proteobacteria]', 0.7, 0.101852, 
                          7.0, 0.00073396297302392652]]
        known_high_01 = [['k__Bacteria; p__Proteobacteria; '\
                          'c__Gammaproteobacteria; o__Enterobacteriales; '\
                          'f__Enterbacteriaceae', 0.002, 7.6e-05, 
                           26.0, 1.4507298345686689e-22], 
                         ['k__Bacteria; p__[Proteobacteria]', 0.7, 0.101852, 
                          7.0, 0.00073396297302392652]]
        known_low_10  = [['k__Bacteria; p__Actinobacteria; c__Actinobacteria;'\
                          ' o__Actinomycetales; f__Dietziaceae; g__', 
                          0.0, 0.044336, 0.0, 0.016163391728614671], 
                         ['k__Bacteria', 
                         0.02, 0.704584, 0.0, 0.051170083967289066], 
                         ['k__Archaea; p__Crenarchaeota; c__Thaumarchaeota;'\
                         ' o__Cenarchaeales; f__Cenarchaeaceae; '\
                         'g__Nitrosopumilus', 
                         0.03, 0.289032, 0.0, 0.065525264907055]]
        known_low_05  = [['k__Bacteria; p__Actinobacteria; c__Actinobacteria;'\
                          ' o__Actinomycetales; f__Dietziaceae; g__', 
                          0.0, 0.044336, 0.0, 0.016163391728614671]]
        known_low_01  = []
        # Sets up test values
        (test_high_10, test_low_10) = calculate_tax_rank_1(sample = self.sample,
                                                        population = self.pop,
                                                        taxa = self.taxa,
                                                        critical_value = 0.1)

        (test_high_05, test_low_05) = calculate_tax_rank_1(sample = self.sample,
                                                         population = self.pop,
                                                         taxa = self.taxa)
        
        (test_high_01, test_low_01) = calculate_tax_rank_1(sample = self.sample,
                                                        population = self.pop,
                                                        taxa = self.taxa,
                                                        critical_value = 0.01)

        self.assertEqual(known_high_10, test_high_10)
        self.assertEqual(known_low_10, test_low_10)
        self.assertEqual(known_high_05, test_high_05)
        self.assertEqual(known_low_05, test_low_05)
        self.assertEqual(known_high_01, test_high_01)
        self.assertEqual(known_low_01, test_low_01)

    def test_convert_taxa(self):
        """Checks that convert_taxa runs sanely"""
        # Sets up test values
        test_value = [['Rose', 9, 10.101, 0.7401, 97.643, 0.6802, 0.2252], 
                      ['Martha', 10, 10.201, 0.3823, 20.027, 0.0023, 0.2302]]
        test_keys = ['%i', '%1.2f', '%1.2f', '%1.1f\\%%', '%1.1f\\%%', 'SKIP']
        test_hund = [False, False, True, False, True, False]              

        # Checks that errors are called correctly
        with self.assertRaises(TypeError):
            convert_taxa('test_value', test_keys, test_hund)
        
        with self.assertRaises(TypeError):
            convert_taxa(test_value.append('Cats'), test_keys, test_hund)
        test_value.pop(test_value.index('Cats'))
        test_value.append(['Cats'])

        with self.assertRaises(ValueError):
            convert_taxa(test_value, test_keys, test_hund)
        test_value.pop(test_value.index(['Cats']))

        with self.assertRaises(TypeError):
            convert_taxa(test_value, set(test_keys), test_hund)
        with self.assertRaises(ValueError):
            convert_taxa(test_value, test_keys[0:1], test_hund)

        with self.assertRaises(TypeError):
            convert_taxa(test_value, set(test_keys), test_hund)
        with self.assertRaises(ValueError):
            convert_taxa(test_value, test_keys[0:1], test_hund)

        with self.assertRaises(TypeError):
            convert_taxa(test_value, test_keys, set(test_hund))
        with self.assertRaises(ValueError):
            convert_taxa(test_value, test_keys, test_hund[0:1])

        # Tests with list values
        test  = convert_taxa(test_value, test_keys, test_hund)
        known = [['Rose', '9', '10.10', '74.01', '97.6\%', '68.0\%'],
                 ['Martha', '10', '10.20', '38.23', '20.0\%', '0.2\%']]
        self.assertEqual(test, known)

        # Tests with string key
        key = '%i'
        test_value = [['Rose', 9, 10.101, 0.7401, 97.643, 0.6802, 0.2252], 
                      ['Martha', 10, 10.201, 0.3823, 20.027, 0.0023, 0.2302]]
        test  = convert_taxa(test_value, key, test_hund)
        known = [['Rose', '9', '10', '74', '97', '68', '0'],
                 ['Martha', '10', '10', '38', '20', '0', '0']]
        self.assertEqual(known, test)

         # Tests with hundred key
        hund = False
        test_value = [['Rose', 9, 10.101, 0.7401, 97.643, 0.6802, 0.2252], 
                      ['Martha', 10, 10.201, 0.3823, 20.027, 0.0023, 0.2302]]
        test  = convert_taxa(test_value, test_keys, hund)
        known = [['Rose', '9', '10.10', '0.74', '97.6\%', '0.7\\%'],
                 ['Martha', '10', '10.20', '0.38', '20.0\%', '0.0\\%']]
        self.assertEqual(known, test)

    def test_convert_taxa_to_list(self):
        format_list = ['BOLD', 'REG', 'REG', 'REG', 'REG', 'REG', 'COLOR']

        known_latex_format = '\\begin{itemize}\n\\item \\textbf{Unclassified '\
        'Kingdom Bacteria}\n\\item Unclassified Phylum Proteobacteria\n\\item '\
        'Unclassified Class Gammaproteobacteria\n\\item Unclassified Order '\
        'Enterobacteriales\n\\item Unclassified Family Enterbacteriaceae\n'\
        '\\item Genus \\textit{Escherichia}\n\\item \\textcolor{red}{\\textit'\
        '{Escherichia coli}}\n\\end{itemize}'

        known_latex_comma = '\\textbf{Unclassified Kingdom Bacteria}, '\
        'Unclassified Phylum Proteobacteria, Unclassified Class '\
        'Gammaproteobacteria, Unclassified Order Enterobacteriales, '\
        'Unclassified Family Enterbacteriaceae, Genus \\textit{Escherichia}, '\
        '\\textcolor{blue}{\\textit{Escherichia coli}}'

        known_raw_format = '\n o   *Unclassified Kingdom Bacteria*\n o   '\
        'Unclassified Phylum Proteobacteria\n o   Unclassified Class '\
        'Gammaproteobacteria\n o   Unclassified Order Enterobacteriales\n o   '\
        'Unclassified Family Enterbacteriaceae\n o   Genus Escherichia\n o   '\
        '*Escherichia coli*'

        known_raw_comma = '*Unclassified Kingdom Bacteria*, Unclassified '\
        'Phylum Proteobacteria, Unclassified Class Gammaproteobacteria, '\
        'Unclassified Order Enterobacteriales, Unclassified Family '\
        'Enterbacteriaceae, Genus Escherichia, *Escherichia coli*'

        test_list_latex_format = convert_taxa_to_list(self.test_list, 
                                                      render_mode = 'LATEX',
                                                      tax_format = format_list, 
                                                      comma = False, 
                                                      color = 'red')
        test_list_latex_comma = convert_taxa_to_list(self.test_list, 
                                                      render_mode = 'LATEX',
                                                      tax_format = format_list,
                                                      comma = True, 
                                                      color = 'blue')
        test_list_raw_format = convert_taxa_to_list(self.test_list,
                                                    render_mode = 'RAW',
                                                    tax_format = format_list,
                                                    comma = False)
        test_list_raw_comma = convert_taxa_to_list(self.test_list,
                                                    render_mode = 'RAW',
                                                    tax_format = format_list,
                                                    comma = True)

        self.assertEqual(test_list_latex_format, known_latex_format)
        self.assertEqual(test_list_latex_comma, known_latex_comma)
        self.assertEqual(test_list_raw_format, known_raw_format)
        self.assertEqual(test_list_raw_comma, known_raw_comma)

    def test_clean_otu_string(self):
        known_raw = ['*Kingdom Bacteria*',
                     'Phylum Proteobacteria',
                     'Class Gammaproteobacteria',
                     'Order Enterobacteriales',
                     'Family Enterbacteriaceae',
                     '*Genus Escherichia*',
                     '*Escherichia coli*']
        known_latex = ['\\textbf{Unclassified Kingdom Bacteria}',
                       'Unclassified Phylum Proteobacteria',
                       'Unclassified Class Gammaproteobacteria',
                       'Unclassified Order Enterobacteriales',
                       'Unclassified Family Enterbacteriaceae',
                       '\\textcolor{red}{Genus \\textit{Escherichia}}',
                       '\\textcolor{blue}{\\textit{Escherichia coli}}']

        for idx, taxon in enumerate(self.test_list):
            # Sets the format
            if idx == 0:
                format = 'BOLD'
                color = 'red'
            elif idx == 5:
                format = 'COLOR'                
                color = 'red'
            elif idx == 6:
                format = 'COLOR'
                color = 'blue'
            else:
                format = 'NONE'
                color = 'red'

            # Generates the cleaned string
            test_string_raw = clean_otu_string(taxon, render_mode = 'RAW', \
                format = format, color = color)
            test_string_latex = clean_otu_string(taxon, render_mode = 'LATEX',\
                format = format, color = color, unclassified = True)
            self.assertEqual(test_string_raw, known_raw[idx])
            self.assertEqual(test_string_latex, known_latex[idx])

    def test_render_latex_header(self):
        test_alignment = ['l', 'c', 'c', 'r']

        # Checks that appropriate errors are raised
        with self.assertRaises(ValueError):
            render_latex_header(['Doctor', 'Companion', 'Writer', 'Queen', \
                'Years'], alignment = ['l', 'c'])

        with self.assertRaises(TypeError):
            render_latex_header(self.header, alignment = 2)

        # Checks that the code is loading properly
        test_center_number = render_latex_header(self.header, 
                                                 numbering = True)

        test_align_no_nums = render_latex_header(self.header, \
                                                 alignment = test_alignment,\
                                                 numbering = False)

        known_header_center_number = '\\begin{tabular}{r c c c c }\n\\hline\n '\
        ' & Taxonomy & Doctor & Humans & Fold \\\\\n\\hline\n'

        known_header_align_no_num = '\\begin{tabular}{ l c c r }\n\\hline\n'\
        'Taxonomy & Doctor & Humans & Fold \\\\\n\\hline\n'

        self.assertEqual(test_center_number, known_header_center_number)

        self.assertEqual(test_align_no_nums, known_header_align_no_num)

    def test_render_raw_header(self):
        # Sets constants for table formatting
        test_header = ['Doctor', 'Companion(s)', 'Writer']
        category_length = 23        
        HEADER_BAR = "--------------------------------------------------------"\
        "-------------------"
        SPACER = '                                    '
        TAX_SPACE = 27
        
        # Sets up known values
        known_header_numbering = '--------------------------------------------'\
        '------------------------------------\n     Doctor                    '\
        ' Companion(s)           Writer                 \n--------------------'\
        '-------------------------------------------------------\n'
        known_header_no_number = '--------------------------------------------'\
        '-------------------------------\nDoctor                     Companion'\
        '(s)           Writer                 \n------------------------------'\
        '---------------------------------------------\n'

        # Checks the function
        test_header_numbering = render_raw_header(test_header,
                                                  numbering = True, 
                                                  header_bar = HEADER_BAR, 
                                                  spacer = SPACER,
                                                  tax_len = TAX_SPACE, 
                                                  cat_len = category_length)
        
        test_header_no_number = render_raw_header(test_header,
                                                  numbering = False, 
                                                  header_bar = HEADER_BAR, 
                                                  spacer = SPACER, 
                                                  tax_len = TAX_SPACE, 
                                                  cat_len = category_length)


        self.assertEqual(known_header_numbering, test_header_numbering)
        self.assertEqual(known_header_no_number, test_header_no_number)

    def test_render_latex_macro(self):
        # Sets up variable for testing        
        known_macro = '\\def\\TaxonomyA{Genus \\textit{Escherichia}}\n'\
                      '\\def\\DoctorA{10.00%}\n'\
                      '\\def\\HumansA{1.00%}\n'\
                      '\\def\\FoldA{10}\n'\
                      '\\def\\TaxonomyB{\\textit{Escherichia coli}}\n'\
                      '\\def\\DoctorB{0.10%}\n'\
                      '\\def\\HumansB{1.00%}\n'\
                      '\\def\\FoldB{0.10}'

        # Calculates test variable
        test_macro = generate_latex_macro(self.test_table, self.header)

        # Checks the outputs match
        self.assertEqual(known_macro, test_macro)

    def test_convert_taxa_table(self):
        # Sets up inputs
        header_code = '<THIS IS HEADER CODE FOR LATEX>'
        alignment_vary = ['l', 'c', 'c', 'r']

        # Sets up known values
        default_known = '-----------------------------------------------------'\
                        '---------------------------\n'\
                        '     Taxonomy                   Doctor         Humans'\
                        '         Fold           \n'\
                        '-----------------------------------------------------'\
                        '----------------------\n\n'\
                        '( 1) Genus Escherichia          10.00%         1.00% '\
                        '         10             \n'\
                        '( 2) Escherichia coli           0.10%          1.00% '\
                        '         0.10           \n'\
                        '-----------------------------------------------------'\
                        '---------------------------'

        latex_default_known = '\\begin{tabular}{r c c c c }\n\\hline\n'\
                              '  & Taxonomy & Doctor & Humans & Fold \\\\\n'\
                              '\\hline\n'\
                              '\\\\\n'\
                              '( 1)  & Genus \\textit{Escherichia} & 10.00% &'\
                              ' 1.00% & 10\\\\\n'\
                              '( 2)  & \\textit{Escherichia coli} & 0.10% & '\
                              '1.00% & 0.10\\\\\n\\hline\n\\end{tabular}'

        raw_no_numbers_known = '----------------------------------------------'\
                               '-----------------------------\n'\
                               'Taxonomy                   Doctor         '\
                               'Humans         Fold           \n'\
                               '----------------------------------------------'\
                               '-----------------------------\n\n'\
                               'Genus Escherichia          10.00%         '\
                               '1.00%          10             \n'\
                               'Escherichia coli           0.10%          '\
                               '1.00%          0.10           \n'\
                               '----------------------------------------------'\
                               '-----------------------------'

        latex_no_number_known = '\\begin{tabular}{ c c c c }\n'\
                                '\\hline\n'\
                                'Taxonomy & Doctor & Humans & Fold \\\\\n'\
                                '\\hline\n'\
                                '\\\\\n'\
                                'Genus \\textit{Escherichia} & 10.00% & 1.00% '\
                                '& 10\\\\\n'\
                                '\\textit{Escherichia coli} & 0.10% & 1.00% & '\
                                '0.10\\\\\n'\
                                '\\hline\n'\
                                '\\end{tabular}'

        raw_alignment_known = '-----------------------------------------------'\
                              '---------------------------------\n'\
                              '     Taxonomy                   Doctor         '\
                              'Humans         Fold           \n'\
                              '-----------------------------------------------'\
                              '----------------------------\n'\
                              '\n'\
                              '( 1) Genus Escherichia          10.00%         '\
                              '1.00%          10             \n'\
                              '( 2) Escherichia coli           0.10%          '\
                              '1.00%          0.10           \n'\
                              '-----------------------------------------------'\
                              '---------------------------------'

        latex_alignment_known = '\\begin{tabular}{r l c c r }\n'\
                                '\\hline\n'\
                                '  & Taxonomy & Doctor & Humans & Fold \\\\\n'\
                                '\\hline\n'\
                                '\\\\\n'\
                                '( 1)  & Genus \\textit{Escherichia} & 10.00% '\
                                '& 1.00% & 10\\\\\n'\
                                '( 2)  & \\textit{Escherichia coli} & 0.10% & '\
                                '1.00% & 0.10\\\\\n'\
                                '\\hline\n'\
                                '\\end{tabular}'

        raw_code_known = '----------------------------------------------------'\
                         '----------------------------\n'\
                         '     Taxonomy                   Doctor         '\
                         'Humans         Fold           \n'\
                         '----------------------------------------------------'\
                         '-----------------------\n'\
                         '\n'\
                         '( 1) Genus Escherichia          10.00%         1.00%'\
                         '          10             \n'\
                         '( 2) Escherichia coli           0.10%          1.00%'\
                         '          0.10           \n'\
                         '----------------------------------------------------'\
                         '----------------------------'

        latex_code_known = '<THIS IS HEADER CODE FOR LATEX>\\\\\n'\
                           '( 1)  & Genus \\textit{Escherichia} & 10.00% & '\
                           '1.00% & 10\\\\\n'\
                           '( 2)  & \\textit{Escherichia coli} & 0.10% & 1.00%'\
                           ' & 0.10\\\\\n'\
                           '\\hline\n'\
                           '\\end{tabular}'

        # Generates tables and compares to the knowns
        default_test = convert_taxa_to_table(corr_taxa = self.test_table, 
                                             header = self.header)
        self.assertEqual(default_test, default_known)

        latex_default_test = convert_taxa_to_table(corr_taxa = self.test_table,
                                                   header = self.header,
                                                   render_mode = 'LATEX')
        self.assertEqual(latex_default_test, latex_default_known)

        raw_no_nums_test = convert_taxa_to_table(corr_taxa = self.test_table, 
                                                 header = self.header,
                                                 render_mode = 'RAW',
                                                 numbering = False)
        self.assertEqual(raw_no_nums_test, raw_no_numbers_known)

        latex_no_nums_test = convert_taxa_to_table(corr_taxa = self.test_table, 
                                                   header = self.header,
                                                   render_mode = 'LATEX',
                                                   numbering = False)
        self.assertEqual(latex_no_nums_test, latex_no_number_known)
        
        raw_align_test = convert_taxa_to_table(corr_taxa = self.test_table, 
                                               header = self.header,
                                               render_mode = 'RAW',
                                               numbering = True,
                                               alignment = alignment_vary)
        self.assertEqual(raw_align_test, raw_alignment_known)

        latex_align_test = convert_taxa_to_table(corr_taxa = self.test_table, 
                                                 header = self.header,
                                                 render_mode = 'LATEX',
                                                 numbering = True,
                                                 alignment = alignment_vary)
        self.assertEqual(latex_align_test, latex_alignment_known)

        raw_header_test = convert_taxa_to_table(corr_taxa = self.test_table, 
                                                header = self.header,
                                                render_mode = 'RAW',
                                                numbering = True,
                                                alignment = alignment_vary,
                                                header_code = header_code)
        self.assertEqual(raw_header_test, raw_code_known)

        latex_header_test = convert_taxa_to_table(corr_taxa = self.test_table, 
                                                  header = self.header,
                                                  render_mode = 'LATEX',
                                                  numbering = True,
                                                  alignment = alignment_vary,
                                                  header_code = header_code)
        self.assertEqual(latex_header_test, latex_code_known)
        

if __name__ == '__main__':
    main()