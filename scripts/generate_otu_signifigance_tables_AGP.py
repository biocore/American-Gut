#!/usr/bin/env python

from argparse import ArgumentParser
from os import mkdir
from os.path import isfile, exists, join as pjoin
from americangut.generate_otu_signifigance_tables import (taxa_importer,
									          calculate_tax_rank_1,
									          convert_taxa,
									          convert_taxa_to_list,
									          generate_latex_macro)

def main(taxa, table, sample_ids, output_dir, samples_to_analyze = None):
    """Creates LaTeX formatted significant OTU lists

    INPUTS:
        taxa -- a numpy vector with greengenes taxonomy strings

        tax_table -- a numpy array with the relative frequencies of taxonomies
            (rows) for each give sample (column)

        samples_to_analyze -- a numpy vector of sample ids associated with the 
            tax_table values

        output_dir -- a directory where the final files should be saved.

        samples_to_analyze -- a list of samples_to_analyze which should be used 
                    to generate data. If this is left empty, all the samples in 
                    the table will be used.

    OUTPUTS:
        Generates text files containing LaTex encoded strings which creates a 
        LaTeX macro dictionary with the information for creating a table of most
        abundant taxa, most enriched taxa, and rare and unique taxa. Rare 
        defined as present in less than 10% of the total population. The unique 
        taxa are bolded in the lists. 
    """
    # Sets table constants
    RENDERING = "LATEX"
    FORMAT_SIGNIFIGANCE = ["VAL_100", "VAL_100", "VAL_INT", "SKIP"]
    FORMAT_ABUNDANCE = ["VAL_100"]    
    MACRO_CATS_SIGNIFICANCE = ['enrichTaxon','enrichSampl', 'enrichPopul', 'enrichFoldd']
    MACRO_CATS_ABUNDANCE = ['abundTaxon', 'abundSampl']

    FILE_PRECURSER = 'macros_'
    FILE_EXTENSION = '.tex'

    # Number of taxa shown is an indexing value, it is one less than what is 
    # actually shown.
    NUMBER_OF_TAXA_SHOWN = 5

    # Sets up samples for which tables are being generated    
    if samples_to_analyze == None:
        samples_to_test = sample_ids
    else:
        samples_to_test = samples_to_analyze

    for idx, sample_id in enumerate(sample_ids):
        if sample_id in samples_to_test:
            sample = table[:,idx]

            population = delete(table, idx, 1)

            # Calculates tax rank tables
            (unique, rare, low, high, absent, abundance) = \
                calculate_tax_rank_1(sample, population, taxa)

            # Generates formatted enriched table
            formatted_high = convert_taxa(high[0:NUMBER_OF_TAXA_SHOWN],
                                          render_mode = RENDERING, 
                                          formatting_keys = FORMAT_SIGNIFIGANCE)

            high_formatted = generate_latex_macro(formatted_high, \
                categories = MACRO_CATS_SIGNIFICANCE)

            # Generates formatted abundance table
            formatted_abundance = convert_taxa(abundance[0:NUMBER_OF_TAXA_SHOWN],
                                            render_mode = RENDERING,
                                            formatting_keys = FORMAT_ABUNDANCE)
            abundance_formatted = generate_latex_macro(formatted_abundance, \
                categories = MACRO_CATS_ABUNDANCE)
        
            # Generates formatted list
            rare_format = []
            rare_combined = []
            for taxon in unique:
                rare_combined.append(taxon)
                rare_format.append('COLOR')
            for taxon in rare:
                rare_combined.append(taxon)
                rare_format.append('REG')

            number_rare_tax = len(rare_combined)

            if number_rare_tax > NUMBER_OF_TAXA_SHOWN + 1 and len(unique) == 0:
                rare_formatted = ["Your sample contained %i rare "\
                "taxa, including the following: " % number_rare_tax]
                rare_formatted.append(convert_taxa_to_list(\
                    rare_combined[:NUMBER_OF_TAXA_SHOWN ], 
                    tax_format = rare_format,
                    render_mode = RENDERING, 
                    comma = True))
                rare_formatted = ''.join(rare_formatted)                

            elif number_rare_tax > NUMBER_OF_TAXA_SHOWN + 1:
                rare_formatted = ["This sample contained %i rare and " \
                     "\\textcolor{red}{%i unique} taxa, including "\
                     "the following: " % (len(rare), len(unique))]
                rare_formatted.append(convert_taxa_to_list(\
                    rare_combined[:NUMBER_OF_TAXA_SHOWN ], 
                    tax_format = rare_format,
                    render_mode = RENDERING, 
                    comma = True))
                rare_formatted = ''.join(rare_formatted)

            elif number_rare_tax > 0 and len(unique) == 0:
                rare_formatted = ['This sample included the following rare taxa: ']
                rare_formatted.append(convert_taxa_to_list(rare_combined, 
                                                    tax_format = rare_format,
                                                    render_mode = RENDERING, 
                                                    comma = True))
                rare_formatted = ''.join(rare_formatted)
    
            elif number_rare_tax > 0 and len(unique) > 0:
                rare_formatted = ['This sample included the following rare or'
                ' \\textcolor{red}{unique} taxa: ']
                rare_formatted.append(convert_taxa_to_list(rare_combined, 
                                                    tax_format = rare_format,
                                                    render_mode = RENDERING, 
                                                    comma = True))
                rare_formatted = ''.join(rare_formatted)
            
            else:
                rare_formatted = "There were no rare or unique taxa found"\
                             " in this sample."

            # Saves the file
            file_name = pjoin(output_dir, '%s%s%s' % (FILE_PRECURSER, sample_id, 
                FILE_EXTENSION))

            file_for_editing = open(file_name, 'w')
            # file_for_editing.write('% Participant Name\n\\def\\yourname'\
            #     '{Michael Pollan or longer name}\n\n')
            file_for_editing.write('%% Abundance Table\n%s\n\n\n' \
                % abundance_formatted)
            file_for_editing.write('%% Enrichment Table\n%s\n\n\n' \
                % high_formatted)
            file_for_editing.write('%% Rare List\n\\def\\rareList{%s}\n' \
                % rare_formatted)
            file_for_editing.close()

# Sets up command line parsing
parser = ArgumentParser(description = "Creates lists and tables of enriched, abundance and rare taxa")

parser.add_argument('-i', '--input', help = 'Path to taxonomy table [REQUIRED]')
parser.add_argument('-o', '--output', \
                    help = 'Path to the output directory [REQUIRED]')
parser.add_argument('-s', '--samples', default = None, \
                    help = 'Sample IDs to be analyzed. If no value is '\
                    'specified, all samples in the taxonomy file will be'\
                    ' analyzed.')

if __name__ == '__main__':

    args = parser.parse_args()

    # Checks the tax table file path is sane and loads it.
    if not args.input:
        parser.error('An input taxonomy table is required')
    elif not isfile(args.input):
        parser.error("The supplied taxonomy file does not exist in the path.")
    else:
        (taxa, table, sample_ids) = taxa_importer(args.input)

    # Checks the output directory is sane.
    if not args.output:
        parser.error('An output directory must be supplied.')
    elif not exists(args.output):
        mkdir(args.output)        

    output_dir = args.output
    

    # Parses the sample IDs as a list
    if args.samples:
        samples_to_analyze = []
        for sample in args.samples.split(','): 
            samples_to_analyze.append(sample)
    else:
        samples_to_analyze = None

    main(taxa = taxa, table = table, \
        sample_ids = sample_ids, output_dir = output_dir, \
        samples_to_analyze = samples_to_analyze)
