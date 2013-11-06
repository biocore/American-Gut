#!/usr/bin/env python

from os import mkdir
from os.path import isfile, exists, join as pjoin
from numpy import loadtxt, delete, mean, shape, argsort, sort
from scipy.stats import ttest_1samp
from argparse import ArgumentParser

__author__ = "Justine Debelius"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Justine Debelius"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Justine Debelius"
__email__ = "j.debelius@gmail.com" 

def taxa_importer(taxa_table_file):
    """Loads text taxonomy files as numpy arrays.

    INPUT:
        taxa_table_file -- the contents of the open taxonomy file

    OUTPUTS:
        taxonomy -- a numpy vector with greengenes taxonomy strings

        tax_table -- a numpy array with the relative frequencies of taxonomies
            (rows) for each give sample (column)

        sample_ids -- a numpy vector of sample ids associated with the 
            tax_table values
    """
    # Loads the file as strings to determine the shape, then pulls of 
    # taxonomic designations and Sample IDs.
    tax_table = loadtxt(taxa_table_file, dtype = "string", \
        comments = "'", delimiter = "\t")
    (num_rows, num_cols) = tax_table.shape

    # The sample Ids are taken as the first row
    sample_ids = tax_table[0,1:]    
    taxonomy = tax_table[1:,0]       
    tax_table = delete(tax_table, 0, 0)
    tax_table = delete(tax_table, 0, 1)
    tax_table = tax_table.astype('float')
    # Returns the result
    return taxonomy, tax_table, sample_ids

def calculate_tax_rank_1(sample, population, taxa, rare_threshold=0.1, \
    abundance_threshold=0.95):
    """Identifies unique and rare samples in the population and preforms a 
    case 1 t-test on common samples.

    INPUTS:
        sample -- a one dimensional numpy array containing the taxonomic
                    frequency values for a single sample

        population -- a numpy array containing containing taxonomic frequency
                    values. Samples are columns, taxa are rows.

        taxa -- an array of greengenes ids associated the sample and 
                    population frequencies

        rare_threshhold -- a value between 0 and 1 indicating the maximum 
                    fraction of the population than can have an OTU for that 
                    OTU to be considered rare.

        abundance_threshhold -- a value between 0 and 1 indicating the minimum 
                    fraction of a sample to be represented by the most abundant 
                    OTUs. 

    OUTPUTS:
        unique -- a list of greengenes taxonomy strings that are found
                    only in the sample

        rare -- a list of greengenes taxonomy strings that are found in
                    the sample and ten percent of the population

        high -- a list of lists with greengenes strings, sample frequency, 
                    average population frequency, the ratio of values, and the 
                    p-value

        low -- a list of lists with greengenes strings, sample frequency, 
                    average population frequency, the ratio of values, and the 
                    p-value

        absent -- a list of greengenes taxonomy strings that are not found 
                    in the sample but are present in the population.

        abundant -- a list of lists of greenegenes taxonomy strings and the
                    frequencies representing the most abundant taxa in the 
                    sample.
    """
    # Rare taxa are defined as appearing in less than 10% of the samples

    (num_taxa, num_samples) = shape(population)

    # Calculates binary and count matrices
    sample_bin = sample > 0
    population_bin = population > 0
    population_count = population_bin.sum(1)

    # Identifies absent taxonomies in the sample
    absent = taxa[((sample_bin == 0) * (population_count != 0)) == 1]

    # Sorts the sample by abundance
    abundance_data = sort(sample)
    abundance_rank = argsort(sample)
    abundance_taxa = taxa[abundance_rank]
   
    # Identifies the taxonomy up to the abundance threshold    
    abundance_watch = 0
    abundant = []

    for idx, frequency in enumerate(reversed(abundance_data)):
        tax_idx = len(abundance_data) - (idx + 1)
        abundance_watch = abundance_watch + frequency
        abundant.append([abundance_taxa[tax_idx], frequency])
        if abundance_watch > abundance_threshold:
            break

    # Identifies unique taxa and removes them from the table    
    unique = []
    rare = []
    remove_index = []    

    # Identifies rare and unique taxa
    for idx, taxon in enumerate(taxa):
        if sample_bin[idx] == 1 and population_count[idx] == 0:
            unique.append(taxon)
            remove_index.append(idx)

        elif sample_bin[idx] == 1 and \
             population_count[idx] < num_samples*rare_threshold:
            rare.append(taxon)
            remove_index.append(idx)

    # Removes taxa identified as unique from the available set
    taxa = delete(taxa, remove_index)
    sample_bin = delete(sample_bin, remove_index)
    sample = delete(sample, remove_index)
    population_count = delete(population_count, remove_index)
    population = delete(population, remove_index, 0)

    # Identifies taxa that are significantly enriched or depleted in the 
    # population
    high = []
    low = []
    # Determines the ratio 
    population_mean = mean(population,1)
    ratio = sample.astype(float) / population_mean.astype(float)
    # preforms a case 1 t-test comparing the sample and population
    (t_stat, p_stat) = ttest_1samp(population, sample, 1)
    # Preforms a bonferroni correction on the p values
    p_stat = p_stat*num_samples
    
    # Determines list position based on the smallest p values.
    p_order = argsort(p_stat)

    # Goes through the p values and determines if they are enriched or depleted
    for index in p_order:
        if p_stat[index] < 0.05 and ratio[index] > 1:
            high.append([taxa[index], sample[index], population_mean[index], \
                ratio[index], p_stat[index]])

        elif p_stat[index] < 0.05 and ratio[index] > 0 and ratio[index] < 1:
            low.append([taxa[index], sample[index], population_mean[index], \
                ratio[index], p_stat[index]])
   
    return unique, rare, low, high, absent, abundant

def convert_taxa(rough_taxa, render_mode, formatting_keys):
    """Takes a dictionary of taxonomy and corresponding values and formats
    for inclusion in an output table.

    INPUTS:

        rough_taxa -- a dictionary of greengenes taxonomy strings keyed to a
                    list of numeric values.

        render_mode -- a string describing the format for the table: "RAW",
                    "HTML" or "LATEX".

        formatting_keys -- a string for converting values to strings. "VAL_INT"
                    converts the value to a string with an integer, "VAL_FLOAT" 
                    gives a truncated floating value as a string, "VAL_PER" adds
                    a percent symbol to the end of the string, and "100_PER" 
                    multiplies the value by 100 percent and adds a percent 
                    symbol at the end.

    OUTPUTS:

        formatted_taxa -- a list of string with formatting for the final table. 
    """
    formatted_taxa = []

    for element in rough_taxa:
        taxon = element[0]
        element.pop(0)
        new_element = [taxon]
        for idx, item in enumerate(element):
            if formatting_keys[idx] == "VAL_INT":
                new_element.append("%i" % item)

            elif formatting_keys[idx] == "VAL_FLOAT":
                new_element.append("%1.1f" % item)

            elif formatting_keys[idx] == 'VAL_100':
                new_element.append('%1.1f' % (item*100))

            elif formatting_keys[idx] == 'VAL_100_DEC_ALIGN':
                seperate = '%1.1f' % (item*100)
                seperate = seperate.split('.')
                new_element.append(' & '.join(seperate))

            elif formatting_keys[idx] == "VAL_PER" and render_mode == "LATEX":
                new_element.append("%1.1f\\%%" % item)

            elif formatting_keys[idx] == "100_PER" and render_mode == "LATEX":
                new_element.append("%1.1f\\%%" % (item*100))

            elif formatting_keys[idx] == "VAL_PER":
                new_element.append("%1.1f%%" % item)

            elif formatting_keys[idx] == "100_PER":
                new_element.append("%1.1f%%" % (item*100))

        formatted_taxa.append(new_element)

    return formatted_taxa

def convert_taxa_to_list(raw_taxa, tax_format, render_mode, comma = False):
    """Converts a list of greengenes strings to a text encoded list for printing

    INPUTS:
        raw_taxa -- a python list object containing taxonomy strings from 
                    greengenes to be included in the final, formated output.

        tax_format -- a list specifiying if an argument should be bolded 
                        (denoted by "BOLD") or left alone ("REG")

        render_mode -- a python string describing the way the out should be 
                    formatted. Options are LATEX, corresponding to LaTex code,
                    or RAW. LATEX will give a string which encodes a table. 
                    RAW will give a text file suitable for viewing, although 
                    RAW formats do not include italics.

        comma -- a binary value indicating whether the list should be single 
                    line comma separated list (TRUE), or a list format with each
                    item on its own line (FALSE).

    OUTPUT:
        format_list -- a python string formatted to give a list of taxa 
                    according to the supplied formatting mode."""

    # Sets up precurser text
    if render_mode == "LATEX":
        prelist = '\\begin{itemize}'
        antelist = '\n\\end{itemize}'
        preitem = '\n\\item '
        anteitem = ''
    else:
        prelist = ''
        antelist = ''
        preitem = '\n o   '
        anteitem = ''

    # Creates the list
    format_list = []
    if comma == True:
        for idx, taxon in enumerate(raw_taxa): 
            format_list.append(clean_otu_string(taxon, render_mode, \
                tax_format[idx].upper()))
        format_list = ', '.join(format_list)
    else:
        format_list.append(prelist)
        for idx, taxon in enumerate(raw_taxa):
            format_list.append('%s%s%s' % (preitem, clean_otu_string(taxon, \
                render_mode, tax_format[idx].upper()), anteitem))
        format_list.append(antelist)
        format_list = ''.join(format_list)

    return format_list
    
def clean_otu_string(greengenes_string, render_mode, format=False):
    """Distills a greengenes string to its high taxonomic resolution

    INPUTS:
        greengenes_string -- a greengenes string describing taxonomy

        render_mode -- a string ("LATEX", "HTML" or "RAW") which describes 
                    the way the table will be formatted. LATEX or HTML gives a
                    string containing formatting code.
        bold -- a binary value indication if the output string should be bolded.
                    In raw text, *bold* is render with *.
    OUTPUTS:
        cleaned_taxon -- a formatted string describing taxonomic information

    """
    # Preallocates level descriptors
    TAX_DES = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 
               'Species']

    # Sets up text formatting strings
    if render_mode == "LATEX":
        italic_before = '\\textit{'
        italic_after = '}'
        bold_before = '\\textbf{'
        bold_after = '}'
        color_before = '\\textcolor{red}{'
        color_after = '}'

    else:
        italic_before = ''
        italic_after = ''
        bold_before = '*'
        bold_after = '*'
        color_before = '*'
        color_after = '*'

    # Splits the taxonomy at the ; and removes the designation header. 
    split_tax = [i.split('__',1)[-1] for i in \
        greengenes_string.strip().split('; ')]
    
    # Identifies the highest level of resolution at which taxonomy is defined
    for id_, level in enumerate(split_tax):
        if level != '':                
            no_levels = id_

    # Sets up taxonomy string
    if no_levels < 5:
        cleaned_taxon = '%s %s' % (TAX_DES[no_levels], split_tax[no_levels])
    elif no_levels == 5:
        cleaned_taxon = '%s %s%s%s' % (TAX_DES[no_levels], italic_before, \
            split_tax[no_levels], italic_after)
    elif no_levels == 6:
        cleaned_taxon = '%s%s %s%s' % (italic_before, split_tax[no_levels-1],
            split_tax[no_levels], italic_after)
    else:
        cleaned_taxon = 'Kingdom %s' % split_tax

    cleaned_taxon = cleaned_taxon.replace('[', '').replace(']', '')
    cleaned_taxon = cleaned_taxon.replace('_','-')

    # Bolds taxon if necessary
    if format == 'BOLD':
        cleaned_taxon = ''.join([bold_before, cleaned_taxon, bold_after])
    elif format == 'COLOR':
        cleaned_taxon = ''.join([color_before, cleaned_taxon, color_after])

    return cleaned_taxon

def render_latex_header(header, numbering = True, alignment = 'c'):
    """Creates the header text for a LaTeX-encoded list

    INPUTS:
        header -- a python list of strings that describe the columns in the 
                    table. There can be no more than 5 header columns in a 
                    table.

        numbering -- a binary value that will add numbers along the side of 
                    the table if true.

        alignment -- a string of list of strings describing how the latex 
                    columns should be aligned. If a string (single value) is 
                    provided, it will be used for all columns. A list must have 
                    the same number of elements as the header list. Options for 
                    this string include 'l' (left), 'r' (right), and 'c' 
                    (center). For more information on alignment in LaTex 
                    headers, see <http://en.wikibooks.org/wiki/LaTeX/Tables>.

    OUTPUT:
        table_header -- a string encoded to begin a LaTeX formatted table
    """
    
    # Sets up a column for numbering if desired
    if numbering == 1:
        format_elements = ['{r']
        format_header = [' ']
    else:
        format_elements = ['{']
        format_header = []

    # Sets up alignment string
    if isinstance(alignment, (str, unicode)):
        for element in header:
            format_elements.append(alignment)

    elif isinstance(alignment, list) and len(header) != len(alignment):
        raise ValueError, 'The header list and alignment list must have the'\
            ' same number of elements.'

    elif isinstance(alignment, list):
        format_elements.extend(alignment)       

    else:
        raise TypeError, 'Alignment must be a string or list type.'

    format_elements.append('}')

    # Sets up header
    format_header.extend(header)

    # Combines into a properly spaced header
    table_header = '\\begin{tabular}%s\n\\hline\n%s \\\\\n\\hline\n' \
        % (' '.join(format_elements), ' & '.join(format_header))

    return table_header

def render_raw_header(header, numbering, header_bar, spacer, tax_len, cat_len):
    """Creates the header text for a raw text encoded table
    """
    
    # Initializes the table row with the header bar
    if numbering:
        header_elements = ['-----%s\n     ' % header_bar]
    else:
        header_elements = ['%s\n' % header_bar]

    # Creates the header elements for the table
    for idx, category in enumerate(header):
        clean_cat = '%s%s' % (category, spacer)
        if idx == 0:
            header_elements.append(clean_cat[:tax_len])
        else:
            header_elements.append(clean_cat[:cat_len])

    # Terminates the table header with a horizontal bar
    header_elements.append('\n%s\n' % header_bar)

    table_header = ''.join(header_elements)

    return table_header

def generate_latex_macro(corr_taxa, categories):
    """Generates a LaTeX macro for use in a template
    
    INPUTS:
        corr_taxa -- a list of lists where the first item is a greengenes 
                    taxonomy string and  subsequent items are associated 
                    descriptions.
        categories -- a list of strings describing the macro data to be 
                    substituted. i.e. Name, Sample, etc. 

    OUTPUT:
        A LaTeX macro with the definitions provided in categories


    """
    # Preallocates the an indexing variable
    ALPHABET = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 
                'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 
                'Y', 'Z']
    RENDER = 'LATEX'

    format_table = []

    # Combines categories with data and mapping index
    for idx, taxon_description in enumerate(corr_taxa):
        for id_, cat in enumerate(categories):
            if id_ == 0:
                format_table.append('\\def\\%s%s{%s}' % (cat, ALPHABET[idx], 
                                  clean_otu_string(taxon_description[0], 
                                  render_mode = RENDER)))
            else:
                format_table.append('\\def\\%s%s{%s}' % (cat, ALPHABET[idx], 
                                  taxon_description[id_]))

    # Inserts line breaks
    table = '\n'.join(format_table)


    return table

def convert_taxa_to_table(corr_taxa, header, render_mode = "RAW", \
    numbering = True, alignment = 'c', header_code=None):
    """Creates a text-encoded table

    INPUTS:
        corr_taxa -- a list of lists where the first item is a greengenes 
                    taxonomy string and  subsequent items are associated 
                    descriptions.

        header -- a python list of strings that describe the columns in 
                    the table. For raw text tables, the width of the columns is
                    set by the number of entries in the header. RAW tables can
                    have no more than 5 columns (corresponding to 15 characters
                    in each column). 
        
        render_mode -- a string ("LATEX", "HTML",  or "RAW") which describes 
                    the way the table will be formatted. LATEX or HTML gives a
                    string containing formatting code. 

        numbering -- a binary value that will add numbers along the side of the
                    table when True

        alignment -- a list or list of descriptors used to column alignment in 
                    LaTex encoded tables. This is only used when a table is 
                    formatted in LaTeX.
                    
        header_code -- a fully formatted LaTeX header description for custom 
                    formatting of the table
    OUTPUT:
        table -- a string encoding a table.

    """

    # Sets constants for table formatting
    HEADER_BAR = "------------------------------------------------------------"\
    "---------------"
    SPACER = '                                    '
    TAX_SPACE = 27

    # Category lengths are set up for an 80 character in a raw table
    CATEGORY_LEN_2 = 47
    CATEGORY_LEN_3 = 23
    CATEGORY_LEN_4 = 15
    CATEGORY_LEN_5 = 11

    # Sets up the the formatting key
    FORMAT_KEY = False

    # Sets up preformatted text by table type
    if render_mode == 'LATEX':  
        # Creates the table header
        if header_code == None:
            table_header = render_latex_header(header, numbering, \
                alignment = alignment)
        else:
            table_header = header_code

        # Creates formatting text around the table row
        anterow = '\\\\\n'  
        row_seperator = ' & '

        # Creates the text to termiante the table
        table_end = '\\hline\n\\end{tabular}'
        
    else:
        # Determines the appropriate category length for each row
        header_len = len(header)

        if header_len == 2:
            category_len = CATEGORY_LEN_2
        elif header_len == 3:
            category_len = CATEGORY_LEN_3
        elif header_len == 4:
            category_len = CATEGORY_LEN_4
        elif header_len == 5:
            category_len = CATEGORY_LEN_5
        else:
            raise ValueError, "There cannot be more than 5 header categories"

        # Creates the table header
        table_header = render_raw_header(header, numbering, \
            header_bar = HEADER_BAR, spacer = SPACER, tax_len = TAX_SPACE, \
            cat_len = category_len)

        # Creates the text around the table  
        anterow = '\n'
        row_seperator = ''

        # Sets up text to termiante the table
        if numbering == 1:
            table_end = '-----%s' % HEADER_BAR
        else:
            table_end = HEADER_BAR

    # Starts creating the table
    table_code = [table_header]

    # Creates the table rows
    for idx, taxon_description in enumerate(corr_taxa):
        # Sets up the row text
        table_row = []
        # Pulls out the taxon and descriptor, separates them and formats them
        taxon = taxon_description[0]
        clean_taxon = clean_otu_string(taxon, render_mode = render_mode, \
            bold = FORMAT_KEY)
        description = taxon_description[1:]

        # Adds numbering to the begining of the row if appropriate
        if numbering == True and idx < 10:
            table_row.append('( %d) ' % (idx + 1))
        elif numbering == True:
            table_row.append('(%d) ' % (idx + 1))
        
        # Pads the raw row if necessary
        if render_mode == "RAW":
            clean_taxon = '%s%s' % (clean_taxon, SPACER)
            clean_taxon = clean_taxon[:TAX_SPACE]

        # Adds the taxonomy data to the row
        table_row.append('%s' % (clean_taxon))

        # Sets up the data for the rows
        for item in description:
            # Pads the row if necessary
            if render_mode == "RAW":
                item = '%s%s' % (item, SPACER)
                item = item[:category_len]              

            # Adds the item to the table
            table_row.append('%s' % (item))

        # Inserts the correct spacer and adds to the table body
        table_code.append(row_seperator.join(table_row))

    table_code.append(table_end)
    
    # Returns an encoded string
    table = anterow.join(table_code)

    return table

def generate_otu_signifigance_tables_AGP(taxa, table, sample_ids, output_dir, \
    samples_to_analyze = None):
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
                rare_formatted = ["This sample contained %i rare "\
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
                rare_formatted = ['This sample included the follow rare taxa: ']
                rare_formatted.append(convert_taxa_to_list(rare_combined, 
                                                    tax_format = rare_format,
                                                    render_mode = RENDERING, 
                                                    comma = True))
                rare_formatted = ''.join(rare_formatted)
    
            elif number_rare_tax > 0 and len(unique) > 0:
                rare_formatted = ['This sample included the follow rare or'
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
parser = ArgumentParser(description = "Creates LaTeX formatted significant '\
                        'OTU lists and tables")

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

    generate_otu_signifigance_tables_AGP(taxa = taxa, table = table, \
        sample_ids = sample_ids, output_dir = output_dir, \
        samples_to_analyze = samples_to_analyze)
