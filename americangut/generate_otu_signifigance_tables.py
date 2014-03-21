#!/usr/bin/env python

from numpy import mean, shape, argsort, sort, sum as nsum, delete, seterr
from scipy.stats import ttest_1samp

__author__ = "Justine Debelius"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Justine Debelius"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Justine Debelius"
__email__ = "j.debelius@gmail.com" 

def calculate_abundance(sample, taxa, sum_min = 0.95):
    """Ranks taxa in a sample in order of abundance

    INPUTS:
        sample -- a one dimensional numpy array or list containing taxonomic 
                    frequencies in a single sample

        population -- a numpy array containing containing taxonomic frequency
                    values. Samples are columns, taxa are rows.

        taxa -- a one dimensional numpy array or list of greengenes ids 
                    associated the sample

        sum_min -- a value between 0 and 1 indicating the minimum fraction of a 
                    sample to be represented by the sum of the most abundant 
                    OTUs.

    OUTPUTS:
        abundant -- a list of lists of greenegenes taxonomy strings and the
                    frequencies representing the most abundant taxa in the 
                    sample.
    """ 

    if len(sample) != len(taxa):
        raise ValueError, 'The number of enteries in samples and taxa must be'\
            ' equal.'

    # Sorts the sample by abundance
    abundance_data = sort(sample)[::-1]
    abundance_rank = argsort(sample)[::-1]

    # List comprehension; faster. Also possible in dictionaries?
    abundance_taxa = [taxa[rank] for rank in abundance_rank]
   
    # Identifies the taxonomy up to the abundance threshold    
    abundance_watch = 0
    abundant = []
    for idx, frequency in enumerate(abundance_data):
        abundance_watch = abundance_watch + frequency
        abundant.append([abundance_taxa[idx], round(frequency, 6)])
        if abundance_watch > sum_min:
            break

    return abundant

def calculate_tax_rank_1(sample, population, taxa, critical_value = 0.05):
    """Identifies unique and rare samples in the population and preforms a 
    case 1 t-test on common samples.

    INPUTS:
        sample -- a one dimensional numpy array containing the taxonomic
                    frequency values for a single sample

        population -- a numpy array containing containing taxonomic frequency
                    values. Samples are columns, taxa are rows.

        taxa -- an array of greengenes ids associated the sample and 
                    population frequencies

        critical_value -- the alpha for use in the t-test


    OUTPUTS:       
        high -- a list of lists with greengenes strings, sample frequency, 
                    average population frequency, the ratio of values, and the 
                    p-value

        low -- a list of lists with greengenes strings, sample frequency, 
                    average population frequency, the ratio of values, and the 
                    p-value
    """

    # Rare taxa are defined as appearing in less than 10% of the samples

    (num_taxa, num_samples) = shape(population)
    
    if num_taxa != len(taxa):
        raise ValueError, 'The number of entries in samples and taxa must'\
            ' be equal.'

    # Identifies taxa that are significantly enriched or depleted in the 
    # population
    high = []
    low = []

    # Identifies taxa which are not populated
    population_count = nsum(population > 0, axis = 1)

    for idx, count in enumerate(population_count):
        # Removes any line which is equal to zero
        if count == 0:
            population = delete(population, idx, 0)           
            sample = delete(sample, idx)            
            taxa = delete(taxa, idx)           

    seterr(all='raise')
    # Determines the ratio 
    population_mean = mean(population,1)
    ratio = sample.astype(float) / population_mean.astype(float)
    # preforms a case 1 t-test comparing the sample and population
    t_stat = []
    p_stat = []
    # Could potentially use qiime functions
    (t_stat, p_stat) = ttest_1samp(population, sample, 1)

    # Preforms a bonferroni correction on the p values
    p_stat = p_stat*num_taxa

    # Determines list position based on the smallest p values.
    p_order = argsort(p_stat)

    # Goes through the p values and determines if they are enriched or depleted
    for index in p_order:
        if p_stat[index] >= critical_value:
            continue

        list_value = [taxa[index], 
                      round(sample[index], 6), 
                      round(population_mean[index],6),
                      round(ratio[index], 0), 
                      p_stat[index]]        
        if ratio[index] > 1:
            
            high.append(list_value)

        else:
            low.append(list_value)
   
    return high, low

def convert_taxa(rough_taxa, formatting_keys = '%1.2f', hundredx = False):
    """Takes a dictionary of taxonomy and corresponding values and formats
    for inclusion in an output table.

    INPUTS:

        rough_taxa -- a list of lists with a descriptor string followed by
                   a list of corresponding values

        formatting_keys --  a list describing the way values following the 
                    taxonomy string should be formatting using string formats. 
                    For example, %1.2f, %2d, %i. A value of 'SKIP' will ignore 
                    that value and remove it from the output list.

    OUTPUTS:

        formatted_taxa -- a list of string with formatting for the final table. 
    """
    
    print rough_taxa
    print rough_taxa.__class__
    
    # Checks the rough_taxa argument is sane
    if not isinstance(rough_taxa, list):
        raise TypeError('rough_taxa must be a list of lists')
    num_ent = len(rough_taxa[0])
    for entry in rough_taxa:
        if not isinstance(entry, list):
            raise TypeError('rough_taxa must be a list of lists')
        if not len(entry) == num_ent:
            raise ValueError('list size is inconsistant')
    num_rough = num_ent-1


    if isinstance(formatting_keys, list):
        num_keys = len(formatting_keys)
    else:
        num_keys = 1

    if isinstance(hundredx, list):
        num_hund = len(hundredx)
    else:
        num_hund = 1

    if not isinstance(formatting_keys, (list, str)):
        raise TypeError, 'formatting_keys must be a list or string.'
    if not num_rough == num_keys and isinstance(formatting_keys, list):
        raise ValueError('The number of elements in rough_taxa (%i) and the '
                         'number of elements in formatting_keys (%i) must be '
                         'equal.' %(num_rough, num_keys))

    elif not isinstance(hundredx, (list, bool)):
        raise TypeError, 'hundredx must be a list or bool.'
    elif not num_rough == num_hund and isinstance(hundredx, list):
        raise ValueError('The number of elements in rough_taxa(%i) and the '
                         'number of elements in hundredx(%i) must be equal.' 
                         %(num_rough, num_hund))

    # Converts formatting keys and hundredx to lists
    if isinstance(formatting_keys, str):
        formatting_keys = [formatting_keys]*num_rough

    if isinstance(hundredx, bool):
        hundredx = [hundredx]*num_rough

    # Creates formatted list
    formatted_taxa = []

    for element in rough_taxa:
        taxon = element[0]
        element.pop(0)
        new_element = [taxon]
        for idx, item in enumerate(element):

            if formatting_keys[idx] == 'SKIP':
                continue

            if hundredx[idx] == True:
                item = item * 100
            new_element.append(formatting_keys[idx] % item)

        formatted_taxa.append(new_element)

    return formatted_taxa

def convert_taxa_to_list(raw_taxa, tax_format, render_mode, comma = False, \
    color = 'red'):
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
            format_list.append(clean_otu_string(taxon, 
                                            render_mode = render_mode, \
                                            format = tax_format[idx].upper(), 
                                            unclassified = True,
                                            color = color))

        format_list = ', '.join(format_list)
    else:
        format_list.append(prelist)
        for idx, taxon in enumerate(raw_taxa):
            format_list.append('%s%s%s' % (preitem, 
                                clean_otu_string(taxon, \
                                render_mode, format = tax_format[idx].upper(),
                                unclassified = True,
                                color = color), 
                                anteitem))

        format_list.append(antelist)
        format_list = ''.join(format_list)

    return format_list
    
def clean_otu_string(greengenes_string, render_mode, format=False, \
    unclassified = False, color = 'red'):
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
        color_before = '\\textcolor{%s}{' % color
        color_after = '}'

    else:
        italic_before = ''
        italic_after = ''
        bold_before = '*'
        bold_after = '*'
        color_before = '*'
        color_after = '*'

    if unclassified == True:
        classified = 'Unclassified '
    else:
        classified = ''

    # Splits the taxonomy at the ; and removes the designation header. 
    split_tax = [i.split('__',1)[-1] for i in \
        greengenes_string.strip().split('; ')]
    
    # Identifies the highest level of resolution at which taxonomy is defined
    for id_, level in enumerate(split_tax):
        if level != '':                
            no_levels = id_

    # Sets up taxonomy string
    if no_levels < 5:
        cleaned_taxon = '%s%s %s' % (classified, TAX_DES[no_levels], \
            split_tax[no_levels])
    elif no_levels == 5:
        cleaned_taxon = '%s %s%s%s' % (TAX_DES[no_levels], italic_before, \
            split_tax[no_levels], italic_after)
    elif no_levels == 6:
        cleaned_taxon = '%s%s %s%s' % (italic_before, split_tax[no_levels-1],
            split_tax[no_levels], italic_after)
    else:
        cleaned_taxon = '%sKingdom %s' % (classified, split_tax)

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
    # Preallocates an indexing variable
    ALPHABET = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 
                'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 
                'Y', 'Z']
    # Rendering must be LaTeX for a LaTex macro
    RENDER = 'LATEX'

    format_table = []

    # Combines categories with data and mapping index
    for idx, taxon_description in enumerate(corr_taxa):
        if taxon_description[0] == '':
            for cat in categories:
                format_table.append('\\def\\%s%s{}' % (cat, ALPHABET[idx]))
                
        else:
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
            format = FORMAT_KEY)
        description = taxon_description[1:]

        # Adds numbering to the begining of the row if appropriate
        # Use string formatting in f
        if numbering:
            table_row.append('(%2d) ' % (idx + 1))
        
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
