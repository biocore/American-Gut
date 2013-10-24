#!/usr/bin/env python

from os import mkdir
from os.path import isfile, exists
from numpy import (loadtxt, array, delete, mean, shape, argsort, sort)
from scipy.stats import ttest_1samp
from argparse import ArgumentParser

__author__ = "Justine Debelius"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Justine Debelius"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Justine Debelius"
__email__ = "j.debelius@gmail.com"

def taxa_importer(taxa_table_fp):
    """Loads text taxonomy files as numpy arrays.

    INPUT:
        taxa_table_fp -- a string describing the location of the taxonomy file

    OUTPUTS:
        taxonomy -- a numpy vector with greengenes taxonomy strings

        tax_table -- a numpy array with the relative frequencies of taxonomies
            (rows) for each give sample (column)

        sample_ids -- a numpy vector of sample ids associated with the 
            tax_table values
    """
    # Loads the file as strings to determine the shape, then pulls of 
    # taxonomic designations and Sample IDs.
    tax_table = loadtxt(taxa_table_fp, dtype = "string", \
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

def calculate_tax_rank_1(sample, population, taxa):
    """Identifies unique and rare samples in the population and preforms a 
    case 1 t-test on common samples.

    INPUTS:
        sample -- a one dimensional numpy array containing the taxonomic
                    frequency values for a single sample

        population -- a numpy array containing containing taxonomic frequency
                    values. Samples are columns, taxa are rows.

        taxa -- an array of greengenes ids associated the sample and 
                    population frequencies

    OUTPUTS:
        unique_taxa -- a list of greengenes taxonomy strings that are found
                    only in the sample

        rare_taxa -- a list of greengenes taxonomy strings that are found in
                    the sample and ten percent of the population

        high_taxa -- a dictionary of greengenes strings, keyed to the sample 
                    frequency, average population frequency, the ratio of 
                    values, and the p-value

        low_taxa -- a dictionary of greengene strings keyed to the sample 
                    frequency, average population frequency, the ratio of 
                    values, and the p-value
    """
    # Rare taxa are defined as appearing in less than 10% of the samples
    RARE_THRESHHOLD = 0.1

    (num_taxa, num_samples) = shape(population)
    sample_shape = sample.shape

    # Calculates binary and count matrices
    sample_bin = sample > 0
    population_bin = population > 0
    population_count = population_bin.sum(1)
   
    # Identifies unique taxa and removes them from the table    
    unique = []
    rare = []
    new_taxa = []
    remove_index = []    

    # Identifies rare and unique taxa
    for idx, taxon in enumerate(taxa):
        if sample_bin[idx] == 1 and population_count[idx] == 0:
            unique.append(taxon)
            remove_index.append(idx)

        elif sample_bin[idx] == 1 and \
             population_count[idx] < num_samples*RARE_THRESHHOLD:
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
    ratio = sample / population_mean
    # preforms a case 1 t-test comparing the sample and population
    (t_stat, p_stat) = ttest_1samp(population, sample, 1)
    # Preforms a bonferroni correction on the p values
    p_stat = p_stat*num_samples
    
    # Determines list position based on the smallest p values.
    p_sort = sort(p_stat)
    p_order = argsort(p_stat)

    # Goes through the p values and determines if they are enriched or depleted
    for index in p_order:
        if p_stat[index] < 0.05 and ratio[index] > 1:
            high.append([taxa[index], sample[index], population_mean[index], \
                ratio[index], p_stat[index]])

        elif p_stat[index] < 0.05 and ratio[index] > 0 and ratio[index] < 1:
            low.append([taxa[index], sample[index], population_mean[index], \
                ratio[index], p_stat[index]])
   
    return unique, rare, low, high

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

        formatted_taxa -- a string of strings with formatting for the final 
                    table. 
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
                new_element.append("%1.2f" % item)

            elif formatting_keys[idx] == "VAL_PER" and render_mode == "LATEX":
                new_element.append("%1.2f\\%%" % item)

            elif formatting_keys[idx] == "100_PER" and render_mode == "LATEX":
                new_element.append("%1.2f\\%%" % (item*100))

            elif formatting_keys[idx] == "VAL_PER":
                new_element.append("%1.2f%%" % item)

            elif formatting_keys[idx] == "100_PER":
                new_element.append("%1.2f%%" % (item*100))

        formatted_taxa.append(new_element)

    return formatted_taxa

def taxa_to_table(corr_taxa, header, render_mode="RAW", numbering = 1):
    """taxa_to_table takes a set of greengenes strings and corresponding 
    numeric values and converts them into a string or formatted table.

    INPUTS:
        corr_taxa -- a python dictonary that relates the green genes taxonomy 
                    string to the list of corresponding taxonomy values. 
                    Ideally, these should be strings already formatted for use
                    in the table.

        table_header -- a python list of strings that describe the volumns in 
                    the table. For raw text tables, the width of the columns is
                    set by the number of entries in the header. RAW tables can
                    have no more than 5 columns (corresponding to 15 characters
                    in each column).

        render_mode -- a string ("LATEX" or "RAW") which describes 
                    the way the table will be formatted. LATEX or HTML gives a
                    string containing formatting code. 
    
        numbering -- a binary value that will add numbers along the side of the
                    table if true. Automatically FALSE for TSV rendering.

    OUTPUTS:
        format_table -- a python string formatted to give a table of taxa when
                    rendered in the program specified by render_mode."""

    # Sets up the constant string designations, describing the phylogentic 
    # levels 
    TAX_DES = ['kingdom', 'phylum', 'class', 'order', 'family', \
    'genus', 'species']
    
    if render_mode == "LATEX":
        format_table = format_latex_table(corr_taxa, header, numbering, TAX_DES)
    else:
        format_table = format_raw_table(corr_taxa, header, numbering, TAX_DES)

    # Returns the formatted table string
    return format_table

def format_raw_table(corr_taxa, header, numbering, tax_des):
    """converts a greengenes ids and frequency to a raw text formatted table

    INPUTS:
        corr_taxa -- a list of lists that relates the green genes 
                    taxonomy string to the list of corresponding taxonomy 
                    values. Ideally, these should be strings already formatted 
                    for use in the table.

        table_header -- a python list of strings that describe the volumns 
                    in the table. For raw text tables, the width of the columns 
                    is set by the number of entries in the header. RAW tables 
                    can have no more than 5 columns (corresponding to 15 
                    characters in each column).

        render_mode -- a string ("LATEX", "HTML",  or "RAW") which describes 
                    the way the table will be formatted. LATEX or HTML gives
                     a string containing formatting code. 
    
        numbering -- a binary value that will add numbers along the side of 
                        the table if true. Automatically FALSE for TSV 
                        rendering.

        tax_des -- a list encoding phylogenetic levels

    OUTPUTS:
        format_table -- a python string formatted to give a table of taxa 
                    when rendered in the program specified by render_mode.
        """


    # Sets constants
    HEADER_BAR = "--------------------------------------------------------"\
                 "-------------------"
    SPACER = '                                    ' 
    TAX_SPACE = 27
    # Category lengths are set up for an 80 character table
    CATEGORY_LEN_2 = 47
    CATEGORY_LEN_3 = 23
    CATEGORY_LEN_4 = 15
    CATEGORY_LEN_5 = 11

    # Initializes the table
    header_len = len(header)

    # Sets up the category length based on the number of header categories
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

    # Creates initial table row
    if numbering == 1:
        header_element = ['-----%s\n      ' % HEADER_BAR]
    else:
        header_element = ['%s\n' % HEADER_BAR]

    # Adds the header description
    for counter, category in enumerate(header):
        clean_cat = "%s%s" % (category, spacer)
        if counter == 0:
            header_element.append('%s' % clean_cat[0:TAX_SPACE])
        else:
            header_element.append(clean_cat[0:category_len])

    # Terminates the table header
    if numbering == 1:
        header_element.append('\n-----%s' % HEADER_BAR)
    else:
        header_element.append('\n%s' % HEADER_BAR)

    format_table = [''.join(header_element)]

    for idx, element in enumerate(corr_taxa):
        otu = element[0]            
        otu_values = element[1:]
        
        # Cleans up otus for use in the table
        split_tax = [i.split('__',1)[-1] for i in otu.split('; ')]
        no_levels = len(split_tax)
        if no_levels < 7:
            clean_otu = "%s %s%s" \
            % (tax_des[no_levels - 1], \
                split_tax[no_levels - 1], SPACER)
        elif no_levels == 7:
            clean_otu = "%s %s%s" \
            % (split_tax[no_levels - 2], \
                split_tax[no_levels - 1], SPACER)
        else:
            clean_otu = "kingdom %s phylum Other%s" \
                % (split_tax, SPACER)
        
        # Sets up the first column
        if numbering == 1 and count < 10:
            table_row = ['\n ( %d) %s' % (count, clean_otu[0:27])]
        elif numbering == 1:
            table_row = ['\n (%d) %s' % (count, clean_otu[0:TAX_SPACE])]
        else:
            table_row = ['\n%s' % (clean_otu[0:TAX_SPACE])]
        
        # Adds row information to the table
        for value in otu_values:
            val_expand = "%s%s" % (value, SPACER)
            table_row.append("\t%s" % (val_expand[0:(category_len)]))
        
        table_row = ''.join(table_row)   
        format_table.append(table_row)

    # Terminates the table
    if numbering == 1:
        format_table.append('\n-----%s' % HEADER_BAR)
    else:
        format_table.append('\n%s' % HEADER_BAR)

    format_table = ''.join(format_table)

    return format_table

def format_latex_table(corr_taxa, header, numbering, tax_des):
    """converts a greengenes ids and frequency to a raw text formatted table
    
    INPUTS:
        corr_taxa -- a python dictonary that relates the green genes 
                    taxonomy string to the list of corresponding taxonomy 
                    values. Ideally, these should be strings already formatted 
                    for use in the table.

        table_header -- a python list of strings that describe the volumns 
                    in the table. For raw text tables, the width of the columns 
                    is set by the number of entries in the header. RAW tables 
                    can have no more than 5 columns (corresponding to 15 
                    characters in each column).

        render_mode -- a string ("LATEX", "HTML",  or "RAW") which describes 
                    the way the table will be formatted. LATEX or HTML gives
                    a string containing formatting code. 
    
        numbering -- a binary value that will add numbers along the side of 
                    the table if true. Automatically FALSE for TSV 
                    rendering.

        tax_des -- a list encoding phylogenetic levels

    OUTPUTS:
        format_table -- a python string formatted to give a table of taxa 
                        when rendered in the program specified by render_mode.
        """


    # Initializes the table
    if numbering == 1:
        format_elements = ["{r"]
        header_elements = [" & "]
    else:
        format_elements = ["{"]
        header_elements = [""]

    # Creates formatting designation and the header
    for category in header: 
        format_elements.append(" c")
        header_elements.append("%s & " % category)

    format_elements.append("}")        
    format_code = ''.join(format_elements)
    header_code = ''.join(header_elements)
    header_code = header_code[0:(len(header_code)-3)]
    # Sets up the intialization code
    format_table = ["\\begin{tabular}%s\n\\hline\n%s \\\\ \n\\hline\n" \
        %(format_code, header_code)]

    for idx, element in enumerate(corr_taxa):
        count = idx + 1
        otu = element[0]
        otu_values = element[1:]            
            
        # Cleans up the otu for use in the table
        split_tax = [i.split('__',1)[-1] for i in otu.strip().split('; ')]
        for id_, level in enumerate(split_tax):
            if level != '':                
                no_levels = id_

        if no_levels < 6:
            clean_otu = "%s %s" % (tax_des[no_levels - 1],  \
                split_tax[no_levels - 1])

        elif no_levels == 6:
            clean_otu = "%s \\textit{%s}" % (tax_des[no_levels - 1], \
                split_tax[no_levels - 1])
            
        elif no_levels == 7:
            clean_otu = "\\textit{%s %s}" \
            % (split_tax[no_levels - 2], \
               split_tax[no_levels - 1])
            
        else:
            clean_otu = "kingdom \\textit{%s} phylum Other" \
            % split_tax

        # Sets up the table entry
        table_row = []
        if numbering == 1 and count < 10:
            table_row.append("( %d) & %s" % (count, clean_otu))
            
        elif numbering == 1:
            table_row.append("(%d) & %s" % (count, clean_otu))
            
        else:
            table_row.append("%s" % clean_otu)
            
        # Adds numeric data to the table
        for value in otu_values:table_row.append(" & %s" % value)
                
        # Appends to create a line break
        table_row.append("\\\\\n")

        format_table.append(''.join(table_row))

    # Terminate the table
    format_table.append('\\hline\n\\end{tabular}')

    format_table = ''.join(format_table)

    return format_table

def render_latex_list(raw_taxa, tax_format, tax_des):
    """Creates a series of items for a LaTex encoded list
    INPUTS:
        raw_taxa -- a python list object containing taxonomy strings from 
                        greengenes to be included in the final, formated output.

        tax_format -- a list specifiying if an argument should be bolded 
                    (denoted by "BOLD") or left alone ("REG")

    OUTPUTS:
        format_list_items -- Latex encoded list items"""
    
    format_list_items = []       

    format_list_items = []

    for (idx, element) in enumerate(raw_taxa):
        split_tax = [i.split('__',1)[-1] for i in element.strip().split('; ')]
        for id_, level in enumerate(split_tax):
            if level != '':                
                no_levels = id_
        list_item = ['\n\\item ']           
        if tax_format[idx] == 'BOLD':
            if no_levels < 6:
                list_item.append("\\textbf{%s %s}" % (tax_des[no_levels], \
                    split_tax[no_levels]))
            elif no_levels == 6:
                list_item.append("\\textbf{%s \\textit{%s}}" \
                    % (tax_des[no_levels], \
                    split_tax[no_levels]))
            elif no_levels == 7:
                list_item.append("\\textbf{\\textit{%s %s}}" \
                    % (split_tax[no_levels - 1], \
                    split_tax[no_levels]))
            elif no_levels > 7:
                list_item.append("\\textbf{kingdom %s}" \
                    % split_tax)

        else:
            if no_levels < 6:
                list_item.append("%s %s" % (tax_des[no_levels], \
                    split_tax[no_levels]))
            elif no_levels == 6:
                list_item.append("%s \\textit{%s}" \
                    % (tax_des[no_levels], \
                    split_tax[no_levels]))
            elif no_levels == 7:
                list_item.append("\\textit{%s %s}" \
                    % (split_tax[no_levels - 1], \
                    split_tax[no_levels]))
            elif no_levels > 7:
                list_item.append("kingdom %s" % split_tax)
        
        list_item = ''.join(list_item)
        list_item = list_item.strip('[').strip(']')
        format_list_items.append(list_item)

    format_list_items = ''.join(format_list_items)

    return format_list_items
   
def render_raw_list(raw_taxa, tax_format, tax_des):
    """Creates a series of items for a raw text list
    INPUTS:
        raw_taxa -- a python list object containing taxonomy strings from 
                    greengenes to be included in the final, formated output.

        tax_format -- a list specifiying if an argument should be bolded 
                    (denoted by "BOLD") or left alone ("REG")

    OUTPUTS:
        format_list_items -- raw text encoded list items"""

    format_list_items = []

    for idx, element in enumerate(raw_taxa):
        split_tax = [i.split('__',1)[-1] for i in element.strip().split('; ')]
        for id_, level in enumerate(split_tax):
            if level != '':                
                no_levels = id_

        list_item = ['\n     o  ']
        if tax_format(idx) == 'BOLD':
            if no_levels < 7:
                list_item.append("*%s %s*" \
                    % (tax_des[no_levels - 1], \
                    split_tax[no_levels - 1]))
            elif no_levels == 7:
                list_item.append("*%s %s*" \
                    % (split_tax[no_levels - 2], \
                        split_tax[no_levels - 1]))
            else:
                list_item.append("*kingdom %s*" \
                    % split_tax)
        else:
            if no_levels < 7:
                list_item.append("%s %s" % (tax_des[no_levels - 1], \
                        split_tax[no_levels - 1]))
            elif no_levels == 7:
                list_item.append("%s %s" \
                    % (split_tax[no_levels - 2], \
                    split_tax[no_levels - 1]))
            else:
                list_item.append("kingdom %s" \
                    % split_tax)
            
        list_item = ''.join(list_item)    
        format_list_items.append(list_item)

    format_list_items = ''.join(format_list_items)

    return format_list_items

def taxa_to_list(raw_taxa, tax_format, render_mode="RAW"):
    """taxa_to_list takes a list of greengenes taxonomy strings and converts
    it to a text string that can be printed to a document.

    INPUTS:
        raw_taxa -- a python list object containing taxonomy strings from 
                    greengenes to be included in the final, formated output.

        tax_format -- a list specifiying if an argument should be bolded 
                        (denoted by "BOLD") or left alone ("REG")

        render_mode -- a python string describing the way the out should be 
                    formatted. Options are LATEX, corresponding to LaTex code,
                    HTML, or RAW. LaTEX or HTML will give a string of code for
                    inclusion in the corresponding document. RAW will give a 
                    text file suitable for viewing, although RAW formats do not
                    include italics.
    
    OUTPUTS:
        format_list -- a python string formatted to give a list of taxa 
                    according to the supplied formatting mode."""

    # Sets up the constant string designations, describing the phylogentic levels 
    TAX_DES = ['kingdom', 'phylum', 'class', 'order', 'family',\
    'genus', 'species']

    # Sets up precurser formatting text
    if render_mode == "LATEX":
        format_list = ["\\begin{itemize}"]
        format_list.append(render_latex_list(raw_taxa, tax_format, TAX_DES))
        format_list.append('\n\\end{itemize}')

    else:
        format_list = []
        format_list.append(render_raw_list(raw_taxa, tax_format, TAX_DES))
        format_list.append('\n')

    format_list = ''.join(format_list)

    # Returns formatted string
    return format_list

def generate_otu_signifigance_tables_AGP(taxa, table, samples, output_dir, \
    sample_ids = None):
    """Creates LaTeX formatted significant OTU lists

    INPUTS:
        taxa -- a numpy vector with greengenes taxonomy strings

        tax_table -- a numpy array with the relative frequencies of taxonomies
            (rows) for each give sample (column)

        sample_ids -- a numpy vector of sample ids associated with the 
            tax_table values

        output_dir -- a directory where the final files should be saved.

        sample_ids -- a list of sample_ids which should be used to generate 
                    data. If this is left empty, all the samples in the table 
                    will be used.

    OUTPUTS:
        Generates text files containing LaTex encoded strings which creates a 
        formatted table of taxa enriched in a single sample 
        (Table_<SAMPLE_ID>.txt) and a list of rare and unique samples 
        (List_<SAMPLE_ID>). Rare defined as present in less than 10% of the 
        total population. The unique taxa are bolded in the lists. 
    """
    # Sets table constants
    RENDERING = "LATEX"
    FORMAT_KEYS = ["100_PER", "100_PER", "VAL_INT", "SKIP"]
    TABLE_HEADER = ['Taxonomy', 'Sample', 'Population', 'Fold Difference']
    # Number of taxa shown is an indexing value, it is one less than what is 
    # actually shown.
    NUMBER_OF_TAXA_SHOWN = 4

    # Checks the output directory is sane

    # Sets up samples for which tables are being generated
    if sample_ids == None:
        samples_to_test = samples
    else:
        samples_to_test = sample_ids

    for idx, sample_id in enumerate(samples_to_test):
        # Sets up the sample and population sets
        population = delete(table, idx, 1)
        sample = table[:,idx]


        # Calculates tax rank tables
        (unique, rare, low, high) = calculate_tax_rank_1(sample, population, \
            taxa)

        # Generates formatted table
        formatted_high = convert_taxa(high[0:NUMBER_OF_TAXA_SHOWN], \
            render_mode = RENDERING, formatting_keys = FORMAT_KEYS)
        high_formatted = taxa_to_table(formatted_high, TABLE_HEADER, \
            render_mode = RENDERING)

        # Generates formatted list
        rare_format = []
        rare_combined = []
        for taxon in unique:
            rare_combined.append(taxon)
            rare_format.append('BOLD')
        for taxon in rare:
            rare_combined.append(taxon)
            rare_format.append('REG')

        number_rare_tax = len(rare_combined)

        if number_rare_tax > NUMBER_OF_TAXA_SHOWN + 1:
            rare_formatted = ["This sample contained %i rare or unique taxa,"\
                              " including the following.\\\n" % number_rare_tax]
            rare_formatted.append(taxa_to_list(\
                rare_combined[:NUMBER_OF_TAXA_SHOWN ], rare_format, RENDERING))
            rare_formatted = ''.join(rare_formatted)
    
        elif number_rare_tax > 0:
            rare_formatted = taxa_to_list(rare_combined, rare_format, \
                RENDERING)

        else:
            rare_formatted = "There were no rare or unique samples found"\
                         " in this sample."

        # Saves the file
        file_table_name = "%s/Table_%s.txt" % (output_dir, sample_id)
        file_list_name = "%s/List_%s.txt" % (output_dir, sample_id)

        file_table = open(file_table_name, 'w')
        file_table.write(high_formatted)
        file_table.close()

        file_list = open(file_list_name, 'w')
        file_list.write(rare_formatted)
        file_list.close()

#american_gut_fp = "/Users/jwdebelius/Desktop/FecesSplit/L6.txt"
#output_dir = "/Users/jwdebelius/Desktop/TestOut/"
#sample_ids = ['000007117.1075649', '000005634.1053886', '000005637.1053909']
#generate_otu_signifigance_tables_AGP(american_gut_fp, output_dir, \
#    sample_ids = sample_ids)

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
        raise ValueError, "The supplied taxonomy file does not exist in the path."
    else:
        (taxa, table, sample_ids) = taxa_importer(args.input)

    # Checks the output directory is sane.
    if not args.output:
        parser.error('An output directory must be supplied.')
    elif not exists(args.output):
        mkdir(args.output)
        output_dir = args.output

    if output_dir[-1] != "/":
        temp_dir_name = [output_dir]
        temp_dir_name.append('/')
        output_dir = ''.join(temp_dir_name)

    # Parses the sample IDs as a list
    if args.samples:
        samples_to_analyze = []
        for sample in args.samples.split(','): 
            samples_to_analyze.append(sample)
    else:
        samples_to_analyze = None

    generate_otu_signifigance_tables_AGP(taxa = taxa, table = table, \
        samples = sample_ids, output_dir = output_dir, \
        sample_ids = samples_to_analyze)
