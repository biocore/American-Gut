#!/usr/bin/env python

from argparse import ArgumentParser
from os import mkdir
from numpy import array, delete
from biom.parse import parse_biom_table
from os.path import isfile, exists, join as pjoin
from americangut.generate_otu_signifigance_tables import (calculate_abundance,
                                                          calculate_tax_rank_1,
                                                          convert_taxa,
                                                          convert_taxa_to_list,
                                                          clean_greengenes_string,
                                                          build_latex_macro,
                                                          format_date)
from americangut.taxtree import build_tree_from_taxontable, sample_rare_unique
from americangut.make_phyla_plots import map_to_2D_dict

__author__ = "Justine Debelius"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Justine Debelius", "Daniel McDonald"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Justine Debelius"
__email__ = "Justine.Debelius@colorado.edu"


def main(taxa_table, output_dir, mapping=None, samples_to_analyze=None):
    """Creates LaTeX formatted significant OTU lists

    INPUTS:
        tax_table -- a numpy array with the relative frequencies of taxonomies
            (rows) for each give sample (column)

        output_dir -- a directory where the final files should be saved.

        mapping -- a 2D dictionary of mapping data where the sample id is keyed
                    to a dictionary of metadata.

        samples_to_analyze -- a list of samples_to_analyze which should be used
                    to generate data. If None, all the samples will be used.
                    DEFAULT: None

    OUTPUTS:
        Generates text files containing LaTex encoded strings which creates a
        LaTeX macro dictionary with the information for creating a table of
        most abundant taxa, most enriched taxa, and rare and unique taxa. Rare
        defined as present in less than 10% of the total population. The unique
        taxa are bolded in the lists.
    """
     # Sets up the way samples should be converted
    SAMPLE_CONVERTER = {'feces': 'fecal',
                        'oral_cavity': 'oral',
                        'skin': 'skin'}

    DUMMY = ['', '', '', '']
    COUNT = [0, 1, 2, 3, 4, 5, 6, 7]

     # Sets up the way samples should be converted
    SAMPLE_CONVERTER = {'feces': 'fecal',
                        'oral_cavity': 'oral',
                        'skin': 'skin'}

    DUMMY = ['', '', '', '']
    COUNT = [0, 1, 2, 3, 4, 5, 6, 7]
    # Sets table constants
    RENDERING = "LATEX"
    RARE_THRESH = 0.1

    SUM_MIN = 1

    FORMAT_SIGNIFIGANCE = ['%1.2f', "%1.2f", "%i", "SKIP"]
    SIGNIFIGANCE_HUNDRED = [True, True, False, False]
    MACRO_CATS_SIGNIFICANCE = ['enrichTaxon', 'enrichSampl', 'enrichPopul',
                               'enrichFold']
    MACRO_FORM_SIGNIFICANCE = [lambda x: clean_greengenes_string(x,
                               render_mode='LATEX'),
                               lambda x: x,
                               lambda x: x,
                               lambda x: x]

    DUMMY = ['', '', '', '']
    COUNT = [0, 1, 2, 3, 4, 5, 6, 7]

    FORMAT_ABUNDANCE = ["%1.1f"]
    ABUNDANCE_HUNDRED = [True]
    MACRO_CATS_ABUNDANCE = ['abundTaxon', 'abundSampl']
    MACRO_FORM_ABUNDANCE = [lambda x: clean_greengenes_string(x,
                            render_mode='LATEX'), lambda x: x]

    FILE_PRECURSER = 'macros_'
    FILE_EXTENSION = '.tex'

    DATE_FIELD = 'COLLECTION_DATE'
    DATE_FORMAT_SHORT = '%m/%d/%y'
    DATE_FORMAT_LONG = '%m/%d/%Y'

    UNKNOWNS = set(['None', 'NONE', 'none', 'NA', 'na', 'UNKNOWN', 'unknown'])
    DATE_OUT = '%B %d, %Y'
    TIME_FIELD = 'SAMPLE_TIME'

    # Number of taxa shown is an indexing value, it is one less than what is
    # actually shown.
    NUM_TAXA_SHOW = 5

    # Builds the the taxomnomy tree for the table and identifies the
    # rare/unique taxa in each sample
    tree, all_taxa = build_tree_from_taxontable(taxa_table)

    # Sets up samples for which tables are being generated
    if not samples_to_analyze is None:
        samples_to_test = samples_to_analyze
    else:
        samples_to_test = all_taxa.keys()

    if samples_to_test:
        samples_to_test = set(samples_to_test)
        tmp = {k: v for k, v in all_taxa.items() if k in samples_to_test}
        all_taxa = tmp
        if not samples_to_test:
            raise ValueError("No samples!")

    # Generates lists and tables for each sample
    for samp, filtered_table, rare, unique in sample_rare_unique(tree,
                                                                 tax_table,
                                                                 all_taxa,
                                                                 RARE_THRESH):
        # Sets up filename
        file_name = pjoin(output_dir, '%s%s%s' % (FILE_PRECURSER, samp,
                          FILE_EXTENSION))

        filt_fun = lambda v, i, md: v.sum() > 0
        filtered_table = filtered_table.filterObservations(filt_fun)
        abund_table = tax_table.filterObservations(filt_fun)

        # Gets sample information for the whole table
        abund_sample = abund_table.sampleData(samp)
        abund_taxa = abund_table.ObservationIds

        # Gets sample information for other filtered samples
        filt_taxa = filtered_table.ObservationIds
        population = array([filtered_table.observationData(i) for i in
                            filtered_table.ObservationIds])

        sample_position = filtered_table.getSampleIndex(samp)
        filt_sample = filtered_table.sampleData(samp)

        population = delete(population, sample_position, 1)

        # Converts the lists into greengenes strings for later processing
        greengenes_rare = []
        greengenes_unique = []
        for taxon in rare:
            greengenes_rare.append(';'.join(taxon))
        for taxon in unique:
            greengenes_unique.append(';'.join(taxon))

        # Formats the rare and unique lists
        rare_format = []
        rare_combined = []
        for taxon in greengenes_unique:
            rare_combined.append(taxon)
            rare_format.append('COLOR')
        for taxon in greengenes_rare:
            rare_combined.append(taxon)
            rare_format.append('REG')

        number_rare_tax = len(rare_combined)
        num_rare = len(rare)
        num_unique = len(unique)

        rare_formatted = \
            convert_taxa_to_list(rare_combined[0:NUM_TAXA_SHOW],
                                 tax_format=rare_format,
                                 render_mode=RENDERING,
                                 comma=True)

        if num_unique > 0:
            unique_string = ' and \\textcolor{red}{%i unique}' % num_unique
        else:
            unique_string = ''

        if number_rare_tax == 0:
            rare_formatted = "There were no rare or unique taxa found in "\
                "your sample."

        elif 0 < number_rare_tax <= NUM_TAXA_SHOW:
            rare_formatted = 'Your sample contained the following rare%s '\
                'taxa: %s.' % (unique_string, rare_formatted)

        else:
            rare_formatted = 'Your sample contained %i rare%s taxa, '\
                'including the following: %s.' \
                % (num_rare, unique_string,
                   rare_formatted)

        # Calculates abundance rank
        (abundance) = calculate_abundance(abund_sample, abund_taxa,
                                          sum_min=SUM_MIN)

        # Generates formatted abundance table
        formatted_abundance = convert_taxa(abundance[0:NUM_TAXA_SHOW],
                                           formatting_keys=FORMAT_ABUNDANCE,
                                           hundredx=ABUNDANCE_HUNDRED)

        abundance_formatted = \
            build_latex_macro(formatted_abundance,
                              categories=MACRO_CATS_ABUNDANCE,
                              format=MACRO_FORM_ABUNDANCE)

        (high, low) = calculate_tax_rank_1(sample=filt_sample,
                                           population=population,
                                           taxa=filt_taxa,
                                           critical_value=0.05)

        if len(high) == 0:
            formatted_high = [['', '', '', '']]*NUM_TAXA_SHOW

        if len(high) < NUM_TAXA_SHOW:
            # Formats the known high taxa
            formatted_high = \
                convert_taxa(high[0:NUM_TAXA_SHOW],
                             formatting_keys=FORMAT_SIGNIFIGANCE,
                             hundredx=SIGNIFIGANCE_HUNDRED)

            # Adds the dummy list to the end
            for idx in COUNT:
                if idx == (NUM_TAXA_SHOW - len(high)):
                    break
                formatted_high.append(DUMMY)

        else:
            formatted_high = convert_taxa(high[0:NUM_TAXA_SHOW],
                                          formatting_keys=FORMAT_SIGNIFIGANCE,
                                          hundredx=SIGNIFIGANCE_HUNDRED)

        high_formatted = build_latex_macro(formatted_high,
                                           categories=MACRO_CATS_SIGNIFICANCE,
                                           format=MACRO_FORM_SIGNIFICANCE)

       # Handles date parsing
        if mapping is not None and mapping[samp][DATE_FIELD] not in UNKNOWNS:
            try:
                sample_date = format_date(mapping[samp],
                                          date_field=DATE_FIELD,
                                          d_form_in=DATE_FORMAT_SHORT,
                                          format_out=DATE_OUT)
            except:
                sample_date = format_date(mapping[samp],
                                          date_field=DATE_FIELD,
                                          d_form_in=DATE_FORMAT_LONG,
                                          format_out=DATE_OUT)
        else:
            sample_date = 'unknown'

        # Removes a zero character from the date
        if ',' in sample_date and sample_date[sample_date.index(',')-2] == '0':
                zero_pos = sample_date.index(',')-2
                sample_date = ''.join([sample_date[:zero_pos],
                                       sample_date[zero_pos+1:]])

        else:
            sample_date = 'unknown'

        # Handles sample parsing
        if mapping is not None and mapping[samp][TIME_FIELD] not in UNKNOWNS:
            sample_time = mapping[samp][TIME_FIELD].lower()
        else:
            sample_time = 'unknown'

        if mapping is not None:
            sample_type_prelim = mapping[samp]['BODY_HABITAT'].split(':')[1]
            if sample_type_prelim in SAMPLE_CONVERTER:
                sample_type = SAMPLE_CONVERTER[sample_type_prelim]
            elif sample_type in UNKNOWNS:
                sample_time = 'unknown'
            else:
                sample_type = sample_type_prelim.lower()
        else:
            sample_type = 'unknown'

        # Saves the file
        file_for_editing = open(file_name, 'w')
        file_for_editing.write('%% Barcode\n\\def\\barcode{%s}\n\n'
                               % samp.split('.')[0])
        file_for_editing.write('%% Sample Type\n\\def\\sampletype{%s}\n\n'
                               % sample_type)
        file_for_editing.write('%% Sample Date\n\\def\\sampledate{%s}\n'
                               '\\def\\sampletime{%s}\n\n\n'
                               % (sample_date, sample_time))
        file_for_editing.write('%% Abundance Table\n%s\n\n\n'
                               % abundance_formatted)
        file_for_editing.write('%% Enrichment Table\n%s\n\n\n'
                               % high_formatted)
        file_for_editing.write('%% Rare List\n\\def\\rareList{%s}\n'
                               % rare_formatted)
        file_for_editing.close()

# Sets up command line parsing
parser = ArgumentParser(description="Creates lists and tables of enriched, "
                        "abundance and rare taxa")

parser.add_argument('-i', '--input',
                    help='Path to taxonomy table [REQUIRED]')
parser.add_argument('-o', '--output',
                    help='Path to the output directory [REQUIRED]')
parser.add_argument('-s', '--samples',
                    default=None,
                    help='Sample IDs to be analyzed. If no value is '
                    'specified, all samples in the taxonomy file will be'
                    ' analyzed.')
parser.add_argument('-m', '--mapping',
                    default=None,
                    help='Path to a mapping file')

if __name__ == '__main__':

    args = parser.parse_args()

    # Checks the tax table file path is sane and loads it.
    if not args.input:
        parser.error('An input taxonomy table is required')
    elif not isfile(args.input):
        parser.error("The supplied taxonomy file does not exist in the path.")
    else:
        tax_table = parse_biom_table(open(args.input))

    # Checks the output directory is sane.
    if not args.output:
        parser.error('An output directory must be supplied.')
    elif not exists(args.output):
        mkdir(args.output)

    output_dir = args.output

    if args.mapping and not isfile(args.mapping):
        parser.error('The supplied mapping file does not exist in the path.')
    elif args.mapping:
        mapping = map_to_2D_dict(open(args.mapping, 'U'))
    else:
        mapping = args.mapping

    # Parses the sample IDs as a list
    if args.samples:
        samples_to_analyze = []
        for sample in args.samples.split(','):
            samples_to_analyze.append(sample)
    else:
        samples_to_analyze = None

    main(taxa_table=tax_table,
         output_dir=output_dir,
         mapping=mapping,
         samples_to_analyze=samples_to_analyze)
