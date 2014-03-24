#!/usr/bin/env python

__author__ = "Sam Way"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Sam Way"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Sam Way"
__email__ = "samuel.way@colorado.edu"

from numpy import array, argsort, zeros

def get_filtered_taxa_summary(mapping_file, taxa_summary_file, metadata_category, \
    metadata_value, top_n_taxa=7, select_taxa=None):
    """ Get a simplified taxonomy table.

        Inputs:
        mapping_file - Input mapping file (sample ids match taxa file)
        taxa_summary_file - Taxonomy summary file
        metadata_category - Used to select specific samples from the taxa file
        metadata_value - Value used to select specific samples from the taxa file
        top_n_taxa - If taxonomy groups aren't specified use the top N most abundant
        select_taxa - List of desired taxonomic groups

        Outputs:
        filtered_sample_ids - selected sample ids
        taxa_labels - taxonomic labels (including "Other")
        collapsed_taxa_table - simplified taxonomy table
    """

    mapping_fp = open(mapping_file, 'rU')
    mapping_dict, comments = parse_mapping_file_to_dict(mapping_fp)
    taxa_fp = open(taxa_summary_file, 'rU')
    sample_ids, taxa_ids, taxa_table = parse_taxa_summary_table(taxa_fp)
    taxa_ids = [ taxa_id.split('__')[-1] for taxa_id in taxa_ids ]

    selected_ids = [ key for key in mapping_dict.keys() if \
        mapping_dict[key][metadata_category] == metadata_value ]

    if len(selected_ids) < 1:
        raise ValueError('No sample ids match metadata_value="%s" in metadata_category="%s"' \
            % (metadata_value, metadata_category))

    sample_id_indices = [ i for i in xrange(len(sample_ids)) if sample_ids[i] in selected_ids ]
    filtered_taxa_table = taxa_table[:, sample_id_indices]
    filtered_sample_ids = [ sample_ids[idx] for idx in sample_id_indices ]

    if select_taxa is None:
        # If select_taxa is None, take the top N most abundant
        totals = filtered_taxa_table.sum(axis=1)
        taxa_indices = argsort(-totals)
        top_taxa = taxa_indices[:top_n_taxa]
        other_taxa = taxa_indices[top_n_taxa:]
        taxa_labels = [ taxa_ids[idx] for idx in top_taxa ]
    else:
        # List of taxa was supplied, use those
        top_taxa = [ taxa_ids.index(x) for x in select_taxa ]
        other_taxa = [ t for t in xrange(len(taxa_ids)) if t not in top_taxa ]
        taxa_labels = select_taxa

    taxa_labels.append('Other')
    N = len(taxa_labels) # Number of classes/labels
    M = filtered_taxa_table.shape[1] # Number of samples after filtering

    # Sort samples by most_abundant_taxa
    sort_sample_indices = argsort(-filtered_taxa_table[top_taxa[0], :])
    filtered_taxa_table = filtered_taxa_table[:, sort_sample_indices]
    filtered_sample_ids = [ filtered_sample_ids[idx] for idx in sort_sample_indices ]

    # Collapse "Others" rows into single row
    collapsed_taxa_table = zeros((N, filtered_taxa_table.shape[1]))
    collapsed_taxa_table[:-1, :] = filtered_taxa_table[top_taxa, :]
    collapsed_taxa_table[-1, :] = filtered_taxa_table[other_taxa, :].sum(axis=0)
    total = collapsed_taxa_table.sum(axis=0)
    collapsed_taxa_table = collapsed_taxa_table / total

    return filtered_sample_ids, taxa_labels, collapsed_taxa_table

def parse_taxa_summary_table(taxa_summary_fp, cast_as=float):
    """ Parse a taxa summary table 
        Returns tuple: sample_ids, otu_ids, matrix of OTUs(rows) by samples(cols)
    """ 
    sample_ids = []
    taxa_ids = []
    taxa_table = [] 

    header_line = taxa_summary_fp.readline()
    sample_ids = header_line.strip().split('\t')[1:]
    num_samples = len(sample_ids)
   
    for line in taxa_summary_fp:
        line = line.strip()
        if not line:
            continue
        line_pieces = line.split('\t')
        if len(line_pieces[1:]) != num_samples:
            raise ValueError("Error in taxa summary file - number of values does not " \
                "match the number of samples")
        taxa_table.append(array(map(cast_as, line_pieces[1:])))
        taxa_ids.append(line_pieces[0])

    return sample_ids, taxa_ids, array(taxa_table)

def parse_mapping_file_to_dict(mapping_file_fp):
    """ Takes an open mapping file and parses it to a dictionary structure """ 

    metadata_dict = {} 
    metadata_categories = [] 
    comments = [] 

    for line in mapping_file_fp:
        line = line.strip()
        if not line:
            continue 

        if line.startswith('#'):
            line = line[1:] 
            if not metadata_categories:
                metadata_categories = line.split('\t')[1:]  
                num_categories = len(metadata_categories)
            else:
                comments.append(line)
        else: 
            line_pieces = line.split('\t')
            if len(line_pieces[1:]) != num_categories:
                raise ValueError("Error in mapping file - number of metadata values does not " \
                    "match the number of metadata categories") 

            sample_id = line_pieces[0]
            metadata_values = line_pieces[1:]
            metadata_dict[sample_id] = { key:value for key, value in \
                zip(metadata_categories, metadata_values) }

    return metadata_dict, comments
    
