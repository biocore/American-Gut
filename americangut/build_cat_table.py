#!/usr/bin/env python 
# build_cat_table.py

from numpy import (array, arange, zeros, mean, median, concatenate, std)
from itertools import chain
from TaxPlot import check_data_array


__author__ = "Justine Debelius"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Justine Debelius"]
__license__ = "BSD" 
__version__ = "unversioned"
__maintainer__ = "Justine Debelius"
__email__ = "Justine.Debelius@colorado.edu"


def check_category(mapping, category):
    """Checks that a category exists in a mapping file"""
    # Preforms a sanity check on the inputs
    if not isinstance(mapping, dict) or not \
        isinstance(mapping[mapping.keys()[0]], dict):
        raise TypeError('mapping must be a 2D dict.')
    if not isinstance(category, str):
        raise TypeError('category must be a string')

    # Checks the categories exist for all instances of the metadata
    possible_cats = set(mapping[mapping.keys()[0]].keys())
    for (key, val) in mapping.iteritems():
        level_cat = set(val.keys())
        if not level_cat == possible_cats:
            raise ValueError('The metadata categories are not consistent '
                             'between samples.')

    # Checks the category is possible
    if category not in possible_cats:
        raise ValueError('%s is not a category in the metadata file.' 
                         % category)

def identify_groups(mapping, category):
    """Determines the group identigies in a mapping category
    INPUTS:
        mapping -- a 2D dictionary of sample IDs keyed to mapping categories. 
                    The first key in each sample dictionary should be #SampleId
                    or #SampleID.
        category -- a category in the mapping dictionary
    
    OUTPUTS:
        groups -- a set of possible values for the mapping category
    """
    # Checks the category is sane
    check_category(mapping, category)

    # Assigns the groups to a set
    groups = set()
    for (sampleid, meta) in mapping.iteritems():
        groups = groups.union([meta[category]])
    
    return groups

def get_category_position(all_cats, cat_descr):
    """Determines the position of a category in a list of categories
    INPUTS:
        all_cats -- a list of category descriptions in the form of 
                    strings or tuples (all enteries must be of the 
                    same type).
                
        cat_descr -- a string that matches part or all of the desired 
                    category. This description must be an exact match 
                    and must be unique to the category.
                    
    OUTPUTS:
        cat_pos -- the position of the target category in all the 
                    categories
    """
    
    # Looks for the target category in the list of categories
    watch_pos = []
    for (pos, cat) in enumerate(all_cats):
        if isinstance(cat, tuple):
            cat_bool = [cat_descr in i for i in chain(cat)]
            if any(cat_bool):
                watch_pos.append(pos)
        elif isinstance(cat, str) and cat_descr in cat:
             watch_pos.append(pos)
        
    # Preforms a sanity check on the position
    num_pos = len(watch_pos)
    if not num_pos == 1:
        raise ValueError('A unique position cannot be identified.')
        
    return watch_pos[0]

def sort_alphabetically(unord):
    """Sorts the list in alphabetical order
    INPUTS:
        unordered -- a list of strings to be sorted

    OUTPUTS:
        alpha -- the list in alphabetical order
        a_order -- the order of the items in the list
    """
    pos_asign = [(val, loc) for (loc, val) in enumerate(unord)]
    pos_order = sorted(pos_asign)
    [alpha, a_order] = zip(*pos_order)
    alpha = list(alpha)
    a_order = list(a_order)
    print alpha
    print a_order

    return list(alpha), list(a_order)

def fuzzy_match(target, possible_matches):
    """Determines if the target phrase is in a unique entry in possible_matches

    INPUTS:
        target -- a string or integer which uniquely identifies the group.

        possible_matches -- the list of values to be identified

    OUTPUTS:
        pos -- the location of the match in the possible_matches

        match -- the actual value of the matched data
    """

    # Identifies the type of data to be matched
    watch_pos = []
    for pos, possible in enumerate(possible_matches):
        if not isinstance(possible, (str, tuple, list)):
            raise TypeError('fuzzy_match is for strings and iterable classes.')
        if isinstance(possible, str) and target in possible:
            watch_pos.append(pos)
            match = possible
        elif isinstance(possible, tuple) or isinstance(possible, list):
            pos_check = [target in chunk for chunk in possible]            
            if any(pos_check):
                watch_pos.append(pos)
                match = possible
        
    # Checks the match is sane
    if len(watch_pos) == 0:
        raise ValueError('No match could be found for the data.')
    elif len(watch_pos) > 1:
        raise ValueError('The matches were not unique.')
        
    pos = watch_pos[0]
    
    return pos, match

def identify_sample_group(sample_ids, mapping, category):
    """Determines which samples belong to each group of the category.

    INPUTS:
        sample_ids -- A list of the Sample IDs corresponding to the columns in 
                    data

        mapping -- a 2D dictionary of sample IDs keyed to mapping categories. 

        category -- a category in the mapping dictionary

    OUTPUTS:
        group_assign -- a dictionary keying the group name to the sample ids 
                in the group.

    """

    # Identifies the groups associated with the category
    groups = identify_groups(mapping, category)
    group_assign = {group:[] for group in groups}

    # Adds the sample id to the appropriate list
    for samp in sample_ids:
        try:
            samp_group = mapping[samp][category]
        except:
            raise ValueError('Not all the sample ids have mapping data.')
        group_assign[samp_group].append(samp)

    return group_assign

def build_sub_table(data, sample_ids, target_ids):
    """Gets a subsample of the data based on the sample ids

    INPUTS:
        data -- A numpy array of the abundance values. Rows correspond to the 
                    metadata category, and columns correspond to sample ids.

        sample_ids -- A list of the Sample IDs corresponding to the columns in 
                    data

    OUTPUTS:
    """
     # Preforms a sanity check on the inputs
    [num_rows, num_cols] = data.shape
    num_ids = len(sample_ids)
    if isinstance(target_ids, str):
        target_ids = [target_ids]
        num_targets = 1
    elif isinstance(target_ids, list):
        num_targets = len(target_ids)
    else:
        raise TypeError('target_ids must be a string or a list of strings')
        
    if not num_cols == num_ids:
        raise ValueError('Each column in data must have a corresponding sample'
            ' id')
  
    new_data = zeros((num_rows, num_targets))
    
    for (idx, t_id) in enumerate(target_ids):
        # Target ids must exist in sample ids
        if t_id not in sample_ids:
            raise ValueError('%s is not a valid id.' %t_id)
        sample_pos = sample_ids.index(t_id)
        data_vec = data[:,sample_pos]
        new_data[:,idx] = data_vec
        
    return new_data

def sort_samples(data, sample_ids, sort_key=None, reverse=False):
    """Sorts data columns and sample ids using the supplied key
    
    INPUTS:
        data -- A numpy array of the abundance values. Rows correspond to the 
                    metadata category, and columns correspond to sample ids.
        
        sample_ids -- A list of the Sample IDs corresponding to the columns in 
                    data
        
        sort_key -- A numeric value or indicating the location in the categories
                    to use for sorting the data, or a None value. None indicates
                    data should be sorted by the sample ids. Common Categories 
                    must be provided to sort numerically. 
                    DEFAULT: None 
                    
        reverse -- A binary value indicating if data should be sorted in 
                    ascending (False) or descending (True) order. 
                    DEFAULT: False
                    
    OUTPUTS:
        sorted_ids -- an ordered list of sample ids
        
        sorted_data -- an ordered numpy array of the data        
    """
    
    # Preforms a sanity check
    [num_rows, num_cols] = data.shape
    num_samp = len(sample_ids)
    if not num_cols == num_samp:
        raise ValueError('Each sample observation must have a corresponding '
            'sample id.')
    if sort_key is not None and sort_key > (num_rows-1):
        raise ValueError('The sort_key must be a row in data.')
        
    # Pulls out item to be sorted
    if sort_key is None:
        sort_vec = sample_ids
    else:
        sort_vec = list(data[sort_key,:])
        
    # Identifies the order for samples
    order = [loc for (val, loc) in sorted((val, loc) 
             for (loc, val) in enumerate(sort_vec))]
    
    if reverse:
        order = order[::-1]
            
    # Orders the sample_ids and data
    sorted_ids = [sample_ids[loc] for loc in order]
    sorted_data = data[:,order]
    
    return sorted_ids, sorted_data            

def sort_categories(data, categories, first_cat, sort_method='ALPHA_WRAP', 
    category_order=None, reverse=False):
    """Sorts data rows and categories using the supplied information
    INPUTS:
        data -- A numpy array of the abundance values. Rows correspond to the 
                    metadata category, and columns correspond to sample ids.
                    
        categories -- a list of the categories corresponding to the 
                    current order in the data table.
        
        first_cat -- the category which should appear first in the list of 
                    categories when randomly ordered. This may be a string, 
                    tuple, or list.

        sort_method -- a string describing how the data should be sorted. 
                    "ALPHA" says start at the first_cat value, and follow with 
                        the categories alphabetically, excluding the first 
                        category. For example, if the groups are DOG, CAT, BIRD,
                        SMALL, and FISH, and the first_cat is DOG, the order 
                        would be DOG, BIRD, CAT, FISH, SMALL.
                    "ALPHA_WRAP" starts at the first_cat value and sorts 
                        alphabetically following the first category. Using the 
                        same set of parameters as above, the final order would 
                        be DOG, FISH, SMALL, BIRD, CAT.
                    "RETAIN" starts at the first category, and then keeps the 
                        same order which was input before. (DOG, CAT, BIRD, 
                        SMALL, FISH).
                    "CUSTOM" uses the category_order provided.
                    DEFAULT: 'ALPHA_WRAP'
        
        category_order -- a list giving the values in the custom desired order.

        reverse -- A binary value indicating if data should be sorted in 
                    ascending (False) or descending (True) order. 
                    DEFAULT: False
    
    OUTPUTS:
        sorted_cats -- a list of the categories sorted according to the 
                    parameters and identifying the rows in data

        sorted_data -- a numpy array where the rows are sorted to correspond 
                    with sorted_cats
    """
    
    METHODS = set(['ALPHA', 'ALPHA_WRAP', 'RETAIN', 'CUSTOM'])

    # Preforms a sanity check
    [num_rows, num_cols] = data.shape
    num_cats = len(categories)
    if not num_rows == num_cats:
        raise ValueError('Each category must correspond to an observation in '
            'the data.')

    if sort_method not in METHODS:
        raise ValueError('The sort method is not supported')

    if sort_method is 'CUSTOM' and category_order is None:
        raise ValueError('A custom category order must be supplied.')

    if category_order is not None and len(category_order) > num_cats:
        raise ValueError('We cannot sort more categories than exist in the '
            'data table.')

        
    # Determines if alphabetical sorting is required       
    if 'ALPHA' in sort_method:
       [alphabetical, cat_order] = sort_alphabetically(categories)
    else:
        cat_order = [loc for (loc, val) in enumerate(categories)] 

    # Determines the position of the first_cat
    first_pos = get_category_position(categories, first_cat)
    first_ord = cat_order.index(first_pos)
    
    # Determines the order for the categories
    if sort_method == 'ALPHA' or sort_method == 'RETAIN':
        cat_order.remove(first_pos)
        # Inserts the ordered element at the first position
        cat_order.insert(0, first_pos)
        order = cat_order        
        
    elif sort_method == 'ALPHA_WRAP':
        order = cat_order[first_ord:]        
        order.extend(cat_order[:(first_ord)])


    elif sort_method == 'CUSTOM':
        for cat in category_order:            
            order = [get_category_position(categories, cat) for cat in category_order]
    
    else:
        raise ValueError('The sorting method cannot be determined.')
    
    # Reverses the order if necessary
    if reverse:
        order = order[::-1]

    # Orders the sample data
    sorted_cats = [categories[loc] for loc in order]
    # print sorted_cats
    sorted_data = data[order,:]

    return sorted_cats, sorted_data

