#!/usr/bin/env python

from biom.parse import parse_biom_table

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Daniel McDonald", "Justine Debelius"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"


def create_node(name):
    """Create and return a new node"""
    return {'name': name, 'popcount': 0, 'children': []}


def add_node(cur_node, node):
    """Add a node as a child to an existing node"""
    cur_node['children'].append(node)


def get_node(cur_node, name):
    """Fetch a node from an existing node, return None if not found"""
    for c in cur_node['children']:
        if c['name'] == name:
            return c
    return None


def update_tree(tree, taxa_by_sample):
    """Updates a nightmare nested tree structure

    all nodes have the following keys:
        name     : node name
        popcount : number of times observed in the population
        children : list of nodes
    """
    if tree is None:
        tree = create_node('root')

    for list_of_taxa in taxa_by_sample:
        tree['popcount'] += 1
        seen = set([])
        for tax_string in list_of_taxa:
            cur_node = tree

            for taxon in tax_string:
                if taxon.endswith('__') or '[' in taxon:
                    break

                node = get_node(cur_node, taxon)

                if node is None:
                    node = create_node(taxon)
                    add_node(cur_node, node)

                if taxon not in seen:
                    node['popcount'] += 1
                    seen.add(taxon)

                cur_node = node
    return tree


def get_rare_unique(tree, sample_taxa, rare_threshold):
    """Returns the rare and unique taxa in a sample

    sample_taxa    : a list of taxon strings
    rare_threshold : the level at which a taxon is considered rare
    """
    unique = []
    rare = []

    popsize = float(tree['popcount'])
    for tax_string in sample_taxa:
        cur_node = tree

        for taxon in tax_string:
            if taxon.endswith('__') or '[' in taxon:
                break

            node = get_node(cur_node, taxon)

            if node is None:
                raise ValueError("%s doesn't exist!" % taxon)

            cur_node = node

        if cur_node['popcount'] == 1:
            unique.append(tax_string)
        elif cur_node['popcount'] / popsize <= rare_threshold:
            rare.append(tax_string)

    return (rare, unique)


def traverse(node):
    """Post-order traversal of the full tree"""
    for c in node['children']:
        for gc in traverse(c):
            yield gc
    yield node


def build_tree_from_taxontable(table):
    """Construct a tree from a taxon table

    returns (tree, sample_taxa_lookup) where sample_taxa_lookup is:
    {sample_id:[['k__foo','p__bar',...]]}
    """
    sample_taxa_lookup = {}
    for taxa_freqs, sample_id, sample_md in table.iterSamples():
        sample_taxa = []
        for taxon, freq in zip(table.ObservationIds, taxa_freqs):
            if freq > 0:
                sample_taxa.append(taxon.split('; '))
        sample_taxa_lookup[sample_id] = sample_taxa

    tree = update_tree(None, sample_taxa_lookup.values())

    return tree, sample_taxa_lookup


def sample_rare_unique(tree, table, all_sample_taxa, rare_threshold):
    """Get the rare and unique taxa per sample

    returns (sample_id, biom table w/o rare and uniques, rare, unique)
    """
    def make_filter_f(r, u):
        r = ['; '.join(i) for i in r]
        u = ['; '.join(i) for i in u]
        to_remove = set(r).union(set(u))

        def f(v, i, md):
            if i in to_remove:
                return False
            return True
        return f

    for sample_id, sample_taxa in all_sample_taxa.iteritems():
        rare, unique = get_rare_unique(tree, sample_taxa, rare_threshold)
        filter_f = make_filter_f(rare, unique)

        if table is None:
            yield (sample_id, None, rare, unique)
        else:
            filtered = table.filterObservations(filter_f)
            yield (sample_id, filtered, rare, unique)

if __name__ == '__main__':
    from sys import argv
    table = parse_biom_table(open(argv[1]))
    tree, all_taxa = build_tree_from_taxontable(table)
    print "#SampleID\ttype\ttaxon"
    for samp, ft, rare, uniq in sample_rare_unique(tree, None, all_taxa, 0.01):
        for r in rare:
            print "%s\t%s\t%s" % (samp, 'rare', r)
        for u in uniq:
            print "%s\t%s\t%s" % (samp, 'unique', u)
