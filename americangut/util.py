#!/usr/bin/env python

import os
import urllib2
import gzip
import time
from itertools import izip
from StringIO import StringIO
from lxml import etree

from skbio.parse.sequences import parse_fastq


__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Daniel McDonald", "Adam Robbins-Pianka"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"


def pick_rarifaction_level(id_, lookups):
    """Determine which lookup has the appropriate key

    id_ is a barcode, e.g., '000001000'
    lookups is a list of tuples, e.g., [('10k',{'000001000':'000001000.123'})]

    The order of the lookups matters. The first lookup found with the key will
    be returned.

    None is returned if the key is not found
    """
    for name, lookup in lookups:
        if id_ in lookup:
            return name
    return None

def parse_mapping_file(open_file):
    """return (header, [(sample_id, all_other_fields)])

    """
    header = open_file.readline().strip()
    res = []

    for l in open_file:
        res.append(l.strip().split('\t',1))

    return (header, res)

def verify_subset(table, mapping):
    """Returns True/False if the table is a subset"""
    ids = set([i[0] for i in mapping])
    t_ids = set(table.SampleIds)

    return t_ids.issubset(ids)

def slice_mapping_file(table, mapping):
    """Returns a new mapping corresponding to just the ids in the table"""
    t_ids = set(table.SampleIds)
    res = []

    for id_, l in mapping:
        if id_ in t_ids:
            res.append('\t'.join([id_, l]))

    return res

def check_file(f):
    """Verify a file (or directory) exists"""
    if not os.path.exists(f):
        raise IOError("Cannot continue! The file %s does not exist!" % f)

def trim_fasta(input_fasta, output_fasta, length):
    """Trim FASTA sequences to a given length

    input_fasta: should be an open file. Every two lines should compose a
                 complete FASTA record (header, sequence)
    output_fasta: should be an open file ready for writing
    length: what length to trim the sequences to. Sequences shorter than
            length will not be modified.
    """
    # reads the FASTA file two lines at at a time
    # Assumptions: 1) each FASTA record is two lines
    #              2) There are no incomplete FASTA records
    for header, sequence in izip(input_fasta, input_fasta):
        header = header.strip()
        sequence = sequence.strip()[:length]
        output_fasta.write("%s\n%s\n" % (header, sequence))

def concatenate_files(input_files, output_file, read_chunk=10000):
    """Concatenate all input files and produce an output file

    input_fps is a list of open files
    output_fp is an open file ready for writing
    """
    for infile in input_files:
        chunk = infile.read(read_chunk)
        while chunk:
            output_file.write(chunk)
            chunk = infile.read(read_chunk)

def fetch_study_details(accession):
    """Fetch secondary accession and FASTQ details

    yields [(secondary_accession, fastq_url)]
    """
    url_fmt = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=%(accession)s&result=read_run&fields=secondary_sample_accession,submitted_ftp"
    res = fetch_url(url_fmt % {'accession': accession})

    return [tuple(l.strip().split('\t')) for l in res.readlines()[1:]]

def fetch_url(url):
    """Return an open file handle"""
    # really should use requests instead of urllib2
    attempts = 0
    res = None

    while attempts < 5:
        attempts += 1
        try:
            res = urllib2.urlopen(url)
        except urllib2.HTTPError as e:
            if e.code == 500:
                time.sleep(5)
                continue
            else:
                raise

    if res is None:
        raise ValueError("Failed at fetching %s" % url)

    return StringIO(res.read())

def fetch_seqs_fastq(url):
    """Fetch a FTP item"""
    # not using a url_fmt here as the directory structure has potential to
    # be different between studies
    if not url.startswith('ftp://'):
        url = "ftp://%s" % url

    res = fetch_url(url)

    return gzip.GzipFile(fileobj=res)

def fetch_metadata_xml(accession):
    """Fetch sample metadata"""
    url_fmt = "http://www.ebi.ac.uk/ena/data/view/%(accession)s&display=xml"
    res = fetch_url(url_fmt % {'accession': accession})

    root = etree.parse(res).getroot()
    sample = root.getchildren()[0]
    attributes = sample.find('SAMPLE_ATTRIBUTES')

    metadata = {}
    for node in attributes.iterfind('SAMPLE_ATTRIBUTE'):
        tag, value = node.getchildren()
        metadata[tag.text.strip('" ')] = value.text.strip('" ')
    return metadata

def fetch_study(accession, metadata_path, fasta_path):
    """Fetch and dump a full study

    Grab and dump a full study
    """
    all_md = {}
    all_cols = set([])
    md_f = open(metadata_path, 'w')
    fasta_path = open(fasta_path, 'w')
    for sample, fastq_url in fetch_study_details(accession):
        # in the form seqs_000007123.1075697.fastq.gz
        # and unfortunately, the suffix (1075697) is missing and parts of the
        # current results processing depend on the suffix.
        fastq_filename = fastq_url.rsplit('/')[-1]
        qiimedb_samplename = fastq_filename.split('_')[-1].rsplit('.', 2)[0]

        md = fetch_metadata_xml(sample)
        all_md[qiimedb_samplename] = md
        all_cols.update(md)

        # write out fasta
        for id_, seq, qual in parse_fastq(fetch_seqs_fastq(fastq_url)):
            fasta_path.write(">%s\n%s\n" % (id_, seq))

    header = list(all_cols)
    md_f.write('#SampleID\t')
    md_f.write('\t'.join(header))
    md_f.write('\n')
    for sampleid, values in all_md.iteritems():
        to_write = [values.get(k, "no_data") for k in header]
        to_write.insert(0, sampleid)
        md_f.write('\t'.join(to_write))
        md_f.write('\n')

    md_f.close()
    fasta_path.close()

