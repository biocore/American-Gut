#!/usr/bin/env python

import os
import urllib2
import gzip
import time
from itertools import izip
from StringIO import StringIO
from collections import defaultdict

from lxml import etree
from skbio.parse.sequences import parse_fastq, parse_fasta

import americangut as ag


__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Daniel McDonald", "Adam Robbins-Pianka"]
__license__ = "BSD"
__version__ = "unversioned"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"


def get_path(path):
    """Get a relative path to the working directory

    Parameters
    ----------
    path : str
        A path

    Notes
    -----
    This method does not care if the path exists or not

    Returns
    -------
    str
        The filepath
    """
    return os.path.join(ag.working_dir, path)


def get_new_path(path):
    """Get a new relative path to the working directory

    Parameters
    ----------
    path : str
        A path that does not exist

    Notes
    -----
    It is only assured that the path does not exist at the time of function
    evaluation.

    Raises
    ------
    IOError
        If the path exists

    Returns
    -------
    str
        The filepath
    """
    path = get_path(path)

    if os.path.exists(path):
        raise IOError('%s already exists.' % path)

    return path


def get_existing_path(path):
    """Get an existing relative path to the working directory

    Parameters
    ----------
    path : str
        A path that exists

    Notes
    -----
    It is only assured that the path exists at the time of function evaluation

    Raises
    ------
    IOError
        If the path does not exist

    Returns
    -------
    str
        The filepath
    """
    path = get_path(path)

    if not os.path.exists(path):
        raise IOError('%s does not exist.' % path)

    return path


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
        res.append(l.strip().split('\t', 1))

    return (header, res)


def verify_subset(table, mapping):
    """Returns True/False if the table is a subset"""
    ids = set([i[0] for i in mapping])
    t_ids = set(table.ids())

    return t_ids.issubset(ids)


def slice_mapping_file(table, mapping):
    """Returns a new mapping corresponding to just the ids in the table"""
    t_ids = set(table.ids())
    res = []

    for id_, l in mapping:
        if id_ in t_ids:
            res.append('\t'.join([id_, l]))

    return res


def check_file(f, e=IOError):
    """Verify a file (or directory) exists"""
    if not os.path.exists(f):
        raise e("Cannot continue! The file %s does not exist!" % f)


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
    url_fmt = "http://www.ebi.ac.uk/ena/data/warehouse/" \
              "filereport?accession=%(accession)s&result=read_run&" \
              "fields=secondary_sample_accession,submitted_ftp"
    res = fetch_url(url_fmt % {'accession': accession})

    for line in res.readlines()[1:]:
        if 'ERA371447' in line:
            # Corrupt sequence files were uploaded to EBI for one of the AG
            # rounds. Ignoring entries associated with this accession works
            # around the corruption
            continue

        parts = line.strip().split('\t')
        if len(parts) != 2:
            continue
        else:
            yield tuple(parts)


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
        if value.text is None:
            metadata[tag.text.strip('" ').upper()] = 'no_data'
        else:
            metadata[tag.text.strip('" ').upper()] = value.text.strip('" ')

    description = sample.find('DESCRIPTION')
    metadata['Description'] = description.text.strip('" ')

    return metadata


def fetch_study(accession, base_dir):
    """Fetch and dump a full study

    Grab and dump a full study
    """
    metadata_path = os.path.join(base_dir, '%s.txt' % accession)
    fasta_path = os.path.join(base_dir, '%s.fna' % accession)

    if os.path.exists(fasta_path) and os.path.exists(metadata_path):
        # it appears we already have the accession, so short circuit
        return

    all_md = {}
    all_cols = set(['BarcodeSequence', 'LinkerPrimerSequence'])
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
        try:
            for id_, seq, qual in parse_fastq(fetch_seqs_fastq(fastq_url)):
                fasta_path.write(">%s\n%s\n" % (id_, seq))
        except:
            continue

    header = list(all_cols)
    md_f.write('#SampleID\t')
    md_f.write('\t'.join(header))
    md_f.write('\n')
    for sampleid, values in all_md.iteritems():
        to_write = [values.get(k, "no_data").encode('utf-8') for k in header]
        to_write.insert(0, sampleid)
        md_f.write('\t'.join(to_write))
        md_f.write('\n')

    md_f.close()
    fasta_path.close()


def count_seqs(seqs_fp, subset=None):
    """Could the number of FASTA records"""
    if subset is None:
        return sum(1 for line in seqs_fp if line.startswith(">"))
    else:
        subset = set(subset)
        count = 0
        for id_, seq in parse_fasta(seqs_fp):
            parts = id_.split()

            # check if the ID is there, and handle the qiimedb suffix case
            if parts[0] in subset:
                count += 1
            elif parts[0].split('.')[0] in subset:
                count += 1
        return count


def count_unique_participants(metadata_fp, criteria=None):
    """Count the number of unique participants"""
    if criteria is None:
        criteria = {}

    header = {k: i for i, k in enumerate(
              metadata_fp.next().strip().split('\t'))}

    count = set()
    for line in metadata_fp:
        line = line.strip().split('\t')
        keep = True
        for crit, val in criteria.items():
            if line[header[crit]] != val:
                keep = False
        if keep:
            count.add(line[header['HOST_SUBJECT_ID']])

    return len(count)


def count_samples(metadata_fp, criteria=None):
    """Count the number of samples

    criteria : dict
        Header keys and values to restrict by
    """
    if criteria is None:
        criteria = {}

    header = {k: i for i, k in enumerate(
              metadata_fp.next().strip().split('\t'))}

    count = 0
    for line in metadata_fp:
        line = line.strip().split('\t')
        keep = True
        for crit, val in criteria.items():
            if line[header[crit]] != val:
                keep = False
        if keep:
            count += 1

    return count


simple_matter_map = {
    'feces': 'FECAL',
    'sebum': 'SKIN',
    'tongue': 'ORAL',
    'skin': 'SKIN',
    'mouth': 'ORAL',
    'gingiva': 'ORAL',
    'gingival epithelium': 'ORAL',
    'nares': 'SKIN',
    'skin of hand': 'SKIN',
    'hand': 'SKIN',
    'skin of head': 'SKIN',
    'hand skin': 'SKIN',
    'throat': 'ORAL',
    'auricular region zone of skin': 'SKIN',
    'mucosa of tongue': 'ORAL',
    'mucosa of vagina': 'SKIN',
    'palatine tonsil': 'ORAL',
    'hard palate': 'ORAL',
    'saliva': 'ORAL',
    'stool': 'FECAL',
    'vagina': 'SKIN',
    'fossa': 'SKIN',
    'buccal mucosa': 'ORAL',
    'vaginal fornix': 'SKIN',
    'hair follicle': 'SKIN',
    'nostril': 'SKIN'
}


def clean_and_reformat_mapping(in_fp, out_fp, body_site_column_name,
                               exp_acronym):
    """Simplify the mapping file for use in figures

    in_fp : input file-like object
    out_fp : output file-like object
    body_site_column_name : specify the column name for body
    exp_acronym : short name for the study

    Returns a dict containing a description of any unprocessed samples.
    """
    errors = defaultdict(list)

    mapping_lines = [l.strip('\n').split('\t') for l in in_fp]

    header = mapping_lines[0]
    header_low = [x.lower() for x in header]

    bodysite_idx = header_low.index(body_site_column_name.lower())
    country_idx = header_low.index('country')

    new_mapping_lines = [header[:]]
    new_mapping_lines[0].append('SIMPLE_BODY_SITE')
    new_mapping_lines[0].append('TITLE_ACRONYM')
    new_mapping_lines[0].append('TITLE_BODY_SITE')
    new_mapping_lines[0].append('HMP_SITE')

    for l in mapping_lines[1:]:
        new_line = l[:]
        sample_id = new_line[0]
        body_site = new_line[bodysite_idx]
        country = new_line[country_idx]

        # grab the body site
        if body_site.startswith('UBERON_'):
            body_site = body_site.split('_', 1)[-1].replace("_", " ")
        elif body_site.startswith('UBERON:'):
            body_site = body_site.split(':', 1)[-1]
        elif body_site in ['NA', 'unknown', '', 'no_data']:
            errors[('unspecified_bodysite', body_site)].append(sample_id)
            continue
        else:
            raise ValueError("Cannot process: %s, %s" % (sample_id, body_site))

        # remap the body site
        if body_site.lower() not in simple_matter_map:
            errors[('unknown_bodysite', body_site)].append(sample_id)
            continue
        else:
            body_site = simple_matter_map[body_site.lower()]

        if exp_acronym == 'HMP':
            hmp_site = 'HMP-%s' % body_site
        else:
            hmp_site = body_site

        # simplify the country
        if country.startswith('GAZ:'):
            new_line[country_idx] = country.split(':', 1)[-1]
        else:
            errors[('unknown_country', country)].append(sample_id)
            continue

        new_line.append(body_site)
        new_line.append(exp_acronym)
        new_line.append("%s-%s" % (exp_acronym, body_site))
        new_line.append(hmp_site)

        new_mapping_lines.append(new_line)

    out_fp.write('\n'.join(['\t'.join(l) for l in new_mapping_lines]))
    out_fp.write('\n')

    return errors
