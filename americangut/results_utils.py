#!/usr/bin/env python

import re
import os
import shutil
import zipfile
from functools import partial
from americangut.util import check_file
from biom.parse import parse_biom_table
from itertools import izip
from collections import defaultdict


# These are the data files in the American-Gut repository that are used for
# results processing
_data_files = [
        # (directory, filename)
        ('AG', 'AG_100nt.biom.gz'),
        ('AG', 'AG_100nt.txt'),
        ('PGP', 'PGP_100nt.biom.gz'),
        ('PGP', 'PGP_100nt.txt'),
        ('HMP', 'HMPv35_100nt.biom.gz'),
        ('HMP', 'HMPv35_100nt.txt'),
        ('GG', 'GG_100nt.biom.gz'),
        ('GG', 'GG_100nt.txt')
        ]


# These are the Latex templates for the different results types
_templates = {
        'fecal': ('template_gut.tex', 'macros_gut.tex'),
        'oralskin': ('template_oralskin.tex', 'macros_oralskin.tex')
        }


def get_path(d, f):
    """Check and get a path, or throw IOError"""
    path = os.path.join(d, f)
    check_file(path)
    return path


def get_repository_dir():
    """Get the root of the American-Gut repository"""
    expected = os.path.abspath(__file__).rsplit('/', 2)[0]

    # get_path verifies the existance of these directories
    _ = get_path(expected, 'data')
    _ = get_path(expected, 'latex')

    return expected


def get_repository_data():
    """Get the path to the data"""
    return get_path(get_repository_dir(), 'data')


def get_repository_latex():
    """Get the path to the latex directory"""
    return get_path(get_repository_dir(), 'latex')


def get_repository_latex_pdfs(sample_type):
    """Get the Latex static PDFs directory"""
    latex_dir = get_repository_latex()

    if sample_type == 'oralskin':
        pdfs_dir = get_path(latex_dir, 'pdfs-oralskin')
    elif sample_type == 'fecal':
        pdfs_dir = get_path(latex_dir, 'pdfs-gut')
    else:
        raise ValueError("Unknown sample type: %s" % sample_type)

    check_file(pdfs_dir)

    return pdfs_dir


def _stage_static_latex(sample_type, working_dir):
    latex_dir = get_repository_latex()

    for item in _templates[sample_type]:
        src = get_path(latex_dir, item)
        shutil.copy(src, working_dir)


def _stage_static_pdfs(sample_type, working_dir):
    pdfs_dir = get_repository_latex_pdfs(sample_type)

    for f in os.listdir(pdfs_dir):
        src = get_path(pdfs_dir, f)
        shutil.copy(src, working_dir)


def _stage_static_data(working_dir, debug):
    data_dir = get_repository_data()

    for d, f in _data_files:
        if debug:
            d = d + '_debug'
        src = get_path(get_path(data_dir, d), f)
        shutil.copy(src, working_dir)


def stage_static_files(sample_type, working_dir, debug=False):
    """Stage static files in the current working directory"""
    _stage_static_data(working_dir, debug)
    _stage_static_latex(sample_type, working_dir)
    _stage_static_pdfs(sample_type, working_dir)


# use participant names only if the data are available.
# NOTE: these data are _not_ part of the github repository for
#       privacy reasons.
def parse_identifying_data(path, passwd, embedded_file='participants.txt'):
    """Process identifying data if available

    The expected format of the file is a passworded zipfile that contains
    an embedded, tab delimited file. The format of the tab delimited file
    is expected to be barcode TAB participant name
    """
    if path is not None:
        zf = zipfile.ZipFile(path)
        zf.setpassword(passwd)

        participants = {}
        for l in zf.read(embedded_file).splitlines():
            if l.startswith('#'):
                continue

            bc, name = l.strip().split('\t')[:2]
            participants[bc] = name.replace(",", "")

        print "Using identified data!"
    else:
        participants = None

    return participants


def parse_previously_printed(path):
    """Returns the set of previously printed barcodes

    The format of the file to be parsed is a single column of sample barcodes
    """
    if path is not None:
        prev_printed = set([l.strip() for l in open(path)])
    else:
        prev_printed = set([])
    return prev_printed


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


def massage_mapping(in_fp, out_fp, body_site_column_name, exp_acronym):
    """Simplify the mapping file for use in figures

    in_fp : input file-like object
    out_fp : output file-like object
    body_site_column_name : specify the column name for body
    exp_acronym : short name for the study

    Returns False on failure, True on success
    """
    def err_msg(issue, id_):
        print "SampleID: %s, %s" % (id_, issue)


    age_cat_map = [(0,2,'Baby'),
                   (2,13,'Child'),
                   (13,20,'Teen'),
                   (20,30,'20s'),
                   (30,40,'30s'),
                   (40,50,'40s'),
                   (50,60,'50s'),
                   (60,70,'60s'),
                   (70,80,'70s'),
                   (80,99999,'Older than 80')]
    bmi_cat_map = [(0, 18.5,'Underweight'),
                   (18.5, 25,'Normal'),
                   (25, 30,'Overweight'),
                   (30, 35,'Moderately obese'),
                   (35, 40,'Severely obese'),
                   (40, 99999,'Very severely obese')]

    mapping_lines = [l.strip().split('\t') for l in in_fp]

    header = mapping_lines[0]
    header_low = [x.lower() for x in header]

    bodysite_idx = header_low.index(body_site_column_name.lower())
    country_idx = header_low.index('country')

    try:
        age_idx = header_low.index('age')
    except ValueError:
        age_idx = None

    try:
        bmi_idx = header_low.index('bmi')
    except ValueError:
        bmi_idx = None

    new_mapping_lines = [header[:]]
    new_mapping_lines[0].append('SIMPLE_BODY_SITE')
    new_mapping_lines[0].append('TITLE_ACRONYM')
    new_mapping_lines[0].append('TITLE_BODY_SITE')
    new_mapping_lines[0].append('HMP_SITE')

    if age_idx is not None:
        new_mapping_lines[0].append('AGE_CATEGORY')
    if bmi_idx is not None:
        new_mapping_lines[0].append('BMI_CATEGORY')

    for l in mapping_lines[1:]:
        new_line = l[:]
        body_site = new_line[bodysite_idx]
        country = new_line[country_idx]

        # grab the body site
        if body_site.startswith('UBERON_'):
            body_site = body_site.split('_',1)[-1].replace("_"," ")
        elif body_site.startswith('UBERON:'):
            body_site = body_site.split(':',1)[-1]
        elif body_site in ['NA', 'unknown']:
            # controls, environmental, etc
            continue
        else:
            err_msg("Unknown body site: %s" % body_site, new_line[0])
            continue

        # remap the body site
        if body_site.lower() not in simple_matter_map:
            err_msg("Could not remap: %s" % body_site, new_line[0])
            continue
        else:
            body_site = simple_matter_map[body_site.lower()]

        if exp_acronym == 'HMP':
            hmp_site = 'HMP-%s' % body_site
        else:
            hmp_site = body_site

        # simplify the country
        if country.startswith('GAZ:'):
            country = country.split(':',1)[-1]
        else:
            err_msg("Could not parse country %s" % country, new_line[0])
            continue

        if age_idx is not None:
            age_cat = None
            if new_line[age_idx] in ['NA','None']:
                age_cat = 'Unknown'
            else:
                try:
                    # PGP is currently in age ranges, ignoring those for now
                    age = float(new_line[age_idx])
                except ValueError:
                    age_cat = 'Unknown'

            if age_cat is not 'Unknown':
                for low,high,cat in age_cat_map:
                    if low <= age < high:
                        age_cat = cat
                        break
                if age_cat is None:
                    err_msg("Unknown age: %f", new_line[0])
                    continue

        if bmi_idx is not None:
            if new_line[bmi_idx] in ['NA','', 'None']:
                bmi_cat = 'Unknown'
            else:
                bmi = float(new_line[bmi_idx])
                bmi_cat = None
                for low,high,cat in bmi_cat_map:
                    if low <= bmi < high:
                        bmi_cat = cat
                        break
                if bmi_cat is None:
                    err_msg("Unknown BMI: %f" % bmi, new_line[0])

        new_line.append(body_site)
        new_line.append(exp_acronym)
        new_line.append("%s-%s" % (exp_acronym, body_site))
        new_line[country_idx] = country
        new_line.append(hmp_site)

        if age_idx is not None:
            new_line.append(age_cat)

        if bmi_idx is not None:
            new_line.append(bmi_cat)

        new_mapping_lines.append(new_line)

    out_fp.write('\n'.join(['\t'.join(l) for l in new_mapping_lines]))
    out_fp.write('\n')


def filter_mapping_file(in_fp, out_fp, columns_to_keep):
    """Filter out columns in a mapping file

    in_fp : the input file-like object
    out_fp : the output file-like object
    columns_to_keep : a dict of the columns to keep, valued by specific category
        value if desired to filter out samples that don't meet a given
        criteria. In other words, a row is retained if the function associated
        with the key "foo" returns True, or the row is retained if the value
        associated with "foo" is None.
    """
    lines = [l.strip().split('\t') for l in in_fp]
    header = lines[0][:]
    header_lower = [x.lower() for x in header]

    # ensure SampleID is always first
    new_header = ["#SampleID"]
    indices = [0] # always keep SampleID
    for c in columns_to_keep:
        if c.lower() not in header_lower:
            raise ValueError("Cannot find %s!" % c)

        indices.append(header_lower.index(c.lower()))
        new_header.append(c)
    columns_to_keep['#SampleID'] = None # add for consistency

    new_lines = [new_header]
    for l in lines[1:]:
        new_line = []

        keep = True
        # fetch values from specific columns
        for column, index in zip(new_header, indices):
            value = l[index]
            if columns_to_keep[column] is None:
                new_line.append(value)
            elif not columns_to_keep[column](value):
                keep = False
                break
            else:
                new_line.append(value)

        if keep:
            new_lines.append(new_line)

    out_fp.write('\n'.join(['\t'.join(l) for l in new_lines]))
    out_fp.write('\n')


def construct_svg_smash_commands(files, ids, cmd_format, cmd_args):
    """Format the SVG smashing commands

    files : list of files
    ids : set of ids
    cmd_format : a string to format
    cmd_args : a dict of strings that can get filled into cmd_format
    """
    commands = []
    for f in files:
        if not f.startswith('Figure'):
            continue

        prefix, remainder = f.split('.', 1)

        try:
            id_, remainder = remainder.rsplit('_', 1)
        except:
            # GLOBAL SVG for each figure
            assert remainder == 'GLOBAL'
            continue

        # ignore svgs for non-AG points
        if id_ not in ids:
            continue

        args = cmd_args.copy()
        args['sample_id'] = id_
        args['prefix'] = prefix
        commands.append(cmd_format % args)
    return commands


def chunk_list(items, chunk_size=25):
    """Chunk up a list of items"""
    start = 0
    for end in range(chunk_size, len(items) + chunk_size, chunk_size):
        chunk = items[start:end]
        start = end
        yield chunk


def construct_phyla_plots_cmds(sample_ids, cmd_format, cmd_args):
    """Constuct the phlya plots commands"""
    commands = []
    for chunk in chunk_list(sample_ids):
        args = cmd_args.copy()
        args['samples'] = ','.join(chunk)
        commands.append(cmd_format % args)
    return commands


def per_sample_taxa_summaries(open_table, output_format):
    """Write out per-sample taxonomy summaries

    open_table : an open file-like biom table
    output_format : a path that supports a string format, eg:
        foo/bar_%s.txt
    """
    t = parse_biom_table(open_table)
    header = "#taxon\trelative_abundance\n"

    for v, id_, md in t.iterSamples():
        with open(output_format % id_, 'w') as f:
            f.write(header)

            for sorted_v, taxa in sorted(zip(v, t.ObservationIds))[::-1]:
                if sorted_v:
                    f.write("%s\t%f\n" % (taxa, sorted_v))


class MissingFigure(Exception):
    pass


def bootstrap_result(rel_existing_path, static_paths, base_cmd_fmt,
                     to_pdf_fmt, sample_id, name):
    """Stage for results

    sample_id : an id
    name : None or str
    rel_existing_path : a function that gets an existing path
    static_paths : a dict of paths
    base_cmd_fmt : base format for the commands to execute
    to_pdf_fmt : base format for the call to construct the latex PDF
    """
    if name is None:
        unidentified = rel_existing_path('unidentified')
        bootstrap_path = lambda x: os.path.join(unidentified, x)
    else:
        identified = rel_existing_path('identified')
        bootstrap_path = lambda x: os.path.join(identified, x)

    template_path = rel_existing_path('template_files')
    indiv_dir = bootstrap_path(sample_id)
    pdf_dir = os.path.join(indiv_dir, 'pdfs-gut')
    tex_path = lambda x: os.path.join(indiv_dir, x)
    fig_pdf_path = lambda x: os.path.join(pdf_dir, x)
    template_files_path = lambda x: os.path.join(template_path, x)

    fig1_src = template_files_path("Figure_1.%s_huge.pdf" % sample_id)
    fig2_src = template_files_path("Figure_2.%s_huge.pdf" % sample_id)
    fig3_src = template_files_path("Figure_3.%s_huge.pdf" % sample_id)
    fig4_src = template_files_path("Figure_4_%s.pdf" % sample_id)
    fig6_src = template_files_path("Figure_6_%s.txt" % sample_id)
    macros_src = template_files_path("macros_%s.tex" % sample_id)

    fig1_dst = fig_pdf_path("figure1.pdf")
    fig2_dst = fig_pdf_path("figure2.pdf")
    fig3_dst = fig_pdf_path("figure3.pdf")
    fig4_dst = fig_pdf_path("figure4.pdf")
    fig6_dst = tex_path("%s_taxa.txt" % sample_id)
    macros_dst = tex_path("macros_gut.tex")
    template_dst = tex_path('%s.tex' % sample_id)

    check_file(fig1_src, e=MissingFigure)
    check_file(fig2_src, e=MissingFigure)
    check_file(fig3_src, e=MissingFigure)
    check_file(fig4_src, e=MissingFigure)
    check_file(fig6_src, e=MissingFigure)
    check_file(macros_src, e=MissingFigure)

    cmds = []
    cmds.append('mkdir -p %s' % pdf_dir)
    cmds.append('cp %s %s' % (fig1_src, fig1_dst))
    cmds.append('cp %s %s' % (fig2_src, fig2_dst))
    cmds.append('cp %s %s' % (fig3_src, fig3_dst))
    cmds.append('cp %s %s' % (fig4_src, fig4_dst))
    cmds.append('cp %s %s' % (fig6_src, fig6_dst))
    cmds.append('cp %s %s' % (macros_src, macros_dst))
    cmds.append('cp %s %s' % (static_paths['template'], template_dst))
    cmds.append('cp %s %s' % (static_paths['aglogo'], pdf_dir))
    cmds.append('cp %s %s' % (static_paths['fig1_legend'], pdf_dir))
    cmds.append('cp %s %s' % (static_paths['fig2_legend'], pdf_dir))
    cmds.append('cp %s %s' % (static_paths['fig2_2ndlegend'], pdf_dir))
    cmds.append('cp %s %s' % (static_paths['fig3_legend'], pdf_dir))
    cmds.append('cp %s %s' % (static_paths['fig4_overlay'], pdf_dir))
    cmds.append('cp %s %s' % (static_paths['fig1_ovals'], pdf_dir))
    cmds.append('cp %s %s' % (static_paths['fig2_ovals'], pdf_dir))
    cmds.append('cp %s %s' % (static_paths['ball_legend'], pdf_dir))
    cmds.append('cp %s %s' % (static_paths['title'], pdf_dir))

    name_fmt = "echo '\n\def\yourname{%s}\n' >> %s"
    if name is None:
        cmds.append(name_fmt % ("unidentified", macros_dst))
    else:
        cmds.append(name_fmt % (name, macros_dst))

    indiv_cmd = base_cmd_fmt % (static_paths['working_dir'], '; '.join(cmds))
    latex_cmd = to_pdf_fmt % {'path': indiv_dir, 'input': template_dst}

    return (indiv_cmd, latex_cmd)


def construct_bootstrap_and_latex_commands(ids, participants, rel_existing_path,
                                           static_paths, base_cmd_fmt,
                                           to_pdf_fmt):
    """Construct the commands to bootstrap results and latex generation

    ids : an iterable of ids
    participants : None or a dict mapping barcodes to participant names
    rel_existing_path : a function that gets an existing path
    static_paths : a dict of paths
    base_cmd_fmt : base format for the commands to execute
    to_pdf_fmt : base format for the call to construct the latex PDF
    """
    bs_f = partial(bootstrap_result, rel_existing_path, static_paths,
                   base_cmd_fmt, to_pdf_fmt)
    indiv_cmds = []
    latex_cmds = []
    missing = []
    for sample_id in ids:
        name = None
        if participants is not None:
            bc = sample_id.split('.')[0]
            if bc in participants:
                name = participants[bc]

        # unidentified
        try:
            indiv_cmd, latex_cmd = bs_f(sample_id, None)
        except MissingFigure:
            missing.append(sample_id)
            continue

        indiv_cmds.append(indiv_cmd)
        latex_cmds.append(latex_cmd)

        # identified
        if name:
            indiv_cmd, latex_cmd = bs_f(sample_id, name)
            indiv_cmds.append(indiv_cmd)
            latex_cmds.append(latex_cmd)

    return (indiv_cmds, latex_cmds, missing)


def harvest(path):
    """harvest PDFs"""
    harvest_path = os.path.join(path, 'harvested')
    if not os.path.exists(harvest_path):
        os.mkdir(harvest_path)

    for dirpath, dirnames, filenames in os.walk(path):
        try:
            dirnames.remove('harvested')
        except ValueError:
            pass

        try:
            dirnames.remove('pdfs')
        except ValueError:
            pass

        sample_suffix = re.search('\d+\.\d+', dirpath)
        if sample_suffix is None:
            continue
        else:
            sample = sample_suffix.group().split('.')[0]

        pdf = os.path.join(path, dirpath, "%s.pdf" % sample)
        if os.path.exists(pdf):
            yield "mv %s %s" % (pdf, harvest_path)


def pdf_smash(path, tag, pdf_smash_fmt, n_per_result=30,
              previously_printed=None):
    """Combine sets of PDFs into single documents

    path : a path to where the PDFs are
    tag : some tag to put on the file names
    pdf_smash_fmt : command format to use
    n_per_result : number of PDFs to smash together
    previously_printed : set of previously printed barcodes or None
    """
    if previously_printed is None:
        previously_printed = set([])

    files = []
    for f in os.listdir(path):
        bc, extension = os.path.splitext(f)
        if extension != '.pdf':
            continue
        if bc in previously_printed:
            continue
        files.append(os.path.join(path, f))

    result_path = os.path.join(path, 'pdf_smash')
    if not os.path.exists(result_path):
        os.mkdir(result_path)

    # sort and then filter out. filtering is by barcode without prefix
    sort_key = lambda x: int(x.rsplit('/')[-1].split('.')[0])
    files_ordered = sorted(files, key=sort_key)

    smash_set = []
    barcode_set = []
    bc_f = lambda ch: '\n'.join([f.rsplit('/')[-1].split('.')[0] for f in ch])

    for chunk in chunk_list(files_ordered, n_per_result):
        smash_set.append(' '.join(chunk))
        barcode_set.append(bc_f(chunk))

    smash_cmds = []
    smash_basename = os.path.join(result_path, "%s_smashset_%d")
    for set_number, (pdfs, barcodes) in enumerate(zip(smash_set, barcode_set)):
        filename_base = smash_basename % (tag, set_number)
        filename_pdf = filename_base + '.pdf'
        filename_txt = filename_base + '.txt'

        with open(filename_txt, 'w') as f:
            f.write(barcodes)
            f.write('\n')

        smash_cmds.append(pdf_smash_fmt % {'output': filename_pdf,
                                           'pdfs': pdfs})

    ordered_barcodes_path = os.path.join(result_path, 'ordered_barcodes.txt')

    with open(ordered_barcodes_path, 'w') as ordered_barcodes:
        ordered_barcodes.write('\n'.join(barcode_set))
        ordered_barcodes.write('\n')

    return smash_cmds


def count_unique_sequences_per_otu(otu_ids, otu_map_file, input_seqs_file):
    """Counts unique sequences per-OTU for a given set of OTUs

    otu_ids: a set of OTU IDs
    otu_map_file: file-like object in the format of an OTU map
    input_seqs_file: FASTA containing sequences that were used to generate
                     the otu_map_file

    Returns a nested dict structure: {otu_id: {sequence: count}}
    """
    # This will hold the OTU map for the OTUs in otu_ids
    otu_map = {x: set() for x in otu_ids}

    # go through the otu map and save the lines of interest to the otu_map
    # data structure above
    print "Reading OTU map..."
    for line in otu_map_file:
        otu_id, seq_ids = line.strip().split('\t', 1)
        if otu_id in otu_ids:
            otu_map[otu_id] = set(seq_ids.split('\t'))

    # this will hold, for each OTU in otus, counts of each unique sequence
    # observed in that OTU
    unique_counts = {x: defaultdict(int) for x in otu_ids}

    # go through input fasta file TWO LINES AT A TIME, counting unique
    # sequences in each OTU of intrest
    print "Reading FASTA file and counting unique sequences..."
    for header, sequence in izip(input_seqs_file, input_seqs_file):
        header = header.strip()
        sequence = sequence.strip()
        seq_id = header.split(' ', 1)[0][1:]
        for otu_id in otu_ids:
            if seq_id in otu_map[otu_id]:
                unique_counts[otu_id][sequence] += 1
                break
    print "Done."

    return unique_counts


def write_contaminant_fasta(unique_counts, output_file, abundance_threshold):
    """Writes a FASTA file of sequences determined to be contaminants

    If one unique sequences composes more than abundance_threshold of the OTU,
    that sequences is marked as a contaminant and written to output_file.

    unique_counts: a nested dict of the form {otu_id: {sequence: count}}
                   E.g., the output of count_unique_sequences_per_otu
    output_file: a file-like object ready for writing
    abundance_threshold: If a sequence composes more than this percent of an
                         OTU, then it is marked as a contaminant
    """
    for otu_id, otu_counts in unique_counts.iteritems():
        otu_total_count = sum([count for seq, count in otu_counts.iteritems()])

        counter = 0
        for seq, count in sorted(otu_counts.items(), key=lambda x:x[1],
                reverse=True):
            counter += 1
            if 1.0*count/otu_total_count > abundance_threshold:
                output_file.write('>%s_%d\n%s\n' % (otu_id, counter, seq))
