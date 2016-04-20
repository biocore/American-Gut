#!/usr/bin/env python

import os
import shutil
import zipfile
from itertools import izip
from collections import defaultdict
from functools import partial

import matplotlib.pyplot as plt
import seaborn as sn
from biom.parse import parse_biom_table

from americangut.util import check_file

# These are the data files in the American-Gut repository that are used for
# results processing
_data_files = [
    # (directory, filename)
    ('PGP', 'PGP_100nt.biom.gz'),
    ('PGP', 'PGP_100nt.txt'),
    ('HMP', 'HMPv35_100nt.biom.gz'),
    ('HMP', 'HMPv35_100nt.txt'),
    ('GG', 'GG_100nt.biom.gz'),
    ('GG', 'GG_100nt.txt'),
    ('AG', 'pgp_agp_barcodes.txt'),
    ('AG', 'BLOOM.fasta')
]


# These are the Latex templates for the different results types
_templates = {
    'fecal': 'template_gut.tex',
    'oral': 'template_oralskin.tex',
    'skin': 'template_oralskin.tex'
}


_mod1_bits = ['metadata_charts.json',
              'mod1_main.tex',
              'stacked_plots_key_taxa.txt',
              'fecal_identified.txt',
              'skin_identified.txt',
              'oral_identified.txt']


def get_path(d, f):
    """Check and get a path, or throw IOError"""
    path = os.path.join(d, f)
    check_file(path)
    return path


def get_repository_dir():
    """Get the root of the American-Gut repository"""
    expected = os.path.abspath(__file__).rsplit('/', 2)[0]

    # get_path verifies the existance of these directories
    get_path(expected, 'data')
    get_path(expected, 'latex')

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
    sample_type = sample_type.lower()

    if sample_type == 'skin':
        pdfs_dir = get_path(latex_dir, 'pdfs-oralskin')
    elif sample_type == 'oral':
        pdfs_dir = get_path(latex_dir, 'pdfs-oralskin')
    elif sample_type == 'fecal':
        pdfs_dir = get_path(latex_dir, 'pdfs-gut')
    elif sample_type == 'mod1':
        pdfs_dir = get_path(latex_dir, 'pdfs-mod1')
    else:
        raise ValueError("Unknown sample type: %s" % sample_type)

    check_file(pdfs_dir)

    return pdfs_dir


def stage_static_latex(sample_type, working_dir):
    latex_dir = get_repository_latex()

    item = _templates[sample_type.lower()]
    src = get_path(latex_dir, item)
    shutil.copy(src, working_dir)


def stage_static_pdfs(sample_type, working_dir):
    pdfs_dir = get_repository_latex_pdfs(sample_type)

    for f in os.listdir(pdfs_dir):
        if f.endswith('.pdf'):
            src = get_path(pdfs_dir, f)
            shutil.copy(src, working_dir)


def _stage_static_data(working_dir, debug):
    data_dir = get_repository_data()

    for d, f in _data_files:
        if debug:
            d = d + '_debug'
        src = get_path(get_path(data_dir, d), f)
        shutil.copy(src, working_dir)


def _stage_static_mod1(working_dir):
    latex_dir = get_repository_latex()
    statics = get_repository_latex_pdfs('mod1')

    for f in _mod1_bits:
        src = get_path(latex_dir, f)
        shutil.copy(src, working_dir)

    for f in os.listdir(statics):
        src = get_path(statics, f)
        shutil.copy(src, working_dir)


def stage_static_files(sample_type, working_dir, debug=False):
    """Stage static files in the current working directory"""
    stage_static_latex(sample_type, working_dir)
    stage_static_pdfs(sample_type, working_dir)


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


def filter_mapping_file(in_fp, out_fp, columns_to_keep):
    """Filter out columns in a mapping file

    in_fp : the input file-like object
    out_fp : the output file-like object
    columns_to_keep : a dict of the columns to keep, valued by specific
        category value if desired to filter out samples that don't meet a given
        criteria. In other words, a row is retained if the function associated
        with the key "foo" returns True, or the row is retained if the value
        associated with "foo" is None.
    """
    lines = [l.strip().split('\t') for l in in_fp]
    header = lines[0][:]
    header_lower = [x.lower() for x in header]

    # ensure SampleID is always first
    new_header = ["#SampleID"]
    indices = [0]  # always keep SampleID
    for c in columns_to_keep:
        if c.lower() not in header_lower:
            raise ValueError("Cannot find %s!" % c)

        indices.append(header_lower.index(c.lower()))
        new_header.append(c)
    columns_to_keep['#SampleID'] = None  # add for consistency

    new_lines = [new_header]
    for l in lines[1:]:
        new_line = []

        keep = True
        # fetch values from specific columns
        for column, index in zip(new_header, indices):
            try:
                value = l[index]
            except:
                raise
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

    for v, id_, md in t.iter():
        with open(output_format % id_, 'w') as f:
            f.write(header)

            for sorted_v, taxa in \
                    sorted(zip(v, t.ids(axis='observation')))[::-1]:
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

        def bootstrap_path(x):
            return os.path.join(unidentified, x)
    else:
        identified = rel_existing_path('identified')

        def bootstrap_path(x):
            return os.path.join(identified, x)

    template_path = rel_existing_path('template_files')
    indiv_dir = bootstrap_path(sample_id)
    pdf_dir = os.path.join(indiv_dir, 'pdfs-gut')

    def tex_path(x):
        return os.path.join(indiv_dir, x)

    def fig_pdf_path(x):
        return os.path.join(pdf_dir, x)

    def template_files_path(x):
        return os.path.join(template_path, x)

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


def construct_bootstrap_and_latex_commands(ids, participants,
                                           rel_existing_path,
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
    """harvest PDFs and taxa summaries"""
    harvest_path = os.path.join(path, 'harvested')
    if not os.path.exists(harvest_path):
        os.mkdir(harvest_path)

    for name in os.listdir(path):
        if 'harvested' in name:
            continue
        if 'pdfs' in name:
            continue
        try:
            float(name)  # should work for 000001111.123123 and 000001111
        except ValueError:
            continue

        dst_name = name.split('.')[0] if '.' in name else name

        src = os.path.join(path, name, "%s.pdf" % name)
        dst = os.path.join(harvest_path, "%s.pdf" % dst_name)

        if os.path.exists(src):
            yield "mv %s %s" % (src, dst)

        src = os.path.join(path, name, "%s_taxa.txt" % name)
        dst = os.path.join(harvest_path, "%s.txt" % dst_name)

        if os.path.exists(src):
            yield "mv %s %s" % (src, dst)


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
    def sort_key(x):
        return int(x.rsplit('/')[-1].split('.')[0])
    files_ordered = sorted(files, key=sort_key)

    smash_set = []
    barcode_set = []

    def bc_f(ch):
        return '\n'.join([f.rsplit('/')[-1].split('.')[0] for f in ch])

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

    return unique_counts


def write_bloom_fasta(unique_counts, output_file, abundance_threshold):
    """Writes a FASTA file of sequences determined to be the result of a bloom

    If one unique sequences composes more than abundance_threshold of the OTU,
    that sequences is marked as a bloom sequence and written to output_file.

    unique_counts: a nested dict of the form {otu_id: {sequence: count}}
                   E.g., the output of count_unique_sequences_per_otu
    output_file: a file-like object ready for writing
    abundance_threshold: If a sequence composes more than this percent of an
                         OTU, then it is marked as a bloom sequence
    """
    for otu_id, otu_counts in unique_counts.iteritems():
        otu_total_count = sum([count for seq, count in otu_counts.iteritems()])

        counter = 0
        for seq, count in sorted(otu_counts.items(), key=lambda x: x[1],
                                 reverse=True):
            counter += 1
            if 1.0*count/otu_total_count > abundance_threshold:
                output_file.write('>%s_%d\n%s\n' % (otu_id, counter, seq))


def plot_alpha(sample, alpha_map, alpha_field, group_field='SIMPLE_BODY_SITE',
               output_dir=None, xlabel=None, fp=None, sample_color='#525252',
               highlight_range=None, categorical=False, debug=False):
    """Generates a distrbution plot for the data

    Parameters
    ----------
    sample : str
        The sample ID to be plotted
    alpha_map_fp : pandas DataFrame
        A pandas dataframe containing the sample metadata. The sample ID
        should be given in the `'#SampleID'` column, a column with
        the name given by `alpha_field` contains alpha diversity values,
        and the `group_field` column specifying the groups which should be
        used to seperate the data for making the distribution plot.
    alpha_field : str
        The name of the column in `alpha_map` which includes the alpha
        diversity values.
    group_field : str
        Default is 'SIMPLE_BODY_SITE'. The name of the column in `alpha_map`
        which provides the grouping for generating distribution plots.
    output_dir : str
        The location where the alpha diversity figures should be saved.
    xlabel : str
        Text describing the quantity on the x-axis.
    sample_color : str, optional
        Hex color for the sample line. Default '#525252'
    highlight_range: list of [float, float], optional
        If given, the range to highlight under the curve. Will be lighter shade
        of sample_color. Default None.
    categorical : bool, optional
        Whether the text should be added categorically or per sample.
        Default False (add per sample)

    Returns
    -------
    If the sample is present, a matplotlib figure with the alpha diversity
    distribution and a line indicating the sample value is returned. If a
    file path is specified, the figure will be saved at the filepath instead
    of returning.

    If debug is passed, the following parameters are returned:
        group : str
            The value of the `group_field` for the sample
        group_alpha : ndarray
            The alpha diversity values associated with the group
        sample_alpha : float
            The alpha diversity for the sample
        xlabel : str
            The label used for the x-axis of the plot.

    """

    # Explicitly casts the alpha diversity to a float
    alpha_map[alpha_field] = alpha_map[alpha_field].astype(float)

    # Draws the observations and group
    group = alpha_map.loc[sample, group_field]
    sample_alpha = alpha_map.loc[sample, alpha_field]
    group_alpha = alpha_map.loc[alpha_map[group_field] == group, alpha_field]
    if categorical:
        # Keep only barcoded items, which can be converted to float.
        # Categories will be strings in the index
        def test(val):
            try:
                float(val)
                return True
            except ValueError:
                return False
        group_alpha = group_alpha.loc[[test(v) for v in group_alpha.index]]

    if xlabel is None:
        xlabel = '%sdiversity' % alpha_field.split('1')[0].replace('_', ' ')

    if debug:
        return group, group_alpha, sample_alpha, xlabel

    # Defines the group color. This is currently hardcoded, although the
    # longer term plan is to substitute in function which will define the color
    # based on the relationship between the sample and a yet to be written
    # predicted value.
    group_color = '#1f78b4'

    with sn.axes_style('ticks', {'axes.facecolor': 'none'}):
        # Sets up the axis for plotting
        ax = plt.axes()

        # Plots the distribution
        sn.kdeplot(group_alpha,
                   ax=ax,
                   legend=False,
                   color=group_color)
        ylim = ax.get_ylim()

        # Plots the individual line
        ax.plot([sample_alpha, sample_alpha], [-1, 1], color=sample_color)

        # Returns the y-limits to the original value and removes the ticks
        ax.set_ylim(ylim)
        ax.set_yticks([])
        # Removes the spine
        sn.despine(offset=5, trim=True, top=True, left=True, right=True)
        # Updates the xticks to match the correct font
        ax.set_xticklabels(map(int, ax.get_xticks()), size=11)
        ax.set_xlabel(xlabel, size=13)

        # Highlight range, if given
        if highlight_range is not None:
            xticks = ax.get_xticks()
            x, y = ax.lines[0].get_xydata().T
            m = (highlight_range[0] < x) & (x < highlight_range[1])
            ax.fill_between(x[m], 0, y[m], alpha=.15, color=group_color)
            # Fill between changes axis ever so slightly, so need to reset
            # back to what we want.
            ax.set_xticks(xticks)

        # Adds text describing the sample or category
        if categorical:
            text = '\n\n%s:\t%1.1f' % (sample, sample_alpha)
        else:
            text = 'Your Sample:\t%1.1f\nAverage:\t%1.1f' % (
                sample_alpha, group_alpha.mean())
        ax.text(x=ax.get_xticks().max(),
                y=ax.get_ylim()[1]*0.85,
                s=text,
                ha='right',
                size=11,
                )

    # Sets the figure size
    fig = ax.figure
    fig.set_size_inches((5, 2.5))
    ax.set_position((0.125, 0.375, 0.75, 0.5))

    if fp is None:
        return fig
    else:
        fig.savefig(fp, dpi=300)
        fig.clear()
