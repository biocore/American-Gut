#!/usr/bin/env python

import os
import shutil
import zipfile
from americangut.util import check_file


# These are the data files in the American-Gut repository that are used for
# results processing
_data_files = [
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


def _stage_static_data(working_dir):
    data_dir = get_repository_data()

    for d, f in _data_files:
        src = get_path(get_path(data_dir, d), f)
        shutil.copy(src, working_dir)


def stage_static_files(sample_type, working_dir):
    """Stage static files in the current working directory"""
    _stage_static_data(working_dir)
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
        'feces':'FECAL',
        'sebum':'SKIN',
        'tongue':'ORAL',
        'skin':'SKIN',
        'mouth':'ORAL',
        'gingiva':'ORAL',
        'gingival epithelium':'ORAL',
        'nares':'SKIN',
        'skin of hand':'SKIN',
        'hand':'SKIN',
        'skin of head':'SKIN',
        'hand skin':'SKIN',
        'throat':'ORAL',
        'auricular region zone of skin':'SKIN',
        'mucosa of tongue':'ORAL',
        'mucosa of vagina':'SKIN',
        'palatine tonsil':'ORAL',
        'hard palate':'ORAL',
        'saliva':'ORAL',
        'stool':'FECAL',
        'vagina':'SKIN',
        'fossa':'SKIN',
        'buccal mucosa':'ORAL',
        'vaginal fornix':'SKIN',
        'hair follicle':'SKIN',
        'nostril':'SKIN'
        }

def massage_mapping(in_fp, out_fp, body_site_column_name, exp_acronym):
    """Simplify the mapping file for use in figures

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

    mapping_lines = [l.strip().split('\t') for l in open(in_fp)]

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
            err_msg("Count not parse country %s" % country, new_line[0])
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

    output = open(out_fp,'w')
    output.write('\n'.join(['\t'.join(l) for l in new_mapping_lines]))
    output.write('\n')
    output.close()


def test_massage_mapping(ag_100nt):
    """Exercise the massage mapping code, verify expected results"""
    # output test file
    ag_100nt_m_TEST_fp = 'AG_100nt_TEST.txt'

    massage_mapping(ag_100nt, ag_100nt_m_TEST_fp, 'body_site', 'AGP')

    # verify the resulting header structure
    test_mapping = [l.strip().split('\t') for l in open(ag_100nt_m_TEST_fp)]
    test_header = test_mapping[0]
    test_header_length = len(test_header)
    assert test_header[-6:] == ['SIMPLE_BODY_SITE', 'TITLE_ACRONYM',
                                'TITLE_BODY_SITE', 'HMP_SITE',
                                'AGE_CATEGORY', 'BMI_CATEGORY']

    # verify each line in the test file
    for l in test_mapping[1:]:
        assert l[-6] in ['FECAL', 'SKIN', 'ORAL']
        assert l[-5] == 'AGP'
        acro, site = l[-4].split('-')
        assert acro == 'AGP'
        assert site in ['FECAL', 'SKIN', 'ORAL']
        assert l[-3] in ['FECAL', 'SKIN', 'ORAL']
        assert len(l) == test_header_length


def filter_mapping_file(in_fp, out_fp, columns_to_keep):
    """Filter out columns in a mapping file

    in_fp : the input file path
    out_fp : the output file path
    columns_to_keep : a dict of the columns to keep, valued by specific category
        value if desired to filter out samples that don't meet a given
        criteria. In other words, a row is retained if the function associated
        with the key "foo" returns True, or the row is retained if the value
        associated with "foo" is None.
    """
    lines = [l.strip().split('\t') for l in open(in_fp, 'U')]
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

    with open(out_fp, 'w') as out:
        out.write('\n'.join(['\t'.join(l) for l in new_lines]))
        out.write('\n')


def test_filter_mapping_file(ag_gg_mm_fp):
    # test output file
    fmf_TEST_fp = 'AG_GG_fmf_TEST.txt'

    # filter to just fecal samples, keep the age and title_acronym columns
    criteria = {'SIMPLE_BODY_SITE': lambda x: x == 'FECAL',
                'AGE': None,
                'TITLE_ACRONYM': None}
    filter_mapping_file(ag_gg_mm_fp, fmf_TEST_fp, criteria)

    # parse output
    test_mapping = [l.strip().split('\t') for l in open(fmf_TEST_fp)]

    # fish header, verify sanity of it
    test_header = test_mapping[0]
    assert len(test_header) == 4
    assert test_header[0] == '#SampleID'
    assert sorted(test_header) == sorted(['#SampleID', 'SIMPLE_BODY_SITE',
                                          'AGE', 'TITLE_ACRONYM'])

    # check each record
    test_sbs = test_header.index('SIMPLE_BODY_SITE')
    test_acro = test_header.index('TITLE_ACRONYM')
    for l in test_mapping[1:]:
        assert len(l) == 4
        assert l[test_sbs] == 'FECAL'
        assert l[test_acro] in ['AGP', 'GG']


def construct_svg_smash_commands(files, ids, cmd_format, cmd_args):
    """Format the SVG smashing commands"""
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
