# establish any minimals for the notebook environment
import os
import shutil
from distutils.spawn import find_executable

import qiime
import qiime_default_reference as qdr

import americangut as ag
from americangut.results_utils import get_repository_dir
from americangut.util import get_existing_path


_EBI_ACCESSIONS = ['ERP012803']
_TEST_ACCESSIONS = ['ag_testing']


# essential paths relative to the working_dir to be used between notebooks

paths = {
    # raw files
    'raw': {
        'sequences': '01-raw/sequences.fna',
        'metadata': '01-raw/metadata.txt',
        },

    # sequences filtered for blooms
    'filtered': {
        'sequences-notrim': '02-filtered/sequences-notrim.fna',
        'sequences-100nt': '02-filtered/sequences-100nt.fna',

        # only fecal sequences (for filtering for blooms)
        'fecal-sequences': '02-filtered/fecal-sequences.fna',

        # observed bloom sequences in samples
        'observed-blooms': '02-filtered/observed-blooms',
        'observed-blooms-biom': '02-filtered/observed-blooms/otu_table.biom',
        'observed-blooms-otu-map':
            ('02-filtered/observed-blooms/sortmerna_picked_otus/'
             'fecal-sequences_otus.txt'),
        },

    # resulting OTU data
    'otus': {
        'notrim': {
            'ag': '03-otus/notrim/gg-13_8-97-percent',
            'ag-biom': '03-otus/notrim/gg-13_8-97-percent/otu_table.biom'
            },
        '100nt': {
            'ag': '03-otus/100nt/gg-13_8-97-percent',
            'ag-biom': '03-otus/100nt/gg-13_8-97-percent/otu_table.biom'
        }
    },

    # merged files for diversity analyses
    'meta': {
        'ag-gg-100nt-biom': '04-meta/ag-gg-100nt.biom',
        'pgp-hmp-100nt-biom': '04-meta/pgp-hmp-100nt.biom',
        'ag-pgp-hmp-gg-100nt-biom': '04-meta/ag-pgp-hmp-gg-100nt.biom',
        'ag-cleaned-md': '04-meta/ag-cleaned.txt',
        'gg-cleaned-md': '04-meta/gg-cleaned.txt',
        'pgp-cleaned-md': '04-meta/pgp-cleaned.txt',
        'hmp-cleaned-md': '04-meta/hmp-cleaned.txt',
        'ag-gg-cleaned-md': '04-meta/ag-gg-cleaned.txt',
        'pgp-hmp-cleaned-md': '04-meta/pgp-hmp-cleaned.txt',
        'ag-pgp-hmp-gg-cleaned-md': '04-meta/ag-pgp-hmp-gg-cleaned.txt',
        },

    # alpha diversity analysis files
    'alpha': {
        '1k': {
            'ag-notrim-multiple': '05-alpha/1k/ag-notrim-multiple',
            'ag-notrim': '05-alpha/1k/ag-notrim',
            'ag-notrim-pd': '05-alpha/1k/ag-notrim/PD_whole_tree.txt',
            'ag-notrim-chao1': '05-alpha/1k/ag-notrim/chao1.txt',
            'ag-notrim-observedotus': ('05-alpha/1k/ag-notrim/'
                                       'observed_otus.txt'),
            'ag-notrim-shannon': '05-alpha/1k/ag-notrim/shannon.txt',

            'ag-pgp-hmp-gg-100nt-multiple': ('05-alpha/1k/'
                                             'ag-pgp-hmp-gg-100nt-multiple'),
            'ag-pgp-hmp-gg-100nt': '05-alpha/1k/ag-pgp-hmp-gg-100nt',
            'ag-pgp-hmp-gg-100nt-pd':
                '05-alpha/1k/ag-pgp-hmp-gg-100nt/PD_whole_tree.txt',
            'ag-pgp-hmp-gg-100nt-chao1':
                '05-alpha/1k/ag-pgp-hmp-gg-100nt/chao1.txt',
            'ag-pgp-hmp-gg-100nt-observedotus':
                '05-alpha/1k/ag-pgp-hmp-gg-100nt/observed_otus.txt',
            'ag-pgp-hmp-gg-100nt-shannon':
                '05-alpha/1k/ag-pgp-hmp-gg-100nt/shannon.txt',
            },
        '10k': {
            'ag-notrim-multiple': '05-alpha/10k/ag-notrim-multiple',
            'ag-notrim': '05-alpha/10k/ag-notrim',
            'ag-notrim-pd': '05-alpha/10k/ag-notrim/PD_whole_tree.txt',
            'ag-notrim-chao1': '05-alpha/10k/ag-notrim/chao1.txt',
            'ag-notrim-observedotus': ('05-alpha/10k/ag-notrim/'
                                       'observed_otus.txt'),
            'ag-notrim-shannon': '05-alpha/10k/ag-notrim/shannon.txt',

            'ag-pgp-hmp-gg-100nt-multiple': ('05-alpha/10k/'
                                             'ag-pgp-hmp-gg-100nt-multiple'),
            'ag-pgp-hmp-gg-100nt': '05-alpha/10k/ag-pgp-hmp-gg-100nt',
            'ag-pgp-hmp-gg-100nt-pd':
                '05-alpha/10k/ag-pgp-hmp-gg-100nt/PD_whole_tree.txt',
            'ag-pgp-hmp-gg-100nt-chao1':
                '05-alpha/10k/ag-pgp-hmp-gg-100nt/chao1.txt',
            'ag-pgp-hmp-gg-100nt-observedotus':
                '05-alpha/10k/ag-pgp-hmp-gg-100nt/observed_otus.txt',
            'ag-pgp-hmp-gg-100nt-shannon':
                '05-alpha/10k/ag-pgp-hmp-gg-100nt/shannon.txt',
            }
        },

    # beta diversity analysis files
    'beta': {
        'notrim': {
            '1k': {
                'ag-biom': '06-beta/notrim/1k/ag.biom',
                'ag': '06-beta/notrim/1k/ag',

                'ag-unifrac':
                    '06-beta/notrim/1k/ag/unweighted_unifrac_ag.txt',
                'ag-unifrac-pc':
                    '06-beta/notrim/1k/ag/unweighted_unifrac_ag-pc.txt',

                'ag-oral-unifrac':
                    '06-beta/notrim/1k/ag/unweighted_unifrac_ag-oral.txt',
                'ag-oral-unifrac-pc':
                    ('06-beta/notrim/1k/ag/'
                     'unweighted_unifrac_ag-oral-pc.txt'),

                'ag-skin-unifrac':
                    '06-beta/notrim/1k/ag/unweighted_unifrac_ag-skin.txt',
                'ag-skin-unifrac-pc':
                    ('06-beta/notrim/1k/ag/'
                     'unweighted_unifrac_ag-skin-pc.txt'),

                'ag-fecal-unifrac':
                    '06-beta/notrim/1k/ag/unweighted_unifrac_ag-fecal.txt',
                'ag-fecal-unifrac-pc':
                    ('06-beta/notrim/1k/ag/'
                     'unweighted_unifrac_ag-fecal-pc.txt'),
                'ag-wunifrac':
                    '06-beta/notrim/1k/ag/weighted_unifrac_ag.txt',
                'ag-wunifrac-pc':
                    '06-beta/notrim/1k/ag/weighted_unifrac_ag-pc.txt',

                'ag-oral-wunifrac':
                    '06-beta/notrim/1k/ag/weighted_unifrac_ag-oral.txt',
                'ag-oral-wunifrac-pc':
                    '06-beta/notrim/1k/ag/weighted_unifrac_ag-oral-pc.txt',

                'ag-skin-wunifrac':
                    '06-beta/notrim/1k/ag/weighted_unifrac_ag-skin.txt',
                'ag-skin-wunifrac-pc':
                    '06-beta/notrim/1k/ag/weighted_unifrac_ag-skin-pc.txt',

                'ag-fecal-wunifrac':
                    '06-beta/notrim/1k/ag/weighted_unifrac_ag-fecal.txt',
                'ag-fecal-wunifrac-pc':
                    '06-beta/notrim/1k/ag/weighted_unifrac_ag-fecal-pc.txt',
                    },

            '10k': {
                'ag-biom': '06-beta/notrim/10k/ag.biom',
                'ag': '06-beta/notrim/10k/ag',

                'ag-unifrac':
                    '06-beta/notrim/10k/ag/unweighted_unifrac_ag.txt',
                'ag-unifrac-pc':
                    '06-beta/notrim/10k/ag/unweighted_unifrac_ag-pc.txt',

                'ag-oral-unifrac':
                    '06-beta/notrim/10k/ag/unweighted_unifrac_ag-oral.txt',
                'ag-oral-unifrac-pc':
                    ('06-beta/notrim/10k/ag/'
                     'unweighted_unifrac_ag-oral-pc.txt'),

                'ag-skin-unifrac':
                    '06-beta/notrim/10k/ag/unweighted_unifrac_ag-skin.txt',
                'ag-skin-unifrac-pc':
                    ('06-beta/notrim/10k/ag/'
                     'unweighted_unifrac_ag-skin-pc.txt'),

                'ag-fecal-unifrac':
                    '06-beta/notrim/10k/ag/unweighted_unifrac_ag-fecal.txt',
                'ag-fecal-unifrac-pc':
                    ('06-beta/notrim/10k/ag/'
                     'unweighted_unifrac_ag-fecal-pc.txt'),
                'ag-wunifrac':
                    '06-beta/notrim/10k/ag/weighted_unifrac_ag.txt',
                'ag-wunifrac-pc':
                    '06-beta/notrim/10k/ag/weighted_unifrac_ag-pc.txt',

                'ag-oral-wunifrac':
                    '06-beta/notrim/10k/ag/weighted_unifrac_ag-oral.txt',
                'ag-oral-wunifrac-pc':
                    '06-beta/notrim/10k/ag/weighted_unifrac_ag-oral-pc.txt',

                'ag-skin-wunifrac':
                    '06-beta/notrim/10k/ag/weighted_unifrac_ag-skin.txt',
                'ag-skin-wunifrac-pc':
                    '06-beta/notrim/10k/ag/weighted_unifrac_ag-skin-pc.txt',

                'ag-fecal-wunifrac':
                    '06-beta/notrim/10k/ag/weighted_unifrac_ag-fecal.txt',
                'ag-fecal-wunifrac-pc':
                    ('06-beta/notrim/10k/ag/'
                     'weighted_unifrac_ag-fecal-pc.txt'),
                    },
            },
        '100nt': {
            '1k': {
                'ag-biom': '06-beta/100nt/1k/ag.biom',
                'ag': '06-beta/100nt/1k/ag',

                'ag-pgp-hmp-gg-biom': '06-beta/100nt/1k/ag-pgp-hmp-gg.biom',
                'ag-pgp-hmp-gg': '06-beta/100nt/1k/ag-pgp-hmp-gg',

                'ag-unifrac':
                    ('06-beta/100nt/1k/ag/'
                     'unweighted_unifrac_ag.txt'),
                'ag-unifrac-pc':
                    ('06-beta/100nt/1k/ag/'
                     'unweighted_unifrac_ag-pc.txt'),

                'ag-oral-unifrac':
                    ('06-beta/100nt/1k/ag/'
                     'unweighted_unifrac_ag-oral.txt'),
                'ag-oral-unifrac-pc':
                    ('06-beta/100nt/1k/ag/'
                     'unweighted_unifrac_ag-oral-pc.txt'),

                'ag-fecal-unifrac':
                    ('06-beta/100nt/1k/ag/'
                     'unweighted_unifrac_ag-fecal.txt'),
                'ag-fecal-unifrac-pc':
                    ('06-beta/100nt/1k/ag/'
                     'unweighted_unifrac_ag-fecal-pc.txt'),

                'ag-skin-unifrac':
                    ('06-beta/100nt/1k/ag/'
                     'unweighted_unifrac_ag-skin.txt'),
                'ag-skin-unifrac-pc':
                    ('06-beta/100nt/1k/ag/'
                     'unweighted_unifrac_ag-skin-pc.txt'),

                'ag-wunifrac':
                    ('06-beta/100nt/1k/ag/'
                     'weighted_unifrac_ag.txt'),
                'ag-wunifrac-pc':
                    ('06-beta/100nt/1k/ag/'
                     'weighted_unifrac_ag-pc.txt'),

                'ag-oral-wunifrac':
                    ('06-beta/100nt/1k/ag/'
                     'weighted_unifrac_ag-oral.txt'),
                'ag-oral-wunifrac-pc':
                    ('06-beta/100nt/1k/ag/'
                     'weighted_unifrac_ag-oral-pc.txt'),

                'ag-fecal-wunifrac':
                    ('06-beta/100nt/1k/ag/'
                     'weighted_unifrac_ag-fecal.txt'),
                'ag-fecal-wunifrac-pc':
                    ('06-beta/100nt/1k/ag/'
                     'weighted_unifrac_ag-fecal-pc.txt'),

                'ag-skin-wunifrac':
                    ('06-beta/100nt/1k/ag/'
                     'weighted_unifrac_ag-skin.txt'),
                'ag-skin-wunifrac-pc':
                    ('06-beta/100nt/1k/ag/'
                     'weighted_unifrac_ag-skin-pc.txt'),

                'ag-pgp-hmp-gg-unifrac':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-pgp-hmp-gg.txt'),
                'ag-pgp-hmp-gg-unifrac-pc':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-pgp-hmp-gg-pc.txt'),

                'ag-pgp-hmp-gg-skin-unifrac':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-pgp-hmp-gg-skin.txt'),
                'ag-pgp-hmp-gg-skin-unifrac-pc':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-pgp-hmp-gg-skin-pc.txt'),

                'ag-pgp-hmp-gg-oral-unifrac':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-pgp-hmp-gg-oral.txt'),
                'ag-pgp-hmp-gg-oral-unifrac-pc':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-pgp-hmp-gg-oral-pc.txt'),

                'ag-pgp-hmp-gg-fecal-unifrac':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-pgp-hmp-gg-fecal.txt'),
                'ag-pgp-hmp-gg-fecal-unifrac-pc':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-pgp-hmp-gg-fecal-pc.txt'),

                'ag-pgp-hmp-gg-skin-wunifrac':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-pgp-hmp-gg-skin.txt'),
                'ag-pgp-hmp-gg-skin-wunifrac-pc':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-pgp-hmp-gg-skin-pc.txt'),

                'ag-pgp-hmp-gg-oral-wunifrac':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-pgp-hmp-gg-oral.txt'),
                'ag-pgp-hmp-gg-oral-wunifrac-pc':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-pgp-hmp-gg-oral-pc.txt'),

                'ag-pgp-hmp-gg-fecal-wunifrac':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-pgp-hmp-gg-fecal.txt'),
                'ag-pgp-hmp-gg-fecal-wunifrac-pc':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-pgp-hmp-gg-fecal-pc.txt'),

                'ag-pgp-hmp-gg-wunifrac':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-pgp-hmp-gg.txt'),
                'ag-pgp-hmp-gg-wunifrac-pc':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-pgp-hmp-gg-pc.txt'),

                'ag-gg-unifrac':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-gg.txt'),
                'ag-gg-unifrac-pc':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-gg-pc.txt'),

                'ag-gg-skin-unifrac':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-gg-skin.txt'),
                'ag-gg-skin-unifrac-pc':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-gg-skin-pc.txt'),

                'ag-gg-oral-unifrac':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-gg-oral.txt'),
                'ag-gg-oral-unifrac-pc':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-gg-oral-pc.txt'),

                'ag-gg-fecal-unifrac':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-gg-fecal.txt'),
                'ag-gg-fecal-unifrac-pc':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-gg-fecal-pc.txt'),

                'ag-gg-wunifrac':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-gg.txt'),
                'ag-gg-wunifrac-pc':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-gg-pc.txt'),

                'ag-gg-skin-wunifrac':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-gg-skin.txt'),
                'ag-gg-skin-wunifrac-pc':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-gg-skin-pc.txt'),

                'ag-gg-oral-wunifrac':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-gg-oral.txt'),
                'ag-gg-oral-wunifrac-pc':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-gg-oral-pc.txt'),

                'ag-gg-fecal-wunifrac':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-gg-fecal.txt'),
                'ag-gg-fecal-wunifrac-pc':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-gg-fecal-pc.txt'),

                'ag-gg-subsampled-unifrac':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-gg-subsampled.txt'),
                'ag-gg-subsampled-unifrac-pc':
                    ('06-beta/100nt/1k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-gg-subsampled-pc.txt'),
            },
            '10k': {
                'ag-biom': '06-beta/100nt/10k/ag.biom',
                'ag': '06-beta/100nt/10k/ag',

                'ag-pgp-hmp-gg-biom': '06-beta/100nt/10k/ag-pgp-hmp-gg.biom',
                'ag-pgp-hmp-gg': '06-beta/100nt/10k/ag-pgp-hmp-gg',

                'ag-unifrac':
                    ('06-beta/100nt/10k/ag/'
                     'unweighted_unifrac_ag.txt'),
                'ag-unifrac-pc':
                    ('06-beta/100nt/10k/ag/'
                     'unweighted_unifrac_ag-pc.txt'),

                'ag-oral-unifrac':
                    ('06-beta/100nt/10k/ag/'
                     'unweighted_unifrac_ag-oral.txt'),
                'ag-oral-unifrac-pc':
                    ('06-beta/100nt/10k/ag/'
                     'unweighted_unifrac_ag-oral-pc.txt'),

                'ag-fecal-unifrac':
                    ('06-beta/100nt/10k/ag/'
                     'unweighted_unifrac_ag-fecal.txt'),
                'ag-fecal-unifrac-pc':
                    ('06-beta/100nt/10k/ag/'
                     'unweighted_unifrac_ag-fecal-pc.txt'),

                'ag-skin-unifrac':
                    ('06-beta/100nt/10k/ag/'
                     'unweighted_unifrac_ag-skin.txt'),
                'ag-skin-unifrac-pc':
                    ('06-beta/100nt/10k/ag/'
                     'unweighted_unifrac_ag-skin-pc.txt'),

                'ag-wunifrac':
                    ('06-beta/100nt/10k/ag/'
                     'weighted_unifrac_ag.txt'),
                'ag-wunifrac-pc':
                    ('06-beta/100nt/10k/ag/'
                     'weighted_unifrac_ag-pc.txt'),

                'ag-oral-wunifrac':
                    ('06-beta/100nt/10k/ag/'
                     'weighted_unifrac_ag-oral.txt'),
                'ag-oral-wunifrac-pc':
                    ('06-beta/100nt/10k/ag/'
                     'weighted_unifrac_ag-oral-pc.txt'),

                'ag-fecal-wunifrac':
                    ('06-beta/100nt/10k/ag/'
                     'weighted_unifrac_ag-fecal.txt'),
                'ag-fecal-wunifrac-pc':
                    ('06-beta/100nt/10k/ag/'
                     'weighted_unifrac_ag-fecal-pc.txt'),

                'ag-skin-wunifrac':
                    ('06-beta/100nt/10k/ag/'
                     'weighted_unifrac_ag-skin.txt'),
                'ag-skin-wunifrac-pc':
                    ('06-beta/100nt/10k/ag/'
                     'weighted_unifrac_ag-skin-pc.txt'),

                'ag-pgp-hmp-gg-unifrac':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-pgp-hmp-gg.txt'),
                'ag-pgp-hmp-gg-unifrac-pc':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-pgp-hmp-gg-pc.txt'),

                'ag-pgp-hmp-gg-skin-unifrac':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-pgp-hmp-gg-skin.txt'),
                'ag-pgp-hmp-gg-skin-unifrac-pc':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-pgp-hmp-gg-skin-pc.txt'),

                'ag-pgp-hmp-gg-oral-unifrac':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-pgp-hmp-gg-oral.txt'),
                'ag-pgp-hmp-gg-oral-unifrac-pc':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-pgp-hmp-gg-oral-pc.txt'),

                'ag-pgp-hmp-gg-fecal-unifrac':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-pgp-hmp-gg-fecal.txt'),
                'ag-pgp-hmp-gg-fecal-unifrac-pc':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-pgp-hmp-gg-fecal-pc.txt'),

                'ag-pgp-hmp-gg-skin-wunifrac':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-pgp-hmp-gg-skin.txt'),
                'ag-pgp-hmp-gg-skin-wunifrac-pc':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-pgp-hmp-gg-skin-pc.txt'),

                'ag-pgp-hmp-gg-oral-wunifrac':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-pgp-hmp-gg-oral.txt'),
                'ag-pgp-hmp-gg-oral-wunifrac-pc':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-pgp-hmp-gg-oral-pc.txt'),

                'ag-pgp-hmp-gg-fecal-wunifrac':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-pgp-hmp-gg-fecal.txt'),
                'ag-pgp-hmp-gg-fecal-wunifrac-pc':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-pgp-hmp-gg-fecal-pc.txt'),

                'ag-pgp-hmp-gg-wunifrac':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-pgp-hmp-gg.txt'),
                'ag-pgp-hmp-gg-wunifrac-pc':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-pgp-hmp-gg-pc.txt'),

                'ag-gg-unifrac':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-gg.txt'),
                'ag-gg-unifrac-pc':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-gg-pc.txt'),

                'ag-gg-skin-unifrac':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-gg-skin.txt'),
                'ag-gg-skin-unifrac-pc':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-gg-skin-pc.txt'),

                'ag-gg-oral-unifrac':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-gg-oral.txt'),
                'ag-gg-oral-unifrac-pc':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-gg-oral-pc.txt'),

                'ag-gg-fecal-unifrac':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-gg-fecal.txt'),
                'ag-gg-fecal-unifrac-pc':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-gg-fecal-pc.txt'),

                'ag-gg-subsampled-unifrac':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-gg-subsampled.txt'),
                'ag-gg-subsampled-unifrac-pc':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'unweighted_unifrac_ag-gg-subsampled-pc.txt'),

                'ag-gg-wunifrac':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-gg.txt'),
                'ag-gg-wunifrac-pc':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-gg-pc.txt'),

                'ag-gg-skin-wunifrac':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-gg-skin.txt'),
                'ag-gg-skin-wunifrac-pc':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-gg-skin-pc.txt'),

                'ag-gg-oral-wunifrac':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-gg-oral.txt'),
                'ag-gg-oral-wunifrac-pc':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-gg-oral-pc.txt'),

                'ag-gg-fecal-wunifrac':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-gg-fecal.txt'),
                'ag-gg-fecal-wunifrac-pc':
                    ('06-beta/100nt/10k/ag-pgp-hmp-gg/'
                     'weighted_unifrac_ag-gg-fecal-pc.txt'),
            }
        },
    },
    # taxonomy summaries
    'taxa': {
        '100nt': {
            'base': '07-taxa/100nt',
            'L2': {
                'ag-tsv': '07-taxa/100nt/otu_table_L2.txt',
                'ag-biom': '07-taxa/100nt/otu_table_L2.biom',
                'ag-md': '07-taxa/100nt/ag-cleaned_L2.txt',
                'ag-skin-biom': '07-taxa/100nt/otu_table_skin_L2.biom',
                'ag-oral-biom': '07-taxa/100nt/otu_table_oral_L2.biom',
                'ag-fecal-biom': '07-taxa/100nt/otu_table_fecal_L2.biom',
                },
            'L3': {
                'ag-tsv': '07-taxa/100nt/otu_table_L3.txt',
                'ag-biom': '07-taxa/100nt/otu_table_L3.biom',
                'ag-md': '07-taxa/100nt/ag-cleaned_L3.txt',
                'ag-skin-biom': '07-taxa/100nt/otu_table_skin_L3.biom',
                'ag-oral-biom': '07-taxa/100nt/otu_table_oral_L3.biom',
                'ag-fecal-biom': '07-taxa/100nt/otu_table_fecal_L3.biom',
                },
            'L6': {
                'ag-tsv': '07-taxa/100nt/otu_table_L6.txt',
                'ag-biom': '07-taxa/100nt/otu_table_L6.biom',
                'ag-md': '07-taxa/100nt/ag-cleaned_L6.txt',
                'ag-skin-biom': '07-taxa/100nt/otu_table_skin_L6.biom',
                'ag-oral-biom': '07-taxa/100nt/otu_table_oral_L6.biom',
                'ag-fecal-biom': '07-taxa/100nt/otu_table_fecal_L6.biom',
                }
            },
        'notrim': {
            'base': '07-taxa/notrim',
            'L2': {
                'ag-tsv': '07-taxa/notrim/otu_table_L2.txt',
                'ag-biom': '07-taxa/notrim/otu_table_L2.biom',
                'ag-md': '07-taxa/notrim/ag-cleaned_L2.txt',
                'ag-skin-biom': '07-taxa/notrim/otu_table_skin_L2.biom',
                'ag-oral-biom': '07-taxa/notrim/otu_table_oral_L2.biom',
                'ag-fecal-biom': '07-taxa/notrim/otu_table_fecal_L2.biom',
                },
            'L3': {
                'ag-tsv': '07-taxa/notrim/otu_table_L3.txt',
                'ag-biom': '07-taxa/notrim/otu_table_L3.biom',
                'ag-md': '07-taxa/notrim/ag-cleaned_L3.txt',
                'ag-skin-biom': '07-taxa/notrim/otu_table_skin_L3.biom',
                'ag-oral-biom': '07-taxa/notrim/otu_table_oral_L3.biom',
                'ag-fecal-biom': '07-taxa/notrim/otu_table_fecal_L3.biom',
                },
            'L6': {
                'ag-tsv': '07-taxa/notrim/otu_table_L6.txt',
                'ag-biom': '07-taxa/notrim/otu_table_L6.biom',
                'ag-md': '07-taxa/notrim/ag-cleaned_L6.txt',
                'ag-skin-biom': '07-taxa/notrim/otu_table_skin_L6.biom',
                'ag-oral-biom': '07-taxa/notrim/otu_table_oral_L6.biom',
                'ag-fecal-biom': '07-taxa/notrim/otu_table_fecal_L6.biom',
                }
            }
        },

    # collapsed samples
    'collapsed': {
        '100nt': {
            'alpha-map': '08-collapsed/100nt/alpha_map.txt',
            '1k': {
                'ag-biom': '08-collapsed/100nt/1k/ag.biom',
                'ag-fecal-biom': '08-collapsed/100nt/1k/ag-fecal.biom',
                'ag-skin-biom': '08-collapsed/100nt/1k/ag-oral.biom',
                'ag-oral-biom': '08-collapsed/100nt/1k/ag-skin.biom',

                'ag-fecal-sex-biom':  ('08-collapsed/100nt/1k/'
                                       'ag-fecal-sex.biom'),
                'ag-fecal-diet-biom': ('08-collapsed/100nt/1k/'
                                       'ag-fecal-diet.biom'),
                'ag-fecal-age-biom':  ('08-collapsed/100nt/1k/'
                                       'ag-fecal-age.biom'),
                'ag-fecal-bmi-biom':  ('08-collapsed/100nt/1k/'
                                       'ag-fecal-bmi.biom'),

                'ag-oral-sex-biom':  '08-collapsed/100nt/1k/ag-oral-sex.biom',
                'ag-oral-diet-biom': '08-collapsed/100nt/1k/ag-oral-diet.biom',
                'ag-oral-age-biom':  '08-collapsed/100nt/1k/ag-oral-age.biom',
                'ag-oral-flossing-biom':  ('08-collapsed/100nt/1k/'
                                           'ag-oral-flossing.biom'),

                'ag-skin-sex-biom':  '08-collapsed/100nt/1k/ag-skin-sex.biom',
                'ag-skin-cosmetics-biom': ('08-collapsed/100nt/1k/'
                                           'ag-skin-cosmetics.biom'),
                'ag-skin-age-biom': '08-collapsed/100nt/1k/ag-skin-age.biom',
                'ag-skin-hand-biom': ('08-collapsed/100nt/1k/'
                                      'ag-skin-hand.biom'),
                },
            '10k': {
                'ag-biom': '08-collapsed/100nt/10k/ag.biom',
                'ag-fecal-biom': '08-collapsed/100nt/10k/ag-fecal.biom',
                'ag-skin-biom': '08-collapsed/100nt/10k/ag-oral.biom',
                'ag-oral-biom': '08-collapsed/100nt/10k/ag-skin.biom',

                'ag-fecal-sex-biom': ('08-collapsed/100nt/10k/'
                                      'ag-fecal-sex.biom'),
                'ag-fecal-diet-biom': ('08-collapsed/100nt/10k/'
                                       'ag-fecal-diet.biom'),
                'ag-fecal-age-biom': ('08-collapsed/100nt/10k/'
                                      'ag-fecal-age.biom'),
                'ag-fecal-bmi-biom': ('08-collapsed/100nt/10k/'
                                      'ag-fecal-bmi.biom'),

                'ag-oral-sex-biom':  '08-collapsed/100nt/10k/ag-oral-sex.biom',
                'ag-oral-diet-biom': ('08-collapsed/100nt/10k/'
                                      'ag-oral-diet.biom'),
                'ag-oral-age-biom': '08-collapsed/100nt/10k/ag-oral-age.biom',
                'ag-oral-flossing-biom': ('08-collapsed/100nt/10k/'
                                          'ag-oral-flossing.biom'),

                'ag-skin-sex-biom': '08-collapsed/100nt/10k/ag-skin-sex.biom',
                'ag-skin-cosmetics-biom': ('08-collapsed/100nt/10k/'
                                           'ag-skin-cosmetics.biom'),
                'ag-skin-age-biom': '08-collapsed/100nt/10k/ag-skin-age.biom',
                'ag-skin-hand-biom': ('08-collapsed/100nt/10k/'
                                      'ag-skin-hand.biom'),
                }
            },
        'notrim': {
            'alpha-map': '08-collapsed/notrim/alpha-map.txt',
            '1k': {
                'ag-biom': '08-collapsed/notrim/1k/ag.biom',
                'ag-fecal-biom': '08-collapsed/notrim/1k/ag-fecal.biom',
                'ag-skin-biom': '08-collapsed/notrim/1k/ag-oral.biom',
                'ag-oral-biom': '08-collapsed/notrim/1k/ag-skin.biom',

                'ag-fecal-sex-biom': ('08-collapsed/notrim/1k/'
                                      'ag-fecal-sex.biom'),
                'ag-fecal-diet-biom': ('08-collapsed/notrim/1k/'
                                       'ag-fecal-diet.biom'),
                'ag-fecal-age-biom': ('08-collapsed/notrim/1k/'
                                      'ag-fecal-age.biom'),
                'ag-fecal-bmi-biom': ('08-collapsed/notrim/1k/'
                                      'ag-fecal-bmi.biom'),

                'ag-oral-sex-biom':  '08-collapsed/notrim/1k/ag-oral-sex.biom',
                'ag-oral-diet-biom': ('08-collapsed/notrim/1k/'
                                      'ag-oral-diet.biom'),
                'ag-oral-age-biom':  '08-collapsed/notrim/1k/ag-oral-age.biom',
                'ag-oral-flossing-biom': ('08-collapsed/notrim/1k/'
                                          'ag-oral-flossing.biom'),

                'ag-skin-sex-biom':  '08-collapsed/notrim/1k/ag-skin-sex.biom',
                'ag-skin-cosmetics-biom': ('08-collapsed/notrim/1k/'
                                           'ag-skin-cosmetics.biom'),
                'ag-skin-age-biom':  '08-collapsed/notrim/1k/ag-skin-age.biom',
                'ag-skin-hand-biom': ('08-collapsed/notrim/1k/'
                                      'ag-skin-hand.biom'),
                },
            '10k': {
                'ag-biom': '08-collapsed/notrim/10k/ag.biom',
                'ag-fecal-biom': '08-collapsed/notrim/10k/ag-fecal.biom',
                'ag-skin-biom': '08-collapsed/notrim/10k/ag-oral.biom',
                'ag-oral-biom': '08-collapsed/notrim/10k/ag-skin.biom',

                'ag-fecal-sex-biom': ('08-collapsed/notrim/10k/'
                                      'ag-fecal-sex.biom'),
                'ag-fecal-diet-biom': ('08-collapsed/notrim/10k/'
                                       'ag-fecal-diet.biom'),
                'ag-fecal-age-biom':  ('08-collapsed/notrim/10k/'
                                       'ag-fecal-age.biom'),
                'ag-fecal-bmi-biom':  ('08-collapsed/notrim/10k/'
                                       'ag-fecal-bmi.biom'),

                'ag-oral-sex-biom': '08-collapsed/notrim/10k/ag-oral-sex.biom',
                'ag-oral-diet-biom': ('08-collapsed/notrim/10k/'
                                      'ag-oral-diet.biom'),
                'ag-oral-age-biom': '08-collapsed/notrim/10k/ag-oral-age.biom',
                'ag-oral-flossing-biom': ('08-collapsed/notrim/10k/'
                                          'ag-oral-flossing.biom'),

                'ag-skin-sex-biom': '08-collapsed/notrim/10k/ag-skin-sex.biom',
                'ag-skin-cosmetics-biom': ('08-collapsed/notrim/10k/'
                                           'ag-skin-cosmetics.biom'),
                'ag-skin-age-biom': '08-collapsed/notrim/10k/ag-skin-age.biom',
                'ag-skin-hand-biom':  ('08-collapsed/notrim/10k/'
                                       'ag-skin-hand.biom'),
                }
            }
        },


    # per-sample results
    'per-sample': {
        'successful-ids': '09-per-sample/successful_ids.txt',
        'unsuccessful-ids': '09-per-sample/unsuccessful_ids.txt',
        'results': '09-per-sample/results',
        'statics-fecal': '09-per-sample/statics-fecal',
        'statics-oral': '09-per-sample/statics-oral',
        'statics-skin': '09-per-sample/statics-skin',
    },

    'populated-templates': {
        'result-pdfs': '10-populated-templates/pdfs/',
        'result-taxa': '10-populated-templates/taxa/',
        'successful-pdfs': '10-populated-templates/successful_ids.txt',
        'unsuccessful-pdfs': '10-populated-templates/unsuccessful_ids.txt',
    },

    'demux': {},
    'otu-map': {},

    'package': {
        'split': {
            'notrim-map': '01-raw/metadata.txt',
            'notrim-raw-otu': ('03-otus/notrim/gg-13_8-97-percent/'
                               'otu_table.biom'),
            'notrim-raw-dir': '11-packaged/split/notrim/raw/',
            'notrim-1k-otu': '06-beta/notrim/1k/ag.biom',
            'notrim-1k-dir': '11-packaged/split/notrim/1k/',
            'notrim-10k-otu': '06-beta/notrim/10k/ag.biom',
            'notrim-10k-dir': '11-packaged/split/notrim/10k/',

            '100nt-map': '01-raw/metadata.txt',
            '100nt-raw-otu': '03-otus/100nt/gg-13_8-97-percent/otu_table.biom',
            '100nt-raw-dir': '11-packaged/split/100nt/raw/',
            '100nt-1k-otu': '06-beta/100nt/1k/ag.biom',
            '100nt-1k-dir': '11-packaged/split/100nt/1k/',
            '100nt-10k-otu': '06-beta/100nt/10k/ag.biom',
            '100nt-10k-dir': '11-packaged/split/100nt/10k/',
            },

        'single_ids': {
            'fecal-1k': '11-packaged/fecal/single_ids_1k.txt',
            'fecal-10k': '11-packaged/fecal/single_ids_10k.txt',
            'fecal-unrare': '11-packaged/fecal/single_ids_unrarefied.txt',

            'oral-1k': '11-packaged/oral/single_ids_1k.txt',
            'oral-10k': '11-packaged/oral/single_ids_10k.txt',
            'oral-unrare': '11-packaged/oral/single_ids_unrarefied.txt',

            'skin-1k': '11-packaged/skin/single_ids_1k.txt',
            'skin-10k': '11-packaged/skin/single_ids_10k.txt',
            'skin-unrare': '11-packaged/skin/single_ids_unrarefied.txt',
            },

        'all_participants_all_samples': {
            'fecal': {
                'notrim': {
                    'source-unrare-otu': ('11-packaged/split/notrim/raw/'
                                          'otu_table__BODY_HABITAT_UBERON:'
                                          'feces__.biom'),
                    'source-unrare-map': ('11-packaged/split/notrim/raw/'
                                          'metadata__BODY_HABITAT_UBERON:'
                                          'feces__.txt'),
                    'sink-unrare-dir': ('11-packaged/fecal/notrim'),
                    'sink-unrare-otu': ('11-packaged/fecal/notrim/ag_fecal.'
                                        'biom'),
                    'sink-unrare-map': ('11-packaged/fecal/notrim/ag_fecal.'
                                        'txt'),

                    'source-1k-otu': ('11-packaged/split/notrim/1k/'
                                      'ag__BODY_HABITAT_UBERON:feces__.biom'),
                    'source-1k-map': ('11-packaged/split/notrim/1k/'
                                      'metadata__BODY_HABITAT_UBERON:'
                                      'feces__.txt'),
                    'source-1k-pd': '05-alpha/1k/ag-notrim/PD_whole_tree.txt',
                    'source-1k-chao1': '05-alpha/1k/ag-notrim/chao1.txt',
                    'source-1k-observedotus': ('05-alpha/1k/ag-notrim/'
                                               'observed_otus.txt'),
                    'source-1k-shannon': '05-alpha/1k/ag-notrim/shannon.txt',
                    'source-1k-unweighted-unifrac': ('06-beta/notrim/1k/ag/'
                                                     'unweighted_unifrac_ag-'
                                                     'fecal.txt'),
                    'source-1k-weighted-unifrac': ('06-beta/notrim/1k/ag/'
                                                   'weighted_unifrac_ag-'
                                                   'fecal.txt'),

                    'source-10k-otu': ('11-packaged/split/notrim/10k/'
                                       'ag__BODY_HABITAT_UBERON:feces__.biom'),
                    'source-10k-map': ('11-packaged/split/notrim/10k/'
                                       'metadata__BODY_HABITAT_UBERON:feces__.'
                                       'txt'),
                    'source-10k-pd': ('05-alpha/10k/ag-notrim/PD_whole_tree.'
                                      'txt'),
                    'source-10k-chao1': ('05-alpha/10k/ag-notrim/chao1.'
                                         'txt'),
                    'source-10k-observedotus': ('05-alpha/10k/ag-notrim/'
                                                'observed_otus.txt'),
                    'source-10k-shannon': '05-alpha/10k/ag-notrim/shannon.txt',
                    'source-10k-unweighted-unifrac': ('06-beta/notrim/10k/ag/'
                                                      'unweighted_unifrac_ag'
                                                      '-fecal.txt'),
                    'source-10k-weighted-unifrac': ('06-beta/notrim/10k/ag/'
                                                    'weighted_unifrac_ag-'
                                                    'fecal.txt'),

                    'sink-1k-dir': ('11-packaged/fecal/notrim/'
                                    'all_participants/all_samples/1k'),
                    'sink-1k-otu': ('11-packaged/fecal/notrim/'
                                    'all_participants/all_samples/1k/'
                                    'ag_1k_fecal.biom'),
                    'sink-1k-map': ('11-packaged/fecal/notrim/'
                                    'all_participants/all_samples/1k/'
                                    'ag_1k_fecal.txt'),
                    'sink-1k-unweighted-unifrac': ('11-packaged/fecal/notrim/'
                                                   'all_participants/'
                                                   'all_samples/1k/'
                                                   'unweighted_unifrac_ag_1k_'
                                                   'fecal.txt'),
                    'sink-1k-weighted-unifrac': ('11-packaged/fecal/notrim/'
                                                 'all_participants/'
                                                 'all_samples/1k/'
                                                 'weighted_unifrac_ag_1k_'
                                                 'fecal.txt'),

                    'sink-10k-dir': ('11-packaged/fecal/notrim/'
                                     'all_participants/all_samples/10k'),
                    'sink-10k-otu': ('11-packaged/fecal/notrim/'
                                     'all_participants/all_samples/10k/'
                                     'ag_10k_fecal.biom'),
                    'sink-10k-map': ('11-packaged/fecal/notrim/'
                                     'all_participants/all_samples/10k/'
                                     'ag_10k_fecal.txt'),
                    'sink-10k-unweighted-unifrac': ('11-packaged/fecal/notrim/'
                                                    'all_participants/'
                                                    'all_samples/10k/'
                                                    'unweighted_unifrac_ag_'
                                                    '10k_fecal.txt'),
                    'sink-10k-weighted-unifrac': ('11-packaged/fecal/notrim/'
                                                  'all_participants/'
                                                  'all_samples/10k/'
                                                  'weighted_unifrac_ag_10k_'
                                                  'fecal.txt'),
                    },
                '100nt': {
                    'source-unrare-otu': ('11-packaged/split/100nt/raw/'
                                          'otu_table__BODY_HABITAT_UBERON:'
                                          'feces__.biom'),
                    'source-unrare-map': ('11-packaged/split/100nt/raw/'
                                          'metadata__BODY_HABITAT_UBERON:'
                                          'feces__.txt'),
                    'sink-unrare-dir': ('11-packaged/fecal/100nt'),
                    'sink-unrare-otu': ('11-packaged/fecal/100nt/ag_fecal.'
                                        'biom'),
                    'sink-unrare-map': ('11-packaged/fecal/100nt/ag_fecal.'
                                        'txt'),

                    'source-1k-otu': ('11-packaged/split/100nt/1k/'
                                      'ag__BODY_HABITAT_UBERON:feces__.biom'),
                    'source-1k-map': ('11-packaged/split/100nt/1k/'
                                      'metadata__BODY_HABITAT_UBERON:'
                                      'feces__.txt'),
                    'source-1k-pd': ('05-alpha/1k/ag-pgp-hmp-gg-100nt/'
                                     'PD_whole_tree.txt'),
                    'source-1k-chao1': ('05-alpha/1k/ag-pgp-hmp-gg-100nt/'
                                        'chao1.txt'),
                    'source-1k-observedotus': ('05-alpha/1k/'
                                               'ag-pgp-hmp-gg-100nt/'
                                               'observed_otus.txt'),
                    'source-1k-shannon': ('05-alpha/1k/ag-pgp-hmp-gg-100nt/'
                                          'shannon.txt'),
                    'source-1k-unweighted-unifrac': ('06-beta/100nt/1k/ag/'
                                                     'unweighted_unifrac_ag-'
                                                     'fecal.txt'),
                    'source-1k-weighted-unifrac': ('06-beta/100nt/1k/ag/'
                                                   'weighted_unifrac_ag-'
                                                   'fecal.txt'),

                    'source-10k-otu': ('11-packaged/split/100nt/10k/'
                                       'ag__BODY_HABITAT_UBERON:feces__.biom'),
                    'source-10k-map': ('11-packaged/split/100nt/10k/metadata'
                                       '__BODY_HABITAT_UBERON:feces__.txt'),
                    'source-10k-pd': ('05-alpha/10k/ag-pgp-hmp-gg-100nt/'
                                      'PD_whole_tree.txt'),
                    'source-10k-chao1': ('05-alpha/10k/ag-pgp-hmp-gg-100nt/'
                                         'chao1.txt'),
                    'source-10k-observedotus': ('05-alpha/10k/'
                                                'ag-pgp-hmp-gg-100nt/'
                                                'observed_otus.txt'),
                    'source-10k-shannon': ('05-alpha/10k/ag-pgp-hmp-gg-100nt/'
                                           'shannon.txt'),
                    'source-10k-unweighted-unifrac': ('06-beta/100nt/10k/ag/'
                                                      'unweighted_unifrac_ag-'
                                                      'fecal.txt'),
                    'source-10k-weighted-unifrac': ('06-beta/100nt/10k/ag/'
                                                    'weighted_unifrac_ag-'
                                                    'fecal.txt'),

                    'sink-1k-dir': ('11-packaged/fecal/100nt/all_participants/'
                                    'all_samples/1k'),
                    'sink-1k-otu': ('11-packaged/fecal/100nt/all_participants/'
                                    'all_samples/1k/ag_1k_fecal.biom'),
                    'sink-1k-map': ('11-packaged/fecal/100nt/all_participants/'
                                    'all_samples/1k/ag_1k_fecal.txt'),
                    'sink-1k-unweighted-unifrac': ('11-packaged/fecal/100nt/'
                                                   'all_participants/'
                                                   'all_samples/1k/'
                                                   'unweighted_unifrac_ag_1k_'
                                                   'fecal.txt'),
                    'sink-1k-weighted-unifrac': ('11-packaged/fecal/100nt/'
                                                 'all_participants/'
                                                 'all_samples/1k/'
                                                 'weighted_unifrac_ag_1k_'
                                                 'fecal.txt'),

                    'sink-10k-dir': ('11-packaged/fecal/100nt/'
                                     'all_participants/all_samples/10k'),
                    'sink-10k-otu': ('11-packaged/fecal/100nt/'
                                     'all_participants/all_samples/10k/'
                                     'ag_10k_fecal.biom'),
                    'sink-10k-map': ('11-packaged/fecal/100nt/'
                                     'all_participants/all_samples/10k/'
                                     'ag_10k_fecal.txt'),
                    'sink-10k-unweighted-unifrac': ('11-packaged/fecal/100nt/'
                                                    'all_participants/'
                                                    'all_samples/10k/'
                                                    'unweighted_unifrac_ag_'
                                                    '10k_fecal.txt'),
                    'sink-10k-weighted-unifrac': ('11-packaged/fecal/100nt/'
                                                  'all_participants/'
                                                  'all_samples/10k/'
                                                  'weighted_unifrac_ag_10k_'
                                                  'fecal.txt'),
                    },
                },
            'oral': {
                'notrim': {
                    'source-unrare-otu': ('11-packaged/split/notrim/raw/'
                                          'otu_table__BODY_HABITAT_UBERON:'
                                          'oral cavity__.biom'),
                    'source-unrare-map': ('11-packaged/split/notrim/raw/'
                                          'metadata__BODY_HABITAT_UBERON:'
                                          'oral cavity__.txt'),
                    'sink-unrare-dir': ('11-packaged/oral/notrim'),
                    'sink-unrare-otu': ('11-packaged/oral/notrim/ag_oral.'
                                        'biom'),
                    'sink-unrare-map': ('11-packaged/oral/notrim/ag_oral.'
                                        'txt'),

                    'source-1k-otu': ('11-packaged/split/notrim/1k/'
                                      'ag__BODY_HABITAT_UBERON:oral '
                                      'cavity__.biom'),
                    'source-1k-map': ('11-packaged/split/notrim/1k/'
                                      'metadata__BODY_HABITAT_UBERON:oral '
                                      'cavity__.txt'),
                    'source-1k-pd': '05-alpha/1k/ag-notrim/PD_whole_tree.txt',
                    'source-1k-chao1': ('05-alpha/1k/ag-notrim/chao1.txt'),
                    'source-1k-observedotus': ('05-alpha/1k/ag-notrim/'
                                               'observed_otus.txt'),
                    'source-1k-shannon': '05-alpha/1k/ag-notrim/shannon.txt',
                    'source-1k-unweighted-unifrac': ('06-beta/notrim/1k/ag/'
                                                     'unweighted_unifrac_ag-'
                                                     'oral.txt'),
                    'source-1k-weighted-unifrac': ('06-beta/notrim/1k/ag/'
                                                   'weighted_unifrac_ag-oral.'
                                                   'txt'),

                    'source-10k-otu': ('11-packaged/split/notrim/10k/'
                                       'ag__BODY_HABITAT_UBERON:oral '
                                       'cavity__.biom'),
                    'source-10k-map': ('11-packaged/split/notrim/10k/'
                                       'metadata__BODY_HABITAT_UBERON:oral '
                                       'cavity__.txt'),
                    'source-10k-pd': ('05-alpha/10k/ag-notrim/'
                                      'PD_whole_tree.txt'),
                    'source-10k-chao1': '05-alpha/10k/ag-notrim/chao1.txt',
                    'source-10k-observedotus': ('05-alpha/10k/ag-notrim/'
                                                'observed_otus.txt'),
                    'source-10k-shannon': '05-alpha/10k/ag-notrim/shannon.txt',
                    'source-10k-unweighted-unifrac': ('06-beta/notrim/10k/ag/'
                                                      'unweighted_unifrac_ag'
                                                      '-oral.txt'),
                    'source-10k-weighted-unifrac': ('06-beta/notrim/10k/ag/'
                                                    'weighted_unifrac_ag-'
                                                    'oral.txt'),

                    'sink-1k-dir': ('11-packaged/oral/notrim/all_participants/'
                                    'all_samples/1k'),
                    'sink-1k-otu': ('11-packaged/oral/notrim/all_participants/'
                                    'all_samples/1k/ag_1k_oral.biom'),
                    'sink-1k-map': ('11-packaged/oral/notrim/all_participants/'
                                    'all_samples/1k/ag_1k_oral.txt'),
                    'sink-1k-unweighted-unifrac': ('11-packaged/oral/notrim/'
                                                   'all_participants/'
                                                   'all_samples/1k/'
                                                   'unweighted_unifrac_ag_1k_'
                                                   'oral.txt'),
                    'sink-1k-weighted-unifrac': ('11-packaged/oral/notrim/'
                                                 'all_participants/'
                                                 'all_samples/1k/'
                                                 'weighted_unifrac_ag_1k_'
                                                 'oral.txt'),

                    'sink-10k-dir': ('11-packaged/oral/notrim/'
                                     'all_participants/all_samples/10k'),
                    'sink-10k-otu': ('11-packaged/oral/notrim/'
                                     'all_participants/all_samples/10k/'
                                     'ag_10k_oral.biom'),
                    'sink-10k-map': ('11-packaged/oral/notrim/'
                                     'all_participants/all_samples/10k/'
                                     'ag_10k_oral.txt'),
                    'sink-10k-unweighted-unifrac': ('11-packaged/oral/notrim/'
                                                    'all_participants/'
                                                    'all_samples/10k/'
                                                    'unweighted_unifrac_ag_'
                                                    '10k_oral.txt'),
                    'sink-10k-weighted-unifrac': ('11-packaged/oral/notrim/'
                                                  'all_participants/'
                                                  'all_samples/10k/'
                                                  'weighted_unifrac_ag_10k_'
                                                  'oral.txt'),
                    },
                '100nt': {
                    'source-unrare-otu': ('11-packaged/split/100nt/raw/'
                                          'otu_table__BODY_HABITAT_UBERON:'
                                          'oral cavity__.biom'),
                    'source-unrare-map': ('11-packaged/split/100nt/raw/'
                                          'metadata__BODY_HABITAT_UBERON:'
                                          'oral cavity__.txt'),
                    'sink-unrare-dir': ('11-packaged/oral/100nt'),
                    'sink-unrare-otu': ('11-packaged/oral/100nt/ag_oral.'
                                        'biom'),
                    'sink-unrare-map': ('11-packaged/oral/100nt/ag_oral.'
                                        'txt'),

                    'source-1k-otu': ('11-packaged/split/100nt/1k/'
                                      'ag__BODY_HABITAT_UBERON:'
                                      'oral cavity__.biom'),
                    'source-1k-map': ('11-packaged/split/100nt/1k/'
                                      'metadata__BODY_HABITAT_UBERON:'
                                      'oral cavity__.txt'),
                    'source-1k-pd': ('05-alpha/1k/ag-pgp-hmp-gg-100nt/'
                                     'PD_whole_tree.txt'),
                    'source-1k-chao1': ('05-alpha/1k/ag-pgp-hmp-gg-100nt/'
                                        'chao1.txt'),
                    'source-1k-observedotus': ('05-alpha/1k/ag-pgp-hmp-gg-'
                                               '100nt/observed_otus.txt'),
                    'source-1k-shannon': ('05-alpha/1k/ag-pgp-hmp-gg-100nt/'
                                          'shannon.txt'),
                    'source-1k-unweighted-unifrac': ('06-beta/100nt/1k/ag/'
                                                     'unweighted_unifrac_ag-'
                                                     'oral.txt'),
                    'source-1k-weighted-unifrac': ('06-beta/100nt/1k/ag/'
                                                   'weighted_unifrac_ag-oral'
                                                   '.txt'),

                    'source-10k-otu': ('11-packaged/split/100nt/10k/'
                                       'ag__BODY_HABITAT_UBERON:oral '
                                       'cavity__.biom'),
                    'source-10k-map': ('11-packaged/split/100nt/10k/'
                                       'metadata__BODY_HABITAT_UBERON:oral '
                                       'cavity__.txt'),
                    'source-10k-pd': ('05-alpha/10k/ag-pgp-hmp-gg-100nt/'
                                      'PD_whole_tree.txt'),
                    'source-10k-chao1': ('05-alpha/10k/ag-pgp-hmp-gg-100nt/'
                                         'chao1.txt'),
                    'source-10k-observedotus': ('05-alpha/10k/ag-pgp-hmp-gg-'
                                                '100nt/observed_otus.txt'),
                    'source-10k-shannon': ('05-alpha/10k/ag-pgp-hmp-gg-100nt/'
                                           'shannon.txt'),
                    'source-10k-unweighted-unifrac': ('06-beta/100nt/10k/ag/'
                                                      'unweighted_unifrac_ag-'
                                                      'oral.txt'),
                    'source-10k-weighted-unifrac': ('06-beta/100nt/10k/ag/'
                                                    'weighted_unifrac_ag-oral'
                                                    '.txt'),

                    'sink-1k-dir': ('11-packaged/oral/100nt/all_participants/'
                                    'all_samples/1k'),
                    'sink-1k-otu': ('11-packaged/oral/100nt/all_participants/'
                                    'all_samples/1k/ag_1k_oral.biom'),
                    'sink-1k-map': ('11-packaged/oral/100nt/all_participants/'
                                    'all_samples/1k/ag_1k_oral.txt'),
                    'sink-1k-unweighted-unifrac': ('11-packaged/oral/100nt/'
                                                   'all_participants/'
                                                   'all_samples/1k/'
                                                   'unweighted_unifrac_ag_1k_'
                                                   'oral.txt'),
                    'sink-1k-weighted-unifrac': ('11-packaged/oral/100nt/'
                                                 'all_participants/'
                                                 'all_samples/1k/'
                                                 'weighted_unifrac_ag_1k_'
                                                 'oral.txt'),

                    'sink-10k-dir': ('11-packaged/oral/100nt/all_participants/'
                                     'all_samples/10k'),
                    'sink-10k-otu': ('11-packaged/oral/100nt/all_participants/'
                                     'all_samples/10k/ag_10k_oral.biom'),
                    'sink-10k-map': ('11-packaged/oral/100nt/all_participants/'
                                     'all_samples/10k/ag_10k_oral.txt'),
                    'sink-10k-unweighted-unifrac': ('11-packaged/oral/100nt/'
                                                    'all_participants/'
                                                    'all_samples/10k/'
                                                    'unweighted_unifrac_ag_10k'
                                                    '_oral.txt'),
                    'sink-10k-weighted-unifrac': ('11-packaged/oral/100nt/'
                                                  'all_participants/'
                                                  'all_samples/10k/'
                                                  'weighted_unifrac_ag_10k_'
                                                  'oral.txt'),
                    },
                },
            'skin': {
                'notrim': {
                    'source-unrare-otu': ('11-packaged/split/notrim/raw/'
                                          'otu_table__BODY_HABITAT_UBERON:'
                                          'skin__.biom'),
                    'source-unrare-map': ('11-packaged/split/notrim/raw/'
                                          'metadata__BODY_HABITAT_UBERON:'
                                          'skin__.txt'),
                    'sink-unrare-dir': ('11-packaged/skin/notrim'),
                    'sink-unrare-otu': ('11-packaged/skin/notrim/ag_skin.'
                                        'biom'),
                    'sink-unrare-map': ('11-packaged/skin/notrim/ag_skin.'
                                        'txt'),

                    'source-1k-otu': ('11-packaged/split/notrim/1k/'
                                      'ag__BODY_HABITAT_UBERON:skin__.biom'),
                    'source-1k-map': ('11-packaged/split/notrim/1k/'
                                      'metadata__BODY_HABITAT_UBERON:'
                                      'skin__.txt'),
                    'source-1k-pd': '05-alpha/1k/ag-notrim/PD_whole_tree.txt',
                    'source-1k-chao1': '05-alpha/1k/ag-notrim/chao1.txt',
                    'source-1k-observedotus': ('05-alpha/1k/ag-notrim/'
                                               'observed_otus.txt'),
                    'source-1k-shannon': '05-alpha/1k/ag-notrim/shannon.txt',
                    'source-1k-unweighted-unifrac': ('06-beta/notrim/1k/ag/'
                                                     'unweighted_unifrac_ag-'
                                                     'skin.txt'),
                    'source-1k-weighted-unifrac': ('06-beta/notrim/1k/ag/'
                                                   'weighted_unifrac_ag-'
                                                   'skin.txt'),

                    'source-10k-otu': ('11-packaged/split/notrim/10k/'
                                       'ag__BODY_HABITAT_UBERON:skin__.biom'),
                    'source-10k-map': ('11-packaged/split/notrim/10k/'
                                       'metadata__BODY_HABITAT_UBERON:'
                                       'skin__.txt'),
                    'source-10k-pd': ('05-alpha/10k/ag-notrim/PD_whole_'
                                      'tree.txt'),
                    'source-10k-chao1': '05-alpha/10k/ag-notrim/chao1.txt',
                    'source-10k-observedotus': ('05-alpha/10k/ag-notrim/'
                                                'observed_otus.txt'),
                    'source-10k-shannon': '05-alpha/10k/ag-notrim/shannon.txt',
                    'source-10k-unweighted-unifrac': ('06-beta/notrim/10k/ag/'
                                                      'unweighted_unifrac_ag-'
                                                      'skin.txt'),
                    'source-10k-weighted-unifrac': ('06-beta/notrim/10k/ag/'
                                                    'weighted_unifrac_ag-'
                                                    'skin.txt'),

                    'sink-1k-dir': ('11-packaged/skin/notrim/all_participants/'
                                    'all_samples/1k'),
                    'sink-1k-otu': ('11-packaged/skin/notrim/all_participants/'
                                    'all_samples/1k/ag_1k_skin.biom'),
                    'sink-1k-map': ('11-packaged/skin/notrim/all_participants/'
                                    'all_samples/1k/ag_1k_skin.txt'),
                    'sink-1k-unweighted-unifrac': ('11-packaged/skin/notrim/'
                                                   'all_participants/'
                                                   'all_samples/1k/'
                                                   'unweighted_unifrac_ag_1k_'
                                                   'skin.txt'),
                    'sink-1k-weighted-unifrac': ('11-packaged/skin/notrim/'
                                                 'all_participants/'
                                                 'all_samples/1k/'
                                                 'weighted_unifrac_ag_1k_'
                                                 'skin.txt'),

                    'sink-10k-dir': ('11-packaged/skin/notrim/'
                                     'all_participants/all_samples/10k'),
                    'sink-10k-otu': ('11-packaged/skin/notrim/'
                                     'all_participants/all_samples/10k/'
                                     'ag_10k_skin.biom'),
                    'sink-10k-map': ('11-packaged/skin/notrim/'
                                     'all_participants/all_samples/10k/'
                                     'ag_10k_skin.txt'),
                    'sink-10k-unweighted-unifrac': ('11-packaged/skin/notrim/'
                                                    'all_participants/'
                                                    'all_samples/10k/'
                                                    'unweighted_unifrac_ag_'
                                                    '10k_skin.txt'),
                    'sink-10k-weighted-unifrac': ('11-packaged/skin/notrim/'
                                                  'all_participants/'
                                                  'all_samples/10k/'
                                                  'weighted_unifrac_ag_10k_'
                                                  'skin.txt'),
                    },
                '100nt': {
                    'source-unrare-otu': ('11-packaged/split/100nt/raw/'
                                          'otu_table__BODY_HABITAT_UBERON:skin'
                                          '__.biom'),
                    'source-unrare-map': ('11-packaged/split/100nt/raw/'
                                          'metadata__BODY_HABITAT_UBERON:'
                                          'skin__.txt'),
                    'sink-unrare-dir': ('11-packaged/skin/100nt'),
                    'sink-unrare-otu': ('11-packaged/skin/100nt/ag_skin.'
                                        'biom'),
                    'sink-unrare-map': ('11-packaged/skin/100nt/ag_skin.'
                                        'txt'),

                    'source-1k-otu': ('11-packaged/split/100nt/1k/'
                                      'ag__BODY_HABITAT_UBERON:skin__.biom'),
                    'source-1k-map': ('11-packaged/split/100nt/1k/'
                                      'metadata__BODY_HABITAT_UBERON:'
                                      'skin__.txt'),
                    'source-1k-pd': ('05-alpha/1k/ag-pgp-hmp-gg-100nt/'
                                     'PD_whole_tree.txt'),
                    'source-1k-chao1': ('05-alpha/1k/ag-pgp-hmp-gg-100nt/'
                                        'chao1.txt'),
                    'source-1k-observedotus': ('05-alpha/1k/'
                                               'ag-pgp-hmp-gg-100nt/'
                                               'observed_otus.txt'),
                    'source-1k-shannon': ('05-alpha/1k/ag-pgp-hmp-gg-100nt/'
                                          'shannon.txt'),
                    'source-1k-unweighted-unifrac': ('06-beta/100nt/1k/ag/'
                                                     'unweighted_unifrac_ag-'
                                                     'skin.txt'),
                    'source-1k-weighted-unifrac': ('06-beta/100nt/1k/ag/'
                                                   'weighted_unifrac_ag-'
                                                   'skin.txt'),

                    'source-10k-otu': ('11-packaged/split/100nt/10k/'
                                       'ag__BODY_HABITAT_UBERON:skin__.biom'),
                    'source-10k-map': ('11-packaged/split/100nt/10k/'
                                       'metadata__BODY_HABITAT_UBERON:'
                                       'skin__.txt'),
                    'source-10k-pd': ('05-alpha/10k/ag-pgp-hmp-gg-100nt/'
                                      'PD_whole_tree.txt'),
                    'source-10k-chao1': ('05-alpha/10k/ag-pgp-hmp-gg-100nt/'
                                         'chao1.txt'),
                    'source-10k-observedotus': ('05-alpha/10k/'
                                                'ag-pgp-hmp-gg-100nt/'
                                                'observed_otus.txt'),
                    'source-10k-shannon': ('05-alpha/10k/'
                                           'ag-pgp-hmp-gg-100nt/shannon.txt'),
                    'source-10k-unweighted-unifrac': ('06-beta/100nt/10k/ag/'
                                                      'unweighted_unifrac_ag-'
                                                      'skin.txt'),
                    'source-10k-weighted-unifrac': ('06-beta/100nt/10k/ag/'
                                                    'weighted_unifrac_ag-'
                                                    'skin.txt'),

                    'sink-1k-dir': ('11-packaged/skin/100nt/all_participants/'
                                    'all_samples/1k'),
                    'sink-1k-otu': ('11-packaged/skin/100nt/all_participants/'
                                    'all_samples/1k/ag_1k_skin.biom'),
                    'sink-1k-map': ('11-packaged/skin/100nt/all_participants/'
                                    'all_samples/1k/ag_1k_skin.txt'),
                    'sink-1k-unweighted-unifrac': ('11-packaged/skin/100nt/'
                                                   'all_participants/'
                                                   'all_samples/1k/'
                                                   'unweighted_unifrac_ag_'
                                                   '1k_skin.txt'),
                    'sink-1k-weighted-unifrac': ('11-packaged/skin/100nt/'
                                                 'all_participants/'
                                                 'all_samples/1k/weighted_'
                                                 'unifrac_ag_1k_skin.txt'),

                    'sink-10k-dir': ('11-packaged/skin/100nt/all_participants/'
                                     'all_samples/10k'),
                    'sink-10k-otu': ('11-packaged/skin/100nt/all_participants/'
                                     'all_samples/10k/ag_10k_skin.biom'),
                    'sink-10k-map': ('11-packaged/skin/100nt/all_participants/'
                                     'all_samples/10k/ag_10k_skin.txt'),
                    'sink-10k-unweighted-unifrac': ('11-packaged/skin/100nt/'
                                                    'all_participants/'
                                                    'all_samples/10k/'
                                                    'unweighted_unifrac_ag_'
                                                    '10k_skin.txt'),
                    'sink-10k-weighted-unifrac': ('11-packaged/skin/100nt/'
                                                  'all_participants/'
                                                  'all_samples/10k/'
                                                  'weighted_unifrac_ag_10k_'
                                                  'skin.txt'),
                    },
                },
            },

        'all_participants_one_sample': {
            'fecal': {
                'notrim': {

                    'source-1k-otu': ('11-packaged/fecal/notrim/'
                                      'all_participants/all_samples/1k/'
                                      'ag_1k_fecal.biom'),
                    'source-1k-map': ('11-packaged/fecal/notrim/'
                                      'all_participants/all_samples/1k/'
                                      'ag_1k_fecal.txt'),
                    'source-1k-unweighted-unifrac': ('11-packaged/fecal/'
                                                     'notrim/all_participants/'
                                                     'all_samples/1k/'
                                                     'unweighted_unifrac_ag_1k'
                                                     '_fecal.txt'),
                    'source-1k-weighted-unifrac': ('11-packaged/fecal/'
                                                   'notrim/all_participants/'
                                                   'all_samples/1k/'
                                                   'weighted_unifrac_ag_1k_'
                                                   'fecal.txt'),

                    'source-10k-otu': ('11-packaged/fecal/notrim/'
                                       'all_participants/all_samples/10k/'
                                       'ag_10k_fecal.biom'),
                    'source-10k-map': ('11-packaged/fecal/notrim/'
                                       'all_participants/all_samples/10k/'
                                       'ag_10k_fecal.txt'),
                    'source-10k-unweighted-unifrac': ('11-packaged/fecal/'
                                                      'notrim/'
                                                      'all_participants/'
                                                      'all_samples/10k/'
                                                      'unweighted_unifrac_ag_'
                                                      '10k_fecal.txt'),
                    'source-10k-weighted-unifrac': ('11-packaged/fecal/'
                                                    'notrim/all_participants/'
                                                    'all_samples/10k/'
                                                    'weighted_unifrac_ag_'
                                                    '10k_fecal.txt'),

                    'sink-1k-dir': ('11-packaged/fecal/notrim/'
                                    'all_participants/one_sample/1k'),
                    'sink-1k-otu': ('11-packaged/fecal/notrim/'
                                    'all_participants/one_sample/1k/'
                                    'ag_1k_fecal.biom'),
                    'sink-1k-map': ('11-packaged/fecal/notrim/'
                                    'all_participants/one_sample/1k/'
                                    'ag_1k_fecal.txt'),
                    'sink-1k-unweighted-unifrac': ('11-packaged/fecal/notrim/'
                                                   'all_participants/'
                                                   'one_sample/1k/'
                                                   'unweighted_unifrac_ag_1k_'
                                                   'fecal.txt'),
                    'sink-1k-weighted-unifrac': ('11-packaged/fecal/notrim/'
                                                 'all_participants/one_sample/'
                                                 '1k/weighted_unifrac_ag_1k_'
                                                 'fecal.txt'),

                    'sink-10k-dir': ('11-packaged/fecal/notrim/'
                                     'all_participants/one_sample/10k'),
                    'sink-10k-otu': ('11-packaged/fecal/notrim/'
                                     'all_participants/one_sample/10k/'
                                     'ag_10k_fecal.biom'),
                    'sink-10k-map': ('11-packaged/fecal/notrim/'
                                     'all_participants/one_sample/10k/'
                                     'ag_10k_fecal.txt'),
                    'sink-10k-unweighted-unifrac': ('11-packaged/fecal/notrim/'
                                                    'all_participants/'
                                                    'one_sample/10k/'
                                                    'unweighted_unifrac_'
                                                    'ag_10k_fecal.txt'),
                    'sink-10k-weighted-unifrac': ('11-packaged/fecal/notrim/'
                                                  'all_participants/'
                                                  'one_sample/10k/weighted_'
                                                  'unifrac_ag_10k_fecal.txt'),
                    },
                '100nt': {

                    'source-1k-otu': ('11-packaged/fecal/100nt/'
                                      'all_participants/all_samples/1k/'
                                      'ag_1k_fecal.biom'),
                    'source-1k-map': ('11-packaged/fecal/100nt/'
                                      'all_participants/all_samples/1k/'
                                      'ag_1k_fecal.txt'),
                    'source-1k-unweighted-unifrac': ('11-packaged/fecal/100nt/'
                                                     'all_participants/'
                                                     'all_samples/1k/'
                                                     'unweighted_unifrac_ag_'
                                                     '1k_fecal.txt'),
                    'source-1k-weighted-unifrac': ('11-packaged/fecal/100nt/'
                                                   'all_participants/'
                                                   'all_samples/1k/'
                                                   'weighted_unifrac_ag_1k_'
                                                   'fecal.txt'),

                    'source-10k-otu': ('11-packaged/fecal/100nt/'
                                       'all_participants/all_samples/10k/'
                                       'ag_10k_fecal.biom'),
                    'source-10k-map': ('11-packaged/fecal/100nt/'
                                       'all_participants/all_samples/10k/'
                                       'ag_10k_fecal.txt'),
                    'source-10k-unweighted-unifrac': ('11-packaged/fecal/'
                                                      '100nt/all_participants/'
                                                      'all_samples/10k/'
                                                      'unweighted_unifrac_ag_'
                                                      '10k_fecal.txt'),
                    'source-10k-weighted-unifrac': ('11-packaged/fecal/100nt/'
                                                    'all_participants/'
                                                    'all_samples/10k/'
                                                    'weighted_unifrac_'
                                                    'ag_10k_fecal.txt'),

                    'sink-1k-dir': ('11-packaged/fecal/100nt/'
                                    'all_participants/one_sample/1k'),
                    'sink-1k-otu': ('11-packaged/fecal/100nt/'
                                    'all_participants/one_sample/1k/'
                                    'ag_1k_fecal.biom'),
                    'sink-1k-map': ('11-packaged/fecal/100nt/'
                                    'all_participants/one_sample/1k/'
                                    'ag_1k_fecal.txt'),
                    'sink-1k-unweighted-unifrac': ('11-packaged/fecal/100nt/'
                                                   'all_participants/'
                                                   'one_sample/1k/'
                                                   'unweighted_unifrac_ag_1k_'
                                                   'fecal.txt'),
                    'sink-1k-weighted-unifrac': ('11-packaged/fecal/100nt/'
                                                 'all_participants/one_sample/'
                                                 '1k/weighted_unifrac_ag_1k_'
                                                 'fecal.txt'),

                    'sink-10k-dir': ('11-packaged/fecal/100nt/'
                                     'all_participants/one_sample/10k'),
                    'sink-10k-otu': ('11-packaged/fecal/100nt/'
                                     'all_participants/one_sample/10k/'
                                     'ag_10k_fecal.biom'),
                    'sink-10k-map': ('11-packaged/fecal/100nt/'
                                     'all_participants/one_sample/10k/'
                                     'ag_10k_fecal.txt'),
                    'sink-10k-unweighted-unifrac': ('11-packaged/fecal/100nt/'
                                                    'all_participants/'
                                                    'one_sample/10k/'
                                                    'unweighted_unifrac_'
                                                    'ag_10k_fecal.txt'),
                    'sink-10k-weighted-unifrac': ('11-packaged/fecal/100nt/'
                                                  'all_participants/'
                                                  'one_sample/10k/weighted_'
                                                  'unifrac_ag_10k_fecal.txt'),
                    },
                },
            'oral': {
                'notrim': {

                    'source-1k-otu': ('11-packaged/oral/notrim/'
                                      'all_participants/all_samples/1k/'
                                      'ag_1k_oral.biom'),
                    'source-1k-map': ('11-packaged/oral/notrim/'
                                      'all_participants/all_samples/1k/'
                                      'ag_1k_oral.txt'),
                    'source-1k-unweighted-unifrac': ('11-packaged/oral/notrim/'
                                                     'all_participants/'
                                                     'all_samples/1k/'
                                                     'unweighted_unifrac_ag_'
                                                     '1k_oral.txt'),
                    'source-1k-weighted-unifrac': ('11-packaged/oral/notrim/'
                                                   'all_participants/'
                                                   'all_samples/1k/'
                                                   'weighted_unifrac_ag_1k_'
                                                   'oral.txt'),

                    'source-10k-otu': ('11-packaged/oral/notrim/'
                                       'all_participants/'
                                       'all_samples/10k/'
                                       'ag_10k_oral.biom'),
                    'source-10k-map': ('11-packaged/oral/notrim/'
                                       'all_participants/'
                                       'all_samples/10k/'
                                       'ag_10k_oral.txt'),
                    'source-10k-unweighted-unifrac': ('11-packaged/oral/'
                                                      'notrim/'
                                                      'all_participants/'
                                                      'all_samples/10k/'
                                                      'unweighted_unifrac_'
                                                      'ag_10k_oral.txt'),
                    'source-10k-weighted-unifrac': ('11-packaged/oral/'
                                                    'notrim/all_participants/'
                                                    'all_samples/10k/weighted'
                                                    '_unifrac_ag_10k_'
                                                    'oral.txt'),

                    'sink-1k-dir': ('11-packaged/oral/notrim/all_participants/'
                                    'one_sample/1k'),
                    'sink-1k-otu': ('11-packaged/oral/notrim/all_participants/'
                                    'one_sample/1k/ag_1k_oral.biom'),
                    'sink-1k-map': ('11-packaged/oral/notrim/all_participants/'
                                    'one_sample/1k/ag_1k_oral.txt'),
                    'sink-1k-unweighted-unifrac': ('11-packaged/oral/notrim/'
                                                   'all_participants/'
                                                   'one_sample/1k/'
                                                   'unweighted_unifrac_ag_1k_'
                                                   'oral.txt'),
                    'sink-1k-weighted-unifrac': ('11-packaged/oral/notrim/'
                                                 'all_participants/one_sample/'
                                                 '1k/weighted_unifrac_ag_1k_'
                                                 'oral.txt'),

                    'sink-10k-dir': ('11-packaged/oral/notrim/'
                                     'all_participants/one_sample/10k'),
                    'sink-10k-otu': ('11-packaged/oral/notrim/'
                                     'all_participants/one_sample/10k/'
                                     'ag_10k_oral.biom'),
                    'sink-10k-map': ('11-packaged/oral/notrim/'
                                     'all_participants/one_sample/10k/'
                                     'ag_10k_oral.txt'),
                    'sink-10k-unweighted-unifrac': ('11-packaged/oral/notrim/'
                                                    'all_participants/'
                                                    'one_sample/10k/'
                                                    'unweighted_unifrac_ag_'
                                                    '10k_oral.txt'),
                    'sink-10k-weighted-unifrac': ('11-packaged/oral/notrim/'
                                                  'all_participants/'
                                                  'one_sample/10k/'
                                                  'weighted_unifrac_ag_10k_'
                                                  'oral.txt'),
                    },
                '100nt': {

                    'source-1k-otu': ('11-packaged/oral/100nt/'
                                      'all_participants/all_samples/1k/'
                                      'ag_1k_oral.biom'),
                    'source-1k-map': ('11-packaged/oral/100nt/'
                                      'all_participants/all_samples/1k/'
                                      'ag_1k_oral.txt'),
                    'source-1k-unweighted-unifrac': ('11-packaged/oral/100nt/'
                                                     'all_participants/'
                                                     'all_samples/1k/'
                                                     'unweighted_unifrac_ag_'
                                                     '1k_oral.txt'),
                    'source-1k-weighted-unifrac': ('11-packaged/oral/100nt/'
                                                   'all_participants/'
                                                   'all_samples/1k/'
                                                   'weighted_unifrac_ag_'
                                                   '1k_oral.txt'),

                    'source-10k-otu': ('11-packaged/oral/100nt/'
                                       'all_participants/all_samples/10k/'
                                       'ag_10k_oral.biom'),
                    'source-10k-map': ('11-packaged/oral/100nt/'
                                       'all_participants/all_samples/10k/'
                                       'ag_10k_oral.txt'),
                    'source-10k-unweighted-unifrac': ('11-packaged/oral/100nt/'
                                                      'all_participants/'
                                                      'all_samples/10k/'
                                                      'unweighted_unifrac_'
                                                      'ag_10k_oral.txt'),
                    'source-10k-weighted-unifrac': ('11-packaged/oral/100nt/'
                                                    'all_participants/'
                                                    'all_samples/10k/'
                                                    'weighted_unifrac_'
                                                    'ag_10k_oral.txt'),

                    'sink-1k-dir': ('11-packaged/oral/100nt/all_participants/'
                                    'one_sample/1k'),
                    'sink-1k-otu': ('11-packaged/oral/100nt/all_participants/'
                                    'one_sample/1k/ag_1k_oral.biom'),
                    'sink-1k-map': ('11-packaged/oral/100nt/all_participants/'
                                    'one_sample/1k/ag_1k_oral.txt'),
                    'sink-1k-unweighted-unifrac': ('11-packaged/oral/100nt/'
                                                   'all_participants/'
                                                   'one_sample/1k/unweighted_'
                                                   'unifrac_ag_1k_oral.txt'),
                    'sink-1k-weighted-unifrac': ('11-packaged/oral/100nt/'
                                                 'all_participants/one_sample/'
                                                 '1k/weighted_unifrac_'
                                                 'ag_1k_oral.txt'),

                    'sink-10k-dir': ('11-packaged/oral/100nt/all_participants/'
                                     'one_sample/10k'),
                    'sink-10k-otu': ('11-packaged/oral/100nt/all_participants/'
                                     'one_sample/10k/ag_10k_oral.biom'),
                    'sink-10k-map': ('11-packaged/oral/100nt/all_participants/'
                                     'one_sample/10k/ag_10k_oral.txt'),
                    'sink-10k-unweighted-unifrac': ('11-packaged/oral/100nt/'
                                                    'all_participants/'
                                                    'one_sample/10k/'
                                                    'unweighted_unifrac_'
                                                    'ag_10k_oral.txt'),
                    'sink-10k-weighted-unifrac': ('11-packaged/oral/100nt/'
                                                  'all_participants/'
                                                  'one_sample/10k/weighted_'
                                                  'unifrac_ag_10k_oral.txt'),
                    },
                },
            'skin': {
                'notrim': {
                    'source-1k-otu': ('11-packaged/skin/notrim/'
                                      'all_participants/all_samples/1k/'
                                      'ag_1k_skin.biom'),
                    'source-1k-map': ('11-packaged/skin/notrim/'
                                      'all_participants/all_samples/1k/'
                                      'ag_1k_skin.txt'),
                    'source-1k-unweighted-unifrac': ('11-packaged/skin/notrim/'
                                                     'all_participants/'
                                                     'all_samples/1k/'
                                                     'unweighted_unifrac_ag_'
                                                     '1k_skin.txt'),
                    'source-1k-weighted-unifrac': ('11-packaged/skin/notrim/'
                                                   'all_participants/'
                                                   'all_samples/1k/'
                                                   'weighted_unifrac_ag_1k_'
                                                   'skin.txt'),

                    'source-10k-otu': ('11-packaged/skin/notrim/'
                                       'all_participants/all_samples/10k/'
                                       'ag_10k_skin.biom'),
                    'source-10k-map': ('11-packaged/skin/notrim/'
                                       'all_participants/all_samples/10k/'
                                       'ag_10k_skin.txt'),
                    'source-10k-unweighted-unifrac': ('11-packaged/skin/'
                                                      'notrim/'
                                                      'all_participants/'
                                                      'all_samples/10k/'
                                                      'unweighted_unifrac_ag_'
                                                      '10k_skin.txt'),
                    'source-10k-weighted-unifrac': ('11-packaged/skin/notrim/'
                                                    'all_participants/'
                                                    'all_samples/10k/'
                                                    'weighted_unifrac_ag_10k_'
                                                    'skin.txt'),

                    'sink-1k-dir': ('11-packaged/skin/notrim/all_participants/'
                                    'one_sample/1k'),
                    'sink-1k-otu': ('11-packaged/skin/notrim/all_participants/'
                                    'one_sample/1k/ag_1k_skin.biom'),
                    'sink-1k-map': ('11-packaged/skin/notrim/all_participants/'
                                    'one_sample/1k/ag_1k_skin.txt'),
                    'sink-1k-unweighted-unifrac': ('11-packaged/skin/notrim/'
                                                   'all_participants/'
                                                   'one_sample/1k/'
                                                   'unweighted_unifrac_ag_1k_'
                                                   'skin.txt'),
                    'sink-1k-weighted-unifrac': ('11-packaged/skin/notrim/'
                                                 'all_participants/one_sample/'
                                                 '1k/weighted_unifrac_ag_'
                                                 '1k_skin.txt'),

                    'sink-10k-dir': ('11-packaged/skin/notrim/'
                                     'all_participants/one_sample/10k'),
                    'sink-10k-otu': ('11-packaged/skin/notrim/'
                                     'all_participants/one_sample/10k/'
                                     'ag_10k_skin.biom'),
                    'sink-10k-map': ('11-packaged/skin/notrim/'
                                     'all_participants/one_sample/10k/'
                                     'ag_10k_skin.txt'),
                    'sink-10k-unweighted-unifrac': ('11-packaged/skin/notrim/'
                                                    'all_participants/'
                                                    'one_sample/10k/'
                                                    'unweighted_unifrac_ag_'
                                                    '10k_skin.txt'),
                    'sink-10k-weighted-unifrac': ('11-packaged/skin/notrim/'
                                                  'all_participants/'
                                                  'one_sample/10k/'
                                                  'weighted_unifrac_ag_10k_'
                                                  'skin.txt'),
                    },
                '100nt': {
                    'source-1k-otu': ('11-packaged/skin/100nt/'
                                      'all_participants/all_samples/1k/'
                                      'ag_1k_skin.biom'),
                    'source-1k-map': ('11-packaged/skin/100nt/'
                                      'all_participants/all_samples/1k/'
                                      'ag_1k_skin.txt'),
                    'source-1k-unweighted-unifrac': ('11-packaged/skin/100nt/'
                                                     'all_participants/'
                                                     'all_samples/1k/'
                                                     'unweighted_unifrac_'
                                                     'ag_1k_skin.txt'),
                    'source-1k-weighted-unifrac': ('11-packaged/skin/100nt/'
                                                   'all_participants/'
                                                   'all_samples/1k/'
                                                   'weighted_unifrac_'
                                                   'ag_1k_skin.txt'),

                    'source-10k-otu': ('11-packaged/skin/100nt/'
                                       'all_participants/all_samples/10k/'
                                       'ag_10k_skin.biom'),
                    'source-10k-map': ('11-packaged/skin/100nt/'
                                       'all_participants/all_samples/10k/'
                                       'ag_10k_skin.txt'),
                    'source-10k-unweighted-unifrac': ('11-packaged/skin/100nt/'
                                                      'all_participants/'
                                                      'all_samples/10k/'
                                                      'unweighted_unifrac_ag_'
                                                      '10k_skin.txt'),
                    'source-10k-weighted-unifrac': ('11-packaged/skin/100nt/'
                                                    'all_participants/'
                                                    'all_samples/10k/'
                                                    'weighted_unifrac_ag_'
                                                    '10k_skin.txt'),

                    'sink-1k-dir': ('11-packaged/skin/100nt/all_participants/'
                                    'one_sample/1k'),
                    'sink-1k-otu': ('11-packaged/skin/100nt/all_participants/'
                                    'one_sample/1k/ag_1k_skin.biom'),
                    'sink-1k-map': ('11-packaged/skin/100nt/all_participants/'
                                    'one_sample/1k/ag_1k_skin.txt'),
                    'sink-1k-unweighted-unifrac': ('11-packaged/skin/100nt/'
                                                   'all_participants/'
                                                   'one_sample/1k/'
                                                   'unweighted_unifrac_ag_'
                                                   '1k_skin.txt'),
                    'sink-1k-weighted-unifrac': ('11-packaged/skin/100nt/'
                                                 'all_participants/'
                                                 'one_sample/1k/'
                                                 'weighted_unifrac_ag_'
                                                 '1k_skin.txt'),

                    'sink-10k-dir': ('11-packaged/skin/100nt/all_participants/'
                                     'one_sample/10k'),
                    'sink-10k-otu': ('11-packaged/skin/100nt/all_participants/'
                                     'one_sample/10k/ag_10k_skin.biom'),
                    'sink-10k-map': ('11-packaged/skin/100nt/all_participants/'
                                     'one_sample/10k/ag_10k_skin.txt'),
                    'sink-10k-unweighted-unifrac': ('11-packaged/skin/100nt/'
                                                    'all_participants/'
                                                    'one_sample/10k/'
                                                    'unweighted_unifrac_'
                                                    'ag_10k_skin.txt'),
                    'sink-10k-weighted-unifrac': ('11-packaged/skin/100nt/'
                                                  'all_participants/'
                                                  'one_sample/10k/'
                                                  'weighted_unifrac_ag_'
                                                  '10k_skin.txt'),
                    },
                },
            },

        'sub_participants_all_samples': {
            'fecal': {
                'notrim': {
                    'source-1k-otu': ('11-packaged/fecal/notrim/'
                                      'all_participants/all_samples/1k/'
                                      'ag_1k_fecal.biom'),
                    'source-1k-map': ('11-packaged/fecal/notrim/'
                                      'all_participants/all_samples/1k/'
                                      'ag_1k_fecal.txt'),
                    'source-1k-unweighted-unifrac': ('11-packaged/fecal/'
                                                     'notrim/all_participants/'
                                                     'all_samples/1k/'
                                                     'unweighted_unifrac_ag_'
                                                     '1k_fecal.txt'),
                    'source-1k-weighted-unifrac': ('11-packaged/fecal/notrim/'
                                                   'all_participants/'
                                                   'all_samples/1k/'
                                                   'weighted_unifrac_ag_1k_'
                                                   'fecal.txt'),

                    'source-10k-otu': ('11-packaged/fecal/notrim/'
                                       'all_participants/all_samples/10k/'
                                       'ag_10k_fecal.biom'),
                    'source-10k-map': ('11-packaged/fecal/notrim/'
                                       'all_participants/all_samples/10k/'
                                       'ag_10k_fecal.txt'),
                    'source-10k-unweighted-unifrac': ('11-packaged/fecal/'
                                                      'notrim/'
                                                      'all_participants/'
                                                      'all_samples/10k/'
                                                      'unweighted_unifrac_'
                                                      'ag_10k_fecal.txt'),
                    'source-10k-weighted-unifrac': ('11-packaged/fecal/notrim/'
                                                    'all_participants/'
                                                    'all_samples/10k/'
                                                    'weighted_unifrac_ag_'
                                                    '10k_fecal.txt'),

                    'sink-1k-dir': ('11-packaged/fecal/notrim/'
                                    'sub_participants/all_samples/1k'),
                    'sink-1k-otu': ('11-packaged/fecal/notrim/'
                                    'sub_participants/all_samples/1k/'
                                    'ag_1k_fecal.biom'),
                    'sink-1k-map': ('11-packaged/fecal/notrim/'
                                    'sub_participants/all_samples/1k/'
                                    'ag_1k_fecal.txt'),
                    'sink-1k-unweighted-unifrac': ('11-packaged/fecal/notrim/'
                                                   'sub_participants/'
                                                   'all_samples/1k/'
                                                   'unweighted_unifrac_ag_1k_'
                                                   'fecal.txt'),
                    'sink-1k-weighted-unifrac': ('11-packaged/fecal/notrim/'
                                                 'sub_participants/'
                                                 'all_samples/1k/'
                                                 'weighted_unifrac_ag_1k_'
                                                 'fecal.txt'),

                    'sink-10k-dir': ('11-packaged/fecal/notrim/'
                                     'sub_participants/all_samples/10k'),
                    'sink-10k-otu': ('11-packaged/fecal/notrim/'
                                     'sub_participants/all_samples/10k/'
                                     'ag_10k_fecal.biom'),
                    'sink-10k-map': ('11-packaged/fecal/notrim/'
                                     'sub_participants/all_samples/10k/'
                                     'ag_10k_fecal.txt'),
                    'sink-10k-unweighted-unifrac': ('11-packaged/fecal/'
                                                    'notrim/sub_participants/'
                                                    'all_samples/10k/'
                                                    'unweighted_unifrac_ag_'
                                                    '10k_fecal.txt'),
                    'sink-10k-weighted-unifrac': ('11-packaged/fecal/notrim/'
                                                  'sub_participants/'
                                                  'all_samples/10k/'
                                                  'weighted_unifrac_ag_10k_'
                                                  'fecal.txt'),
                    },
                '100nt': {
                    'source-1k-otu': ('11-packaged/fecal/100nt/'
                                      'all_participants/all_samples/1k/'
                                      'ag_1k_fecal.biom'),
                    'source-1k-map': ('11-packaged/fecal/100nt/'
                                      'all_participants/all_samples/1k/'
                                      'ag_1k_fecal.txt'),
                    'source-1k-unweighted-unifrac': ('11-packaged/fecal/100nt/'
                                                     'all_participants/'
                                                     'all_samples/1k/'
                                                     'unweighted_unifrac_ag_'
                                                     '1k_fecal.txt'),
                    'source-1k-weighted-unifrac': ('11-packaged/fecal/100nt/'
                                                   'all_participants/'
                                                   'all_samples/1k/'
                                                   'weighted_unifrac_ag_1k_'
                                                   'fecal.txt'),

                    'source-10k-otu': ('11-packaged/fecal/100nt/'
                                       'all_participants/all_samples/10k/'
                                       'ag_10k_fecal.biom'),
                    'source-10k-map': ('11-packaged/fecal/100nt/'
                                       'all_participants/all_samples/10k/'
                                       'ag_10k_fecal.txt'),
                    'source-10k-unweighted-unifrac': ('11-packaged/fecal/'
                                                      '100nt/all_participants/'
                                                      'all_samples/10k/'
                                                      'unweighted_unifrac_ag_'
                                                      '10k_fecal.txt'),
                    'source-10k-weighted-unifrac': ('11-packaged/fecal/'
                                                    '100nt/all_participants/'
                                                    'all_samples/10k/'
                                                    'weighted_unifrac_ag_10k_'
                                                    'fecal.txt'),

                    'sink-1k-dir': ('11-packaged/fecal/100nt/sub_participants/'
                                    'all_samples/1k'),
                    'sink-1k-otu': ('11-packaged/fecal/100nt/sub_participants/'
                                    'all_samples/1k/ag_1k_fecal.biom'),
                    'sink-1k-map': ('11-packaged/fecal/100nt/sub_participants/'
                                    'all_samples/1k/ag_1k_fecal.txt'),
                    'sink-1k-unweighted-unifrac': ('11-packaged/fecal/100nt/'
                                                   'sub_participants/'
                                                   'all_samples/1k/'
                                                   'unweighted_unifrac_ag_1k_'
                                                   'fecal.txt'),
                    'sink-1k-weighted-unifrac': ('11-packaged/fecal/100nt/'
                                                 'sub_participants/'
                                                 'all_samples/1k/'
                                                 'weighted_unifrac_ag_1k_'
                                                 'fecal.txt'),

                    'sink-10k-dir': ('11-packaged/fecal/100nt/'
                                     'sub_participants/all_samples/10k'),
                    'sink-10k-otu': ('11-packaged/fecal/100nt/'
                                     'sub_participants/all_samples/10k/'
                                     'ag_10k_fecal.biom'),
                    'sink-10k-map': ('11-packaged/fecal/100nt/'
                                     'sub_participants/all_samples/10k/'
                                     'ag_10k_fecal.txt'),
                    'sink-10k-unweighted-unifrac': ('11-packaged/fecal/100nt/'
                                                    'sub_participants/'
                                                    'all_samples/10k/'
                                                    'unweighted_unifrac_ag_'
                                                    '10k_fecal.txt'),
                    'sink-10k-weighted-unifrac': ('11-packaged/fecal/100nt/'
                                                  'sub_participants/'
                                                  'all_samples/10k/'
                                                  'weighted_unifrac_ag_10k_'
                                                  'fecal.txt'),
                    },
                },
            },

        'sub_participants_one_sample': {
            'fecal': {
                'notrim': {

                    'source-1k-otu': ('11-packaged/fecal/notrim/'
                                      'all_participants/one_sample/1k/'
                                      'ag_1k_fecal.biom'),
                    'source-1k-map': ('11-packaged/fecal/notrim/'
                                      'all_participants/one_sample/1k/'
                                      'ag_1k_fecal.txt'),
                    'source-1k-unweighted-unifrac': ('11-packaged/fecal/'
                                                     'notrim/all_participants/'
                                                     'one_sample/1k/'
                                                     'unweighted_unifrac_ag_'
                                                     '1k_fecal.txt'),
                    'source-1k-weighted-unifrac': ('11-packaged/fecal/'
                                                   'notrim/all_participants/'
                                                   'one_sample/1k/weighted_'
                                                   'unifrac_ag_1k_fecal.txt'),

                    'source-10k-otu': ('11-packaged/fecal/notrim/'
                                       'all_participants/one_sample/10k/'
                                       'ag_10k_fecal.biom'),
                    'source-10k-map': ('11-packaged/fecal/notrim/'
                                       'all_participants/one_sample/10k/'
                                       'ag_10k_fecal.txt'),
                    'source-10k-unweighted-unifrac': ('11-packaged/fecal/'
                                                      'notrim/'
                                                      'all_participants/'
                                                      'one_sample/10k/'
                                                      'unweighted_unifrac_ag_'
                                                      '10k_fecal.txt'),
                    'source-10k-weighted-unifrac': ('11-packaged/fecal/'
                                                    'notrim/all_participants/'
                                                    'one_sample/10k/'
                                                    'weighted_unifrac_ag_10k_'
                                                    'fecal.txt'),

                    'sink-1k-dir': ('11-packaged/fecal/notrim/'
                                    'sub_participants/one_sample/1k'),
                    'sink-1k-otu': ('11-packaged/fecal/notrim/'
                                    'sub_participants/one_sample/1k/'
                                    'ag_1k_fecal.biom'),
                    'sink-1k-map': ('11-packaged/fecal/notrim/'
                                    'sub_participants/one_sample/1k/'
                                    'ag_1k_fecal.txt'),
                    'sink-1k-unweighted-unifrac': ('11-packaged/fecal/'
                                                   'notrim/sub_participants/'
                                                   'one_sample/1k/'
                                                   'unweighted_unifrac_ag_1k_'
                                                   'fecal.txt'),
                    'sink-1k-weighted-unifrac': ('11-packaged/fecal/notrim/'
                                                 'sub_participants/one_sample/'
                                                 '1k/weighted_unifrac_ag_1k_'
                                                 'fecal.txt'),

                    'sink-10k-dir': ('11-packaged/fecal/notrim/'
                                     'sub_participants/one_sample/10k'),
                    'sink-10k-otu': ('11-packaged/fecal/notrim/'
                                     'sub_participants/one_sample/10k/'
                                     'ag_10k_fecal.biom'),
                    'sink-10k-map': ('11-packaged/fecal/notrim/'
                                     'sub_participants/one_sample/10k/'
                                     'ag_10k_fecal.txt'),
                    'sink-10k-unweighted-unifrac': ('11-packaged/fecal/'
                                                    'notrim/sub_participants/'
                                                    'one_sample/10k/'
                                                    'unweighted_unifrac_ag_'
                                                    '10k_fecal.txt'),
                    'sink-10k-weighted-unifrac': ('11-packaged/fecal/notrim/'
                                                  'sub_participants/'
                                                  'one_sample/10k/'
                                                  'weighted_unifrac_ag_10k_'
                                                  'fecal.txt'),
                    },
                '100nt': {

                    'source-1k-otu': ('11-packaged/fecal/100nt/'
                                      'all_participants/one_sample/1k/ag_1k_'
                                      'fecal.biom'),
                    'source-1k-map': ('11-packaged/fecal/100nt/'
                                      'all_participants/one_sample/1k/ag_1k_'
                                      'fecal.txt'),
                    'source-1k-unweighted-unifrac': ('11-packaged/fecal/100nt'
                                                     '/all_participants/'
                                                     'one_sample/1k/'
                                                     'unweighted_unifrac_ag_'
                                                     '1k_fecal.txt'),
                    'source-1k-weighted-unifrac': ('11-packaged/fecal/100nt/'
                                                   'all_participants/'
                                                   'one_sample/1k/'
                                                   'weighted_unifrac_ag_1k_'
                                                   'fecal.txt'),

                    'source-10k-otu': ('11-packaged/fecal/100nt/'
                                       'all_participants/one_sample/10k/'
                                       'ag_10k_fecal.biom'),
                    'source-10k-map': ('11-packaged/fecal/100nt/'
                                       'all_participants/one_sample/10k/'
                                       'ag_10k_fecal.txt'),
                    'source-10k-unweighted-unifrac': ('11-packaged/fecal/'
                                                      '100nt/all_participants/'
                                                      'one_sample/10k/'
                                                      'unweighted_unifrac_ag_'
                                                      '10k_fecal.txt'),
                    'source-10k-weighted-unifrac': ('11-packaged/fecal/100nt/'
                                                    'all_participants/'
                                                    'one_sample/10k/'
                                                    'weighted_unifrac_ag_10k_'
                                                    'fecal.txt'),

                    'sink-1k-dir': ('11-packaged/fecal/100nt/'
                                    'sub_participants/one_sample/1k'),
                    'sink-1k-otu': ('11-packaged/fecal/100nt/'
                                    'sub_participants/one_sample/1k/'
                                    'ag_1k_fecal.biom'),
                    'sink-1k-map': ('11-packaged/fecal/100nt/'
                                    'sub_participants/one_sample/1k/'
                                    'ag_1k_fecal.txt'),
                    'sink-1k-unweighted-unifrac': ('11-packaged/fecal/100nt/'
                                                   'sub_participants/'
                                                   'one_sample/1k/'
                                                   'unweighted_unifrac_ag_1k_'
                                                   'fecal.txt'),
                    'sink-1k-weighted-unifrac': ('11-packaged/fecal/100nt/'
                                                 'sub_participants/one_sample/'
                                                 '1k/weighted_unifrac_ag_1k_'
                                                 'fecal.txt'),

                    'sink-10k-dir': ('11-packaged/fecal/100nt/'
                                     'sub_participants/one_sample/10k'),
                    'sink-10k-otu': ('11-packaged/fecal/100nt/'
                                     'sub_participants/one_sample/10k/'
                                     'ag_10k_fecal.biom'),
                    'sink-10k-map': ('11-packaged/fecal/100nt/'
                                     'sub_participants/one_sample/10k/'
                                     'ag_10k_fecal.txt'),
                    'sink-10k-unweighted-unifrac': ('11-packaged/fecal/100nt/'
                                                    'sub_participants/'
                                                    'one_sample/10k/'
                                                    'unweighted_unifrac_ag_'
                                                    '10k_fecal.txt'),
                    'sink-10k-weighted-unifrac': ('11-packaged/fecal/100nt/'
                                                  'sub_participants/'
                                                  'one_sample/10k/weighted_'
                                                  'unifrac_ag_10k_fecal.txt'),
                    },
                },
            },
        },
    }

alpha_metrics = {'pd': 'PD_whole_tree',
                 'chao1': 'chao1',
                 'shannon': 'shannon',
                 'observedotus': 'observed_otus'
                 }


def _assert_environment():
    if qiime.__version__ != '1.9.1':
        obs_version = qiime.__version__
        raise ImportError("QIIME 1.9.1 is not in the environment, found "
                          "%s." % obs_version)

    if find_executable('print_qiime_config.py') is None:
        raise EnvironmentError("The QIIME executables are not in $PATH.")

    if find_executable('mod2_pcoa.py') is None:
        raise EnvironmentError("The AG scripts are not in $PATH.")
_assert_environment()


def activate(chp):
    """Activate a chapter

    Parameters
    ----------
    chp : str
        The chapter.

    Returns
    -------
    str
        The path to the activated directory

    Notes
    -----
    Activation creates the chapter processing directory if it does
    not already exist.
    """
    path = os.path.join(ag.WORKING_DIR, chp)
    if not os.path.exists(path):
        os.mkdir(path)
    return path


def get_sortmerna_index():
    """Return the absolute path of a SortMeRNA index if available"""
    return os.environ.get('AG_SMR_INDEX')


def get_rarefaction_depth():
    """Return the rarefaction depth to use"""
    if ag.is_test_env():
        return ("50", "100")
    else:
        return ("1000", "10000")


def get_reference_set():
    """Get the reference set to use for OTU picking

    Returns
    -------
    str
        The file path to the reference sequences.
    str
        The file path to the reference taxonomy.
    """
    if ag.is_test_env():
        repo = get_repository_dir()
        ref_seqs = os.path.join(repo, 'tests/data/otus.fna')
        ref_tax = os.path.join(repo, 'tests/data/otus.txt')
        return ref_seqs, ref_tax
    else:
        return qdr.get_reference_sequences(), qdr.get_reference_taxonomy()


def get_cpu_count():
    """Get the CPU count to use

    Returns
    -------
    int
        The number of CPUs available for use.

    Notes
    -----
    This method provides a layer of indirection in the event that a) a user
    does not want to allow all cores to be used or b) for testing purposes
    as Travis CI appears to be unreliable under multiprocessing and QIIME.
    """
    if os.environ.get('AG_CPU_COUNT') is not None:
        return int(os.environ.get('AG_CPU_COUNT'))
    else:
        # import is scoped as far its use is contained here.
        import multiprocessing
        return multiprocessing.cpu_count()


def get_hmp():
    """Get the HMP 100nt table and mapping"""
    return _get_data('HMP', 'HMPv35_100nt')


def get_pgp():
    """Get the PGP 100nt table and mapping"""
    return _get_data('PGP', 'PGP_100nt')


def get_global_gut():
    """Get the Global Gut table and mapping"""
    return _get_data('GG', 'GG_100nt')


def _get_data(data_dir, tag):
    """Get a non-AG table and mapping file

    Parameters
    ----------
    data_dir : str
        The base data path
    tag : str
        The filetag (e.g., HMPv35_100nt)

    Returns
    -------
    (str, str)
        The filepath to the table, and the filepath to the mapping file.

    Notes
    -----
    If $AG_TESTING == 'True', then the data returned will correspond to the
    test dataset.

    Raises
    ------
    IOError
        If the filepaths are not accessible
    """
    repo = get_repository_dir()
    data = 'tests/data' if ag.is_test_env() else 'data'
    base = os.path.join(repo, data)

    table = os.path.join(base, data_dir, '%s.biom' % tag)
    mapping = os.path.join(base, data_dir, '%s.txt' % tag)

    if not os.path.exists(table):
        raise IOError("Unable to access: %s" % table)
    if not os.path.exists(mapping):
        raise IOError("Unable to access: %s" % table)

    return table, mapping


def get_study_accessions():
    """Get the accessions to use, or redirect to test data

    Returns
    -------
    list of str
        The accessions, which are expected to be basenames for the actual data.
        For instance, the accession "foo" would have sequences as "foo.fna" and
        metadata as "foo.txt".

    Notes
    -----
    If $AG_TESTING == 'True', then the accessions returned will
    correspond to the test dataset.
    """
    if ag.is_test_env():
        _stage_test_accessions()
        return _TEST_ACCESSIONS[:]
    else:
        return _EBI_ACCESSIONS[:]


def get_files(rootdir, suffix):
    """Get the filepaths with a given suffix

    Parameters
    ----------
    rootdir : str
        The root directory to look under
    suffix : str
        The file suffix of the files to keep

    Returns
    -------
    fps : list, str
        List of file paths for all of the
        sample fasta files

    Note
    ----
    This only looks at the directory names under the
    root directory.  This assumes that the sample names
    correspond to the folders within the base folder
    """
    fps = []
    for root, dirs, files in os.walk(rootdir, followlinks=True):
        for _file in files:
            if _file.endswith(".%s" % suffix):
                fps.append(os.path.join(root, _file))
    return fps


def get_bloom_sequences():
    """Get the filepath to the bloom sequences

    Returns
    -------
    str
        The filepath to the bloom sequences

    Raises
    ------
    IOError
        If the path does not exist
    """
    repo = get_repository_dir()
    return get_existing_path(os.path.join(repo, 'data/AG/BLOOM.fasta'))


def _stage_test_accessions():
    """Stage test data

    Notes
    -----
    Staging copies the test dataset into the working directory. This "tricks"
    the fetch_study mechanism as it'll appear that the data have already been
    sourced from EBI.
    """
    repo = get_repository_dir()
    for acc in _TEST_ACCESSIONS:
        src = os.path.join(repo, 'tests/data/%s' % acc)
        dst = os.path.join(ag.WORKING_DIR, '01-raw/%s' % acc)

        if not os.path.exists(os.path.join(ag.WORKING_DIR, '01-raw')):
            os.mkdir('01-raw')

        shutil.copytree(src, dst)


def write_readme(chp_path):
    """Writes the readme file for the processed data"""

    with open(os.path.join(chp_path, 'readme.txt'), 'w') as f_:
        f_.write('\n'.join(readme))


readme = ["American Gut File Sets",
          "----------------------",
          "",
          "The American Gut data has been packaged into directories for easy"
          " access.",
          "",
          "Datasets",
          "========",
          "",
          "American Gut data was primarily divided by sequence trim length, "
          "bodysite and rarefaction depth. Trimmed sequences were used to "
          "facilitate meta analysis with other datasets. Participant results"
          " are returned at a lower rarefaction depth than are used in "
          "analysis.",
          "",
          "We've chosen to organize the data by bodysite, trim length, and "
          "then by rarefaction depth. For each bodysite and trim length, an "
          "unrarefied mapping file and biom table are provided. ",
          "",
          "Data is additionally partitioned into all the samples from all "
          "participants at a particiular bodysite, and a single sample per "
          "indiviudual at each bodysite. The single sample was selected at "
          "random from all the samples from an individual which met the "
          "rarefaction criteria. This sample is used across all datasets. If "
          "no sample met the rarefaction critieria, the level is stepped down"
          " and the proceedure is repeated. This way, the sample selected for"
          " each individual is consistent across all the datasets while "
          "maximizing the number of samples for analysis.",
          "",
          "A list of the barcodes in each single sample dataset are provided "
          "for each bodysite (i.e. `single_1k.txt`).",
          "",
          "For the fecal samples, we additionally filtered the data to include"
          " samples from individuals in a healthy subset of adults. The "
          "criteria for a participant to be included in this group included a"
          " reported age between 20 and 69, a BMI between 18.5 and 30, and no"
          " reported history of Inflammatory Bowel Disease, Diabetes, or "
          "antibiotic use in the past year. A single sample per individual is"
          " provided in each subset.",
          "",
          "Data Dictionary",
          "===============",
          "A data dictionary describing all the base columns in the mapping "
          "file is provided as the `data_dictionary.csv` in the parent "
          "partition directory.",
          "",
          "",
          "Files",
          "=====",
          "Within each dataset directory, the following files are provided:",
          "",
          "Metadata file",
          "+++++++++++++",
          "The mapping file downloaded from EBI. Alpha diversity (PD whole "
          "tree, shannon, choa1, and observed OTUs) for the rarefaction "
          "depth, and every depth lower have been added.",
          "",
          "OTU table",
          "+++++++++",
          "A rarefied biom table.",
          "",
          "Distance Matrices",
          "+++++++++++++++++",
          "The weighted and unweighted UniFrac distance associated with the "
          "set of samples.",
          "",
          "",
          "Analysis Directory Structure",
          "============================",
          "fecal/",
          "   single_ids_1k.txt",
          "   single_ids_10k.txt",
          "   single_ids_unrarefied.txt",
          "   100nt/",
          "       ag_fecal.biom",
          "       ag_fecal.txt",
          "       all_participants/",
          "           all_samples/",
          "               1k/",
          "                   ag_1k_fecal.biom",
          "                   ag_1k_fecal.txt",
          "                   unweighted_unifrac_ag_1k_fecal.txt",
          "                   weighted_unifrac_ag_1k_fecal.txt",
          "               10k/",
          "                   ag_10k_fecal.biom",
          "                   ag_10k_fecal.txt",
          "                   unweighted_unifrac_ag_10k_fecal.txt",
          "                   weighted_unifrac_ag_10k_fecal.txt",
          "           one_sample/",
          "               1k/",
          "                   ag_1k_fecal.biom",
          "                   ag_1k_fecal.txt",
          "                   unweighted_unifrac_ag_1k_fecal.txt",
          "                   weighted_unifrac_ag_1k_fecal.txt",
          "               10k/",
          "                   ag_10k_fecal.biom",
          "                   ag_10k_fecal.txt",
          "                   unweighted_unifrac_ag_10k_fecal.txt",
          "                   weighted_unifrac_ag_10k_fecal.txt",
          "       sub_participants/",
          "           all_samples/",
          "               1k/",
          "                   ag_1k_fecal.biom",
          "                   ag_1k_fecal.txt",
          "                   unweighted_unifrac_ag_1k_fecal.txt",
          "                   weighted_unifrac_ag_1k_fecal.txt",
          "               10k/",
          "                   ag_10k_fecal.biom",
          "                   ag_10k_fecal.txt",
          "                   unweighted_unifrac_ag_10k_fecal.txt",
          "                   weighted_unifrac_ag_10k_fecal.txt",
          "           one_sample/",
          "               1k/",
          "                   ag_1k_fecal.biom",
          "                   ag_1k_fecal.txt",
          "                   unweighted_unifrac_ag_1k_fecal.txt",
          "                   weighted_unifrac_ag_1k_fecal.txt",
          "               10k/",
          "                   ag_10k_fecal.biom",
          "                   ag_10k_fecal.txt",
          "                   unweighted_unifrac_ag_10k_fecal.txt",
          "                   weighted_unifrac_ag_10k_fecal.txt",
          "   notrim/",
          "       100nt/",
          "       ag_fecal.biom",
          "       ag_fecal.txt",
          "       all_participants/",
          "           all_samples/",
          "               1k/",
          "                   ag_1k_fecal.biom",
          "                   ag_1k_fecal.txt",
          "                   unweighted_unifrac_ag_1k_fecal.txt",
          "                   weighted_unifrac_ag_1k_fecal.txt",
          "               10k/",
          "                   ag_10k_fecal.biom",
          "                   ag_10k_fecal.txt",
          "                   unweighted_unifrac_ag_10k_fecal.txt",
          "                   weighted_unifrac_ag_10k_fecal.txt",
          "           one_sample/",
          "               1k/",
          "                   ag_1k_fecal.biom",
          "                   ag_1k_fecal.txt",
          "                   unweighted_unifrac_ag_1k_fecal.txt",
          "                   weighted_unifrac_ag_1k_fecal.txt",
          "               10k/",
          "                   ag_10k_fecal.biom",
          "                   ag_10k_fecal.txt",
          "                   unweighted_unifrac_ag_10k_fecal.txt",
          "                   weighted_unifrac_ag_10k_fecal.txt",
          "       sub_participants/",
          "           all_samples/",
          "               1k/",
          "                   ag_1k_fecal.biom",
          "                   ag_1k_fecal.txt",
          "                   unweighted_unifrac_ag_1k_fecal.txt",
          "                   weighted_unifrac_ag_1k_fecal.txt",
          "               10k/",
          "                   ag_10k_fecal.biom",
          "                   ag_10k_fecal.txt",
          "                   unweighted_unifrac_ag_10k_fecal.txt",
          "                   weighted_unifrac_ag_10k_fecal.txt",
          "           one_sample/",
          "               1k/",
          "                   ag_1k_fecal.biom",
          "                   ag_1k_fecal.txt",
          "                   unweighted_unifrac_ag_1k_fecal.txt",
          "                   weighted_unifrac_ag_1k_fecal.txt",
          "               10k/",
          "                   ag_10k_fecal.biom",
          "                   ag_10k_fecal.txt",
          "                   unweighted_unifrac_ag_10k_fecal.txt",
          "                   weighted_unifrac_ag_10k_fecal.txt",
          "oral/",
          "   single_ids_1k.txt",
          "   single_ids_10k.txt",
          "   single_ids_unrarefied.txt",
          "   100nt/",
          "       ag_oral.biom",
          "       ag_oral.txt",
          "       all_participants/",
          "           all_samples/",
          "               1k/",
          "                   ag_1k_oral.biom",
          "                   ag_1k_oral.txt",
          "                   unweighted_unifrac_ag_1k_oral.txt",
          "                   weighted_unifrac_ag_1k_oral.txt",
          "               10k/",
          "                   ag_10k_oral.biom",
          "                   ag_10k_oral.txt",
          "                   unweighted_unifrac_ag_10k_oral.txt",
          "                   weighted_unifrac_ag_10k_oral.txt",
          "           one_sample/",
          "               1k/",
          "                   ag_1k_oral.biom",
          "                   ag_1k_oral.txt",
          "                   unweighted_unifrac_ag_1k_oral.txt",
          "                   weighted_unifrac_ag_1k_oral.txt",
          "               10k/",
          "                   ag_10k_oral.biom",
          "                   ag_10k_oral.txt",
          "                   unweighted_unifrac_ag_10k_oral.txt",
          "                   weighted_unifrac_ag_10k_oral.txt",
          "   notrim/",
          "       100nt/",
          "       ag_oral.biom",
          "       ag_oral.txt",
          "       all_participants/",
          "           all_samples/",
          "               1k/",
          "                   ag_1k_oral.biom",
          "                   ag_1k_oral.txt",
          "                   unweighted_unifrac_ag_1k_oral.txt",
          "                   weighted_unifrac_ag_1k_oral.txt",
          "               10k/",
          "                   ag_10k_oral.biom",
          "                   ag_10k_oral.txt",
          "                   unweighted_unifrac_ag_10k_oral.txt",
          "                   weighted_unifrac_ag_10k_oral.txt",
          "           one_sample/",
          "               1k/",
          "                   ag_1k_oral.biom",
          "                   ag_1k_oral.txt",
          "                   unweighted_unifrac_ag_1k_oral.txt",
          "                   weighted_unifrac_ag_1k_oral.txt",
          "               10k/",
          "                   ag_10k_oral.biom",
          "                   ag_10k_oral.txt",
          "                   unweighted_unifrac_ag_10k_oral.txt",
          "                   weighted_unifrac_ag_10k_oral.txt",
          "skin/",
          "   single_ids_1k.txt",
          "   single_ids_10k.txt",
          "   single_ids_unrarefied.txt",
          "   100nt/",
          "       ag_skin.biom",
          "       ag_skin.txt",
          "       all_participants/",
          "           all_samples/",
          "               1k/",
          "                   ag_1k_skin.biom",
          "                   ag_1k_skin.txt",
          "                   unweighted_unifrac_ag_1k_skin.txt",
          "                   weighted_unifrac_ag_1k_skin.txt",
          "               10k/",
          "                   ag_10k_skin.biom",
          "                   ag_10k_skin.txt",
          "                   unweighted_unifrac_ag_10k_skin.txt",
          "                   weighted_unifrac_ag_10k_skin.txt",
          "           one_sample/",
          "               1k/",
          "                   ag_1k_skin.biom",
          "                   ag_1k_skin.txt",
          "                   unweighted_unifrac_ag_1k_skin.txt",
          "                   weighted_unifrac_ag_1k_skin.txt",
          "               10k/",
          "                   ag_10k_skin.biom",
          "                   ag_10k_skin.txt",
          "                   unweighted_unifrac_ag_10k_skin.txt",
          "                   weighted_unifrac_ag_10k_skin.txt",
          "   notrim/",
          "       100nt/",
          "       ag_skin.biom",
          "       ag_skin.txt",
          "       all_participants/",
          "           all_samples/",
          "               1k/",
          "                   ag_1k_skin.biom",
          "                   ag_1k_skin.txt",
          "                   unweighted_unifrac_ag_1k_skin.txt",
          "                   weighted_unifrac_ag_1k_skin.txt",
          "               10k/",
          "                   ag_10k_skin.biom",
          "                   ag_10k_skin.txt",
          "                   unweighted_unifrac_ag_10k_skin.txt",
          "                   weighted_unifrac_ag_10k_skin.txt",
          "           one_sample/",
          "               1k/",
          "                   ag_1k_skin.biom",
          "                   ag_1k_skin.txt",
          "                   unweighted_unifrac_ag_1k_skin.txt",
          "                   weighted_unifrac_ag_1k_skin.txt",
          "               10k/",
          "                   ag_10k_skin.biom",
          "                   ag_10k_skin.txt",
          "                   unweighted_unifrac_ag_10k_skin.txt",
          "                   weighted_unifrac_ag_10k_skin.txt",
          ]
