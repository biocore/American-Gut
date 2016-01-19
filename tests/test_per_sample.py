import os
import shutil
from unittest import TestCase, main

import numpy as np
import pandas as pd
import numpy.testing as npt

import americangut as ag
import americangut.notebook_environment as agenv
import americangut.per_sample as agps


class PerSampleTests(TestCase):
    def test_create_opts(self):
        obs = agps.create_opts('fecal', 'somepath', gradient_color_by='foo',
                               barchart_categories=('sex', 'age'))
        self.assertEqual(obs['gradient_color_by'], 'foo')
        self.assertIn('AGE_CAT', obs['barchart_categories'])
        self.assertIn('SEX', obs['barchart_categories'])
        self.assertEqual(obs['sample_type'], 'fecal')

    def test_sample_type_processor(self):
        funcs = [lambda a, b: {c: c for c in b},
                 lambda a, b: {c: None for c in b},
                 lambda a, b: {c: a[c] for c in b if c != 'c'}]
        ids = ['a', 'b', 'c', 'd']
        opts = {'a': 'opt-a', 'b': 'opt-b', 'c': 'opt-c', 'd': 'opt-d'}
        exp = {'a': ['a', 'opt-a'],
               'b': ['b', 'opt-b'],
               'c': ['c'],
               'd': ['d', 'opt-d']}
        obs = agps.sample_type_processor(funcs, opts, ids)
        self.assertEqual(obs, exp)

    def test_result_path(self):
        self.assertEqual(agps._result_path(agenv.paths, 'a'),
                         agenv.paths['per-sample']['results'] + '/a')

    def test_base_barcode(self):
        self.assertEqual(agps._base_barcode('asd.foo'), 'foo')

    def test_partition_samples_by_bodysite(self):
        df = pd.DataFrame()
        df['SIMPLE_BODY_SITE'] = ['FECAL', 'FECAL', 'SKIN', 'ORAL', 'ORAL',
                                  'NA', 'FECAL']
        site_to_functions = [('FECAL', 'a fecal function'),
                             ('SKIN', 'a skin function'),
                             ('ORAL', 'a oral function')]
        parts = list(agps.partition_samples_by_bodysite(df, site_to_functions))
        fecal_parts, skin_parts, oral_parts = parts
        self.assertEqual(fecal_parts[0], 'a fecal function')
        self.assertEqual(oral_parts[0], 'a oral function')
        self.assertEqual(skin_parts[0], 'a skin function')
        npt.assert_array_equal(fecal_parts[1], [0, 1, 6])
        npt.assert_array_equal(oral_parts[1], [3, 4])
        npt.assert_array_equal(skin_parts[1], [2])

    def test_merge_error_reports(self):
        report_a = {'a': None, 'b': None, 'c': 'foo'}
        report_b = {'a': None, 'b': 'bar', 'c': 'bar'}
        report_c = {'a': None, 'b': None, 'c': 'baz'}
        exp = {'a': [], 'b': ['bar'], 'c': ['foo', 'bar', 'baz']}
        obs = agps.merge_error_reports(*[report_a, report_b, report_c])
        self.assertEqual(obs, exp)

    def test_iter_ids_over_system_call(self):
        cmd_fmt = "ls %(result_path)s %(id)s"
        ids = ['/', '/asdasd', '/usr']
        opts = {'per-sample': {'results': ''}}

        obs = agps._iter_ids_over_system_call(cmd_fmt, ids, opts)

        # ls behaves differently on BSD and Linux
        self.assertEqual(obs['/'], None)
        self.assertTrue(obs['/asdasd'].startswith('FAILED (ls: '))
        self.assertTrue(obs['/asdasd'].endswith('ls /asdasd /asdasd'))
        self.assertEqual(obs['/usr'], None)

    # When the unit test suite is run, we cannot assume that the expected
    # inputs to these methods are available. The intent of these next
    # tests are to verify that the expected commands for system calls are
    # being formed as expected indirectly via forcing failures.
    def test_taxon_significance(self):
        exp = {'test': ('FAILED (generate_otu_signifigance_tables_AGP.py: '
                        'error: The supplied taxonomy file does not exist '
                        'in the path.): '
                        'generate_otu_signifigance_tables_AGP.py '
                        '-i bar.biom -o foo/test -s test -m baz')}
        ids = ['test']
        opts = {'per-sample': {'results': 'foo'},
                'taxa': {'notrim': {'L6': {'ag-bar-biom': 'bar.biom'}}},
                'meta': {'ag-cleaned-md': 'baz'},
                'sample_type': 'bar'}

        obs = agps.taxon_significance(opts, ids)
        self.assertEqual(obs, exp)

    def test_body_site_pcoa(self):
        exp = {'test': ('FAILED (Error: Invalid value for "--coords": Path '
                        '"foo" does not exist.): mod2_pcoa.py body_site '
                        '--coords foo --mapping_file bar --output baz/test '
                        '--filename figure1.pdf --sample test')}
        ids = ['test']
        opts = {'per-sample': {'results': 'baz'},
                'beta':
                    {'100nt':
                        {'1k':
                            {'ag-pgp-hmp-gg-unifrac-pc': 'foo'}}},
                'meta': {'ag-pgp-hmp-gg-cleaned-md': 'bar'}}

        obs = agps.body_site_pcoa(opts, ids)
        self.assertEqual(obs, exp)

    def test_countly_pcoa(self):
        exp = {'test': ('FAILED (Error: Invalid value for "--distmat": Path '
                        '"foo" does not exist.): mod2_pcoa.py country '
                        '--distmat foo --coords bar --mapping_file baz '
                        '--output foobar/test --filename figure2.pdf '
                        '--sample test')}
        ids = ['test']
        opts = {'per-sample': {'results': 'foobar'},
                'beta':
                    {'100nt':
                        {'1k':
                            {'ag-gg-subsampled-unifrac-pc': 'bar',
                             'ag-gg-unifrac': 'foo'}}},
                'meta': {'ag-gg-cleaned-md': 'baz'}}

        obs = agps.country_pcoa(opts, ids)
        self.assertEqual(obs, exp)

    def test_gradient_pcoa(self):
        exp = {'test': ('FAILED (Error: Invalid value for "--coords": Path '
                        '"foo" does not exist.): mod2_pcoa.py gradient '
                        '--coords foo --mapping_file bar --output baz/test '
                        '--filename figure3.pdf --color foobar '
                        '--sample test')}
        ids = ['test']
        opts = {'per-sample': {'results': 'baz'},
                'beta': {'100nt': {'1k': {'ag-what-unifrac-pc': 'foo'}}},
                'taxa': {'notrim': {'L2': {'ag-md': 'bar'}}},
                'sample_type': 'what',
                'gradient_color_by': 'foobar'}

        obs = agps.gradient_pcoa(opts, ids)
        self.assertEqual(obs, exp)

    def test_pie_plot(self):
        exp = {'test': ('FAILED (make_pie_plot_AGP.py: error: The supplied '
                        'taxonomy file does not exist in the path.): '
                        'make_pie_plot_AGP.py -i foo -o bar/test -s test')}
        ids = ['test']
        opts = {'per-sample': {'results': 'bar'},
                'taxa': {'notrim': {'L3': {'ag-tsv': 'foo'}}}}

        obs = agps.pie_plot(opts, ids)
        self.assertEqual(obs, exp)

    def test_bar_chart(self):
        exp = {'test': ('FAILED (make_phyla_plots_AGP.py: error: The supplied '
                        'biom table does not exist in the path.): '
                        'make_phyla_plots_AGP.py -i foo -m baz -o bar/test '
                        '-c stuff -t what -s test')}
        ids = ['test']
        opts = {'per-sample': {'results': 'bar'},
                'collapsed':
                    {'notrim':
                        {'1k': {'ag-what-biom': 'foo'}}},
                'meta': {'ag-cleaned-md': 'baz'},
                'barchart_categories': 'stuff',
                'sample_type': 'what'}

        obs = agps.bar_chart(opts, ids)
        self.assertEqual(obs, exp)

    def test_taxa_summaries(self):
        ids = ['USygt45.M.418662', 'missing']
        exp = {'USygt45.M.418662': None, 'missing': "ID not found"}
        opts = {'per-sample': {'results': 'bar'},
                'taxa': {'notrim': {'L6': {'ag-fecal-biom':
                                           agenv.get_global_gut()[0]}}},
                'sample_type': 'fecal'}

        os.mkdir('bar')
        os.mkdir('bar/USygt45.M.418662')

        obs = agps.taxa_summaries(opts, ids)
        self.assertEqual(obs, exp)

        shutil.rmtree('bar')

    def test_per_sample_directory(self):
        ids = ['test']
        opts = {'per-sample': {'results': 'bar'}}
        os.mkdir('bar')
        agps.per_sample_directory(opts, ids)

        self.assertTrue(os.path.exists('bar/test'))
        os.rmdir('bar/test')
        os.rmdir('bar')

    def test_stage_per_sample_specific_statics(self):
        opts = {'chp-path': '.',
                'per-sample': {'statics-fecal': 'stuff',
                               'results': ''},
                'sample_type': 'fecal'}
        ids = ['foo']
        exp = {'foo': "Cannot stage template."}
        obs = agps.stage_per_sample_specific_statics(opts, ids)
        self.assertEqual(obs, exp)

    def test_plot_alpha(self):
        map_ = pd.DataFrame(
            data=np.array([
                ['skin', '1990', 'female', 'Verity', 'US', 12.5],
                ['fecal', '1990', 'female', 'Verity', 'US', 8.6],
                ['fecal', '1987', 'male', 'Alex', 'US', 7.9],
                ['fecal', '1993', 'female', 'Annie', 'US', 7.5],
                ['skin', '1989', 'male', 'Dominic', 'UK', 14.0],
                ['fecal', '1986', 'female', 'Sarah', 'US', 15.0],
                ['oral', '1988', 'female', 'Shelby', 'AUS', 4.2],
                ]),
            index=['VeP0', 'VeP1', 'AxP0', 'AnP0', 'DoD0', 'SaZ0', 'ShT0'],
            columns=['SIMPLE_BODY_SITE', 'BIRTH_YEAR', 'SEX',
                     'HOST_SUBJECT_ID', 'NATIONALITY', 'alpha'],
            )
        sample = 'VeP0'
        kgroup = 'skin'
        kgalpha = pd.Series([12.5, 14.0],
                            index=['VeP0', 'DoD0'],
                            name='alpha',
                            dtype=object)
        ksalpha = 12.5
        kxlabel = 'alphadiversity'

        (tgroup, tgalpha, tsalpha, txlabel) = \
            agps._plot_alpha(sample=sample,
                             alpha_map=map_,
                             alpha_field='alpha',
                             debug=True)

        self.assertEqual(tgroup, kgroup)
        npt.assert_array_equal(tgalpha.values, kgalpha.values)
        self.assertEqual(tsalpha, ksalpha)
        self.assertEqual(txlabel, kxlabel)

    def test_alpha_plot_metric_error_field_error(self):
        opts = {'collapsed': {
                    '100nt': {
                        'alpha-map':
                        os.path.join(ag.WORKING_DIR.split('American-Gut')[0],
                                     'American-Gut/tests/files/'
                                     'test_mapping.txt')}},
                'sample_type': 'fecal'}
        with self.assertRaises(ValueError):
            agps.alpha_plot(opts, ['sample_a', 'sample_b'])

    def test_alpha_plot_group_field_error(self):
        opts = {'collapsed': {
                    '100nt': {
                        'alpha-map':
                        os.path.join(ag.WORKING_DIR.split('American-Gut')[0],
                                     'American-Gut/tests/files/'
                                     'test_mapping_alpha.txt')}},
                'sample_type': 'fecal'}
        with self.assertRaises(ValueError):
            agps.alpha_plot(opts, ['sample_a', 'sample_b'])

if __name__ == '__main__':
    main()
