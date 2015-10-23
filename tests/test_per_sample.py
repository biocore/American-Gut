from unittest import TestCase, main

import pandas as pd
import numpy.testing as npt

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
                         agenv.paths['per-sample-results'] + '/a')

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
        opts = {'per-sample-results': ''}

        # the command for /asdasd is duplicated as that is the effect of the
        # formatting with "result_path" as it encodes the ID as well.
        exp = {'/': None, '/asdasd': ('FAILED (ls: /asdasd: No such file '
                                      'or directory): ls /asdasd /asdasd'),
               '/usr': None}
        obs = agps._iter_ids_over_system_call(cmd_fmt, ids, opts)
        self.assertEqual(obs, exp)

    # When the unit test suite is run, we cannot assume that the expected
    # inputs to these methods are available. The intent of these next
    # tests are to verify that the expected commands for system calls are
    # being formed as expected indirectly via forcing failures.
    def test_taxon_significance(self):
        exp = {'test': ('FAILED (generate_otu_signifigance_tables_AGP.py: '
                        'error: The supplied taxonomy file does not exist '
                        'in the path.): '
                        'generate_otu_signifigance_tables_AGP.py '
                        '-i bar.biom -o foo/test -s test')}
        ids = ['test']
        opts = {'per-sample-results': 'foo',
                'ag-L6-taxa-bar-biom': 'bar.biom',
                'sample_type': 'bar'}

        obs = agps.taxon_significance(opts, ids)
        self.assertEqual(obs, exp)

    def test_body_site_pcoa(self):
        exp = {'test': ('FAILED (Error: Invalid value for "--coords": Path '
                        '"foo" does not exist.): mod2_pcoa.py body_site '
                        '--coords foo --mapping_file bar --output baz/test '
                        '--prefix Figure_1 --samples test')}
        ids = ['test']
        opts = {'per-sample-results': 'baz',
                'ag-pgp-hmp-gg-100nt-1k-unifrac-pc': 'foo',
                'ag-cleaned-md': 'bar'}

        obs = agps.body_site_pcoa(opts, ids)
        self.assertEqual(obs, exp)

    def test_countly_pcoa(self):
        exp = {'test': ('FAILED (Error: Invalid value for "--distmat": Path '
                        '"foo" does not exist.): mod2_pcoa.py country '
                        '--distmat foo --coords bar --mapping_file baz '
                        '--output foobar/test --prefix Figure_2 '
                        '--samples test')}
        ids = ['test']
        opts = {'per-sample-results': 'foobar',
                'ag-gg-100nt-1k-subsampled-unifrac-pc': 'bar',
                'ag-gg-100nt-1k-bdiv-unifrac': 'foo',
                'ag-gg-cleaned-md': 'baz'}

        obs = agps.country_pcoa(opts, ids)
        self.assertEqual(obs, exp)

    def test_gradient_pcoa(self):
        exp = {'test': ('FAILED (Error: Invalid value for "--coords": Path '
                        '"foo" does not exist.): mod2_pcoa.py gradient '
                        '--coords foo --mapping_file bar --output baz/test '
                        '--prefix Figure_3 --color foobar '
                        '--samples test')}
        ids = ['test']
        opts = {'per-sample-results': 'baz',
                'ag-100nt-what-1k-unifrac-pc': 'foo',
                'ag-L2-taxa-md': 'bar',
                'sample_type': 'what',
                'gradient-color-by': 'foobar'}

        obs = agps.gradient_pcoa(opts, ids)
        self.assertEqual(obs, exp)

    def test_pie_plot(self):
        exp = {'test': ('FAILED (make_pie_plot_AGP.py: error: The supplied '
                        'taxonomy file does not exist in the path.): '
                        'make_pie_plot_AGP.py -i foo -o bar/test -s test')}
        ids = ['test']
        opts = {'per-sample-results': 'bar',
                'ag-L3-taxa-tsv': 'foo'}

        obs = agps.pie_plot(opts, ids)
        self.assertEqual(obs, exp)

    def test_bar_chart(self):
        self.fail()

    def test_taxa_summaries(self):
        self.fail()
    def test_per_sample_directory(self):
        self.fail()
    def test_stage_per_sample_specific_statics(self):
        self.fail()

if __name__ == '__main__':
    main()
