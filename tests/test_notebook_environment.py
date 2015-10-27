from unittest import TestCase, main

import qiime

import americangut as ag
import americangut.notebook_environment as agenv


class NotebookEnvironmentTests(TestCase):
    def test_assert_environment(self):
        self.assertEqual(agenv._assert_environment(), None)
        old = qiime.__version__
        qiime.__version__ = 'foo'
        with self.assertRaises(ImportError):
            agenv._assert_environment()
        qiime.__version__ = old

    def test_get_study_accessions(self):
        ag._TEST_ENV = ''
        self.assertEqual(agenv.get_study_accessions(), agenv._EBI_ACCESSIONS)
        ag._TEST_ENV = 'True'
        self.assertEqual(agenv.get_study_accessions(), agenv._TEST_ACCESSIONS)

    def test_get_bloom_sequences(self):
        self.assertEqual(agenv.get_bloom_sequences().rsplit('/')[-1],
                         'BLOOM.fasta')


if __name__ == '__main__':
    main()
