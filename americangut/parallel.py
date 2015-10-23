import multiprocessing as mp

import americangut.results_utils as agru
import americangut.notebook_environment as agenv


def dispatcher(success_fp, fail_fp, partitions):
    """Dispatch execution over a pool of processors

    Parameters
    ----------
    success_fp : file-like object
        A file-like object to write a list of successful sample IDs too
    fail_fp : file-like object
        A file-like object to write a list of unsuccessful sample IDs too,
        and any associated error information
    partitions : iterable
        Yields a function and an iterable of IDs. It is expected that this
        function will return a dict of {str: list}.
    """
    pool = mp.Pool(processes=agenv.get_cpu_count())

    success_fp.write('%s\n' % '#SampleID')
    fail_fp.write('%s\t%s\n' % ('#SampleID', 'Error(s)'))

    for func, ids in partitions:
        for success_details in pool.map(func, list(agru.chunk_list(ids))):
            for id_, detail in success_details.items():
                if detail:
                    fail_fp.write("%s\t%s\n" % (id_, '\t'.join(detail)))
                else:
                    success_fp.write("%s\n" % id_)
