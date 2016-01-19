import multiprocessing as mp
import logging
import traceback
import sys
from functools import partial

import matplotlib
matplotlib.use('Agg')

import americangut as ag
import americangut.results_utils as agru
import americangut.notebook_environment as agenv


def run_functor(functor, *args, **kwargs):
    """
    Given a functor, run it and return its result. We can use this with
    multiprocessing.map and map it over a list of job functors to do them.

    Handles getting more than multiprocessing's pitiful exception output

    This function was derived from:
    http://stackoverflow.com/a/16618842/19741
    """
    try:
        # This is where you do your actual work
        return functor(*args, **kwargs)
    except:
        # Put all exception text into an exception and raise that
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))


def dispatcher(success_fp, fail_fp, partitions):
    """Dispatch execution over a pool of processors

    Parameters
    ----------
    success_fp : file-like object
        A file-like object to write a list of successful sample IDs too
    fail_fp : file-like object
        A file-like object to write a list of unsuccessful sample IDs too,
        and any associated error information
    partitions : Iterable of (function, Iterable of str)
        Yields a function and an iterable of IDs. It is expected that the
        functions yielded will have the following signature:

        {str: list} <- function(list of str)
    """
    if ag.is_test_env():
        logger = mp.log_to_stderr()
        logger.setLevel(logging.INFO)

    pool = mp.Pool(processes=agenv.get_cpu_count())

    success_fp.write('%s\n' % '#SampleID')
    fail_fp.write('%s\t%s\n' % ('#SampleID', 'Error(s)'))

    for func, ids in partitions:
        functor = partial(run_functor, func)
        for success_details in pool.map(functor, list(agru.chunk_list(ids))):
            for id_, detail in success_details.items():
                if detail:
                    fail_fp.write("%s\t%s\n" % (id_, '\t'.join(detail)))
                else:
                    success_fp.write("%s\n" % id_)
