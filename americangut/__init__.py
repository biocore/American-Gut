import os

WORKING_DIR = os.path.join(os.path.abspath('.'), 'agp_processing')
if not os.path.exists(WORKING_DIR):
    os.mkdir(WORKING_DIR)

_TEST_ENV = os.environ.get('AG_TESTING') == 'True'


def is_test_env():
    return _TEST_ENV

__version__ = "0.0.1"

__all__ = ['WORKING_DIR']
