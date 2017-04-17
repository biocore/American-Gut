import os

BASE_DIR = os.path.join(os.path.abspath(__file__).split('American-Gut')[0],
                        'American-Gut')
WORKING_DIR = os.path.join(BASE_DIR, 'ipynb/primary-processing/agp_processing')
preprocessing_dir = os.environ.get('AG_WORKING_DIR', 'agp_processing')
if not os.path.exists(WORKING_DIR):
    os.mkdir(WORKING_DIR)

_TEST_ENV = os.environ.get('AG_TESTING') == 'True'


def is_test_env():
    return _TEST_ENV


def staged_raw_data():
    """Get the path of staged raw data or None"""
    return os.environ.get('AG_RAW_DATA')

__version__ = "0.0.1"

__all__ = ['WORKING_DIR']
