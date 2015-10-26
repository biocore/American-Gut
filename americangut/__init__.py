import os

WORKING_DIR = os.path.join(os.path.abspath('.'), 'agp_processing')
if not os.path.exists(WORKING_DIR):
    os.mkdir(WORKING_DIR)

__version__ = "0.0.1"

__all__ = ['WORKING_DIR']
