import os

WORKING_DIR = 'agp_processing'
if not os.path.exists(WORKING_DIR):
    os.mkdir(WORKING_DIR)


__all__ = ['WORKING_DIR']
