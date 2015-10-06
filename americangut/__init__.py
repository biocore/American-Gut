import os

working_dir = 'agp_processing'
if not os.path.exists(working_dir):
    os.mkdir(working_dir)


__all__ = ['working_dir']
