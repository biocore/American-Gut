American-Gut
============

American Gut open-access data and IPython notebooks

# INSTALL

## Basics

American-Gut repository is intended to be used as a project/repo
meaning there is no need to install it (ignore `setup.py` at the moment).

After cloning the repository and before using the scripts user should
install necessary dependencies. Two approaches are supported at the moment.

### Conda based

If you're choice of package manager is conda dependencies can be
installed with

```bash
$ conda install --file ./conda_requirements.txt
$ pip install -r ./pip_requirements.txt
```

If you would like to install dependencies within a conda environment be
sure to change to the appropriate environment prior to the installation
of dependencies.

***Note***: Be aware that with pip some libraries will have to be compiled
from source so appropriate system libraries should be installed
prior to running the pip command. For more details take a look
at Supported Systems section.


### PIP based
```bash
$ pip install numpy==1.9.2
$ pip install -r ./pip_requirements.txt
```

If you would like to install dependencies within a virtualenv environment be
sure to change to the appropriate environment prior to the installation
of dependencies.
	
***Note***: Be aware that with pip some libraries will have to be compiled
from source so appropriate system libraries should be installed
prior to running the pip command. For more details take a look
at Supported Systems section.

## Supported Operating Systems / Distributions

### Debian 8

Tested with Debian 8.3.0 (amd64).

To compile dependencies from source appropriate libraries can be installed
(as root/sudo) with
```bash
(root/sudo)$ aptitude install pkg-config libxslt1-dev libxml2 libfreetype6 \
    build-essential python-pip python-dev liblapack-dev liblapack3 \
    libfreetype6-dev libblas-dev libblas3 gfortran libhdf5-serial-dev libsm6
```

# RUN

## Basics

Although American-Gut repo provides separate scripts (`scripts` folder)
and a package (`americangut` folder) it is primarily intended to be used
through notebooks (`ipynb` folder).

There are a few environment variable that can be used to customize the run:

- AG\_TESTING: if set to `True` scripts will not download AmericanGut
  EBI data (ERP012803) bit instead work with test data (subset of the original
  EBI data). This is useful for testing.
- AG\_CPU\_COUNT: Number of process to use when parallelizing code (defaults to
  the number of cores)

To generate reports (pdfs) a TeX distribution should be installed on the system.

## Adjusting environment on POSIX systems

Since American-Gut repo contains scripts and packages we need to adjust
PYTHONPATH and PATH to reflect this. Therefore, prior to working with notebooks
execute the following from within the American-Gut repo:

```bash
REPO=`pwd`
$ export PYTHONPATH=$REPO/:$PYTHONPATH
$ export PATH=$REPO/scripts:$PATH
```

If needed adjust `AG_*` environment variables from Basics section.

## Run notebooks

Notebooks are written in two formats and therefore require
different profiles.

### Markdown based notebooks

Markdown based notebooks can be found in `./ipynb/primary-processing/` folder
and have extension `md`.
To use these notebooks we first need to create a profile for `ag_ipymd`
with

```bash
$ ipython profile create ag_ipymd
```

and adjust newly created `/path/to/.ipython/profile_ag_ipymd/ipython_notebook_config.py'
by addding
```
#------------------------
# ipymd
#------------------------
c.NotebookApp.contents_manager_class = 'ipymd.IPymdContentsManager'
```
to the end of the file.

Now, we can start ipython with
```bash
$ ipython notebook --profile=ag_ipymd
```
and visit the newly started notebook server by going to http://localhost:8888

### Markdown based notebooks

Notebooks in native notebook format (ipynb) in `./ipynb/` folder
and have the extension `ipynb`.
To use these notebooks we first need to create a profile for `ag_default`
with

```bash
$ ipython profile create ag_default
```

Now, we can start ipython with
```bash
$ ipython --profile=ag_default notebook
```
and visit the newly started notebook server by going to http://localhost:8888
