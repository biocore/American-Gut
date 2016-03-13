American-Gut
============

American Gut open-access data and IPython notebooks

# INSTALL

## Basics

Amering-Gut repository is intended to be used as a project/repo
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

### PIP based
```bash
$ pip install numpy==1.9.2
$ pip install -r ./pip_requirements.txt
```

If you would like to install dependencies within a virualenv environment be
sure to change to the appropriate environment prior to the installation
of dependencies.
	
*Note*: Be aware that with pip some libraries will have to be compiled
from source so having appropriate system libraries should be installed
prior to running the pip command. For more details take a look
at Supported Systems chapter.

## Supported Operating Systems / Distributions

### Debian 8

Tested with Debian 8.3.0 (amd64).

To compile dependencies from source appropriate libraries can be installed
(as root/sudo) with
```bash
(root/sudo)$ aptitude install pkg_config libxslt1-dev libxml2 libfreetype6 \
    buildessential python-pip python-dev liblapack-dev liblapack3 \
    libblas-dev libblas3 gfortran libhdf5-serial-dev
```
