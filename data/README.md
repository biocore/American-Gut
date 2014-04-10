Data
====

OTU Tables and QIIME-compliant mapping files used to generate figures and statistics for the American Gut project. All data are de-identified.

Data sources
------------

The following studies are being used to provide context for the American Gut data:

* [Human Microbiome Project](http://www.ncbi.nlm.nih.gov/pubmed/22699609), v35 reads. The sequence data were generated using 454 instruments.
* [Global Gut](http://www.ncbi.nlm.nih.gov/pubmed/22699611). The sequence data were generated using HiSeq instruments.
* [Personal Genome Project](http://personalgenomes.org) microbiome samples, unpublished. The sequence data were generated using MiSeq instruments.

Data prefixes
-------------

Each study used is described by an acronym:

* AG, American Gut
* HMP, Human Microbiome Project
* GG, Global Gut
* PGP, Personal Genome Project

File tags
---------

The provided [BIOM](http://biom-format.org) tables have a few different tags in the filenames to describe the included data.

* 100nt - The sequences were trimmed to 100 nucleotides prior to OTU picking
* even1k - The full table was rarified to 1000 sequences per sample
* even10k - The full table was rarified to 10000 sequences per sample

The trimming is necessary when combining data from studies in which different sequences technologies were used (e.g., HiSeq vs. MiSeq).

Debug data
----------

The debug data files are sourced from the main data files, but are 10% random subsets (by sample) of what is in them main files. The purpose of the debug files is to reduce processing load on the results framework for testing purposes. 
