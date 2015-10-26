# Contributing to American Gut

The [American Gut Project](http://americangut.org) is the largest a crowd-funded science project looking to map the topology of the human superorganisms. One of the goals of the American Gut Project is transparency about data processing and technique development. You can find source code used in American Gut analyses under public revision in the American Gut Repository on [Github](https://github.com/biocore/American-Gut).

This document covers what you should do to get started contributing to American Gut. You should read this whole document before you consider submitting code to American Gut. This will save time for both you and the American Gut developers.

## Types of Submissions
We are looking for submission of new analyses (although it may be a good idea to discuss your analysis with the development team before submitting your pull request), bug fixes, documentation updates, additions, and fixes.

When considering submitting a new analysis to American Gut, you should begin by posting an issue on [the American Gut Issue Tracker](https://github.com/biocore/American-Gut/issues). The information which needs to be included in this post will differ based on the type of contribution. Your contribution will also need to be tested (discussed below).

* For new features, describe why the functionality being proposed is relevant. The functionality must be demonstrated as relevant to other users, or other analyses. If appropriate, you may be encouraged to push the functionality to other biocore packages, such as [QIIME](https://github.com/biocore/qiime) or [scikit-bio](https://github.com/biocore/scikit-bio).

* For new analyses, you’ll want to describe the new approach. Explain why you have selected this approach, citing as necessary. Your analysis should be presented as an Jupyter Notebook (see below).

* For bug fixes, you should provide a detailed description of the bug so other developers can reproduce it. Bugs may relate to errors in code, documentation, test, or the data hosted in the repository.
<br/>You should the following information in your bug report:
    * The exact command or function call that can be run to reproduce the bug.
    * A link to all necessary input files for reproducing the bug, or a list of samples which produce the bug. This is *extremely* useful to other developers, and it is likely that if you don’t provide this information, you will not get the response you’re asking for. 
<br>Often, this process will help you understand the bug, as well.

* For documentation issues, please first post an issue describing what you propose to add, where you’d like to add it, and a description of why it is an important addition/ For documentation improvements and fixes, you should post an issue of what is currently wrong or missing, and how you propose to address it. 

When you post your issue, the American Gut developers will respond to let you know if we agree with the addition or change. It’s important to go through this step to avoid wasting time working on a feature that will not be included.

## Code Review
When you submit code in American Gut, it will be reviewed by one or more American Gut developer. These reviews are intended to confirm a few points:
* Your code is sufficiently well tested (see Testing Guidelines below)
* Your analyses adheres to submission guidelines (see Jupyter Notebooks below)
* Your code and analysis are sufficiently well documented
* Your update provides relevant changes or additions.

This process is designed to ensure the quality of American Gut submissions, and can provide useful experience for new developers.

For big changes, if you’d like feedback on your code as you work, you should request help in the issue that you created, and one of the American Gut developers will work with you to perform regular code review. 

## Submitting Code to American Gut
American Gut is hosted on [GitHub](http://www.github.com), and we use GitHub's [Pull Request](https://help.github.com/articles/using-pull-requests) mechanism for accepting submissions. You should go through the following steps to submit code or analyses to American Gut.

* Begin by [creating an issue](https://github.com/biocore/american-gut/issues) describing your proposed analysis or change. This should include a description of the proposed change, and a note in the issue documentation that you'd like to work on it. Once you hear back from a maintainer that it's okay to make the change, we will assign the issue to you.
* [Fork](https://help.github.com/articles/fork-a-repo) the American Gut repository on the Github Website to your Github account.
* Clone your forked repository to the system where you'll be developing using `git clone`
* Ensure that you have the latest version of all the files (especially important if you cloned a long time ago). You should do this by adding American Gut as a remote repository and then pulling from that repository. You'll only need to run the `git remote` step one time:
```
git checkout master
git remote add upstream master http://github.com/biocore/american-gut.git
git pull upstream master
```
* Create a new topic branch that you will use to make your changes with `git checkout -b`:<br>
```git checkout -b my-topic-branch```

* Make changes to your branch. You can add them using `git add` and `git commit`. Don't forget to update the associated scripts and tests. You should make incremental commits, rather than waiting to make one massive commit. Write a descriptive message for each commit.
* When you think you're ready to submit your contribution, again insure that you have the latest version of the repository, incase something changed while you were working on your edits.
* Test your code, to make sure nothing is unexpectedly broken.
* Once the tests past, you should push your changes to your forked GitHub repository. This can be done through the command line, using the command,<br>
```
git push origin my-topic-branch
```

* Issue a [pull request](https://help.github.com/articles/using-pull-requests) on the GitHub website to requset that we merge your changes. One of the American Gut developers will review your code. If we request changes (which is highly probable), *do not issue a new pull request*. Your pull request will be updated automatically.

## Coding Guidelines

We adhere to the [PEP 8](http://www.python.org/dev/peps/pep-0008/) python coding guidelines for code and documentation standards. Before submitting any code to scikit-bio, you should read these carefully and apply the guidelines in your code.

## Testing Guidelines
Code submitted to American Gut should be unit tested as much as possible. Tests should be added to the *test* directory. 

# Jupyter Notebooks for Analysis
We use [Jupyter Notebooks](http://jupyter.org) to document and describe analyses performed on American Gut data. There are two primary goals of these notebooks. First, they create a way to reproducibly generate results. They also provide an opportunity for scientific communication and public outreach. Jupyter has the advantage over previous IPython notebooks in that they can accept a variety of backend kernels, including R and Julia.

## Audience
Jupyter Notebooks are intended to be read by a broad audience. For instance, they may be useful to individual participants, interested in a particular topic, or members of the scientific community looking for more information about the analysis being performed in a paper.

Notebooks represent a *sample calculation* for analysis. If a similar set of analysis techniques is applied to the same set of data, multiple notebooks do not need to be produced for each analysis. The results of the notebooks can then be hosted at another location.

We may publish analysis links on the [American Gut blog](http://americangut.org/?page_id=202). We recommend providing both documentation for the code used, and logic behind the analysis. A lengthy introduction section, or text describing variables may also be useful.


## Index Notebook

Please update the [index notebook](https://github.com/biocore/American-Gut/blob/master/ipynb/index.ipynb) when a new analysis is added.

## Sections
Jupyter Notebooks should more or less follow the outline of a paper.

### Credits and License
American Gut analyses are distributed under a BSD license. Please also provide a data of analysis. Typical text is 

<blockquote><strong>License</strong>: BSD<br>
<strong>Copyright</strong>: Copyright American Gut Project, 2015<br>
</blockquote>


You may also add the following optional information:

* Author name
* Author contact (i.e. email, twiter)
* Last update


### Table of Contents
Describe the major sections in the notebook. We recomend using the document tools from the [ICalico Project](http://calicoproject.org/ICalico) to create a table of contents from the headers in the notebook.

### Introduction
A brief introduction to the topic addressed in the notebook. This may include theory behind the analysis technique, a brief review of relevant literature and/or a description of the data being analyzed. Scientific citations may be included in this section. The introduction may also be subdivided, as appropriate.

### Notebook Requirements
This provides a listing of the current versions of software used to generate the notebook. It is designed to provide a consistent environment. Typically, software compatible with the current version of QIIME is recommended.

### Function Imports
Libraries may be imported into Jupyter. It is recommended most code be passed to the notebook, rather than writing functions in the notebook. Separating code from presentation makes it easier to test. It also increases the probability that the code will be reproduced in other environments.

### Parameter and File Path Definition
We suggest setting analysis parameters and file paths before the analysis is performed. This clarifies the information, and can help avoid bugs if parameters of file paths need to be changed.

### Analysis Steps
The analysis pipeline should then be detailed. We suggest using markdown cells to provide relevant theory about the analysis step, as well as comments in the code cells.

### Discussion or Conclusions
The notebook should end with a discussion of the relevant results. Depending on the style of the notebook, it may be more appropriate to do result-by-result discussion. In that case, a short conclusion or prospectus should be provided.

<!--

### References
A numbered list of references used in the document should be provided. The in-text references should be hyperlinked to the bibliography, at the bottom of the document. We recommend using the ICalico Bibliography plugin, although this currently only allows the use of MLA citation. Alternatively, we recommend a numbered-list of hyperlinked references. Numbered references should following the following formats:

* **Journal Articles** 
<blockquote> <ol><li>last_1, i1; last_2, i2, ..., and last_last, ilast. (year). "Title Hyperlinked to the article." *Journal*. **Volume**: inclusive page numbers.</li>
<li>Lozupone, C.; and Knight, R. (2005). “[UniFrac: a new phylogenetic method for comparing microbial communities](http://www.ncbi.nlm.nih.gov/pubmed/16332807).” *Appl Enviro Microbiol.* **71**: 8228-8235.
</li></ol> </blockquote>

* **Book**<blockquote>
	1. last\_1, i1; last\_2, i2, ..., last_last, ilast. (Year). *Title, hyperlinked if possible*. Edition. Publisher City: Publisher. pp. pages.
	2. Zar, J. (1999) *Biostatistical Analysis*. Fourth Ed. Upper Saddle River: Prentice Hall. pp 185.</blockquote>

* **Book Chapter**
<blockquote>
<ol><li> last\_1, i1; last\_2, i2, ..., last_last, ilast. (Year). "Chapter Title with hyperlink, if available." *Italicized Book Title*. Ch chapter. Edition. Publisher City: Publisher. pp. pages.
</li><li>Cohen, J. (1988) “<a href="http://www.lrdc.pitt.edu/schneider/p2465/Readings/Cohen,%201988%20(Statistical%20Power,%20273-406).pdf">The Analysis of Variance</a>”. *Statistical Power Analysis for the Behavioral Sciences*. Ch 8. Second Ed. Hillsdale: Lawrence Erlbaum and Associates. pp. 273 - 288.
</li></ol>
</blockquote>

* **Websites**
<blockquote>
	<ol><li>last\_1, i1; last\_2, i2, ..., last_last, ilast. (Year). "Webpage name with hyperlink". *Website title*. Sponsoring organization.
	2. Sponsoring Organization. (Year). "Webpage name with hyperlink." *Web Page title*.
	2. National Park Service. (2015). “[Mammal Checklist](http://www.nps.gov/yell/learn/nature/mammalscheck.htm).” *Yellowstone National Park*.
</blockquote>

# Getting help with git

If you're new to ``git``, you'll probably find [gitref.org](http://gitref.org/) helpful.
-->
