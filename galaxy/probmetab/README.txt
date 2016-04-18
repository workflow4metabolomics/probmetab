------------------
ProbMetab Workflow
------------------

-----------
Description
-----------

This is a Galaxy Wrapper version="1.0.0" for the ProbMetab tool.

ProbMetab, an R package that promotes substantial improvement in automatic probabilistic liquid chromatography-mass spectrometry-based metabolome annotation. The inference engine core is based on a Bayesian model implemented to (i) allow diverse source of experimental data and metadata to be systematically incorporated into the model with alternative ways to calculate the likelihood function and (ii) allow sensitive selection of biologically meaningful biochemical reaction databases as Dirichlet-categorical prior distribution. Additionally, to ensure result interpretation by system biologists, we display the annotation in a network where observed mass peaks are connected if their candidate metabolites are substrate/product of known biochemical reactions. This graph can be overlaid with other graph-based analysis, such as partial correlation networks, in a visualization scheme exported to Cytoscap, , with web and stand-alone versions.

-----------
Please cite
-----------

Authors Ricardo R. Silva et al. (2013)

Bioinformatics. 2014 May 1;30(9):1336-7. doi: 10.1093/bioinformatics/btu019. Epub 2014 Jan 17. "ProbMetab: an R package for Bayesian probabilistic annotation of LC-MS-based metabolomics".

--------------
Galaxy Wrapper
--------------

ABIMS TEAM Misharl Monsoor mmonsoor@sb-roscoff.fr

Contributors Yann Guitton yann.guitton@univ-nantes.fr and Ricardo R. Silva.


------------
Requirements
------------

You need to install the R package library ProbMetab.

Download the package here:

http://labpib.fmrp.usp.br/methods/probmetab/resources/ProbMetab_1.0.tar.gz

--------------------------------------------------------------------
Instructions for integration of the this tool into the workflow-system
Galaxy (http://getgalaxy.org)
--------------------------------------------------------------------

For installing the tool into your Galaxy installation, please do the following:

        - First of all, you need to have an account on the ABIMS server. Follow the link to abims website, where it is explained how to get an account on W4M:
http://abims.sb-roscoff.fr/faq

        - For a manual installation, copy the image from the  static/images/ into the  "galaxy-dist/static/images/"



Last but not least, restart Galaxy.


