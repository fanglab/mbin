====
mBin
====


.. image:: https://img.shields.io/pypi/v/mbin.svg
        :target: https://pypi.python.org/pypi/mbin

.. image:: https://img.shields.io/travis/fanglab/mbin.svg
        :target: https://travis-ci.org/fanglab/mbin

.. image:: https://readthedocs.org/projects/mbin/badge/?version=latest
        :target: https://mbin.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://pyup.io/repos/github/fanglab/mbin/shield.svg
     :target: https://pyup.io/repos/github/fanglab/mbin/
     :alt: Updates


mBin: a methylation-based binning framework for metagenomic SMRT sequencing reads

The mBin pipeline is designed to discover the unique signals of DNA methylation in metagenomic SMRT sequencing reads and leverage them for organism binning of assembled contigs or unassembled reads.


* Free software: BSD license
* Documentation: https://mbin.readthedocs.io. <-- IN PROGRESS

For basic usage instructions, run "mbin --help"


Features
--------
mBin can take as input either unaligned single molecule real-time (SMRT) reads from a PacBio instrument or contigs assembled from SMRT reads. Methylation scores are calculated from individual inter-pulse duration (IPD) metrics embedded in each sequencing read that record the polymerase kinetics during sequencing and indicate the presence or absence of DNA methylation at the level of individual nucleotides.

By aggregating these IPD metrics across multiple sites for a motif and across multiple reads aligned to a contig, mBin generates methylation scores for motifs and uses these to construct methylation profiles for reads and contigs. Methylation profiles can then be used as epigenetic features for unsupervised metagenomic binning. mBin can also generate methylation scores for contigs that are given binning assignments by other binning tools (with the --cross_cov_bins option).

* TODO


Installation
---------
For a comprehensive guide on how to install mBin and its dependencies, see :doc:`installation`


Citations
---------
TBA

Credits
---------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage

