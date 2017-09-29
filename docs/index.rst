mBin documentation
==================
mBin: a methylation-based binning framework for metagenomic SMRT sequencing reads

The mBin pipeline is designed to discover the unique signals of DNA methylation in metagenomic SMRT sequencing reads and leverage them for organism binning of assembled contigs or unassembled reads.


* Free software: BSD license
* Documentation: https://mbin.readthedocs.io. <-- IN PROGRESS
* Please report any issues through the issue tracker on `GitHub <https://github.com/fanglab/mbin/issues>`__

For basic usage instructions and list of optional parameters, run ``mbin --help``

Features
^^^^^^^^
mBin can take as input either unaligned single molecule real-time (SMRT) reads from a PacBio instrument or contigs assembled from SMRT reads. Methylation scores are calculated from individual inter-pulse duration (IPD) metrics embedded in each sequencing read that record the polymerase kinetics during sequencing and indicate the presence or absence of DNA methylation at the level of individual nucleotides.

By aggregating these IPD metrics across multiple sites for a motif and across multiple reads aligned to a contig, mBin generates methylation scores for motifs and uses these to construct methylation profiles for reads and contigs. Methylation profiles can then be used as epigenetic features for unsupervised metagenomic binning. mBin can also generate methylation scores for contigs that are given binning assignments by other binning tools (with the --cross_cov_bins option).


Installation
^^^^^^^^^^^^
For a comprehensive guide on how to install mBin and its dependencies, see :doc:`installation`


Contribute
^^^^^^^^^^
* Issue tracker: `GitHub <https://github.com/fanglab/mbin/issues>`__
* Source code: `GitHub <https://github.com/fanglab/mbin>`__


Contents
^^^^^^^^
.. toctree::
   :maxdepth: 2

   readme
   installation
   usage
   modules
   contributing

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
