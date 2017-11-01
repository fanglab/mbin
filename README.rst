=============
mBin overview
=============


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

The mBin pipeline is designed to discover the unique signals of DNA methylation in metagenomic SMRT sequencing reads and leverage them for organism binning of assembled contigs or unassembled reads. Because *all* cellular DNA is modified by the same methyltransferase, DNA methylation signals can be used for binning not just chromosomal sequences, but also extrachromosomal mobile genetic elements like plasmids.

The pipeline consists of four routines:

1. *buildcontrols*: Gets unmethylated IPD values for motifs from whole-genome amplified (WGA) sequencing 
2. *filtermotifs*: Identifies methylated motifs in native metagenomic sequencing
3. *methylprofiles*: Creates methylation profiles for sequences using speified motifs 
4. *mapfeatures*: Visualizes landscape of methylation features across all sequences

Documentation
-------------
For a comprehensive guide on how to install and run mBin, please see the full `documentation <https://mbin.readthedocs.io/en/latest/>`__.


Citations
---------
Beaulaurier J, Zhu S, Deikus G, Mogno I, Zhang XS, Davis-Richardson A, Canepa R, Triplett E, Faith J, Sebra R, Schadt EE. Metagenomic binning and association of plasmids with bacterial host genomes using DNA methylation. *In press at Nature Biotechnology*.

Credits
---------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage

