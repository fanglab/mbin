=====
Usage
=====

Extract control IPDs from WGA sequencing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
All mBin analyses require a one-time step of creating a set of control IPD values using SMRT data from whole-genome
amplified (WGA) sequencing. This consists of control (unmethylated) IPD values for all motifs that will be queried during the process of motif discovery and methylation profile construction. The IPD for each unmethylated motif is very dependant on the sequencing chemistry used and therefore best results are obtained by ensuring that the same chemistry kit (e.g. P6-C4) is used for both the WGA and native DNA sequencing runs.

This is a time- and resource-intensive process, but can be omitted for subsequent analyses once the control IPD values are extracted. However, the control IPD values should be re-extracted from WGA sequencing data whenever the sequencing chemistry is updated.

IN PROGRESS

Detect methylated motifs
^^^^^^^^^^^^^^^^^^^^^^^^^^^
IN PROGRESS

Build methylation profiles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
IN PROGRESS

Visualize results
^^^^^^^^^^^^^^^^^^^
IN PROGRESS
