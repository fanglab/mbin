=====
Usage
=====

Extract control IPDs from WGA sequencing with *buildcontrols*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. code-block:: console

	Usage: buildcontrols [--help] [options] wga_aligned_reads.cmp.h5

	Example:

	buildcontrols -i --procs=4 --control_pkl_name=control_means.pkl wga_aligned_reads.cmp.h5

	Options:
	  -h, --help                              Show this help message and exit
	  -d, --debug                             Increase verbosity of logging
	  -i, --info                              Add basic logging
	  --logFile=LOGFILE                       Write logging to file [log.controls]
	  --subreadlength_min=SUBREADLENGTH_MIN   Minimum subread length to include for analysis [100]
	  --readlength_min=READLENGTH_MIN         Minimum read length to include for analysis [100]
	  --min_kmer=MIN_KMER                     Minimum motif size to scan (contiguous motifs) [4]
	  --max_kmer=MAX_KMER                     Maximum motif size to scan (contiguous motifs) [6]
	  --no_bipartite                          Omit bipartite motifs [False]
	  --mod_bases=MOD_BASES                   String containing bases to query for mods. Changing this is not recommended ['A']
	  --minAcc=MINACC                         Min subread accuracy of read [0.8]
	  --minMapQV=MINMAPQV                     Min mapping QV of aligned read [240]
	  --procs=PROCS                           Number of cores to use [4]
	  --N_reads=N_READS                       Number of qualifying reads to include in analysis [1000000000]
	  --min_motif_count=MIN_MOTIF_COUNT       Number of motif sites required in WGA data to be included in controls dictionary [10]
	  --control_pkl_name=CONTROL_PKL_NAME     Filename to save control IPD data from WGA sequencing [control_ipds.pkl]

All mBin analyses require a one-time step of creating a set of control IPD values using SMRT data from whole-genome
amplified (WGA) sequencing. This WGA sequencing can be obtained from any bacterial genomic sequencing and does not have to be metagenomic. 

The control (unmethylated) IPD values for all motifs that will be queried during the process of motif discovery and methylation profile construction. The IPD for each unmethylated motif is very dependant on the sequencing chemistry used and therefore best results are obtained by ensuring that the same chemistry kit (e.g. P6-C4) is used for both the WGA and native DNA sequencing runs.

This procedure collects control IPD information for all possible motifs in the defined query space as specified with the options ``--min_kmer``, ``--max_kmer``, ``--bipart_first``, ``--bipart_Ns``, and ``--bipart_second``. Bipartite motifs describe those motifs (e.g. ACCTNNNNNCTT) that begin with a set of determinate bases (ACCT), followed by a set of indeterminate bases (NNNNN), followed by a second set of determinate bases (CTT). Accomodating every single possible bipartite motif configuration would cause the query motif space to balloon exponentially to a size that is not feasible for analysis. Therefore, constraints are placed on the acceptable lengths of each of these individual components of a bipartite motif using the options ``--bipart_first``, ``--bipart_Ns``, and ``--bipart_second``.

This is a time- and resource-intensive process, but can be omitted for subsequent analyses once the control IPD values are extracted. However, the control IPD values should be re-extracted from WGA sequencing data whenever the sequencing chemistry is updated.

A WGA aligned_reads.cmp.h5 file containing SMRT read alignments is required for this step.

.. code-block:: console
	
	$ buildcontrols -i --procs=4 --control_pkl_name=control_means.pkl wga_aligned_reads.cmp.h5

This process can be expedited in three ways:

1) Increasing ``--procs``
2) Using ``--N_reads`` to only use a subset of the aligned reads in the WGA aligned_reads.cmp.h5 file, rather than using **all** aligned reads.
3) Using ``--no_bipartite`` to omit bipartite motifs (e.g. AGCNNNNNNGTCT) from the control dictionary, only creating control values for contiguous motifs (e.g. CTGCAG). However, this will serious hamper the ability to assess the full richness of methylated motifs in the native sequencing data.

The output of this step is a `pickled <https://docs.python.org/2/library/pickle.html>`__ file (called *control_means.pkl* in the above example) containing a dictionary of the mean IPD values for all queried motifs in the WGA data. These can be manually inspected in a Python interactive shell using the following commands:

.. code-block:: console
	
	>import pickle
	>control_means = pickle.load(open("control_means.pkl", "rb"))
	>for motif,ipd in control_means.iteritems():
	>    print motif, ipd

Which should produce something similar to the following:

.. code-block:: console

	TAAGGA-5 0.22476
	GGCAAG-4 -1.08934782609
	ATANNNNNNTGCA-2 -0.306076923077
	CTGATC-3 0.344641025641
	ATANNNNNNTGCA-0 1.19992307692
	ATTCGG-0 0.000526315789474
	GTCTA-4 0.151090909091
	.....

It is important that the WGA sequencing data used for this step is of at least moderate depth and sequence complexity in order to provide sufficient control data points across the full spectrum of possible motifs. In subsequent analyses, any motifs lacking control IPD values will be discarded from the analysis, so try to include all motifs in the control data if possible.

Detect methylated motifs with *filtermotifs*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. code-block:: console

	Usage: filtermotifs [--help] [options] aligned_reads.cmp.h5

	Examples:

	Using a cmp.h5 file of aligned reads as input (recommended):
	filtermotifs -i --procs=4 --contigs=reference.fasta --control_pkl_name=control_means.pkl aligned_reads.cmp.h5

	Using a bas.h5 file of unaligned reads as input (not recommended):
	filtermotifs -i --procs=4 --control_pkl_name=control_means.pkl m12345.bas.h5

	Using a FOFN file of containing multiple files of bas.h5 unaligned reads as input (not recommended):
	filtermotifs -i --procs=4 --control_pkl_name=control_means.pkl bas.h5.fofn

	Options:
	  -h, --help                             Show this help message and exit
	  -d, --debug                            Increase verbosity of logging
	  -i, --info                             Add basic logging
	  --logFile=LOGFILE                      Write logging to file [log.controls]
	  --procs=PROCS                          Number of cores to use [4]
	  --contigs=CONTIGS                      Fasta file containing entries for the assembled contigs [None]
	  --control_pkl_name=CONTROL_PKL_NAME    Filename of control IPD data from WGA sequencing, generated using buildcontrols 
	                                         [control_ipds.pkl]
	  --motifs_fn=MOTIFS_FN                  Filename to save output filtered motifs [motifs.txt]
	  --N_reads=N_READS                      Number of reads to include for motif filtering [20000]
	  --tmp=TMP                              Directory where numerous temporary files will be written [filter_tmp]
	  --minAcc=MINACC                        Min subread accuracy of read [0.8]
	  --minMapQV=MINMAPQV                    Min mapping QV of aligned read [240]
	  --minReadScore=MINREADSCORE            Min read score of an unaligned read [0.0]
	  --maxPausiness=MAXPAUSINESS            Max pausiness value of an unaligned read [1000]
	  --subreadlength_min=SUBREADLENGTH_MIN  Minimum subread length to include for analysis [100]
	  --readlength_min=READLENGTH_MIN        Minimum read length to include for analysis [100]
	  --readlength_max=READLENGTH_MAX        Maximum read length to include for analysis [10000000]
	  --minQV=MINQV                          If base has QV < minQV, do not include [0]
	  --min_kmer=MIN_KMER                    Minimum motif size to scan (contiguous motifs) [4]
	  --max_kmer=MAX_KMER                    Maximum motif size to scan (contiguous motifs) [6]
	  --no_bipartite                         Omit bipartite motifs [False]
	  --bipart_first=BIPART_FIRST            Bipartite motif configuration: acceptable length of first determinate component 
	                                         (comma-separated string of integers) [3,4]
	  --bipart_Ns=BIPART_NS                  Bipartite motif configuration: acceptable length of middle indeterminate component 
	                                         (comma-separated string of integers) [5,6]
	  --bipart_second=BIPART_SECOND          Bipartite motif configuration: acceptable length of second determinate component 
	                                         (comma-separated string of integers) [3,4]
	  --mod_bases=MOD_BASES                  String containing bases to query for mods ['A']
	  --minMotifIPD=MINMOTIFIPD              Min motif contig IPD for inclusion of motif in final set [1.7]
	  --min_motif_reads=MIN_MOTIF_READS      Min number of reads with motif hits to keep for motif filtering (only if using 
	                                         unaligned reads as input) [20]
	  --min_motif_N=MIN_MOTIF_N              Min number of motif IPD values required to keep for motif filtering [20]
	  --cross_cov_bins=CROSS_COV_BINS        Path to file containing binning results from CONCOCT. Will use to improve motif 
	                                         discovery. Only works with contig-level analysis (cmp.h5 input) inputs. File format 
	                                         should be '<contig_name>,<bin_id>' [None]

After the control IPD values have been tabulated and stored, methylated motifs can then detected in the HDF5 files of native, metagenomic sequencing data. Both aligned reads (cmp.h5) and unaligned reads (bas.h5 or FOFN containing multiple bas.h5 files) are supported as input, but the use of aligned reads is strongly recommended for motif filtering. Unaligned reads contain significant sequencing errors that introduce noise into the IPD signals for motifs, making it difficult to detect truly methylated motifs from unaligned reads.

.. code-block:: console

	filtermotifs -i --procs=4 --contigs=reference.fasta --control_pkl_name=control_means.pkl native_aligned_reads.cmp.h5

Unless otherwise specified (using ``--min_kmer``, ``--max_kmer``, ``--bipart_first``, ``--bipart_Ns``, ``--bipart_second``, or ``--no_bipartite`` options), this procedure starts with motifs for which control data exists (``--control_pkl_name``) and discards all motifs that do not have a methylation score (native IPD - control IPD) greater than the value defined by ``--minMotifIPD``.

The original query space of motifs is very large and can include several hundred thousand motifs if bipartite motifs are included (recommended). To ease computational demands for storing IPD data for all motifs in the query space, only a subset of reads (``--N_reads``) are examined from the input HDF5 files. This value of ``--N_reads`` can be modified according to available computational resources and acceptable running time.

The filtered motifs are written to the output file specified by ``--motifs_fn``. Motifs are listed with the sequence string, followed by the 0-based index of the methylated base. For example, GATC-1 indicates that the A position in the motif is methylated. This ``--motifs_fn`` output file serves as input to *methylprofiles*, which constructs methylation profiles across the filtered motifs.


Build methylation profiles with *methylprofiles*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. code-block:: console

	Usage: methylprofiles [--help] [options] input_hdf5 motifs.txt

	Example:

	Using a cmp.h5 file of aligned reads as input:
	methylprofiles -i --procs=4 --control_pkl_name=control_means.pkl --contigs=reference.fasta aligned_reads.cmp.h5 motifs.txt

	Using a bas.h5 file of unaligned reads as input:
	methylprofiles -i --procs=4 --control_pkl_name=control_means.pkl m12345.bas.h5

	Using a FOFN file of containing multiple files of bas.h5 unaligned reads as input:
	methylprofiles -i --procs=4 --control_pkl_name=control_means.pkl bas.h5.fofn

	Options:
	  -h, --help                              show this help message and exit
	  -d, --debug                             Increase verbosity of logging
	  -i, --info                              Add basic logging
	  --logFile=LOGFILE                       Write logging to file [log.controls]
	  --prefix=PREFIX                         Prefix to use for output files [None]
	  --tmp=TMP                               Directory where numerous temporary files will be written [profiles_tmp]
	  --contigs=CONTIGS                       Fasta file containing entries for the assembled contigs [None]
	  --minReadScore=MINREADSCORE             Min read score of an unaligned read [0.0]
	  --maxPausiness=MAXPAUSINESS             Max pausiness value of an unaligned read [1000]
	  --minQV=MINQV                           If base has QV < minQV, do not include [0]
	  --subreadlength_min=SUBREADLENGTH_MIN   Minimum subread length to include for analysis [100]
	  --readlength_min=READLENGTH_MIN         Minimum read length to include for analysis [100]
	  --minContigLength=MINCONTIGLENGTH       Min length of contig to consider [10000]
	  --comp_kmer=COMP_KMER                   Kmer size to use for sequence composition measurements [5]
	  --aligned_read_barcodes                 Also output features for individual aligned reads, not just 
	                                          contigs (requires cmp.h5 input) [False]
	  --minAcc=MINACC                         Min subread accuracy of read [0.8]
	  --minMapQV=MINMAPQV                     Min mapping QV of aligned read [240]
	  --procs=PROCS                           Number of processors to use [4]
	  --N_reads=N_READS                       Number of qualifying reads to include (from each bas.h5 if input is FOFN 
	                                          of bas.h5 files) in analysis [1000000000]
	  --control_pkl_name=CONTROL_PKL_NAME     Filename to save control IPD data from WGA sequencing [control_ipds.pkl]
	  --subtract_control=SUBTRACT_CONTROL     Subtract control IPDs in final calculations [True]
	  --cross_cov_bins=CROSS_COV_BINS         Path to file containing binning results from CONCOCT. Will use to improve 
	                                          motif discovery. Only works with contig-level analysis (cmp.h5 input) inputs. 
	                                          File format should be '<contig_name>,<bin_id>' [None]

*methylprofiles* compiles the methylation profiles across the motifs specified in the arguments. The methylation profiles can be constructed using either native contigs (\*.cmp.h5) or unaligned reads (\*.bas.h5), the latter of which can be supplied as a single bas.h5 file or a FOFN containing multiple \*.bas.h5 files (each *bas.h5 file on a new line). 

The output consists of three separate files containing methylation features, as well as other relevant features for binning:

1. **Methylation features**: <prefix>_<seq>_methyl_features.txt
	* Column 1: <seq> id
	* Column 2: <seq> length
	* Columns 3-M: Methylation scores (for M motifs)
2. **Motif counts**: <prefix>_<seq>_motif_counts.txt
	* Column 1: <seq> id
	* Column 2: <seq> length
	* Columns 3-M: Motif counts (for M motifs)
3. **Alternative features**: <prefix>_<seq>_other_features.txt
	* Column 1: <seq> id
	* Column 2: <seq> length
	* For dtype = read or align: 
		* Columns 3-N: k-mer frequencies (for N k-mers)
	* For dtype = contig: 
		* Column 3: Contig coverage
		* Columns 4-N: k-mer frequencies (for N k-mers)

Where <prefix> is defined by ``--prefix`` and <seq> is the sequence data type: contig, align, or read. When inputting \*.cmp.h5 reads for methylation profiling, *methylprofiles* will always generate contig level features and will optionally generate align level features for the reads comprising each contig (using ``--aligned_read_barcodes``). When \*.bas.h5 reads are input, only read level features will be output.

For example, the following command with generate methylation (and other) profiles for a set of contigs contained in aligned_reads.cmp.h5:

.. code-block:: console
	
	methylprofiles -i --procs=4 --prefix=test --control_pkl_name=control_means.pkl --contigs=reference.fasta aligned_reads.cmp.h5 motifs.txt

These profiles will be contained in the following output files:

1. test_contig_methyl_features.txt
2. test_contig_motif_counts.txt
3. test_contig_other_features.txt

Visualize feature landscape with *mapfeatures*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: console

	Usage: mapfeatures [--help] [options] <SEQ>_methyl_features.txt <SEQ>_other_features.txt

	Examples:

	mapfeatures -i --labels=<LABELS.txt> --size_markers contig_methyl_features.txt contig_other_features.txt

	mapfeatures -i --labels=<LABELS.txt> --l_min=500 align_methyl_features.txt align_other_features.txt

	mapfeatures -i --l_min=10000 --n_seqs=1000 read_methyl_features.txt read_other_features.txt

	Options:
	  -h, --help                Show this help message and exit
	  -d, --debug               Increase verbosity of logging
	  -i, --info                Add basic logging
	  --logFile=LOGFILE         Write logging to file [log.controls]
	  --prefix=PREFIX           Prefix to use for output files [None]
	  --size_markers            Adjust marker size in plot according to sequence length [False]
	  --dim_reduce=DIM_REDUCE   Dimensionality reduction algorithm to apply (bhtsne or pca) [bhtsne]
	  --labels=LABELS           Tab-delimited file (no header) of sequence labels (seq_name\tlabel_name) [None]
	  --l_min=L_MIN             Minimum read length to include for analysis [0]
	  --n_seqs=N_SEQS           Number of sequences to subsample [all]
	  --n_dims=N_DIMS           Number of dimensions to reduce to for visualization (only n_dims=2 will be plotted) [2]
	  --n_iters=N_ITERS         Number of iterations to use for BH-tSNE [500]

*mapfeatures* visualizes the landscape of high-dimensional sequence features using the Barnes Hut approximation of t-SNE (*PCA support coming soon*). The sequence features that are output from *methylprofiles* are often high-dimensional (>3D), making it difficult to visualize the sequences. To ease this visualization for resolution of discrete sequence clusters in the feature space, t-SNE is used to reduce the dimensionality of the methylation, composition, and coverage features to 2D. The resulting 2D maps, which can be overlaid with sequence annotation labels (generated with Kraken_, for instance), often reveals sequence clustering in the 2D feature space representing distinct taxonomical groups for binning.

.. _Kraken: https://ccb.jhu.edu/software/kraken/

Two files from *methylprofiles* serve as input to *mapfeatures*: ``<prefix>_<seq>_methyl_features.txt`` and ``<prefix>_<seq>_other_features.txt``, which contain the high-dimensional methylation score features and coverage & composition features, respectively. As before, <prefix> is defined by ``--prefix`` and <seq> is the sequence data type: contig, align, or read.

For example, the command:

.. code-block:: console

	mapfeatures -i --labels=<LABELS.txt> --size_markers contig_methyl_features.txt contig_other_features.txt

Will generate three 2D t-SNE maps representing the contig feature space:

1. **Methylation t-SNE map**: contig_methyl_features.bhtsne.png
2. **Composition t-SNE map**: contig_other_features.comp.bhtsne.png
3. **Coverage+composition t-SNE map**: contig_other_features.covcomp.bhtsne.png

Each file has an accompanying tab-delimited \*.txt file containing the 2D coordinates from each t-SNE map that can be used for additional plotting purposes by the user. The file supplied to ``--labels`` must contain tab-delimited sequence labels for **all** sequences listed in the first column of the <seq>_methyl_features.txt and <seq>_other_features.txt files that are output from *methylprofiles*. Any sequences lacking annotation must be included using the label 'Unlabeled', for example:

+------------+-------------+
| contig_1   | Bacteroides |
+------------+-------------+
| contig_2   | Clostridium |
+------------+-------------+
| contig_4   | Escherichia |
+------------+-------------+
| contig_5   | Bacteroides |
+------------+-------------+
| contig_6   | Unlabeled   |
+------------+-------------+
