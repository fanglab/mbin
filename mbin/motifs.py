import os,sys
import optparse
import logging
from pbcore.io.align.CmpH5IO import CmpH5Reader
from pbcore.io.BasH5IO import BasH5Reader
import glob
import numpy as np
import logging
import shutil
import pickle
import math
import mbin
import motif_tools
from Bio import SeqIO
import multiprocessing
import subprocess
import read_scanner
from collections import defaultdict
import copy
from itertools import izip
import unicodedata
import re

def launch():
	opts,h5_files = __parseArgs()
	__initLog(opts)

	opts = get_control_motifs(opts)

	filter_motifs(opts, h5_files)

	print >> sys.stderr, "mBin motif filtering has finished running. See log for details."

def filter_motifs(opts, h5_files):
	"""
	Samples N=<N_reads> reads and looks for evidence of methylation
	in all the motifs in the query space.
	"""
	if os.path.exists(opts.tmp):
		shutil.rmtree(opts.tmp)
	os.mkdir(opts.tmp)

	##########################################################
	# Check which motifs are in the dictionary of control IPDs
	##########################################################
	logging.info("")
	logging.info("Assessing the motifs for which we have control IPD values...")
	controls_d = pickle.load(open(opts.control_pkl_name, "r"))
	all_m      = controls_d.keys()
	
	logging.info(" -- Contiguous (e.g. CATG):")
	for l in range(1,20):
		hits = len([m for m in all_m if len(m.split("-")[0])==l and m.find("N")<0])
		if hits>0:
			logging.info("    Found %s motifs of length %s" % (hits,l))
	
	logging.info(" -- Bipartite (e.g. ACGNNNNNCTT):")
	bi_hits = len([m for m in all_m if m.find("N")>-1])
	logging.info("    Found %s motifs" % bi_hits)
	logging.info("Found %s total motifs" % len(all_m))
	logging.info("")

	mbinRunner     = mbin.mbinRunner(opts)
	##################################################
	# Launch analysis of <N_reads> for motif filtering
	##################################################
	for i,h5_file in enumerate(h5_files):
		logging.info("Creating %s barcodes (%s motifs) from %s..." % (opts.N_reads, (len(opts.motifs)+len(opts.bi_motifs)), h5_file))
		mbinRunner.launch_data_loader( h5_file, opts.N_reads, i, opts )

	if opts.h5_type=="bas":
		# Combine subread data across multiple movies
		logging.info("Combining subread data across all movies...")
		results = mbinRunner.combine_subread_data_across_bas_movies()
		logging.info("Done.")
		# Combine movie-merged subreads to get read-level barcodes
		logging.info("Combining subreads to get read-level barcodes...")
		results = mbinRunner.bas_combine_subreads_for_read_level( tmp_run=True )
		logging.info("Done.")

	filter_runner = FilterRunner(opts, h5_files)
	filter_runner.run(mbinRunner)

	logging.info("Cleaning up temp files from motif filtering...")
	shutil.rmtree(opts.tmp)

def get_motif_scores_from_read( args ):
	j                   = args[0]
	control_means       = args[1]
	reads_ipds_fn       = args[2]
	reads_ipds_N_fn     = args[3]
	reads_ipds_kmers_fn = args[4]
	opts                = args[5]
	
	incr                = 1000
	if j%incr==0:
		logging.info("   Read %s" % j)
	kmers               = np.loadtxt(reads_ipds_kmers_fn, dtype="str")
	motifs_max_score    = {}
	files               = [reads_ipds_fn    +".trans", reads_ipds_N_fn  +".trans"]
	for i,(l1,l2) in enumerate(izip( open(files[0]), open(files[1]))):
		if i==j:
			motif = kmers[i]
			if not motifs_max_score.get(motif):
				motifs_max_score[motif] = -10
			n_reads_passing = 0
			for read_i in range(len(l1.split())):
				case_read_N = int(l2.split()[read_i])
				if case_read_N>=opts.min_motif_N:
					case_mean            = float(l1.split()[read_i])
					score                = case_mean - control_means[motif]
					if score > opts.minMotifIPD:
						n_reads_passing += 1
	if n_reads_passing >= opts.min_motif_reads:
		return motif

def process_contig_chunk( args ):
	chunk_id      = args[0]
	control_pkl   = args[1]
	cut_CMDs      = args[2]
	kmers         = args[3]
	cols_chunk    = args[4]
	contig_id     = args[5]
	n_chunks      = args[6]
	n_contigs     = args[7]
	opts          = args[8]
	logging.info("  - Contig %s/%s: chunk %s/%s" % ((contig_id+1), n_contigs, (chunk_id+1), (n_chunks+1)))
	control_means = pickle.load(open(control_pkl, "rb"))
	contig_motifs = {}
	case_motif_Ns = {}
	
	for cut_CMD in cut_CMDs:
		sts,stdOutErr = mbin.run_OS_command( cut_CMD )
	
	fns                = map(lambda x: x.split("> ")[-1], cut_CMDs)
	contig_ipds_sub    = np.loadtxt(fns[0], dtype="float")
	contig_ipds_N_sub  = np.loadtxt(fns[1], dtype="int")
	# If there is only one row (read) for this contig, still treat as
	# a 2d matrix of many reads
	contig_ipds_sub    = np.atleast_2d(contig_ipds_sub)
	contig_ipds_N_sub  = np.atleast_2d(contig_ipds_N_sub)
	for j in range(len(cols_chunk)):
		motif = kmers[cols_chunk[j]]
		case_contig_Ns    = contig_ipds_N_sub[:,j]
		if control_means.get(motif):
			case_contig_means = contig_ipds_sub[:,j]
			if np.sum(case_contig_Ns)>0:
				case_mean = np.dot(case_contig_means, case_contig_Ns) / np.sum(case_contig_Ns)
			else:
				case_mean = 0
			score                = case_mean - control_means[motif]
			contig_motifs[motif] = score
			case_motif_Ns[motif] = np.sum(case_contig_Ns)
	return contig_motifs,case_motif_Ns

def chunk_case_control_files( opts, contig, j, contigs_N ):
	logging.info("   ...chunking contig %s (%s/%s)..." % (contig,(j+1),contigs_N))

	contig_SCp           = {}
	contig_SCp_N         = {}
	keeper_motifs        = set()
	control_means        = pickle.load(open(opts.control_pkl_name, "rb"))
	contig_ipds_fn       = os.path.join( opts.tmp, "%s_ipds.tmp"      % contig)
	contig_ipds_N_fn     = os.path.join( opts.tmp, "%s_ipdsN.tmp"     % contig)
	contig_ipds_kmers_fn = os.path.join( opts.tmp, "%s_ipdskmers.tmp" % contig)
	kmers                = np.loadtxt(contig_ipds_kmers_fn, dtype="str")
	fns                  = [contig_ipds_fn, contig_ipds_N_fn]
	
	n_chunks    = 99
	chunksize   = int(math.ceil(float( len(kmers)/n_chunks )))
	cols_chunks = list(chunks( range(len(kmers)), chunksize ))
	args        = []
	for i,cols_chunk in enumerate(cols_chunks):
		cut_CMDs = []
		for fn in fns:
			cut_cols = "%s-%s" % ((cols_chunk[0]+1), (cols_chunk[-1]+1))
			in_fn    = fn
			out_fn   = fn+".sub.%s" % i
			cut_CMD  = "cut -d$\'\\t\' -f%s %s > %s" % (cut_cols, in_fn, out_fn)
			cut_CMDs.append(cut_CMD)
		args.append( (i, opts.control_pkl_name, cut_CMDs, kmers, cols_chunk, j, n_chunks, contigs_N, opts) )
	results = mbin.launch_pool(opts.procs, process_contig_chunk, args)
	for i,result in enumerate(results):
		for motif in result[0].keys(): # contig_SCp,contig_SCp_N
			contig_SCp[motif]   = result[0][motif]
			contig_SCp_N[motif] = result[1][motif]

	return control_means, contig_SCp, contig_SCp_N, contig

def get_control_motifs(opts):
	"""
	We can only filter motifs from the set of motifs for which we
	have control values, so see what motifs are present in the 
	control dictionary.
	"""
	control_means  = pickle.load(open(opts.control_pkl_name, "rb"))
	opts.motifs    = set([m for m in control_means.keys() if m.find("N")<0])
	opts.bi_motifs = set([m for m in control_means.keys() if m.find("N")>0])

	return opts

def transpose_file( fn ):
	def run_OS( CMD ):
		p         = subprocess.Popen(CMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdOutErr = p.communicate()
		sts       = p.returncode
		if sts != 0:
			logging.warning("Failed command: %s" % CMD)
		return sts, stdOutErr

	transpose_dir  = os.path.dirname(os.path.realpath(read_scanner.__file__))
	trans_script   = os.path.join(transpose_dir, "transpose.awk")
	trans_CMD      = "awk -f %s %s > %s.trans" % (trans_script, fn, fn)
	logging.info(trans_CMD)
	sts, stdOutErr = run_OS(trans_CMD)

def transpose_contig_matrix( args ):
	contig  = args[0]
	opts    = args[1]
	logging.info("  Transposing %s" % contig)
	contig_ipds_fn       = os.path.join( opts.tmp, "%s_ipds.tmp"       % contig)
	contig_ipds_kmers_fn = os.path.join( opts.tmp, "%s_ipdskmers.tmp"  % contig)
	contig_ipds_N_fn     = os.path.join( opts.tmp, "%s_ipdsN.tmp"      % contig)
	contig_ipds          = np.loadtxt(contig_ipds_fn,       dtype="float")
	contig_ipds_kmers    = np.loadtxt(contig_ipds_kmers_fn, dtype="str")
	contig_ipds_N        = np.loadtxt(contig_ipds_N_fn,     dtype="int")
	if len(contig_ipds.shape)==1:
		contig_ipds   = contig_ipds.reshape(1,contig_ipds.shape[0])
		contig_ipds_N = contig_ipds_N.reshape(1,contig_ipds_N.shape[0])

	contig_ipds    = contig_ipds.T
	contig_ipds_N  = contig_ipds_N.T
	np.savetxt(contig_ipds_fn+".trans",   contig_ipds,   fmt="%.4f", delimiter="\t")
	np.savetxt(contig_ipds_N_fn+".trans", contig_ipds_N, fmt="%s",   delimiter="\t")
	return None

def build_bin_dict( contig, bin_id, contig_SCp_dict, bin_SCp_dict ):
	"""
	For each contig, incorporate its methylation score values into 
	the running bin methylation score values.
	"""
	logging.info("Adding methylation scores from contig %s to bin %s..." % (contig, bin_id))
	for motif in contig_SCp_dict["SCp"].keys():
		"""
		If the bin doesn't yet have values for a motif, initialize
		them with zeroes.
		"""
		if not bin_SCp_dict["SCp"].get(motif):
			bin_SCp_dict["SCp"][motif]   = 0.0
			bin_SCp_dict["SCp_N"][motif] = 0

		# Compute contig and bin score sums for a motif
		contig_sum    = contig_SCp_dict["SCp"][motif] * contig_SCp_dict["SCp_N"][motif]
		bin_sum       = bin_SCp_dict["SCp"][motif]    * bin_SCp_dict["SCp_N"][motif]
		
		# Get new sum for the bin and divide by new count to get new score
		new_bin_sum   = bin_sum + contig_sum
		new_bin_N     = bin_SCp_dict["SCp_N"][motif] + contig_SCp_dict["SCp_N"][motif]
		
		if new_bin_N>0:
			new_bin_score = new_bin_sum / new_bin_N
		else:
			new_bin_score = 0

		bin_SCp_dict["SCp"][motif]   = new_bin_score
		bin_SCp_dict["SCp_N"][motif] = new_bin_N

	return bin_SCp_dict

def stream_case_control_files( tup ):
	opts                 = tup[0]
	contig               = tup[1]
	j                    = tup[2]
	contigs_N            = tup[3]

	logging.info("   ...streaming contig %s (%s/%s)..." % (contig,(j+1),contigs_N) )

	control_means        = pickle.load(open(opts.control_pkl_name, "rb"))
	contig_SCp           = {}
	contig_SCp_N         = {}
	keeper_motifs        = set()
	contig_ipds_fn       = os.path.join( opts.tmp, "%s_ipds.tmp"      % contig)
	contig_ipds_kmers_fn = os.path.join( opts.tmp, "%s_ipdskmers.tmp" % contig)
	contig_ipds_N_fn     = os.path.join( opts.tmp, "%s_ipdsN.tmp"     % contig)
	kmers = np.loadtxt(contig_ipds_kmers_fn, dtype="str")
	files = [contig_ipds_fn+".trans", contig_ipds_N_fn+".trans"]
	for i,(l1,l2) in enumerate(izip( open(files[0]), open(files[1]))):
		motif          = kmers[i]
		case_contig_Ns = map(lambda x: int(x), l2.split())
		if control_means.get(motif):
			case_contig_means    = map(lambda x: float(x), l1.split())
			if np.sum(case_contig_Ns)>0:
				case_mean = np.dot(case_contig_means, case_contig_Ns) / np.sum(case_contig_Ns)
			else:
				case_mean = 0
			score                = case_mean - control_means[motif]
			contig_SCp[motif]    = score
			contig_SCp_N[motif]  = np.sum(case_contig_Ns)
	
	return contig_SCp, contig_SCp_N, contig

def chunks( l, n ):
	"""
	Yield successive n-sized chunks from l.
	"""
	for i in xrange(0, len(l), n):
		yield l[i:i+n]

def simplify_motifs( tup ):
	"""
	Call functions to find the shortest version of a motif
	and look for possible degenerate versions of the motif.
	The "units" referred to are contigs by default, but are
	bins of contigs if --cross_cov_bins=<FILE> is used to 
	group contigs based on bin assignments.
	"""
	j                = tup[0]
	unit_SCp         = tup[1]
	unit_SCp_N       = tup[2]
	opts             = tup[3]
	units_N          = tup[4]
	unit_name        = tup[5]
	unit_type        = tup[6]

	# Keep only the motifs with scores higher than the minMotifIPD value
	highscore_motifs = dict( [(motif,SCp) for motif,SCp in unit_SCp.items() if SCp>=opts.minMotifIPD and unit_SCp_N[motif]>=opts.min_motif_N] )

	# Keep only the shortest version of the high scoring motifs (reduces redundancy)
	keeper_motifs = motif_tools.shorten_motifs( unit_SCp, highscore_motifs )
	logging.info("  - %s %s %s/%s after motif shortening: %s" % (unit_type, unit_name, (j+1), units_N, ";".join(list(keeper_motifs))))

	# See if any sets of motifs can be simplified to a degenerate motif
	keeper_motifs, unit_SCp, unit_SCp_N, degen_members = motif_tools.refine_degen_motifs( keeper_motifs, \
																						  unit_SCp, \
																						  unit_SCp_N )
	logging.info("  - %s %s %s/%s after degenerate motif calling: %s" % (unit_type, unit_name, (j+1), units_N, ";".join(list(keeper_motifs))))

	rpt_str = ""
	for motif in list(keeper_motifs):
		rpt_str += "%s (%.3f ; %s)," % (motif, unit_SCp[motif], unit_SCp_N[motif])
	logging.info("  - %s %s %s/%s detected motifs: %s" % (unit_type, unit_name, (j+1), units_N, rpt_str))

	# Pull the control IPD values for the keeper motifs
	control_means = pickle.load(open(opts.control_pkl_name, "rb"))
	# Add control values for new degenerate bases (mean of the specified motifs)
	for d_motif,o_motifs in degen_members.iteritems():
		degen_mean = np.mean([control_means[o_m] for o_m in o_motifs])
		control_means[d_motif] = degen_mean
	
	# Delete all control mean entries for motifs not in the keeper set
	# for motif in control_means.keys():
	# 	if motif not in keeper_motifs:
	# 		del control_means[motif]

	return keeper_motifs, control_means

def __parseArgs():
	"""Handle command line argument parsing"""

	usage = """%prog [--help] [options] input_hdf5

	Examples:

	Using a cmp.h5 file of aligned reads as input (recommended):
	filtermotifs -i --procs=4 --contigs=reference.fasta --control_pkl_name=control_means.pkl aligned_reads.cmp.h5

	Using a bas.h5 file of unaligned reads as input (not recommended):
	filtermotifs -i --procs=4 --contigs=reference.fasta --control_pkl_name=control_means.pkl m12345.bas.h5

	Using a FOFN file of containing multiple files of bas.h5 unaligned reads as input (not recommended):
	filtermotifs -i --procs=4 --contigs=reference.fasta --control_pkl_name=control_means.pkl bas.h5.fofn
	"""

	parser = optparse.OptionParser( usage=usage, description=__doc__ )

	parser.add_option( "-d", "--debug", action="store_true", help="Increase verbosity of logging" )

	parser.add_option( "-i", "--info", action="store_true", help="Add basic logging" )

	parser.add_option( "--logFile", type="str", help="Write logging to file [log.controls]" )

	parser.add_option( "--procs", type="int", help="Number of cores to use [4]" )
	
	parser.add_option( "--contigs", type="str", help="Fasta file containing entries for the assembled contigs [None]" )

	parser.add_option( "--control_pkl_name", type="str", help="Filename of control IPD data from WGA sequencing, generated using buildcontrols [control_ipds.pkl]" )
	
	parser.add_option( "--motifs_fn", type="str", help="Filename to save output filtered motifs [motifs.txt]" )
	
	parser.add_option( "--N_reads", type="int", help="Number of reads to include for motif filtering [20000]" )
	
	parser.add_option( "--tmp", type="str", help="Directory where numerous temporary files will be written [filter_tmp]" )

	parser.add_option( "--minAcc", type="float", help="Min subread accuracy of read [0.8]" )

	parser.add_option( "--minMapQV", type="float", help="Min mapping QV of aligned read [240]" )

	parser.add_option( "--minReadScore", type="float", help="Min read score of an unaligned read [0.0]" )

	parser.add_option( "--maxPausiness", type="float", help="Max pausiness value of an unaligned read [1000]" )

	parser.add_option( "--subreadlength_min", type="int", help="Minimum subread length to include for analysis [100]" )

	parser.add_option( "--readlength_min", type="int", help="Minimum read length to include for analysis [100]" )

	parser.add_option( "--readlength_max", type="int", help="Maximum read length to include for analysis [10000000]" )

	parser.add_option( "--minQV", type="float", help="If base has QV < minQV, do not include [0]" )

	parser.add_option( "--min_kmer", type="int", help="Minimum motif size to scan (contiguous motifs) [4]" )

	parser.add_option( "--max_kmer", type="int", help="Maximum motif size to scan (contiguous motifs) [6]" )

	parser.add_option( "--no_bipartite", action="store_true", help="Omit bipartite motifs [False]" )
	
	parser.add_option( "--bipart_first", type="str", help="Bipartite motif configuration: acceptable length of first determinate component (comma-separated string of integers) [3,4]" )
	
	parser.add_option( "--bipart_Ns", type="str", help="Bipartite motif configuration: acceptable length of middle indeterminate component (comma-separated string of integers) [5,6]" )
	
	parser.add_option( "--bipart_second", type="str", help="Bipartite motif configuration: acceptable length of second determinate component (comma-separated string of integers) [3,4]" )

	parser.add_option( "--mod_bases", type="str", help="String containing bases to query for mods ['A']" )
	
	parser.add_option( "--minMotifIPD", type="float", help="Min motif contig IPD for inclusion of motif in final set [1.7]" )

	parser.add_option( "--min_motif_reads", type="int", help="Min number of reads with motif hits to keep for motif filtering (only if using unaligned reads as input) [20]" )
	
	parser.add_option( "--min_motif_N", type="int", help="Min number of motif IPD values required to keep for motif filtering [20]" )
	
	parser.add_option( "--cross_cov_bins", type="str", help="Path to file containing binning results from CONCOCT. Will use to improve motif discovery. Only works with contig-level analysis (cmp.h5 input) inputs. File format should be '<contig_name>,<bin_id>' [None]" )
	
	parser.set_defaults( logFile="log.filtermotifs",           \
						 info=False,                           \
						 debug=False,                          \
						 procs=4,                              \
						 contigs=None,                         \
						 control_pkl_name="control_means.pkl", \
						 motifs_fn="motifs.txt",               \
						 N_reads=20000,                        \
						 tmp="filter_tmp",                     \
						 minAcc=0.8,                           \
						 minMapQV=240,                         \
						 minReadScore=0.0,                     \
						 maxPausiness=1000,                    \
						 subreadlength_min=100,                \
						 readlength_min=100,                   \
						 readlength_max=10000000,              \
						 minQV=0.0,                            \
						 min_kmer=4,                           \
						 max_kmer=6,                           \
						 mod_bases="A",                        \
						 no_bipartite=False,                   \
						 bipart_first="3,4",                   \
						 bipart_Ns="5,6",                      \
						 bipart_second="3,4",                  \
						 minMotifIPD=1.7,                      \
						 min_motif_N=20,                       \
						 min_motif_reads=20,                   \
						 cross_cov_bins=None,                  \
						 )

	opts, args = parser.parse_args( )

	h5_files   = __check_input( opts, args, parser )

	if opts.no_bipartite:
		opts.bipartite = False
	else:
		opts.bipartite = True

	############################################
	# Define the types of bipartite motifs to 
	# include in the analysis. This describes the
	# acceptable sizes of the three components of
	# bipartite motifs. For example, the motif 
	# ACCT/NNNNN/CTT (first/Ns/last) would be 
	# described by 4/5/3.
	
	first  = map(lambda x: int(x), opts.bipart_first.split(","))
	middle = map(lambda x: int(x), opts.bipart_Ns.split(","))
	second = map(lambda x: int(x), opts.bipart_second.split(","))
	opts.bipart_config   = [(first), (middle), (second)]
	
	# As set using default values, acceptible bipartite 
	# motifs would have the following component lengths:
	# First: 3 or 4 ACGT bases
	# Ns: 5 or 6 unspecified N bases
	# Last: 3 or 4 ACGT bases
	############################################
	
	opts.minContigLength = 0
	opts.comp_kmer       = 5
	opts.cross_cov_bins  = None
	opts.sam             = None
	opts.motifs_file     = None
	opts.skip_motifs     = None
	opts.control_run     = False
	opts.bas_whitelist   = None
	
	opts.control_pkl_name = os.path.abspath(opts.control_pkl_name)

	return opts,h5_files

def __initLog( opts ):
	"""Sets up logging based on command line arguments. Allows for three levels of logging:
	logging.error( ): always emitted
	logging.info( ) : emitted with --info or --debug
	logging.debug( ): only with --debug"""

	if os.path.exists(opts.logFile):
		os.remove(opts.logFile)

	logLevel = logging.DEBUG if opts.debug \
				else logging.INFO if opts.info \
				else logging.ERROR

	logger = logging.getLogger("")
	logger.setLevel(logLevel)
	
	# create file handler which logs even debug messages
	fh = logging.FileHandler(opts.logFile)
	fh.setLevel(logLevel)
	
	# create console handler with a higher log level
	ch = logging.StreamHandler()
	ch.setLevel(logLevel)
	
	# create formatter and add it to the handlers
	logFormat = "%(asctime)s [%(levelname)s] %(message)s"
	formatter = logging.Formatter(logFormat, "%Y-%m-%d %H:%M:%S")
	ch.setFormatter(formatter)
	fh.setFormatter(formatter)
	
	# add the handlers to logger
	logger.addHandler(ch)
	logger.addHandler(fh)

def __check_input( opts, args, parser ):
	"""
	Make sure the input is in the form of either a cmp.h5 file of aligned reads
	or a FOFN of unaligned bas.h5 files. Also make sure that a reference fasta 
	file is specified if 
	"""
	arg            = args[0]
	h5_files       = []
	opts.h5_labels = {}

	if arg[-6:]=="cmp.h5":
		print "Found cmp.h5 of aligned reads:"

		opts.h5_type                = "cmp"
		opts.cmph5_contig_lens      = {}
		opts.cmph5_contig_lens[arg] = {}

		h5_files.append(arg)
		print "  -- %s" % arg
		print "Getting contig information from %s..." % arg
		reader = CmpH5Reader(arg)
		for entry in reader.referenceInfoTable:
			name                                   = entry[3]
			length                                 = entry[4]
			slug_name                              = mbin.slugify(name)
			opts.cmph5_contig_lens[arg][slug_name] = length
			opts.h5_labels[arg]                    = "remove"
		reader.close()

	elif arg[-6:]=="bas.h5":
		print "Found bas.h5 of unaligned reads:"
		opts.h5_type        = "bas"
		h5_files.append(arg)
		opts.h5_labels[arg] = "remove"
		print "  -- %s" % arg

	elif arg[-5:]==".fofn":
		print "Found FOFN of bas.h5 files:"
		opts.h5_type = "bas"
		fns          = map(lambda x: x.strip("\n"), np.atleast_1d(open(arg, "r").read()))
		h5_files     = fns
		for fn in fns:
			print "  -- %s" % fn
			opts.h5_labels[fn] = "remove"

	if opts.h5_type=="bas":
		print "*************************************************************"
		print "* Motif filtering using unaligned reads is not recommended. *" 
		print "*         Aligned reads work much better for this!          *"
		print "*************************************************************"
		print ""

	if opts.h5_type=="bas" and opts.cross_cov_bins!=None:
		parser.error("Use of the --cross_cov_bins option is not compatible with bas.h5 inputs!")

	return h5_files

class FilterRunner:
	def __init__(self, opts, h5_files):
		self.opts     = opts
		self.h5_files = h5_files

	def add_degen_motifs( self, motifs, orig_control_means ):
		"""
		If a predetermined set of motifs is input using --motifs_file option,
		create a new entry for the degen motif in the control values dictionary
		by combining the existing data from the various specified versions 
		of the motif.
		"""
		keys_str          = "\n".join(orig_control_means.keys())
		new_control_means = orig_control_means
		n_degen           = 0
		for m in motifs:
			new_m = motif_tools.sub_bases(m)
			if new_m!=m:
				matches              = re.findall(new_m, keys_str)
				degen_mean           = np.mean([orig_control_means[match] for match in matches])
				new_control_means[m] = degen_mean
				logging.info("Adding degenerate motif %s to controls: %s" % (m, degen_mean))
				n_degen += 1

		return new_control_means, n_degen

	def bas_stream_files( self ):
		reads_ipds_fn        = os.path.join( self.opts.tmp, "read_ipds.tmp"      )
		reads_ipds_kmers_fn  = os.path.join( self.opts.tmp, "read_ipdskmers.tmp")
		reads_ipds_N_fn      = os.path.join( self.opts.tmp, "read_ipdsN.tmp"    )
		all_motifs           = defaultdict(list)
		logging.info("Unpickling the control IPDs...")
		control_means        = pickle.load(open(self.opts.control_pkl_name, "r"))
		logging.info("Done.")
		args = []
		for j,line in enumerate(open(reads_ipds_fn+".trans", "r").xreadlines()):
			args.append( (j, copy.copy(control_means), reads_ipds_fn, reads_ipds_N_fn, reads_ipds_kmers_fn, self.opts) )

		results          = mbin.launch_pool( self.opts.procs, get_motif_scores_from_read, args )
		highscore_motifs = [x for x in results if x is not None]

		# Keep only the shortest version of the high scoring motifs (reduces redundancy)
		keeper_motifs    = set()
		if len(highscore_motifs)>0:
			shortest_contiguous = min([len(m.split("-")[0]) for m in highscore_motifs])
			shortest_motifs     = [m for m in highscore_motifs if len(m.split("-")[0])==shortest_contiguous]
			to_del = set()
			for shorty in shortest_motifs:
				shorty_str =     shorty.split("-")[0]
				shorty_idx = int(shorty.split("-")[1])
				for motif in highscore_motifs:
					if motif!=shorty:
						motif_str =     motif.split("-")[0]
						motif_idx = int(motif.split("-")[1])
						match = re.search(shorty_str, motif_str)
						if match != None:
							if (shorty_idx + match.start()) == motif_idx:
								to_del.add( motif )
			for motif in highscore_motifs:
				if motif not in to_del:
					keeper_motifs.add(motif)
		return keeper_motifs

	def run(self, mbinRunner):
		####################################################
		# Filter out motifs without significant signatures
		####################################################
		logging.info("Getting top motifs from each contig...")

		if self.opts.h5_type=="cmp":
			
			self.contig_fasta_lens = {}
			for entry in SeqIO.parse(self.opts.contigs, "fasta"):
				name = entry.id
				if name.find("|quiver")>-1:
					# SMRT assemblies add |quiver to contig names, but this
					# gets dropped from the contig names in the cmp.h5 file.
					name = name.replace("|quiver","")
				self.contig_fasta_lens[mbin.slugify(name)] = len(entry.seq)

			contig_ipds_fns = glob.glob( os.path.join(self.opts.tmp, "*_ipds.tmp") )
			contigs         = map(lambda x: os.path.basename(x).split("_ipds.tmp")[0], contig_ipds_fns)
			ipds_fn_dict    = dict([(os.path.basename(ipds_fn).split("_ipds.tmp")[0],ipds_fn) for ipds_fn in contig_ipds_fns])

			contigs_for_transpose = []
			contigs_for_chunking  = []
			maxsize_for_transpose = 25000000 #25Mb
			for name in contigs:
				fsize = os.path.getsize(ipds_fn_dict[name])
				if fsize < maxsize_for_transpose:
					contigs_for_transpose.append(name)
				else:
					contigs_for_chunking.append(name)

			logging.info("Transposing %s case contigs..." % len(contigs_for_transpose))
			args    = [(contig,self.opts) for contig in contigs_for_transpose]
			results = mbin.launch_pool( self.opts.procs, transpose_contig_matrix, args )

			streamed_contig_dicts  = {}
			if len(contigs_for_transpose)>0:
				logging.info("Streaming through %s contigs..." % len(contigs_for_transpose))
				args    = [(self.opts, contig, i, len(contigs_for_transpose)) for i,contig in enumerate(contigs_for_transpose)]
				results = mbin.launch_pool( self.opts.procs, stream_case_control_files, args)
				
				streamed_contig_SCp    = map(lambda x: x[0], results)
				streamed_contig_SCp_N  = map(lambda x: x[1], results)
				streamed_contigs       = map(lambda x: x[2], results)
				
				for i,contig in enumerate(streamed_contigs):
					streamed_contig_dicts[contig] = {"SCp":streamed_contig_SCp[i], "SCp_N":streamed_contig_SCp_N[i]}

			chunked_contigs_dicts  = {}
			if len(contigs_for_chunking)>0:
				logging.info("Chunking %s contigs..." % len(contigs_for_chunking))
				
				for i,contig in enumerate(contigs_for_chunking):
					control_means,contig_SCp,contig_SCp_N,contig = chunk_case_control_files( self.opts, contig, i, len(contigs_for_chunking) )
					chunked_contigs_dicts[contig]                = {"SCp":contig_SCp, "SCp_N":contig_SCp_N}
				
			# Combine the contig dictionaries from both streaming and chunked paths
			def merge_two_dicts(x, y):
				"""Given two dicts, merge them into a new dict as a shallow copy."""
				z = x.copy()
				z.update(y)
				return z

			contig_dicts = merge_two_dicts( streamed_contig_dicts, chunked_contigs_dicts )

			keeper_control_ipds = {}
			keeper_motifs       = set()

			if self.opts.cross_cov_bins!=None:
				"""
				Using contig<-->bin mappings, collect methylation data from each
				contig and compile them into methylation scores for each bin.
				Then discover motifs based on bin-level scores.
				"""
				bin_map = {}
				for line in open(self.opts.cross_cov_bins, "rb").xreadlines():
					line            =     line.strip()
					contig          =     line.split(",")[0]
					bin_id          = int(line.split(",")[1])
					bin_map[contig] = bin_id

				bin_contig_dicts = {}
				for bin_id in bin_map.values():
					# Initialize the bin-level methylation dictionary
					bin_contig_dicts[bin_id] = {"SCp":{}, "SCp_N":{}}

				for contig,contig_d in contig_dicts.iteritems():
					# Make sure contig is binned
					if bin_map.get(contig):
						bin_id = bin_map[contig]
						bin_contig_dicts[bin_id] = build_bin_dict( contig, bin_id, contig_d, bin_contig_dicts[bin_id] )
					else:
						logging.info("Contig %s not found in cross-coverage binning results." % contig)
				
				bin_ids = list(set(bin_map.values()))
				bin_ids.sort()

				args = []
				for bin_id in bin_ids:
					# For each bin, do motif filtering and refinement
					if len(bin_contig_dicts[bin_id]["SCp"].keys())>0:
						
						bin_copy_contig_dicts = copy.deepcopy(bin_contig_dicts[bin_id])
						
						args.append( (bin_id,                         \
									  bin_copy_contig_dicts["SCp"],   \
									  bin_copy_contig_dicts["SCp_N"], \
									  self.opts,                      \
									  len(bin_ids),                   \
									  bin_id,                         \
									  "bin") )

				results = mbin.launch_pool( self.opts.procs, simplify_motifs, args )

				bin_keeper_motifs_list = map(lambda x: x[0], results)
				control_means_list     = map(lambda x: x[1], results)

				"""
				Add the control means for these bin motifs to the 
				complete set of control means for all detected motifs.
				"""
				for bin_keeper_motifs in bin_keeper_motifs_list:
					keeper_motifs = keeper_motifs | bin_keeper_motifs

				for sub_control_means in control_means_list:
					for motif,score in sub_control_means.iteritems():
						control_means[motif] = score

			else:
				args = []
				for j,(contig,contig_d) in enumerate(contig_dicts.iteritems()):
					# For each contig, do motif filtering and refinement
					copy_contig_dicts = copy.deepcopy(contig_d)

					args.append( (j,                          \
								  copy_contig_dicts["SCp"],   \
								  copy_contig_dicts["SCp_N"], \
								  self.opts,                  \
								  len(contig_dicts.keys()),   \
								  contig,                     \
								  "contig") )

				results = mbin.launch_pool( self.opts.procs, simplify_motifs, args )

				contig_keeper_motifs_list = map(lambda x: x[0], results)
				control_means_list        = map(lambda x: x[1], results)

				"""
				Add the control means for these bin motifs to the 
				complete set of control means for all detected motifs.
				"""
				for contig_keeper_motifs in contig_keeper_motifs_list:
					keeper_motifs = keeper_motifs | contig_keeper_motifs

				for control_means in control_means_list:
					for motif,score in control_means.iteritems():
						keeper_control_ipds[motif] = score

			# Rewrite the control so that it includes the new degenerate motifs.
			control_means,n_degen = self.add_degen_motifs( keeper_motifs, control_means)
			if n_degen>0:
				pickle.dump(control_means, open(self.opts.control_pkl_name, "wb"))
		
		elif self.opts.h5_type=="bas":
			logging.info("Transposing reads...")
			files   = [os.path.join(self.opts.tmp, "read_ipds.tmp"), os.path.join(self.opts.tmp, "read_ipdsN.tmp")]
			results = mbin.launch_pool( len(files), transpose_file, files )
			logging.info("Done.")
			logging.info("Streaming through reads for motif filtering...")
			keeper_motifs = self.bas_stream_files()
			
		logging.info("Keeping %s motifs for further analysis" % len(keeper_motifs))
		self.motifs = list(keeper_motifs)
		n_motifs    = len(keeper_motifs)
		f_motifs    = open(self.opts.motifs_fn, "w")
		for motif in keeper_motifs:
			f_motifs.write("%s\n" % motif)
		f_motifs.close()


if __name__ == "__main__":
	main()