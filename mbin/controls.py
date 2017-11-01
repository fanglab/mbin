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

def launch():
	opts,control_h5 = __parseArgs()
	__initLog(opts)

	extract_controls(opts, control_h5)

	print >> sys.stderr, "mBin control extraction has finished running. See log for details."

def extract_controls(opts, control_h5):
	"""

	"""
	controls       = ControlRunner(control_h5, opts)
	mbinRunner     = mbin.mbinRunner(opts)

	# Pulling the IPD data for each motif from the WGA cmp.h5 file
	motifs,bi_motifs = motif_tools.build_motif_dict(opts)
	opts.motifs      = motifs
	opts.bi_motifs   = bi_motifs

	logging.info("")
	logging.info("Preparing to create new control data in %s" % opts.control_tmp)
	controls.goto_control_output_dir()
	
	opts           = controls.scan_WGA_h5()
	filter_N_reads = opts.N_reads

	mbinRunner.launch_data_loader( control_h5, filter_N_reads, 1, opts )
	
	controls.analyze_WGA_reads()
	logging.info("Done.")
	logging.info("")

	# Building dictionary of mean control IPD values for each motif
	logging.info("Building dictionary of control values for all motifs...")
	logging.info("   * Initial build requires significant time and memory.")
	controls.combine_control_data_from_contigs()

	
	control_means    = controls.build_control_IPD_dict(motifs, bi_motifs)
	controls.return_to_orig_dir()

	logging.info("")
	logging.info("Cleaning up temp files from control data processing...")
	shutil.rmtree(opts.control_tmp)

	# Controls are loaded into control_means, now pickle them for easy
	# passing between parallel processes
	pickle.dump(control_means, open(opts.control_pkl_name, "wb"))

def chunks( l, n ):
	"""
	Yield successive n-sized chunks from l.
	"""
	for i in xrange(0, len(l), n):
		yield l[i:i+n]

def __parseArgs():
	"""Handle command line argument parsing"""

	usage = """%prog [--help] [options] wga_input_hdf5

	Example:

	buildcontrols -i --procs=4 --control_pkl_name=control_means.pkl wga_aligned_reads.cmp.h5

	"""

	parser = optparse.OptionParser( usage=usage, description=__doc__ )

	parser.add_option( "-d", "--debug", action="store_true", help="Increase verbosity of logging" )

	parser.add_option( "-i", "--info", action="store_true", help="Add basic logging" )

	parser.add_option( "--logFile", type="str", help="Write logging to file [log.controls]" )

	parser.add_option( "--subreadlength_min", type="int", help="Minimum subread length to include for analysis [100]" )

	parser.add_option( "--readlength_min", type="int", help="Minimum read length to include for analysis [100]" )

	parser.add_option( "--min_kmer", type="int", help="Minimum motif size to scan (contiguous motifs) [4]" )

	parser.add_option( "--max_kmer", type="int", help="Maximum motif size to scan (contiguous motifs) [6]" )

	parser.add_option( "--no_bipartite", action="store_true", help="Omit bipartite motifs [False]" )
	
	parser.add_option( "--bipart_first", type="str", help="Bipartite motif configuration: acceptable length of first determinate component (comma-separated string of integers) [3,4]" )
	
	parser.add_option( "--bipart_Ns", type="str", help="Bipartite motif configuration: acceptable length of middle indeterminate component (comma-separated string of integers) [5,6]" )
	
	parser.add_option( "--bipart_second", type="str", help="Bipartite motif configuration: acceptable length of second determinate component (comma-separated string of integers) [3,4]" )

	parser.add_option( "--mod_bases", type="str", help="String containing bases to query for mods. Changing this is not recommended ['A'] ['A']" )

	parser.add_option( "--minAcc", type="float", help="Min subread accuracy of read [0.8]" )

	parser.add_option( "--minMapQV", type="float", help="Min mapping QV of aligned read [240]" )

	parser.add_option( "--procs", type="int", help="Number of cores to use [4]" )

	parser.add_option( "--N_reads", type="int", help="Number of qualifying reads to include in analysis [1000000000]" )
	
	parser.add_option( "--min_motif_count", type="int", help="Number of motif sites required in WGA data to be included in controls dictionary [10]" )

	parser.add_option( "--control_pkl_name", type="str", help="Filename to save control IPD data from WGA sequencing [control_ipds.pkl]" )
	
	parser.set_defaults( logFile="log.buildcontrols",        \
						 info=False,                         \
						 debug=False,                        \
						 subreadlength_min=100,              \
						 readlength_min=100,                 \
						 min_kmer=4,                         \
						 max_kmer=6,                         \
						 no_bipartite=False,                 \
						 bipart_first="3,4",                 \
						 bipart_Ns="5,6",                    \
						 bipart_second="3,4",                \
						 mod_bases="A",                      \
						 minAcc=0.8,                         \
						 minMapQV=240,                       \
						 procs=4,                            \
						 N_reads=1000000000,                 \
						 min_motif_count=10,                 \
						 control_pkl_name="control_ipds.pkl")

	opts, args         = parser.parse_args( )
	control_h5         = __check_input( opts, args, parser )

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
	
	# As set, acceptible bipartite motifs would have
	# the following component lengths.
	# First: 3 or 4 ACGT bases
	# Ns: 5 or 6 unspecified N bases
	# Last: 3 or 4 ACGT bases
	############################################

	opts.control_tmp     = "ctrl_tmp"
	opts.tmp             = "tmp"
	opts.minContigLength = 0
	opts.comp_kmer       = 5
	opts.h5_type         = "cmp"
	opts.cross_cov_bins  = None
	opts.sam             = None
	opts.motifs_file     = None
	opts.skip_motifs     = None
	
	opts.control_pkl_name = os.path.abspath(opts.control_pkl_name)

	return opts,control_h5

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
	control_h5 = os.path.abspath(args[0])

	if control_h5[-6:]!="cmp.h5":
		parser.error("Please input a *.cmp.h5 file of aligned reads")

	# if opts.contigs==None:
	# 	parser.error("Please specify the fasta file used for the alignments in %s!" % control_h5)

	if len(args) != 1:
		parser.error( "Expected 1 argument." )

	return control_h5

def process_contig_chunk( args ):
	chunk_id        = args[0]
	cut_CMDs        = args[1]
	kmers           = args[2]
	cols_chunk      = args[3]
	n_chunks        = args[4]
	min_motif_count = args[5]
	logging.info("  - Control data: chunk %s/%s" % ((chunk_id+1), (n_chunks+1)))
	control_means = {}
	
	for cut_CMD in cut_CMDs:
		sts,stdOutErr = mbin.run_OS_command( cut_CMD )
	
	fns                = map(lambda x: x.split("> ")[-1], cut_CMDs)
	control_ipds_sub   = np.loadtxt(fns[0], dtype="float")
	control_ipds_N_sub = np.loadtxt(fns[1], dtype="int")
	# If there is only one row (read) for this contig, still treat as
	# a 2d matrix of many reads
	control_ipds_sub   = np.atleast_2d(control_ipds_sub)
	control_ipds_N_sub = np.atleast_2d(control_ipds_N_sub)
	
	not_found     = 0
	for j in range(len(cols_chunk)):
		motif = kmers[cols_chunk[j]]
		if np.sum(control_ipds_N_sub[:,j])>=min_motif_count:
			if np.sum(control_ipds_N_sub[:,j])>0:
				control_mean = np.dot(control_ipds_sub[:,j], control_ipds_N_sub[:,j]) / np.sum(control_ipds_N_sub[:,j])
			else:
				control_mean = 0
			control_means[motif] = control_mean
		else:
			not_found += 1

	return control_means,not_found

class ControlRunner:
	def __init__( self, wga_h5, opts ):
		"""
		Point to the appropriate WGA sequencing data files
		to generate the control IPD values for the motifs. 
		"""
		self.control_h5 = wga_h5
		self.opts       = opts
		self.orig_dir   = os.getcwd()

	def build_insilico_controls( self ):
		"""
		To be added...
		"""
		pass

	def goto_control_output_dir( self ):
		"""
		Create directory where data from control reads 
		will be gathered and analyzed.
		"""
		# Make control directory
		if os.path.exists(self.opts.control_tmp):
			shutil.rmtree(self.opts.control_tmp)
		os.mkdir(self.opts.control_tmp)
		os.chdir(self.opts.control_tmp)
		
		# Make tmp directory inside control directory
		if os.path.exists(self.opts.tmp):
			shutil.rmtree(self.opts.tmp)
		os.mkdir(self.opts.tmp)

	def return_to_orig_dir( self ):
		"""
		Back out of the control directory.
		"""
		os.chdir(self.orig_dir)

	def scan_WGA_h5( self ):
		"""
		Get some necessary information about the WGA cmp.h5 
		being used to generate the control IPD data.
		"""
		self.opts.h5_labels                          = {}
		self.opts.cmph5_contig_lens                  = {}
		self.opts.h5_labels[self.control_h5]         = "control"
		self.opts.cmph5_contig_lens[self.control_h5] = {}
		
		reader = CmpH5Reader(self.control_h5)
		for entry in reader.referenceInfoTable:
			name      = entry[3]
			length    = entry[4]
			slug_name = mbin.slugify(name)
			self.opts.cmph5_contig_lens[self.control_h5][slug_name] = length
		reader.close()

		return self.opts

	def analyze_WGA_reads( self ):
		"""
		Launch read scanning pipeline for building up
		control IPD values for motifs.
		"""
		control_fns = glob.glob( os.path.join(self.opts.tmp, "*.tmp"))
		ftypes      = set( map(lambda x: "_".join(os.path.basename(x).split("_")[2:]), control_fns) )
		for ftype in ftypes:
			if ftype in ["compkmers.tmp", "ipdskmers.tmp"]:
				first_fn = glob.glob( os.path.join(self.opts.tmp, "unitig_*_%s" % ftype) )[0]
				shutil.copy( first_fn, os.path.join(self.opts.tmp, "control_%s" % ftype) )
				for fn in glob.glob( os.path.join(self.opts.tmp, "unitig_*_%s" % ftype) ):
					os.remove(fn)
			else:
				to_cat = glob.glob( os.path.join(self.opts.tmp, "unitig_*_%s" % ftype) )
				to_cat.sort()
				outname = os.path.join(self.opts.tmp, "control_%s" % ftype )
				mbin.cat_list_of_files(to_cat, outname)
		
	def chunk_control_matrices( self, control_ipds_fn, control_ipds_N_fn, control_kmers_fn ):
		"""

		"""
		kmers       = np.atleast_1d(np.loadtxt(control_kmers_fn, dtype="str"))
		fns         = [control_ipds_fn, control_ipds_N_fn]
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
			args.append( (i, cut_CMDs, kmers, cols_chunk, n_chunks, self.opts.min_motif_count) )
		
		results = mbin.launch_pool(self.opts.procs, process_contig_chunk, args)
		
		logging.info("Combining motifs from all chunks of control data...")
		not_found     = 0
		control_means = {}
		for i,result in enumerate(results):
			not_found += result[1]
			for motif in result[0].keys():
				control_means[motif] = result[0][motif]
		logging.info("Done.")

		return control_means,not_found

	def combine_control_data_from_contigs( self ):
		"""
		If control WGA data contains multiple contigs, this 
		will combine them into one file for each data type
		so that the control IPD dictionary will be generated
		using data from all contigs.
		"""
		contigs_fns    = glob.glob( os.path.join(self.opts.tmp, "control_*.tmp") )
		labels_fns     = []
		strands_fns    = []
		lengths_fns    = []
		readnames_fns  = []
		ipds_fns       = []
		ipds_N_fns     = []
		comp_N_fns     = []
		comp_kmers_fns = []
		ipds_kmers_fns = []
		for fn in contigs_fns:
			if   fn.find("_labels.")>-1:
				labels_fns.append(fn)
			elif fn.find("_strand.")>-1:
				strands_fns.append(fn)
			elif fn.find("_lengths.")>-1:
				lengths_fns.append(fn)
			elif fn.find("_readnames.")>-1:
				readnames_fns.append(fn)
			elif fn.find("_ipds.")>-1:
				ipds_fns.append(fn)
			elif fn.find("_ipdsN.")>-1:
				ipds_N_fns.append(fn)
			elif fn.find("_compN.")>-1:
				comp_N_fns.append(fn)
			elif fn.find("_compkmers.")>-1:
				comp_kmers_fns.append(fn)
			elif fn.find("_ipdskmers.")>-1:
				ipds_kmers_fns.append(fn)

		labels_fns.sort()
		strands_fns.sort()
		lengths_fns.sort()
		readnames_fns.sort()
		ipds_fns.sort()
		ipds_N_fns.sort()
		comp_N_fns.sort()

		mbin.cat_list_of_files(labels_fns,    "control_labels.tmp")
		mbin.cat_list_of_files(strands_fns,   "control_strands.tmp")
		mbin.cat_list_of_files(lengths_fns,   "control_lengths.tmp")
		mbin.cat_list_of_files(readnames_fns, "control_names.tmp")
		mbin.cat_list_of_files(ipds_fns,      "control_ipds.tmp")
		mbin.cat_list_of_files(ipds_N_fns,    "control_ipdsN.tmp")
		mbin.cat_list_of_files(comp_N_fns,    "control_compN.tmp")
		
		shutil.copy(comp_kmers_fns[0], "control_compkmers.tmp")
		shutil.copy(ipds_kmers_fns[0], "control_ipdskmers.tmp")
		x = [os.remove(fn) for fn in comp_kmers_fns]
		x = [os.remove(fn) for fn in ipds_kmers_fns]

	def build_control_IPD_dict( self, motifs, bi_motifs ):
		"""

		"""
		control_ipds_fn   = glob.glob( "control_ipds.tmp" )
		control_ipds_N_fn = glob.glob( "control_ipdsN.tmp")
		control_kmers_fn  = glob.glob( "control_ipdskmers.tmp")

		if (len(control_ipds_fn)>1 or len(control_ipds_N_fn)>1 or len(control_kmers_fn)>1):
			raise Exception("*** Double check the control files. There should not be multiples for a file type.")

		control_means,not_found = self.chunk_control_matrices(control_ipds_fn[0], control_ipds_N_fn[0], control_kmers_fn[0])

		if not_found > 0:
			logging.info("")
			logging.warning("WARNING: could not find sufficient instances (>=%s) for %s motifs (out of %s total) in control data!" % (self.opts.min_motif_count, not_found, (len(motifs)+len(bi_motifs))))
			logging.warning("   * If this is alarming, try reducing --min_motif_count or increasing --N_reads, although you just might not have those motifs in your reference sequence.")
		
		logging.info("")
		logging.info("Writing control data to a pickled file: %s" % self.opts.control_pkl_name)
		pickle.dump( control_means, open( self.opts.control_pkl_name, "wb" ) )

		return control_means

if __name__ == "__main__":
	main()