import os,sys
from collections import Counter,defaultdict
from itertools import groupby
import math
import numpy as np
from pbcore.io.align.CmpH5IO import CmpH5Reader
from pbcore.io.BasH5IO import BasH5Reader
import shutil
import re
import copy
import random
from scipy import stats
import multiproc
import read_scanner
import cmph5_read
import baxh5_read
import optparse
import logging
import multiprocessing
import operator
import subprocess
import plotting
import glob
from itertools import izip
import warnings
import pickle
import pysam
import motif_tools
from control_ipds import ControlRunner

def run_OS_command( CMD ):
	p         = subprocess.Popen(CMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdOutErr = p.communicate()
	sts       = p.returncode
	if sts != 0:
		logging.warning("Failed command: %s" % CMD)
	return sts, stdOutErr

def launch_pool( procs, funct, args ):
	p    = multiprocessing.Pool(processes=procs)
	try:
		results = p.map(funct, args)
		p.close()
		p.join()
	except KeyboardInterrupt:
		p.terminate()
	return results

def process_contig_chunk( args ):
	chunk_id      = args[0]
	control_pkl   = args[1]
	cut_CMDs      = args[2]
	kmers         = args[3]
	cols_chunk    = args[4]
	min_motif_N   = args[5]
	contig_id     = args[6]
	n_chunks      = args[7]
	n_contigs     = args[8]
	opts          = args[9]
	logging.info("  - Contig %s/%s: chunk %s/%s" % ((contig_id+1), n_contigs, (chunk_id+1), (n_chunks+1)))
	control_means = pickle.load(open(control_pkl, "rb"))
	contig_motifs = {}
	case_motif_Ns = {}
	
	for cut_CMD in cut_CMDs:
		sts,stdOutErr = run_OS_command( cut_CMD )
	
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

def cat_contig_files_from_subprocs( tup ):
	contig    = tup[0]
	tmp       = tup[1]
	j         = tup[2]
	N_contigs = tup[3]
	logging.info("...contig %s (%s/%s)" % (contig, j, N_contigs))
	contig_fns = glob.glob( os.path.join(tmp, "chunk_*", "%s_*.tmp" % contig) )
	ftypes = set( map(lambda x: os.path.basename(x).split("_")[-1], contig_fns) )
	for ftype in ftypes:
		if ftype in ["compkmers.tmp", "ipdskmers.tmp"]:
			first_fn = glob.glob( os.path.join(tmp, "chunk_*", "%s_%s" % (contig, ftype)) )[0]
			shutil.copy( first_fn, os.path.join(tmp, "%s_%s" % (contig,ftype)))
		else:
			to_cat  = glob.glob( os.path.join(tmp, "chunk_*", "%s_%s" % (contig, ftype)) )
			to_cat.sort()
			outname = os.path.join(tmp, "%s_%s" % (contig,ftype) )
			cmph5_read.cat_list_of_files(to_cat, outname, del_ins=False)

def cat_subreads_files_from_subprocs( subproc_tmp_files, movie_i, opts ):
	for ftype in subproc_tmp_files.keys():
		if ftype in ["compkmers.tmp", "ipdskmers.tmp"]:
			first_fn = glob.glob( os.path.join(opts.tmp, "chunk_*", "subreads_%s" % ftype) )[0]
			shutil.copy( first_fn, os.path.join(opts.tmp, "subreads_%s" % ftype.replace(".tmp", ".%s.tmp" % movie_i)))
		else:
			ftype_fns = list(subproc_tmp_files[ftype])
			sorted(ftype_fns)
			catted_name = os.path.join(opts.tmp, "subreads_%s" % ftype.replace(".tmp", ".%s.tmp" % movie_i))
			baxh5_read.cat_list_of_files( ftype_fns, catted_name, del_ins=False )

def transpose_file( fn ):
	def run_OS( CMD ):
		p         = subprocess.Popen(CMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdOutErr = p.communicate()
		sts       = p.returncode
		if sts != 0:
			logging.warning("Failed command: %s" % CMD)
		return sts, stdOutErr

	trans_script   = "/hpc/users/beaulj01/gitRepo/eBinning/src/transpose.awk"
	trans_CMD      = "awk -f %s %s > %s.trans" % (trans_script, fn, fn)
	logging.info(trans_CMD)
	sts, stdOutErr = run_OS(trans_CMD)

def transpose_contig_matrix( args ):
	contig  = args[0]
	control = args[1]
	logging.info("  Transposing %s" % contig)
	if not control:
		contig_ipds_fn       = os.path.join( "tmp", "%s_ipds.tmp"       % contig)
		contig_ipds_kmers_fn = os.path.join( "tmp", "%s_ipdskmers.tmp"  % contig)
		contig_ipds_N_fn     = os.path.join( "tmp", "%s_ipdsN.tmp"      % contig)
	else:
		contig_ipds_fn       = os.path.join( "control", "tmp", "%s_ipds.tmp"       % contig)
		contig_ipds_kmers_fn = os.path.join( "control", "tmp", "%s_ipdskmers.tmp"  % contig)
		contig_ipds_N_fn     = os.path.join( "control", "tmp", "%s_ipdsN.tmp"      % contig)
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
	minMotifIPD      = tup[3]
	min_motif_N      = tup[4]
	units_N          = tup[5]
	unit_name        = tup[6]
	control_means_fn = tup[7]
	unit_type        = tup[8]

	# Keep only the motifs with scores higher than the minMotifIPD value
	highscore_motifs = dict( [(motif,SCp) for motif,SCp in unit_SCp.items() if SCp>=minMotifIPD and unit_SCp_N[motif]>=min_motif_N] )

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
	control_means = pickle.load(open(control_means_fn, "rb"))
	# Add control values for new degenerate bases (mean of the specified motifs)
	for d_motif,o_motifs in degen_members.iteritems():
		degen_mean = np.mean([control_means[o_m] for o_m in o_motifs])
		control_means[d_motif] = degen_mean
	
	# Delete all control mean entries for motifs not in the keeper set
	for motif in control_means.keys():
		if motif not in keeper_motifs:
			del control_means[motif]

	return keeper_motifs, control_means

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
	control_means_fn     = tup[0]
	contig               = tup[1]
	j                    = tup[2]
	contigs_N            = tup[3]
	minMotifIPD          = tup[4]
	min_motif_N          = tup[5]
	
	logging.info("   ...streaming contig %s (%s/%s)..." % (contig,(j+1),contigs_N) )

	control_means        = pickle.load(open(control_means_fn, "rb"))
	contig_SCp           = {}
	contig_SCp_N         = {}
	keeper_motifs        = set()
	contig_ipds_fn       = os.path.join( "tmp", "%s_ipds.tmp"      % contig)
	contig_ipds_kmers_fn = os.path.join( "tmp", "%s_ipdskmers.tmp" % contig)
	contig_ipds_N_fn     = os.path.join( "tmp", "%s_ipdsN.tmp"     % contig)
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

def chunk_case_control_files( control_means_fn, contig, j, contigs_N, minMotifIPD, min_motif_N, procs, opts ):
	def chunks( l, n ):
		"""
		Yield successive n-sized chunks from l.
		"""
		for i in xrange(0, len(l), n):
			yield l[i:i+n]

	logging.info("   ...chunking contig %s (%s/%s)..." % (contig,(j+1),contigs_N))

	contig_SCp           = {}
	contig_SCp_N         = {}
	keeper_motifs        = set()
	control_means        = pickle.load(open(control_means_fn, "rb"))
	contig_ipds_fn       = os.path.join( "tmp", "%s_ipds.tmp"      % contig)
	contig_ipds_N_fn     = os.path.join( "tmp", "%s_ipdsN.tmp"     % contig)
	contig_ipds_kmers_fn = os.path.join( "tmp", "%s_ipdskmers.tmp" % contig)
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
		args.append( (i, control_means_fn, cut_CMDs, kmers, cols_chunk, min_motif_N, j, n_chunks, contigs_N, opts) )
	results = launch_pool(procs, process_contig_chunk, args)
	for i,result in enumerate(results):
		for motif in result[0].keys(): # contig_SCp,contig_SCp_N
			contig_SCp[motif]   = result[0][motif]
			contig_SCp_N[motif] = result[1][motif]

	return control_means, contig_SCp, contig_SCp_N, contig

def get_motif_scores_from_read( args ):
	j                   = args[0]
	control_means       = args[1]
	reads_ipds_fn       = args[2]
	reads_ipds_N_fn     = args[3]
	reads_ipds_kmers_fn = args[4]
	minMotifIPD         = args[5]
	min_motif_N         = args[6]
	min_motif_reads     = args[7]
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
				if case_read_N>=min_motif_N:
					case_mean            = float(l1.split()[read_i])
					score                = case_mean - control_means[motif]
					if score > minMotifIPD:
						n_reads_passing += 1
	if n_reads_passing >= min_motif_reads:
		return motif

def combine_subreads_for_read_level( tup ):
	cmph5_file   = tup[0]
	contig       = tup[1]
	tmp          = tup[2]
	h5_labels    = tup[3]
	j            = tup[4]
	N_contigs    = tup[5]
	logging.info("...contig %s (%s/%s)" % (contig, j, N_contigs))
	contig_fns = glob.glob( os.path.join(tmp, "%s_*.tmp" % contig) )
	for fn in contig_fns:
		if   fn.find("_labels.tmp")>-1:
			labels_fn     = fn
		elif fn.find("_lengths.tmp")>-1:
			lengths_fn    = fn
		elif fn.find("_readnames.tmp")>-1:
			readnames_fn  = fn
		elif fn.find("_ipds.tmp")>-1:
			ipds_fn       = fn
		elif fn.find("_ipdsN.tmp")>-1:
			ipds_N_fn     = fn
		elif fn.find("_compN.tmp")>-1:
			comp_N_fn     = fn
		elif fn.find("_compkmers.tmp")>-1:
			comp_kmers_fn = fn
		elif fn.find("_ipdskmers.tmp")>-1:
			ipds_kmers_fn = fn
		elif fn.find("_strand.tmp")>-1:
			strand_fn     = fn

	readname_row_idx = defaultdict(list)
	for i,line in enumerate(open(readnames_fn, "r").xreadlines()):
		readname = line.strip()
		readname_row_idx[readname].append(i)
	n_subreads = i+1

	for line in open(ipds_kmers_fn).xreadlines():
		n_motifs = len(line.strip().split("\t"))

	sublengths = np.loadtxt(lengths_fn,    dtype="int")
	ipds       = np.loadtxt(ipds_fn,       dtype="float")
	ipds_N     = np.loadtxt(ipds_N_fn,     dtype="int")
	ipds_kmers = np.loadtxt(ipds_kmers_fn, dtype="string")
	comp_N     = np.loadtxt(comp_N_fn,     dtype="int")
	comp_kmers = np.loadtxt(comp_kmers_fn, dtype="string")
	strands    = np.loadtxt(strand_fn,     dtype="int")

	if n_subreads>1:
		if n_motifs==1:
			# Only one motif survived motif-filtering
			# Case #3
			ipds       = ipds.reshape(ipds.shape[0],1)
			ipds_N     = ipds_N.reshape(ipds_N.shape[0],1)
			ipds_kmers = ipds_kmers.reshape(1)
		elif n_motifs>1:
			# The ipd data is loaded properly in matrix form
			# Case #1
			pass
		elif n_motifs==0:
			pass
	elif n_subreads==1:
		sublengths = sublengths.reshape(1)
		strands    = strands.reshape(1)
		comp_N     = comp_N.reshape(1,comp_N.shape[0])
		if n_motifs==1:
			# Case #4
			ipds       = ipds.reshape(1,1)
			ipds_N     = ipds_N.reshape(1,1)
			ipds_kmers = ipds_kmers.reshape(1)
		elif n_motifs>1:
			# Case #2
			ipds       = ipds.reshape(1,ipds.shape[0])
			ipds_N     = ipds_N.reshape(1,ipds_N.shape[0])
		elif n_motifs==0:
			pass
		
	reads_names_fn      = os.path.join(tmp, "%s_read_names.tmp"      % contig)
	reads_labels_fn     = os.path.join(tmp, "%s_read_labels.tmp"     % contig)
	reads_lengths_fn    = os.path.join(tmp, "%s_read_lengths.tmp"    % contig)
	reads_contig_fn     = os.path.join(tmp, "%s_read_contig.tmp"     % contig)
	reads_ipds_fn       = os.path.join(tmp, "%s_read_ipds.tmp"       % contig)
	reads_ipds_N_fn     = os.path.join(tmp, "%s_read_ipdsN.tmp"     % contig)
	reads_ipds_kmers_fn = os.path.join(tmp, "%s_read_ipdskmers.tmp" % contig)
	reads_comp_N_fn     = os.path.join(tmp, "%s_read_compN.tmp"     % contig)
	reads_comp_kmers_fn = os.path.join(tmp, "%s_read_compkmers.tmp" % contig)
	reads_strands_fn    = os.path.join(tmp, "%s_read_strands.tmp"    % contig)

	f_reads      = open(reads_names_fn,      "w")
	f_labels     = open(reads_labels_fn,     "w")
	f_lengths    = open(reads_lengths_fn,    "w")
	f_contig     = open(reads_contig_fn,     "w")
	f_ipds       = open(reads_ipds_fn,       "w")
	f_ipds_N     = open(reads_ipds_N_fn,     "w")
	f_comp_N     = open(reads_comp_N_fn,     "w")
	f_strands    = open(reads_strands_fn, "w")
	label        = h5_labels[cmph5_file]
	for readname,row_idx in readname_row_idx.iteritems():
		f_reads.write(  "%s\n" % readname)
		f_labels.write( "%s\n" % label)
		f_contig.write( "%s\n" % contig)
		comp_N_list  = comp_N[row_idx,:].sum(axis=0)
		strands_list = strands[row_idx]
		# ipds_list   = ipds[row_idx,:].mean(axis=0)
		ipds_N_list  = ipds_N[row_idx,:].sum(axis=0)
		ipds_list    = []
		for k in range(n_motifs):
			read_motifs_N = np.sum(ipds_N[row_idx,k])
			read_ipds_sum = np.sum(ipds[row_idx,k] * ipds_N[row_idx,k])
			if read_motifs_N > 0:
				motif_mean = read_ipds_sum / read_motifs_N
			else:
				motif_mean = 0.0
			ipds_list.append(motif_mean)

		# Normalize composition kmer counts
		normed_comp_N_list = map(lambda x: math.log( float(x)/sum(comp_N_list) ), comp_N_list)
		readlength         = sublengths[row_idx].sum()
		f_lengths.write( "%s\n" % readlength)
		f_ipds.write(    "%s\n" % "\t".join(map(lambda x: str(round(x,3)), ipds_list)))
		f_ipds_N.write(  "%s\n" % "\t".join(map(lambda x: str(x),          ipds_N_list)))
		f_comp_N.write(  "%s\n" % "\t".join(map(lambda x: str(round(x,6)), normed_comp_N_list)))
		f_strands.write( "%s\n" %  ",".join(map(lambda x: str(x),          strands_list)))
	
	shutil.copy(ipds_kmers_fn, os.path.join(tmp, "%s_read_ipdskmers.tmp" % contig))
	shutil.copy(comp_kmers_fn, os.path.join(tmp, "%s_read_compkmers.tmp" % contig))
	f_reads.close()
	f_labels.close()
	f_lengths.close()
	f_contig.close()
	f_ipds.close()
	f_ipds_N.close()
	f_comp_N.close()
	f_strands.close()

	# Remove the subread-level barcodes for each contig
	for fn in contig_fns:
		os.remove(fn)

class eBinner:
	def __init__( self ):
		"""
		Parse the options and arguments, then instantiate the logger. 
		"""
		self.__parseArgs( )
		self.__initLog( )

	def chunks( self, l, n ):
		"""
		Yield successive n-sized chunks from l.
		"""
		for i in xrange(0, len(l), n):
			yield l[i:i+n]

	def run_tSNE( self, cmd ):
		sts        = 1
		perplexity = 30
		while sts!=0 and perplexity>1:
			cmd            = cmd.split(" -p ")[0]
			cmd           += " -p %s" % perplexity
			logging.info("Running %s" % cmd)
			sts, stdOutErr = run_OS_command( cmd )
			perplexity -= 1
			if sts!=0:
				logging.info("   ...trying again (p = %s)..." % perplexity)

	def launch_data_loader( self, h5_file, N_reads, movie_i ):
		logging.info("Loading data from %s..." % h5_file)
		results = self.launch_subprocs( h5_file, N_reads )
		if self.opts.h5_type=="cmp":
			contig_tmps = defaultdict(list)
			for p,result in enumerate(results):
				for fn in result:
					contig = "_".join(fn.split("/")[1].split("_")[:-1])
					contig_tmps[contig].append(fn)
			logging.info("Catting contig-specific data from parallel processes...")
			args    = [(contig, self.opts.tmp, i, len(contig_tmps.keys())) for i,contig in enumerate(contig_tmps.keys())]
			results = launch_pool( self.opts.procs, cat_contig_files_from_subprocs, args )
		elif self.opts.h5_type=="bas":
			logging.info("Catting subprocess-specific data from parallel processes...")
			subproc_tmp_files = defaultdict(list)
			for result in results:
				for fn in result:
					ftype = fn.split(os.sep)[-1].split("subreads_")[-1]
					subproc_tmp_files[ftype].append(os.path.join(self.opts.tmp, fn))
			
			results = cat_subreads_files_from_subprocs( subproc_tmp_files, movie_i, self.opts )			
			logging.info("Done.")	

		for chunkdir in glob.glob( os.path.join(self.opts.tmp, "chunk_*")):
			shutil.rmtree(chunkdir)

	def fasta_iter( self, fasta_name ):
		"""
		Given a fasta file, yield tuples of (header, sequence).
		"""
		fh = open(fasta_name)
		# ditch the boolean (x[0]) and just keep the header or sequence since
		# we know they alternate.
		faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
		for header in faiter:
			# drop the ">"
			header = header.next()[1:].strip()
			# join all sequence lines to one.
			seq = "".join(s.strip() for s in faiter.next())
			yield header, seq

	def combine_read_level_barcodes_across_contigs( self ):
		contigs_fns = glob.glob( os.path.join(self.opts.tmp, "*_read_*.tmp") )
		labels_fns     = []
		lengths_fns    = []
		refs_fns       = []
		readnames_fns  = []
		strands_fns    = []
		ipds_fns       = []
		ipds_N_fns     = []
		comp_N_fns     = []
		comp_kmers_fns = []
		ipds_kmers_fns = []
		for fn in contigs_fns:
			if   fn.find("_labels.tmp")>-1:
				labels_fns.append(fn)
			elif fn.find("_lengths.tmp")>-1:
				lengths_fns.append(fn)
			elif fn.find("_contig.tmp")>-1:
				refs_fns.append(fn)
			elif fn.find("_names.tmp")>-1:
				readnames_fns.append(fn)
			elif fn.find("_strands.tmp")>-1:
				strands_fns.append(fn)
			elif fn.find("_ipds.tmp")>-1:
				ipds_fns.append(fn)
			elif fn.find("_ipdsN.tmp")>-1:
				ipds_N_fns.append(fn)
			elif fn.find("_compN.tmp")>-1:
				comp_N_fns.append(fn)
			elif fn.find("_compkmers.tmp")>-1:
				comp_kmers_fns.append(fn)
			elif fn.find("_ipdskmers.tmp")>-1:
				ipds_kmers_fns.append(fn)

		labels_fns.sort()
		lengths_fns.sort()
		refs_fns.sort()
		readnames_fns.sort()
		strands_fns.sort()
		ipds_fns.sort()
		ipds_N_fns.sort()
		comp_N_fns.sort()
		

		cmph5_read.cat_list_of_files(labels_fns,    self.read_labels_fn, del_ins=False)
		cmph5_read.cat_list_of_files(lengths_fns,   self.read_lengths_fn, del_ins=False)
		cmph5_read.cat_list_of_files(refs_fns,      self.read_refs_fn, del_ins=False)
		cmph5_read.cat_list_of_files(readnames_fns, self.read_names_fn, del_ins=False)
		cmph5_read.cat_list_of_files(strands_fns,   self.read_strands_fn, del_ins=False)
		cmph5_read.cat_list_of_files(ipds_fns,      self.read_SMp_fn, del_ins=False)
		cmph5_read.cat_list_of_files(ipds_N_fns,    self.read_SMp_counts_fn, del_ins=False)
		cmph5_read.cat_list_of_files(comp_N_fns,    self.read_comp_fn, del_ins=False)
		shutil.copy(comp_kmers_fns[0],              self.read_comp_kmers_fn)
		shutil.copy(ipds_kmers_fns[0],              self.read_SMp_kmers_fn)

		for fn in comp_kmers_fns:
			os.remove(fn)

		# for fn in ipds_kmers_fns:
		# 	os.remove(fn)

	def combine_subreads_for_contig_level( self, cmph5_file ):
		readnames         = np.loadtxt(self.read_names_fn,       dtype="str")
		read_refs         = np.loadtxt(self.read_refs_fn,        dtype="str")
		reads_SMp         = np.loadtxt(self.read_SMp_fn,         dtype="float")
		reads_SMp_counts  = np.loadtxt(self.read_SMp_counts_fn,  dtype="int")
		reads_lengths     = np.loadtxt(self.read_lengths_fn,     dtype="int")
		
		if len(reads_SMp.shape)==1:
			reads_SMp        = reads_SMp.reshape(reads_SMp.shape[0],1)
			reads_SMp_counts = reads_SMp_counts.reshape(reads_SMp_counts.shape[0],1)
		
		n_motifs  = reads_SMp.shape[1]
		ref_SCp   = {}
		ref_SCp_N = {}
		ref_bases = {}
		for ref in set(read_refs):
			idx            = read_refs==ref
			ref_bases[ref] = reads_lengths[idx].sum()
			ref_SCp[ref]   = []
			ref_SCp_N[ref] = []
			for j in range(n_motifs):
				motif_counts  = reads_SMp_counts[idx,j]
				motif_vals    = reads_SMp[idx,j]
				# Don't include data from reads without a motif instance
				motif_vals    = motif_vals[motif_counts>0]
				motif_counts  = motif_counts[motif_counts>0]
				if sum(motif_counts) > 0:
					contig_SCp   = sum(motif_vals * motif_counts) / sum(motif_counts)
					contig_SCp_N = sum(motif_counts)
				else:
					contig_SCp   = 0.0
					contig_SCp_N = 0
				ref_SCp[ref].append( contig_SCp )
				ref_SCp_N[ref].append( contig_SCp_N )

		def calc_contig_comps( contigs_fasta ):
			contigs_comp_strings = {}
			contig_fasta_lens    = {}
			for i,(name,contig_seq) in enumerate(self.fasta_iter(contigs_fasta)):
				contig_comp  = []
				contig_kmers = read_scanner.kmer_freq( "cmp", contig_seq, 0, self.opts )
				for kmer,count in contig_kmers.iteritems():
					kmer_normed_comp = math.log(float(count) / sum(contig_kmers.values()))
					contig_comp.append(kmer_normed_comp)
				contig_comp_str               = "\t".join(map(lambda x: str(round(x,6)), contig_comp))
				if name.find("|quiver")>-1:
					name = name.replace("|quiver", "") 
				refName                       = cmph5_read.slugify(name)
				contigs_comp_strings[refName] = contig_comp_str
				contig_fasta_lens[refName]    = len(contig_seq)
			contigs_comp_strings["unknown"]   = "\t".join(map(lambda x: str(round(x,3)), [0.000] * len(contig_comp)))
			return contigs_comp_strings,contig_fasta_lens

		# Calculate composition profiles for all contigs in the contigs.fasta file
		contigs_comp_strings,self.contig_fasta_lens = calc_contig_comps( self.opts.contigs )

		f_ref_names       = open(self.contig_names_fn,   "w")
		f_ref_SCp         = open(self.contig_SCp_fn,     "w")
		f_ref_SCp_N       = open(self.contig_SCp_N_fn,   "w")
		f_ref_comp        = open(self.contig_comp_fn,    "w")
		f_ref_labels      = open(self.contig_labels_fn,  "w")
		f_ref_covs        = open(self.contig_cov_fn,     "w")
		f_ref_lens        = open(self.contig_lengths_fn, "w")
		if self.opts.comp_only:
			ref_names         = contigs_comp_strings.keys()
			ref_names_str     = "\n".join(ref_names)
			f_ref_names.write("%s\n" % ref_names_str)
			for ref in contigs_comp_strings.keys():
				ref_SCp_str   = "0.000\t0.000\t0.000"
				ref_SCp_N_str = "0\t0\t0"
				f_ref_SCp.write(   "%s\n" % ref_SCp_str)
				f_ref_SCp_N.write( "%s\n" % ref_SCp_N_str)
				
				# Sometimes the contig name in cmp.h5 trims off the "|quiver" in "<unitig_N>|quiver"
				try:
					f_ref_comp.write(  "%s\n" % contigs_comp_strings[ref])
				except KeyError:
					# quiver_ref = cmph5_read.slugify(ref+"|quiver")
					quiver_ref = cmph5_read.slugify(ref)
					f_ref_comp.write(  "%s\n" % contigs_comp_strings[quiver_ref])
				
				f_ref_labels.write("%s\n" % self.opts.h5_labels[cmph5_file])
				if self.opts.sam!=None:
					if ref=="unknown":
						coverage = 0
					else:
						coverage = float(ref_bases[ref]) / self.contig_fasta_lens[ref]
				else:
					coverage = float(ref_bases[ref]) / self.opts.cmph5_contig_lens[cmph5_file][ref]
				f_ref_covs.write(  "%.2f\n" % coverage)
		else:
			ref_names         = ref_SCp.keys()
			ref_names_str     = "\n".join(ref_names)
			f_ref_names.write("%s\n" % ref_names_str)
			for ref in ref_names:
				ref_SCp_str   = "\t".join(map(lambda x: str(round(x,3)), ref_SCp[ref]))
				ref_SCp_N_str = "\t".join(map(lambda x: str(x), ref_SCp_N[ref]))
				f_ref_SCp.write(   "%s\n" % ref_SCp_str)
				f_ref_SCp_N.write( "%s\n" % ref_SCp_N_str)
				
				# Sometimes the contig name in cmp.h5 trims off the "|quiver" in "<unitig_N>|quiver"
				try:
					f_ref_comp.write(  "%s\n" % contigs_comp_strings[ref])
				except KeyError:
					# quiver_ref = cmph5_read.slugify(ref+"|quiver")
					quiver_ref = cmph5_read.slugify(ref)
					f_ref_comp.write(  "%s\n" % contigs_comp_strings[quiver_ref])
				
				f_ref_labels.write("%s\n" % self.opts.h5_labels[cmph5_file])
				if self.opts.sam!=None:
					if ref=="unknown":
						coverage = 0
					else:
						coverage = float(ref_bases[ref]) / self.contig_fasta_lens[ref]
				else:
					coverage = float(ref_bases[ref]) / self.opts.cmph5_contig_lens[cmph5_file][ref]
				f_ref_covs.write(  "%.2f\n" % coverage)
				f_ref_lens.write(  "%s\n" % self.opts.cmph5_contig_lens[cmph5_file][ref])
		f_ref_names.close()
		f_ref_SCp.close()
		f_ref_SCp_N.close()
		f_ref_comp.close()
		f_ref_labels.close()
		f_ref_covs.close()
		f_ref_lens.close()

	def combine_contigs_for_bin_level( self ):
		"""
		If --cross_cov_bins is specified, aggregate IPD scores from contigs
		into bin-level scores based on their assignments from cross-coverage
		binning results. Recall that not all contigs will always have bin 
		assignments.
		"""

		bin_map = dict(np.loadtxt(self.opts.cross_cov_bins, delimiter=",", dtype="str"))

		contigs_names   = np.loadtxt(self.contig_names_fn,   dtype="str")
		contigs_SCp     = np.loadtxt(self.contig_SCp_fn,     dtype="float")
		contigs_SCp_N   = np.loadtxt(self.contig_SCp_N_fn,   dtype="int")
		contigs_lengths = np.loadtxt(self.contig_lengths_fn, dtype="int")
		
		if len(contigs_SCp.shape)==1:
			contigs_SCp   = contigs_SCp.reshape(contigs_SCp.shape[0],1)
			contigs_SCp_N = contigs_SCp_N.reshape(contigs_SCp_N.shape[0],1)

		# Ignore contigs without bin assignments
		binned_contigs_names   = []
		binned_contigs_SCp     = []
		binned_contigs_SCp_N   = []
		binned_contigs_lengths = []
		binned_contigs_bins    = []
		for i,name in enumerate(contigs_names):
			if bin_map.get(name):
				binned_contigs_names.append(name)
				binned_contigs_SCp.append(contigs_SCp[i,:])
				binned_contigs_SCp_N.append(contigs_SCp_N[i,:])
				binned_contigs_lengths.append(contigs_lengths[i])
				binned_contigs_bins.append(bin_map[name])
			else:
				pass

		binned_contigs_names   = np.array(binned_contigs_names)
		binned_contigs_SCp     = np.array(binned_contigs_SCp)
		binned_contigs_SCp_N   = np.array(binned_contigs_SCp_N)
		binned_contigs_lengths = np.array(binned_contigs_lengths)
		binned_contigs_bins    = np.array(binned_contigs_bins)

		n_motifs  = binned_contigs_SCp.shape[1]
		bins_SBp   = {}
		bins_SBp_N = {}
		bins_sizes = {}
		for bin_id in set(bin_map.values()):
			idx                = binned_contigs_bins==bin_id
			bins_sizes[bin_id] = binned_contigs_lengths[idx].sum()
			bins_SBp[bin_id]   = []
			bins_SBp_N[bin_id] = []
			for j in range(n_motifs):
				motif_counts  = binned_contigs_SCp_N[idx,j]
				motif_vals    = binned_contigs_SCp[idx,j]
				# Don't include data from contigs without a motif instance
				motif_vals    = motif_vals[motif_counts>0]
				motif_counts  = motif_counts[motif_counts>0]
				if sum(motif_counts) > 0:
					motif_SBp   = sum(motif_vals * motif_counts) / sum(motif_counts)
					motif_SBp_N = sum(motif_counts)
				else:
					motif_SBp   = 0.0
					motif_SBp_N = 0
				bins_SBp[bin_id].append( motif_SBp )
				bins_SBp_N[bin_id].append( motif_SBp_N )

		f_bin_names = open(self.bin_names_fn, "wb")
		f_bin_SBp   = open(self.bin_SCp_fn,   "wb")
		f_bin_SBp_N = open(self.bin_SCp_N_fn, "wb")
		f_bin_sizes = open(self.bin_sizes_fn, "wb")
		
		bin_names = bins_SBp.keys()
		bin_names.sort()

		for bin_id in bin_names:
			if bins_SBp.get(bin_id):
				bins_SBp_str   = "\t".join(map(lambda x: str(round(x,3)), bins_SBp[bin_id]))
				bins_SBp_N_str = "\t".join(map(lambda x: str(x), bins_SBp_N[bin_id]))
				f_bin_SBp.write(   "%s\n" % bins_SBp_str)
				f_bin_SBp_N.write( "%s\n" % bins_SBp_N_str)
				f_bin_sizes.write( "%s\n" % bins_sizes[bin_id])
				f_bin_names.write( "%s\n" % bin_id)
			else:
				# Did not collect data for that bin
				pass
		
		f_bin_names.close()
		f_bin_SBp.close()
		f_bin_SBp_N.close()
		f_bin_sizes.close()

	def combine_subread_data_across_bas_movies( self ):
		subreads_fns   = glob.glob( os.path.join(self.opts.tmp, "subreads_*.*.tmp") )
		labels_fns     = []
		lengths_fns    = []
		readnames_fns  = []
		ipds_fns       = []
		ipds_N_fns     = []
		comp_N_fns     = []
		comp_kmers_fns = []
		ipds_kmers_fns = []
		for fn in subreads_fns:
			if   fn.find("_labels.")>-1:
				labels_fns.append(fn)
			elif fn.find("_lengths.")>-1:
				lengths_fns.append(fn)
			elif fn.find("_names.")>-1:
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
		lengths_fns.sort()
		readnames_fns.sort()
		ipds_fns.sort()
		ipds_N_fns.sort()
		comp_N_fns.sort()
		
		cmph5_read.cat_list_of_files(labels_fns,    os.path.join(self.opts.tmp, "subreads_labels.tmp"))
		cmph5_read.cat_list_of_files(lengths_fns,   os.path.join(self.opts.tmp, "subreads_lengths.tmp"))
		cmph5_read.cat_list_of_files(readnames_fns, os.path.join(self.opts.tmp, "subreads_names.tmp"))
		cmph5_read.cat_list_of_files(ipds_fns,      os.path.join(self.opts.tmp, "subreads_ipds.tmp"))
		cmph5_read.cat_list_of_files(ipds_N_fns,    os.path.join(self.opts.tmp, "subreads_ipdsN.tmp"))
		cmph5_read.cat_list_of_files(comp_N_fns,    os.path.join(self.opts.tmp, "subreads_compN.tmp"))
		shutil.copy(comp_kmers_fns[0],              os.path.join(self.opts.tmp, "subreads_compkmers.tmp"))
		shutil.copy(ipds_kmers_fns[0],              os.path.join(self.opts.tmp, "subreads_ipdskmers.tmp"))

		for fn in comp_kmers_fns:
			os.remove(fn)
		for fn in ipds_kmers_fns:
			os.remove(fn)

	def bas_combine_subreads_for_read_level( self, tmp_run=False ):
		subreads_fns = glob.glob( os.path.join(self.opts.tmp, "subreads_*.tmp") )
		for fn in subreads_fns:
			if   fn.find("_labels.tmp")>-1:
				labels_fn     = fn
			elif fn.find("_lengths.tmp")>-1:
				lengths_fn    = fn
			elif fn.find("_names.tmp")>-1:
				readnames_fn  = fn
			elif fn.find("_ipds.tmp")>-1:
				ipds_fn       = fn
			elif fn.find("_ipdsN.tmp")>-1:
				ipds_N_fn     = fn
			elif fn.find("_compN.tmp")>-1:
				comp_N_fn     = fn
			elif fn.find("_compkmers.tmp")>-1:
				comp_kmers_fn = fn
			elif fn.find("_ipdskmers.tmp")>-1:
				ipds_kmers_fn = fn

		readname_row_idx = defaultdict(list)
		for i,line in enumerate(open(readnames_fn, "r").xreadlines()):
			readname = line.strip()
			readname_row_idx[readname].append(i)
		n_subreads = i+1

		for line in open(ipds_kmers_fn).xreadlines():
			n_motifs = len(line.strip().split("\t"))

		sublengths = np.loadtxt(lengths_fn,    dtype="int")
		ipds       = np.loadtxt(ipds_fn,       dtype="float")
		ipds_N     = np.loadtxt(ipds_N_fn,     dtype="int")
		ipds_kmers = np.loadtxt(ipds_kmers_fn, dtype="string")
		comp_N     = np.loadtxt(comp_N_fn,     dtype="int")
		comp_kmers = np.loadtxt(comp_kmers_fn, dtype="string")
		labels     = np.loadtxt(labels_fn,     dtype="string")

		if n_subreads>1:
			if n_motifs==1:
				# Only one motif survived motif-filtering
				# Case #3
				ipds       = ipds.reshape(ipds.shape[0],1)
				ipds_N     = ipds_N.reshape(ipds_N.shape[0],1)
				ipds_kmers = ipds_kmers.reshape(1)
			elif n_motifs>1:
				# The ipd data is loaded properly in matrix form
				# Case #1
				pass
			elif n_motifs==0:
				pass
		elif n_subreads==1:
			sublengths = sublengths.reshape(1)
			comp_N     = comp_N.reshape(1,comp_N.shape[0])
			if n_motifs==1:
				# Case #4
				ipds       = ipds.reshape(1,1)
				ipds_N     = ipds_N.reshape(1,1)
				ipds_kmers = ipds_kmers.reshape(1)
			elif n_motifs>1:
				# Case #2
				ipds       = ipds.reshape(1,ipds.shape[0])
				ipds_N     = ipds_N.reshape(1,ipds_N.shape[0])
			elif n_motifs==0:
				pass
		
		if tmp_run:
			reads_names_fn      = os.path.join(self.opts.tmp, "read_names.tmp"     )
			reads_labels_fn     = os.path.join(self.opts.tmp, "read_labels.tmp"    )
			reads_lengths_fn    = os.path.join(self.opts.tmp, "read_lengths.tmp"   )
			reads_ipds_fn       = os.path.join(self.opts.tmp, "read_ipds.tmp"      )
			reads_ipds_N_fn     = os.path.join(self.opts.tmp, "read_ipdsN.tmp"    )
			reads_ipds_kmers_fn = os.path.join(self.opts.tmp, "read_ipdskmers.tmp")
			reads_comp_N_fn     = os.path.join(self.opts.tmp, "read_compN.tmp"    )
			reads_comp_kmers_fn = os.path.join(self.opts.tmp, "read_compkmers.tmp")
		else:
			reads_names_fn      = self.read_names_fn     
			reads_labels_fn     = self.read_labels_fn    
			reads_lengths_fn    = self.read_lengths_fn   
			reads_ipds_fn       = self.read_SMp_fn      
			reads_ipds_N_fn     = self.read_SMp_counts_fn    
			reads_ipds_kmers_fn = self.read_SMp_kmers_fn
			reads_comp_N_fn     = self.read_comp_fn    
			reads_comp_kmers_fn = self.read_comp_kmers_fn

		f_reads      = open(reads_names_fn,      "w")
		f_labels     = open(reads_labels_fn,     "w")
		f_lengths    = open(reads_lengths_fn,    "w")
		f_ipds       = open(reads_ipds_fn,       "w")
		f_ipds_N     = open(reads_ipds_N_fn,     "w")
		f_comp_N     = open(reads_comp_N_fn,     "w")
		for readname,row_idx in readname_row_idx.iteritems():
			f_reads.write(  "%s\n" % readname)
			f_labels.write( "%s\n" % labels[row_idx][0])
			comp_N_list  = comp_N[row_idx,:].sum(axis=0)
			ipds_list    = ipds[row_idx,:].mean(axis=0)
			ipds_N_list  = ipds_N[row_idx,:].sum(axis=0)

			# Normalize composition kmer counts
			normed_comp_N_list = map(lambda x: math.log( float(x)/sum(comp_N_list) ), comp_N_list)
			readlength         = sublengths[row_idx].sum()
			f_lengths.write( "%s\n" % readlength)
			f_ipds.write(    "%s\n" % "\t".join(map(lambda x: str(round(x,3)), ipds_list)))
			f_ipds_N.write(  "%s\n" % "\t".join(map(lambda x: str(x),          ipds_N_list)))
			f_comp_N.write(  "%s\n" % "\t".join(map(lambda x: str(round(x,6)), normed_comp_N_list)))
		
		shutil.copy(ipds_kmers_fn, reads_ipds_kmers_fn)
		shutil.copy(comp_kmers_fn, reads_comp_kmers_fn)
		f_reads.close()
		f_labels.close()
		f_lengths.close()
		f_ipds.close()
		f_ipds_N.close()
		f_comp_N.close()

		# Remove the subread-level barcodes for each contig
		# for fn in subreads_fns:
		# 	os.remove(fn)

	def bas_stream_files( self ):
		# control_ipds_fn      = os.path.join(self.opts.control_dir, self.opts.tmp, "read_ipds.tmp")
		# control_ipds_N_fn    = os.path.join(self.opts.control_dir, self.opts.tmp, "read_ipdsN.tmp")
		
		reads_ipds_fn        = os.path.join( self.opts.tmp, "read_ipds.tmp"      )
		reads_ipds_kmers_fn  = os.path.join( self.opts.tmp, "read_ipdskmers.tmp")
		reads_ipds_N_fn      = os.path.join( self.opts.tmp, "read_ipdsN.tmp"    )
		all_motifs           = defaultdict(list)
		logging.info("Unpickling the control IPDs...")
		control_means        = pickle.load(open(self.opts.control_pkl, "r"))
		logging.info("Done.")
		args = []
		for j,line in enumerate(open(reads_ipds_fn+".trans", "r").xreadlines()):
			args.append( (j, copy.copy(control_means), reads_ipds_fn, reads_ipds_N_fn, reads_ipds_kmers_fn, self.opts.minMotifIPD, self.opts.min_motif_N, self.opts.min_motif_reads) )

		results          = launch_pool( self.opts.procs, get_motif_scores_from_read, args )
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

	def get_read_refs_from_SAM( self ):
		read_refs = {}
		read_lens = {}
		samfile   = pysam.AlignmentFile(self.opts.sam, "rb")
		for k,aln in enumerate(samfile.fetch()):
			if aln.mapq == 254:
				read            = "/".join(aln.query_name.split("/")[:2])
				read            = "_".join(read.split("_")[1:])
				# ref             = aln.reference_name.split("|")[0]
				ref             = cmph5_read.slugify(aln.reference_name)
				read_refs[read] = ref
				read_lens[read] = len(aln.seq)

		# Write the read_refs file to match the readnames file
		readnames = np.loadtxt(self.read_names_fn, dtype="str")
		f_names   = open(self.read_refs_fn,    "wb")
		f_lens    = open(self.read_lengths_fn, "wb")
		for readname in readnames:
			try:
				f_names.write("%s\n" % read_refs[readname])
			except KeyError:
				f_names.write("unknown\n")
			try:
				f_lens.write("%s\n" % read_lens[readname])
			except KeyError:
				f_lens.write("0\n")
		f_names.close()
		f_lens.close()

	def write_contig_length_labels_file( self, cmph5_file ):
		contigs           = np.atleast_1d(np.loadtxt(self.contig_names_fn, dtype="str"))
		f                 = open(self.contig_lengths_fn, "w")
		for contig in contigs:
			if self.opts.sam!=None:
				if contig=="unknown":
					length = 1000
				else:
					length = self.contig_fasta_lens[contig]
			else:
				length = self.opts.cmph5_contig_lens[cmph5_file][contig]
			f.write("%s\n" % length)
		f.close()

	def write_readname_species_map( self, readname_species_fn ):
		readname_species_map = {}
		if os.path.exists(readname_species_fn):
			for line in open(readname_species_fn, "r").xreadlines():
				line = line.strip()
				spec = line.split("\t")[1]
				read = line.split("\t")[0]
				readname_species_map[read] = spec
		else:
			readnames = np.loadtxt(self.read_names_fn, dtype="str")
			for read in readnames:
				readname_species_map[read] = "mix"
		return readname_species_map

	def write_movie_species_map( self, movie_map_fn ):
		cell_species_map = {}
		if os.path.exists(movie_map_fn):
			for line in open(movie_map_fn, "r").xreadlines():
				line  = line.strip()
				spec  = line.split(":")[0]
				cell  = line.split(":")[1].split(".")[0]
				cell_species_map[cell] = spec
		else:
			readnames = np.loadtxt(self.read_names_fn, dtype="str")
			cellnames = set(map(lambda x: x.split("/")[0], readnames))
			for cellname in cellnames:
				cell_species_map[cellname] = "mix"
		return cell_species_map

	def rewrite_read_labels_synth_mix( self, cell_species_map ):
		os.rename(self.read_labels_fn, self.read_labels_fn+".tmp")
		f = open(self.read_labels_fn, "w")
		for readname in open(self.read_names_fn, "r").xreadlines():
			readname = readname.strip()
			cell     = readname.split("/")[0]
			new_lab  = cell_species_map[cell]
			f.write("%s\n" % new_lab)
		f.close()
		os.remove(self.read_labels_fn+".tmp")

	def rewrite_read_labels_real_mix( self, readname_species_map ):
		os.rename(self.read_labels_fn, self.read_labels_fn+".tmp")
		f = open(self.read_labels_fn, "w")
		for readname in open(self.read_names_fn, "r").xreadlines():
			readname = readname.strip()
			try:
				new_lab = readname_species_map[readname]
			except KeyError:
				new_lab = "unknown"
			f.write("%s\n" % new_lab)
		f.close()
		os.remove(self.read_labels_fn+".tmp")

	def rewrite_contig_labels_synth_mix( self, cell_species_map ):
		readnames    = np.loadtxt(self.read_names_fn,     dtype="str")
		contigs      = np.atleast_1d(np.loadtxt(self.contig_names_fn,   dtype="str"))
		refs         = np.loadtxt(self.read_refs_fn,      dtype="str")
		spec_counter = {}
		for contig in set(refs):
			spec_counter[contig] = Counter()
			idx                  = refs==contig
			ref_reads            = readnames[idx]
			for readname in ref_reads:
				cell   = readname.split("/")[0]
				spec   = cell_species_map[cell]
				spec_counter[contig][spec] += 1
		percents  = {}
		f_mapping = open("contig_species_mappings.out", "w")
		specs     = set()
		for contig,counter in spec_counter.iteritems():
			for spec,N in counter.iteritems():
				specs.add(spec)
		header  = "contig_name\t"
		header += "\t".join(specs)
		f_mapping.write("%s\n" % header)
		for contig,counter in spec_counter.iteritems():
			pcts_str         = contig
			percents[contig] = {}
			for spec in specs:
				if counter.get(spec):
					pct = 100 * float(counter[spec]) / sum(counter.values())
				else:
					pct = 0.0
				percents[contig][spec] = pct
				pcts_str += "\t%.1f" % pct
			f_mapping.write("%s\n" % pcts_str)
		f_mapping.close()
		os.rename(self.contig_labels_fn, self.contig_labels_fn+".tmp")
		f         = open(self.contig_labels_fn, "w")
		for contig in contigs:
			max_spec = max(percents[contig].iteritems(), key=operator.itemgetter(1))[0]
			f.write("%s\n" % max_spec)
		f.close()
		os.remove(self.contig_labels_fn+".tmp")

	def rewrite_contig_labels_real_mix( self, readname_species_map ):
		readnames    = np.loadtxt(self.read_names_fn,   dtype="str")
		contigs      = np.atleast_1d(np.loadtxt(self.contig_names_fn, dtype="str"))
		refs         = np.loadtxt(self.read_refs_fn,    dtype="str")
		spec_counter = {}
		for contig in set(contigs):
			spec_counter[contig] = Counter()
			idx                  = refs==contig
			ref_reads            = readnames[idx]
			for readname in ref_reads:
				readname = readname.strip()
				try:
					spec = readname_species_map[readname]
				except KeyError:
					spec = "unknown"
				spec_counter[contig][spec] += 1
		percents = {}
		for contig,counter in spec_counter.iteritems():
			percents[contig] = {}
			for spec,N in counter.iteritems():
				pct = 100 * float(N) / sum(counter.values())
				percents[contig][spec] = pct
		os.rename(self.contig_labels_fn, self.contig_labels_fn+".tmp")
		f = open(self.contig_labels_fn, "w")
		for contig in contigs:
			if len(percents[contig])==0:
				max_spec = "unknown" 
			else:
				max_spec = max(percents[contig].iteritems(), key=operator.itemgetter(1))[0]
			f.write("%s\n" % max_spec)
		f.close()
		os.remove(self.contig_labels_fn+".tmp")

	def digest_large_contigs_for_tSNE( self ):
		names      = np.atleast_1d(np.loadtxt(self.contig_names_fn,     dtype="str"))
		labels     = np.atleast_1d(np.loadtxt(self.contig_labels_fn,    dtype="str"))
		covs       = np.atleast_1d(np.loadtxt(self.contig_cov_fn,       dtype="float"))
		lengths    = np.atleast_1d(np.loadtxt(self.contig_lengths_fn,   dtype="int"))
		SCp        = np.atleast_1d(np.loadtxt(self.contig_SCp_fn,       dtype="float"))
		SCp_N      = np.atleast_1d(np.loadtxt(self.contig_SCp_N_fn,     dtype="float"))
		comp       = np.atleast_1d(np.loadtxt(self.contig_comp_fn,      dtype="float"))

		# if len(SCp.shape)==1:
		# 	SCp   = SCp.reshape(SCp.shape[0],1)
		# 	SCp_N = SCp_N.reshape(SCp_N.shape[0],1)

		# if len(comp.shape)==1:
		# 	comp   = comp.reshape(comp.shape[0],1)

		max_size          = self.opts.subcontig_size
		ref_comp_fn       = "contigs.frag_comp"
		ref_SCp_fn        = "contigs.frag_SCp"
		ref_SCp_N_fn      = "contigs.frag_SCp_N"
		ref_names_fn      = "contigs.frag_names"
		ref_labels_fn     = "contigs.frag_labels"
		ref_covs_fn       = "contigs.frag_covs"
		ref_lengths_fn    = "contigs.frag_lengths"
		f_ref_comp        = open(ref_comp_fn,    "w")
		f_ref_SCp         = open(ref_SCp_fn,     "w")
		f_ref_SCp_N       = open(ref_SCp_N_fn,   "w")
		f_ref_names       = open(ref_names_fn,   "w")
		f_ref_labels      = open(ref_labels_fn,  "w")
		f_ref_covs        = open(ref_covs_fn,    "w")
		f_ref_lengths     = open(ref_lengths_fn, "w")
		for i,length in enumerate(lengths):
			if length>max_size:
				name            = names[i]
				label           = labels[i]
				cov             = covs[i]
				# Copy the original contig SCp and composition vectors; apply to fragments
				if len(SCp.shape)==1:
					contig_SCp_str    = "\t".join(map(lambda x: str(x), SCp[:]))
					contig_SCp_N_str  = "\t".join(map(lambda x: str(x), SCp_N[:]))
					contig_comp_str   = "\t".join(map(lambda x: str(x), comp[:]))
				else: 
					contig_SCp_str    = "\t".join(map(lambda x: str(x), SCp[i,:]))
					contig_SCp_N_str  = "\t".join(map(lambda x: str(x), SCp_N[i,:]))
					contig_comp_str   = "\t".join(map(lambda x: str(x), comp[i,:]))
				
				for j,(refName,contig_seq) in enumerate(self.fasta_iter(self.opts.contigs)):
				
					if refName.find("|quiver")>-1:
						# SMRT assemblies add |quiver to contig names, but this
						# gets dropped from the contig names in the cmp.h5 file.
						refName = refName.replace("|quiver","")
			
					if cmph5_read.slugify(refName) == name:
						seq_chunks = list(self.chunks(contig_seq, max_size))
						# Omit the chunk containing the remainder sequence
						seq_chunks = seq_chunks[:-1]
						for seq_chunk in seq_chunks:
							contig_comp  = []
							contig_kmers = read_scanner.kmer_freq( "cmp", seq_chunk, 0, self.opts )
							for kmer,count in contig_kmers.iteritems():
								kmer_normed_comp = math.log(float(count) / sum(contig_kmers.values()))
								contig_comp.append(kmer_normed_comp)
							# contig_comp_str = "\t".join(map(lambda x: str(round(x,6)), contig_comp))
							f_ref_comp.write(   "%s\n" %  contig_comp_str)
							f_ref_SCp.write(    "%s\n" %  contig_SCp_str)
							f_ref_SCp_N.write(  "%s\n" %  contig_SCp_N_str)
							f_ref_names.write(  "%s\n" %  name)
							f_ref_labels.write( "%s\n" %  label)
							f_ref_covs.write(   "%s\n" %  cov)
							f_ref_lengths.write("%s\n" %  len(seq_chunk))
		f_ref_comp.close()
		f_ref_SCp.close()
		f_ref_SCp_N.close()
		f_ref_names.close()
		f_ref_labels.close()
		f_ref_covs.close()
		f_ref_lengths.close()

		def append_file( orig, new ):
			cat_CMD = "cat %s %s > tmp.appended" % (orig, new)
			sts, stdOutErr = run_OS_command( cat_CMD )
			os.rename("tmp.appended", orig)
			os.remove(new)
			
		append_file( self.contig_comp_fn,    ref_comp_fn )
		append_file( self.contig_SCp_fn,     ref_SCp_fn )
		append_file( self.contig_SCp_N_fn,   ref_SCp_N_fn )
		append_file( self.contig_names_fn,   ref_names_fn )
		append_file( self.contig_labels_fn,  ref_labels_fn )
		append_file( self.contig_cov_fn,     ref_covs_fn )
		append_file( self.contig_lengths_fn, ref_lengths_fn )

	def run( self ):
		##########################################################
		# Create set of motifs to analyze for methylation profiles
		##########################################################
		if self.opts.motifs_file != None:
			self.motifs, self.bi_motifs = motif_tools.motifs_from_file( self.opts )
		else:
			self.motifs, self.bi_motifs = motif_tools.build_motif_dict( self.opts )

		if os.path.exists(self.opts.tmp):
			shutil.rmtree(self.opts.tmp)
		os.mkdir(self.opts.tmp)

		##########################################################
		# Define output filenames
		##########################################################
		if self.opts.h5_type=="cmp":
			#############################################
			# Define the alignment-level output filenames
			#############################################
			self.read_names_fn           = "reads.names"
			self.read_refs_fn            = "reads.refs"
			self.read_labels_fn          = "reads.labels"
			self.read_lengths_fn         = "reads.lengths"
			self.read_strands_fn         = "reads.strands"
			self.read_comp_kmers_fn      = "reads.comp_kmers"
			self.read_comp_counts_fn     = "reads.comp_counts"
			self.read_comp_fn            = "reads.comp"
			self.read_comp_2D_fn         = "reads.comp.2D"
			self.read_comp_2D_z_fn       = "reads.comp.2D.zscores"
			self.read_SMp_kmers_fn       = "reads.SMp_kmers"
			self.read_SMp_counts_fn      = "reads.SMp_counts"
			self.read_SMp_fn             = "reads.SMp"
			self.read_SMp_2D_fn          = "reads.SMp.2D"
			self.read_SMp_2D_z_fn        = "reads.SMp.2D.zscores"
			self.read_4D_z_fn            = "reads.combined.4D.zscores"
			self.read_combo_2D_fn        = "reads.combined.2D"

			# #############################################
			# # Define the contig-level output filenames
			# #############################################
			self.contig_names_fn         = "contigs.names"
			self.contig_labels_fn        = "contigs.labels"
			self.contig_lengths_fn       = "contigs.lengths"
			self.contig_comp_fn          = "contigs.comp"
			self.contig_comp_2D_fn       = "contigs.comp.2D"
			self.contig_comp_2D_z_fn     = "contigs.comp.2D.zscores"
			self.contig_cov_comp_3D_fn   = "contigs.cov_comp.3D"
			self.contig_cov_comp_2D_fn   = "contigs.cov_comp.2D"
			self.contig_cov_comp_2D_z_fn = "contigs.cov_comp.2D.zscores"
			self.contig_SCp_fn           = "contigs.SCp"
			self.contig_SCp_N_fn         = "contigs.SCp_N"
			self.contig_SCp_2D_fn        = "contigs.SCp.2D"
			self.contig_SCp_2D_z_fn      = "contigs.SCp.2D.zscores"
			# self.contig_SCp_kmers_fn     = "contigs.SCp_kmers"
			self.contig_cov_fn           = "contigs.cov"
			self.contig_4D_z_fn          = "contigs.combined.4D.zscores"
			self.contig_combo_2D_fn      = "contigs.combined.2D"
		
			if self.opts.cross_cov_bins!=None:
				self.bin_names_fn        = "bins.names"
				self.bin_sizes_fn        = "bins.sizes"
				self.bin_SCp_fn          = "bins.SCp"
				self.bin_SCp_N_fn        = "bins.SCp_N"

		elif self.opts.h5_type=="bas":
			
			self.read_names_fn           = "reads.raw.names"
			self.read_labels_fn          = "reads.raw.labels"
			self.read_lengths_fn         = "reads.raw.lengths"
			self.read_comp_kmers_fn      = "reads.raw.comp_kmers"
			self.read_comp_counts_fn     = "reads.raw.comp_counts"
			self.read_comp_fn            = "reads.raw.comp"
			self.read_comp_2D_fn         = "reads.raw.comp.2D"
			self.read_comp_2D_z_fn       = "reads.raw.comp.2D.zscores"
			self.read_SMp_kmers_fn       = "reads.raw.SMp_kmers"
			self.read_SMp_counts_fn      = "reads.raw.SMp_counts"
			self.read_SMp_fn             = "reads.raw.SMp"
			self.read_SMp_2D_fn          = "reads.raw.SMp.2D"
			self.read_SMp_2D_z_fn        = "reads.raw.SMp.2D.zscores"
			self.read_4D_z_fn            = "reads.raw.combined.4D.zscores"
			self.read_combo_2D_fn        = "reads.raw.combined.2D"
			
			if self.opts.sam!=None:
				####################################################
				# Define contig filenames for SAM --> contig pathway
				####################################################
				self.read_refs_fn            = "reads.ccs.refs"
				self.contig_names_fn         = "contigs.names"
				self.contig_labels_fn        = "contigs.labels"
				self.contig_lengths_fn       = "contigs.lengths"
				self.contig_comp_fn          = "contigs.comp"
				self.contig_comp_2D_fn       = "contigs.comp.2D"
				self.contig_comp_2D_z_fn     = "contigs.comp.2D.zscores"
				self.contig_cov_comp_3D_fn   = "contigs.cov_comp.3D"
				self.contig_cov_comp_2D_fn   = "contigs.cov_comp.2D"
				self.contig_cov_comp_2D_z_fn = "contigs.cov_comp.2D.zscores"
				self.contig_SCp_fn           = "contigs.SCp"
				self.contig_SCp_N_fn         = "contigs.SCp_N"
				self.contig_SCp_2D_fn        = "contigs.SCp.2D"
				self.contig_SCp_2D_z_fn      = "contigs.SCp.2D.zscores"
				# self.contig_SCp_kmers_fn     = "contigs.SCp_kmers"
				self.contig_cov_fn           = "contigs.cov"
				self.contig_4D_z_fn          = "contigs.combined.4D.zscores"
				self.contig_combo_2D_fn      = "contigs.combined.2D"

		self.opts.h5_labels         = {}
		self.opts.cmph5_contig_lens = {}

		##########################################################
		# Point to or build the control IPD data using WGA sequencing
		##########################################################
		filter_N_reads = min(self.opts.N_motif_reads, self.opts.N_reads)

		# wga_cmp_h5 = "/hpc/users/beaulj01/projects/ebinning/control/P6_C4_Bb_B31_aligned_reads.cmp.h5"
		wga_cmp_h5 = "/hpc/users/beaulj01/projects/ebinning/control/P6_C4_Cdiff_PID22134_aligned_reads.cmp.h5"
		wga_bas_h5 = "/hpc/users/beaulj01/projects/ebinning/gnoto/isolates/rgnavus/m150430_033851_42163R_c100796292550000001823171709091517_s1_p0.1.bax.h5"
		if self.opts.h5_type=="cmp":
			control_h5 = wga_cmp_h5
		elif self.opts.h5_type=="bas":
			control_h5 = wga_bas_h5

		self.opts.control_run = True

		logging.info("Checking for existing control IPD data...")
		self.opts.control_pkl = os.path.join(self.opts.control_dir, self.opts.tmp, "control_means.pkl")
		controls              = ControlRunner(control_h5, self.opts)
		launch_control        = controls.check_control_file()
		
		if launch_control:
			logging.info("No control data found -- preparing to create new control data...")
			logging.info("   -- Being created in %s" % self.opts.control_dir)
			controls.goto_control_output_dir()
			self.opts = controls.scan_WGA_h5()
			self.launch_data_loader( control_h5, filter_N_reads, 1 )
			controls.analyze_WGA_reads()
			logging.info("Done.")

			logging.info("Building dictionary of control values for all motifs...")
			logging.info("   * Initial build requires significant time and memory.")
			control_means = controls.build_control_IPD_dict(self.motifs, self.bi_motifs)
			logging.info("Done.")
		else:
			logging.info("Control data found at %s" % self.opts.control_dir)
			
			logging.info("Loading pickled control IPD values...")
			control_means = controls.load_control_pickle()
			logging.info("Done.")

		if self.opts.motifs_file!=None:
			# Degen motifs in the motifs_file might not be in the existing
			# dictionary of control values. Need to construct it by combining
			# the existing data from the various specified versions of the motif
			
			control_means = controls.add_degen_motifs( self.opts.motifs_file, control_means)
			local_pkl     = os.path.basename(self.opts.control_pkl)+".wdegen"
			pickle.dump(control_means, open(local_pkl, "wb"))

		controls.return_to_orig_dir()
		self.opts.control_run = False
		
		###############################
		# *.h5 file metadata collection
		###############################
		h5_files                    = []
		for line in open(self.opts.h5_files).xreadlines():
			line                     = line.strip("\n")
			h5_fn                    = line.split()[0]
			label                    = line.split()[1]
			self.opts.h5_labels[h5_fn]         = label
			h5_files.append(h5_fn)
			if self.opts.h5_type=="cmp":
				self.opts.cmph5_contig_lens[h5_fn] = {}
				logging.info("Getting contig information from %s..." % h5_fn)
				reader = CmpH5Reader(h5_fn)
				for entry in reader.referenceInfoTable:
					name                                = entry[3]
					length                              = entry[4]
					slug_name                           = cmph5_read.slugify(name)
					self.opts.cmph5_contig_lens[h5_fn][slug_name] = length
				reader.close()
				logging.info("Done.")

		if self.opts.h5_type=="cmp" and len(h5_files)>1:
			raise Exception("For aligned read analysis, only one cmp.h5 file at a time is supported.")

		if not self.opts.comp_only and self.opts.motifs_file == None:
			#########################################################
			# Launch analysis of <filter_N_reads> for motif filtering
			#########################################################
			# filter_N_reads        = min(5000, self.opts.N_reads)
			filter_N_reads        = min(20000, self.opts.N_reads)
			self.opts.control_run = False
			for i,h5_file in enumerate(h5_files):
				logging.info("Creating %s barcodes (%s motifs) from %s..." % (filter_N_reads, (len(self.motifs)+len(self.bi_motifs)), h5_file))
				self.launch_data_loader( h5_file, filter_N_reads, i )
				logging.info("Done.")

			if self.opts.h5_type=="bas":
				# Combine subread data across multiple movies
				logging.info("Combining subread data across all movies...")
				results = self.combine_subread_data_across_bas_movies()
				logging.info("Done.")
				# Combine movie-merged subreads to get read-level barcodes
				logging.info("Combining subreads to get read-level barcodes...")
				results = self.bas_combine_subreads_for_read_level( tmp_run=True )
				logging.info("Done.")
				
			####################################################
			# Filter out motifs without significant signatures
			####################################################
			logging.info("Getting top motifs from each contig...")
			# if self.opts.h5_type=="cmp":
			# 	# control_ipds_fn   = "/hpc/users/beaulj01/projects/ebinning/control/control_quiver_ipds.tmp"
			# 	# control_ipds_N_fn = "/hpc/users/beaulj01/projects/ebinning/control/control_quiver_ipds_N.tmp"
			# 	control_ipds_fn   = os.path.join(control_dir, self.opts.tmp, "control_ipds.tmp" )
			# 	control_ipds_N_fn = os.path.join(control_dir, self.opts.tmp, "control_ipdsN.tmp")
			# elif self.opts.h5_type=="bas":
			# 	control_ipds_fn   = os.path.join("control", self.opts.tmp, "read_ipds.tmp")
			# 	control_ipds_N_fn = os.path.join("control", self.opts.tmp, "read_ipdsN.tmp")

			if self.opts.h5_type=="cmp":
				
				self.contig_fasta_lens = {}
				for (name,seq) in self.fasta_iter(self.opts.contigs):
					if name.find("|quiver")>-1:
						# SMRT assemblies add |quiver to contig names, but this
						# gets dropped from the contig names in the cmp.h5 file.
						name = name.replace("|quiver","")
					self.contig_fasta_lens[cmph5_read.slugify(name)] = len(seq)

				contig_ipds_fns = glob.glob( os.path.join("tmp", "*_ipds.tmp") )
				contigs         = map(lambda x: os.path.basename(x).split("_ipds.tmp")[0], contig_ipds_fns)
				ipds_fn_dict    = dict([(os.path.basename(ipds_fn).split("_ipds.tmp")[0],ipds_fn) for ipds_fn in contig_ipds_fns])

				contigs_for_transpose = []
				contigs_for_chunking  = []
				# maxsize_for_transpose = 50000000 #50Mb
				maxsize_for_transpose = 25000000 #25Mb
				# maxsize_for_transpose = 10000000 #10Mb
				for name in contigs:
					fsize = os.path.getsize(ipds_fn_dict[name])
					if fsize < maxsize_for_transpose:
						contigs_for_transpose.append(name)
					else:
						contigs_for_chunking.append(name)

				logging.info("Transposing %s case contigs..." % len(contigs_for_transpose))
				args    = [(contig,False) for contig in contigs_for_transpose]
				results = launch_pool( self.opts.procs, transpose_contig_matrix, args )
				logging.info("Done.")

				streamed_contig_dicts  = {}
				if len(contigs_for_transpose)>0:
					logging.info("Streaming through %s contigs..." % len(contigs_for_transpose))
					args    = [(self.opts.control_pkl, contig, i, len(contigs_for_transpose), self.opts.minMotifIPD, self.opts.min_motif_N) for i,contig in enumerate(contigs_for_transpose)]
					results = launch_pool( self.opts.procs, stream_case_control_files, args)
					
					streamed_contig_SCp    = map(lambda x: x[0], results)
					streamed_contig_SCp_N  = map(lambda x: x[1], results)
					streamed_contigs       = map(lambda x: x[2], results)
					
					for i,contig in enumerate(streamed_contigs):
						streamed_contig_dicts[contig] = {"SCp":streamed_contig_SCp[i], "SCp_N":streamed_contig_SCp_N[i]}

					logging.info("Done.")
				
				chunked_contigs_dicts  = {}
				if len(contigs_for_chunking)>0:
					logging.info("Chunking %s contigs..." % len(contigs_for_chunking))
					
					for i,contig in enumerate(contigs_for_chunking):
						control_means,contig_SCp,contig_SCp_N,contig = chunk_case_control_files( self.opts.control_pkl, contig, i, len(contigs_for_chunking), self.opts.minMotifIPD, self.opts.min_motif_N, self.opts.procs, self.opts )
						
						chunked_contigs_dicts[contig] = {"SCp":contig_SCp, "SCp_N":contig_SCp_N}
					
					logging.info("Done.")

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
							
							args.append( (bin_id, \
										  bin_copy_contig_dicts["SCp"],   \
										  bin_copy_contig_dicts["SCp_N"], \
										  self.opts.minMotifIPD,          \
										  self.opts.min_motif_N,          \
										  len(bin_ids),                   \
										  bin_id,                         \
										  self.opts.control_pkl,          \
										  "bin") )

					results = launch_pool( self.opts.procs, simplify_motifs, args )

					bin_keeper_motifs_list = map(lambda x: x[0], results)
					control_means_list     = map(lambda x: x[1], results)

					"""
					Add the control means for these bin motifs to the 
					complete set of control means for all detected motifs.
					"""
					for bin_keeper_motifs in bin_keeper_motifs_list:
						keeper_motifs = keeper_motifs | bin_keeper_motifs

					for control_means in control_means_list:
						for motif,score in control_means.iteritems():
							keeper_control_ipds[motif] = score

				else:

					args = []
					for j,(contig,contig_d) in enumerate(contig_dicts.iteritems()):
						# For each contig, do motif filtering and refinement
					
						copy_contig_dicts = copy.deepcopy(contig_d)

						args.append( (j,                          \
									  copy_contig_dicts["SCp"],   \
									  copy_contig_dicts["SCp_N"], \
									  self.opts.minMotifIPD,      \
									  self.opts.min_motif_N,      \
									  len(contig_dicts.keys()),   \
									  contig,                     \
									  self.opts.control_pkl,      \
									  "contig") )

					results = launch_pool( self.opts.procs, simplify_motifs, args )

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
				# But write the control pickle locally so that it doesn't get overwritten
				# if another job uses the same control directory.
				local_pkl = os.path.basename(self.opts.control_pkl)+".wdegen"
				pickle.dump( keeper_control_ipds, open(local_pkl, "wb" ) )
			
			elif self.opts.h5_type=="bas":
				# if os.path.exists(control_ipds_fn+".trans") and os.path.exists(control_ipds_N_fn+".trans"):
				# 	pass
				# else:
				# 	logging.info("Transposing control matrices....")
				# 	control_files = [control_ipds_fn, control_ipds_N_fn]
				# 	results       = launch_pool( len(control_files), transpose_file, control_files )
				# 	logging.info("Done.")
				
				logging.info("Transposing reads...")
				files   = [os.path.join(self.opts.tmp, "read_ipds.tmp"), os.path.join(self.opts.tmp, "read_ipdsN.tmp")]
				results = launch_pool( len(files), transpose_file, files )
				logging.info("Done.")
				logging.info("Streaming through reads for motif filtering...")
				keeper_motifs = self.bas_stream_files()
				
				keeper_control_ipds = {}
				control_means       = pickle.load(open(self.opts.control_pkl, "r"))
				for motif in list(keeper_motifs):
					keeper_control_ipds[motif] = control_means[motif]
				logging.info("%s" % ",".join(list(keeper_motifs)))
				logging.info("Done.")

				local_pkl = os.path.basename(self.opts.control_pkl)+".filtered"
				pickle.dump( keeper_control_ipds, open(local_pkl, "wb" ) )

			logging.info("Keeping %s motifs for further analysis:" % len(keeper_motifs))
			self.motifs = list(keeper_motifs)
			motifs_fn   = "motifs.txt"
			n_motifs    = len(keeper_motifs)
			f_motifs    = open(motifs_fn, "w")
			if n_motifs == 0 and not self.opts.motif_discov_only:
				keeper_motifs = set(["GATC-1"])
			for motif in keeper_motifs:
				f_motifs.write("%s\n" % motif)
			f_motifs.close()

			if self.opts.motif_discov_only:
				logging.info("Discovered %s motifs -- EXITING!" % len(keeper_motifs))
				sys.exit()

		else:
			self.opts.control_run = False
			if self.opts.motifs_file != None:
				keeper_motifs = set(np.loadtxt(self.opts.motifs_file, dtype="str"))
			else:
				keeper_motifs = set(["GATC-1","CATG-1"])
			self.motifs = list(keeper_motifs)
			motifs_fn   = "motifs.txt"
			n_motifs    = len(keeper_motifs)
			f_motifs    = open(motifs_fn, "w")
			for motif in keeper_motifs:
				f_motifs.write("%s\n" % motif)
			f_motifs.close()





		##################################################
		# Re-run analysis using the filtered set of motifs
		##################################################
		logging.info("Calling the pipeline using only %s filtered motifs..." % len(keeper_motifs))
		to_del = glob.glob( os.path.join(self.opts.tmp, "*") )
		for fn in to_del:
			os.remove(fn)
		self.opts.motifs_file = motifs_fn
		
		for i,h5_file in enumerate(h5_files):
			logging.info("Creating subread-level barcodes (%s motifs) from %s..." % (len(keeper_motifs), h5_file))
			self.launch_data_loader(  h5_file, self.opts.N_reads, i )
			logging.info("Done.")


			if self.opts.h5_type=="cmp":
				logging.info("Combining subread-level barcodes to get read-level barcodes from each contig...")
				contig_labels_fns = glob.glob( os.path.join(self.opts.tmp, "*_labels.tmp") )
				contigs           = map(lambda x: os.path.basename(x).split("_labels.tmp")[0], contig_labels_fns)
				args    = [ (h5_file, contig, self.opts.tmp, self.opts.h5_labels, i, len(contigs)) for i,contig in enumerate(contigs)]
				results = launch_pool( self.opts.procs, combine_subreads_for_read_level, args )


				logging.info("Combining read-level barcodes from all contigs...")
				self.combine_read_level_barcodes_across_contigs()
				logging.info("Done.")


				logging.info("Creating contig-level barcodes (%s motifs) from %s..." % (len(keeper_motifs), h5_file))
				self.combine_subreads_for_contig_level( h5_file )
				logging.info("Done.")
				n_contigs = len(np.atleast_1d(np.loadtxt(self.contig_names_fn, dtype="str")))


				if self.opts.cross_cov_bins!=None:
					logging.info("Creating bin-level barcodes (%s motifs) using %s..." % (len(keeper_motifs), self.opts.cross_cov_bins))
					self.combine_contigs_for_bin_level()
					logging.info("Done.")


		if self.opts.h5_type=="bas":
			logging.info("Combining subread data across all movies...")
			results = self.combine_subread_data_across_bas_movies()
			logging.info("Done.")
			logging.info("Combining subreads to get read-level barcodes...")
			results = self.bas_combine_subreads_for_read_level()
			logging.info("Done.")

			if self.opts.sam!=None:
				logging.info("Writing read-contig assignments based on %s..." % self.opts.sam)
				self.get_read_refs_from_SAM()
				logging.info("Done.")
				for i,h5_file in enumerate(h5_files):
					logging.info("Creating contig-level barcodes (%s motifs) from %s..." % (len(keeper_motifs), h5_file))
					self.combine_subreads_for_contig_level( h5_file )
					logging.info("Done.")
				n_contigs = len(np.atleast_1d(np.loadtxt(self.contig_names_fn, dtype="str")))

		movie_map_fn         = "movie_mapping.txt"
		readname_species_fn  = "read_species_mapping.txt"
		readname_species_map = self.write_readname_species_map( readname_species_fn )
		cell_species_map     = self.write_movie_species_map( movie_map_fn )

		################################
		# Apply species mapping to reads
		################################
		if self.opts.synth_mix:
			self.rewrite_read_labels_synth_mix( cell_species_map )
		elif self.opts.real_mix:
			self.rewrite_read_labels_real_mix( readname_species_map )
		
		if self.opts.h5_type=="cmp" or self.opts.sam!=None:
			##################################
			# Apply species mapping to contigs
			##################################
			if self.opts.synth_mix:
				self.rewrite_contig_labels_synth_mix( cell_species_map )
			elif self.opts.real_mix:
				self.rewrite_contig_labels_real_mix( readname_species_map )
			self.write_contig_length_labels_file( h5_files[0] )
			
			##########################################################
			# Reduce contigs composition dimensionality from ND --> 2D
			##########################################################
			logging.info("Reducing dimensionality of the CONTIGS composition matrix...")
			logging.info("    - digesting large contigs to provide more weight in tSNE")
			self.digest_large_contigs_for_tSNE()
			sne_CMD    = "python ~/gitRepo/eBinning/src/bhtsne.py -v -i %s -l %s -o %s -s contigs.comp.scatter.png -z %s -n contigs_composition_digest" %    \
																																	(self.contig_comp_fn,    \
																																	 self.contig_labels_fn,  \
																																	 self.contig_comp_2D_fn, \
																																	 self.contig_lengths_fn)
			self.run_tSNE(sne_CMD)
			
			##################################################################################
			# Adding contigs coverage column to 2D composition matrix, reducing from 3D --> 2D
			##################################################################################
			logging.info("Adding contig coverage column to 2D composition matrix and reducing to 2D...")
			paste_CMD      = "paste %s %s > %s" % (self.contig_comp_2D_fn, self.contig_cov_fn, self.contig_cov_comp_3D_fn)
			sts, stdOutErr = run_OS_command( paste_CMD )
			sne_CMD  = "python ~/gitRepo/eBinning/src/bhtsne.py -v -i %s -l %s -o %s -s contigs.cov_comp.scatter.png -z %s -n contigs_composition_coverage_digest" %    \
																																	(self.contig_cov_comp_3D_fn, \
																																	 self.contig_labels_fn,  \
																																	 self.contig_cov_comp_2D_fn, \
																																	 self.contig_lengths_fn)
			self.run_tSNE(sne_CMD)

			##################################################
			# Reduce contigs IPD dimensionality from ND --> 2D
			##################################################
			logging.info("Reducing dimensionality of the CONTIGS SCp matrix...")
			sne_CMD = "python ~/gitRepo/eBinning/src/bhtsne.py -v -i %s -l %s -o %s -s contigs.SCp.tSNE.scatter.png -z %s -n contigs_SCp_tSNE_%s_motifs_digest" % \
																																	(self.contig_SCp_fn,          \
																																	 self.contig_labels_fn,       \
																																	 self.contig_SCp_2D_fn,       \
																																	 self.contig_lengths_fn,      \
																																	 n_motifs)
			self.run_tSNE(sne_CMD)
			logging.info("Done.")

			#########################################################################
			# Apply Z-score transformation for combining composition and IPD matrices
			#########################################################################
			logging.info("Z-score tranforming the 2D matrices...")
			# for fn in [self.read_SMp_2D_fn, self.read_comp_2D_fn, self.contig_SCp_2D_fn, self.contig_cov_comp_2D_fn, self.contig_comp_2D_fn]:
			for fn in [self.contig_SCp_2D_fn, self.contig_cov_comp_2D_fn, self.contig_comp_2D_fn]:
				if fn==self.contig_SCp_2D_fn and self.opts.h5_type=="cmp" and self.opts.comp_only:
					shutil.copy(self.contig_SCp_fn, self.contig_SCp_2D_fn)
					shutil.copy(self.contig_SCp_2D_fn, self.contig_SCp_2D_fn+".zscores")
				else:
					m         = np.loadtxt(fn, dtype="float")
					m_std     = stats.mstats.zscore(m, axis=0)
					out_fn    = fn + ".zscores"
					np.savetxt(out_fn, m_std, fmt='%.3f', delimiter="\t")
			logging.info("Done.")

			paste_CMD      = "paste %s %s > %s" % (self.contig_cov_comp_2D_z_fn, self.contig_SCp_2D_z_fn, self.contig_4D_z_fn)
			sts, stdOutErr = run_OS_command( paste_CMD )
			sne_CMD        = "python ~/gitRepo/eBinning/src/bhtsne.py -v -i %s -l %s -o %s -s contigs.combined.scatter.png -z %s -n contigs_combined_%s_motifs" % \
																																	(self.contig_4D_z_fn,         \
																																	 self.contig_labels_fn,       \
																																	 self.contig_combo_2D_fn,     \
																																	 self.contig_lengths_fn,      \
																																	 n_motifs)
			self.run_tSNE(sne_CMD)

			########################################################################
			# Make both digested and undigested versions of the contigs output files
			########################################################################
			def drop_appended_lines( fn, n_to_keep ):
				head_CMD = "head -%s %s > tmp.head" % (n_to_keep, fn)
				sts, stdOutErr = run_OS_command( head_CMD )
				os.rename("tmp.head", fn)

			logging.info("...removing the fragmented contig seeds...")
			shutil.copy( self.contig_comp_fn,          self.contig_comp_fn+".digest" )
			shutil.copy( self.contig_comp_2D_fn,       self.contig_comp_2D_fn+".digest" )
			shutil.copy( self.contig_comp_2D_z_fn,     self.contig_comp_2D_z_fn+".digest" )
			shutil.copy( self.contig_cov_comp_3D_fn,   self.contig_cov_comp_3D_fn+".digest" )
			shutil.copy( self.contig_cov_comp_2D_fn,   self.contig_cov_comp_2D_fn+".digest" )
			shutil.copy( self.contig_cov_comp_2D_z_fn, self.contig_cov_comp_2D_z_fn+".digest" )
			shutil.copy( self.contig_SCp_fn,           self.contig_SCp_fn+".digest" )
			shutil.copy( self.contig_SCp_N_fn,         self.contig_SCp_N_fn+".digest" )
			shutil.copy( self.contig_SCp_2D_fn,        self.contig_SCp_2D_fn+".digest" )
			shutil.copy( self.contig_SCp_2D_z_fn,      self.contig_SCp_2D_z_fn+".digest" )
			shutil.copy( self.contig_4D_z_fn,          self.contig_4D_z_fn+".digest" )
			shutil.copy( self.contig_combo_2D_fn,      self.contig_combo_2D_fn+".digest" )
			shutil.copy( self.contig_names_fn,         self.contig_names_fn+".digest" )
			shutil.copy( self.contig_labels_fn,        self.contig_labels_fn+".digest" )
			shutil.copy( self.contig_cov_fn,           self.contig_cov_fn+".digest" )
			shutil.copy( self.contig_lengths_fn,       self.contig_lengths_fn+".digest" )
			drop_appended_lines( self.contig_comp_fn,          n_contigs )
			drop_appended_lines( self.contig_comp_2D_fn,       n_contigs )
			drop_appended_lines( self.contig_comp_2D_z_fn,     n_contigs )
			drop_appended_lines( self.contig_cov_comp_2D_fn,   n_contigs )
			drop_appended_lines( self.contig_cov_comp_2D_z_fn, n_contigs )
			drop_appended_lines( self.contig_cov_comp_3D_fn,   n_contigs )
			drop_appended_lines( self.contig_SCp_fn,           n_contigs )
			drop_appended_lines( self.contig_SCp_N_fn,         n_contigs )
			drop_appended_lines( self.contig_SCp_2D_fn,        n_contigs )
			drop_appended_lines( self.contig_SCp_2D_z_fn,      n_contigs )
			drop_appended_lines( self.contig_4D_z_fn,          n_contigs )
			drop_appended_lines( self.contig_combo_2D_fn,      n_contigs )
			drop_appended_lines( self.contig_names_fn,         n_contigs )
			drop_appended_lines( self.contig_labels_fn,        n_contigs )
			drop_appended_lines( self.contig_cov_fn,           n_contigs )
			drop_appended_lines( self.contig_lengths_fn,       n_contigs )
			
			#########
			warnings.filterwarnings("ignore", category=DeprecationWarning, module="matplotlib")
			#########

			SCp = np.loadtxt(self.contig_SCp_fn, dtype="float")
			if len(SCp.shape)!=1:
				# There is ore than one contig for plotting
				plotting.scatterplot(np.loadtxt(self.contig_comp_2D_fn,     dtype="float"), np.loadtxt(self.contig_labels_fn, dtype="str"), "contigs.comp.noDigest.scatter.png",     np.loadtxt(self.contig_lengths_fn, dtype="int"), "contigs_comp_noDigest")
				plotting.scatterplot(np.loadtxt(self.contig_cov_comp_2D_fn, dtype="float"), np.loadtxt(self.contig_labels_fn, dtype="str"), "contigs.cov_comp.noDigest.scatter.png", np.loadtxt(self.contig_lengths_fn, dtype="int"), "contigs_cov_comp_noDigest")
				plotting.scatterplot(np.loadtxt(self.contig_SCp_2D_fn,      dtype="float"), np.loadtxt(self.contig_labels_fn, dtype="str"), "contigs.SCp.tSNE.noDigest.scatter.png", np.loadtxt(self.contig_lengths_fn, dtype="int"), "contigs_SCp_noDigest")
				plotting.scatterplot(np.loadtxt(self.contig_combo_2D_fn,    dtype="float"), np.loadtxt(self.contig_labels_fn, dtype="str"), "contigs.combined.noDigest.scatter.png", np.loadtxt(self.contig_lengths_fn, dtype="int"), "contigs_combined_noDigest")
			else:
				pass
		
		else:
			#######################################################
			# Reduce reads composition dimensionality from ND --> 2D
			#######################################################
			logging.info("Reducing dimensionality of the READS composition matrix...")
			sne_CMD = "python ~/gitRepo/eBinning/src/bhtsne.py -v -i %s -l %s -o %s -s reads.comp.scatter.png -n reads_composition" %                     \
																																	(self.read_comp_fn,   \
																																	 self.read_labels_fn, \
																																	 self.read_comp_2D_fn)
			self.run_tSNE(sne_CMD)

			###############################################
			# Reduce reads IPD dimensionality from ND --> 2D
			###############################################
			logging.info("Reducing dimensionality of the READS SMp matrix...")
			sne_CMD = "python ~/gitRepo/eBinning/src/bhtsne.py -v -i %s -l %s -o %s -s reads.SMp.tSNE.scatter.png -n reads_SMp_tSNE_%s_motifs_digest" %   \
																																	(self.read_SMp_fn,    \
																																	 self.read_labels_fn, \
																																	 self.read_SMp_2D_fn, \
																																	 n_motifs)
			self.run_tSNE(sne_CMD)

			# logging.info("Combining and tSNE reducing the separately reduced SMp/SCp and composition matrices...")
			# paste_CMD      = "paste %s %s > %s" % (self.read_comp_2D_z_fn, self.read_SMp_2D_z_fn, self.read_4D_z_fn)
			# sts, stdOutErr = run_OS_command( paste_CMD )
			# sne_CMD        = "python ~/gitRepo/eBinning/src/bhtsne.py -v -i %s -l %s -o %s -s reads.combined.scatter.png -n reads_combined" %        \
			# 																														(self.read_4D_z_fn,   \
			# 																														 self.read_labels_fn, \
			# 																														 self.read_combo_2D_fn)
			# self.run_tSNE(sne_CMD)

		
		# shutil.rmtree( self.opts.tmp )

		f = open("ordered_motifs.txt", "w")
		motifs = np.loadtxt(self.read_SMp_kmers_fn, dtype="str")
		for m in motifs:
			f.write("%s\n" % m)
		f.close()

		logging.info("Pipeline finished.")

	def launch_subprocs( self, h5_file, N_reads ):
		"""
		"""
		logging.debug("Creating tasks...")
		tasks   = multiprocessing.JoinableQueue()
		results = multiprocessing.Queue()
		logging.debug("Done.")

		if self.opts.h5_type=="cmp":
			reader     = CmpH5Reader(h5_file)
			to_check   = reader
			entries    = range(len(to_check))
		elif self.opts.h5_type=="bas":
			reader     = BasH5Reader(h5_file)
			if self.opts.bas_whitelist != None and not self.opts.control_run:
				logging.info("Intersecting with whitelist...")
				bas_whitelist = set(np.loadtxt(self.opts.bas_whitelist, dtype="str"))
				to_check = [z for z in reader if z.zmwName in bas_whitelist]
			else:
				to_check = reader
			# Filter on zmw metrics
			pre      = len(to_check)
			logging.info("Starting with %s reads..." % pre)
			to_check = [z for z in to_check if z.zmwMetric("Productivity")==1 and \
											   z.zmwMetric("ReadScore")>self.opts.minReadScore and \
											   z.zmwMetric("Pausiness")<self.opts.maxPausiness]
			post     = len(to_check)
			logging.info("Dropped %s reads due to poor zmw metrics (%s remain)" % ((pre-post),post))
			# Filter on read length 
			pre      = post
			to_check = [z for z in to_check if np.sum([len(sub) for sub in z.subreads])>=self.opts.readlength_min]
			post     = len(to_check)
			logging.info("Dropped %s reads < %s (%s remain)" % ((pre-post),self.opts.readlength_min, post))
			entries  = np.array([z.holeNumber for z in to_check])

		reader.close()
		if len(entries) <= self.opts.procs * 5:
			procs = 1
		else:
			procs = self.opts.procs

		logging.debug("Starting consumers...")
		consumers     = [ multiproc.Consumer(tasks, results) for i in xrange(procs) ]
		for w in consumers:
			w.start()
		logging.debug("Done.")

		num_jobs       = procs
		N_target_reads = {}
		reads_left     = N_reads
		procs_left     = procs

		for job in range(num_jobs):
			N_target_reads[job] = int(math.ceil(float(reads_left)/procs_left))
			reads_left -= N_target_reads[job]
			procs_left -= 1

		logging.debug("Partitioning %s into %s chunks for analysis..." % (h5_file, num_jobs))
		chunksize      = int(math.ceil(float( len(entries)/procs )))
		logging.info("Querying %s reads using %s chunks of size %s..." % (len(to_check), procs, chunksize))
		entries_chunks = list(self.chunks( entries, chunksize ))
		for chunk_id in range(procs):
			n   = N_target_reads[chunk_id]
			# If --N_reads used, this ensures we touch as many contigs as possible
			idx = entries_chunks[chunk_id]
			np.random.shuffle(idx)
			if self.opts.h5_type=="cmp":
				tasks.put(cmph5_read.subread_motif_processor( h5_file, chunk_id, idx, n, self.motifs, self.bi_motifs, self.opts) )
				logging.debug("...%s (%s alignments)" % (chunk_id, len(entries)))
			elif self.opts.h5_type=="bas":
				tasks.put(baxh5_read.subread_motif_processor( h5_file, chunk_id, idx, n, self.motifs, self.bi_motifs, self.opts) )
				logging.debug("...%s (%s reads)" % (chunk_id, len(entries)))
		logging.debug("Done")
		
		# Add a poison pill for each consumer
		for i in xrange(num_jobs):
			tasks.put(None)

		# Wait for all of the tasks to finish
		tasks.join()
		
		# Start printing results
		logging.info("Combining results data from all chunks...")
		parallel_results = []
		while num_jobs:
			result = results.get()
			parallel_results.append(result)
			num_jobs -= 1
			logging.info("...%s/%s" % ((procs-num_jobs),procs))
		logging.info("Done.")
		return parallel_results

	def __parseArgs( self ):
		"""Handle command line argument parsing"""

		usage = """%prog [--help] [options]

		USAGE INSTRUCTIONS
		"""

		parser = optparse.OptionParser( usage=usage, description=__doc__ )

		parser.add_option( "-d", "--debug", action="store_true", help="Increase verbosity of logging" )

		parser.add_option( "-i", "--info", action="store_true", help="Add basic logging" )

		parser.add_option( "--synth_mix", action="store_true", help="The mixture of aligned reads are synthetically created \
																	 by mixing species-specific samples [False]" )

		parser.add_option( "--real_mix", action="store_true", help="The mixture of aligned reads are from a truly mixed sample \
																	of DNA [False]" )

		parser.add_option( "--logFile", type="str", help="Write logging to file [log.out]" )

		parser.add_option( "--subreadlength_min", type="int", help="Minimum subread length to include for analysis [100]" )

		parser.add_option( "--readlength_min", type="int", help="Minimum read length to include for analysis [100]" )

		parser.add_option( "--readlength_max", type="int", help="Maximum read length to include for analysis [10000000]" )

		parser.add_option( "--minQV", type="float", help="If base has QV < minQV, do not include [0]" )

		parser.add_option( "--min_IPD", type="float", help="If motif IPD is < min_IPD, set to zero [1.0]" )

		parser.add_option( "--min_pct", type="float", help="Remove motifs if they are significant in < min_pct of all reads [5.0]" )

		parser.add_option( "--min_kmer", type="int", help="Minimum motif size to scan (contiguous motifs) [4]" )

		parser.add_option( "--max_kmer", type="int", help="Maximum motif size to scan (contiguous motifs) [6]" )

		parser.add_option( "--mod_bases", type="str", help="String containing bases to query for mods ['AC']" )

		parser.add_option( "--comp_kmer", type="int", help="Kmer size to use for sequence composition measurements [5]" )

		parser.add_option( "--minAcc", type="float", help="Min subread accuracy of read [0.8]" )

		parser.add_option( "--minMapQV", type="float", help="Min mapping QV of aligned read [240]" )

		parser.add_option( "--minReadScore", type="float", help="Min read score of an unaligned read [0.0]" )

		parser.add_option( "--maxPausiness", type="float", help="Max pausiness value of an unaligned read [1000]" )

		parser.add_option( "--bipartite", action="store_true", help="Search bipartite motifs also [False]" )

		parser.add_option( "--comp_only", action="store_true", help="Only run composition binning (no eBinning) [False]" )

		parser.add_option( "--minContigLength", type="int", help="Min length of contig to consider [10000]" )

		parser.add_option( "--minContigForMotifs", type="int", help="Min length of contig to use when getting top motifs [50000]" )

		parser.add_option( "--get_cmph5_stats", action="store_true", help="Generate read stats and exit [False]" )

		parser.add_option( "--skip_motifs", type="str", help="File containing specific motifs to SKIP [None]" )

		parser.add_option( "--motifs_file", type="str", help="File containing specific motifs to include [None]" )

		parser.add_option( "--tmp", type="str", help="Directory where numerous temporary files will be written [tmp]" )

		parser.add_option( "--procs", type="int", help="Number of cores to use [8]" )

		parser.add_option( "--N_reads", type="int", help="Number of qualifying reads to include in analysis [1000000000]" )

		parser.add_option( "--h5_files", type="str", help="Path to file listing the *.h5 files to analyze and their labels [None]" )

		parser.add_option( "--h5_type", type="str", help="cmp (aligned reads) or bax (unaligned) [cmp]" )
		
		parser.add_option( "--contigs", type="str", help="Fasta file containing entries for the assembled contigs [None]" )
		
		parser.add_option( "--bas_whitelist", type="str", help="File containing reads to use in the read-level binning [None]" )
		
		parser.add_option( "--subcontig_size", type="int", help="Size to decompose contigs into for tSNE weighting [50000]" )
		
		parser.add_option( "--min_motif_N", type="int", help="Min number of IPDs from a read/contig to keep for motif filtering [20]" )
		
		parser.add_option( "--min_motif_reads", type="int", help="Min number of reads with motif hits to keep for motif filtering [20]" )
		
		parser.add_option( "--N_motif_reads", type="int", help="Number of reads to scan for motif discovery procedure (will use N_reads \
																if smaller than N_motif_reads) [20000]" )
		
		parser.add_option( "--minMotifIPD", type="float", help="Min motif contig IPD for inclusion of motif in final set [1.7]" )
		
		parser.add_option( "--printIPDs", action="store_true", help="Print out IPD values for SCp ROC curve [False]" )
		
		parser.add_option( "--subtract_control", type="str", help="Subtract control IPDs in final calculations [True]" )
		
		parser.add_option( "--sam", type="str", help="Path to SAM file with CCS alignments for assigning read-level IPDs to contigs [None]" )
		
		parser.add_option( "--control_dir", type="str", help="Path to location to write control IPD data from WGA sequencing [None]" )
		
		parser.add_option( "--cross_cov_bins", type="str", help="Path to file containing binning results from CONCOCT. Will use to \
																 improve motif discovery. Only works with contig-level analysis \
																 (cmp.h5 input) inputs. File format should be '<contig_name>,<bin_id>' \
																 [None]" )
		
		parser.add_option( "--motif_discov_only", action="store_true", help="Quit pipeline after motif discovery [False]" )
		parser.set_defaults( logFile="log.out",                  \
							 info=False,                         \
							 debug=False,                        \
							 synth_mix=False,                    \
							 real_mix=False,                     \
							 printIPDs=False,                    \
							 subtract_control=True,              \
							 subreadlength_min=100,              \
							 readlength_min=100,                 \
							 readlength_max=10000000,            \
							 minQV=0.0,                          \
							 min_IPD=1.0,                        \
							 min_pct=5.0,                        \
							 min_kmer=4,                         \
							 max_kmer=6,                         \
							 comp_kmer=5,                        \
							 mod_bases="A",                      \
							 minAcc=0.8,                         \
							 minMapQV=240,                       \
							 minReadScore=0.0,                   \
							 maxPausiness=1000,                  \
							 minContigLength=10000,              \
							 minContigForMotifs=25000,           \
							 minMotifIPD=1.7,                    \
							 min_motif_N=20,                     \
							 min_motif_reads=20,                 \
							 N_motif_reads=20000,                \
							 get_cmph5_stats=False,              \
							 bipartite=False,                    \
							 comp_only=False,                    \
							 skip_motifs=None,                   \
							 motifs_file=None,                   \
							 tmp="tmp",                          \
							 procs=8,                            \
							 N_reads=1000000000,                 \
							 h5_files=None,                      \
							 h5_type="cmp",                      \
							 bas_whitelist=None,                 \
							 contigs=None,                       \
							 sam=None,                           \
							 control_dir="Not_set",              \
							 cross_cov_bins=None,                \
							 motif_discov_only=False,            \
							 subcontig_size=50000)

		self.opts, args = parser.parse_args( )

		self.opts.bipart_config = [(3,4), (5,6), (3,4)]
		# self.opts.bipart_config = [(3,), (5,), (3,)]

		if len(args) != 0:
			parser.error( "Expected 0 arguments." )		

		if not self.opts.synth_mix and not self.opts.real_mix:
			parser.error( "Specify --real_mix or --synth_mix!" )

		if self.opts.sam!=None and self.opts.motifs_file==None:
			parser.error( "Use of SAM file for read<-->contig mapping only supported with input motifs file using --motifs_file!" )

		if self.opts.control_dir=="Not_set":
			parser.error("Please specify where to build control IPD data (or where an existing build exists)!")

		if self.opts.h5_type=="bas" and self.opts.cross_cov_bins!=None:
			parser.error("Use of the --cross_cov_bins option is not compatible with bas.h5 inputs!")

		self.opts.control_dir = os.path.abspath(self.opts.control_dir)

	def __initLog( self ):
		"""Sets up logging based on command line arguments. Allows for three levels of logging:
		logging.error( ): always emitted
		logging.info( ) : emitted with --info or --debug
		logging.debug( ): only with --debug"""

		if os.path.exists(self.opts.logFile):
			os.remove(self.opts.logFile)

		logLevel = logging.DEBUG if self.opts.debug \
					else logging.INFO if self.opts.info \
					else logging.ERROR

		self.logger = logging.getLogger("")
		self.logger.setLevel(logLevel)
		
		# create file handler which logs even debug messages
		fh = logging.FileHandler(self.opts.logFile)
		fh.setLevel(logLevel)
		
		# create console handler with a higher log level
		ch = logging.StreamHandler()
		ch.setLevel(logLevel)
		
		# create formatter and add it to the handlers
		logFormat = "%(asctime)s [%(levelname)s] %(message)s"
		formatter = logging.Formatter(logFormat)
		ch.setFormatter(formatter)
		fh.setFormatter(formatter)
		
		# add the handlers to logger
		self.logger.addHandler(ch)
		self.logger.addHandler(fh)

if __name__=="__main__":
	app     = eBinner()
	results = app.run()
