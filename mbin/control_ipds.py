import os,sys
import glob
import numpy as np
import logging
from collections import defaultdict
from itertools import product
import re
from pbcore.io.align.CmpH5IO import CmpH5Reader
from pbcore.io.BasH5IO import BasH5Reader
import cmph5_read
import shutil
import pickle
import multiprocessing
import math
import subprocess
import motif_tools

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

def chunks( l, n ):
	"""
	Yield successive n-sized chunks from l.
	"""
	for i in xrange(0, len(l), n):
		yield l[i:i+n]

def process_contig_chunk( args ):
	chunk_id      = args[0]
	cut_CMDs      = args[1]
	kmers         = args[2]
	cols_chunk    = args[3]
	n_chunks      = args[4]
	min_motif_N   = args[5]
	logging.info("  - Control data: chunk %s/%s" % ((chunk_id+1), (n_chunks+1)))
	control_means = {}
	
	for cut_CMD in cut_CMDs:
		sts,stdOutErr = run_OS_command( cut_CMD )
	
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
		if np.sum(control_ipds_N_sub[:,j])>=min_motif_N:
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
		self.control_h5       = wga_h5
		self.opts             = opts
		self.orig_dir = os.getcwd()

	def build_insilico_controls( self ):
		"""
		To be added...
		"""
		pass

	def check_control_file( self ):
		"""
		If control data already exists, point to it and skip.
		Otherwise, specify where to generate the control data.
		"""
		if (self.opts.use_control!=None and os.path.exists(self.opts.use_control)):
			# Control data already exists
			launch_control = False
		else:
			# Will need to run the read scan to build control data
			launch_control = True
		
		return launch_control

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

	def add_degen_motifs( self, motifs, orig_control_means ):
		"""
		If a predetermined set of motifs is input using --motifs_file option,
		create a new entry for the degen motif in the control values dictionary
		by combining the existing data from the various specified versions 
		of the motif.
		"""
		keys_str          = "\n".join(orig_control_means.keys())
		new_control_means = orig_control_means
		for m in motifs:
			new_m = motif_tools.sub_bases(m)
			if new_m!=m:
				matches    = re.findall(new_m, keys_str)
				degen_mean = np.mean([orig_control_means[match] for match in matches])
				new_control_means[m] = degen_mean
				logging.info("Adding degenerate motif %s to controls: %s" % (m, degen_mean))

		return new_control_means

	def scan_WGA_h5( self ):
		"""
		Get some necessary information about the WGA bas/cmp.h5 
		being used to generate the control IPD data.
		"""
		if self.opts.h5_type=="cmp":
			# contig-level analysis
			self.opts.h5_labels[self.control_h5]         = "control"
			self.opts.cmph5_contig_lens[self.control_h5] = {}
			
			reader = CmpH5Reader(self.control_h5)
			for entry in reader.referenceInfoTable:
				name      = entry[3]
				length    = entry[4]
				slug_name = cmph5_read.slugify(name)
				self.opts.cmph5_contig_lens[self.control_h5][slug_name] = length
			reader.close()
		elif self.opts.h5_type=="bas":
			self.opts.h5_labels[self.control_h5] = "control"

		return self.opts

	def analyze_WGA_reads( self ):
		"""
		Launch read scanning pipeline for building up
		control IPD values for motifs.
		"""
		if self.opts.h5_type=="cmp":
			# Using reads aligned to contigs
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
					cmph5_read.cat_list_of_files(to_cat, outname)
		
		elif self.opts.h5_type=="bas":
			# Using unaligned reads
			results     = self.combine_subread_data_across_bas_movies()
			results     = self.bas_combine_subreads_for_read_level( tmp_run=True )
			control_fns = glob.glob( os.path.join(self.opts.tmp, "*.tmp"))

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
			args.append( (i, cut_CMDs, kmers, cols_chunk, n_chunks, self.opts.min_motif_N) )
		
		results = launch_pool(self.opts.procs, process_contig_chunk, args)
		
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
		
		cmph5_read.cat_list_of_files(labels_fns,    os.path.join(self.opts.tmp, "control_labels.tmp"))
		cmph5_read.cat_list_of_files(strands_fns,   os.path.join(self.opts.tmp, "control_strands.tmp"))
		cmph5_read.cat_list_of_files(lengths_fns,   os.path.join(self.opts.tmp, "control_lengths.tmp"))
		cmph5_read.cat_list_of_files(readnames_fns, os.path.join(self.opts.tmp, "control_names.tmp"))
		cmph5_read.cat_list_of_files(ipds_fns,      os.path.join(self.opts.tmp, "control_ipds.tmp"))
		cmph5_read.cat_list_of_files(ipds_N_fns,    os.path.join(self.opts.tmp, "control_ipdsN.tmp"))
		cmph5_read.cat_list_of_files(comp_N_fns,    os.path.join(self.opts.tmp, "control_compN.tmp"))
		
		shutil.copy(comp_kmers_fns[0],              os.path.join(self.opts.tmp, "control_compkmers.tmp"))
		shutil.copy(ipds_kmers_fns[0],              os.path.join(self.opts.tmp, "control_ipdskmers.tmp"))
		x = [os.remove(fn) for fn in comp_kmers_fns]
		x = [os.remove(fn) for fn in ipds_kmers_fns]

	def build_control_IPD_dict( self, motifs, bi_motifs ):
		"""

		"""
		if self.opts.h5_type=="cmp":
			control_ipds_fn   = glob.glob(os.path.join(self.opts.tmp, "control_ipds.tmp" ))
			control_ipds_N_fn = glob.glob(os.path.join(self.opts.tmp, "control_ipdsN.tmp"))
			control_kmers_fn  = glob.glob(os.path.join(self.opts.tmp, "control_ipdskmers.tmp"))
		elif self.opts.h5_type=="bas":
			control_ipds_fn   = glob.glob(os.path.join(self.opts.tmp, "read_ipds.tmp"))
			control_ipds_N_fn = glob.glob(os.path.join(self.opts.tmp, "read_ipdsN.tmp"))
			control_kmers_fn  = glob.glob(os.path.join(self.opts.tmp, "read_ipdskmers.tmp"))

		if (len(control_ipds_fn)>1 or len(control_ipds_N_fn)>1 or len(control_kmers_fn)>1):
			raise Exception("*** Double check the control files. There should not be multiples for a file type.")

		control_means,not_found = self.chunk_control_matrices(control_ipds_fn[0], control_ipds_N_fn[0], control_kmers_fn[0])

		if not_found > 0:
			logging.warning("WARNING: %s/%s motifs did not have sufficient coverage in control data!" % (not_found, (len(motifs)+len(bi_motifs))))
			logging.warning("   * If this is alarming, try increasing --N_motif_reads to get higher coverage and make sure N_reads > N_motif_reads")
		
		logging.info("Writing control data to a pickled file for future re-use: %s" % self.opts.write_control)
		pickle.dump( control_means, open( self.opts.write_control, "wb" ) )
		logging.info("Done.")

		return control_means

	def load_control_pickle( self ):
		"""

		"""
		control_means = pickle.load(open(self.opts.use_control, "rb"))
		return control_means