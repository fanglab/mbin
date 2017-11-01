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
import logging
import multiprocessing
import operator
import subprocess
import glob
from itertools import izip
import pickle
import pysam
import motif_tools
import unicodedata

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

def cat_list_of_files( in_fns, out_fn, del_ins=True ):
	"""
	Given a list of filenames, cat them together into a single file. Then cleans up pre-catted
	single files.
	"""
	if len(in_fns)==0:
		raise Exception("There are no files to cat!")
		
	cat_CMD   = "cat %s > %s" % (" ".join(in_fns), out_fn)
	p         = subprocess.Popen(cat_CMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdOutErr = p.communicate()
	sts       = p.returncode
	if sts != 0:
		raise Exception("Failed cat command: %s" % cat_CMD)
	if del_ins:
		for fn in in_fns:
			try:
				os.remove(fn)
			except OSError:
				pass
	return out_fn

def slugify(text):
	"""Generates an ASCII-only slug."""
	slug = unicodedata.normalize("NFKD",unicode(text)).encode("ascii", "ignore")
	slug = re.sub(r"[^\w]+", " ", slug)
	slug = "-".join(slug.lower().strip().split())
	return slug

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
			cat_list_of_files(to_cat, outname, del_ins=False)

def cat_subreads_files_from_subprocs( subproc_tmp_files, movie_i, opts ):
	for ftype in subproc_tmp_files.keys():
		if ftype in ["compkmers.tmp", "ipdskmers.tmp"]:
			first_fn = glob.glob( os.path.join(opts.tmp, "chunk_*", "subreads_%s" % ftype) )[0]
			shutil.copy( first_fn, os.path.join(opts.tmp, "subreads_%s" % ftype.replace(".tmp", ".%s.tmp" % movie_i)))
		else:
			ftype_fns = list(subproc_tmp_files[ftype])
			sorted(ftype_fns)
			catted_name = os.path.join(opts.tmp, "subreads_%s" % ftype.replace(".tmp", ".%s.tmp" % movie_i))
			cat_list_of_files( ftype_fns, catted_name, del_ins=False )

class mbinRunner:
	def __init__( self, opts ):
		"""
		The options from __main__.py are passed to mbinRunner
		"""
		self.opts = opts
		self.fns  = {}
	
		##########################################################
		# Define output filenames
		##########################################################
		if self.opts.h5_type=="cmp":
			#############################################
			# Define the alignment-level output filenames
			#############################################
			self.fns["read_names"]           = "reads.names"
			self.fns["read_refs"]            = "reads.refs"
			self.fns["read_labels"]          = "reads.labels"
			self.fns["read_lengths"]         = "reads.lengths"
			self.fns["read_strands"]         = "reads.strands"
			self.fns["read_comp_kmers"]      = "reads.comp_kmers"
			self.fns["read_comp_counts"]     = "reads.comp_counts"
			self.fns["read_comp"]            = "reads.comp"
			self.fns["read_comp_2D"]         = "reads.comp.2D"
			self.fns["read_comp_2D_z"]       = "reads.comp.2D.zscores"
			self.fns["read_SMp_kmers"]       = "reads.SMp_kmers"
			self.fns["read_SMp_counts"]      = "reads.SMp_counts"
			self.fns["read_SMp"]             = "reads.SMp"
			self.fns["read_SMp_2D"]          = "reads.SMp.2D"
			self.fns["read_SMp_2D_z"]        = "reads.SMp.2D.zscores"
			self.fns["read_4D_z"]            = "reads.combined.4D.zscores"
			self.fns["read_combo_2D"]        = "reads.combined.2D"

			# #############################################
			# # Define the contig-level output filenames
			# #############################################
			self.fns["contig_names"]         = "contigs.names"
			self.fns["contig_labels"]        = "contigs.labels"
			self.fns["contig_lengths"]       = "contigs.lengths"
			self.fns["contig_comp"]          = "contigs.comp"
			self.fns["contig_comp_2D"]       = "contigs.comp.2D"
			self.fns["contig_comp_2D_z"]     = "contigs.comp.2D.zscores"
			self.fns["contig_cov_comp_3D"]   = "contigs.cov_comp.3D"
			self.fns["contig_cov_comp_2D"]   = "contigs.cov_comp.2D"
			self.fns["contig_cov_comp_2D_z"] = "contigs.cov_comp.2D.zscores"
			self.fns["contig_SCp"]           = "contigs.SCp"
			self.fns["contig_SCp_N"]         = "contigs.SCp_N"
			self.fns["contig_SCp_2D"]        = "contigs.SCp.2D"
			self.fns["contig_SCp_2D_z"]      = "contigs.SCp.2D.zscores"
			self.fns["contig_cov"]           = "contigs.cov"
			self.fns["contig_4D_z"]          = "contigs.combined.4D.zscores"
			self.fns["contig_combo_2D"]      = "contigs.combined.2D"
		
			if self.opts.cross_cov_bins!=None:
				self.fns["bin_names"]        = "bins.names"
				self.fns["bin_sizes"]        = "bins.sizes"
				self.fns["bin_SCp"]          = "bins.SCp"
				self.fns["bin_SCp_N"]        = "bins.SCp_N"

		elif self.opts.h5_type=="bas":
			#############################################
			# Define the unaligned read-level output filenames
			#############################################
			self.fns["read_names"]           = "reads.names"
			self.fns["read_labels"]          = "reads.labels"
			self.fns["read_lengths"]         = "reads.lengths"
			self.fns["read_comp_kmers"]      = "reads.comp_kmers"
			self.fns["read_comp_counts"]     = "reads.comp_counts"
			self.fns["read_comp"]            = "reads.comp"
			self.fns["read_comp_2D"]         = "reads.comp.2D"
			self.fns["read_comp_2D_z"]       = "reads.comp.2D.zscores"
			self.fns["read_SMp_kmers"]       = "reads.SMp_kmers"
			self.fns["read_SMp_counts"]      = "reads.SMp_counts"
			self.fns["read_SMp"]             = "reads.SMp"
			self.fns["read_SMp_2D"]          = "reads.SMp.2D"
			self.fns["read_SMp_2D_z"]        = "reads.SMp.2D.zscores"
			self.fns["read_4D_z"]            = "reads.combined.4D.zscores"
			self.fns["read_combo_2D"]        = "reads.combined.2D"
			
			if self.opts.sam!=None:
				####################################################
				# Define contig filenames when SAM file used to map reads --> contigs
				####################################################
				self.fns["read_refs"]            = "reads.ccs.refs"
				self.fns["contig_names"]         = "contigs.names"
				self.fns["contig_labels"]        = "contigs.labels"
				self.fns["contig_lengths"]       = "contigs.lengths"
				self.fns["contig_comp"]          = "contigs.comp"
				self.fns["contig_comp_2D"]       = "contigs.comp.2D"
				self.fns["contig_comp_2D_z"]     = "contigs.comp.2D.zscores"
				self.fns["contig_cov_comp_3D"]   = "contigs.cov_comp.3D"
				self.fns["contig_cov_comp_2D"]   = "contigs.cov_comp.2D"
				self.fns["contig_cov_comp_2D_z"] = "contigs.cov_comp.2D.zscores"
				self.fns["contig_SCp"]           = "contigs.SCp"
				self.fns["contig_SCp_N"]         = "contigs.SCp_N"
				self.fns["contig_SCp_2D"]        = "contigs.SCp.2D"
				self.fns["contig_SCp_2D_z"]      = "contigs.SCp.2D.zscores"
				self.fns["contig_cov"]           = "contigs.cov"
				self.fns["contig_4D_z"]          = "contigs.combined.4D.zscores"
				self.fns["contig_combo_2D"]      = "contigs.combined.2D"

	def chunks( self, l, n ):
		"""
		Yield successive n-sized chunks from l.
		"""
		for i in xrange(0, len(l), n):
			yield l[i:i+n]

	def launch_data_loader( self, h5_file, N_reads, movie_i, opts ):
		logging.info("Loading data from %s..." % h5_file)
		results = self.launch_subprocs( h5_file, N_reads, opts )
		if opts.h5_type=="cmp":
			contig_tmps = defaultdict(list)
			for p,result in enumerate(results):
				for fn in result:
					contig = "_".join(fn.split("/")[1].split("_")[:-1])
					contig_tmps[contig].append(fn)
			logging.info("Catting contig-specific data from parallel processes...")
			args    = [(contig, opts.tmp, i, len(contig_tmps.keys())) for i,contig in enumerate(contig_tmps.keys())]
			results = launch_pool( opts.procs, cat_contig_files_from_subprocs, args )
		elif opts.h5_type=="bas":
			logging.info("Catting subprocess-specific data from parallel processes...")
			subproc_tmp_files = defaultdict(list)
			for result in results:
				for fn in result:
					ftype = fn.split(os.sep)[-1].split("subreads_")[-1]
					subproc_tmp_files[ftype].append(os.path.join(opts.tmp, fn))
			
			results = cat_subreads_files_from_subprocs( subproc_tmp_files, movie_i, opts )			
			logging.info("Done.")	

		for chunkdir in glob.glob( os.path.join(opts.tmp, "chunk_*")):
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
		

		cat_list_of_files(labels_fns,    os.path.join(self.opts.tmp, self.fns["read_labels"]), del_ins=False)
		cat_list_of_files(lengths_fns,   os.path.join(self.opts.tmp, self.fns["read_lengths"]), del_ins=False)
		cat_list_of_files(refs_fns,      os.path.join(self.opts.tmp, self.fns["read_refs"]), del_ins=False)
		cat_list_of_files(readnames_fns, os.path.join(self.opts.tmp, self.fns["read_names"]), del_ins=False)
		cat_list_of_files(strands_fns,   os.path.join(self.opts.tmp, self.fns["read_strands"]), del_ins=False)
		cat_list_of_files(ipds_fns,      os.path.join(self.opts.tmp, self.fns["read_SMp"]), del_ins=False)
		cat_list_of_files(ipds_N_fns,    os.path.join(self.opts.tmp, self.fns["read_SMp_counts"]), del_ins=False)
		cat_list_of_files(comp_N_fns,    os.path.join(self.opts.tmp, self.fns["read_comp"]), del_ins=False)
		shutil.copy(comp_kmers_fns[0],              os.path.join(self.opts.tmp, self.fns["read_comp_kmers"]))
		shutil.copy(ipds_kmers_fns[0],              os.path.join(self.opts.tmp, self.fns["read_SMp_kmers"]))

		for fn in comp_kmers_fns:
			os.remove(fn)

		for fn in ipds_kmers_fns:
			os.remove(fn)

	def combine_subreads_for_contig_level( self, cmph5_file ):
		readnames         = np.loadtxt(os.path.join(self.opts.tmp, self.fns["read_names"]),       dtype="str")
		read_refs         = np.loadtxt(os.path.join(self.opts.tmp, self.fns["read_refs"]),        dtype="str")
		reads_SMp         = np.loadtxt(os.path.join(self.opts.tmp, self.fns["read_SMp"]),         dtype="float")
		reads_SMp_counts  = np.loadtxt(os.path.join(self.opts.tmp, self.fns["read_SMp_counts"]),  dtype="int")
		reads_lengths     = np.loadtxt(os.path.join(self.opts.tmp, self.fns["read_lengths"]),     dtype="int")
		
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
				contig_comp_str               = "\t".join(map(lambda x: str(round(x,4)), contig_comp))
				if name.find("|quiver")>-1:
					name = name.replace("|quiver", "") 
				refName                       = slugify(name)
				contigs_comp_strings[refName] = contig_comp_str
				contig_fasta_lens[refName]    = len(contig_seq)
			contigs_comp_strings["unknown"]   = "\t".join(map(lambda x: str(round(x,4)), [0.000] * len(contig_comp)))
			return contigs_comp_strings,contig_fasta_lens

		# Calculate composition profiles for all contigs in the contigs.fasta file
		contigs_comp_strings,self.contig_fasta_lens = calc_contig_comps( self.opts.contigs )

		f_ref_names       = open(os.path.join(self.opts.tmp, self.fns["contig_names"]),   "w")
		f_ref_SCp         = open(os.path.join(self.opts.tmp, self.fns["contig_SCp"]),     "w")
		f_ref_SCp_N       = open(os.path.join(self.opts.tmp, self.fns["contig_SCp_N"]),   "w")
		f_ref_comp        = open(os.path.join(self.opts.tmp, self.fns["contig_comp"]),    "w")
		f_ref_labels      = open(os.path.join(self.opts.tmp, self.fns["contig_labels"]),  "w")
		f_ref_covs        = open(os.path.join(self.opts.tmp, self.fns["contig_cov"]),     "w")
		f_ref_lens        = open(os.path.join(self.opts.tmp, self.fns["contig_lengths"]), "w")
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
					# quiver_ref = slugify(ref+"|quiver")
					quiver_ref = slugify(ref)
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
				ref_SCp_str   = "\t".join(map(lambda x: str(round(x,4)), ref_SCp[ref]))
				ref_SCp_N_str = "\t".join(map(lambda x: str(x), ref_SCp_N[ref]))
				f_ref_SCp.write(   "%s\n" % ref_SCp_str)
				f_ref_SCp_N.write( "%s\n" % ref_SCp_N_str)
				
				# Sometimes the contig name in cmp.h5 trims off the "|quiver" in "<unitig_N>|quiver"
				try:
					f_ref_comp.write(  "%s\n" % contigs_comp_strings[ref])
				except KeyError:
					# quiver_ref = slugify(ref+"|quiver")
					quiver_ref = slugify(ref)
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

		contigs_names   = np.loadtxt(os.path.join(self.opts.tmp, self.fns["contig_names"]),   dtype="str")
		contigs_SCp     = np.loadtxt(os.path.join(self.opts.tmp, self.fns["contig_SCp"]),     dtype="float")
		contigs_SCp_N   = np.loadtxt(os.path.join(self.opts.tmp, self.fns["contig_SCp_N"]),   dtype="int")
		contigs_lengths = np.loadtxt(os.path.join(self.opts.tmp, self.fns["contig_lengths"]), dtype="int")
		
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
				bins_SBp_str   = "\t".join(map(lambda x: str(round(x,4)), bins_SBp[bin_id]))
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
		
		cat_list_of_files(labels_fns,    os.path.join(self.opts.tmp, "subreads_labels.tmp"))
		cat_list_of_files(lengths_fns,   os.path.join(self.opts.tmp, "subreads_lengths.tmp"))
		cat_list_of_files(readnames_fns, os.path.join(self.opts.tmp, "subreads_names.tmp"))
		cat_list_of_files(ipds_fns,      os.path.join(self.opts.tmp, "subreads_ipds.tmp"))
		cat_list_of_files(ipds_N_fns,    os.path.join(self.opts.tmp, "subreads_ipdsN.tmp"))
		cat_list_of_files(comp_N_fns,    os.path.join(self.opts.tmp, "subreads_compN.tmp"))
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
			reads_names_fn      = os.path.join(self.opts.tmp, self.fns["read_names"])
			reads_labels_fn     = os.path.join(self.opts.tmp, self.fns["read_labels"])
			reads_lengths_fn    = os.path.join(self.opts.tmp, self.fns["read_lengths"])
			reads_ipds_fn       = os.path.join(self.opts.tmp, self.fns["read_SMp"])
			reads_ipds_N_fn     = os.path.join(self.opts.tmp, self.fns["read_SMp_counts"])
			reads_ipds_kmers_fn = os.path.join(self.opts.tmp, self.fns["read_SMp_kmers"])
			reads_comp_N_fn     = os.path.join(self.opts.tmp, self.fns["read_comp"])
			reads_comp_kmers_fn = os.path.join(self.opts.tmp, self.fns["read_comp_kmers"])

		f_reads   = open(reads_names_fn,   "w")
		f_labels  = open(reads_labels_fn,  "w")
		f_lengths = open(reads_lengths_fn, "w")
		f_ipds    = open(reads_ipds_fn,    "w")
		f_ipds_N  = open(reads_ipds_N_fn,  "w")
		f_comp_N  = open(reads_comp_N_fn,  "w")
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
			f_ipds.write(    "%s\n" % "\t".join(map(lambda x: str(round(x,4)), ipds_list)))
			f_ipds_N.write(  "%s\n" % "\t".join(map(lambda x: str(x),          ipds_N_list)))
			f_comp_N.write(  "%s\n" % "\t".join(map(lambda x: str(round(x,4)), normed_comp_N_list)))
		
		shutil.copy(ipds_kmers_fn, reads_ipds_kmers_fn)
		shutil.copy(comp_kmers_fn, reads_comp_kmers_fn)
		f_reads.close()
		f_labels.close()
		f_lengths.close()
		f_ipds.close()
		f_ipds_N.close()
		f_comp_N.close()

		# Remove the subread-level barcodes for each contig
		for fn in subreads_fns:
			os.remove(fn)

	def get_read_refs_from_SAM( self ):
		read_refs = {}
		read_lens = {}
		samfile   = pysam.AlignmentFile(self.opts.sam, "rb")
		for k,aln in enumerate(samfile.fetch()):
			if aln.mapq == 254:
				read            = "/".join(aln.query_name.split("/")[:2])
				read            = "_".join(read.split("_")[1:])
				# ref             = aln.reference_name.split("|")[0]
				ref             = slugify(aln.reference_name)
				read_refs[read] = ref
				read_lens[read] = len(aln.seq)

		# Write the read_refs file to match the readnames file
		readnames = np.loadtxt(os.path.join( self.opts.tmp, self.fns["read_names"]), dtype="str")
		f_names   = open(os.path.join( self.opts.tmp, self.fns["read_refs"]),    "wb")
		f_lens    = open(os.path.join( self.opts.tmp, self.fns["read_lengths"]), "wb")
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

	def launch_subprocs( self, h5_file, N_reads, opts ):
		"""
		"""
		logging.debug("Creating tasks...")
		tasks   = multiprocessing.JoinableQueue()
		results = multiprocessing.Queue()
		logging.debug("Done.")

		if opts.h5_type=="cmp":
			reader     = CmpH5Reader(h5_file)
			to_check   = reader
			entries    = range(len(to_check))
		elif opts.h5_type=="bas":
			reader     = BasH5Reader(h5_file)
			if opts.bas_whitelist != None and not opts.control_run:
				logging.info("Intersecting with whitelist...")
				bas_whitelist = set(np.loadtxt(opts.bas_whitelist, dtype="str"))
				to_check = [z for z in reader if z.zmwName in bas_whitelist]
			else:
				to_check = reader
			# Filter on zmw metrics
			pre      = len(to_check)
			logging.info("Starting with %s reads..." % pre)
			to_check = [z for z in to_check if z.zmwMetric("Productivity")==1 and \
											   z.zmwMetric("ReadScore")>opts.minReadScore and \
											   z.zmwMetric("Pausiness")<opts.maxPausiness]
			post     = len(to_check)
			logging.info("Dropped %s reads due to poor zmw metrics (%s remain)" % ((pre-post),post))
			# Filter on read length 
			pre      = post
			to_check = [z for z in to_check if np.sum([len(sub) for sub in z.subreads])>=opts.readlength_min]
			post     = len(to_check)
			logging.info("Dropped %s reads < %s (%s remain)" % ((pre-post),opts.readlength_min, post))
			entries  = np.array([z.holeNumber for z in to_check])

		reader.close()
		if len(entries) <= opts.procs * 5:
			procs = 1
		else:
			procs = opts.procs

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
			if opts.h5_type=="cmp":
				tasks.put(cmph5_read.subread_motif_processor( h5_file, chunk_id, idx, n, opts.motifs, opts.bi_motifs, opts) )
				logging.debug("...%s (%s alignments)" % (chunk_id, len(entries)))
			elif opts.h5_type=="bas":
				tasks.put(baxh5_read.subread_motif_processor( h5_file, chunk_id, idx, n, opts.motifs, opts.bi_motifs, opts) )
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


if __name__=="__main__":
	app     = mbinRunner()
	results = app.run()
