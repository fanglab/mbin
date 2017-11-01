import os,sys,glob
import re
import shutil
import math
import multiprocessing
from pbcore.io.align.CmpH5IO import CmpH5Reader
import subprocess
import numpy as np
from collections import OrderedDict,Counter,defaultdict
import logging
import mbin
import read_scanner
from operator import itemgetter
import pickle
import unicodedata

import warnings
warnings.simplefilter("error")

class subread_motif_processor:
	def __init__( self, cmph5, chunk_id, idx, N_target_reads, motifs, bi_motifs, opts ):
		self.cmph5          = cmph5
		self.chunk_id       = chunk_id
		self.idx            = idx
		self.N_target_reads = N_target_reads
		self.motifs         = motifs
		self.bi_motifs      = bi_motifs
		self.opts           = opts

	def __call__( self ):
		
		class ipd_entry:
			def __init__( self, tup ):
				"""
				"""
				self.ref_base  = tup[0]
				self.ipd       = tup[1]
				# self.call      = tup[2]
				# self.read_base = tup[3]
				self.ref_pos   = tup[2]

		class subread:
			def __init__( self, cmph5, alignment, label, opts ):
				leftAnchor   = 1
				rightAnchor  = 1
				self.entries = {}
				self.opts    = opts

				self.subname    = alignment.readName
				movieID         = alignment.movieInfo[0]
				alignedLength   = alignment.referenceSpan
				fps             = alignment.movieInfo[2]
				self.refName    = alignment.referenceInfo[3]
				zmw             = alignment.HoleNumber
				self.mol        = alignment.MoleculeID
				if alignment.isForwardStrand:
					self.strand = 0
				else:
					self.strand = 1
				self.ref_bases  = alignment.reference()
				# self.read_bases = alignment.read()
				
				read_calls      = alignment.transcript()
				ref_pos         = list(alignment.referencePositions())
				IPD             = list(alignment.IPD())
				self.label      = self.opts.h5_labels[cmph5]

				error_mk = []
				for read_call in read_calls:
					# Go through all entries and flag which positions are MM/indels
					if read_call != "M":
						# Mismatch or indel at this position!
						error_mk.append(1)
					else:
						error_mk.append(0)
				
				# Get the indices of all the non-matches
				error_idx = [i for (i,val) in enumerate(error_mk) if val == 1]
				for error_id in error_idx:
					try:
						for j in range(leftAnchor):
							error_mk[error_id - (j+1)] = 1
						for j in range(rightAnchor):
							error_mk[error_id + (j+1)] = 1
					except IndexError:
						pass
				error_mk = np.array(error_mk)

				ipds     = np.array(IPD) / fps
				
				strands  = np.array([self.strand] * len(read_calls))

				self.ref_bases  = np.array(list(self.ref_bases))
				# self.read_bases = np.array(list(self.read_bases))
				self.ref_pos    = np.array(ref_pos)
				read_calls      = np.array(list(read_calls))

				# Mark the error positions, but leave them in the sequence so
				# we can pull out intact motifs from contiguous correct bases
				self.ref_bases[error_mk==1]  = "*"
				# self.read_bases[error_mk==1] = "*"
				read_calls[error_mk==1] = "*"
				ipds[error_mk==1]       = -9
				strands[error_mk==1]    = -9

				# Attach these IPD entries to the subread object
				# for i,tup in enumerate(zip(self.ref_bases, ipds, read_calls, self.read_bases, self.ref_pos)):
				for i,tup in enumerate(zip(self.ref_bases, ipds, self.ref_pos)):
					entry = ipd_entry(tup)
					self.entries[ self.ref_pos[i] ] = ipd_entry(tup)

				# self.cap_outliers()

				self.subread_normalize()

			def cap_outliers( self, max_ipd=10 ):
				"""
				Cap the outlier IPDs at max_ipd seconds.
				"""
				for read_pos,entry in self.entries.iteritems():
					entry.ipd = min(entry.ipd, max_ipd)

			def subread_normalize( self ):
				"""
				Every IPD entry needs to be normalized by the mean IPD of its subread.
				"""
				if len(self.entries) == 0:
					# Nothing to do here.
					return self.entries

				# First populate list of all IPDs per subread. Will use to get normalization factor.
				subread_vals = []
				for entry in self.entries.values():
					# Only do if this IPD is NOT from an error position
					if entry.ipd != -9:
						subread_vals.append(entry.ipd)

				rawIPDs = np.array(map(lambda x: math.log(x + 0.001), subread_vals))
				nfs     = rawIPDs.mean()

				for pos, entry in self.entries.iteritems():
					if entry.ipd == -9:
						newIPD = -9
					else:
						newIPD = math.log(entry.ipd + 0.001) - nfs
					
					entry.ipd = newIPD

			def zip_bases_and_IPDs( self ):
				"""
				Reassemble the read and IPD values using the subread normalized IPDs
				"""
				od        = OrderedDict(sorted(self.entries.items()))
				ref       = []
				ref_pos   = []
				self.ipds = []
				for read_pos, entry in od.items():
					ref.append(entry.ref_base)
					ref_pos.append(entry.ref_pos)
					self.ipds.append(entry.ipd)
				self.ref_str  = "".join(ref)
				self.ref_pos  = ref_pos

		reader = CmpH5Reader(self.cmph5)

		read_refs     = {}
		read_SMp      = {}
		read_SMp_N    = {}
		read_comps    = {}
		read_labs     = {}
		contig_SCp    = {}
		i             = 0
		n_mols        = 0

		cwd = os.getcwd()

		# Periodically (after <chunksize> alignments) write out data to a contig-specific tmp file
		chunksize     = 10
		self.chunkdir = "chunk_%s" % self.chunk_id
		if os.path.exists(os.path.join(self.opts.tmp, self.chunkdir)):
			shutil.rmtree(os.path.join(self.opts.tmp, self.chunkdir))
		os.mkdir(os.path.join(self.opts.tmp, self.chunkdir))
		to_dump = defaultdict(list)

		def dump_data_to_contig_files( refName, to_dump, read_labs ):
			refName           = mbin.slugify(refName)
			ref_subname_fn    = "%s_readnames.tmp"    % refName
			ref_label_fn      = "%s_labels.tmp"       % refName
			ref_length_fn     = "%s_lengths.tmp"      % refName
			ref_ipds_fn       = "%s_ipds.tmp"         % refName
			ref_ipds_N_fn     = "%s_ipdsN.tmp"        % refName
			ref_comp_N_fn     = "%s_compN.tmp"        % refName
			ref_strand_fn     = "%s_strand.tmp"       % refName

			self.tmp_fns.add( os.path.join(self.chunkdir, ref_subname_fn) )
			self.tmp_fns.add( os.path.join(self.chunkdir, ref_label_fn) )
			self.tmp_fns.add( os.path.join(self.chunkdir, ref_length_fn) )
			self.tmp_fns.add( os.path.join(self.chunkdir, ref_ipds_fn) )
			self.tmp_fns.add( os.path.join(self.chunkdir, ref_ipds_N_fn) )
			self.tmp_fns.add( os.path.join(self.chunkdir, ref_comp_N_fn) )
			self.tmp_fns.add( os.path.join(self.chunkdir, ref_strand_fn) )
			f_subnames   = open(os.path.join(self.opts.tmp, self.chunkdir, ref_subname_fn),  "a")
			f_labels     = open(os.path.join(self.opts.tmp, self.chunkdir, ref_label_fn),    "a")
			f_lengths    = open(os.path.join(self.opts.tmp, self.chunkdir, ref_length_fn),   "a")
			f_ipds       = open(os.path.join(self.opts.tmp, self.chunkdir, ref_ipds_fn),     "a")
			f_ipds_N     = open(os.path.join(self.opts.tmp, self.chunkdir, ref_ipds_N_fn),   "a")
			f_comp_N     = open(os.path.join(self.opts.tmp, self.chunkdir, ref_comp_N_fn),   "a")
			f_strand     = open(os.path.join(self.opts.tmp, self.chunkdir, ref_strand_fn),   "a")
			self.tmp_fs.add(f_subnames)
			self.tmp_fs.add(f_labels)
			self.tmp_fs.add(f_ipds)
			self.tmp_fs.add(f_ipds_N)
			self.tmp_fs.add(f_comp_N)
			self.tmp_fs.add(f_strand)
			
			if self.opts.motifs_file!=None and self.opts.subtract_control:
				control_ipds_d = pickle.load( open(self.opts.control_pkl_name,"rb" ) )

			for i,(subread_ipds,subread_comps,readname,subread_length,strand) in enumerate(to_dump[refName]):
				ipd_kmers   = [motif                  for motif in subread_ipds.iterkeys()]
				ipd_means   = [subread_ipds[motif][1] for motif in subread_ipds.iterkeys()]
				ipd_counts  = [subread_ipds[motif][0] for motif in subread_ipds.iterkeys()]

				ipd_means = []
				if self.opts.motifs_file!=None and self.opts.subtract_control:
					for motif in subread_ipds.iterkeys():
						if subread_ipds[motif][1] != 0.0:
							w_control_sub = subread_ipds[motif][1] - control_ipds_d[motif]
							ipd_means.append(w_control_sub)
						else: # Don't subtract control if no ipd values are available (i.e. IPD score == 0.0)
							ipd_means.append(subread_ipds[motif][1])
				else:
					for motif in subread_ipds.iterkeys():
						ipd_means.append(subread_ipds[motif][1])

				comp_kmers  = np.array( [motif   for motif,ipds in subread_comps.items()] )
				comp_counts = np.array( [ipds    for motif,ipds in subread_comps.items()] )
				if i==0 and refName not in self.refName_has_header:
					ref_ipds_kmers_fn = "%s_ipdskmers.tmp" % refName
					ref_comp_kmers_fn = "%s_compkmers.tmp" % refName
					f_ipds_kmers      = open(os.path.join(self.opts.tmp, self.chunkdir, ref_ipds_kmers_fn), "a")
					f_comp_kmers      = open(os.path.join(self.opts.tmp, self.chunkdir, ref_comp_kmers_fn), "a")
					ipds_kmers_str    = "\t".join(ipd_kmers)
					comp_kmers_str    = "\t".join(comp_kmers)
					f_ipds_kmers.write("%s\n" % ipds_kmers_str)
					f_comp_kmers.write("%s\n" % comp_kmers_str)
					f_ipds_kmers.close()
					f_comp_kmers.close()
					self.refName_has_header.add(refName)
				ipds_str        = "\t".join(map(lambda x: str(round(x,3)), ipd_means))
				ipds_N_str      = "\t".join(map(lambda x: str(x),          ipd_counts))
				comp_counts_str = "\t".join(map(lambda x: str(x),          comp_counts))
				f_subnames.write(  "%s\n" % readname)
				f_labels.write(    "%s\n" % read_labs[readname])
				f_lengths.write(   "%s\n" % subread_length)
				f_ipds.write(      "%s\n" % ipds_str)
				f_ipds_N.write(    "%s\n" % ipds_N_str)
				f_comp_N.write(    "%s\n" % comp_counts_str)
				f_strand.write(    "%s\n" % strand)

			for f in self.tmp_fs:
				f.close()
		self.tmp_fs             = set()
		self.tmp_fns            = set()
		self.refName_has_header = set()
		to_check = reader[self.idx]
		for alignment in to_check:
			ref_contig = mbin.slugify(alignment.referenceInfo[3])
			label      = self.opts.h5_labels[self.cmph5]
			ref_len    = self.opts.cmph5_contig_lens[self.cmph5][ref_contig]
			if ref_len >= self.opts.minContigLength and alignment.referenceSpan >= self.opts.readlength_min and alignment.MapQV >= self.opts.minMapQV:
				to_get    = min(self.N_target_reads, len(self.idx))
				incr      = to_get/10
				readname  = "/".join(alignment.readName.split("/")[:-1])
				if len(read_labs.keys())%incr==0 and not read_labs.get(readname):
					logging.info("...chunk %s\t- mol %s/%s (%.1f%%)" % (self.chunk_id, n_mols, to_get, 100*n_mols/to_get))

				read_labs[readname] = label
				read_refs[readname] = ref_contig

				sub = subread( self.cmph5, alignment, label, self.opts )
				sub.zip_bases_and_IPDs()
				subread_ipds,subread_comps = read_scanner.scan_motifs( "cmp",          \
																	   # sub.read_str,   \
																	   sub.ipds,       \
																	   sub.ref_str,    \
																	   sub.strand,     \
																	   self.motifs,    \
																	   self.bi_motifs, \
																	   self.opts )
				
				to_dump[ref_contig].append( (subread_ipds,subread_comps,readname,len(sub.ref_str),sub.strand) )
				# Dump subread IPD and comp data to contig-specific file
				if len(to_dump[ref_contig])%chunksize==0 and len(to_dump[ref_contig])!=0:
					dump_data_to_contig_files( ref_contig, to_dump, read_labs )
					to_dump[ref_contig] = []

				n_mols = len(read_labs.keys())
				i += 1
					
				if n_mols==self.N_target_reads:
					break

		for ref_contig in to_dump.keys():
			dump_data_to_contig_files( ref_contig, to_dump, read_labs )
		for f in self.tmp_fs:
			f.close()
		to_dump = defaultdict(list)

		if i==0:
			logging.info("Chunk %s: no qualifying reads found!" % self.chunk_id)
		
		logging.info("Chunk %s: found %s alignments (%s molecules) > %sbp in %s" % (self.chunk_id, i, len(read_labs.keys()), self.opts.readlength_min, os.path.basename(self.cmph5)))
		reader.close()

		return self.tmp_fns