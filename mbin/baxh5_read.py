import os,sys,glob
import re
import math
import multiprocessing
from pbcore.io.BasH5IO import BasH5Reader
import subprocess
import numpy as np
from collections import OrderedDict,Counter,defaultdict
import logging
import read_scanner
from operator import itemgetter
import pickle

import warnings
warnings.simplefilter("error")

class subread_motif_processor:
	def __init__( self, baxh5, chunk_id, idx, N_target_reads, motifs, bi_motifs, opts ):
		self.baxh5          = baxh5
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
				self.call    = tup[0]
				self.ipd     = tup[1]

		class subread:
			def __init__( self, subname, label, base_calls, IPD, QV, opts ):
				self.subname = subname
				self.label   = label
				self.entries = {}
				self.opts    = opts

				# Pull out read-level information
				base_calls = np.array(list(base_calls))
				QVs        = np.array(QV)
				fps        = 75.0
				ipds       = np.array(IPD) / fps

				if self.opts.minQV != 0:
					# Don't consider basecalls & IPDs with bad QVs
					for i,val in enumerate(QV):
						if val < self.opts.minQV:
							base_calls[i] = "*"
							ipds[i]       = 99999

				# Attach these IPD entries to the molecule object
				for read_pos, tup in enumerate(zip(base_calls, ipds)):
					self.entries[ read_pos ] = ipd_entry(tup)

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
					if entry.ipd != 99999:
						subread_vals.append(entry.ipd)

				nfs = np.mean(map(lambda x: math.log(x + 0.001), subread_vals))

				for read_pos, entry in self.entries.iteritems():
					if entry.ipd != 99999:
						newIPD = math.log(entry.ipd + 0.001) - nfs
						entry.ipd = newIPD

			def zip_bases_and_IPDs( self ):
				"""
				Reassemble the read and IPD values using the subread normalized IPDs
				"""
				od        = OrderedDict(sorted(self.entries.items()))
				read      = []
				self.ipds = []
				for read_pos, entry in od.items():
					read.append(entry.call)
					self.ipds.append(entry.ipd)
				self.read_str = "".join(read)

		def dump_data_to_proc_files( to_dump, read_labs ):
			subname_fn     = "subreads_names.tmp"
			label_fn       = "subreads_labels.tmp"   
			length_fn      = "subreads_lengths.tmp"  
			ipds_fn        = "subreads_ipds.tmp"     
			ipds_N_fn      = "subreads_ipdsN.tmp"   
			comp_N_fn      = "subreads_compN.tmp"

			self.tmp_fns.add( os.path.join(self.chunkdir, subname_fn) )
			self.tmp_fns.add( os.path.join(self.chunkdir, label_fn) )
			self.tmp_fns.add( os.path.join(self.chunkdir, length_fn) )
			self.tmp_fns.add( os.path.join(self.chunkdir, ipds_fn) )
			self.tmp_fns.add( os.path.join(self.chunkdir, ipds_N_fn) )
			self.tmp_fns.add( os.path.join(self.chunkdir, comp_N_fn) )
			f_subnames   = open(os.path.join(self.opts.tmp, self.chunkdir, subname_fn),  "a")
			f_labels     = open(os.path.join(self.opts.tmp, self.chunkdir, label_fn),    "a")
			f_lengths    = open(os.path.join(self.opts.tmp, self.chunkdir, length_fn),   "a")
			f_ipds       = open(os.path.join(self.opts.tmp, self.chunkdir, ipds_fn),     "a")
			f_ipds_N     = open(os.path.join(self.opts.tmp, self.chunkdir, ipds_N_fn),   "a")
			f_comp_N     = open(os.path.join(self.opts.tmp, self.chunkdir, comp_N_fn),   "a")
			self.tmp_fs.add(f_subnames)
			self.tmp_fs.add(f_labels)
			self.tmp_fs.add(f_ipds)
			self.tmp_fs.add(f_ipds_N)
			self.tmp_fs.add(f_comp_N)
			
			if self.opts.motifs_file!=None and self.opts.subtract_control:
				control_ipds_d = pickle.load( open(self.opts.control_pkl_name,"rb" ) )

			for i,(subread_ipds,subread_comps,readname,subread_length) in enumerate(to_dump):
				for motif in subread_ipds.iterkeys():
					if subread_ipds[motif]==[]:
						subread_ipds[motif] = [0.0]
				ipd_kmers   = [motif                                            for motif in subread_ipds.iterkeys()]
				ipd_means   = [subread_ipds[motif][1]                           for motif in subread_ipds.iterkeys()]
				ipd_counts  = [subread_ipds[motif][0]                           for motif in subread_ipds.iterkeys()]

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
				if i==0:
					ipds_kmers_fn = "subreads_ipdskmers.tmp"
					comp_kmers_fn = "subreads_compkmers.tmp"
					self.tmp_fns.add( os.path.join(self.chunkdir, ipds_kmers_fn) )
					self.tmp_fns.add( os.path.join(self.chunkdir, comp_kmers_fn) )
					f_ipds_kmers      = open(os.path.join(self.opts.tmp, self.chunkdir, ipds_kmers_fn), "w")
					f_comp_kmers      = open(os.path.join(self.opts.tmp, self.chunkdir, comp_kmers_fn), "w")
					ipds_kmers_str    = "\t".join(ipd_kmers)
					comp_kmers_str    = "\t".join(comp_kmers)
					f_ipds_kmers.write("%s\n" % ipds_kmers_str)
					f_comp_kmers.write("%s\n" % comp_kmers_str)
					f_ipds_kmers.close()
					f_comp_kmers.close()
				ipds_str        = "\t".join(map(lambda x: str(round(x,3)), ipd_means))
				ipds_N_str      = "\t".join(map(lambda x: str(x),          ipd_counts))
				comp_counts_str = "\t".join(map(lambda x: str(x),          comp_counts))
				f_subnames.write(  "%s\n" % readname)
				f_labels.write(    "%s\n" % read_labs[readname])
				f_lengths.write(   "%s\n" % subread_length)
				f_ipds.write(      "%s\n" % ipds_str)
				f_ipds_N.write(    "%s\n" % ipds_N_str)
				f_comp_N.write(    "%s\n" % comp_counts_str)

			for f in self.tmp_fs:
				f.close()

		reader = BasH5Reader(self.baxh5)

		cwd = os.getcwd()

		# Periodically (after <chunksize> subreads) write out data to a process-specific tmp file
		chunksize     = 10
		self.chunkdir = "chunk_%s" % self.chunk_id
		if os.path.exists(os.path.join(self.opts.tmp, self.chunkdir)):
			shutil.rmtree(os.path.join(self.opts.tmp, self.chunkdir))
		os.mkdir(os.path.join(self.opts.tmp, self.chunkdir))
		to_dump = []

		read_SMp      = {}
		read_SMp_N    = {}
		read_comps    = {}
		read_labs     = {}
		read_lengths  = Counter()
		SMp_kmer_N    = {}

		self.tmp_fs  = set()
		self.tmp_fns = set()

		# if self.opts.bas_whitelist != None and not self.opts.control_run:
		# 	bas_whitelist = set(np.loadtxt(self.opts.bas_whitelist, dtype="str"))
		# 	zmws_to_check = [zmw for zmw in reader[self.idx] if zmw.zmwName in bas_whitelist]
		# else:
		zmws_to_check = reader[self.idx]

		i      = 0
		j      = 0
		to_get = min(self.N_target_reads, len(self.idx))
		incr   = to_get/4

		for y,zmw in enumerate(zmws_to_check):
			has_good_subs = False
			label         = self.opts.h5_labels[self.baxh5]
			subreads = [sub for sub in zmw.subreads]
			# read_length = 0
			# for sub in subreads:
			# 	read_length += len(sub.basecalls())
			
			if i%incr==0:
				logging.info("...chunk %s - ZMW %s/%s (%.1f%%)" % (self.chunk_id, i, to_get, 100*float(i)/to_get))
			
			for w,sub in enumerate(subreads):
				subname             = sub.readName
				subread_length      = len(sub.basecalls())
				hole                = sub.holeNumber
				base_calls          = sub.basecalls()
				IPD                 = sub.IPD()
				QV                  = sub.QualityValue()
				readname            = "/".join(subname.split("/")[:-1])

				if subread_length >= self.opts.subreadlength_min:
					has_good_subs = True
					
					

					sub = subread( subname, label, base_calls, IPD, QV, self.opts )
					sub.zip_bases_and_IPDs()
					subread_ipds,subread_comps = read_scanner.scan_motifs( "bas",          \
																		   sub.ipds,       \
																		   sub.read_str,   \
																		   "None",         \
																		   self.motifs,    \
																		   self.bi_motifs, \
																		   self.opts )
					
					read_labs[readname] = label
					to_dump.append( (subread_ipds,subread_comps,readname,len(sub.read_str)) )
					# Dump subread IPD and comp data to process-specific file
					if len(to_dump)%chunksize==0 and len(to_dump)!=0:
						dump_data_to_proc_files( to_dump, read_labs )
						to_dump = []

					# Combine subread-level IPD info to create read-level eBarcodes
					if not read_SMp.get(readname):
						read_SMp[readname]   = defaultdict(list)
						SMp_kmer_N[readname] = Counter()
					for motif,(N,score) in subread_ipds.iteritems():
						read_vals                    = [score]*N
						read_SMp[readname][motif]   += read_vals
						SMp_kmer_N[readname][motif] += N

					# Combine subread kmer frequency info to get read-level kmer frequencies
					if not read_comps.get(readname):
						read_comps[readname] = Counter()
					for kmer,N in subread_comps.items():
						read_comps[readname][kmer] += N
					j += 1
			if has_good_subs:
				i += 1

			if i==self.N_target_reads:
				break
		# else:
		# 	pass

		# print "N reads:", i, "N subreads:", j
		# print "================="

		dump_data_to_proc_files( to_dump, read_labs )

		for f in self.tmp_fs:
			f.close()
		to_dump = []

		if i==0:
			logging.info("Chunk %s: no qualifying reads found!" % self.chunk_id)

		logging.info("Chunk %s: found %s subreads (%s molecules) > %sbp in %s" % (self.chunk_id, j, len(read_labs.keys()), self.opts.readlength_min, os.path.basename(self.baxh5)))
		reader.close()

		return self.tmp_fns