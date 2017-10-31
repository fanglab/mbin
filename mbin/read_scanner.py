import re
from collections import Counter
from itertools import product
import numpy as np
import logging
import motif_tools

# def find_motif_matches( mode, motif, ref_str, read_str, strand ):
def find_motif_matches( mode, motif, ref_str, strand ):
	if mode == "cmp":
		if strand == 0:
			q_motif = motif_tools.sub_bases( motif_tools.rev_comp_motif(motif) )
		elif strand == 1:
			q_motif = motif_tools.sub_bases( motif_tools.comp_motif(motif) )
		matches_iter = re.finditer(q_motif, ref_str)
	elif mode == "bas":
		q_motif = motif_tools.sub_bases( motif_tools.rev_comp_motif(motif) )
		matches_iter = re.finditer(q_motif, ref_str)

	matches_list = []
	for match in matches_iter:
		matches_list.append(match)
	return matches_list

def kmer_freq ( mode, ref_str, strand, opts ):
	ref_str = ref_str.upper()
	if strand==1:
		ref_str = ref_str[::-1]
	k = opts.comp_kmer
	kmers = []
	for seq in product("ATGC",repeat=k):
		kmers.append( "".join(seq) )

	kmer_counts = Counter()
	for j in range( len(ref_str)-(k-1) ):
		motif    = ref_str[j:j+k]
		kmer_counts[motif] += 1

	# Combine forward and reverse complement motifs into one count
	combined_kmer = Counter()
	for kmer in kmers:
		kmer_rc = motif_tools.rev_comp_motif(kmer)
		if not combined_kmer.get(kmer_rc):
			combined_kmer[kmer] = kmer_counts[kmer] + kmer_counts[kmer_rc] + 1

	return combined_kmer

def scan_motifs( mode, ipds, ref_str, strand, motifs, bi_motifs, opts ):
	"""
	For each motif, find all occurrences in read.
	"""
	read_str     = None
	if mode == "bas":
		strand = 0
	subread_ipds = motif_ipds( mode,      \
							   # read_str,  \
							   ipds,      \
							   ref_str,   \
							   strand,    \
							   motifs,    \
							   bi_motifs, \
							   opts )
	
	def compute_score( ipds ):
		mean = np.mean(ipds)
		return mean

	subread_barcode = {}
	for motif,ipds in subread_ipds.iteritems():
		if len(ipds) == 0:
			score = 0.0
		else:
			score = compute_score( ipds )
		subread_barcode[motif] = ( len(ipds), score )

	subread_kmers = kmer_freq( mode,     \
							   # read_str, \
							   ref_str, \
							   strand,  \
							   opts )
	return subread_barcode, subread_kmers

def walk_over_read( mode, subread_ipds, ref_str, read_ipds, strand, k, opts):
	"""
	Loop over each position in the read string, adding motifs as they
	are encountered in the walk.
	"""
	for j in range( len(ref_str)-3 ):
		seq  = ref_str[j:j+k]
		ipds = read_ipds[j:j+k]

		if seq.find("*") == -1 and seq.find("X") == -1:
			if mode == "cmp":
				if strand == 0:
					q_motif  = motif_tools.rev_comp_motif( seq )
				elif strand == 1:
					q_motif  = motif_tools.comp_motif( seq )
			elif mode == "bas":
				q_motif = motif_tools.rev_comp_motif( ref_str[j:j+k] )
			
			for base in opts.mod_bases:
				ref_indexes = [m.start() for m in re.finditer(base, q_motif)]
				for ref_index in ref_indexes:
					rc_index = len(q_motif) - 1 - ref_index
					
					if mode == "cmp":
						if strand == 0:
							idx = rc_index
						elif strand == 1:
							idx = ref_index
					elif mode == "bas":
						idx = rc_index
					IPD = ipds[idx]
					
					# If a contiguous motif contains an N, skip 
					# if q_motif.find("N")>-1 and not opts.bipartite:
					# 	continue

					ref_motif_str = "%s-%s" % (q_motif, ref_index)
					try:
						subread_ipds[ref_motif_str].append( IPD )
					except KeyError:
						# logging.warning("Motif %s has unexpected characters (N,etc). Skipping..." % ref_motif_str)
						pass
	return subread_ipds

def walk_over_read_bipartite( subread_ipds, ref_str, read_ipds, strand, opts ):
	"""
	Loop over each position in the read string, adding motifs as they
	are encountered in the walk.
	"""
	max_motif_len = 15
	firsts        = opts.bipart_config[0]
	Ns            = opts.bipart_config[1]
	seconds       = opts.bipart_config[2]
	for j in range( len(ref_str)-max_motif_len ):
		for first in firsts:
			last_mod_pos = first-1
			for N in Ns:
				for second in seconds:
					length = first + N + second
					seq    = ref_str[j:j+length]
					ipds   = read_ipds[j:j+length]
					if seq.find("*") == -1 and seq.find("X") == -1:
						if strand == 0:
							q_motif  = motif_tools.rev_comp_motif( seq )
						elif strand == 1:
							q_motif  = motif_tools.comp_motif( seq )

						for base in opts.mod_bases:
							ref_indexes = [m.start() for m in re.finditer(base, q_motif) if m.start() <= last_mod_pos]
							for ref_index in ref_indexes:
								rc_index      = len(q_motif) - 1 - ref_index
								if strand == 0:
									idx = rc_index
								elif strand == 1:
									idx = ref_index
								IPD      = ipds[idx]
								bi_motif = "".join( [q_motif[:first],"N"*N,q_motif[-second:]] )
								ref_motif_str = "%s-%s" % (bi_motif, ref_index)
								
								try:
									subread_ipds[ref_motif_str].append( IPD )
								except KeyError:
									# logging.warning("Motif %s has unexpected characters (N,etc). Skipping..." % ref_motif_str)
									pass
	return subread_ipds

def motif_ipds( mode, read_ipds, ref_str, strand, motifs, bi_motifs, opts ):
	read_id      = np.random.randint(1000000)
	subread_ipds = {}
	if opts.motifs_file != None:
		for motif in motifs:
			if motif.find("-") == -1:
				raise Exception("Specify the position of the modified base in your supplied motifs (0-based)")
			ref_index = int(motif.split("-")[1])
			motif     =     motif.split("-")[0]
			rc_index  = len(motif) - 1 - ref_index

			matches_list = find_motif_matches( mode, motif, ref_str, strand )

			motif_ipds = []
			for match in matches_list:
				motif_start = match.span()[0]
				motif_end   = match.span()[1]
				motif_ipds.append( read_ipds[motif_start:motif_end] )

			if mode == "cmp" and strand == 1:
				ipds = map(lambda x: x[ref_index], motif_ipds)
			else:
				ipds = map(lambda x: x[rc_index], motif_ipds)
			ref_motif_str = "%s-%s" % (motif, ref_index)

			subread_ipds[ref_motif_str] = ipds
	else:
		# Instatiate a motif IPD dict where the default value is 0
		for motif in motifs:
			subread_ipds[motif] = []

		for bi_motif in bi_motifs:
			subread_ipds[bi_motif] = []

		if opts.skip_motifs != None:
			# We want to first remove these bases in the read and ref (and associated IPDs)
			motifs_to_del = []
			for line in open(opts.skip_motifs).readlines():
				motif = line.strip()
				motifs_to_del.append(motif)
				motifs_to_del = map(lambda x: x.split("-")[0], motifs_to_del)
				try:
					del subread_ipds[motif]
				except KeyError:
					pass
			for motif in motifs_to_del:
				# matches_list = find_motif_matches( mode, motif, ref_str, read_str, strand )
				matches_list = find_motif_matches( mode, motif, ref_str, strand )
				motif_ipds = []
				for match in matches_list:
					# Remove the data associated with this motif instance
					motif_start = match.span()[0]
					motif_end   = match.span()[1]
					ref_str     = "".join([ref_str[:motif_start], "X"*len(motif), ref_str[motif_end:]])
					read_ipds   = read_ipds[:motif_start] + ([100.0]*len(motif)) + read_ipds[motif_end:]
					if mode == "cmp":
						ref_str     = "".join([ref_str[:motif_start], "X"*len(motif), ref_str[motif_end:]])

		# Loop over all possible contiguous motif sizes
		for k in range( 4, opts.max_kmer+1 ):
			# subread_ipds = walk_over_read( mode, subread_ipds, read_str, ref_str, read_ipds, strand, k, opts)
			subread_ipds = walk_over_read( mode, subread_ipds, ref_str, read_ipds, strand, k, opts)
		
		# Now scan for valid bipartite motifs
		if opts.bipartite:
			for bi_motif in bi_motifs:
				subread_ipds[bi_motif] = []
			# subread_ipds = walk_over_read_bipartite( subread_ipds, read_str, read_ipds, strand, opts )
			subread_ipds = walk_over_read_bipartite( subread_ipds, ref_str, read_ipds, strand, opts )

	return subread_ipds