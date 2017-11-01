import numpy as np
import logging
from collections import defaultdict
from Bio import motifs as BioMotifs
from itertools import product
import re

def motifs_from_file( opts ):
	"""
	Read through the file containing motifs to query. Should 
	be in the form:
	
	ACGT-0
	GATC-1
	CATG-1

	where the number indicates the 0-based index of the 
	methylated base in the string.
	"""
	motifs    = set()
	bi_motifs = set()
	for i,line in enumerate(open(opts.motifs_file).xreadlines()):
		line = line.strip()
		motifs.add(line)
	logging.info("Added %s motifs from %s" % (i+1, opts.motifs_file))
	return motifs, bi_motifs

def build_motif_dict( opts ):
	"""
	Generate a set of all possible motifs within provided parameters.
	"""
	motifs      = set()
	bi_motifs   = set()
	NUCS        = "ACGT"
	total_kmers = 0
	logging.info("Initiating dictionary of all possible motifs...")
	#########################################
	# Generate all possible contiguous motifs
	#########################################
	for kmer in range(opts.min_kmer,opts.max_kmer+1):
		total_kmers += len(NUCS)**kmer
		logging.info("  - Adding %s %s-mer motifs..." % (len(NUCS)**kmer, kmer))
		for seq in product("ATGC",repeat=kmer):
			string = "".join(seq)
			for base in opts.mod_bases:
				indexes = [m.start() for m in re.finditer(base, string)]
				for index in indexes:
					motif = "%s-%s" % (string, index)							
					motifs.add(motif)
		logging.info("Done: %s possible contiguous motifs\n" % len(motifs))

	if opts.bipartite:
		####################################################################################
		# Generate all possible bipartite motifs. The parameters specifying 
		# the acceptable forms of bipartite motifs are currently hard-coded as:

		# opts.bipart_config = [(3,4), (5,6), (3,4)]
		
		# where the first item in tuple is the possible lengths of the first part 
		# of the motif, the second item contains the possible number of Ns, and 
		# the third item contains the possible lengths of the second part of the motif.
		# Adding to this motif space greatly increases compute time and memory requirements.
		####################################################################################
		logging.info("  - Adding bipartite motifs to search space...")
		firsts  = opts.bipart_config[0]
		Ns      = opts.bipart_config[1]
		seconds = opts.bipart_config[2]
		for bases1 in firsts:
			for num_Ns in Ns:
				for bases2 in seconds:
					last_mod_pos = bases1-1
					for seq1 in product("ATGC",repeat=bases1):
						for seq2 in product("ATGC",repeat=bases2):
							string = "".join(seq1) + ("N"*num_Ns) + "".join(seq2)
							for base in opts.mod_bases:
								indexes = [m.start() for m in re.finditer(base, string) if m.start()<=last_mod_pos]
								for index in indexes:
									motif = "%s-%s" % (string, index)							
									bi_motifs.add(motif)
		logging.info("Done: %s possible bipartite motifs\n" % len(bi_motifs))
	return motifs, bi_motifs

def add_degen_motifs( motifs, orig_control_means ):
	"""
	If a predetermined set of motifs is input using --motifs_file option,
	create a new entry for the degen motif in the control values dictionary
	by combining the existing data from the various specified versions 
	of the motif.
	"""
	keys_str          = "\n".join(orig_control_means.keys())
	new_control_means = orig_control_means
	for m in motifs:
		new_m = sub_bases(m)
		if new_m!=m:
			matches    = re.findall(new_m, keys_str)
			degen_mean = np.mean([orig_control_means[match] for match in matches])
			new_control_means[m] = degen_mean
			logging.info("Adding degenerate motif %s to controls: %s" % (m, degen_mean))

	return new_control_means

def sub_bases( motif ):
	"""
	Return all possible specifications of a motif with degenerate bases.
	"""
	subs = {"W":"[AT]",  \
			"S":"[CG]",  \
			"M":"[AC]",  \
			"K":"[GT]",  \
			"R":"[AG]",  \
			"Y":"[CT]",  \
			"B":"[CGT]", \
			"D":"[AGT]", \
			"H":"[ACT]", \
			"V":"[ACG]", \
			"N":"[ACGTN]"}
	for symbol,sub in subs.iteritems():
		if motif.find(symbol) > -1:
			motif = motif.replace(symbol, sub)
	return motif

def comp_motif( motif ):
	"""
	Return the complement of the input motif.
	"""
	COMP = {"A":"T", \
			"T":"A", \
			"C":"G", \
			"G":"C", \
			"W":"S", \
			"S":"W", \
			"M":"K", \
			"K":"M", \
			"R":"Y", \
			"Y":"R", \
			"B":"V", \
			"D":"H", \
			"H":"D", \
			"V":"B", \
			"N":"N", \
			"X":"X", \
			"*":"*"}
	r_motif = []
	for char in motif:
		r_motif.append( COMP[char] )
	return "".join(r_motif)

def rev_comp_motif( motif ):
	"""
	Return the reverse complement of the input motif.
	"""
	COMP = {"A":"T", \
			"T":"A", \
			"C":"G", \
			"G":"C", \
			"W":"S", \
			"S":"W", \
			"M":"K", \
			"K":"M", \
			"R":"Y", \
			"Y":"R", \
			"B":"V", \
			"D":"H", \
			"H":"D", \
			"V":"B", \
			"N":"N", \
			"X":"X", \
			"*":"*"}
	rc_motif = []
	for char in motif[::-1]:
		rc_motif.append( COMP[char] )
	return "".join(rc_motif)

def shorten_motifs( contig_motifs, highscore_motifs ):
	"""
	Keep only the shortest, most concise version of the high scoring
	motifs (reduces redundancy).
	"""
	keeper_motifs    = set(highscore_motifs.keys())
	if len(highscore_motifs)>0:
		shortest_contiguous = min([len(m.split("-")[0]) for m in highscore_motifs.keys()])
		# (1) Sort by keys; shortest motif to longest
		motifs_s = sorted(highscore_motifs, key=len)
		# (2) For each motif, check if it's contained in a longer version of other motifs
		for m in motifs_s:
			motif_str =     m.split("-")[0]
			motif_idx = int(m.split("-")[1])
			for remaining in list(keeper_motifs):
				remaining_str =     remaining.split("-")[0]
				remaining_idx = int(remaining.split("-")[1])
				match         = re.search(motif_str, remaining_str)
				if match != None and (motif_idx + match.start()) == remaining_idx and len(remaining_str) > len(motif_str):
					# 3. If True, remove the longer version
					keeper_motifs.remove(remaining)
	return keeper_motifs

def wagner_fischer(word_1, word_2):
	"""

	"""
	n = len(word_1) + 1  # counting empty string 
	m = len(word_2) + 1  # counting empty string
 
	# initialize D matrix
	D = np.zeros(shape=(n, m), dtype=np.int)
	D[:,0] = range(n)
	D[0,:] = range(m)
 
	# B is the backtrack matrix. At each index, it contains a triple
	# of booleans, used as flags. if B(i,j) = (1, 1, 0) for example,
	# the distance computed in D(i,j) came from a deletion or a
	# substitution. This is used to compute backtracking later.
	B = np.zeros(shape=(n, m), dtype=[("del", 'b'), 
									  ("sub", 'b'),
									  ("ins", 'b')])
	B[1:,0] = (1, 0, 0) 
	B[0,1:] = (0, 0, 1)
 
	for i, l_1 in enumerate(word_1, start=1):
		for j, l_2 in enumerate(word_2, start=1):
			deletion = D[i-1,j] + 1
			insertion = D[i, j-1] + 1
			substitution = D[i-1,j-1] + (0 if l_1==l_2 else 2)
 
			mo = np.min([deletion, insertion, substitution])
 
			B[i,j] = (deletion==mo, substitution==mo, insertion==mo)
			D[i,j] = mo
	return D, B

def naive_backtrace(B_matrix):
	"""

	"""
	i, j = B_matrix.shape[0]-1, B_matrix.shape[1]-1
	backtrace_idxs = [(i, j)]
 
	while (i, j) != (0, 0):
		if B_matrix[i,j][1]:
			i, j = i-1, j-1
		elif B_matrix[i,j][0]:
			i, j = i-1, j
		elif B_matrix[i,j][2]:
			i, j = i, j-1
		backtrace_idxs.append((i,j))
 
	return backtrace_idxs

def align(word_1, word_2, bt):
	"""

	"""
	aligned_word_1 = []
	aligned_word_2 = []
	operations     = []
 
	backtrace = bt[::-1]  # make it a forward trace
 
	for k in range(len(backtrace) - 1): 
		i_0, j_0 = backtrace[k]
		i_1, j_1 = backtrace[k+1]
 
		w_1_letter = None
		w_2_letter = None
		op = None
 
		if i_1 > i_0 and j_1 > j_0:  # either substitution or no-op
			if word_1[i_0] == word_2[j_0]:  # no-op, same symbol
				w_1_letter = word_1[i_0]
				w_2_letter = word_2[j_0]
				op = " "
			else:  # cost increased: substitution
				w_1_letter = word_1[i_0]
				w_2_letter = word_2[j_0]
				op = "s"
		elif i_0 == i_1:  # insertion
			w_1_letter = " "
			w_2_letter = word_2[j_0]
			op = "i"
		else: #  j_0 == j_1,  deletion
			w_1_letter = word_1[i_0]
			w_2_letter = " "
			op = "d"
 
		aligned_word_1.append(w_1_letter)
		aligned_word_2.append(w_2_letter)
		operations.append(op)
 
	return aligned_word_1, aligned_word_2, operations

def refine_degen_motifs( keeper_motifs, contig_motifs, case_motif_Ns ):
	"""
	Identify redundant instances of degenerate motifs and replace them with the
	most parsimonious representation. 

	'A': 'A', 
	'C': 'C', 
	'G': 'G', 
	'T': 'T', 
	'AC': 'M', 
	'AG': 'R', 
	'AT': 'W', 
	'CG': 'S', 
	'CT': 'Y', 
	'GT': 'K', 
	'ACG': 'V', 
	'ACT': 'H', 
	'AGT': 'D', 
	'CGT': 'B', 
	'ACGT': 'N'

	In order to call a motif consensus to identify degenerate bases, we first
	need to identify those motifs that are likely different specifications
	of the same motif containing degenerate bases. To build this set of related 
	motifs, we filter based on:

	(1) motif length
	(2) methylated position in the motif
	(3) number of central Ns in bipartite motifs
	(4) edit distance between the most representative of the motifs and all
	    others in the set. Most representative determined by all-vs-all edit
	    distance calculations.

	Once this set of related motifs is selected, we call a consensus to get 
	degenerate bases.
	"""
	def edit_distance( word1, word2 ):
		D,B             = wagner_fischer(word1, word2)
		bt              = naive_backtrace(B)
		alignment_table = align(word1, word2, bt)
		n_edits = len([entry for entry in alignment_table[2] if entry!=" "])
		return n_edits

	refined_motifs = []
	degen_members  = defaultdict(list)
	
	# First sort by (1) motif lengths
	lens = map(lambda x: len(x), list(keeper_motifs))
	for size in list(set(lens)):
		size_motifs = [m for m in list(keeper_motifs) if len(m)==size]
		
		# Next, sort by (2) methylated potision in motif
		idxs = set(map(lambda x: x.split("-")[1], size_motifs))
		for idx in list(idxs):
			idx_motifs = [m for m in size_motifs if m.split("-")[1]==idx]
			
			# Count number of Ns in bipartite motifs
			n_N_motifs = defaultdict(list)
			for m in idx_motifs:
				n_N = len([nuc for nuc in m if nuc=="N"])
				n_N_motifs[n_N].append(m)
		
			# Now sort based on (3) number of central Ns in bipartite motifs
			for n_N,motifs in n_N_motifs.iteritems():
				motif_set  = set()
				leftovers  = set()
				
				# Finally, calculate edit distance between remaining motifs.
				# Run all-against-all edit distance to determine the most 
				# representative of the motifs in the set.
				edit_mat = np.zeros([len(motifs), len(motifs)])
				for i in range(len(motifs)):
					for j in range(len(motifs)):
						edit_mat[i,j] = edit_distance(motifs[i], motifs[j])
				
				# Find the motif with the least total edits
				min_idx = np.argmin([np.sum(edit_mat[i,:]) for i in range(edit_mat.shape[0])])
				word1   = motifs[min_idx]

				# Calculate edit distance against all other motifs
				other_motifs = [motif for motif in motifs if motif!=word1]
				for word2 in other_motifs:
					n_edits = edit_distance(word1, word2)
					n_Ns    = len([x for x in word1 if x=="N"])
					if (n_Ns==0 and n_edits<=1) or (n_Ns>0 and n_edits<=2):
						# Close enough edit distance to use in consensus calling
						motif_set.add(word1)
						motif_set.add(word2)
					else:
						# Not close enough edit distance; do not group this motif
						# with the others; will retain and keep separate
						leftovers.add(word2)
				if len(motif_set)==0:
					# No companion motif sets found for consensus. Cannot search for
					# degenerate bases.
					refined_motifs+=motifs
				else:
					# Successfully found companion motifs. Gather information about these
					# motifs and call consensus motif to identify degenerate bases.
					SCp_values    = [contig_motifs[m] for m in list(motif_set)]
					N_values      = [case_motif_Ns[m] for m in list(motif_set)]
					for_consensus = [m.split("-")[0] for m in list(motif_set)]
					# Must treat contiguous and bipartite motifs slightly differently
					if n_N>0:
						# BioMotifs.degenerate_consensus() cannot accept Ns;
						# will temporarily replace Ns with Ts.
						replaced = []
						for j,m in enumerate(for_consensus):
							mk = np.array([i for i,nuc in enumerate(m) if nuc=="N"])
							replaced.append( m.replace("N","T") )
						m = BioMotifs.create(replaced)
						x = list(m.degenerate_consensus)
						for pos in mk:
							x[pos] = "N"
						degen_motif = "".join(x) + "-%s" % idx
					else:
						# No need to replace Ns with Ts
						m = BioMotifs.create(for_consensus)
						degen_motif = str(m.degenerate_consensus) + "-%s" % idx
					
					logging.info("Refining %s motifs to %s" % (len(motif_set), degen_motif))
					for motif in list(motif_set):
						logging.info("   - %s" % motif)
					refined_motifs.append(degen_motif)

					# Calculate SCp based on the mean of all specified motifs it encompasses
					contig_motifs[degen_motif] = np.mean(SCp_values)
					case_motif_Ns[degen_motif] = np.sum(N_values)
					
					# Keep record of which original motifs are now represented by the
					# new degenerate motif.
					for o_motif in list(motif_set):
						degen_members[degen_motif].append(o_motif)

				refined_motifs += list(leftovers)

	return set(refined_motifs), contig_motifs, case_motif_Ns, degen_members