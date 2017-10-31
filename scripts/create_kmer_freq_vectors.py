import os,sys
from itertools import product
from collections import Counter,OrderedDict
from Bio import SeqIO
import numpy as np
import math

######################## USAGE ##########################
#
# python create_kmer_freq_vectors.py <FASTA> <k-mer size>
#
#########################################################

fasta_fn =     sys.argv[1]
k        = int(sys.argv[2])

def kmer_freq ( ref_str, k ):
	"""
	Walk through sequence and return k-mer counts plus
	a pseudocount of 1.
	"""
	ref_str = ref_str.upper()
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
		kmer_rc = rev_comp_motif(kmer)
		if not combined_kmer.get(kmer_rc):
			combined_kmer[kmer] = kmer_counts[kmer] + kmer_counts[kmer_rc] + 1

	return combined_kmer

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

for i,entry in enumerate(SeqIO.parse(fasta_fn, "fasta")):
	# Return k-mer counts for each entry in the fasta
	kmers     = []
	seq_comp  = []
	seq_kmers = kmer_freq( str(entry.seq), k )

	seq_kmers_sort = OrderedDict(sorted(seq_kmers.items()))
	for kmer,count in seq_kmers_sort.iteritems():

		# Normalize each count by total k-mers and log2 transform
		kmer_normed_comp = math.log(float(count) / sum(seq_kmers_sort.values()))
		# kmer_normed_comp = count
		kmers.append(kmer)
		seq_comp.append(kmer_normed_comp)

	kmer_vector = np.array(kmers)
	comp_vector = np.array(seq_comp)

	if i==0:
		# print header
		print "seq_name\t%s" % "\t".join(kmer_vector)

	print "%s\t%s" % (entry.id, "\t".join(map(lambda x: str(round(x,4)), comp_vector)))