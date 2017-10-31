import os,sys
from Bio import SeqIO,Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import glob

fasta   = sys.argv[1]

#########################################
# USAGE: python search_for_circ_contigs.py <CONTIGS.FASTA>
#
# This writes out a set of blasr commands
# specific to the supplied contigs that
# can be executed separately (might
# take awhile to run).
#
#
# Any blasr output files containing 
# alignments (i.e. non-zero file size)
# are possibly circularized and should
# be further investigated.
#########################################

# First preprocess the sequences so that we can align two portions of each sequence to each other
for orig_Seq in SeqIO.parse(fasta, "fasta"):
	overlap = 1000
	overlap = min(overlap, len(orig_Seq.seq))
	
	seq_beg = orig_Seq.seq[:overlap]
	record  = SeqRecord(seq=seq_beg, id=orig_Seq.id+".beg", description=orig_Seq.id+".beg")
	beg_fn  = orig_Seq.id.split("|")[0]+".beg.fa"
	SeqIO.write(record, beg_fn, "fasta")

	seq_rest = orig_Seq.seq[overlap:]
	record   = SeqRecord(seq=seq_rest, id=orig_Seq.id+".rest", description=orig_Seq.id+".rest")
	rest_fn   = orig_Seq.id.split("|")[0]+".rest.fa"
	SeqIO.write(record, rest_fn, "fasta")

	blasr_out = "%s.out" % orig_Seq.id.split("|")[0]
	blasr_CMD = "blasr %s %s -bestn 1 -m 4 -minPctIdentity 98.0 > %s" % (beg_fn, rest_fn, blasr_out)
	
	# Print out each individual blasr command to be executed
	print blasr_CMD