import os,sys
import glob
from Bio import SeqIO
import numpy as np
from collections import defaultdict

####################################
# USAGE
# --------
# While in the directory of data downloaded 
# from https://genome.jgi.doe.gov/pages/dynamicOrganismDownload.jsf?organism=IMG_3300005806,
# run this script:
#
# python write_xiao_MGS_bin_fastas.py
#
# This will create a directory called MGS
# and write a single multi-entry fasta file
# (each entry is a binned contig) for each
# of the 541 MGS bins.

mgs2gene    = np.loadtxt("MmCAG2geneID.txt",     dtype="str", skiprows=1)
gene2contig = np.loadtxt("3300005806.a.map.txt", dtype="str")
contigs     = "3300005806.a.fna"

gene2MGSs = defaultdict(list)
for x in mgs2gene:
	gene2MGSs[x[1]].append(x[0])

contig2gene = {}
for x in gene2contig:
	contig2gene[x[1]] = x[0]

contig2MGSs = {}
for contig,gene in contig2gene.iteritems():
	MGSs                = gene2MGSs[gene]
	contig2MGSs[contig] = MGSs

mgs_contigs = defaultdict(list)
for sequence in SeqIO.parse(contigs,"fasta"):
	contig = sequence.id
	
	if contig2MGSs.get(contig):
		MGSs = contig2MGSs[contig]
		for mgs in MGSs:
			if mgs[:5]=="MmMGS":
				mgs_contigs[mgs].append(sequence)
				# print mgs, sequence.id

for mgs,sequences in mgs_contigs.iteritems():
	mgs_id = mgs.split(":")[1]
	mgs_fn = "MGS/MGS_%s.fasta" % mgs_id 
	SeqIO.write(sequences, mgs_fn, "fasta")