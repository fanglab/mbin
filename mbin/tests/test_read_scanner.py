#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `mbin` package."""

import unittest
import os
from mbin import read_scanner,motif_tools
from mbin import data

class TestReadscanner(unittest.TestCase):
	"""Tests for `mbin` package."""

	def setUp(self):
		"""
		The test sequence will contain instances of non-palindromic and 
		palindromic contiguous motifs (ACCACC-0 and TGATCA-2), as well 
		as a non-palindromic bipartite motif (GACNNNNNACC-1).
		
		              _____________________       _____________________   _________________________________________
		Ref seq:  C   A   C   C   A   C   C   T   T   G   A   T   C   A   G   A   C   T   T   T   T   T   A   C   C
		Methyl:       *                                   *                   *
		Read seq: G   T   G   G   T   G   G   A   A   C   T   A   G   T   C   T   G   A   A   A   A   A   T   G   G
		IPDs:    0.4 2.4 0.3 0.1 0.1 0.2 0.4 0.4 0.2 0.1 2.2 0.4 0.3 0.5 0.1 2.3 0.2 0.6 0.5 0.4 0.3 0.4 0.2 0.1 0.3


		"""
		self.MOTIFS  = ["ACCACC-0","TGATCA-2", "GACNNNNNACC-1"]
		self.REF_STR = "CACCACCTTGATCAGACTTTTTACC"
		self.IPD_VALS = [0.4,2.4,0.3,0.1,0.1,0.2,0.4,0.4,0.2,0.1,2.2,0.4,0.3,0.5,0.1,2.3,0.2,0.6,0.5,0.4,0.3,0.4,0.2,0.1,0.3]

	def tearDown(self):
		"""Tear down test fixtures, if any."""

	def test_find_motif_matches_aligned_strand_f(self):
		mode = "aligned"
		motif = self.MOTIFS[0].split("-")[0] # motif = ACCACC --> rc_motif = GGTGGT
		rc_REF_STR = motif_tools.rev_comp_motif(self.REF_STR)
		strand = 0
		matches_list = read_scanner.find_motif_matches(mode,motif,rc_REF_STR,strand)
		match1 = matches_list[0]
		self.assertEqual((18, 24), match1.span()[:2])

	def test_find_motif_matches_aligned_strand_r(self):
		mode = "aligned"
		motif = self.MOTIFS[0].split("-")[0] # motif = ACCACC --> rc_motif = GGTGGT
		c_REF_STR = motif_tools.comp_motif(self.REF_STR)
		strand = 1
		matches_list = read_scanner.find_motif_matches(mode,motif,c_REF_STR,strand)
		match1 = matches_list[0]
		self.assertEqual((1,7), match1.span()[:2])

	def test_find_motif_matches_unaligned(self):
		mode = "unaligned"
		motif = self.MOTIFS[0].split("-")[0] # motif = ACCACC --> rc_motif = GGTGGT
		rc_REF_STR = motif_tools.rev_comp_motif(self.REF_STR)
		strand = 0
		matches_list = read_scanner.find_motif_matches(mode,motif,rc_REF_STR,strand)
		match1 = matches_list[0]
		self.assertEqual((18, 24), match1.span()[:2])

	def test_kmer_freq_f(self):
		strand = 0
		opts = data.set_opts()
		kmer_counts = read_scanner.kmer_freq(self.REF_STR, strand, opts)
		self.assertTrue(kmer_counts["TGATC"]==3 and kmer_counts["TGGTG"]==2 and kmer_counts["AGACT"]==2 and kmer_counts["TTACC"]==2)

	def test_kmer_freq_r(self):
		strand = 1
		opts = data.set_opts()
		kmer_counts = read_scanner.kmer_freq(self.REF_STR, strand, opts)
		self.assertTrue(kmer_counts["ACTAG"]==3 and kmer_counts["AAATG"]==2 and kmer_counts["TGGTG"]==2 and kmer_counts["AGACT"]==2)

	def test_scan_motifs_aligned_f(self):
		mode = "aligned"
		strand = 0
		rc_REF_STR = motif_tools.rev_comp_motif(self.REF_STR)
		r_IPD_VALS = self.IPD_VALS[::-1]
		motifs = set(self.MOTIFS)
		bi_motifs = set()
		opts = data.set_opts()
		barcode, kmers = read_scanner.scan_motifs(mode, r_IPD_VALS, rc_REF_STR, strand, motifs, bi_motifs, opts)
		self.assertTrue(barcode["ACCACC-0"]==(1, round(2.4,1)) and \
						barcode["TGATCA-2"]==(1, round(2.2,1)) and \
						barcode["GACNNNNNACC-1"]==(1, round(2.3,1)))

	def test_scan_motifs_aligned_r(self):
		mode = "aligned"
		strand = 1
		c_REF_STR = motif_tools.comp_motif(self.REF_STR)
		motifs = set(self.MOTIFS)
		bi_motifs = set()
		opts = data.set_opts()
		barcode, kmers = read_scanner.scan_motifs(mode, self.IPD_VALS, c_REF_STR, strand, motifs, bi_motifs, opts)
		self.assertTrue(barcode["ACCACC-0"]==(1, round(2.4,1)) and \
						barcode["TGATCA-2"]==(1, round(2.2,1)) and \
						barcode["GACNNNNNACC-1"]==(1, round(2.3,1)))

	def test_scan_motifs_unaligned(self):
		mode = "unaligned"
		strand = 0
		rc_REF_STR = motif_tools.rev_comp_motif(self.REF_STR)
		r_IPD_VALS = self.IPD_VALS[::-1]
		motifs = set(self.MOTIFS)
		bi_motifs = set()
		opts = data.set_opts()
		barcode, kmers = read_scanner.scan_motifs(mode, r_IPD_VALS, rc_REF_STR, strand, motifs, bi_motifs, opts)
		self.assertTrue(barcode["ACCACC-0"]==(1, round(2.4,1)) and \
						barcode["TGATCA-2"]==(1, round(2.2,1)) and \
						barcode["GACNNNNNACC-1"]==(1, round(2.3,1)))

	def test_walk_over_read_aligned_f(self):
		mode = "aligned"
		opts = data.set_opts()
		rc_REF_STR = motif_tools.rev_comp_motif(self.REF_STR)
		r_IPD_VALS = self.IPD_VALS[::-1]
		motifs = set(self.MOTIFS[:2])
		subread_ipds = {}
		for motif in motifs:
			subread_ipds[motif] = []
		strand = 0
		k = 6
		subread_ipds = read_scanner.walk_over_read(mode, subread_ipds, rc_REF_STR, r_IPD_VALS, strand, k, opts)
		self.assertTrue(round(subread_ipds["ACCACC-0"][0],1)==2.4 and \
						round(subread_ipds["TGATCA-2"][0],1)==2.2)

	def test_walk_over_read_aligned_r(self):
		mode = "aligned"
		opts = data.set_opts()
		c_REF_STR = motif_tools.comp_motif(self.REF_STR)
		motifs = set(self.MOTIFS[:2])
		subread_ipds = {}
		for motif in motifs:
			subread_ipds[motif] = []
		strand = 1
		k = 6
		subread_ipds = read_scanner.walk_over_read(mode, subread_ipds, c_REF_STR, self.IPD_VALS, strand, k, opts)
		self.assertTrue(round(subread_ipds["ACCACC-0"][0],1)==2.4 and \
						round(subread_ipds["TGATCA-2"][0],1)==2.2)

	def test_walk_over_read_unaligned(self):
		mode = "unaligned"
		opts = data.set_opts()
		rc_REF_STR = motif_tools.rev_comp_motif(self.REF_STR)
		r_IPD_VALS = self.IPD_VALS[::-1]
		motifs = set(self.MOTIFS[:2])
		subread_ipds = {}
		for motif in motifs:
			subread_ipds[motif] = []
		strand = 0
		k = 6
		subread_ipds = read_scanner.walk_over_read(mode, subread_ipds, rc_REF_STR, r_IPD_VALS, strand, k, opts)
		self.assertTrue(round(subread_ipds["ACCACC-0"][0],1)==2.4 and \
						round(subread_ipds["TGATCA-2"][0],1)==2.2)

	def test_walk_over_read_bipartite_aligned_f(self):
		mode = "aligned"
		opts = data.set_opts()
		rc_REF_STR = motif_tools.rev_comp_motif(self.REF_STR)
		r_IPD_VALS = self.IPD_VALS[::-1]
		motifs = set([self.MOTIFS[2]])
		subread_ipds = {}
		for motif in motifs:
			subread_ipds[motif] = []
		strand = 0
		subread_ipds = read_scanner.walk_over_read_bipartite(mode, subread_ipds, rc_REF_STR, r_IPD_VALS, strand, opts)
		self.assertTrue(round(subread_ipds["GACNNNNNACC-1"][0],1)==2.3 )

	def test_walk_over_read_bipartite_aligned_r(self):
		mode = "aligned"
		opts = data.set_opts()
		c_REF_STR = motif_tools.comp_motif(self.REF_STR)
		motifs = set([self.MOTIFS[2]])
		subread_ipds = {}
		for motif in motifs:
			subread_ipds[motif] = []
		strand = 1
		subread_ipds = read_scanner.walk_over_read_bipartite(mode, subread_ipds, c_REF_STR, self.IPD_VALS, strand, opts)
		self.assertTrue(round(subread_ipds["GACNNNNNACC-1"][0],1)==2.3 )

	def test_walk_over_read_bipartite_unaligned(self):
		mode = "unaligned"
		opts = data.set_opts()
		rc_REF_STR = motif_tools.rev_comp_motif(self.REF_STR)
		r_IPD_VALS = self.IPD_VALS[::-1]
		motifs = set([self.MOTIFS[2]])
		subread_ipds = {}
		for motif in motifs:
			subread_ipds[motif] = []
		strand = 0
		subread_ipds = read_scanner.walk_over_read_bipartite(mode, subread_ipds, rc_REF_STR, r_IPD_VALS, strand, opts)
		self.assertTrue(round(subread_ipds["GACNNNNNACC-1"][0],1)==2.3 )

if __name__ == '__main__':
    unittest.main()
