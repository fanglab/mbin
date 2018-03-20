#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `mbin` package."""

import unittest
import os
from mbin import motif_tools
from mbin import data

class TestMbin(unittest.TestCase):
	"""Tests for `mbin` package."""

	def setUp(self):
		"""Set up test fixtures, if any."""

	def tearDown(self):
		"""Tear down test fixtures, if any."""

	def test_cmph5_size(self):
		cmph5_d = data.getCmpH5()
		self.assertEqual(907224, os.path.getsize(cmph5_d["cmph5"]))

	def test_fasta_size(self):
		cmph5_d = data.getCmpH5()
		self.assertEqual(5748, os.path.getsize(cmph5_d["fasta"]))

	def test_bash5_size(self):
		bash5_d = data.getBasH5()
		self.assertEqual(260202, os.path.getsize(bash5_d["bash5"]))

	def test_read_motifs_from_file(self):
		opts = data.set_opts()
		motifs = motif_tools.motifs_from_file(opts)
		self.assertEqual(2,len(motifs))

	def test_build_contiguous_motifs(self):
		opts = data.set_opts()
		opts.bipartite = False
		motifs,bi_motifs = motif_tools.build_motif_sets(opts)
		self.assertEqual(7680,len(motifs))

	def test_build_bipartite_motifs(self):
		opts = data.set_opts()
		opts.min_kmer = 0
		opts.max_kmer = 0
		motifs,bi_motifs = motif_tools.build_motif_sets(opts)
		self.assertEqual(194560,len(bi_motifs))

	def test_add_degen_motifs(self):
		motifs_fn = data.getMotifsFile()
		motifs = set(["GWTC-1"])
		orig_means = {"GATC-1":1.0, "GTTC-1":0.8}
		new_means = motif_tools.add_degen_motifs(motifs,orig_means)
		self.assertEqual(0.90000000000000002,new_means.get("GWTC-1"))

	def test_sub_bases(self):
		motif = "ACHNTYA-0"
		subbed_motif = motif_tools.sub_bases(motif)
		self.assertEqual("AC[ACT][ACGTN]T[CT]A-0",subbed_motif)

	def test_comp_motif(self):
		motif = "ACTTGTTACA"
		comped_motif = motif_tools.comp_motif(motif)
		self.assertEqual("TGAACAATGT",comped_motif)

	def test_rev_comp_motif(self):
		motif = "ACTTGTTACA"
		rev_comped_motif = motif_tools.rev_comp_motif(motif)
		self.assertEqual("TGTAACAAGT",rev_comped_motif)

	def test_shorten_motifs(self):
		highscore_motifs = {"GATCA-1":2.1, "GATCC-1":2.0, "GATCG-1":1.9, "GATCT-1":2.2, "GATC-1":2.0 ,"CATG-1":1.7}
		keeper_motifs = motif_tools.shorten_motifs(highscore_motifs)
		self.assertTrue("GATC-1" in keeper_motifs and len(keeper_motifs)==2)

	# def test_run_mbin_cmph5(self):
	# 	cmph5_fn     = data.getCmpH5()["cmph5"]
	# 	contigs      = data.getCmpH5()["fasta"]
	# 	controls_dir = data.getControl()
	# 	opts         = data.set_opts()
	# 	args         = [cmph5_fn]
	# 	runner       = mbin.mbin_runner(opts, args)
	# 	runner.run()
	# 	self.assertEqual(1,1)

if __name__ == '__main__':
    unittest.main()