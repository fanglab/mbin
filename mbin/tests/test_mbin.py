#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `mbin` package."""


import unittest
import os
from mbin import mbin
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