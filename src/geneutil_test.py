import time, os, random, string, sys, math, traceback, unittest
import geneutil


class test001(unittest.TestCase):
	def test_longest_run(self):
		"""Longest run testcases"""
		self.assertTrue(geneutil.longestRun('AAAAA','A')==5)
		self.assertTrue(geneutil.longestRun('AAATAA','A',1)==6)
		self.assertTrue(geneutil.longestRun('AAATTAA','A',1)==3)
		self.assertTrue(geneutil.longestRun('AAATTAA','A',2)==7)
		self.assertTrue(geneutil.longestRun('TAAATAA','A',1)==6)
		self.assertTrue(geneutil.longestRun('TAAATAAT','A',1)==6)

	def test_longest_run_mult(self):
		"""Longest run testcases with more than one target"""
		self.assertTrue(geneutil.longestRun('QQQQN','QN')==5)
		self.assertTrue(geneutil.longestRun('QQANNQ','QN',1)==6)
		self.assertTrue(geneutil.longestRun('QQNPPQ','QN',1)==3)
		self.assertTrue(geneutil.longestRun('QQQAANN','QN',2)==7)
		self.assertTrue(geneutil.longestRun('ANQNQAN','QN',1)==6)
		self.assertTrue(geneutil.longestRun('ANQNQANP','QN',1)==6)
	
	def test_max_sliding_count(self):
		"""Max Sliding Count testcases"""
		self.assertTrue(geneutil.maxSlidingCount('AAAAA','A')==5)
		self.assertTrue(geneutil.maxSlidingCount('AAAAA','Q')==0)
		self.assertTrue(geneutil.maxSlidingCount('AAATAA','A')==4)
		self.assertTrue(geneutil.maxSlidingCount('AAATTAA','A')==3)
		self.assertTrue(geneutil.maxSlidingCount('MMMMMMMMMMABCABCABCDM','M',10)==10)
		self.assertTrue(geneutil.maxSlidingCount('MMMMMMMMMMABCABCABCDM','C',10)==3)

if __name__=="__main__":
	unittest.main(verbosity=2)
