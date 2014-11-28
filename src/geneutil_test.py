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

class test002(unittest.TestCase):
	def test_gapped_find(self):
		"""Gapped-find testcases"""
		self.assertTrue(geneutil.gappedFind('AAASS--SAA','SSS')==3)
		self.assertTrue(geneutil.gappedFind('AAASS--SAA','SSSS')==-1)
		#print geneutil.gappedFind('AAASS--SAA','SSS',start=False)
		self.assertTrue(geneutil.gappedFind('AAASS--SAA','SSS',start=False)==8)
		self.assertTrue(geneutil.gappedFind('AAASSxxSAA','SSS',start=False)==-1)
		self.assertTrue(geneutil.gappedFind('AAASSxxSAA','SSS',start=False,gap='x')==8)

if __name__=="__main__":
	unittest.main(verbosity=2)
