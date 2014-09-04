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

if __name__=="__main__":
	unittest.main(verbosity=2)
