import time, os, random, string, sys, math, traceback, unittest
import util
from Bio.Phylo import Newick, BaseTree
from Bio import Phylo
import lcr


class test001(unittest.TestCase):
	def test_it(self):
		"""test"""
		self.assertTrue(True)

if __name__=="__main__":
	unittest.main(verbosity=2)
