import time, os, random, string, sys, math, traceback, unittest
from Bio.Phylo import Newick
import phyloutil

def buildTree(depth, breadth, parent):
	for n in range(depth):
		for i in range(breadth): # binary
			node = Newick.Clade()
			node.name = "{}.{}".format(n,i)
			node.parent = parent
			buildTree(depth-1, breadth, node)



class test001(unittest.TestCase):
	def test_search(self):
		"""Find by name"""
		root = Newick.Clade()
		root.parent = None
		buildTree(3,2,root)
		print(root)

		#self.assertTrue(geneutil.longestRun('AAAAA','A')==5)


if __name__=="__main__":
	unittest.main(verbosity=2)
