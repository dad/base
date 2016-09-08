import time, os, random, string, sys, math, traceback, unittest
from Bio.Phylo import Newick, BaseTree
import phyloutil

def buildTree(depth, breadth, parent, id):
	if depth>0:
		for i in range(breadth): # binary
			node = Newick.Clade()
			node.name = "{}.{}".format(id,i+1)
			parent.clades.append(node)
			print("Added node level {}".format(depth))
			buildTree(depth-1, breadth, node, node.name)



class test001(unittest.TestCase):
	def test_search(self):
		"""Find by name"""
		root = Newick.Clade()
		#root.parent = None
		root.name = 'root'
		buildTree(3,2,root,'root')
		for n in BaseTree._level_traverse(root, lambda x: x.clades):
			print(n)
		target_name
		target = phyloutil.findRootByName(target_name)
		self.assertTrue(target.name == target_name)


if __name__=="__main__":
	unittest.main(verbosity=2)
