import time, os, random, string, sys, math, traceback, unittest
import util
from Bio.Phylo import Newick, BaseTree
from Bio import Phylo
import phyloutil

def buildTree(depth, breadth, parent, id):
	if depth>0:
		for i in range(breadth): # binary
			node = Newick.Clade()
			node.name = "{}.{}".format(id,i+1)
			parent.clades.append(node)
			#print("Added node level {}".format(depth))
			buildTree(depth-1, breadth, node, node.name)



class test001(unittest.TestCase):
	def test_search(self):
		"""Find by name"""
		root = Newick.Clade()
		#root.parent = None
		root.name = 'root'
		buildTree(3,2,root,'root')
		#for n in BaseTree._level_traverse(root, lambda x: x.clades):
		#	print(n)
		target_name = "root.1.2.1"
		target = phyloutil.findNodeByName(target_name, root)
		self.assertTrue(target.name == target_name)

	def test_reading_from_class_table(self):
		"""Read table"""
		root = Newick.Clade()
		root.name = "cellular organisms"
		inf = open("./test-phyloutil/test1/Pseudozyma-antarctica-1.txt", 'r')
		table = util.readTable(inf, header=True)
		#print(table)
		phyloutil.parseClassificationTable(root, table)
		inf.close()
		#phyloutil.printTree(root, sys.stdout)
		termlist = list(root.get_terminals())
		self.assertTrue(termlist[0].name=='Moesziomyces antarcticus T-34')

	def test_reading_from_guide_table(self):
		"""Read table"""
		root = Newick.Clade()
		root.name = "cellular organisms"
		#print(root.depths())
		inf = open("./test-phyloutil/test1/Pseudozyma-antarctica-1.txt", 'r')
		table = util.readTable(inf, header=True)
		#print(table)
		phyloutil.parseClassificationTable(root, table)
		inf.close()
		phyloutil.printTree(root, sys.stdout)
		#print(root.depths())
		#Phylo.draw_ascii(root)

	def test_simple_manual_tree(self):
		"""Manually constructed tree"""
		root = Newick.Clade()
		root.name = "root"
		child = Newick.Clade()
		child.name = "kid1"
		root.clades.append(child)
		child = Newick.Clade()
		child.name = "kid2"
		root.clades.append(child)
		phyloutil.printTree(root, sys.stdout)
		#Phylo.draw_ascii(root)

	def test_print(self):
		"""Print tree"""
		root = Newick.Clade()
		#root.parent = None
		root.name = 'root'
		buildTree(3,2,root,'root')
		phyloutil.printTree(root, sys.stdout)
		print(root)

if __name__=="__main__":
	unittest.main(verbosity=2)
