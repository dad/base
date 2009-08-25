#! /usr/local/bin/python

import sys, os, math, string, random, pickle
import newick, copy

if __name__=="__main__":
	tree_str = "((((scer,spar),smik),sbay),scas);"
	tree = newick.tree.parseTree(tree_str)
	node_dict = dict([(x.name, x) for x in tree.nodes])
	tree2 = copy.deepcopy(tree)
	node_dict2 = dict([(x.name, x) for x in tree2.nodes])
	assert str(tree2) == str(tree)

	# Name nodes on whole tree
	if False:
		for n in tree2.nodes:
			if n.name is None:
				leaf_names = [x.name for x in n.leaves]
				leaf_names.sort()
				n.name = '_'.join(leaf_names)
			print n.name
	assert str(tree2) == tree_str
	tree2 = tree2.removeLeaf(node_dict2["sbay"])
	assert str(tree2) == '(scas,((scer,spar),smik));'
	tree2 = tree2.removeLeaf(node_dict2["scas"])
	assert str(tree2) == '((scer,spar),smik);'
	tree2 = tree2.removeLeaf(node_dict2["smik"])
	assert str(tree2) == '(scer,spar);'
	print "** All tests passed **"
	#print tree2
	