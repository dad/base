#! /usr/local/bin/python

import sys, os, math, string, random, pickle
import newick, copy

if __name__=="__main__":
	print "** parsing tree and copying tree"
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
	print "** removing leaves"
	assert str(tree2) == tree_str
	tree2 = tree2.removeLeaf(node_dict2["sbay"])
	assert str(tree2) == '(scas,((scer,spar),smik));'
	tree2 = tree2.removeLeaf(node_dict2["scas"])
	assert str(tree2) == '((scer,spar),smik);'
	tree2 = tree2.removeLeaf(node_dict2["smik"])
	assert str(tree2) == '(scer,spar);'

	print "** isRooted() and unroot()"
	assert tree.isRooted()
	#print tree.root
	leaves = set(tree.leaves)
	#print tree
	tree.unroot()
	#print tree
	assert not tree.isRooted() and set(tree.leaves) == leaves

	t3 = newick.tree.parseTree("(scer,spar);")
	assert not t3.isRooted()

	print "** all nodes have root equal to tree"
	for n in tree.nodes:
		assert n.root == tree

	print "** only one root"
	num_roots = 0
	for n in tree.nodes:
		if n.isRoot():
			num_roots += 1
	assert num_roots == 1

	print "** All tests passed **"
