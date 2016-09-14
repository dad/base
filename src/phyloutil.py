import sys
from Bio import Phylo
from Bio.Phylo import Newick
from Bio.Phylo import NewickIO
from io import StringIO

class TreeError(Exception):
	def __init__(self, value):
		self.value = value

	def __str__(self):
		return repr(self.value)

def printTree(root, out=sys.stdout, indent='', perlevel='  ', iterating=False):
	"""Simple print method for a Phylo tree"""
	if not iterating:
		out.write("\n")
	outstr = ''
	if not root.name is None:
		outstr = root.name
	out.write(indent+outstr+"\n")
	for clade in root.clades:
		printTree(clade, out, indent+perlevel, iterating=True)

def extractSpeciesName(text):
	"""Extract species name from NCBI-style header"""
	res = text.split('[')[1].split(']')[0]
	return res

def cladeNamesEqual(clade1, clade2):
	return clade1.name == clade2.name

def findNodeByName(name, tree):
	"""Find node in tree whose node.name matches name."""
	res = None
	found = tree.find_clades({"name":name})
	try:
		res = next(found)
	except:
		pass
	return res

def readOneTree(stream):
	"""Reads a Newick-formatted tree, permitting lines with comments denoted by leading '#'."""
	tree_string = ""
	lines = stream.readlines()
	for line in lines:
		if not line.strip()[0] == '#':
			tree_string += line.strip()
	trees = NewickIO.parse(StringIO(tree_string))
	tree = next(trees)
	return tree

def buildNodeLookupFromTree(tree, allow_collisions=False):
	    names = {}
	    for clade in tree.find_clades():
	        if clade.name:
	            if not allow_collisions and (clade.name in names):
	                raise ValueError("Duplicate clade name in tree: '{}'".format(clade.name))
	            names[clade.name] = clade
	    return names

def mergeTrees(master, twig, add_to_leaf):
	"""Merge twig into master"""
	twig_clades = []
	cur = twig
	while len(cur.clades)>0:
		twig_clades.append(cur)
		cur = cur.clades[0]
	twig_clades.append(cur)
	#twig_clades = list(Phylo.BaseTree._level_traverse(twig, lambda x: x.clades))
	twig_clades.reverse()
	#for t in twig_clades:
	#	print("\t{}".format(t.name))
	last_twig_clade = None
	added = False
	for clade in twig_clades:
		#print(clade.name)
		# Look up 
		master_nodes = [x for x in master.find_clades({"name": clade.name})]
		if len(master_nodes)>0:
			# Assert that this is an actual match -- just checking!
			master_node = master_nodes[0]
			assert master_node.name == clade.name
			# Assert that we've not already added this clade before
			assert not clade.name in [x.name for x in master_node.clades]
			#print("found node {}, last twig clade {}".format(master_node.name, last_twig_clade.name))
			if not last_twig_clade is None:
				if len(master_node.clades)==0:
					if add_to_leaf:
						master_node.clades.append(last_twig_clade)
						added = True
				else:
					master_node.clades.append(last_twig_clade)
					added = True
				break
		last_twig_clade = clade
	return added

def treeFromClassificationTable(table):
	names = table['name']
	names.reverse()
	cur_child = None
	node = None
	for name in names:
		node = Newick.Clade()
		node.name = name
		if not cur_child is None:
			node.clades.append(cur_child)
		cur_child = node
	return node

def getClassificationEntryByRank(table, rank):
	res = None
	for row in table.dictrows:
		if row['rank'] == rank:
			res = row
	return res
