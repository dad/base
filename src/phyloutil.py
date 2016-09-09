import sys
from Bio import Phylo
from Bio.Phylo import Newick

class TreeError(Exception):
	def __init__(self, value):
		self.value = value

	def __str__(self):
		return repr(self.value)

def printTree(root, out=sys.stdout, indent='', perlevel='  ', iterating=False):
	if not iterating:
		out.write("\n")
	outstr = ''
	if not root.name is None:
		outstr = root.name
	out.write(indent+outstr+"\n")
	for clade in root.clades:
		printTree(clade, out, indent+perlevel, iterating=True)

def cladeNamesEqual(clade1, clade2):
	return clade1.name == clade2.name


def findNodeByName(name, tree):
	"""Find node in tree whose node.name matches name."""
	res = None
	match = Phylo.BaseTree._string_matcher(target=name)
	for n in Phylo.BaseTree._level_traverse(tree, lambda x: x.clades):
		if match(n):
			res = n
			break
	return res

def findNode(node, tree, clades_equal=cladeNamesEqual):
	"""Find node in tree whose node.name matches name."""
	res = None
	for n in Phylo.BaseTree._level_traverse(tree, lambda x: x.clades):
		if clades_equal(node, n):
			res = n
			break
	return res


def mergeTrees(master, twig, clades_equal=cladeNamesEqual):
	"""Merge twig into master"""
	master_node = None
	if not clades_equal(master, twig):
		# Root is not shared
		# Find twig's root in master
		master_node = findNode(twig, master, clades_equal)
	else:
		master_node = master
	if master_node is None:
		# Bail: 
		raise TreeError("Can't find twig root '{}' in master tree".format(twig.name))
	# Now step down twig and master in tandem
	# until twig has information not in master
	cur_master = master_node
	cur_twig = twig
	while clades_equal(cur_master, cur_twig):
		pass #cur_twig = cur_twig.clades[]

def parseClassificationTable(table):
	#node = findNodeByName(name, root)
	names = table['name']
	#print(names)
	names.reverse()
	cur_child = None
	for name in names:
		#print(name)
		node = Newick.Clade()
		node.name = name
		#node_found_in_tree = findNodeByName(node.name, root)
		# DAD: actually need to start from the root.
		if not cur_child is None:
			node.clades.append(cur_child)
		cur_child = node
	return node
'''
	# Now we have constructed a linear twig.
	# Merge this twig into the root tree
	mergeTrees
			#Phylo.draw_ascii(root)
			if not node_found_in_tree is None:
				# Insert the children into the tree, and exit
				node_found_in_tree.clades.append(cur_child)
				#print("Node {} found\n".format(node.name))
				break
			else:
				node.clades.append(cur_child)
		else:
			# This is the leaf node
			#node.
			pass
		cur_child = node
	#Phylo.write(cur_child, sys.stdout, 'newick')
'''
