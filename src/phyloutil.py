import sys
from Bio import Phylo
from Bio.Phylo import Newick


def findNodeByName(name, node):
	"""Find node or child whose node.name matches name."""
	res = None
	match = Phylo.BaseTree._string_matcher(target=name)
	for n in Phylo.BaseTree._level_traverse(node, lambda x: x.clades):
		if match(n):
			res = n
			break
	return res

def parseClassificationTable(root, table):
	#node = findNodeByName(name, root)
	names = table['name']
	#print(names)
	names.reverse()
	cur_child = None
	#Phylo.write(root, sys.stdout, 'newick')
	for name in names:
		node = Newick.Clade()
		node.name = name
		node_found_in_tree = findNodeByName(node.name, root)
		if not cur_child is None:
			if not node_found_in_tree is None:
				# Insert the children into the tree, and exit
				node_found_in_tree.clades.append(cur_child)
				break
			else:
				node.clades.append(cur_child)
		cur_child = node
	#Phylo.write(cur_child, sys.stdout, 'newick')
