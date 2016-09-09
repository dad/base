import sys
from Bio import Phylo
from Bio.Phylo import Newick

def printTree(root, out, indent='', perlevel='  ', iterating=False):
	if not iterating:
		out.write("\n")
	outstr = ''
	if not root.name is None:
		outstr = root.name
	out.write(indent+outstr+"\n")
	for clade in root.clades:
		printTree(clade, out, indent+perlevel, iterating=True)


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
		#print(name)
		node = Newick.Clade()
		node.name = name
		node_found_in_tree = findNodeByName(node.name, root)
		if not cur_child is None:
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
