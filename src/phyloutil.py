from Bio import Phylo
from Bio.Phylo import Newick


def findNodeByName(name, node):
	"""Find node or child whose node.name matches name."""
	res = None
	match = Phylo.BaseTree._string_matcher(target=name)
	for n in node._level_traverse():
		if match(n):
			res = n
			break
	return res

def parseClassificationTable(root, table):
	node = findNodeByName(name, root)
	print(table)
