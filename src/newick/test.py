import sys, os, math, string
sys.path = [os.path.expanduser("~/research/lib")] + sys.path
import newick

def main():
	# A yeast phylogeny
	tree_string = "((((sbay:50.000000,scas:0.000004)8:11.818649,(smik:6.760546,spar:8.815752)9:19.825863)7:0.000000,sklu:0.000004)6);"
	print tree_string
	t = newick.tree.parseTree(tree_string)
	print t
	leaves = t.leaves
	for i in range(len(leaves)-1):
		for j in range(i+1,len(leaves)):
			print leaves[i].name, leaves[j].name, leaves[i].distanceFrom(leaves[j])
		

main()
