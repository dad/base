import sys
from random import shuffle
from Bio import Phylo
from Bio.Phylo import Newick
from Bio.Phylo import NewickIO
from io import StringIO

'''
Tests and utility functions for low-complexity regions (LCRs).
'''

def randomizeResidues(sequences):
	"""Preserve sequence length, randomize states (residues) between sequences, return randomized sequences"""
	# General strategy (no optimization): build list of all residues in all sequences,
	# then sample from list to rebuild sequences.
	all_residues = []
	for s in sequences:
		all_residues += [x for x in s]
	# Randomize the residues
	shuffle(all_residues)
	# Now chip off new sequences, each according to the length of the current sequence
	new_seqs = []
	cur_index = 0
	for s in sequences:
		new_index = cur_index+len(s)
		new_seq = ''.join(all_residues[cur_index:new_index])
		new_seqs.append(new_seq)
		cur_index = new_index
	# Ensure we've used all the residues
	assert new_index == len(all_residues)
	return new_seqs



