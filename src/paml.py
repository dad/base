#!/usr/bin/python
# Begin paml.py
"""Module for running PAML programs.

Original version by Jesse Bloom.
Expanded and maintained by D. Allan Drummond, 2004-2008.
treeDnDsAtCodon and Node class by Claus O. Wilke, 2006-2007.
"""

import shutil, sys, os, re, string, math, random
import codon, newick

class Node:
	def __init__( self, sequence=None, left=None, right=None ):
		self.sequence = sequence
		self.left = left
		self.right = right

	def printTree( self, level=0 ):
		if self.left != None:
			self.left.printTree( level+1 )
		print "      "*level, self.sequence
		if self.right != None:
			self.right.printTree( level+1 )

	def printCodonTree( self, pos, level=0 ):
		if self.left != None:
			self.left.printCodonTree( pos, level+1 )
		print "      "*level, self.sequence[pos:pos+3]
		if self.right != None:
			self.right.printCodonTree( pos, level+1 )

	def allCodonsValid( self, pos ):
		if not codon.codonValid( self.sequence[pos:pos+3] ):
			return False
		else:
			valid = True
			if self.left != None:
				valid = valid and self.right.allCodonsValid( pos )
			if valid:
				if self.right != None:
					return valid and self.left.allCodonsValid( pos )
				else:
					return True
			else:
				return False


def printData( c1, c2, Dn, Ds, N, S ):
	"For debugging purposes only"
	if debug:
		print "%s %s Dn=%g, Ds=%g, N=%g, S=%g, Dn/N=%g, Ds/S=%g" % ( c1, c2, Dn, Ds, N, S, Dn/N, Ds/S )

def treeDnDsAtCodon( node, pos ):
	"Calculate Dn and Ds over entire tree, at position 'pos' only. "
	if node.left == None or node.right == None:
		return ( 0., 0. )
	croot = node.sequence[pos:pos+3]
	cleft = node.left.sequence[pos:pos+3]
	cright = node.right.sequence[pos:pos+3]
	( N, S ) = codon.calcNS( croot )
	#if S==0:
	#	return ( -100000., -100000. )
	( Dn1, Ds1 ) = codon.calcDnDs( croot, cleft )
	#printData( croot, cleft, Dn1, Ds1, N, S )
	( Dn2, Ds2 ) = codon.calcDnDs( croot, cright )
	#printData( croot, cright, Dn2, Ds2, N, S )
	( Dnleft, Dsleft ) = treeDnDsAtCodon( node.left, pos )
	( Dnright, Dsright ) = treeDnDsAtCodon( node.right, pos )
	codon_ds = 0.0
	if S > 0:
		codon_ds = (Ds1 + Ds2)/S
	return ( Dnleft + Dnright + (Dn1 + Dn2)/N, Dsleft + Dsright + codon_ds )

#----------------------------------------------------------------------------------
class PAMLError(Exception):
	"""Error running or processing PAML."""
#-----------------------------------------------------------------------------

def Get_Distance(seq1, seq2, seq_type):
	"""Computes the evolutionary distance(s) between aligned sequences.

	'seq1' and 'seq2' are two aligned sequences, with gaps denoted as '-'.
	Uses the PAML program CodeML to compute the distance(s) between the
    two sequences, and returns it.  The sequences are either aligned
	protein sequences or DNA sequences aligned according to a protein
	alignment.
	  'type' specifies whether we are comparing protein sequences
	  (type = 'protein') or nucleotide codon sequences (type = 'codon').
	  If 'type' is 'protein', returns a single scalar representing the distance.
	  If 'type' is 'codon', returns a 6-tuple with the following elements:
		(maximum likelihood (ML) distance, ML synonomous rate,
		ML non-synonomous rate, Nei and Gojobori (NG) distance,
		NG synonomous rate, NG non-nysynonomous rate)
	Also works if seq1 and seq2 are lists of aligned sequences (they must
	be the same length).  The returned results are then as for the
	single case but are lists of tuples or scalars."""
    #See if we are aligning single sequences or lists of sequences
	#If they are strings, make them lists of length 1
	if isinstance(seq1, str): # sequence is string, so single sequences
		if not isinstance(seq2, str):
			raise PAMLError, "One sequence is a string but the other is not."
		seq1 = [seq1]
		seq2 = [seq2]
	# Now we should have lists
	assert isinstance(seq1, list) and isinstance(seq2, list)
	assert len(seq1) == len(seq2), "Error, aligned sequence lists differ in length."
	if not seq_type in ['protein', 'codon']:
		raise PAMLError, "Distance type of %s is invalid." % seq_type
	cm = CodeML(seq_type)
	phylip_file = 'Get_Distance_phylip.temp'
	codeml_file = 'Get_Distance_codeml.temp'
	tree_file = 'Get_Distance_tree.temp'
	dist = []
	for i in range(len(seq1)):
		s1 = seq1[i]
		s2 = seq2[i]
		assert len(s1) == len(s2), "Error, sequences differ in length for " + str(i)
	 	Write_Phylip(phylip_file, s1, s2)
	 	writeTree(tree_file, s1, s2)
		(syn, ns) = cm.Run(phylip_file, codeml_file, tree_file)
	   	d = cm.Get_Dist()
	 	if d == None:
			print "Error, did not compute protein distance for ", i, ":", d
			return
		dist.append(d)
	osremove(phylip_file)
	osremove(codeml_file)
	assert len(dist) == len(seq1) == len(seq2), "Error, distances differ."
	# Now if we have single sequences return single value,
	# else return list
	if len(dist) == 1:
		return dist[0]
	else:
		return dist

def Get_Distance_NS(seq1, seq2, seq_type):
	"""Computes the evolutionary distance(s) between aligned sequences.

   'seq1' and 'seq2' are two aligned sequences, with gaps denoted as
   '-'.
   Uses the PAML program CodeML to compute the distance(s) between the

    two sequences, and returns it.  The sequences are either aligned
	protein sequences or DNA sequences aligned according to a protein
	alignment.
	  'type' specifies whether we are comparing protein sequences
	  (type = 'protein') or nucleotide codon sequences (type =
	'codon').
	  If 'type' is 'protein', returns a single scalar representing the
	distance.
	  If 'type' is 'codon', returns a 6-tuple with the following
	elements:
	(maximum likelihood (ML) distance, ML synonomous rate,
	ML non-synonomous rate, Nei and Gojobori (NG) distance,
	NG synonomous rate, NG non-nysynonomous rate, number of synonymous sites,
	number of nonsynonymous sites)
	  Also works if seq1 and seq2 are lists of aligned sequences (they
	must
	be the same length).  The returned results are then as for the
	single case but are lists of tuples or scalars."""
    #See if we are aligning single sequences or lists of sequences
	#If they are strings, make them lists of length 1
	if isinstance(seq1, str): # sequence is string, so single sequences
		if not isinstance(seq2, str):
			raise PAMLError, "One sequence is a string but the other is not."
		seq1 = [seq1]
		seq2 = [seq2]
	# Now we should have lists
	assert isinstance(seq1, list) and isinstance(seq2, list)
	assert len(seq1) == len(seq2), "Error, aligned sequence lists differ in length."
	if not seq_type in ['protein', 'codon']:
		raise PAMLError, "Distance type of %s is invalid." % seq_type
	cm = CodeML(seq_type)
	phylip_file = 'Get_Distance_phylip.temp'
	codeml_file = 'Get_Distance_codeml.temp'
	tree_file = 'Get_Distance_tree.temp'
	dist = []
	for i in range(len(seq1)):
		s1 = seq1[i]
		s2 = seq2[i]
		assert len(s1) == len(s2), "Error, sequences differ in length for " + str(i)
	 	Write_Phylip(phylip_file, s1, s2)
	 	writeTree(tree_file, s1, s2)
		tempfilename = cm.Run(phylip_file, codeml_file, tree_file)
		(syn, ns) = cm.Get_NS_Sites(tempfilename)
	   	d = cm.Get_Dist()
	 	if d == None:
			print "Error, did not compute protein distance for ", i, ":", d
			return
		d += (syn, ns)
		dist.append(d)
	osremove(phylip_file)
	osremove(codeml_file)
	assert len(dist) == len(seq1) == len(seq2), "Error, distances differ."
	# Now if we have single sequences return single value,
	# else return list
	if len(dist) == 1:
		return dist[0]
	else:
		return dist

def Get_Distance_NS_Physical(seq1, seq2, seq_labels=None, tree_string=None ):
	"""Computes the evolutionary distance(s) between aligned sequences with
	a specified phylogeny.  Uses physical definition of sites, rather than
	typical mutational-opportunity definition.  See Bierne and Eyre-Walker
	MBE 2003.  Uses FMutSel-F model, see Yang & Nielsen MBE 2008.

	'seqs' is a list of aligned sequences.  'tree' is a phylogeny in the format
	Uses the PAML program CodeML to compute the distance(s) between the
    two sequences, and returns it.  The sequences must be DNA sequences aligned according to a protein
	alignment.

	Returns a four-tuple
		(ML nonsynonymous rate, ML synonomous rate, number of nonsynonymous sites, number of synonymous sites)"""
	seq_type = 'codon'
	if not seq_type in ['protein', 'codon']:
		raise PAMLError, "Distance type of %s is invalid." % seq_type
	cm = CodeML(seq_type)
	codeml_outfile = 'get-distance-tmp.txt'
	codeml_infile = 'get-distance-infile-tmp.txt'
	tree_file = 'get-distance-tree-tmp.txt'

	seqs = [seq1,seq2]
	if not seq_labels:
		seq_labels = ['seq_%d' % (i+1,) for i in range(len(seqs))]
	if not tree_string:
		tree_string = "(" + ','.join(seq_labels) + ")"

 	writeMultipleTree(tree_file, seqs, tree_string)
 	writeMultipleSequences(codeml_infile, seqs, seq_labels)
 	option_dict = dict([("model","0"),
 		("NSsites","0"),  # one omega for tree
 		("runmode","0"),  # user tree
		("RateAncestor","0"),  # don't reconstruct the ancestral states
		("noisy","9"),  # maximize output
		("CodonFreq","7"),  # FMutSel model
		("estFreq","0"),  # -F model (no ML estimation of codon frequencies; compute them from data)
		("seqtype","1")  # lots of output
		])
	tempfilename = cm.Run(codeml_infile, codeml_outfile, tree_file, option_dict)
   	(ignore1, ignore2, nn, ns) = cm.getMultipleDistPhysical(tempfilename)
   	rst_fname = "rst1"
   	(dn, ds) = cm.getPhysicalDistRST(rst_fname)
	#osremove(codeml_infile)
	#osremove(codeml_outfile)
	#osremove(tree_file)
	#osremove(tempfilename)
	#osremove(rst_fname)
	return (dn, ds, nn, ns)

def Get_Tree_Distance_NS(seqs, seq_labels=None, tree_string=None, options=[] ):
	"""Computes the evolutionary distance(s) between aligned sequences with
	a specified phylogeny.

	'seqs' is a list of aligned DNA sequences.  'tree' is a phylogeny in Newick format
	Uses the PAML program CodeML to compute the distance(s) between the
    two sequences, and returns it.  The sequences must be DNA sequences aligned according to a protein
	alignment.

	Returns a four-tuple:
		(maximum likelihood (ML) nonsynonymous rate, ML synonymous rate,
		number of nonsynonymous sites, number of synonymous sites, transition-transversion ratio (kappa)"""
	seq_type = 'codon'
	if not seq_type in ['protein', 'codon']:
		raise PAMLError, "Distance type of %s is invalid." % seq_type
	cm = CodeML(seq_type)
	codeml_outfile = 'get-distance-tmp.txt'
	codeml_infile = 'get-distance-infile-tmp.txt'
	tree_file = 'get-distance-tree-tmp.txt'

	if not seq_labels:
		seq_labels = ['seq_%d' % (i+1,) for i in range(len(seqs))]
	if not tree_string:
		tree_string = "(" + ','.join(seq_labels) + ")"

 	writeMultipleTree(tree_file, seqs, tree_string)
 	writeMultipleSequences(codeml_infile, seqs, seq_labels)
 	option_dict = dict([("model","0"), ("NSsites","1"), ("runmode","0"), ("RateAncestor","0"), \
						("noisy","9"), ("CodonFreq","2"), ("seqtype","1")]+options)
	tempfilename = cm.Run(codeml_infile, codeml_outfile, tree_file, option_dict)

	#'2ML.t', '2ML.dS', '2ML.dN', '2NG.t', '2NG.dS', '2NG.dN'
	(ns, ns) = cm.Get_NS_Sites(tempfilename)
	d = cm.Get_Dist()
	(t, ds, dn, ngt, dsng, dnng) = d
	osremove(codeml_infile)
	osremove(codeml_outfile)
	osremove(tree_file)
	return (dn, ds, nn, ns)


def Get_Tree_Distance_Physical(seqs, seq_labels=None, tree_string=None, options=[] ):
	"""Computes the evolutionary distance(s) between aligned sequences with
	a specified phylogeny.  Uses physical definition of sites, rather than
	typical mutational-opportunity definition.  See Bierne and Eyre-Walker
	MBE 2003.

	'seqs' is a list of aligned sequences.  'tree' is a phylogeny in Newick format.
	Uses the PAML program CodeML to compute the distance(s) between the
    two sequences, and returns it.  The sequences must be DNA sequences aligned according to a protein
	alignment.
	Returns a four-tuple:
		(maximum likelihood (ML) nonsynonymous rate, ML synonymous rate,
		number of nonsynonymous sites, number of synonymous sites)"""
	seq_type = 'codon'
	if not seq_type in ['protein', 'codon']:
		raise PAMLError, "Distance type of %s is invalid." % seq_type
	cm = CodeML(seq_type)
	codeml_outfile = 'get-distance-tmp.txt'
	codeml_infile = 'get-distance-infile-tmp.txt'
	tree_file = 'get-distance-tree-tmp.txt'
	if not seq_labels:
		seq_labels = ['seq_%d' % (i+1,) for i in range(len(seqs))]
	if not tree_string:
		tree_string = "(" + ','.join(seq_labels) + ")"

 	writeMultipleTree(tree_file, seqs, tree_string)
 	writeMultipleSequences(codeml_infile, seqs, seq_labels)
 	option_dict = dict([("model","0"), ("NSsites","1"), ("runmode","0"), ("RateAncestor","0"), \
						("noisy","9"), ("CodonFreq","2"), ("seqtype","1")]+options)
	tempfilename = cm.Run(codeml_infile, codeml_outfile, tree_file, option_dict)
   	(dn, ds, nn, ns) = cm.Get_Multiple_Dist_Physical(tempfilename)
	osremove(codeml_infile)
	osremove(codeml_outfile)
	osremove(tree_file)
	return (dn, ds, nn, ns)

def Get_Tree_Distance_Physical_Kappa(seqs, seq_labels=None, tree_string=None, options=[] ):
	return getTreeDistancePhysicalKappa(seqs, seq_labels, tree_string, options)

def getTreeDistancePhysicalKappa(seqs, seq_labels=None, tree_string=None, options=[] ):
	"""Computes the evolutionary distance(s) between aligned sequences with
	a specified phylogeny.  Uses physical definition of sites, rather than
	typical mutational-opportunity definition.  See Bierne and Eyre-Walker
	MBE 2003.

	'seqs' is a list of aligned sequences.  'tree' is a phylogeny in Newick format.
	Uses the PAML program CodeML to compute the distance(s) between the
    two sequences, and returns it.  The sequences must be DNA sequences aligned according to a protein
	alignment.
	Returns a five-tuple:
		(maximum likelihood (ML) nonsynonymous rate, ML synonymous rate,
		number of nonsynonymous sites, number of synonymous sites, transition-transversion rate ratio (kappa)"""
	seq_type = 'codon'
	if not seq_type in ['protein', 'codon']:
		raise PAMLError, "Distance type of %s is invalid." % seq_type
	cm = CodeML(seq_type)
	codeml_outfile = 'get-distance-tmp.txt'
	codeml_infile = 'get-distance-infile-tmp.txt'
	tree_file = 'get-distance-tree-tmp.txt'

	if not seq_labels:
		seq_labels = ['seq_%d' % (i+1,) for i in range(len(seqs))]
	if not tree_string:
		tree_string = "(" + ','.join(seq_labels) + ")"

 	writeMultipleTree(tree_file, seqs, tree_string)
 	writeMultipleSequences(codeml_infile, seqs, seq_labels)
 	option_dict = dict([
 		("model","0"),
 		("NSsites","0"),      # one omega for tree
 		("runmode","0"),      # user tree
		("RateAncestor","0"), # don't reconstruct the ancestral states
		("noisy","9"),        # maximize output
		("CodonFreq","7"),    # FMutSel model
		("estFreq","0"),      # -F model (no ML estimation of codon frequencies; compute them from data)
		("fix_kappa","0"),    # don't fix kappa, estimate it
		("seqtype","1")       # codons
		])
	# Add user options
	for (k,v) in options:
		option_dict[k] = v
	tempfilename = cm.Run(codeml_infile, codeml_outfile, tree_file, option_dict)
   	(ignore1, ignore2, nn, ns) = cm.getMultipleDistPhysical(tempfilename)
   	rst_fname = "rst1"
   	(dn, ds) = cm.getPhysicalDistRST(rst_fname)
	kappa = cm.getKappa(codeml_outfile)
	if False:
		osremove(codeml_infile)
		osremove(codeml_outfile)
		osremove(tree_file)
	#print (dn, ds, nn, ns, kappa)
	return (dn, ds, nn, ns, kappa)

def Get_Site_Rates_And_Ancestor(seqs, seq_labels=None, tree_string=None ):
	"""Reconstructs an ancestor and per-site substitutions from an alignment.

	'seqs' is a list of aligned sequences.  'tree_string' is a phylogeny in Newick format needed by PAML,
	e.g. '((seq<n>,seq<n+1>),(seq<n+2>,seq<n+3>))'
	Uses the PAML program codeml to reconstruct the ancestral sequence by maximum likelihood
	and returns the ancestor and the per-site total synonymous and nonsynonymous
	changes in the tree, in a list of (syn,nsyn) pairs, according to the Suzuki and Gojobori MBE 16(10) (1999) method."""
	#if not seq_type in ['protein', 'codon']:
	#	raise PAMLError, "Distance type of %s is invalid." % seq_type
	seq_type = "codon"

	if len(seqs)<2:
		ancestor = seqs[0]
		rates = [(0.0,0.0)]*len(seqs[0])
		return (rates, ancestor)

	cm = CodeML(seq_type, control_file = "codeml.ctl")
	codeml_file = 'get-rate-ancestor-tmp.txt'
	codeml_infile = 'get-rate-ancestor-infile-tmp.txt'
	tree_file = 'get-rate-ancestor-tree-tmp.txt'

	if not seq_labels:
		seq_labels = ['seq_%d' % (i+1,) for i in range(len(seqs))]
	if not tree_string:
		tree_string = "(" + ','.join(seq_labels) + ")"
	# Remove gaps in aligned seqs before reconstruction
	n_seqs = len(seqs)
	new_seqs = ['']*n_seqs
	ngaps = 0
	for i in range(len(seqs[0])/3):
		codons = [s[3*i:3*i+3] for s in seqs]
		if not ('---' in codons):
			for j in range(n_seqs):
				new_seqs[j] += codons[j]
		else:
			ngaps += 1
	# Perhaps there are no remaining residues.  PAML chokes if this is so.
	if len(new_seqs[0]) == 0:
		raise PAMLError, "No ungapped residues in input sequences!"
 	writeMultipleTree(tree_file, new_seqs, tree_string)
 	writeMultipleSequences(codeml_infile, new_seqs, seq_labels)
 	option_dict = dict([("model","0"), ("NSsites","0"), ("runmode","0"),\
						("RateAncestor","1"),("Small_Diff","1e-6"),("verbose","2"),("seqtype","1")])
	tempfilename = cm.Run(codeml_infile, codeml_file, tree_file, option_dict)
	(rawrates, rawancestor) = cm.Get_Site_Rates_And_Ancestor()

	# Add the gaps back in
	rates = []
	ancestor = ''
	j = 0
	for i in range(len(seqs[0])/3):
		codons = [s[3*i:3*i+3] for s in seqs]
		if not ('---' in codons):
			rates.append(rawrates[j])
			ancestor += rawancestor[3*j:3*j+3]
			j += 1
		else:
			rates.append((-1,-1))
			ancestor += "---"

	if False:
		osremove(codeml_infile)
		osremove(codeml_file)
		osremove(tempfilename)
	return (rates, ancestor)

def Get_Site_Rates_Ancestors_And_Probabilities(seqs, seq_labels=None, tree_string=None ):
	"""Reconstructs an ancestor and per-site substitutions from an alignment.

	'seqs' is a list of aligned sequences.  'tree_string' is a phylogeny in the format needed by PAML,
	e.g. '((seq<n>,seq<n+1>),(seq<n+2>,seq<n+3>))'
	Uses the PAML program codeml to reconstruct the ancestral sequence by maximum likelihood
	and returns the ancestor and the per-site total synonymous and nonsynonymous
	changes in the tree, in a list of (syn,nsyn) pairs, according to the Suzuki and Gojobori MBE 16(10) (1999) method."""
	#if not seq_type in ['protein', 'codon']:
	#	raise PAMLError, "Distance type of %s is invalid." % seq_type
	seq_type = "codon"
	if len(seqs)<2:
		ancestor = seqs[0]
		rates = [(0.0,0.0)]*len(seqs[0])
		return (rates, ancestor)

	cm = CodeML(seq_type, control_file = "codeml.ctl")
	codeml_file = 'get-rate-ancestor-tmp.txt'
	codeml_infile = 'get-rate-ancestor-infile-tmp.txt'
	tree_file = 'get-rate-ancestor-tree-tmp.txt'
	if not seq_labels:
		seq_labels = ['seq_%d' % (i+1,) for i in range(len(seqs))]
	if not tree_string:
		tree_string = "(" + ','.join(seq_labels) + ")"
	# Remove gaps in aligned seqs before reconstruction
	n_seqs = len(seqs)
	new_seqs = ['']*n_seqs
	ngaps = 0
	for i in range(len(seqs[0])/3):
		codons = [s[3*i:3*i+3] for s in seqs]
		if not ('---' in codons):
			for j in range(n_seqs):
				new_seqs[j] += codons[j]
		else:
			ngaps += 1
	# Perhaps there are no remaining residues.  PAML chokes if this is so.
	if len(new_seqs[0]) == 0:
		raise PAMLError, "No aligned residues without any gaps in input sequences!"
 	writeMultipleTree(tree_file, new_seqs, tree_string)
 	writeMultipleSequences(codeml_infile, new_seqs, seq_labels)
 	option_dict = dict([("model","0"), ("NSsites","0"), ("runmode","0"),\
						("RateAncestor","1"),("Small_Diff","1e-6"),("verbose","2"),("seqtype","1")])
	#tempfilename = cm.Run(codeml_infile, codeml_file, tree_file, option_dict)
	(rawrates, rawsitetypes, rawcounts, ancestor_ids, rawancestors, rawnodeprobs, reconstructed_tree) = cm.Get_Site_Rates_Ancestors_And_Probabilities()

	# Add the gaps back in
	rates = []
	ancestors = ['']*len(rawancestors)
	# Probability of each amino acid at each node of the reconstructed tree
	nodeprobs = [[] for ps in rawnodeprobs]
	# Number of synonymous and nonsynonymous sites
	sitetypes = []
	counts = []
	j = 0
	for i in range(len(seqs[0])/3):
		codons = [s[3*i:3*i+3] for s in seqs]
		if not ('---' in codons):
			# Add in data
			rates.append(rawrates[j])
			sitetypes.append(rawsitetypes[j])
			counts.append(rawcounts[j])
			for ai in range(len(rawancestors)):
				ancestors[ai] += rawancestors[ai][3*j:3*j+3]
			for ni in range(len(rawnodeprobs)):
				nodeprobs[ni] += [rawnodeprobs[ni][j]]
			j += 1
		else:
			# Add in gaps
			rates.append((-1,-1))
			sitetypes.append((-1,-1))
			counts.append((-1,-1))
			for ai in range(len(rawancestors)):
				ancestors[ai] += '---'
			for ni in range(len(rawnodeprobs)):
				nodeprobs[ni] += [-1]

	if False:
		osremove(codeml_infile)
		osremove(codeml_file)
		osremove(tempfilename)
	return (rates, sitetypes, counts, ancestor_ids, ancestors, nodeprobs, reconstructed_tree)

def makeCumulativeProbabilities(abs_probs):
	psum = sum([p for (d,p) in abs_probs])
	# Normalize sums to ensure that sum is 1.0
	new_probs = [(p/psum, d) for (d,p) in abs_probs]
	psum2 = sum([p for (p,d) in new_probs])
	assert abs(psum2-1) < 1e-5, "psum = %1.6f" % psum2

	new_probs
	new_probs.sort()
	new_probs.reverse()
	cum_p = 0.0
	for i in range(len(new_probs)):
		(prob,data) = new_probs[i]
		# Take care of inevitable rounding errors, and limit list size: if
		# probabilities have gone to zero, don't add them, and set cumulative
		# to 1.0.
		if prob == 0.0:
			new_probs[i-1] = (1.0, data)
			return new_probs[0:i]
		cum_p += prob
		# turn into cumulative prob
		new_probs[i] = (cum_p, data)
		if cum_p == 1.0:
			return new_probs[0:i+1]
	# If all probabilities are nonzero, still enforce final cumulative
	# probability of 1.0
	new_probs[i] = (1.0, data)
	return new_probs

def sampleCumulative(cum_probs):
	rand = random.random()
	for (p,codon) in cum_probs:
		if rand <= p:
			return codon

def sampleAncestor(cumprob_list):
	seq = ''
	for site in range(len(cumprob_list)):
		codon = sampleCumulative(cumprob_list[site])
		seq += codon
	return seq

def Get_dNdS_Per_Codon(seqs, seq_labels, tree_string, options=[], num_samples=100):
	debugging = False
	for seq in seqs:
		if len(seq)%3 != 0:
			raise PAMLError, "Sequence length %d is not a multiple of three -- only codon sequences allowed for per-site rates" % len(seq)

	seq_type = "codon"
	if len(seqs)<2:
		ancestor = seqs[0]
		dn_rates = [(0.0)]*len(seqs[0])
		ds_rates = [(0.0)]*len(seqs[0])
		return dn_rates, ds_rates

	cm = CodeML(seq_type, control_file = "codeml.ctl")
	codeml_file = 'get-rate-ancestor-tmp.txt'
	codeml_infile = 'get-rate-ancestor-infile-tmp.txt'
	tree_file = 'get-rate-ancestor-tree-tmp.txt'
	if not seq_labels:
		seq_labels = ['seq_%d' % (i+1,) for i in range(len(seqs))]
	if not tree_string:
		tree_string = "(" + ','.join(seq_labels) + ")"
	# Remove gaps in aligned seqs before reconstruction
	n_seqs = len(seqs)
	new_seqs = ['']*n_seqs
	ngaps = 0
	for i in range(len(seqs[0])/3):
		codons = [s[3*i:3*i+3] for s in seqs]
		if not ('---' in codons):
			for j in range(n_seqs):
				new_seqs[j] += codons[j]
		else:
			ngaps += 1

	for seq in new_seqs:
		assert "---" not in seq
	# Perhaps there are no remaining residues.  PAML chokes if this is so.
	if len(new_seqs[0]) == 0:
		raise PAMLError, "No aligned residues without any gaps in input sequences!"
 	writeMultipleTree(tree_file, new_seqs, tree_string)
 	writeMultipleSequences(codeml_infile, new_seqs, seq_labels)
	#("fix_kappa","1"),("kappa","0.1"),
 	option_dict = dict([("model","0"), ("NSsites","0"), ("runmode","0"), \
						("RateAncestor","1"),("Small_Diff","1e-6"),("verbose","2"),("seqtype","1")]+options)
	if debugging:
		print "# Running PAML"
	tempfilename = cm.Run(codeml_infile, codeml_file, tree_file, option_dict)
	if debugging:
		print "# Fetching results"
	(rawrates, rawsitetypes, rawcounts, ancestor_ids, rawancestors, rawnodeprobs_list, reconstructed_tree) = cm.Get_Site_Rates_Ancestors_And_Probabilities()
	if debugging:
		print "# Starting sampling"
	full_tree = newick.tree.parseTree(reconstructed_tree)
	full_tree_node_dict = dict([(n.name,n) for n in full_tree.tree_nodes()])
	# Sample sequences from the node probabilities
	dn_ds_per_site = []
	prot_length = len(new_seqs[0])/3
	assert len(rawrates) == prot_length
	assert len(rawsitetypes) == prot_length

	# Convert node probabilities -- which are, for each site, (codon, P) pairs -- into easily sampleable
	# lists, namely cumulative probabilities
	cumprobs_list = []
	for node_num in range(len(rawnodeprobs_list)):
		nprobs = rawnodeprobs_list[node_num]
		cumprobs_list.append([])
		for site in range(prot_length):
			if False and debugging:
				print "nprobs[",site,"]", nprobs[site]
				print "abs:", ancestor_ids[node_num], site, ["%s:%1.3f" % (c,p) for (c,p) in nprobs[site] if p>0]
			cum_prob = makeCumulativeProbabilities(nprobs[site])
			if debugging and False:
				print "cum:", ancestor_ids[node_num], site, ["%s:%1.3f" % (c,p) for (p,c) in cum_prob if p>0]
			cumprobs_list[node_num].append(cum_prob)
	for snum in range(num_samples):
		if debugging:
			print "# Sample", snum
		# Sample sequences for nodes and make ancestors
		sample_ancestors = []
		i = 0
		for nprobs in cumprobs_list:
			s_anc = sampleAncestor(nprobs)
			if debugging:
				print ancestor_ids[i], s_anc, rawancestors[i]
			i += 1
			sample_ancestors.append(s_anc)
		# Assemble sequences into tree
		all_seqs = new_seqs + sample_ancestors
		all_ids = seq_labels + ancestor_ids
		# full_tree_node_dict contain
		for (id, seq) in zip(all_ids, all_seqs):
			full_tree_node_dict[id].sequence = seq
		dn_ds_per_site.append([])
		for i in range( prot_length ):
			( Dn, Ds ) = treeDnDsAtCodon( full_tree, 3*i )
			dn_ds_per_site[snum].append((Dn,Ds))
	# Average over samples.
	avg_dns = [0.0]*prot_length
	avg_dss = [0.0]*prot_length
	for site in range(prot_length):
		dns = []
		dss = []
		for snum in range(num_samples):
			(dn,ds) = dn_ds_per_site[snum][site]
			dns.append(dn)
			dss.append(ds)
		avg_dns[site] = sum(dns)/float(len(dns))
		avg_dss[site] = sum(dss)/float(len(dss))

	# Add the gaps back in
	gap_avg_dns = [0.0]*(len(seqs[0])/3)
	gap_avg_dss = [0.0]*(len(seqs[0])/3)
	j = 0
	for i in range(len(seqs[0])/3):
		codons = [s[3*i:3*i+3] for s in seqs]
		if not ('---' in codons):
			# Add in data
			gap_avg_dns[i] = avg_dns[j]
			gap_avg_dss[i] = avg_dss[j]
			j += 1
		else:
			# Add in gaps
			gap_avg_dns[i] = None
			gap_avg_dss[i] = None
	assert len(gap_avg_dns) == len(seqs[0])/3
	return gap_avg_dns, gap_avg_dss

def osremove(fname):
	try:
		os.remove(fname)
	except OSError, ose:
		pass

#----------------------------------------------------------------------------------
class CodeML:
	"Class for running the PAML program codeml."
	#-----------------------------------------------------------------------
	def __init__(self, seq_type, dir = '~/', bindir="bin", prog = 'codeml', \
		control_file = 'codeml.ctl', control_aa = 'codeml.ctl', \
		control_codon = 'codeml.ctl'):
	  	"""Setup to run the PAML program CodeML for amino acids or nucleotides.

		'type' specifies if we are looking at amino acids (type = 'protein') or
		nucleotide codon distances (type = 'codon').
		'dir' is the directory where the program and control files are located.
		'prog' is the name of the program in 'dir'.
	   	'control_file' is the name of the control file for the program.
		'control_aa' is the template control file for all runs looking at protein
		sequences.
		'control_codon' is the template control file for all runs looking at codon
		acid sequences."""
		self.dir = os.path.expanduser(dir)
		assert os.path.isdir(self.dir), "Error, directory " + self.dir + " does not exist."
		if seq_type == 'protein':
			self.seq_type = 'protein'
			self.control_aa = os.path.join(self.dir,control_aa)
			assert os.path.isfile(self.control_aa), "Control file " + self.control_aa + " does not exist."
		elif seq_type == 'codon':
			self.seq_type = 'codon'
			self.control_codon = os.path.join(self.dir,control_codon)
			assert os.path.isfile(self.control_codon), "Control file " + self.control_codon + " does not exist."
		else:
			raise PAMLError, "Type of %s is invalid." % seq_type
		self.prog = os.path.join(self.dir,bindir,prog)
		assert os.path.isfile(self.prog)
		self.control_file = control_file

	#-----------------------------------------------------------------------
	def Run(self, infile, outfile, treefile = "codeml-tree-tmp.txt", option_dict=None):
		"""Runs PAML program CodeML.
		'infile' is the input file in Phylip format
		'treefile' is the name of the tree file.
		'outfile' is the name of the output file.
		"""

		# Read in lines from the template control file
		if self.seq_type == 'protein':
			controlfile = open(self.control_aa, 'r')
		elif self.seq_type == 'codon':
			controlfile = open(self.control_codon, 'r')
		else:
			raise PAMLError, "No valid type set: %s." % self.seq_type
		val_dict = ReadPAMLControlFile(controlfile)
		controlfile.close()
		# Write to the control file
		controlfile = file(self.control_file, 'w')
		#print val_dict
		for (k,v) in val_dict.items():
			if k == 'runmode':
				v = '-2'  # Compute pairwise maximum-likelihood values of dN and dS
				#v = '0'  # User tree
			elif k == 'noisy':
				#v = '0'  # Minimize output
				v = '9'  # Maximize output
			elif k == 'verbose':
				v = '0'  # Minimize output
			elif k == 'CodonFreq':
				v = '2' # F3X4, 9 free parameters
			elif k == 'estFreq':
				v = '0' # Compute codon frequencies from the data; no ML estimation
			elif k == 'NSsites':
				v = '0' # One w (dN/dS) for whole tree
			elif k == 'model':
				v = '0' #

		# Set options from user
		if option_dict:
			for (k,v) in option_dict.items():
				val_dict[k] = v
		# Set files
		val_dict['seqfile'] = infile
		val_dict['treefile'] = treefile
		val_dict['outfile'] = outfile

		# Write out new control file
		for (k,v) in val_dict.items():
			controlfile.write('%s = %s\n' % (k,v))
		controlfile.close()
		# run program
		tempfilename = 'codeml-outfile-tmp.txt'
		cmd = "%s %s > %s" % (self.prog, self.control_file, tempfilename)
		code = os.system(cmd)
		#code = os.spawnv(os.P_WAIT, self.prog, [x for x in cmd.split()])
		if code != 0:
			raise PAMLError, "Error running PAML; code was %s." % code
		return tempfilename

	def getKappa(self, out_filename):
		kappa = -1.0
		for line in file(out_filename,'r').readlines():
			if "kappa (ts/tv) =" in line: # found our line
				kappa = float(line.strip().split("=")[1])
		return kappa
	#-----------------------------------------------------------------------
	# (rawrates, rawsitetypes, rawcounts, rawancestors, rawnodeprobs, reconstructed_tree) = cm.Get_Site_Rates_Ancestors_And_Probabilities()
	def Get_Site_Rates_Ancestors_And_Probabilities(self):
		debugging = True
		rates = []
		counts = []
		sitetypes = []
		ancestors = {}
		nodeprobs = {} # Dictionary keyed by ancestral sequence id, yielding a list (len = # of sites) of probability lists (len = 61 sense codons)

		# Read lines from rst file
		rst_file = file('rst','r')
		found = False
		rst_lines = rst_file.readlines()
		rst_file.close()
		# Get reconstructed tree
		index_tree_string = ""
		master_tree_ind = -1
		index_tree_ind = -1
		internal_ind = -1
		for li in range(len(rst_lines)):
			line = rst_lines[li]
			if master_tree_ind>0:
				reconstructed_tree = rst_lines[master_tree_ind].strip()
				index_tree_string = rst_lines[index_tree_ind].strip()
				internal_nodes = rst_lines[internal_ind].strip().split()
				break
			if "Ancestral reconstruction by CODONML." in line:
				master_tree_ind = li+2
				index_tree_ind = li+4
				internal_ind = li+6
		# Process tree to extract names and numbers of nodes
		tree = newick.tree.parseTree(reconstructed_tree)
		#node_names = [node.name for node in tree.tree_nodes()]

		index_tree = newick.tree.parseTree(index_tree_string)
		internal_node_kid_dict = dict([(x.split("..")[1], x.split("..")[0]) for x in internal_nodes])

		all_nodes_labeled = False
		while not all_nodes_labeled:
			for node in index_tree.tree_nodes():
				if node.name: # If not None
					continue
				kid_names = [k.name for k in node.children()]
				for k in kid_names:
					try:
						par = internal_node_kid_dict[k]
						node.name = par
						break
					except KeyError:
						continue
			an = [x for x in index_tree.tree_nodes() if not x.name]
			all_nodes_labeled = (len(an) == 0)
		ind_node_names = [node.name for node in index_tree.tree_nodes()]

		internal_nodes = []
		for (n1, n2) in zip(tree.tree_nodes(), index_tree.tree_nodes()):
			if not n1.name:
				n1.name = n2.name
				internal_nodes.append(n2.name)
		node_names = [node.name for node in tree.tree_nodes()]
		reconstructed_tree = tree.__str__()
		num_sequences = len(node_names)
		num_internal_nodes = len(internal_nodes)
		#print node_names

		# Get node probs from rst
		for node_id in internal_nodes:
			nodeprobs[node_id] = []
			target_line = "Prob distribution at node %s, by site" % node_id
			li = 0
			done = False
			while not done:
				line = rst_lines[li]
				if target_line in line:
					li += 4
					line = rst_lines[li]
					while line.strip() != '':
						rawprobs = line.strip().split(":")[1].split()
						probs = [(x.split("(")[0], float(x.split("(")[1].split(")")[0])) for x in rawprobs]
						nodeprobs[node_id].append(probs)
						counter = len(nodeprobs[node_id])
						li += 1
						line = rst_lines[li]
					done = True
				else:
					li += 1

		#for k in nodeprobs.keys():
		#	print k, len(nodeprobs[k]), len(nodeprobs[k][0])
		# Get ancestors from rst
		li = 0
		anc_ind = 1
		while li < len(rst_lines):
			line = rst_lines[li]
			found = ("List of extant and reconstructed sequences" in line)
			if found:
				li += 4
				internal_node_names = nodeprobs.keys()
				internal_node_names.sort()
				ancestor_lines = rst_lines[li:li+num_sequences]
				break
			li += 1
		recon_seqs = dict([(l.split("#")[1].split()[0], ''.join(l.split()[2:])) for l in ancestor_lines[num_sequences-num_internal_nodes:]])
		for id in internal_node_names:
			ancestors[id] = recon_seqs[id]
		# Get rates from rst1
		rst1_file = file('rst1','r')
        # Line format:
		# 129 S N:   0.359  2.641 Sd Nd:    0.2   2.8
		for line in rst1_file.readlines()[4:]:
			if line.strip() == '':
				break
			flds = line.strip().split()
			ssites = float(flds[3])
			nsites = float(flds[4])
			sitetypes.append((ssites, nsites))
			sdiffs = float(flds[7])
			ndiffs = float(flds[8])
			counts.append((sdiffs,ndiffs))
			srate = 0.0
			if ssites > 0:
				srate = sdiffs/ssites
			nrate = 0.0
			if nsites > 0:
				nrate = ndiffs/nsites
			rates.append((srate, nrate))
		rst1_file.close()

		# Prepare final data
		ancestor_ids = internal_node_names
		ancestors = [ancestors[id] for id in ancestor_ids]
		nodeprobs_list = [nodeprobs[id] for id in ancestor_ids]
		return (rates, sitetypes, counts, ancestor_ids, ancestors, nodeprobs_list, reconstructed_tree)
	#-----------------------------------------------------------------------
	def Get_Site_Rates_And_Ancestor(self):
		rates = []
		ancestor = ''
		# Get ancestor from rst
		rst_file = file('rst','r')
		found = False
		line = rst_file.readline()
		while line:
			if not found:
				found = ("List of extant and reconstructed sequences" in line)
				if found:
					# Read three irrelevant lines
					line = rst_file.readline()
					line = rst_file.readline()
					line = rst_file.readline()
			else:
				# Find last reconstructed sequence.
				if line.strip()=='':  # Blank line means we're done
					break
				flds = line.strip().split()
				# First two fields are sequence identifier (e.g. "node #3"); next N-2 are reconstructed sequence
				ancestor = ''.join(flds[2:])
			line = rst_file.readline()
		rst_file.close()

		# Get rates from rst1
		rst1_file = file('rst1','r')
		found = False
		line = rst1_file.readline()
        # Line format:
		# 129 S N:   0.359  2.641 Sd Nd:    0.2   2.8
		for line in rst1_file.readlines()[3:]:
			if line.strip() == '':
				break
			flds = line.strip().split()
			ssites = float(flds[3])
			nsites = float(flds[4])
			sdiffs = float(flds[7])
			ndiffs = float(flds[8])
			srate = 0.0
			if ssites > 0:
				srate = sdiffs/ssites
			nrate = 0.0
			if nsites > 0:
				nrate = ndiffs/nsites
			rates.append((srate, nrate))
		rst1_file.close()
		return rates, ancestor
	#-----------------------------------------------------------------------
	def Get_NS_Sites(self, tempfile):
		for line in file(tempfile, 'r').readlines():
			flds = line.split()
			if len(flds) < 3:
				continue
			if flds[1] == '1:Sites':
				line_end = ''.join(flds[2:])
				lhs = line_end.split('=')[0]
				(syn, ns) = lhs.split('+')
				syn = float(syn)
				ns = float(ns)
		return syn, ns
	#-----------------------------------------------------------------------
	def Get_Dist(self):
		"""Returns the distance(s) between two sequences after a run of CodeML.

		If the run was with type 'protein', returns a single scalar representing
		the distance between the amino acid sequences.  If the run was with
		type 'codon', returns a 6-tuple, with the entries as:
		'(maximum likelihood (ML) distance, ML synonymous rate,
		ML nonsynonymous rate, Nei and Gojobori 1986 (NG) distance,
		NG synonymous rate, NG nonsynonymous rate)'
		Removes the distance files after they are read.
		Returns 'None' if there is a problem."""
		if self.seq_type == 'protein':
			if not os.path.isfile('2AA.t'):
				raise PAMLError, "Cannot find 2AA.t."
			else:
				file = open('2AA.t', 'r')
				lines = file.readlines()
				file.close()
				assert lines[0].strip() == '2'
				assert len(lines) == 3, "Error, 2AA.t did not have 3 lines."
				entries = lines[2].split()
				assert len(entries) >= 2
				osremove('2AA.t')
				return float(entries[len(entries) - 1]) # return protein distance
		elif self.seq_type == 'codon':
			distances = []
			file_list = ['2ML.t', '2ML.dS', '2ML.dN', '2NG.t', '2NG.dS', '2NG.dN']
			for file_name in file_list:
				if not os.path.isfile(file_name):
					raise PAMLError, "Cannot find distance file %s." % file_name
				else:
					file = open(file_name, 'r')
					lines = file.readlines()
					file.close()
					assert lines[0].strip() == '2'
					assert len(lines) == 3
					entries = lines[2].split()
					assert len(entries) >= 2
					# get the distance
					try:
						dist = float(entries[len(entries) - 1])
					except ValueError:
						dist = -1
					distances.append(dist)
					osremove(file_name)
			assert len(distances) == 6
			return tuple(distances)
		else:
			raise PAMLError, "Type of %s is invalid." % self.seq_type
	#---------------------------------------------------------------------
	def Get_Multiple_Dist(self):
		"""Returns the distance(s) between two sequences after a run of CodeML.

		If the run was with type 'protein', returns a single scalar representing
		the distance between the amino acid sequences.  If the run was with
		type 'codon', returns a 3-tuple, with the entries as:
		'(maximum likelihood (ML) distance, ML synonomous rate,
		ML non-synonomous rate, Nei and Gojobori 1986 (NG) distance,
		NG synonomous rate, NG non-nysynonomous rate)'
		Removes the distance files after they are read.
		Returns 'None' if there is a problem."""
		if self.seq_type == 'protein':
			raise PAMLError, "Protein distances not supported for multiple sequences"
		elif self.seq_type == 'codon':
			distances = []
			file_list = ['rst1']
			for file_name in file_list:
				if not os.path.isfile(file_name):
					raise PAMLError, "Cannot find distance file %s." % file_name
				else:
					file = open(file_name, 'r')
					lines = file.readlines()
					file.close()
					flds = lines[0].strip().split()
					num_branches = (len(flds)-3)/3
					ds_sum = 0.0
					dn_sum = 0.0
					dd_sum = 0.0
					for i in range(num_branches+3, len(flds)-1,2):
						ds_sum += float(flds[i])
					for i in range(num_branches+3+1, len(flds),2):
						dn_sum += float(flds[i])
					for i in range(num_branches):
						dd_sum += float(flds[i])
					# get the distance
					#print dd_sum, ds_sum, dn_sum
			#		osremove(file_name)
			return dd_sum, ds_sum, dn_sum
		else:
			raise PAMLError, "Type of %s is invalid." % self.seq_type

	#---------------------------------------------------------------------
	def Get_Multiple_Dist_Physical(self, tempfilename):
		return self.getMultipleDistPhysical(tempfilename)

	def getMultipleDistPhysical(self, tempfilename):
		"""Returns the distance(s) between two sequences after a run of CodeML.

		If the run was with type 'protein', returns a single scalar representing
		the distance between the amino acid sequences.  If the run was with
		type 'codon', returns a 6-tuple, with the entries as:
		'(maximum likelihood (ML) distance, ML synonomous rate,
		ML non-synonomous rate, Nei and Gojobori 1986 (NG) distance,
		NG synonomous rate, NG non-nysynonomous rate)'
		Removes the distance files after they are read.
		Returns 'None' if there is a problem."""

		debugging = False
		if self.seq_type == 'protein':
			raise PAMLError, "Protein distances not supported for multiple sequences"
		elif self.seq_type == 'codon':
			distances = []

			if not os.path.isfile(tempfilename):
				raise PAMLError, "Cannot find distance-containing temporary file %s." % tempfilename
			else:
				f = file(tempfilename, 'r')
				lines = f.readlines()
				f.close()
				# For NSsites = 0, there will be one line of distances
				# For NSsites > 0, there will be multiple distances which must be weighted
				#    by the proportion of sites with a given w (=dN/dS).
				site_rates = []
				for i in range(len(lines)):
					if 'four-fold sites)' in lines[i]:
						# Target is line i+1
						line = lines[i+1]
						assert "dN*" in line and "dS*" in line
						#             dN*=  0.00736 dS*=  0.15299 S* = 200.65 N* = 657.35
						flds = line.split('=')
						dn = float(flds[1].split()[0])
						ds = float(flds[2].split()[0])
						nsyn = float(flds[3].split()[0])
						nns = float(flds[4].split()[0])
						site_rates.append((dn,ds,nns,nsyn))
				if len(site_rates) == 1:
					# Only one set of rates -- just return
					return (dn, ds, nns, nsyn)
				# Multiple site types: weight the results
				#(probs, omegas) = self.Get_Site_Proportions()
				probs = self.Get_Site_Proportions()
				# No (or not enough) probabilities; just average.
				if len(probs) < len(site_rates):
					n_rates = len(site_rates)
					avg_dn = sum([dn for (dn,ds,nns,nsyn) in site_rates])/n_rates
					avg_ds = sum([ds for (dn,ds,nns,nsyn) in site_rates])/n_rates
					avg_nns = sum([nns for (dn,ds,nns,nsyn) in site_rates])/n_rates
					avg_nsyn = sum([nsyn for (dn,ds,nns,nsyn) in site_rates])/n_rates
					return (avg_dn, avg_ds, avg_nns, avg_nsyn)
				(dn,ds,nns,nsyn) = (0,0,0,0)
				for xi in range(len(probs)):
					p = probs[xi]
					#w = omegas[xi]
					(xdn,xds,xnns,xnsyn) = site_rates[xi]
					if debugging:
						print xdn,xds,xnns,xnsyn, p
					dn += p*xdn
					ds += p*xds
					nns += p*xnns
					nsyn += p*xnsyn
					#osremove(file_name)
			return (dn, ds, nns, nsyn)
		else:
			raise PAMLError, "PAML sequence type of %s is invalid." % self.seq_type

	#---------------------------------------------------------------------
	def getPhysicalDistRST(self, rstfilename):
		"""Returns the physical-sites distances between two sequences after a run of CodeML.

		If the run was with type 'protein', returns a single scalar representing
		the distance between the amino acid sequences.  If the run was with
		type 'codon', returns a 6-tuple, with the entries as:
		'(maximum likelihood (ML) distance, ML synonomous rate,
		ML non-synonomous rate, Nei and Gojobori 1986 (NG) distance,
		NG synonomous rate, NG non-nysynonomous rate)'
		Returns 'None' if there is a problem."""

		debugging = False
		if self.seq_type == 'protein':
			raise PAMLError, "Protein distances not supported for multiple sequences"
		elif self.seq_type == 'codon':
			distances = []

			if not os.path.isfile(rstfilename):
				raise PAMLError, "Cannot find distance-containing temporary file %s." % rstfilename
			else:
				lines = file(rstfilename, 'r').readlines()
				# For
				flds = lines[0].strip().split("\t")
				dn = float(flds[-4])
				ds = float(flds[-3])
			return (dn,ds)
		else:
			raise PAMLError, "PAML sequence type of %s is invalid." % self.seq_type

	#------------------------------------------------------------------------
	def Get_Site_Proportions(self):
		"""Retrieves posterior probabilities of each dN/dS class after a CodeML run.

		"""
		if not os.path.isfile("rst"):
			raise PAMLError, "Cannot find distance-class-probability-containing temporary file %s." % "rst"
		else:
			f = file("rst",'r')
			lines = f.readlines()
			f.close()
			probs = []
			for li in range(len(lines)):
				line = lines[li]
				if "dN/dS for site class" in line:
					n_classes = int(line.strip().split("=")[1].split(")")[0])
					line = lines[li+2]
					probs = [float(x) for x in line.strip().split()[1:]]
					#line = lines[li+3]
					#omegas = [float(x) for x in line.strip().split()[1:]]
					#assert len(omegas) == n_classes
					assert len(probs) == n_classes
		return probs

#-------------------------------------------------------------------------------
def Write_Phylip(filename, seq1, seq2):
	"""Writes the two sequences 'seq1' and 'seq2' to the file 'filename'.

	'seq1' and 'seq2' must be of the same length.
	'filename' is created if it does not exist, and overwritten if it exists.
	The two sequences are given the generic names 'seq_1' and 'seq_2'."""
	assert len(seq1) == len(seq2)
	length = len(seq1)
	file = open(filename, 'w')
	file.write("2 %d\n" % len(seq1))
	file.write("seq_1\n%s\nseq_2\n%s\n" % (seq1, seq2))
	file.close()
#-------------------------------------------------------------------------------
def Write_Sequence_Pair(filename, seq1, seq2):
	"""Writes the two sequences 'seq1' and 'seq2' to the file 'filename'.

	'seq1' and 'seq2' must be of the same length.
	'filename' is created if it does not exist, and overwritten if it exists.
	The two sequences are given the generic names 'seq_1' and 'seq_2'."""
	assert len(seq1) == len(seq2)
	length = len(seq1)
	file = open(filename, 'w')
	file.write("2 %d\n" % len(seq1))
	file.write("seq_1\n%s\nseq_2\n%s\n" % (seq1, seq2))
	file.close()
#-------------------------------------------------------------------------------
def writeMultiplePhylip(filename, seqs):
	"""Writes the sequences 'seqs' to the file 'filename'.

	'filename' is created if it does not exist, and overwritten if it exists.
	Sequences are given the generic names 'seq_1', 'seq_2'... if labels are
	not provided."""
	length = len(seqs[0])
	file = open(filename, 'w')
	file.write("%d %d\n" % (len(seqs), length))
	for i in range(len(seqs)):
		file.write("seq_%d\n%s\n" % (i+1, seqs[i]))
	file.close()
#-------------------------------------------------------------------------------
def writeMultipleSequences(filename, seqs, labels=None):
	"""Writes the sequences 'seqs' to the file 'filename'.

	'filename' is created if it does not exist, and overwritten if it exists.
	Sequences are given the generic names 'seq_1', 'seq_2'... if labels are
	not provided."""
	if not labels:
		labels = ['seq_%d' % (i+1,) for i in range(len(seqs))]

	length = len(seqs[0])
	file = open(filename, 'w')
	file.write("%d %d\n" % (len(seqs), length))
	for i in range(len(seqs)):
		file.write("%s\n%s\n" % (labels[i], seqs[i]))
	file.close()
#-------------------------------------------------------------------------------
def writeTree(filename, seq1, seq2):
	"""Writes tree for two sequences 'seq1' and 'seq2' to the file 'filename'.

	'filename' is created if it does not exist, and overwritten if it exists.
	The two sequences are given the generic names 'seq_1' and 'seq_2'."""
	assert len(seq1) == len(seq2)
	length = len(seq1)
	file = open(filename, 'w')
	file.write("2 1\n") # 2 species, 1 tree
	file.write("(seq_1,seq_2);\n")
	file.close()
#------------------------------------------------------------------------------
def writeMultipleTree(filename, seqs, tree_string):
	"""Writes tree for seqs into the file 'filename'.

	'filename' is created if it does not exist, and overwritten if it exists."""
	# Warning: PAML will silently fail if the tree string is not terminated with a semicolon ';'
	if not tree_string[-1] == ';':
		tree_string += ';'
	num_seqs = len(seqs)
	file = open(filename, 'w')
	file.write("%d 1\n" % num_seqs) # n species, 1 tree
	file.write("%s\n" % tree_string)
	file.close()
#------------------------------------------------------------------------------
def ReadPAMLControlFile(f):
	lines = f.readlines()
	vals = {}
	for line in lines:
		if line.strip() == '':
			continue
		notcomment = line.split('*')[0].strip()
		flds = notcomment.split('=')
		if len(flds) < 2:
			continue
		vals[flds[0].strip()] = flds[1].strip()
	return vals
# End my_paml.py
