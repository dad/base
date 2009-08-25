#! python

"""Module for running PAML programs.

Original version by Jesse Bloom, 2004
Expanded, rewritten and maintained by D. Allan Drummond, 2004-2009.
treeDnDsAtCodon and Node class by Claus O. Wilke, 2006-2007.
"""

import shutil, sys, os, re, string, math, random
import codon, newick

#----------------------------------------------------------------------------------
class PAMLError(Exception):
	"""Error running or processing PAML."""


def getPairwiseRates(seq1, seq2, options={} ):
	"""Computes the evolutionary distance(s) between aligned sequences.
	Returns (dn, ds, nn, nsyn, kappa)
	"""
	cm = CodeML('codon', options)
	cm.loadSequences([seq1, seq2])
	cm.run()
	(dn_ml, ds_ml) = cm.getPairwiseRates()
	(nn, ns) = cm.getSites()
	kappa = cm.getKappa()
	return (dn_ml, ds_ml, nn, ns, kappa)

def getPairwisePhysicalRates(seq1, seq2, options={} ):
	"""Computes the evolutionary distance(s) between aligned sequences.
	Returns (dn, ds, nn, nsyn, kappa)
	"""
	cm = CodeML('codon', options)
	cm.loadSequences([seq1, seq2])
	cm.run()
	(dn, ds, nn, ns) = cm.getPairwisePhysicalRates()
	kappa = cm.getKappa()
	return (dn, ds, nn, ns, kappa)

def getTreeDistancePhysicalKappa(seqs, seq_labels=None, tree_string=None, options={} ):
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

	cm = CodeML('codon', options)
	cm.loadSequences(seqs, seq_labels, tree_string)
	cm.run()
   	(dn, ds, nn, ns) = cm.getMultipleDistPhysical()
   	#(dn, ds) = cm.getPhysicalDistRST()
	kappa = cm.getKappa()
	return (dn, ds, nn, ns, kappa)

def getPairwiseDistancePhysicalKappa(seqs, seq_labels=None, tree_string=None, options=[] ):
	"""Computes the evolutionary distance(s) between two aligned sequences.
	Uses physical definition of sites, rather than
	typical mutational-opportunity definition.  See Bierne and Eyre-Walker
	MBE 2003.

	'seqs' is a list of aligned sequences.  'tree' is a phylogeny in Newick format.
	Uses the PAML program CodeML to compute the distance(s) between the
    two sequences, and returns it.  The sequences must be DNA sequences aligned according to a protein
	alignment.
	Returns a five-tuple:
		(maximum likelihood (ML) nonsynonymous rate, ML synonymous rate,
		number of nonsynonymous sites, number of synonymous sites, transition-transversion rate ratio (kappa)"""

	options_dict = dict(options)
	# Ensure that we're not constraining kappa
	options_dict['fix_kappa'] = '0'
	options_dict['kappa'] = '3'
	# Pairwise comparison
	options_dict['runmode'] = '-2'

	cm = CodeML('codon', options_dict)
	#print cm.options
	#print options_dict
	cm.loadSequences(seqs, seq_labels, tree_string)
	cm.run()
   	(dn, ds, nn, ns) = cm.getMultipleDistPhysical()
	kappa = cm.getKappa()
	return (dn, ds, nn, ns, kappa)

def getCodonSelectionCoefficients(seqs, seq_labels=None, tree_string=None, options={}):
	"""Computes the selection coefficients for codon substitutions according to the model
	of Yang and Nielsen MBE 2008.

	'seqs' is a list of aligned sequences.  'tree' is a phylogeny in Newick format.
	Uses the PAML program CodeML to compute the distance(s) between the
    two sequences, and returns it.  The sequences must be DNA sequences aligned according to a protein
	alignment.

	Returns a list of objects containing the codon selection coefficients.
	"""

	options_dict = options
	# Ensure user tree; won't work with runmode=-2
	options_dict['runmode'] = '0'
	cm = CodeML('codon', options_dict)
	cm.loadSequences(seqs, seq_labels, tree_string)
	cm.run()
	sel_coef_dict = cm.getCodonSelectionCoefficients()
	return sel_coef_dict

# Decorate tree with branch rates.  Each non-root node gets a
# RateResult annotation showing the rates from that node to its parent.
def getBranchRates(seqs, seq_labels, tree_string, options={}):
	options_dict = options
	# Ensure user tree; won't work with runmode=-2
	# DAD: make this work with pairwise comparison...
	# DAD: this all should be independent of seqs...more like a helper function
	# to turn branch_rates into a more usable structure.
	# The goal is to be able to take any pair of nodes on the tree and
	# instantly look up the distance between them.
	options_dict['runmode'] = '0'
	cm = CodeML('codon', options_dict)
	tree = newick.tree.parseTree(tree_string)
	leaf_dict = dict([(x.name,x) for x in tree.leaves])
	assert set(leaf_dict.keys()) == set(seq_labels)
	cm.loadSequences(seqs, seq_labels=seq_labels, tree_string=tree_string)
	cm.run()
	branch_rates = cm.getBranchRates()
	# Now map branch rates onto the tree
	# With n = the number of sequences, the first 1:n node numbers are in the order passed in.
	# The remaining node numbers can be inferred from their branch IDs.
	# Each branch can be uniquely identified by its root and tip.
	# The branch IDs are integer..integer; map the integers to the corresponding subtree strings.
	id_node_map = {}
	for ni in range(len(seq_labels)):
		id_node_map["%d"%(ni+1,)] = leaf_dict[seq_labels[ni]]

	# Sort branch rates by their leaves
	branch_rate_list = []
	for (k,v) in branch_rates.items():
		bid = k.split('..')
		branch_rate_list.append((bid[1], (k,v)))
	branch_rate_list.sort()

	branch_rate_dict = {}
	for (k,v) in [y for (x,y) in branch_rate_list]:
		ids = k.split("..")
		# These are integer IDs
		child_id = ids[1]
		parent_id = ids[0]
		# Should already have the child in our name map
		child_node = id_node_map[child_id]
		if not id_node_map.has_key(parent_id):
			# Get the parent node
			parent_node = child_node.parent
			id_node_map[parent_id] = parent_node
		# Now assign the branch rate entry to the child
		child_node.branch_rate = v
	return tree


def osremove(fname):
	try:
		os.remove(fname)
	except OSError, ose:
		pass




##-------------------------------------------------------
# Classes for processing raw PAML output
##-------------------------------------------------------

class RateResult:
	"""Class for storing results of an evolutionary rate calculation"""
	dn = None
	ds = None
	num_syn_sites = None
	num_ns_sites = None
	omega = None
	kappa = None
	alpha = None
	tree_length = None

	def __init__(self):
		return


class CodonSelectionResult:
	"""Class for storing results of FMutSel estimation of codon-specific selection and mutation parameters"""
	codon_I = None      # codon I
	codon_J = None      # codon J
	FI_FJ = None        # 2Ns_IJ, the population-scaled fitness difference between codon I and J
	pMut_IJ = None      # probability of a mutation from codon I to J
	pSub_IJ = None      # probability of a substitution from codon I to J
	FJ_FI = None        # 2Ns_IJ, the population-scaled fitness difference between codon J and I
	pMut_JI = None      # probability of a mutation from codon J to I
	pSub_JI = None      # probability of a substitution from codon J to I

	def __init__(self):
		return

	def __str__(self):
		if self.codon_I < self.codon_J:
			line = "%s\t%s\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f" % \
				   (self.codon_I, self.codon_J, self.FI_FJ, self.pMut_IJ, self.pSub_IJ, self.FJ_FI, self.pMut_JI, self.pSub_JI)
		else:
			line = "%s\t%s\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f" % \
				   (self.codon_J, self.codon_I, self.FJ_FI, self.pMut_JI, self.pSub_JI, self.FI_FJ, self.pMut_IJ, self.pSub_IJ)
		return line

#----------------------------------------------------------------------------------
class CodeML:
	"Class for running the PAML program codeml."

	## Options
	# All of these options are used by the PAML control file.
	# Do not store non-control-file options here!
	default_options = {
		'runmode':'0',     # 0=user tree, 1=semi-automatic, 2=automatic, 3=StepwiseAddition, (4,5)=PerturbationNNI, -2=pairwise
		'noisy':'9',       # 9=maximum output
		'verbose':'0',     # 0=minimize screen output
		'seqtype':'1',     # 1=codons; 2=aas; 3=codons-->aas
		'CodonFreq':'7',   # 0=1/61 each, 1=F1X4, 2=F3X4 (9 free params), 3=codon table (61 free params), 7=FMutSel
		'estFreq':'0',     # 0=estimate codon freq. from data; 1=ML estimate of codon freqs
		'model':'0',       # 0=one omega for all branches
		'NSsites':'0',     # 0=one omega for all sites, 1=neutral, 2=selection, 3=discrete, 4=freqs, 5=gamma, 6=2gamma, 7=beta, 8=beta&w, 9=beta&gamma, 10=beta&gamma+1, 11=beta&normal>1, 12=0&2normal>1, 13=3normal>0
		'fix_alpha':'1',   # 0=estimate alpha; 1=use fixed value
		'alpha':'0',       # fixed or initial value for alpha; 0=infinity (constant rate)
		'fix_kappa':'0',   # 0=estimate kappa; 1=use fixed value
		'kappa':'3',       # fixed or initial value for alpha
		'fix_rho':'1',     # 0=estimate rho; 1=use fixed value
		'rho':'0.',        # fixed or initial value for rho=correlation parameter of the
							# auto-discrete-gamma model (site-to-site correlation); use fix_rho=1 and rho=0. to implement independent and constant site rates
		'icode':'0',       # 0=standard genetic code
		'fix_omega':'0',   # 0=estimate omega; 1=use fixed value
		'omega':'0.1',     # fixed or initial value for omega=dN/dS
		'RateAncestor':'0',# 0=don't estimate ancestral states; 1=estimate ancestral states
		'Small_Diff':'1e-6', #  a small value used in the difference approximation of derivatives, typically 1e-6
		'cleandata':'1',   # 1=remove ambiguous characters, gaps, etc. from sequences; 0=don't touch sequences
		'method':'0',      # 0=old method of simultaneously updating all parameters; 1=new method of iteratively updating branches
		'clock':'0',       # 0=no clock, rates free to vary from branch to branch
		'getSE':'0',	   # 0=don't get standard errors on estimates, 1=get standard errors
		# PAML file options
		'seqfile':  'codeml-seqfile-tmp.txt',
		'treefile': 'codeml-treefile-tmp.txt',
		'outfile':  'codeml-outfile-tmp.txt'
	}
	FMutSel_F_options = {'CodonFreq':'7', 'estFreq':'0', 'model':'0', 'runmode':'0', 'NSsites':'0'}
	FMutSel_options = {'CodonFreq':'7', 'estFreq':'1', 'model':'0', 'runmode':'0', 'NSsites':'0'}
	F3x4_options = {'CodonFreq':'2', 'model':'0', 'runmode':'0', 'NSsites':'0'}


	model_table = {"FMutSel-F":FMutSel_F_options, "FMutSel":FMutSel_options, "F3x4":F3x4_options}

	#-----------------------------------------------------------------------
	def __init__(self, seq_type='codon', options_dict={}, prog = 'codeml'):
	  	"""Setup to run the PAML program CodeML for amino acids or nucleotides.

		'type' specifies if we are looking at amino acids (type = 'protein') or
		nucleotide codon distances (type = 'codon').
		'prog' is the name of the program in 'dir'.
		"""
		if seq_type == 'protein':
			self.seq_type = 'protein'
		elif seq_type == 'codon':
			self.seq_type = 'codon'
		else:
			raise PAMLError, "Sequence type of %s is invalid; valid options are 'codon' or 'protein'." % self.seq_type
		self.prog = prog
		#if !os.path.isfile(self.prog):
		#	raise PAMLError, "Cannot find PAML program %s" % (self.prog,)

		self.options = self.default_options

		# Pass sequence type through to options
		if self.seq_type == 'protein':
			self.options['seqtype'] = '2'
		elif self.seq_type == 'codon':
			self.options['seqtype'] = '1'

		# Replace default options with user options
		for (k,v) in options_dict.items():
			self.options[k] = v

		# Files
		self.controlfile = 'codeml-controlfile-tmp.txt'
		self.tmpfile = 'codeml-tmpfile-tmp.txt'
		self.treefile = self.options['treefile']
		self.outfile = self.options['outfile']
		self.seqfile = self.options['seqfile']

		# Sequence information
		self.num_sequences = -1
		self.sequence_length = -1

	def getModelOptions(self, model_key):
		opts = self.default_options
		try:
			new_opts = self.model_table[model_key]
			for (k,v) in new_opts.items():
				opts[k] = v
		except KeyError:
			raise PAMLError, "Can't find model options for '%s'" % model_key
		return opts

	#-----------------------------------------------------------------------

	def loadSequences(self, seqs, seq_labels=None, tree_string=None):
		if not seq_labels:
			seq_labels = ['seq_%d' % (i+1,) for i in range(len(seqs))]
		if not tree_string:
			# This generates invalid (star) trees unless given only 2 sequences
			assert len(seqs) == 2, "Auto-generated tree, but there are %d sequences instead of 2" % len(seqs)
			tree_string = "(" + ','.join(seq_labels) + ")"

		writeMultipleTree(self.treefile, seqs, tree_string)
		writeMultipleSequences(self.seqfile, seqs, seq_labels)
		self.num_sequences = len(seqs)
		self.sequence_length = len(seqs[0])
		return seq_labels

	def run(self):
		"""Runs PAML program CodeML.
		"""

		# Create control file
		ctrl_file = file(self.controlfile,'w')

		# Write out new control file
		for (k,v) in self.options.items():
			ctrl_file.write('%s = %s\n' % (k,v))
		ctrl_file.close()
		# run program
		cmd = "%s %s > %s" % (self.prog, self.controlfile, self.tmpfile)
		code = os.system(cmd)
		#code = os.spawnv(os.P_WAIT, self.prog, [x for x in cmd.split()])
		if code != 0:
			raise PAMLError, "Error running PAML; code was %s." % code

	def getKappa(self):
		kappa = -1.0
		if self.options['fix_kappa'] == '1':
			# Kappa is fixed -- just return the value passed in.
			kappa = float(self.options['kappa'])
		else:
			if self.options['runmode'] == '-2':
				# rst line for fix_kappa=0, fix_alpha=1, fix_omega=0
				#seq seq        N       S       dN       dS     dN/dS   tree t   kappa     w       -lnL
				#  2   1   1789.7    670.3   0.0316   1.2116   0.0261   1.0593   2.1322   0.0261 -4213.851
				#
				# rst1 line for:
				#         N?      S?      ?       ?       ?       kappa  ?       t?       ?       w       -lnL
				# 2       820     199     0.5455  0.522   0.568   2.177   0.611   1.475   1.651   0.024   -4211.232
				lines = file('rst','r').readlines()
				i = 0
				found = False
				while not found and i < len(lines):
					if "pairwise comparison (Goldman & Yang 1994)" in lines[i]:
						params_line = lines[i+2].strip().split()
						assert self.options['fix_kappa'] == '0'
						assert self.options['fix_omega'] == '0'
						assert self.options['fix_alpha'] == '1'
						kappa = float(params_line[-3])
						found = True
					i += 1
			elif self.options['runmode'] == '0':
				lines = file(self.outfile,'r').readlines()
				for line in lines:
					if "kappa (ts/tv) =" in line:
						kappa = float(line.strip().split()[-1])
						break
			else:
				raise PAMLError, "Kappa extraction for runmode %s not implemented" % self.options['runmode']
		return kappa
	#-----------------------------------------------------------------------
	def getSites(self):
		for line in file(self.tmpfile, 'r').readlines():
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
	def putBranchRatesOnTree(self, seq_labels, rate_tree):
		# Returns a newick.tree with each node n an n.branch_rate entry of type RateResult
		# giving the evolutionary rate (and other data) separating node n from its parent.
		# The goal is to be able to take any pair of nodes on the tree and
		# instantly look up the distance between them.
		# seq_labels define the order in which the original sequences were passed in;
		# these will be sequentially numbered 1:n
		# rate_tree defines the tree, which must include all seq_labels as leaves and no others.
		leaf_dict = dict([(x.name,x) for x in rate_tree.leaves])
		assert set(leaf_dict.keys()) == set(seq_labels)
		branch_rates = self.getBranchRates()
		# Now map branch rates onto the tree
		# With n = the number of sequences, the first 1:n node numbers are in the order passed in.
		# The remaining node numbers can be inferred from their branch IDs.
		# Each branch can be uniquely identified by its root and tip.
		# The branch IDs are integer..integer; map the integers to the corresponding subtree strings.
		id_node_map = {}
		for ni in range(len(seq_labels)):
			id_node_map["%d"%(ni+1,)] = leaf_dict[seq_labels[ni]]

		# Sort branch rates by the numerical value of the second (child) integer identifier
		branch_rate_list = []
		for (branch_id, rate_result) in branch_rates.items():
			# branch ID is m..n where m > n
			child_id = branch_id.split('..')[1]
			branch_rate_list.append((child_id, (branch_id,rate_result)))
		branch_rate_list.sort()

		# Now map the branch rates onto the nodes of rate_tree
		for (k,v) in [y for (x,y) in branch_rate_list]:
			ids = k.split("..")
			# These are integer IDs
			child_id = ids[1]
			parent_id = ids[0]
			# Should already have the child in our name map
			child_node = id_node_map[child_id]
			if not id_node_map.has_key(parent_id):
				# Get the parent node
				parent_node = child_node.parent
				id_node_map[parent_id] = parent_node
				if parent_node.name is None:
					parent_node.name = parent_id
			# Now assign the branch rate entry to the child
			child_node.branch_rate = v

	#---------------------------------------------------------------------
	def getBranchRates(self):
		"""Returns branch lengths for a set of sequences after a run of CodeML.

		If the run was with type 'protein', returns a number representing
		the distance between the amino acid sequences.  If the run was with
		type 'codon', returns a 2-tuple, with the entries as:
		    (ML nonsynonymous rate, ML synonymous rate)
		"""
		if self.options['runmode'] == '-2':
			raise PAMLError, "CodeML.getBranchRates() not implemented for pairwise comparison (runmode=-2) yet."
		elif self.options['runmode'] == '0':
			if self.seq_type != "codon":
				raise PAMLError, "CodeML.getBranchRates() not implemented for seq_type = %s." % self.seq_type
			# Read the standard output file
			lines = file(self.outfile,'r').readlines()
			# Find the point at which branch values are enumerated
			for li in range(len(lines)):
				line = lines[li]
				# Looking for this line:  branch           t        N        S    dN/dS       dN       dS   N*dN   S*dS
				#print line.strip()
				if "branchtNSdN/dS" in ''.join(line.strip().split()):
					target_line = li + 2
			line = lines[target_line]
			li = target_line
			branch_entries = {}
			# Go until we hit a blank line
			while line.strip() != '':
				flds = line.strip().split()
				(branch_id, t, N, S, omega, dN, dS, NxdN, SxdS) = \
					(flds[0], float(flds[1]), float(flds[2]), float(flds[3]), float(flds[4]), float(flds[5]), float(flds[6]), float(flds[7]), float(flds[8]))
				res = RateResult()
				res.dn = dN
				res.ds = dS
				res.omega = omega
				res.tree_length = t
				res.num_syn_sites = S
				res.num_ns_sites = N
				branch_entries[branch_id] = res
				li = li+1
				line = lines[li]
			return branch_entries

	#---------------------------------------------------------------------
	def putAncestralSequencesOnTree(self, seq_labels, tree):
		anc_seq_dict = self.getAncestralSequences()
		node_dict = dict([(x.name, x) for x in tree.nodes])
		for sid in anc_seq_dict.keys():
			if node_dict.has_key(sid):
				node_dict[sid].sequence = anc_seq_dict[sid]
			else:
				# PAML lists ancestral nodes as "node #n" where we may have nodes named "n" in the tree.
				short_sid = sid[6:]
				try:
					node_dict[short_sid].sequence = anc_seq_dict[sid]
				except KeyError, ke:
					print "# Can't find a node named %s or %s in the tree %s" % (sid, short_sid, str(tree))

	#---------------------------------------------------------------------
	def getAncestralSequences(self):
		"""Returns ancestral sequences after a run of CodeML.

		If the run was with type 'protein', returns a number representing
		the distance between the amino acid sequences.  If the run was with
		type 'codon', returns a 2-tuple, with the entries as:
		    (ML nonsynonymous rate, ML synonymous rate)
		"""
		if self.options['runmode'] == '-2':
			raise PAMLError, "CodeML.getBranchRates() not implemented for pairwise comparison (runmode=-2) yet."
		elif self.options['runmode'] == '0':
			if self.seq_type != "codon":
				raise PAMLError, "CodeML.getBranchRates() not implemented for seq_type = %s." % self.seq_type
			if self.options["RateAncestor"] != "1":
				raise PAMLError, "Can't retrieve ancestral sequences because RateAncestor is set to %s instead of 1" % self.options["RateAncestor"]
			# Read the rst file
			lines = file("rst",'r').readlines()
			# Find the point at which branch values are enumerated
			target_line = -1
			for li in range(len(lines)):
				line = lines[li]
				# Looking for this line:  List of extant and reconstructed sequences
				if line.strip() == "List of extant and reconstructed sequences":
					# Two lines down is a pair of integers, n m, indicating n=number of sequences and m=number of nucleotide sites
					flds = lines[li+2].strip().split()
					num_seqs = int(flds[0])
					seq_len = int(flds[1])
					# Actual start of the sequences is four lines down
					target_line = li + 4
					break
			if target_line < 0:
				print "# Could not find ancestral sequences; are you sure you set RateAncestor=1?"
				return
			line = lines[target_line]
			li = target_line
			ancestral_sequences = {}
			# Read the specified number of sequences
			for ai in range(num_seqs):
				flds = line.strip().split('  ')
				seq_id = flds[0].strip()
				seq = flds[-1].replace(" ","")
				ancestral_sequences[seq_id] = seq
				li = li+1
				line = lines[li]
			return ancestral_sequences

	#---------------------------------------------------------------------
	def getMultipleRates(self):
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
	def getPairwiseGYRates(self):
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
	def getPairwisePhysicalRates(self):
		return self.getMultiplePhysicalRates()

	#---------------------------------------------------------------------
	def getMultiplePhysicalRates(self):
		"""Returns the physical distance(s) between two sequences after a run of CodeML.

		Returns a 4-tuple, with the entries:
		(ML nonsynonymous rate per site, ML synonymous rate per site, ML number of nonsynonymous sites, ML number of synonymous sites)
		Raises PAMLError if there is a problem."""

		debugging = False
		if self.seq_type == 'protein':
			raise PAMLError, "Protein distances not supported for multiple sequences"
		elif self.seq_type == 'codon':
			distances = []

			if not os.path.isfile(self.tmpfile):
				raise PAMLError, "Cannot find distance-containing temporary file %s." % self.tmpfile
			else:
				f = file(self.tmpfile, 'r')
				lines = f.readlines()
				f.close()
				# For runmode = -2 (pairwise):
				#    We must have NSsites = 0, and there will be one line of distances (NSsites > 0 is a PAML error)
				# For runmode = 0 (user tree):
				#    For NSsites = 0, there will be two lines of distances which can be added
				#    For NSsites > 0, there will be multiple distances which must be weighted
				#    by the proportion of sites with a given w (=dN/dS).
				if self.options['runmode'] == '-2':
					# One line only
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
							break
				elif self.options['runmode'] == '0':
					if self.options['NSsites']=='0':
						# Sum the branch lengths
						(dn, ds, nsyn, nns) = (0.0, 0.0, 0.0, 0.0)
						for i in range(len(lines)):
							if 'four-fold sites)' in lines[i]:
								# Target is line i+1
								line = lines[i+1]
								# Line should look like this:
								#             dN*=  0.00736 dS*=  0.15299 S* = 200.65 N* = 657.35
								assert "dN*" in line and "dS*" in line
								flds = line.split('=')
								dn += float(flds[1].split()[0])
								ds += float(flds[2].split()[0])
								nsyn = float(flds[3].split()[0])
								nns = float(flds[4].split()[0])
					else: # NSsites > 0
						raise PAMLError, "Rate retrieval for NSsites > 0 not implemented yet"
						# Multiple site types: weight the results
						site_rates = []
						probs = self.getSiteProportions()
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
				else:
					raise PAMLError, "Rate retrieval for runmode %s not implemented yet" % (self.options['runmode'],)

				return (dn, ds, nns, nsyn)
		else:
			raise PAMLError, "PAML sequence type of %s is invalid." % self.seq_type

	#---------------------------------------------------------------------
	def getPhysicalDistRST(self):
		"""Returns the physical-sites distances between two sequences after a run of CodeML.

		Returns (nonsynonymous rate per physical nonsynonymous site, synonymous rate per physical synonymous site)
		"""
		# rst line for fix_kappa=0, fix_alpha=1, fix_omega=0
		#seq seq        N       S       dN       dS     dN/dS   tree t   kappa     w       -lnL
		#  2   1   1789.7    670.3   0.0316   1.2116   0.0261   1.0593   2.1322   0.0261 -4213.851
		debugging = False
		if self.seq_type == 'protein':
			raise PAMLError, "Protein distances not supported for multiple sequences"
		elif self.seq_type == 'codon':
			distances = []
			rstfilename = 'rst1'
			if not os.path.isfile(rstfilename):
				raise PAMLError, "Cannot find distance-containing temporary file %s." % os.path.join(os.getcwd(),rstfilename)
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
	def getSiteProportions(self):
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

	def getCodonSelectionCoefficients(self):
		"""Retrieves codon substitution selection coefficients
		"""
		if not os.path.isfile("rst"):
			raise PAMLError, "Cannot find distance-class-probability-containing temporary file %s." % "rst"
		f = file("rst",'r')
		lines = f.readlines()
		f.close()

		# Find first line of codon table results
		for i in range(len(lines)):
			if (''.join(lines[i].strip().split())).startswith('IJij2'):
				start_line = i+2
				break

		# Headers
		# I       J       ij      2Ns_IJ  pMut_IJ pSub_IJ 2Ns_JI  pMut_JI pSub_JI
		codon_selection_dict = {}
		for line in lines[start_line:]:
			if line.strip() == '':
				break
			#example line:
			#TTC     TTT     CT      -1.91356        0.22056 0.07306 1.91356 0.03255 0.07306 0.07306 0.07306
			flds = line.strip().split()
			res = CodonSelectionResult()
			res.codon_I = flds[0]
			res.codon_J = flds[1]
			res.FI_FJ = float(flds[3])
			res.pMut_IJ = float(flds[4])
			res.pSub_IJ = float(flds[5])
			res.FJ_FI = float(flds[6])
			res.pMut_JI = float(flds[7])
			res.pSub_JI = float(flds[8])
			key = "%s-%s" % (res.codon_I, res.codon_J)
			codon_selection_dict[key] = res
		return codon_selection_dict


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

