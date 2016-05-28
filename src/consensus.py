#! /usr/bin/python
import sys, os, math, string, random
import translate, stats

_genetic_code = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L',
		'CTA':'L', 'CTG':'L', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'GTT':'V',
		'GTC':'V', 'GTA':'V', 'GTG':'V', 'TCT':'S', 'TCC':'S', 'TCA':'S',
		'TCG':'S', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'ACT':'T',
		'ACC':'T', 'ACA':'T', 'ACG':'T', 'GCT':'A', 'GCC':'A', 'GCA':'A',
		'GCG':'A', 'TAT':'Y', 'TAC':'Y', 'TAA':'STOP', 'TAG':'STOP',
		'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'AAT':'N', 'AAC':'N',
		'AAA':'K', 'AAG':'K', 'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
		'TGT':'C', 'TGC':'C', 'TGA':'STOP', 'TGG':'W', 'CGT':'R',
		'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGT':'S', 'AGC':'S', 'AGA':'R',
		'AGG':'R', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}

def get_distance(x,y):
	assert(len(x)==len(y))
	d = 0
	for (a,b) in zip(x,y):
		if a != b:
			d += 1
	return d

def get_consensus_residue(aas):
	# aas represents a column in a multiple alignment.
	freqs = [(aas.count(x),x) for x in 'ACDEFGHIKLMNPQRSTVWXY-']
	# Now figure out the most frequent amino acid.
	return sorted(freqs,reverse=True)[0]

def make_consensus_protein_from_genes(genes):
	prots = []
	for gene in genes:
		p = translate.Translate(gene)
		prots.append(p)
	return make_consensus_protein(prots)

def make_consensus_protein(prots):
	aligns = prots #align_and_remove_gaps(prots)
	(cons_prot, cons_freq) = make_consensus_protein_from_alignment(aligns)
	return cons_prot, cons_freq, aligns

def make_consensus_protein_from_alignment(aligns):
	# Create a consensus protein, and for each position a frequency
	# of the consensus amino acid.
	cons_prot = []
	cons_freq = []
	for i in range(len(aligns[0])):
		(aa, freq) = get_consensus_residue([p[i] for p in aligns])
		cons_prot.append(aa)
		cons_freq.append(freq/float(len(aligns)))
	cons_prot = ''.join(cons_prot)
	return cons_prot, cons_freq

def reverse_translate(prot):
	# Generate one gene for this protein
	rev_genetic_code = dict([(v, k) for (k, v) in _genetic_code.items()])
	gene = []
	for aa in prot:
		gene.append(rev_genetic_code[aa])
	return ''.join(gene)

def make_sequence_from_frequencies(site_cum_probabilities):
	"""Make a sequence from frequencies alone.
	site_cum_probabilities = ..."""
	n = len(site_cum_probabilities)
	new_sequence = ''
	for site in range(n):
		# Ensure that the site really gets changed.
		new_aa = pick_residue(site_cum_probabilities[site])
		new_sequence += new_aa
	return new_sequence
def make_sequence_at_fixed_distance_from_frequencies(distance, starting_sequence, site_cum_probabilities):
	"""Make a sequence differing at a fixed number of sites from a starting sequence.
	distance = number of desired differences.
	starting_sequence = sequence to which changes are made.
	site_cum_probabilities = ..."""
	all_sites = range(len(starting_sequence))
	# Now remove sites for which there is no variation possible given the probabilities,
	# since we must change the sites ultimately chosen to achieve the desired distance
	sites = []
	for site in all_sites:
		(cum_prob, aa) = site_cum_probabilities[site][0]
		if cum_prob < 1.0:
			sites.append(site)
	if distance > len(sites):
		raise Exception, "Cannot create sequence with %d differences, because only %d sites allow variation." % (distance, len(sites))
	sites_to_change = random.sample(sites, distance) # get sites to change
	new_sequence = [x for x in starting_sequence]
	for site in sites_to_change:
		new_aa = pick_residue(site_cum_probabilities[site], new_sequence[site])
		# Ensure that the site really gets changed.
		while new_aa == new_sequence[site]:
			new_aa = pick_residue(site_cum_probabilities[site], new_sequence[site])
		new_sequence[site] = new_aa
	return ''.join(new_sequence)

def make_sequence_at_fixed_distance(distance, starting_sequence, letters='ACDEFGHIKLMNPQRSTVWY'):
	"""Make a sequence differing at a fixed number of sites from a starting sequence.
	distance = number of desired differences.
	starting_sequence = sequence to which changes are made."""
	sites = range(len(starting_sequence))
	if distance > len(sites):
		raise Exception, "Cannot create sequence with %d differences, because only %d sites allow variation." % (distance, len(sites))
	sites_to_change = random.sample(sites, distance) # get sites to change
	new_sequence = [x for x in starting_sequence]
	for site in sites_to_change:
		new_aa = random.choice(letters)
		# Ensure that the site really gets changed.
		while new_aa == new_sequence[site]:
			new_aa = random.choice(letters)
		new_sequence[site] = new_aa
	return ''.join(new_sequence)
def pick_residue(site_cum_probabilities, avoid_residue = ''):
	"""Choose a residue based on its frequency
	site_cum_probabilities stores the cumulative probabilities."""
	r = random.random()
	for (freq, aa) in site_cum_probabilities:
		if freq > r and not (aa==avoid_residue):
			break
	return aa

def build_site_abs_probabilities(alignment, weights=None, letters='ACDEFGHIKLMNPQRSTVWY-'):
	"""Create a list of probability distributions and associated amino acids
	from the specified alignment.
	"""
	nseq = len(alignment)
	if weights is None:
		weights = [1.0/nseq]*nseq
	site_probabilities_list = []
	for site in range(len(alignment[0])): # index over sites in each sequence
		aas = [y[site] for y in alignment]
		probs = dict([(x, 0.0) for x in letters])
		# Accumulate weights for this site
		for iseq in range(nseq):
			aa = aas[iseq]
			wt = weights[iseq]
			probs[aa] += wt
		#n = float(len(aas))
		#for aa in aas:
		#	probs[aa] += 1.0/n
		prob_list = [(p, aa) for (aa, p) in probs.items()]
		#print prob_list
		site_probabilities_list.append(prob_list)
	return site_probabilities_list

def build_site_abs_probabilities_aa(alignment, weights=None, letters='ACDEFGHIKLMNPQRSTVWY-'):
	"""Create a list of probability distributions and associated amino acids
	from the specified alignment.
	Optionally apply weights for each site of each sequence.
	"""
	nseq = len(alignment)
	nsites = len(alignment[0])
	if weights is None:
		weights = [1.0/nseq]*nseq
	site_probabilities_list = []
	for site in range(nsites): # index over sites in each sequence
		aas = [y[site] for y in alignment]
		probs = dict([(x, 0.0) for x in letters])
		# Accumulate weights for this site
		for iseq in range(nseq):
			aa = aas[iseq]
			wt = weights[iseq]
			probs[aa] += wt
		site_probabilities_list.append(probs.items())
	return site_probabilities_list

def build_site_abs_probabilities_aa_biased(alignment, biases, weights=None, letters='ACDEFGHIKLMNPQRSTVWY-'):
	"""Create a list of probability distributions and associated amino acids
	from the specified alignment.
	Optionally apply weights for each site of each sequence.
	Correct for biases in the aa distribution.
	"""
	nseq = len(alignment)
	nsites = len(alignment[0])
	if weights is None:
		weights = [1.0/nseq]*nseq
	site_probabilities_list = []
	for site in range(nsites): # index over sites in each sequence
		aas = [y[site] for y in alignment]
		probs = dict([(x, 0.0) for x in letters])
		# Accumulate weights for this site
		for iseq in range(nseq):
			aa = aas[iseq]
			wt = weights[iseq]
			bias = biases[site][aa]
			probs[aa] += wt/bias
		site_probabilities_list.append(probs.items())
	return site_probabilities_list

def make_weighted_consensus(alignment, weights=None):
	nseq = len(alignment)
	if weights is None:
		weights = [1.0/nseq]*nseq
	consensus_prot = ""
	probs_list = build_site_abs_probabilities_aa(alignment, weights)
	for site in range(len(alignment[0])): # index over sites in each sequence
		probs = probs_list[site]
		sorted_plist = [(p, aa) for (aa,p) in probs]
		sorted_plist.sort()
		consensus_prot += sorted_plist[-1][1]
	return consensus_prot

def cumulative_probabilities(abs_probs):
	prob_list = [(p, aa) for (aa, p) in abs_probs]
	prob_list.sort() #cmp=lambda x, y: y<x)  # sort in order of decreasing probability
	prob_list.reverse()
	#print prob_list
	cum_probs = []
	cum = 0.0
	for (p, aa) in prob_list:
		if p > 0.0:
			cum += p
		cum_probs.append((cum, aa))
		if cum >= 1.0:
			break
	return cum_probs

def build_site_probabilities(alignment):
	"""Create a list of cumulative probability distributions and associated amino acids
	from the specified alignment.
	"""
	site_cum_probabilities_list = []
	for site in range(len(alignment[0])): # index over sites in each sequence
		aas = [y[site] for y in alignment]
		probs = dict([(x, 0.0) for x in 'ACDEFGHIKLMNPQRSTVWY'])
		# Count all the amino acids in the column and divide by number of column entries
		n = float(len(aas))
		for aa in aas:
			if aa != '-':
				probs[aa] += 1.0/n
		cum_probs = cumulative_probabilities(abs_probs.items())
		site_cum_probabilities_list.append(cum_probs)
	return site_cum_probabilities_list
def build_stability_change_table_from_queries(query_sequences, alignment, reference_sequence):
	# For each site, find the stability change in going from all possible amino acids in the query
	# sequences to the reference sequence.
	# Store stability_changes for retrieval based on (site, from, to) in a dictionary
	stability_change_dict = {}
	for site in range(len(reference_sequence)):
		from_aas = dict([(seq[site], seq[site]) for seq in query_sequences]).keys()
		to_aa = reference_sequence[site]
		for from_aa in from_aas:
			if from_aa == to_aa:
				ddG = 0.0
			else:
				ddG = estimate_stability_change(from_aa, to_aa, site, alignment, reference_sequence)
			stability_change_dict[(site, from_aa, to_aa)] = ddG
			#print (site, from_aa, to_aa), ddG
	return stability_change_dict

def build_mean_stability_change_table(alignment, reference_seqs, letters='ACDEFGHIKLMNPQRSTVWY'):
	# For each site, find the stability change in going from all amino acids in letters
	# to each sequence in reference_seqs.
	# Store stability_changes for retrieval based on (site, from, to) in a dictionary
	# Return the mean stability change for each substitution over all reference seqs.
	stability_changes_dict = {}
	for ref_seq in reference_seqs:
		stab_dict = build_stability_change_table(alignment, ref_seq, letters)
		for (k, v) in stab_dict.items():
			if stability_changes_dict.has_key(k):
				stability_changes_dict[k].append(v)
			else:
				stability_changes_dict[k] = [v]
	mean_ddg_dict = {}
	for (k, ddgs) in stability_changes_dict.items():
		mean_ddg_dict[k] = stats.Mean(ddgs)
	return mean_ddg_dict

def build_stability_change_table(alignment, reference_sequence, letters='ACDEFGHIKLMNPQRSTVWY'):
	# For each site, find the stability change in going from all amino acids in letters
	# to the reference sequence.
	# Store stability_changes for retrieval based on (site, from, to) in a dictionary
	stability_change_dict = {}
	seqIDs = []
	for seq in alignment:
		seqIDs.append(get_distance(seq, reference_sequence)/float(len(seq)))
	for site in range(len(reference_sequence)):
		from_aas = [x for x in letters]
		to_aa = reference_sequence[site]
		for from_aa in from_aas:
			if from_aa == to_aa:
				ddG = 0.0
			else:
				ddG = estimate_stability_change(from_aa, to_aa, site, alignment, reference_sequence, seqIDs, 0.00198, 300, 1.0)
			stability_change_dict[(site, from_aa, to_aa)] = ddG
			#print (site, from_aa, to_aa), ddG
	return stability_change_dict

def build_uncorrected_stability_change_table_with_pseudocounts(alignment, letters='ACDEFGHIKLMNPQRSTVWY'):
	# For each site, find the stability change in going from all amino acids in letters
	# to the reference sequence.  Add 1 to all counts.
	# Store stability_changes for retrieval based on (site, from, to) in a dictionary
	stability_change_dict = {}
	for site in range(len(alignment[0])):
		from_aas = [x for x in letters]
		to_aas = [x for x in letters]
		for from_aa in from_aas:
			for to_aa in to_aas:
				if from_aa == to_aa:
					ddG = 0.0
				else:
					ddG = estimate_uncorrected_stability_change(from_aa, to_aa, site, alignment, 0.00198, 300, 1.0)
				stability_change_dict[(site, from_aa, to_aa)] = ddG
				#print (site, from_aa, to_aa), ddG
	return stability_change_dict

def build_stability_change_table_with_pseudocounts(alignment, reference_sequence, letters='ACDEFGHIKLMNPQRSTVWY'):
	# For each site, find the stability change in going from all amino acids in letters
	# to the reference sequence.  Add 1 to all counts.
	# Store stability_changes for retrieval based on (site, from, to) in a dictionary
	stability_change_dict = {}
	seqIDs = []
	for seq in alignment:
		seqIDs.append(get_distance(seq, reference_sequence)/float(len(seq)))
	mean_seqID = stats.Mean(seqIDs)
	for site in range(len(reference_sequence)):
		from_aas = [x for x in letters]
		to_aas = [x for x in letters]
		for from_aa in from_aas:
			for to_aa in to_aas:
				if from_aa == to_aa:
					ddG = 0.0
				else:
					ddG = estimate_stability_change(from_aa, to_aa, site, alignment, reference_sequence, seqIDs, 0.00198, 300, 1.0) #mean_seqID)
				stability_change_dict[(site, from_aa, to_aa)] = ddG
				#print (site, from_aa, to_aa), ddG
	return stability_change_dict

def estimate_uncorrected_stability_change(from_aa, to_aa, site, alignment, R=0.00198, T=300, pseudo_count_addition=0.0):
	"""Estimate the stability change of an amino acid substitution using the pseudo-equilibrium hypothesis.
	No correction for phylogenetic relationships is performed.
	"""
	Z_from_aa = get_uncorrected_partition_sum(from_aa, site, alignment) + pseudo_count_addition
	Z_to_aa = get_uncorrected_partition_sum(to_aa, site, alignment) + pseudo_count_addition
	ddG = 0.0
	if Z_from_aa == 0:  # Never happens if pseudocount > 0.0
		# What to do if nothing about from_aa is known?
		# Presumably some value smaller than Z_to_aa should  be used.
		# A simple solution is just to assume that the residue actually appears once
		# in the MSA (rather than zero times), and contributes the full weight (1.0).
		#print Z_to_aa, Z_from_aa
		Z_from_aa = 1.0
		#print Z_to_aa, Z_from_aa, ddG
	try:
		ddG = R*T*math.log(Z_to_aa/Z_from_aa)
	except OverflowError, oe:
		print Z_from_aa, Z_to_aa
	return ddG

def get_uncorrected_partition_sum(aa, site, alignment):
	Z = 0.0
	for seq in alignment:
		weight = 1.0
		if seq[site] != aa:
			weight = 0.0
		Z += weight
	return Z

def estimate_stability_change(from_aa, to_aa, site, alignment, reference_sequence, seqIDs, R=0.00198, T=300, pseudo_count_addition=0.0):
	"""Estimate the stability change of an amino acid substitution using the pseudo-equilibrium hypothesis.
	We use the method of Godoy-Ruiz et al., J. Mol. Biol. 336:313-318 (2004).
	"""
	Z_from_aa = get_partition_sum(from_aa, site, alignment, reference_sequence, seqIDs) + pseudo_count_addition
	Z_to_aa = get_partition_sum(to_aa, site, alignment, reference_sequence, seqIDs) + pseudo_count_addition
	ddG = 0.0
	if Z_from_aa == 0:  # Never happens if pseudocount > 0.0
		# What to do if nothing about from_aa is known?
		# Presumably some value smaller than Z_to_aa should  be used.
		# A simple solution is just to assume that the residue actually appears once
		# in the MSA (rather than zero times), and contributes the full weight (1.0).
		#print Z_to_aa, Z_from_aa
		Z_from_aa = 1.0
		#print Z_to_aa, Z_from_aa, ddG
	try:
		ddG = R*T*math.log(Z_to_aa/Z_from_aa)
	except OverflowError, oe:
		print Z_from_aa, Z_to_aa
	return ddG

def get_partition_sum(aa, site, alignment, reference_sequence, seqIDs):
	Z = 0.0
	for a in range(len(alignment)):
		seq = alignment[a]
		ref_i = reference_sequence[site]
		seq_i = seq[site]
		seqID = seqIDs[a]
		weight = partition_weight(seq_i, ref_i, aa, seqID)
		Z += weight
	return Z

def partition_weight(s_i, ref_i, aa, seqID):
	if s_i != aa:
		res = 0.0
	elif ref_i != aa:
		res = seqID
	else:
		res = 1.0-seqID
	return res

def get_protein_stability_difference(from_seq, to_seq, stability_change_dict):
	stability_diff = 0.0
	for site in range(len(from_seq)):
		from_aa = from_seq[site]
		to_aa = to_seq[site]
		if to_aa != from_aa:
			ddG = stability_change_dict[(site, from_aa, to_aa)]
			stability_diff += ddG
	return stability_diff
def get_protein_stability_differences(from_seq, to_seq, ref_seq, stability_change_dict):
	stability_diff = 0.0
	for site in range(len(from_seq)):
		from_aa = from_seq[site]
		to_aa = to_seq[site]
		ref_aa = ref_seq[site]
		if to_aa != from_aa:
			ddG_from = 0.0
			if from_aa != ref_aa:
				ddG_from = stability_change_dict[(site, from_aa, ref_aa)]
			ddG_to = 0.0
			if to_aa != ref_aa:
				ddG_to = stability_change_dict[(site, to_aa, ref_aa)]
			stability_diff += (ddG_to - ddG_from)
	return stability_diff
