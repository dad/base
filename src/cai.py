#!/usr/bin/python
"""Module for manipulating aspects of codon bias.

By D. Allan Drummond, 2004-2012, incorporating material by Jesse Bloom, 2004.
"""

import re, os, sys, string, math, random
import translate, stats

class BioUtilsError(Exception):
	"""Error using one of the bio utils."""

def frequencyToProbability(codon_freq, nucleotide_freq):
	# Turn frequencies into probabilities
	total_codons = float(sum(codon_freq.values()))
	codon_prob = dict([(codon, codon_freq[codon]/total_codons) for codon in codon_freq.keys()])
	total_nucleotides = float(sum(nucleotide_freq.values()))
	codon_prob_from_nucleotide = {}
	for codon in codon_freq.keys():
		cprob = math.exp(sum([math.log(nucleotide_freq[nt]/total_nucleotides) for nt in codon]))
		codon_prob_from_nucleotide[codon] = cprob
	return (codon_prob, codon_prob_from_nucleotide)

def getCodonProbabilities(seq, pseudocount):
	codon_freq = dict([(c,pseudocount) for c in translate.AADNACodons()])
	nucleotide_freq = dict([(nt,pseudocount) for nt in 'ATGC'])
	_accumulateCodonFrequencies(seq, codon_freq, nucleotide_freq)
	(codon_prob, codon_prob_from_nucleotide) = frequencyToProbability(codon_freq, nucleotide_freq)
	return (codon_prob, codon_prob_from_nucleotide, codon_freq, nucleotide_freq)

def _accumulateCodonFrequencies(seq, codon_freq, nucleotide_freq):
	codons = split(seq)
	for codon in codons:
		try:
			codon_freq[codon]+=1.0
		except KeyError,ke:
			pass
	for nt in seq:
		try:
			nucleotide_freq[nt] += 1.0
		except KeyError,ke:
			pass

def getCodonProbabilitiesForMultipleSequences(seqs, pseudocount):
	codon_freq = dict([(c,pseudocount) for c in translate.AADNACodons()])
	nucleotide_freq = dict([(nt,pseudocount) for nt in 'ATGC'])
	# The only reason to do this, versus just _accumulateCodonFrequencies(''.join(seqs)...), is a hedge against
	# the possibility that some sequences have bad lengths, and we want to keep everything in frame; this method
	# resets the reading frame every gene.
	for seq in seqs:
		_accumulateCodonFrequencies(seq, codon_freq, nucleotide_freq)
	(codon_prob, codon_prob_from_nucleotide) = frequencyToProbability(codon_freq, nucleotide_freq)
	return (codon_prob, codon_prob_from_nucleotide, codon_freq, nucleotide_freq)

def getSYangNielsen(seq, codon_prob, codon_prob_from_nt, pseudocount):
	gc = translate.geneticCode(rna=False)
	# Compute normalizations for conditional probabilities of a codon given an amino acid
	alt_codons = {}
	cond_normalization = {}
	for aa in translate.AAs():
		codons = translate.getCodonsForAA(aa, rna=False)
		alt_codons[aa] = codons
		cond_normalization[aa] = sum([codon_prob[c] for c in codons])
	sum_sc = 0.0
	gene_codons = split(seq)
	for to_codon in gene_codons:
		aa = gc[to_codon]
		if not aa == '*':
			sum_sc_i = 0.0
			# Go over all alternative codons and compute the average selection coefficient for moving from that codon to this one
			for from_codon in alt_codons[aa]:
				p_fcod_given_aa = codon_prob[from_codon]/cond_normalization[aa]
				s_from_to = math.log(codon_prob[to_codon]/codon_prob[from_codon]) - math.log(codon_prob_from_nt[to_codon]/codon_prob_from_nt[from_codon])
				sum_sc_i += p_fcod_given_aa*s_from_to
			sum_sc += sum_sc_i
	return sum_sc/len(gene_codons)

def getSYangNielsenNoMutBias(seq, codon_prob, codon_prob_from_nt, pseudocount):
	gc = translate.geneticCode(rna=False)
	# Compute normalizations for conditional probabilities of a codon given an amino acid
	alt_codons = {}
	cond_normalization = {}
	for aa in translate.AAs():
		codons = translate.getCodonsForAA(aa, rna=False)
		alt_codons[aa] = codons
		cond_normalization[aa] = sum([codon_prob[c] for c in codons])
	sum_sc = 0.0
	gene_codons = split(seq)
	for to_codon in gene_codons:
		aa = gc[to_codon]
		if not aa == '*':
			sum_sc_i = 0.0
			# Go over all alternative codons and compute the average selection coefficient for moving from that codon to this one
			for from_codon in alt_codons[aa]:
				p_fcod_given_aa = codon_prob[from_codon]/cond_normalization[aa]
				# Eliminate the mutational bias term
				s_from_to = math.log(codon_prob[to_codon]/codon_prob[from_codon]) # - math.log(codon_prob_from_nt[to_codon]/codon_prob_from_nt[from_codon])
				sum_sc_i += p_fcod_given_aa*s_from_to
			sum_sc += sum_sc_i
	return sum_sc/len(gene_codons)

class CodingFrequencies(object):
	def __init__(self, pseudocount=0):
		# Track the frequency of each codon relative to its synonyms
		self.codon_freqs = {}
		self.aa_counts = {}
		self.total_aa_count = 0
		# Track the frequency of each nucleotide relative to all nucleotides
		self.nucleotide_freqs = {}
		self.nucleotide_counts = 0
		self.gc = translate.geneticCode(rna=False)
		self.pseudocount = pseudocount

		for nt in 'ACGT':
			self.nucleotide_freqs[nt] = 0
		for aa in translate.AAsAndStop():
			codons = translate.getCodonsForAA(aa, rna=False)
			for codon in codons:
				self.codon_freqs[codon] = 0
				self.aa_counts[aa] = 0

	def getPseudocount(self):
		return self.pseudocount

	def countCodon(self, codon):
		try:
			self.codon_freqs[codon] += 1
			self.aa_counts[self.gc[codon]] += 1
			self.total_aa_count += 1
			for i in range(3):
				self.countNucleotide(codon[i])
		except KeyError:
			# Couldn't find the codon...continue
			pass

	def countNucleotide(self, nucleotide):
		self.nucleotide_freqs[nucleotide] += 1
		self.nucleotide_counts += 1

	def addCodons(self, codons):
		for codon in codons:
			self.countCodon(codon)
	
	def addGene(self, gene):
		self.addCodons(splitByFrame(gene,0))
	
	def addGenes(self, genes):
		for g in genes:
			self.addGene(g)

	def getNucleotideProportion(self, nucleotide):
		nt_count = self.nucleotide_counts
		res = 0.0
		if nt_count > 0:
			res = (self.getNucleotideCount(nucleotide) + self.pseudocount)/(float(nt_count) + self.pseudocount)
		return res

	def estimateNucleotideProportion(self, nucleotide):
		nt_count = self.nucleotide_counts
		res = 0.0
		if nt_count > 0:
			res = (self.getNucleotideCount(nucleotide) + self.pseudocount)/(float(nt_count) + self.pseudocount)
		return res

	def getCodonProportion(self, codon):
		aa = self.gc[codon]
		aa_count = self.aa_counts[aa]
		res = 0.0
		if aa_count > 0:
			res = (self.codon_freqs[codon] + self.pseudocount)/(float(aa_count) + self.pseudocount)
		return res

	def estimateCodonProportion(self, codon):
		aa = self.gc[codon]
		aa_count = self.aa_counts[aa]
		res = 1.0/self.total_aa_count # If we've never seen this aa before
		if aa_count > 0:
			res = (self.codon_freqs[codon] + self.pseudocount)/(float(aa_count) + self.pseudocount)
		return res

	def getNucleotideCount(self, nt):
		return self.nucleotide_freqs[nt]

	def getCodonCount(self, codon):
		return self.codon_freqs[codon]

	def getAACount(self, aa):
		return self.aa_counts[aa]
	
	# Relative synonymous codon usage as defined in Sharp and Li NAR 1987.
	# RSCU(codon) = n(codon)*degeneracy(aa)/n(aa)
	def getRelativeSynonymousCodonUsage(self):
		# RSCU is a dictionary keyed by codon
		rscu = {}
		gc = translate.geneticCode(rna=False)
		aa_codons = translate.DNACodons()
		for codon in aa_codons:
			ncodon = self.getCodonCount(codon) + self.pseudocount
			naa = self.getAACount(gc[codon]) + self.pseudocount
			deg = len(translate.getSynonyms(codon,rna=False))
			rscu[codon] = ncodon*deg/float(naa)
		return rscu

	def getRelativeAdaptiveness(self):
		rscus = self.getRelativeSynonymousCodonUsage()
		return self.getRelativeAdaptivenessFromRSCUs(rscus, rna=False)

	def getRelativeAdaptivenessFromRSCUs(self, rscu_dict):
		relad = {}
		for aa in translate.AAsAndStop():
			codons = translate.getCodonsForAA(aa, rna=False)
			# Normalize for maximum RSCU in this family
			max_rscu = max([rscu_dict[c] for c in codons])
			for c in codons:
				relad[c] = rscu_dict[c]/max_rscu
		return relad

	def getSelectionCoefficient(self, from_codon, to_codon):
		sc = None
		gc = translate.geneticCode(rna=False)
		if gc[from_codon] == gc[to_codon]:
			# Y = to, X = from
			# s_X->Y = ln p_X/p_Y - sum_i ln p_X_i/p_Y_i
			cp = self.estimateCodonProportion
			np = self.estimateNucleotideProportion
			from_p = cp(from_codon)
			to_p = cp(to_codon)
			assert to_p > 0.0
			sc = math.log(to_p/from_p) - sum([math.log(np(to_codon[i])/np(from_codon[i])) for i in range(3)])
		return sc
	
	def write(self, outstream):
		for aa in translate.AAsAndStop():
			outstream.write('{0}\t{1}\n'.format(aa, self.getAACount(aa)))
		for aa in translate.AAsAndStop():
			for codon in sorted(translate.getCodonsForAA(aa)):
				outstream.write('{0}\t{1}\t{2}\n'.format(aa, codon, self.getCodonCount(codon)))

def estimateSelectionCoefficientForCodon(codon, cons_cf, var_cf):
	aa = translate.translate(codon)
	ln_p_I_cons = math.log(cons_cf.estimateCodonProportion(codon))
	ln_p_I_var = math.log(var_cf.estimateCodonProportion(codon))
	ln_p_i_sum = 0.0
	for i in range(3):
		ln_p_i_cons = math.log(cons_cf.estimateNucleotideProportion(codon[i]))
		ln_p_i_var = math.log(var_cf.estimateNucleotideProportion(codon[i]))
		ln_p_i_sum += ln_p_i_cons - ln_p_i_var
	scaled_s = ln_p_I_cons - ln_p_I_var - ln_p_i_sum
	n_aas = cons_cf.getAACount(aa) + var_cf.getAACount(aa)
	return scaled_s, n_aas

def estimateSelectionCoefficients(cons_cf, var_cf):
	"""Takes two cai.CodonFrequencies objects; returns estimated selection coefficients for each codon."""
	# 2Ns_var-->cons = ln p_I(cons)/p_I(var) - sum_i ln [p_i(cons)/p_i(var)]
	codon_sc_dict = {}
	for aa in translate.degenerateAAs():
		for codon in translate.getCodonsForAA(aa, rna=False):
			(scaled_s, n_aas) = estimateSelectionCoefficientForCodon(codon, cons_cf, var_cf)
			codon_sc_dict[codon] = (scaled_s, n_aas)
	return codon_sc_dict

def test_estimateSelectionCoefficients():
	gene = 'ATGGATTATACCTACGACTACACTTATGATTATACCTAC'
	prot = translate.translate(gene)
	other_genes = ['ATGGATTATACCTACGACTACACTTATGACTACACTTAT']
	other_prots = [translate.translate(g) for g in other_genes]
	(cons_cdna, var_cdna) = splitByConservation(conservedAAconservedCodon, gene, prot, other_genes, other_prots, 0)
	assert cons_cdna == 'GATTATACCTACGACTACACTTAT'
	assert var_cdna == 'GATTATACCTAC'
	# MDYTYDYTYDYTY
	cons_cf = CodingFrequencies(pseudocount=1)
	var_cf = CodingFrequencies(pseudocount=1)
	cons_cf.addCodons(split(cons_cdna))
	var_cf.addCodons(split(var_cdna))
	fwd_sc_dict = estimateSelectionCoefficients(cons_cf, var_cf)
	rev_sc_dict = estimateSelectionCoefficients(var_cf, cons_cf)
	for codon in split(var_cdna):
		(fwd_sc, fwd_n) = fwd_sc_dict[codon]
		(rev_sc, rev_n) = rev_sc_dict[codon]
		assert fwd_sc == -rev_sc
	print "# test_estimateSelectionCoefficients passed"

def randomizeConservationCategoryForAA(aa, cons_cf, var_cf):
	cons_pseudo = cons_cf.getPseudocount()
	var_pseudo = var_cf.getPseudocount()
	codons = []
	for codon in translate.getCodonsForAA(aa, rna=False):
		n_cons = cons_cf.getCodonCount(codon)
		n_var = var_cf.getCodonCount(codon)
		codons += [codon]*(n_cons + n_var - cons_pseudo - var_pseudo)
		total_n_cons += n_cons - cons_pseudo
		total_n_var += n_var - var_pseudo
		random.shuffle(codons)
		# Now portion them out
		rand_cons_cf.addCodons(codons[0:total_n_cons])
		rand_var_cf.addCodons(codons[total_n_cons:])

def randomizeConservationCategoryForCodon(codon, cons_cf, var_cf):
	cons_pseudo = cons_cf.getPseudocount()
	var_pseudo = var_cf.getPseudocount()
	codons = []
	n_cons = cons_cf.getCodonCount(codon)
	n_var = var_cf.getCodonCount(codon)
	codons += [codon]*(n_cons + n_var - cons_pseudo - var_pseudo)
	total_n_cons += n_cons - cons_pseudo
	total_n_var += n_var - var_pseudo
	random.shuffle(codons)
	# Now portion them out
	rand_cons_cf.addCodons(codons[0:total_n_cons])
	rand_var_cf.addCodons(codons[total_n_cons:])

def randomizeConservationCategory(cons_cf, var_cf):
	"""Takes two cai.CodonFrequencies objects; returns two objects in which the codon counts are preserved but the category of each codon is randomized."""
	cons_pseudo = cons_cf.getPseudocount()
	var_pseudo = var_cf.getPseudocount()
	rand_cons_cf = CodingFrequencies(cons_pseudo)
	rand_var_cf = CodingFrequencies(var_pseudo)
	# Make merged lists of codons for each amino acid
	# If var has n codons, select n at random from the merged list and add
	for aa in translate.degenerateAAs():
		# Assemble the merged list of codons
		codons = []
		total_n_cons = 0
		total_n_var = 0
		for codon in translate.getCodonsForAA(aa, rna=False):
			n_cons = cons_cf.getCodonCount(codon)
			n_var = var_cf.getCodonCount(codon)
			codons += [codon]*(n_cons + n_var) # - cons_pseudo - var_pseudo)
			total_n_cons += n_cons #- cons_pseudo
			total_n_var += n_var #- var_pseudo
		random.shuffle(codons)
		# Now portion them out
		rand_cons_cf.addCodons(codons[0:total_n_cons])
		rand_var_cf.addCodons(codons[total_n_cons:])
		assert len(codons) == total_n_cons + total_n_var
		assert rand_cons_cf.getAACount(aa) == cons_cf.getAACount(aa)
		assert rand_var_cf.getAACount(aa) == var_cf.getAACount(aa)
	return rand_cons_cf, rand_var_cf

def test_randomizeConservationCategory():
	gene = 'ATGGATTATACCTACGACTACACTTATGATTATACCTACGATTATACCTACGACTACACTTATGATTATACCTAC'
	prot = translate.translate(gene)
	other_genes = ['ATGGATTATACCTACGACTACACTTATGACTACACTTATGATTATACCTACGACTACACTTATGACTACACTTAC']
	other_prots = [translate.translate(g) for g in other_genes]
	(cons_cdna, var_cdna) = splitByConservation(conservedAAconservedCodon, gene, prot, other_genes, other_prots, 0)
	assert cons_cdna == 'GATTATACCTACGACTACACTTATGATTATACCTACGACTACACTTATTAC'
	assert var_cdna == 'GATTATACCTACGATTATACC'
	# MDYTYDYTYDYTY
	cons_cf = CodingFrequencies(pseudocount=1)
	var_cf = CodingFrequencies(pseudocount=1)
	cons_cf.addCodons(split(cons_cdna))
	var_cf.addCodons(split(var_cdna))
	sc_dict = estimateSelectionCoefficients(cons_cf, var_cf)
	target_codon = 'TAC'
	rscs = []
	(sc, nc) = sc_dict[target_codon]
	for n in range(10):
		(rc, rv) = randomizeConservationCategory(cons_cf, var_cf)
		r_sc = estimateSelectionCoefficients(rc,rv)
		(rsc, nc) = r_sc[target_codon]
		rscs.append(rsc)
	rscs.sort()
	#print rscs[0], sc, rscs[-1]
	assert rscs[0] <= sc
	assert rscs[-1] >= sc
	#(n, m, sd, se) = stats.statsSummary(rscs)
	#print (sc - m)/sd
	print "# test_randomizeConservationCategory passed"

# Goal: compute p_I and p_i for i=A
def getEmpiricalFrequencies(cdna, pseudocount=0):
	cf = CodingFrequencies(pseudocount)
	#cdna = cdna.upper().replace('U','T')
	codons = split(cdna)
	for codon in codons:
		cf.countCodon(codon)
	return cf

def test_getEmpiricalFrequencies():
	gene = 'ATGGATTATACCTAC'
	cf = getEmpiricalFrequencies(gene)
	assert cf.getNucleotideProportion('A') == 5.0/15
	assert cf.getCodonProportion('TAC') == 1.0/2
	assert cf.getCodonProportion('GGG') == 0.0
	cf = getEmpiricalFrequencies(gene, 1)
	assert cf.getNucleotideProportion('A') == 6.0/16
	assert cf.getCodonProportion('TAC') == 2.0/3
	assert cf.getCodonProportion('GGG') == 0.0
	print "# test_getEmpiricalFrequencies passed"

# Conservation functions...
def conservedAAconservedCodon(site, master_prot, other_prots, master_codon, other_codons):
	other_aas = [p[site] for p in other_prots]
	res = None
	# If no gaps and no amino-acid differences...
	if (not ('-' in other_aas)) and (other_aas.count(master_prot[site]) == len(other_aas)):
		# ...return True if all codons identical, False otherwise
		res = other_codons.count(master_codon) == len(other_codons)
	return res

def conservedAA(site, master_prot, other_prots, master_codon, other_codons):
	other_aas = [p[site] for p in other_prots]
	res = None
	# If no gaps...
	if not '-' in other_aas:
		# ...return True if all amino acids identical, False otherwise
		res = other_aas.count(master_prot[site]) == len(other_aas)
	return res

def conservedCodon(site, master_prot, other_prots, master_codon, other_codons):
	other_aas = [p[site] for p in other_prots]
	res = None
	# If no gaps...
	if not ('-' in other_aas):
		# ...return True if all codons identical, False otherwise
		# Note that this cannot be true unless all the amino acids are also identical.
		res = other_codons.count(master_codon) == len(other_codons)
	return res

def notVariableAA(site, master_prot, other_prots, master_codon, other_codons):
	other_aas = [p[site] for p in other_prots]
	res = None
	# If no gaps...
	if not '-' in other_aas:
		# ...return False if all amino acids are different, True otherwise
		res = not (len(set([master_prot[site]] + other_aas)) == len(other_aas)+1)
	return res

def allConservedAllVariableAA(site, master_prot, other_prots, master_codon, other_codons):
	other_aas = [p[site] for p in other_prots]
	res = None
	# If no gaps...
	if not '-' in other_aas:
		# ...return True if all amino acids identical, False if all different, None otherwise
		s = set([master_prot[site]] + other_aas)
		n = len(s)
		if n == 1:
			res = True
		elif n == (len(other_aas)+1):
			res = False
	return res

def conservedAAunconservedCodon(site, master_prot, other_prots, master_codon, other_codons):
	other_aas = [p[site] for p in other_prots]
	res = None
	# If no gaps, and no amino-acid differences, but at least one codon difference...
	if (not ('-' in other_aas)) and (other_aas.count(master_prot[site]) == len(other_aas)):
		# ...return False if all codons identical, True otherwise
		res = other_codons.count(master_codon) < len(other_codons)
	return res

def splitByConservation(conservationFxn, master_cdna, master_prot, other_cdnas, other_prots, n_terminal_start=0):
	""" Divide cDNA into conserved and variable portions based on conservationFxn.
	"""
	gc = translate.geneticCode(rna=False)
	all_tables = []
	# One table per amino acid class per gene
	dna_codon_map = {}
	for codon in translate.AADNACodons():
		dna_codon_map[codon] = True

	# Split into codons and replace U with T
	other_codon_seqs = [split(s.upper().replace('U','T')) for s in other_cdnas]
	master_codon_seq = split(master_cdna.upper().replace('U','T'))

	# Determine the point at which we begin counting: that's n_terminal_start amino acids in.
	# n_terminal_start is also the number of amino acids to skip; if n_terminal_start = 1,
	# then we skip one amino acid.
	naa = 0
	i = 0
	while naa < n_terminal_start:
		if master_prot[i] != '-':
			naa += 1
		i += 1
	gapped_n_terminal_start = i
	# Can't be any closer to the beginning than n_terminal_start
	assert gapped_n_terminal_start >= n_terminal_start

	cons_cdna = ""
	var_cdna = ""

	for i in range(gapped_n_terminal_start, len(master_prot)):
		master_codon = master_codon_seq[i]
		if master_codon in ['---', 'ATG','TGG','TAG','TGA','TAA']:
			# skip gaps, M, W and stops
			continue
		try:
			aa = gc[master_codon]
			if aa in translate.degenerateAAs(): # skip M, W and stops
				other_codons = [s[i] for s in other_codon_seqs]
				conserved = conservationFxn(i, master_prot, other_prots, master_codon, other_codons)
				if conserved is None:
					continue
				if conserved:
					cons_cdna += master_codon
				else:
					var_cdna += master_codon
		except KeyError:
			# If codon not in genetic code
			continue
	return cons_cdna, var_cdna

def test_splitByConservation():
	gene = 'ATGGATTATACCTAC'
	prot = 'MDYTY'
	other_genes = ['ATGGATTATACCTAC','ATGGACTATACCTAT']
	other_prots = [translate.translate(g) for g in other_genes]
	(cons_cdna, var_cdna) = splitByConservation(conservedAAconservedCodon, gene, prot, other_genes, other_prots, 0)
	assert cons_cdna == 'TATACC'
	assert var_cdna == 'GATTAC'

	gene = 'TGGGATTATACCTAC'
	prot = translate.translate(gene) # 'WDYTY'
	other_genes = ['TGGGATTATAGCTAC','TGGGATTATACCTAC']
	other_prots = [translate.translate(g) for g in other_genes]
	(cons_cdna, var_cdna) = splitByConservation(conservedAA, gene, prot, other_genes, other_prots, 0)
	assert cons_cdna == 'GATTATTAC'
	assert var_cdna == 'ACC'
	print "# test_splitByConservation passed"


# The functions getCodonCounts() and getAkashi2x2TablesForORF() are for inverse-Akashi analyses.
def getCodonCounts(conservationFxn, master_cdna, master_prot, other_cdnas, other_prots, n_terminal_start=0):
	gc = translate.geneticCode(rna=False)
	all_tables = []
	# One table per amino acid class per gene
	conserved_codon_counts = {}
	variable_codon_counts = {}
	dna_codon_map = {}
	for codon in translate.AADNACodons():
		conserved_codon_counts[codon] = 0
		variable_codon_counts[codon] = 0
		dna_codon_map[codon] = True

	# Split into codons and replace U with T
	other_codon_seqs = [split(s.replace('U','T')) for s in other_cdnas]
	master_codon_seq = split(master_cdna.replace('U','T'))

	# Determine the point at which we begin counting: that's n_terminal_start amino acids in.
	# n_terminal_start is also the number of amino acids to skip; if n_terminal_start = 1,
	# then we skip one amino acid.
	naa = 0
	i = 0
	while naa < n_terminal_start:
		if master_prot[i] != '-':
			naa += 1
		i += 1
	gapped_n_terminal_start = i
	# Can't be any closer to the beginning than n_terminal_start
	assert gapped_n_terminal_start >= n_terminal_start

	for i in range(gapped_n_terminal_start, len(master_prot)):
		master_codon = master_codon_seq[i]
		if master_codon in ['---', 'ATG','TGG','TAG','TGA','TAA']:
			# skip gaps, M, W and stops
			continue
		try:
			other_codons = [s[i] for s in other_codon_seqs]
			found = dna_codon_map[master_codon]  # Make sure we know about this codon
			conserved = conservationFxn(i, master_prot, other_prots, master_codon, other_codons)
			if conserved is None:
				continue
			if conserved:
				conserved_codon_counts[master_codon] += 1
			else:
				variable_codon_counts[master_codon] += 1
		except KeyError:
			continue
	return conserved_codon_counts, variable_codon_counts

def getAkashi2x2TablesForORF(conservationFxn, aligned_cdna, aligned_prot, other_aligned_cdnas, other_aligned_prots, pseudocount, n_terminal_start=0):
	# Now build the tables for this gene
	# Compute conserved--preferred association using each codon as preferred in turn
	gc = translate.geneticCode(rna=False)
	(conserved_codon_counts, variable_codon_counts) = getCodonCounts(conservationFxn, aligned_cdna, aligned_prot, other_aligned_cdnas, other_aligned_prots, n_terminal_start)
	gene_codon_tables = {}
	for codon in translate.AADNACodons():
		gene_codon_tables[codon] = []
	#for aa in translate.degenerateAAs():
	for aa in translate.AAs():
		codons = translate.getCodonsForAA(aa, rna=False)
		for codon in codons:
			# Get contingency table for each codon
			# cons-pref, cons-un, var-pref, var-un
			cons_not_codon = sum([conserved_codon_counts[x] for x in conserved_codon_counts.keys() if (gc[x] == aa) and (x != codon)])
			var_not_codon = sum([variable_codon_counts[x] for x in variable_codon_counts.keys() if (gc[x] == aa) and (x != codon)])
			# Add the pseudocount to each entry
			table = [conserved_codon_counts[codon]+pseudocount, cons_not_codon+pseudocount, variable_codon_counts[codon]+pseudocount, var_not_codon+pseudocount]
			gene_codon_tables[codon].append(tuple(table))
	return gene_codon_tables

def getAkashi2x2TablesForORFRefCodon(conservationFxn, reference_codon_dict, aligned_cdna, aligned_prot, other_aligned_cdnas, other_aligned_prots, pseudocount, n_terminal_start=0):
	# Now build the tables for this gene
	# Compute conserved--preferred association using each codon as preferred in turn
	gc = translate.geneticCode(rna=False)
	(conserved_codon_counts, variable_codon_counts) = getCodonCounts(conservationFxn, aligned_cdna, aligned_prot, other_aligned_cdnas, other_aligned_prots, n_terminal_start)
	gene_codon_tables = {}
	for codon in translate.AADNACodons():
		gene_codon_tables[codon] = []
	for aa in translate.degenerateAAs():
		codons = translate.getCodonsForAA(aa, rna=False)
		for codon in codons:
			# Get contingency table for each codon
			# cons-pref, cons-un, var-pref, var-un
			# Get reference codon
			ref_codon = reference_codon_dict[codon]
			# Add the pseudocount to each entry
			table = [conserved_codon_counts[codon]+pseudocount,
					 conserved_codon_counts[ref_codon]+pseudocount,
					 variable_codon_counts[codon]+pseudocount,
					 variable_codon_counts[ref_codon]+pseudocount]
			gene_codon_tables[codon].append(tuple(table))
	return gene_codon_tables

def getAkashi2x2TablesForORFAllCodonPairs(conservationFxn, codon_pair_list, aligned_cdna, aligned_prot, other_aligned_cdnas, other_aligned_prots, pseudocount, n_terminal_start=0):
	# Now build the tables for this gene
	# Compute conserved--preferred association using each codon as preferred in turn
	gc = translate.geneticCode(rna=False)
	(conserved_codon_counts, variable_codon_counts) = getCodonCounts(conservationFxn, aligned_cdna, aligned_prot, other_aligned_cdnas, other_aligned_prots, n_terminal_start)
	# Collect counts by codon pair
	gene_codon_tables = {}
	for k in codon_pair_list:
		gene_codon_tables[k] = []
	for aa in translate.degenerateAAs():
		for (codon, ref_codon) in codon_pair_list:
			# Get contingency table for each codon/ref-codon pair
			# cons-pref, cons-un, var-pref, var-un
			# Add the pseudocount to each entry
			table = [conserved_codon_counts[codon]+pseudocount,
					 conserved_codon_counts[ref_codon]+pseudocount,
					 variable_codon_counts[codon]+pseudocount,
					 variable_codon_counts[ref_codon]+pseudocount]
			gene_codon_tables[(codon,ref_codon)].append(tuple(table))
	return gene_codon_tables



def getFractionRare(gene, rel_adapt, cutoff):
	"""Returns the proportion of genes with codons below a cutoff of relative adaptiveness"""
	n_rare = 0
	n_tot = 0
	for i in range(0, len(gene), 3):
		codon = gene[i:i+3]
		# Only check codons which can be optimized, which
		# excludes Met, Trp, and the 3 stop codons
		if codon not in ['ATG','TGG','TGA','TAA','TAG']:
			n_tot += 1
			relad = rel_adapt[codon]
			if relad <= cutoff:
				n_rare += 1
	return float(n_rare)/n_tot

def getCAI(gene, ln_rel_adapt):
	"""Returns the CAI of a gene, or 'None' if that gene is invalid.

	Reference: Sharp and Li, Nucleic Acids Res, 15:1281-1295, 1987."""
	if gene == []:
		raise BioUtilsError, "Empty gene."
	try:
		log_sum = n = 0
		codons = splitByFrame(gene,0)
		for codon in codons:
			if codon in ['ATG','TGG','AUG','UGG']:
				continue # only one of these codons
			try:
				log_sum += ln_rel_adapt[codon]
				n += 1
			except KeyError:
				pass
				#if i + 3 == len(gene) and codon in ['TAA', 'TAG', 'TGA','UAA', 'UAG', 'UGA']:
				#	pass # this is last codon, so don't worry if stop codon
				#else:
				#	pass
				#raise BioUtilsError, "Codon %s is invalid." % codon
		CAI = math.exp(log_sum / float(n))
		return CAI
	except ZeroDivisionError:
		raise BioUtilsError, "Empty gene."
	except IndexError:
		raise BioUtilsError, "Gene is of invalid length."


def getFop(gene, optimal_list, pseudocount=0.0):
	(nOpt, nTot) = getNumOptimalCodons(gene.replace("U","T"), optimal_list)
	fop = 0.0
	if nTot + pseudocount > 0:
		fop = (float(nOpt) + pseudocount)/nTot
	return fop

def isOptimal(codon, optimal_list):
	return codon in optimal_list

def getNumOptimalCodons(gene, optimal_list):
	nOpt = 0
	nTot = 0
	for i in range(0, len(gene), 3):
		codon = gene[i:i+3]
		# Only check codons which can be optimized, which
		# excludes Met, Trp, and the 3 stop codons
		if codon not in ['ATG','TGG','TGA','TAA','TAG']:
			nTot += 1
			if codon in optimal_list:
				nOpt += 1
	return nOpt, nTot

def split(nts):
	return splitByFrame(nts,0)

def split_by_frame(nts, frame_index):
	return splitByFrame(nts, frame_index)

def splitByFrame(nts, frame_index):
	n = int(math.floor(len(nts)/3.0))
	codon_list = [nts[3*i+frame_index:3*i+3+frame_index] for i in range(n)]
	if 3*n+frame_index < len(nts):
		codon_list += [nts[3*n+frame_index:len(nts)]]
	return codon_list

def test_splitByFrame():
	# Tests the split_by_frame function
	gene1 = 'ATGGATTAGAC'
	gene2 = 'GATACCCAG'
	print split_by_frame(gene1, 0)
	print split_by_frame(gene1, 1)
	print split_by_frame(gene1, 2)
	print split_by_frame(gene2, 0)
	print split_by_frame(gene2, 1)
	print split_by_frame(gene2, 2)

def getGC(gene, pseudocount = 0.0):
	gene = gene.upper()
	return (gene.count('G') + gene.count('C') + pseudocount)/float(len(gene))

def getGC3(gene, pseudocount = 0.0):
	codons = split_by_frame(gene.upper(), 0)
	return (len([x for x in codons if len(x)==3 and (x[2]=='G' or x[2]=='C')]) + pseudocount)/float(len(codons))

def getGCi(gene, indices, pseudocount = 0.0):
	codons = split_by_frame(gene.upper(), 0)
	sum_gc = 0.0
	for i in indices:
		assert i>0 and i<=3 # codon boundaries; one-based!
		sum_gc += len([x for x in codons if len(x)==3 and (x[i-1] in 'GC')])
	return (sum_gc + pseudocount)/(float(len(indices)*len(codons)))

def getContent(gene, nts, pseudocount=0.0):
	cont = 0
	for x in gene:
		if x in nts:
			cont += 1
	return (cont+pseudocount)/float(len(gene))


def test_getGC3():
	gene1 = 'ATGGATTAGAC'
	gene2 = 'GATACCCAG'
	print getGC3(gene2) == 2.0/3.0
	print getGC(gene2) == 5.0/9

def test_getGCi():
	gene1 = 'ATGGATTAGAC'
	gene2 = 'GATACCCAG'
	print getGCi(gene2,[3]) == getGC3(gene2)
	print getGCi(gene2,[1,2,3]) == getGC(gene2)
	print getGCi(gene1,[1,2]) == 1.0/8
	print getGCi(gene2,[1,2]) == 3.0/6

def getFopGCEndingOld(gene, optimal_list):
	codons = split_by_frame(gene.upper(), 0)
	gc_ending = [x for x in codons if x[2] in "GC" and not x in ['ATG','TGG','TGA','TAA','TAG']]
	gc_ending_and_optimal = [x for x in codons if x[2] in "GC" and not x in ['ATG','TGG','TGA','TAA','TAG'] and x in optimal_list]
	return float(len(gc_ending_and_optimal))/len(gc_ending)

def getFopGCEnding(gene, optimal_list):
	prot = translate.TranslateRaw(gene.upper())
	four_plus_six_fold = 'LVSPTARG'
	codons = split_by_frame(gene.upper(), 0)
	assert len(codons) >= len(prot)-1  # possibility of stop codon.
	gc_codons = [codons[i] for i in range(len(codons)) if prot[i] in four_plus_six_fold and codons[i][2] in 'GC']
	return getFop(''.join(gc_codons), optimal_list)

def getDinucleotideIndex(gene, dn):
	'''Compute the frequency of a dinucleotide relative to the product of the frequencies of its constituent nucleotides.	'''
	dn = dn.upper()
	nt1 = dn[0]
	nt2 = dn[1]
	L = len(gene)
	nt1_freq = gene.count(nt1)/float(L)
	nt2_freq = gene.count(nt2)/float(L)
	if nt1 == nt2:
		# Must special-case this because .count() does not work properly
		frame1 = [x[i:i+2] for i in range(0,L-3,2)]
		frame2 = [x[i:i+2] for i in range(1,L-2,2)]
		# Q: what if length is odd vs. even?
		dn_freq = (frame1.count(dn)+frame2.count(dn))/float(L-1)
	else:
		dn_freq = gene.count(dn)/float(L-1)
	index = 0.0
	if nt1_freq>0 and nt2_freq>0:
		index = dn_freq/(nt1_freq*nt2_freq)
	return index

# Per Lavner and Kotlar, Gene 245 (2005)
def getFopNucleotideCorrected(gene, optimal_list, noncoding_sequence_fragments):
	# For each codon, compute expected frequency in noncoding region
	# Use simple definition from Lavner and Kotlar of E_i^(nc)(g), in which
	# the noncoding frequency of each nucleotide in the third position is
	# used to weight each codon.
	nts = 'ATGCU'
	ntfreqs = {}
	noncoding_seq = noncoding_sequence_fragments
	if isinstance(noncoding_sequence_fragments, list):
		noncoding_seq = ''.join(noncoding_sequence_fragments)
	len_noncoding = float(len(noncoding_seq))
	for nt in nts:
		ntf = noncoding_seq.count(nt)
		# Since zeros will lead to divide-by-zero errors, substitute
		# one nucleotide count if, by chance, the noncoding
		# region doesn't have any of this nucleotide.  Bad!
		if ntf == 0:
			ntf = 1
		ntfreqs[nt] = ntf/len_noncoding

	fop = 0.0
	ntot = 0
	for codon in split_by_frame(gene,0):
		if codon in optimal_list:
			fop += 1.0/ntfreqs[codon[2]]
			ntot += 1
	return fop/ntot


def getSum(table):
	(a,b,c,d) = table
	return a+b+c+d

def getMHOptimalityTables(left, right, nopFxn):
	# table:
	#              Left |  Right
	#   Optimal     A   |    B
	#  --------------------------------
	# Non-optimal   C   |    D
	ignore = ['ATG', 'TGG', 'TAA', 'TAG' 'TGA']
	tables = {}
	for c in 'ACDEFGHIKLNPQRSTVY':
		tables[c] = [0,0,0,0]

	# Left side of table
	if left != '':
		prot = translate.Translate(left)
		#print noncons, prot
		for i in range(len(prot)):
			codon = left[3*i:3*i+3]
			if not codon in ignore:
				(opt, tot) = nopFxn(codon)
				[a,b,c,d] = tables[prot[i]]
				if opt > 0: # optimal codon on left
					a += 1
				else: # non-optimal codon
					c += 1
				tables[prot[i]] = [a,b,c,d]

	# Right side of table
	if right != '':
		prot = translate.Translate(right)
		#print cons, prot
		for i in range(len(prot)):
			codon = right[3*i:3*i+3]
			if not codon in ignore:
				(opt, tot) = nopFxn(codon)
				[a,b,c,d] = tables[prot[i]]
				if opt > 0:
					# optimal codon on right
					b += 1
				else:
					d += 1
				tables[prot[i]] = [a,b,c,d]

	return [x for x in tables.values() if (getSum(x) > 0)]

genetic_code = translate.geneticCode(rna=False)
'''
genetic_code = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L',
		'CTA':'L', 'CTG':'L', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'GTT':'V',
		'GTC':'V', 'GTA':'V', 'GTG':'V', 'TCT':'S', 'TCC':'S', 'TCA':'S',
		'TCG':'S', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'ACT':'T',
		'ACC':'T', 'ACA':'T', 'ACG':'T', 'GCT':'A', 'GCC':'A', 'GCA':'A',
		'GCG':'A', 'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
		'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'AAT':'N', 'AAC':'N',
		'AAA':'K', 'AAG':'K', 'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
		'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W', 'CGT':'R',
		'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGT':'S', 'AGC':'S', 'AGA':'R',
		'AGG':'R', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}
'''
# DAD: should use translate content here
_aa_codons = [k for k in genetic_code.keys() if genetic_code[k] != "*"]
_syn_codons = dict([(v,k) for (k,v) in genetic_code.items() if v != "*"]).values() + ['TGA','TAG','TAA']
_degeneracy = dict([(c, genetic_code.values().count(genetic_code[c])) for c in genetic_code.keys()])
_aas = 'ACDEFGHIKLMNPQRSTVWY'
_deg_aas = 'ACDEFGHIKLNPQRSTVY'


# Relative synonymous codon usage as defined in Sharp and Li NAR 1987.
# RSCU(codon) = n(codon)*degeneracy(aa)/n(aa)
def getRelativeSynonymousCodonUsage(gene, pseudocount=0, ):
	prot = translate.TranslateRaw(gene)
	codons = splitByFrame(gene,0)
	#base_usage = dict([(aa, prot.count(aa)) for aa in _aas])
	rscu = {}
	gc = translate.geneticCode(rna=False)
	for codon in _aa_codons:
		ncodon = codons.count(codon)
		aa = gc[codon]
		naa = prot.count(aa)
		deg = _degeneracy[codon]
		rscu[codon] = ncodon*deg/float(naa)
	return rscu

def getRelativeAdaptiveness(gene):
	rscus = getRelativeSynonymousCodonUsage(gene)
	return getRelativeAdaptivenessFromRSCUs(rscus, rna=False)

def getRelativeAdaptivenessFromRSCUs(rscus, rna=False):
	relad = {}
	for aa in _aas:
		codons = translate.get_codons_for_aa(aa, rna)
		max_rscu = max([rscus[c] for c in codons])
		for c in codons:
			relad[c] = rscus[c]/max_rscu
	return relad

def logRelativeAdaptiveness(relad):
	ln_relad = {}
	for codon in relad.keys():
		r = relad[codon]
		if r > 0:
			ln_relad[codon] = math.log(r)
		else:
			ln_relad[codon] = math.log(0.001) # A small number, but not too small...?
	return ln_relad

# Compute the RSCU for each family, weighted by species-specific
# weights obtained by regressing a target variable (say,
# expression level) on RSCU.
def wRSCU(gene, weights):
	freqs = getRelativeSynonymousCodonUsage(gene)
	wts = [weights[codon]*freqs[codon] for codon in freqs.keys()]
	wrscu_val = sum(wts)/len(wts)
	return wrscu_val

# From inverse-Akashi analysis, unpublished, Drummond and Marx
_m_extorquens_optimal_codons = ['GCC','TGC','GAC','GAG','TTC','GGC','CAT','ATC','AAG','CTG','AAC','CCG','CAG','CGC','TCG','ACC','GTG','TAC']

# From John Peden thesis, citing Sharp et al 1990
_b_subtilis_optimal_codons = ['TTC', 'CTT', 'ATC', 'GTT','GTA', 'TAC',
'CAA', 'AAC', 'AAA', 'GAC', 'GAA', 'TCT', 'CCT', 'CCA', 'ACT', 'GCT',
'CGT', 'CGC', 'GGT']

# From Ikemura 1985
_e_coli_optimal_codons_ikemura = ['AAA', 'AAC', 'ACC', 'ACT', 'ATC', 'CAG', 'CCG', 'CGC', 'CGT', 'CTG', 'GAA', 'GCA', 'GCG',        'GCT', 'GGC', 'GGT', 'GTT', 'TAC',               'TTC']
# From Sharp and Li NAR 1987
_e_coli_optimal_codons_sharp = ['AAA', 'AAC', 'ACC', 'ATC', 'CAC', 'CAG', 'CCG',  'CGT', 'CTG', 'GAA', 'GAC', 'GCT', 'GGT', 'GTT', 'TAC', 'TCT', 'TGC', 'TTC']
_e_coli_optimal_codons_drummond = ['AAA', 'AAC', 'ACC', 'ATT', 'CAG', 'CAT', 'CCG', 'CGT', 'CTG', 'GAA', 'GAT', 'GCC', 'GGT', 'GTG', 'TAC', 'TCG', 'TGT', 'TTT']
_e_coli_optimal_codons = _e_coli_optimal_codons_sharp

# From Stenico et al 1994
_c_elegans_optimal_codons_stenico = ['TTC','CTT','CTC','ATC','GTC',
'TAC','CAC','CAG','AAC','AAG','GAC','GAG','TCC','CCA','ACC','GCT',
'GCC','TGC','CGT','CGC','GGA']
# From Sharp, C.elegans II, http://www.ncbi.nlm.nih.gov/books/bv.fcgi?rid=ce2
_c_elegans_optimal_codons_sharp = ['TTC', 'CTT', 'CTC', 'ATC', 'GTC',
'TAC','CAC','CAG','AAC','GAC','GAG','TCT','TCC','CCA','ACC','GCT',
'GCC','TGC','CGT','CGC','GGA']
_c_elegans_optimal_codons = _c_elegans_optimal_codons_sharp

# From Sharp, P.M. and Cowe "Synonymous codon usage in Saccharomyces
# cerevisiae," Yeast 7:657-678 (1991)
_s_cerevisiae_optimal_codons_sharp = \
['AAC','AAG','ACC','ACT','AGA','ATC','ATT','CAA','CAC','CCA','GAA','GAC','GCT','GGT','GTC','GTT','TAC','TCC','TCT','TGT','TTC','TTG']
# Only lowest-expression genes from S. cer., inverse-Akashi analysis
_s_cerevisiae_optimal_codons_drummond_ia = \
['AAC','AAG','ACA','AGA','ATT','CAG','CAT','CCC','GAG','GAT','GCT','GGT','GTC','TAC','TCT','TGT','TTC','TTG']
# Best codons by SYN analysis, Drummond unpublished 2011.
_s_cerevisiae_optimal_codons_drummond_syn = \
['AAG','AAT','ACT','AGA','ATT','CAA','CAT','CCA','GAA','GAT','GCT','GGT','GTT','TAC','TCT','TGT','TTG','TTT']
# From Kliman et al. "Selection conflicts, gene expression, and codon usage trends in yeast," J Mol Evol 57:98-109 (2003) Table 1, "pref" designation
_s_cerevisiae_optimal_codons_kliman = \
['AAC','AAG','ACC','ACT','AGA','ATC','ATT','CAA','CAC','CCA','CGT','GAA','GAC','GCC','GCT','GGT','GTC','GTT','TAC','TCC','TCT','TGT','TTC','TTG']
_s_cerevisiae_optimal_codons = _s_cerevisiae_optimal_codons_kliman

# Derived by taking the tRNA gene count for each anticodon.
#_a_thaliana_optimal_codons = ['GCT', 'TGC', 'GAC', 'GAG','TTC','GGC',
#'CAC','ATT','AAG','CTT','AAC','CCA','CAG','AGA','CGT','TCT','ACT','GTT',
#'TAC']

# From Chiapello et al 1998
_a_thaliana_optimal_codons = ['AAG','AAC','ATC','ACC','CGT','AGC','TCC','TAC','TTC', 'TGC','CTC', 'CAG', 'CAC','CCC','GAG','GAC','GTC','GCC','GGT','GGC']
# From inverse-Akashi analysis
_a_thaliana_optimal_codons_drummond = ['AAC','AAG','ACA','ATC','CAC','CAG','CCT','CGA','CTT','GAG','GAT','GCT','GGA','GTG','TAC','TCT','TGC','TTC']
# From inverse-Akashi analysis
_a_lyrata_optimal_codons_drummond = ['AAC','AAG','ACA','ATC','CAC','CAG','CCT','CGA','CTC','GAG','GAT','GCT','GGA','GTG','TAC','TCT','TGC','TTC']

# From Akashi H, Genetics 1994
#_d_melanogaster_optimal_codons = ['TTC','CTC','CTG','ATC','GTC','GTG','TCC','TCG',
#'CCC','ACC','GCC','TAC','CAC','CAG','AAC','AAG','GAC','GAG','TGC','CGT','CGC','GGC']
# From Duret and Mouchiroud, PNAS 1999
_d_melanogaster_optimal_codons_duret99 = \
['CGC','CGT','CTC','CTG','TCC','TCG','ACC','CCC','GCC','GGC','GTC','GTG','AAG','AAC','CAG','CAC','GAG','GAC','TAC','TGC','TTC','ATC']
## From inverse-Akashi analysis
_d_melanogaster_optimal_codons_drummond = \
['TGC','GAC','GAG','TTC','CAC','AAG','AAC','CAG','TAC','ATC','GCC','GGC','CCC','ACC','GTG','TTG','CTG','AGG','CGC','AGC','TCC','TCG']

#['CGC','CGT','CTC','CTG','TCC','TCG','ACG','CCC','GCC','GGA','GTC','GTG','AAG','AAC','CAG','CAC','GAG','GAC','TAC','TGC','TTC','ATC']
_d_melanogaster_optimal_codons = _d_melanogaster_optimal_codons_duret99 #drummond

# From Drummond and Wilke Cell 2008
#_sim_optimal_codons = _s_cerevisiae_optimal_codons #['AAC','AAG','ACC','ACT','AGA','ATC','ATT','CAA','CAC','CCA','GAA','GAC','GCT','GGT','GTC','GTT','TAC','TCC','TCT','TGT','TTC','TTG']
_sim_optimal_codons = ['AAC','AAG','ACC','ACT','AGA','ATC','ATT','CAA','CAC','CCA','GAA','GAC','GCT','GGT','GTC','GTT','TAC','TCC','TCT','TGT','TTC','TTG']

# From Comeron Genetics 2004
_h_sapiens_optimal_codons_comeron = ['AAC','AAG','ACC','ATC','CAC','CAG','CCC','CGC','CTG','GAC','GAG','GCC','GGC','GTG','TAC','TGC','TTC']
# From Lavner and Kotlar Gene 2005
_h_sapiens_optimal_codons_lavner = ['AAC','AAG','ACA','ATT','CAC','CAG','CCT','CGT','CTT','GAA','GAC','GCT','GGC','GTT','TAC','TCT','TGC','TTC']
# From Lavner and Kotlar BMC Genomics 2006
_h_sapiens_optimal_codons_lavner06 = ['AAC','ACA','ACT','ATT','CAC','CAG','CCA','CCT','CTG','CTT','GAC','GCT','GGA','GGC','GTG','GTT','TAC','TCT','TGC','TTC']
_optimal_codons_archetti = ['AGA','TTA','TCA','ACA','CCA','GCA','GGA','GTA','ATA','AAA','','','','','','','','']

# From tRNA gene counts, Waterston et al. Nature 2002, Initial sequencing and analysis of the mouse genome.
_h_sapiens_optimal_codons_waterston = ['AAC','AAG','ACT','ATT','CAC','CAG','CCT','CGT','CTT','GAA','GAC','GCT','GGC','GTT','TAC','TCT','TGC','TTC']
_h_sapiens_optimal_codons_waterston_modified = ['AAC','AAG','ACC','ATC','CAC','CAG','CCC','CGC','CTC','GAA','GAC','GCC','GGC','GTC','GTG','TAC','TCC','TGC','TTC']
# From Drummond et al. unpub. 2006, based on Mantel-Haenszel test associating codons with conserved sites in mouse, rat, and dog
_h_sapiens_optimal_codons_drummond = \
 [      'AAA', 'AGA',  'AAT','ACA','ATT','CAC','CAG',   'CCA','CTG',   'GAT','GAA','GCC','GGA','GTG','TAC','TCC','TGC','TTT']
#['AAC','AAG','ACT',            'ATT','CAC','CAG','CCT','CGT','CTT','GAA','GAC','GCT','GGC','GTT','TAC','TCT','TGC','TTC'] # waterston
'''
http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=Homo+sapiens+[gbpri]
UUU 17.5(676381)  UCU 15.1(585967)  UAU 12.1(470083)  UGU 10.5(407020)
UUC 20.4(789374)  UCC 17.7(684663)  UAC 15.3(592163)  UGC 12.6(487907)
UUA  7.6(294684)  UCA 12.2(471469)  UAA  1.0( 38222)  UGA  1.5( 59528)
UUG 12.9(498920)  UCG  4.4(171428)  UAG  0.8( 30104)  UGG 13.2(510256)

CUU 13.1(508151)  CCU 17.5(676401)  CAU 10.8(419726)  CGU  4.6(176691)
CUC 19.6(759527)  CCC 19.8(767793)  CAC 15.1(583620)  CGC 10.5(405748)
CUA  7.2(276799)  CCA 16.9(653281)  CAA 12.2(473648)  CGA  6.2(239573)
CUG 39.8(1539118) CCG  6.9(268884)  CAG 34.2(1323614) CGG 11.5(443753)

AUU 15.9(615699)  ACU 13.1(506277)  AAU 16.9(653566)  AGU 12.1(469641)
AUC 20.9(808306)  ACC 18.9(732313)  AAC 19.1(739007)  AGC 19.5(753597)
AUA  7.4(288118)  ACA 15.0(580580)  AAA 24.3(940312)  AGA 12.1(466435)
AUG 22.1(853648)  ACG  6.1(234532)  AAG 31.9(1236148) AGG 11.9(461676)

GUU 11.0(426252)  GCU 18.5(715079)  GAU 21.8(842504)  GGU 10.8(416131)
GUC 14.5(562086)  GCC 27.9(1079491) GAC 25.2(973377)  GGC 22.3(862557)
GUA  7.1(273515)  GCA 15.9(614754)  GAA 28.8(1116000) GGA 16.5(637120)
GUG 28.2(1091853) GCG  7.4(286975)  GAG 39.6(1532589) GGG 16.4(636457)
'''
#_h_sapiens_optimal_codons_usage = ['TTC','CTG','ATC','GTG','TCC','CCC','ACC','GCC','TAC','','','','','','','','','','','','']
#_h_sapiens_optimal_codons = _h_sapiens_optimal_codons_waterston_modified
#_h_sapiens_optimal_codons = _h_sapiens_optimal_codons_drummond
_h_sapiens_optimal_codons = _h_sapiens_optimal_codons_comeron

_m_mulatta_optimal_codons = _h_sapiens_optimal_codons

# From tRNA gene counts at http://lowelab.ucsc.edu/GtRNAdb/Mmusc/
_m_musculus_optimal_codons_lowe = ['AAC','AAG',      'ACT','AGA','AGC','ATT','CAC','CAG','CCA','CCT','CGT','CTG',            'GAC','GAG',      'GGC','GTG',      'TAC','TCT','TGC','TTC']
# From tRNA gene counts, Waterston et al. Nature 2002, Initial sequencing and analysis of the mouse genome.
_m_musculus_optimal_codons_waterston = \
 ['AAC','AAG','ACT','ATT','CAC','CAG','CCT','CGT','CTG','GAC','GAG','GCT','GGC','GTG','TAC','TCT','TGC','TTC'] # break ties
#['AAA','AAC','AAG','ACT','AGC','ATT','CAC','CAG','CCT','CGT','CTG','GAA','GAC','GAG','GCT','GGC','GTG','TAC','TCT','TGC','TTC'] # with ties
#['AAG','AAT','ACA','ATT','CAC','CAG','CCC','CGC','CTG','GAC','GAG','GCC','GGA','GTG','TAT','TCC','TGT','TTT']
# From Drummond et al. unpub. 2006, based on Mantel-Haenszel test associating codons with conserved sites in dog, rat and human
_m_musculus_optimal_codons_drummond = \
 ['AAC','AAG','ACC','ATC','CAC','CAG','CCT','CGC','CTG','GAC','GAG','GCC','GGA','GTG','TAT','TCC','TGC','TTT']
#['AAC','AAG','ACC','ATC','CAC','CAG','CCC','CGC','CTG','GAC','GAG','GCC','GGC','GTG','TAC','TGC','TTC'] # comeron
#['AAC','AAG','ACT','ATT','CAC','CAG','CCT','CGT','CTG','GAC','GAG','GCT','GGC','GTG','TAC','TCT','TGC','TTC'] # waterston
#['AAC','AAG','ACC','ATT','CAG','CAT','CCC','CGA','CTG','GAA','GAC','GCC','GGA','GTG','TAT','TCT','TGC','TTC'] # Z=26.69

# Optimal codons from Waterston, but with inosine-modification assumed (T-ending cognates -> C-ending).
_m_musculus_optimal_codons_waterston_modified = \
['AAA','AAC','AAG','ACC','AGC','ATC','CAC','CAG','CCC','CGC','CTG','GAA','GAC','GAG','GCC','GGC','GTG','TAC','TCC','TGC','TTC'] # with ties
#['AAC','AAG','ACC','ATC','CAC','CAG','CCC','CGC','CTG','GAC','GAG','GCC','GGC','GTG','TAC','TCC','TGC','TTC'] # break ties

_m_musculus_optimal_codons = _m_musculus_optimal_codons_waterston_modified
#_m_musculus_optimal_codons = _h_sapiens_optimal_codons
#_m_musculus_optimal_codons = _m_musculus_optimal_codons_drummond
#_h_sapiens_optimal_codons = _m_musculus_optimal_codons_waterston

_t_denticola_optimal_codons_drummond = \
['AAA','AAT','ACA','AGA','ATT','CAG','CAT','CCT','CTT','GAA','GAT','GCT','GGA','GTT','TAT','TCA','TGC','TTT'] # ribosomal
['AAC','AAG','ACC','ATA','CAC','CAG','CCC','CGA','GAC','GAG','GCA','GGC','GTT','TAC','TCG','TGC','TTC','TTG'] # accuracy
_t_denticola_optimal_codons = _t_denticola_optimal_codons_drummond

_t_pallidum_optimal_codons_drummond = \
['AAA','AAC','ACT','ATT','CAA','CAT','CCC','CGT','GAG','GAT','GCA','GGT','GTC','TAC','TCT','TGT','TTG','TTT'] # accuracy
['AAA','AAC','ACA','ATC','CAA','CAC','CCA','CGT','GAG','GAT','GCT','GGC','GTT','TAT','TCA','TGC','TTG','TTT'] # regression on dN
['AAG','AAT','ACT','AGC','AGG','ATC','CAA','CAT','CCG','GAG','GAT','GCT','GGT','GTT','TAT','TGC','TTA','TTC'] # regression on dN/dS
['AAG','AAT','ACG','ATT','CAC','CAG','CCG','CGT','CTT','GAG','GAT','GCG','GGT','GTG','TAC','TCT','TGT','TTT'] # slow dN/dS
['AAC','AAG','ACG','ATT','CAC','CAG','CCG','CGC','CTT','GAG','GAT','GCG','GGT','GTG','TAC','TCT','TGT','TTT'] # slow dN
#['AAG','AAT','ACT','ATT','CAA','CAT','CCT','CGT','GAG','GAT','GCA','GGT','GTG','TAT','TCT','TGT','TTA','TTT']
#['ACT','ATT','CAA','CAT','CCT','CGT','GAG','GAT','GCA','GTG','TCT','TGT','TTA']
#['ACT','ATT','CAA','CAT','CCG','CGC','GAA','GAT','GCA','GTA','TCT','TGT','TTG']
_t_pallidum_optimal_codons = _t_pallidum_optimal_codons_drummond

_h_pylori_optimal_codons_drummond = \
['AAC', 'AAG', 'ACC', 'AGG', 'ATT', 'CAC', 'CAG', 'CCC', 'GAA', 'GAT', 'GCT', 'GGG', 'GTG', 'TAT', 'TCA', 'TGC', 'TTA', 'TTC']
#['AAA','AAC','ACC','AGT','ATT','CAC','CAG','CCG','CGC','GAA','GAT','GCC','GGG','GTG','TAT','TGC','TTA','TTC'] # accuracy
#['AAA','AAT','ACC','AGA','AGC','ATT','CAA','CAT','CCT','GAA','GAT','GCT','GGC','GTG','TAT','TGC','TTA','TTT'] # ribosomal
#['AAA','AAT','ACC','AGC','AGG','ATT','CAA','CAT','CCT','GAA','GAT','GCT','GGC','GTG','TAT','TGC','TTA','TTT'] # flagellar
#['GCC','GCA','CGC','CGA','AGA','AGG','AAC','GAC','TGC','CAA','GAA','GAA','GGC','GGA','CAC','ATC','TTA','TTG','CTC','CTA','AAA','TTC','CCC','CCA','TCA','TCC','AGC','ACC','ACA','TAC','GTA','GTC']
_h_pylori_optimal_codons_ribosomal = \
['AAA','AAT','ACC','AGA','AGC','ATT','CAA','CAT','CCT','GAA','GAT','GCT','GGC','GTG','TAT','TGC','TTA','TTT']
_h_pylori_optimal_codons = _h_pylori_optimal_codons_drummond

_b_burgdorferi_optimal_codons_drummond = ['AAA','AAT','ACA','AGA','ATT','CAA','CAT','CCA','GAA','GAT','GCT','GGA','GTG','TAT','TCT','TGC','TTA','TTT']
#['GGA','TGC','TTC','TTC','CTC','ACC','ATC','GCA','AGC','TCA','CGC','TCC','GAA','AAA','AAG','TTG','CAC','CGA','CTA','GGC','AAC','GAC','GTA','ACA','TAC','CCA','AGA','TTA','CAA']
_b_burgdorferi_optimal_codons_ribosomal = ['AAA','AAT','ACT','AGA','ATT','CAA','CAT','CCT','GAA','GAT','GCT','GGT','GTT','TAT','TCT','TGT','TTA','TTT']
_b_burgdorferi_optimal_codons = _b_burgdorferi_optimal_codons_drummond
_b_burgdorferi_relative_adaptiveness = {'AAA':1.0000, 'AAC':0.2553, 'AAG':0.2481, 'AAT':1.0000, 'ACA':1.0000, 'ACC':0.2651, 'ACG':0.1181, 'ACT':0.7610, 'AGA':1.0000, 'AGC':0.3782, 'AGG':0.3136, 'AGT':0.5642, 'ATA':0.7169, 'ATC':0.1485, 'ATT':1.0000, 'CAA':1.0000, 'CAC':0.3750, 'CAG':0.2150, 'CAT':1.0000, 'CCA':0.8804, 'CCC':0.3636, 'CCG':0.1268, 'CCT':1.0000, 'CGA':0.0866, 'CGC':0.0328, 'CGG':0.0131, 'CGT':0.0787, 'CTA':0.2186, 'CTC':0.0600, 'CTG':0.0619, 'CTT':0.7682, 'GAA':1.0000, 'GAC':0.2265, 'GAG':0.3046, 'GAT':1.0000, 'GCA':0.9523, 'GCC':0.2477, 'GCG':0.1338, 'GCT':1.0000, 'GGA':1.0000, 'GGC':0.3255, 'GGG':0.3159, 'GGT':0.5112, 'GTA':0.4667, 'GTC':0.0753, 'GTG':0.1961, 'GTT':1.0000, 'TAC':0.2745, 'TAT':1.0000, 'TCA':0.7246, 'TCC':0.1295, 'TCG':0.0997, 'TCT':1.0000, 'TGC':0.4677, 'TGT':1.0000, 'TTA':1.0000, 'TTC':0.1092, 'TTG':0.4371, 'TTT':1.0000}
_ln_b_burgdorferi_relative_adaptiveness = logRelativeAdaptiveness(_b_burgdorferi_relative_adaptiveness)

_c_crescentus_optimal_codons = ['GCC','TGC','GAC','GAG','TTC','GGC','CAC','ATC','AAG','CTG','AAC','CCG','CAG','CGC','TCG','ACC','GTG','TAC']
#['AAG','AAT','ACC','ATC','CAC','CAG','CCG','CGG','CTG','GAC','GAG','GCC','GGC','GTG','TAC','TCG','TGC','TTC']


#-----------------------------------------------------------------------------------
# From Sharp and Li, NAR 15(3) 1987, Table 1
# dictionary of E. coli relative adaptiveness, used by 'E_coli_CAI'
_e_coli_relative_adaptiveness_sharp = {
     'TTT' : 0.296, 'TTC' : 1.000, 'TTA' : 0.020, 'TTG' : 0.020,
     'CTT' : 0.042, 'CTC' : 0.037, 'CTA' : 0.007, 'CTG' : 1.000,
     'ATT' : 0.185, 'ATC' : 1.000, 'ATA' : 0.003, 'ATG' : 1.000,
     'GTT' : 1.000, 'GTC' : 0.066, 'GTA' : 0.495, 'GTG' : 0.221,
     'TAT' : 0.239, 'TAC' : 1.000,
     'CAT' : 0.291, 'CAC' : 1.000, 'CAA' : 0.124, 'CAG' : 1.000,
     'AAT' : 0.051, 'AAC' : 1.000, 'AAA' : 1.000, 'AAG' : 0.253,
     'GAT' : 0.434, 'GAC' : 1.000, 'GAA' : 1.000, 'GAG' : 0.259,
     'TCT' : 1.000, 'TCC' : 0.744, 'TCA' : 0.077, 'TCG' : 0.017,
     'CCT' : 0.070, 'CCC' : 0.012, 'CCA' : 0.135, 'CCG' : 1.000,
     'ACT' : 0.965, 'ACC' : 1.000, 'ACA' : 0.076, 'ACG' : 0.099,
     'GCT' : 1.000, 'GCC' : 0.122, 'GCA' : 0.586, 'GCG' : 0.424,
     'TGT' : 0.500, 'TGC' : 1.000,                'TGG' : 1.000,
     'CGT' : 1.000, 'CGC' : 0.356, 'CGA' : 0.004, 'CGG' : 0.004,
     'AGT' : 0.085, 'AGC' : 0.410, 'AGA' : 0.004, 'AGG' : 0.002,
     'GGT' : 1.000, 'GGC' : 0.724, 'GGA' : 0.010, 'GGG' : 0.019}

# From Drummond, unpublished; RSCUs from expression-weighted counts
# Expression from Covert et al. NBT 2004, mean of aerobic wild-type expression

_e_coli_relative_adaptiveness_drummond = {'AAA':1.0000, 'AAC':1.0000, 'AAG':0.3035, 'AAT':0.6361, 'ACA':0.2362, 'ACC':1.0000, 'ACG':0.5237, 'ACT':0.4152, 'AGA':0.0571, 'AGC':1.0000, 'AGG':0.0337, 'AGT':0.4679, 'ATA':0.1004, 'ATC':0.9769, 'ATT':1.0000, 'CAA':0.4664, 'CAC':0.9023, 'CAG':1.0000, 'CAT':1.0000, 'CCA':0.3233, 'CCC':0.1723, 'CCG':1.0000, 'CCT':0.2590, 'CGA':0.1104, 'CGC':0.8954, 'CGG':0.1648, 'CGT':1.0000, 'CTA':0.0541, 'CTC':0.1785, 'CTG':1.0000, 'CTT':0.1733, 'GAA':1.0000, 'GAC':0.6871, 'GAG':0.4196, 'GAT':1.0000, 'GCA':0.6087, 'GCC':0.6960, 'GCG':1.0000, 'GCT':0.5199, 'GGA':0.2055, 'GGC':1.0000, 'GGG':0.3054, 'GGT':0.9104, 'GTA':0.4580, 'GTC':0.5604, 'GTG':1.0000, 'GTT':0.8070, 'TAC':0.9040, 'TAT':1.0000, 'TCA':0.3758, 'TCC':0.6116, 'TCG':0.5102, 'TCT':0.6358, 'TGC':1.0000, 'TGT':0.7557, 'TTA':0.1997, 'TTC':0.9316, 'TTG':0.2155, 'TTT':1.0000}


_e_coli_relative_adaptiveness = _e_coli_relative_adaptiveness_drummond


# compute natural logs of relative adaptivenesses
_ln_e_coli_relative_adaptiveness = logRelativeAdaptiveness(_e_coli_relative_adaptiveness)

#---------------------------------------------------------------------------------
# From Sharp and Cowe, Yeast 1991
# Table 4, for high-bias genes, relative synonymous codon usage (RSCU)
# optimal codons:
# ['TTC','TTG','ATT','ATC','GTT','GTC','TAC','CAC','CAA','AAC','AAG','GAC','GAA','TCT','TCC','CCA','ACT','ACC','GCT','TGT','AGA','GGT']

_yeast_rscu = {
     'TTT' : 0.22, 'TTC' : 1.78, 'TTA' : 0.64, 'TTG' : 5.11,
     'CTT' : 0.03, 'CTC' : 0.00, 'CTA' : 0.19, 'CTG' : 0.03,
     'ATT' : 1.51, 'ATC' : 1.49, 'ATA' : 0.00, 'ATG' : 1.00,
     'GTT' : 2.19, 'GTC' : 1.78, 'GTA' : 0.00, 'GTG' : 0.03,
     'TCT' : 3.41, 'TCC' : 2.27, 'TCA' : 0.08, 'TCG' : 0.00,
     'CCT' : 0.25, 'CCC' : 0.02, 'CCA' : 3.73, 'CCG' : 0.00,
     'ACT' : 1.91, 'ACC' : 2.05, 'ACA' : 0.03, 'ACG' : 0.00,
     'GCT' : 3.07, 'GCC' : 0.91, 'GCA' : 0.02, 'GCG' : 0.00,
     'TAT' : 0.12, 'TAC' : 1.88,
     'CAT' : 0.39, 'CAC' : 1.61, 'CAA' : 1.97, 'CAG' : 0.03,
     'AAT' : 0.12, 'AAC' : 1.88, 'AAA' : 0.20, 'AAG' : 1.80,
     'GAT' : 0.77, 'GAC' : 1.23, 'GAA' : 1.96, 'GAG' : 0.04,
     'TGT' : 1.86, 'TGC' : 0.14,               'TGG' : 1.00,
     'CGT' : 0.53, 'CGC' : 0.00, 'CGA' : 0.00, 'CGG' : 0.00,
     'AGT' : 0.10, 'AGC' : 0.14, 'AGA' : 5.45, 'AGG' : 0.01,
     'GGT' : 3.92, 'GGC' : 0.06, 'GGA' : 0.02, 'GGG' : 0.00}

#_yeast_relative_adaptiveness = getRelativeAdaptivenessFromRSCUs(_yeast_rscu)

# From Sharp and Cowe, Yeast 1991
# Table 4, for high-bias genes, relative adaptiveness
_yeast_relative_adaptiveness_sharp_cowe = {
	'AAA':0.111, 'AAC':1.000, 'AAG':1.000, 'AAT':0.064,
	'ACA':0.015, 'ACC':1.000, 'ACG':0.000, 'ACT':0.932,
	'AGA':1.000, 'AGC':0.041, 'AGG':0.002, 'AGT':0.029,
	'ATA':0.000, 'ATC':0.987, 'ATG':1.000, 'ATT':1.000,
	'CAA':1.000, 'CAC':1.000, 'CAG':0.015, 'CAT':0.242,
	'CCA':1.000, 'CCC':0.005, 'CCG':0.000, 'CCT':0.067,
	'CGA':0.000, 'CGC':0.000, 'CGG':0.000, 'CGT':0.097,
	'CTA':0.037, 'CTC':0.000, 'CTG':0.006, 'CTT':0.006,
	'GAA':1.000, 'GAC':1.000, 'GAG':0.020, 'GAT':0.626,
	'GCA':0.007, 'GCC':0.296, 'GCG':0.000, 'GCT':1.000,
	'GGA':0.005, 'GGC':0.015, 'GGG':0.000, 'GGT':1.000,
	'GTA':0.000, 'GTC':0.813, 'GTG':0.014, 'GTT':1.000,
	'TAC':1.000, 'TAT':0.064,
	'TCA':0.023, 'TCC':0.666, 'TCG':0.000, 'TCT':1.000,
				 'TGC':0.075, 'TGG':1.000, 'TGT':1.000,
	'TTA':0.125, 'TTC':1.000, 'TTG':1.000, 'TTT':0.124}

## From Drummond unpublished 2011
## Top 5% of genes by expression, integrating over 14 mRNA datasets
_s_cerevisiae_relative_adaptiveness = {
 	'AAA':0.3814, 'AAC':1.0000, 'AAG':1.0000, 'AAT':0.3759,
 	'ACA':0.2352, 'ACC':0.7998, 'ACG':0.0424, 'ACT':1.0000,
 	'AGA':1.0000, 'AGC':0.1135, 'AGG':0.0415, 'AGT':0.1607,
 	'ATA':0.1409, 'ATC':0.9370, 'ATT':1.0000, 'ATG':1.0000,
 	'CAA':1.0000, 'CAC':1.0000, 'CAG':0.0961, 'CAT':0.6014,
 	'CCA':1.0000, 'CCC':0.0463, 'CCG':0.0410, 'CCT':0.2576,
 	'CGA':0.0173, 'CGC':0.0184, 'CGG':0.0071, 'CGT':0.2366,
 	'CTA':0.1732, 'CTC':0.0288, 'CTG':0.0479, 'CTT':0.0893,
 	'GAA':1.0000, 'GAC':1.0000, 'GAG':0.1188, 'GAT':0.9319,
 	'GCA':0.1202, 'GCC':0.4108, 'GCG':0.0229, 'GCT':1.0000,
 	'GGA':0.0517, 'GGC':0.0905, 'GGG':0.0193, 'GGT':1.0000,
 	'GTA':0.0976, 'GTC':0.7005, 'GTG':0.0824, 'GTT':1.0000,
 	'TAA':1.0000, 'TAC':1.0000, 'TAG':0.1229, 'TAT':0.3537,
 	'TCA':0.2262, 'TCC':0.6525, 'TCG':0.0616, 'TCT':1.0000,
 	'TGA':0.1102, 'TGC':0.1905, 'TGG':1.0000, 'TGT':1.0000,
 	'TTA':0.3499, 'TTC':1.0000, 'TTG':1.0000, 'TTT':0.4590}

# compute natural logs of relative adaptivenesses
_ln_s_cerevisiae_relative_adaptiveness = logRelativeAdaptiveness(_s_cerevisiae_relative_adaptiveness)
# From Stenico, Lloyd and Sharp, NAR 1994, Table 3, highly biased genes
_c_elegans_rscu = {
     'TTT' : 0.07, 'TTC' : 1.93, 'TTA' : 0.03, 'TTG' : 0.85,
     'CTT' : 2.14, 'CTC' : 2.88, 'CTA' : 0.00, 'CTG' : 0.10,
     'ATT' : 0.65, 'ATC' : 2.33, 'ATA' : 0.01, 'ATG' : 1.00,
     'GTT' : 1.21, 'GTC' : 2.41, 'GTA' : 0.12, 'GTG' : 0.27,
     'TCT' : 1.61, 'TCC' : 3.08, 'TCA' : 0.41, 'TCG' : 0.42,
     'CCT' : 0.13, 'CCC' : 0.09, 'CCA' : 3.73, 'CCG' : 0.05,
     'ACT' : 1.01, 'ACC' : 2.84, 'ACA' : 0.10, 'ACG' : 0.05,
     'GCT' : 1.53, 'GCC' : 2.24, 'GCA' : 0.20, 'GCG' : 0.03,
     'TAT' : 0.23, 'TAC' : 1.77,
     'CAT' : 0.52, 'CAC' : 1.48, 'CAA' : 1.29, 'CAG' : 0.71,
     'AAT' : 0.28, 'AAC' : 1.72, 'AAA' : 0.10, 'AAG' : 1.90,
     'GAT' : 0.82, 'GAC' : 1.18, 'GAA' : 0.61, 'GAG' : 1.39,
     'TGT' : 0.36, 'TGC' : 1.64,               'TGG' : 1.00,
     'CGT' : 2.77, 'CGC' : 2.02, 'CGA' : 0.04, 'CGG' : 0.02,
     'AGT' : 0.13, 'AGC' : 0.35, 'AGA' : 1.08, 'AGG' : 0.07,
     'GGT' : 0.45, 'GGC' : 0.13, 'GGA' : 3.40, 'GGG' : 0.02}


_c_elegans_relative_adaptiveness = getRelativeAdaptivenessFromRSCUs(_c_elegans_rscu)
_ln_c_elegans_relative_adaptiveness = logRelativeAdaptiveness(_c_elegans_relative_adaptiveness)

# Drummond unpub 2010, top 5% of genes by expr (geometric mean expression, 727 high-expression genes by Hill et al. Science 2001.)
_c_elegans_relative_adaptiveness = {'AAA':0.3996, 'AAC':1.0000, 'AAG':1.0000, 'AAT':0.6049, 'ACA':0.5189, 'ACC':1.0000, 'ACG':0.1680, 'ACT':0.9139, 'AGA':0.6655, 'AGC':0.4014, 'AGG':0.0524, 'AGT':0.2839, 'ATA':0.0350, 'ATC':1.0000, 'ATT':0.7588, 'CAA':1.0000, 'CAC':1.0000, 'CAG':0.4923, 'CAT':0.7774, 'CCA':1.0000, 'CCC':0.0393, 'CCG':0.0951, 'CCT':0.0869, 'CGA':0.2080, 'CGC':0.5413, 'CGG':0.0629, 'CGT':1.0000, 'CTA':0.0786, 'CTC':0.9564, 'CTG':0.2577, 'CTT':1.0000, 'GAA':0.8727, 'GAC':0.7583, 'GAG':1.0000, 'GAT':1.0000, 'GCA':0.3478, 'GCC':0.7917, 'GCG':0.0989, 'GCT':1.0000, 'GGA':1.0000, 'GGC':0.0824, 'GGG':0.0324, 'GGT':0.1718, 'GTA':0.1736, 'GTC':0.9209, 'GTG':0.3524, 'GTT':1.0000, 'TAC':1.0000, 'TAT':0.5316, 'TCA':0.6633, 'TCC':0.9184, 'TCG':0.5068, 'TCT':1.0000, 'TGC':1.0000, 'TGT':0.5119, 'TTA':0.0771, 'TTC':1.0000, 'TTG':0.5896, 'TTT':0.2780, 'ATG':1.0000, 'TGG':1.0000}
#_ln_c_elegans_relative_adaptiveness = dict([(codon, math.log(ra)) for (codon,ra) in _c_elegans_relative_adaptiveness.items()])
# Drummond unpub 2010, top 5% of genes (411) by spectral-count abundance (from C. von Mering)
_c_elegans_relative_adaptiveness = {'AAA':0.5888, 'AAC':1.0000, 'AAG':1.0000, 'AAT':0.8607, 'ACA':0.6015, 'ACC':0.8226, 'ACG':0.2429, 'ACT':1.0000, 'AGA':0.6592, 'AGC':0.4302, 'AGG':0.0519, 'AGT':0.3817, 'ATA':0.0609, 'ATC':0.9421, 'ATT':1.0000, 'CAA':1.0000, 'CAC':0.9739, 'CAG':0.4735, 'CAT':1.0000, 'CCA':1.0000, 'CCC':0.0784, 'CCG':0.1865, 'CCT':0.1525, 'CGA':0.2837, 'CGC':0.4545, 'CGG':0.0991, 'CGT':1.0000, 'CTA':0.1131, 'CTC':0.8038, 'CTG':0.3011, 'CTT':1.0000, 'GAA':1.0000, 'GAC':0.6615, 'GAG':0.9167, 'GAT':1.0000, 'GCA':0.4238, 'GCC':0.7373, 'GCG':0.1677, 'GCT':1.0000, 'GGA':1.0000, 'GGC':0.1258, 'GGG':0.0559, 'GGT':0.2383, 'GTA':0.2152, 'GTC':0.7704, 'GTG':0.4285, 'GTT':1.0000, 'TAC':1.0000, 'TAT':0.6773, 'TCA':0.7257, 'TCC':0.6963, 'TCG':0.5472, 'TCT':1.0000, 'TGC':1.0000, 'TGT':0.7515, 'TTA':0.0891, 'TTC':1.0000, 'TTG':0.6570, 'TTT':0.3978, 'ATG':1.0000, 'TGG':1.0000}
_ln_c_elegans_relative_adaptiveness = dict([(codon, math.log(ra)) for (codon,ra) in _c_elegans_relative_adaptiveness.items()])


## From Drummond 2009 unpublished.
_d_melanogaster_relative_adaptiveness = {
	'AAA':0.2928, 'AAC':1.0000, 'AAG':1.0000, 'AAT':0.6137,
	'ACA':0.3839, 'ACC':1.0000, 'ACG':0.5351, 'ACT':0.3788,
	'AGA':0.1685, 'AGC':0.9209, 'AGG':0.2363, 'AGT':0.4493,
	'ATA':0.2676, 'ATC':1.0000, 'ATT':0.6287, 'CAA':0.3617,
	'CAC':1.0000, 'CAG':1.0000, 'CAT':0.5830, 'CCA':0.6041,
	'CCC':1.0000, 'CCG':0.7232, 'CCT':0.3252, 'CGA':0.3128,
	'CGC':1.0000, 'CGG':0.3028, 'CGT':0.5296, 'CTA':0.1664,
	'CTC':0.3227, 'CTG':1.0000, 'CTT':0.1987, 'GAA':0.3980,
	'GAC':0.9859, 'GAG':1.0000, 'GAT':1.0000, 'GCA':0.2841,
	'GCC':1.0000, 'GCG':0.3246, 'GCT':0.4138, 'GGA':0.5904,
	'GGC':1.0000, 'GGG':0.1233, 'GGT':0.5163, 'GTA':0.1895,
	'GTC':0.5314, 'GTG':1.0000, 'GTT':0.3633, 'TAC':1.0000,
	'TAT':0.4524, 'TCA':0.3279, 'TCC':1.0000, 'TCG':0.8279,
	'TCT':0.3471, 'TGC':1.0000, 'TGT':0.3337, 'TTA':0.0832,
	'TTC':1.0000, 'TTG':0.3648, 'TTT':0.4520}

#---------------------------------------------------------------------------------
# From Shields, D.C. and Sharp, P.M., NAR 15(19):8023-8040 (1987)
_b_subtilis_relative_adaptiveness_shields = {
	'TTT':0.571, 'TTC':1.000, 'TTA':1.000, 'TTG':0.036,
	'CTT':0.857, 'CTC':0.143, 'CTA':0.500, 'CTG':0.071,
	'ATT':0.500, 'ATC':1.000, 'ATA':0.071, 'ATG':1.000,
	'GTT':1.000, 'GTC':0.188, 'GTA':0.750, 'GTG':0.438,
	'TCT':1.000, 'TCC':0.021, 'TCA':0.458, 'TCG':0.021,
	'CCT':1.000, 'CCC':0.071, 'CCA':0.714, 'CCG':0.143,
	'ACT':1.000, 'ACC':0.033, 'ACA':0.867, 'ACG':0.200,
	'GCT':1.000, 'GCC':0.025, 'GCA':0.275, 'GCG':0.125,
	'TAT':0.500, 'TAC':1.000,
	'CAT':1.000, 'CAC':0.083, 'CAA':1.000, 'CAG':0.214,
	'AAT':0.417, 'AAC':1.000, 'AAA':1.000, 'AAG':0.097,
	'GAT':0.417, 'GAC':1.000, 'GAA':1.000, 'GAG':0.412,
	'TGT':1.000, 'TGC':1.000, 'TGG':1.000,
	'CGT':1.000, 'CGC':0.609, 'CGA':0.022, 'CGG':0.043,
	'AGT':0.125, 'AGC':0.208, 'AGA':0.435, 'AGG':0.022,
	'GGT':0.955, 'GGC':0.773, 'GGA':1.000, 'GGG':0.045}

_b_subtilis_relative_adaptiveness_drummond = {'AAA':1.0000, 'AAC':0.9709, 'AAG':0.3981, 'AAT':1.0000, 'ACA':1.0000, 'ACC':0.3056, 'ACG':0.6558, 'ACT':0.4438, 'AGA':0.8998, 'AGC':1.0000, 'AGG':0.2598, 'AGT':0.3931, 'ATA':0.1986, 'ATC':0.8279, 'ATT':1.0000, 'CAA':1.0000, 'CAC':0.6087, 'CAG':0.8602, 'CAT':1.0000, 'CCA':0.4244, 'CCC':0.1754, 'CCG':1.0000, 'CCT':0.6587, 'CGA':0.3148, 'CGC':0.9839, 'CGG':0.5805, 'CGT':1.0000, 'CTA':0.1799, 'CTC':0.4233, 'CTG':0.8590, 'CTT':1.0000, 'GAA':1.0000, 'GAC':0.6481, 'GAG':0.4666, 'GAT':1.0000, 'GCA':0.9896, 'GCC':0.6693, 'GCG':1.0000, 'GCT':0.9466, 'GGA':0.8784, 'GGC':1.0000, 'GGG':0.4240, 'GGT':0.5969, 'GTA':0.6679, 'GTC':0.8022, 'GTG':0.8042, 'GTT':1.0000, 'TAC':0.6044, 'TAT':1.0000, 'TCA':0.9220, 'TCC':0.5426, 'TCG':0.3580, 'TCT':0.8474, 'TGC':1.0000, 'TGT':0.6646, 'TTA':0.5984, 'TTC':0.5567, 'TTG':0.5787, 'TTT':1.0000}

_b_subtilis_relative_adaptiveness = _b_subtilis_relative_adaptiveness_drummond
_b_subtilis_optimal_codons = [x for (x,v) in _b_subtilis_relative_adaptiveness_drummond.items() if v>=1.0]
# compute natural logs of relative adaptivenesses
_ln_b_subtilis_relative_adaptiveness = {}
for codon in _b_subtilis_relative_adaptiveness:
    _ln_b_subtilis_relative_adaptiveness[codon] = \
			 math.log(_b_subtilis_relative_adaptiveness[codon])


#-----------------------------------------------------------------------------------
# From Lafay et al., Microbiology 146 (2000), Table 1, highly expressed RSCU
# dictionary of H. pylori RSCU, used by 'H_pylori_CAI'
_h_pylori_rscu = {
     'TTT' : 1.49, 'TTC' : 0.51, 'TTA' : 1.99, 'TTG' : 1.80,
     'CTT' : 1.00, 'CTC' : 0.54, 'CTA' : 0.43, 'CTG' : 0.24,
     'ATT' : 1.61, 'ATC' : 1.22, 'ATA' : 0.17, 'ATG' : 1.00,
     'GTT' : 0.93, 'GTC' : 0.47, 'GTA' : 0.55, 'GTG' : 2.05,
     'TAT' : 1.34, 'TAC' : 0.66,
     'CAT' : 1.31, 'CAC' : 0.69, 'CAA' : 1.63, 'CAG' : 0.38,
     'AAT' : 1.18, 'AAC' : 0.82, 'AAA' : 1.43, 'AAG' : 0.57,
     'GAT' : 1.38, 'GAC' : 0.62, 'GAA' : 1.50, 'GAG' : 0.50,
     'TCT' : 1.56, 'TCC' : 0.48, 'TCA' : 0.63, 'TCG' : 0.31,
     'CCT' : 2.01, 'CCC' : 0.56, 'CCA' : 0.91, 'CCG' : 0.52,
     'ACT' : 1.40, 'ACC' : 1.32, 'ACA' : 0.56, 'ACG' : 0.71,
     'GCT' : 1.63, 'GCC' : 0.65, 'GCA' : 0.57, 'GCG' : 1.15,
     'TGT' : 0.65, 'TGC' : 1.35,                'TGG' : 1.00,
     'CGT' : 0.80, 'CGC' : 1.48, 'CGA' : 0.24, 'CGG' : 0.08,
     'AGT' : 0.69, 'AGC' : 2.32, 'AGA' : 2.16, 'AGG' : 1.25,
     'GGT' : 1.01, 'GGC' : 1.53, 'GGA' : 0.36, 'GGG' : 1.10}

_h_pylori_relative_adaptiveness = getRelativeAdaptivenessFromRSCUs(_h_pylori_rscu)
_ln_h_pylori_relative_adaptiveness = logRelativeAdaptiveness(_h_pylori_relative_adaptiveness)

_d_melanogaster_relative_adaptiveness = {'AAA':0.2734, 'AAC':1.0000, 'AAG':1.0000, 'AAT':0.5066, 'ACA':0.2820, 'ACC':1.0000, 'ACG':0.3460, 'ACT':0.3851, 'AGA':0.1277, 'AGC':0.6741, 'AGG':0.2028, 'AGT':0.2759, 'ATA':0.1695, 'ATC':1.0000, 'ATT':0.5920, 'CAA':0.3418, 'CAC':1.0000, 'CAG':1.0000, 'CAT':0.5262, 'CCA':0.5105, 'CCC':1.0000, 'CCG':0.4876, 'CCT':0.3101, 'CGA':0.2186, 'CGC':1.0000, 'CGG':0.1722, 'CGT':0.6485, 'CTA':0.1212, 'CTC':0.2980, 'CTG':1.0000, 'CTT':0.1960, 'GAA':0.3675, 'GAC':0.9898, 'GAG':1.0000, 'GAT':1.0000, 'GCA':0.2253, 'GCC':1.0000, 'GCG':0.2091, 'GCT':0.4556, 'GGA':0.6337, 'GGC':1.0000, 'GGG':0.0839, 'GGT':0.5962, 'GTA':0.1826, 'GTC':0.5982, 'GTG':1.0000, 'GTT':0.4320, 'TAC':1.0000, 'TAT':0.3843, 'TCA':0.2609, 'TCC':1.0000, 'TCG':0.6577, 'TCT':0.4154, 'TGC':1.0000, 'TGT':0.2955, 'TTA':0.0853, 'TTC':1.0000, 'TTG':0.3514, 'TTT':0.3370, 'ATG':1.0000, 'TGG':1.0000}
_ln_d_melanogaster_relative_adaptiveness = dict([(codon, math.log(ra)) for (codon,ra) in _d_melanogaster_relative_adaptiveness.items()])

def E_coli_CAI(gene):
	return getCAI(gene, _ln_e_coli_relative_adaptiveness)
def Yeast_CAI(gene):
	return getCAI(gene, _ln_s_cerevisiae_relative_adaptiveness)
def B_subtilis_CAI(gene):
	return getCAI(gene, _ln_b_subtilis_relative_adaptiveness)
def H_pylori_CAI(gene):
	return getCAI(gene, _ln_h_pylori_relative_adaptiveness)
def C_elegans_CAI(gene):
	return getCAI(gene, _ln_c_elegans_relative_adaptiveness)
def D_melanogaster_CAI(gene):
	return getCAI(gene, _ln_d_melanogaster_relative_adaptiveness)

def E_coli_Fop(gene):
	return getFop(gene, _e_coli_optimal_codons)
def B_subtilis_Fop(gene):
	return getFop(gene, _b_subtilis_optimal_codons)
def C_elegans_Fop(gene):
	return getFop(gene, _c_elegans_optimal_codons)
def Yeast_Fop(gene):
	return getFop(gene, _s_cerevisiae_optimal_codons)
def A_thaliana_Fop(gene):
	return getFop(gene, _a_thaliana_optimal_codons)
def Melanogaster_Fop(gene):
	return getFop(gene, _d_melanogaster_optimal_codons)
def Mouse_Fop(gene):
	return getFop(gene, _m_musculus_optimal_codons)
def Human_Fop(gene):
	return getFop(gene, _h_sapiens_optimal_codons)
def Human_Fop_Nucleotide_Corrected(gene, noncoding_nucleotides):
	return getFopNucleotideCorrected(gene, _h_sapiens_optimal_codons, noncoding_nucleotides)
def Sim_Fop(gene):
	return getFop(gene, _sim_optimal_codons)
def T_pallidum_Fop(gene):
	return getFop(gene, _t_pallidum_optimal_codons)
def H_pylori_Fop(gene):
	return getFop(gene, _h_pylori_optimal_codons)

def E_coli_Nop(gene):
	return getNumOptimalCodons(gene, _e_coli_optimal_codons)
def B_subtilis_Nop(gene):
	return getNumOptimalCodons(gene, _b_subtilis_optimal_codons)
def C_elegans_Nop(gene):
	return getNumOptimalCodons(gene, _c_elegans_optimal_codons)
def Yeast_Nop(gene):
	return getNumOptimalCodons(gene, _s_cerevisiae_optimal_codons)
def A_thaliana_Nop(gene):
	return getNumOptimalCodons(gene, _a_thaliana_optimal_codons)


def get_18_random_optimal_codons():
	candidate_codons = [(c,aa) for (c,aa) in translate._genetic_code.items() if (not 'U' in c and not aa in 'MW*')]
	target_aas = _aas
	done = False
	opt_codons = []
	for aa in target_aas:
		codons = [xc for (xc,xaa) in candidate_codons if xaa==aa]
		if len(codons)>1:
			opt_codons += random.sample(codons,1)
	return opt_codons

def getCAIFunction(master_species):
	if master_species == 'scer':
		fxn = Yeast_CAI
	elif master_species == 'cele' or master_species == 'celegans':
		fxn = C_elegans_CAI
	elif master_species == 'ecoli':
		fxn = E_coli_CAI
	elif master_species == 'dmel':
		fxn = D_melanogaster_CAI
	elif master_species == 'hpylori':
		fxn = H_pylori_CAI
	elif master_species == 'bsubtilis':
		fxn = B_subtilis_CAI
	else:
		raise KeyError, "Relative adaptivness values for species %s not found" % master_species
	return fxn

def getRelativeAdaptivenessValues(master_species):
	vals = None
	if master_species == 'scer':
		vals = _s_cerevisiae_relative_adaptiveness
	return vals

def getOptimalCodons(master_species):
	if master_species == 'human' or master_species == 'hsapiens':
		opt_codons = _h_sapiens_optimal_codons
	elif master_species == 'mouse':
		opt_codons = _m_musculus_optimal_codons
	elif master_species == 'macaque':
		opt_codons = _m_mulatta_optimal_codons
	elif master_species == 'scer':
		opt_codons = _s_cerevisiae_optimal_codons
	elif master_species == 'cele' or master_species == 'celegans':
		opt_codons = _c_elegans_optimal_codons
	elif master_species == 'melanogaster':
		opt_codons = _d_melanogaster_optimal_codons
	elif master_species == 'dmel':
		opt_codons = _d_melanogaster_optimal_codons
	elif master_species == 'sim':
		opt_codons = _sim_optimal_codons
	elif master_species == 'ecoli':
		opt_codons = _e_coli_optimal_codons
	elif master_species == 'tdenticola':
		opt_codons = _t_denticola_optimal_codons
	elif master_species == 'tpallidum':
		opt_codons = _t_pallidum_optimal_codons
	elif master_species == 'hpylori':
		opt_codons = _h_pylori_optimal_codons
	elif master_species == 'bburgdorferi':
		opt_codons = _b_burgdorferi_optimal_codons
	elif master_species == 'ccrescentus':
		opt_codons = _c_crescentus_optimal_codons
	elif master_species == 'mextorq':
		opt_codons = _m_extorquens_optimal_codons
	elif master_species == 'athaliana':
		opt_codons = _a_thaliana_optimal_codons_drummond
	elif master_species == 'bsubtilis':
		opt_codons = _b_subtilis_optimal_codons
	else:
		raise KeyError, "Optimal codons for species %s not found" % master_species
	opt_codons.sort()
	return opt_codons

if __name__=='__main__':
	test_getEmpiricalFrequencies()
	test_splitByConservation()
	test_randomizeConservationCategory()
	test_estimateSelectionCoefficients()

