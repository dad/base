# Goal: use Yang and Nielsen 2008 to infer selection coefficients for codons
# Want to get both the mutation rates -- to get the stationary frequencies for each codon,
# and the 2NS_ij = F_i - F_j for each codon pair.
# Ultimate objective is to get

class SYNCalculator:
	def __init__(self):
		self.codon_prob = None
		self.codon_prob_from_nucleotide = None
		self.cond_normalization = None
		self.codon_freq = None
		self.nucleotide_freq = None
		self.codon_prob_given_aa = None
		self.pseudocount = 0
	
	def initializeFromSequences(self, seqs, pseudocount):
		self.pseudocount = pseudocount
		
		self.codon_freq = dict([(c,self.pseudocount) for c in translate.AADNACodons()])
		self.nucleotide_freq = dict([(nt,self.pseudocount) for nt in 'ATGC'])
		# The only reason to do this, versus just _accumulateCodonFrequencies(''.join(seqs)...), is a hedge against
		# the possibility that some sequences have bad lengths, and we want to keep everything in frame; this method
		# resets the reading frame every gene.
		for seq in seqs:
			if len(seq) % 3 == 0:
				self._accumulateCodonFrequencies(seq, self.codon_freq, self.nucleotide_freq)
		
		# Turn frequencies into probabilities
		total_codons = float(sum(self.codon_freq.values()))
		self.codon_prob = dict([(codon, self.codon_freq[codon]/total_codons) for codon in self.codon_freq.keys()])
		total_nucleotides = float(sum(self.nucleotide_freq.values()))
		self.codon_prob_from_nucleotide = {}
		for codon in self.codon_freq.keys():
			cprob = math.exp(sum([math.log(self.nucleotide_freq[nt]/total_nucleotides) for nt in codon]))
			self.codon_prob_from_nucleotide[codon] = cprob

		self.codon_prob_given_aa = {}
		for aa in translate.AAs():
			codons = translate.getCodonsForAA(aa, rna=False)
			alt_codons[aa] = codons
			marginal_prob = sum([codon_prob[c] for c in codons])
			for codon in codons:
				# Compute the conditional probability of a codon, given that the amino acid is specified
				self.codon_prob_given_aa[codon] = self.codon_prob[codon]/marginal_prob
	
	def initializeFromFrequencies(self, codon_freqs, nucleotide_freqs):
		self.codon_freq = dict(codon_freqs.items())
		self.nucleotide_freq = dict(nucleotide_freq.items())
	
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
		return (codon_prob, codon_prob_from_nucleotide, codon_freq, nucleotide_freq)

	def getCodonProbabilityGivenAA(self, codon):
		return self.codon_prob_given_aa[codon]

	def getCodonProbability(self, codon):
	
	def getCodonProbabilityFromNucleotideComposition(self, codon):
		return self.codon_prob_from_nucleotide[codon]
	
	def getCodonRawFrequency(self, codon):
		return self.codon_freq[codon]
		
	def getNucleotideRawFrequency(self, nt):
		return self.nucleotide_freq[nt]
		
	def getSYangNielsen(self, seq):
		gc = translate.geneticCode()
		# Compute normalizations for conditional probabilities of a codon given an amino acid
		sum_sc = 0.0
		gene_codons = split(seq)
		for to_codon in gene_codons:
			
			aa = gc[to_codon]
			if not aa == '*':
				sum_sc_i = 0.0
				# Go over all alternative codons and compute the average selection coefficient for moving from that codon to this one
				for from_codon in alt_codons[aa]:
					s_from_to = math.log(self.codon_prob[to_codon]/self.codon_prob[from_codon]) - math.log(self.codon_prob_from_nt[to_codon]/self.codon_prob_from_nt[from_codon])
					sum_sc_i += self.codon_prob_given_aa[from_codon]*s_from_to
				sum_sc += sum_sc_i
		return sum_sc/len(gene_codons)
		

if __name__=='__main__':
	
	