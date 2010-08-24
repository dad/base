import sys, os, math, string, pickle
import biofile, translate, cai, util
from optparse import OptionParser

# Goal: use Yang and Nielsen 2008 to infer selection coefficients for codons
# Want to get both the mutation rates -- to get the stationary frequencies for each codon,
# and the 2NS_ij = F_i - F_j for each codon pair.
# Ultimate objective is to get

class Calculator:
	def __init__(self):
		self.codon_prob = None
		self.codon_prob_from_nucleotide = None
		self.cond_normalization = None
		self.codon_freq = None
		self.nucleotide_freq = None
		self.codon_prob_given_aa = None
		self.pseudocount = 0
		self.codon_syn_scores = None

	def initializeFromSequences(self, seqs, pseudocount):
		self.pseudocount = pseudocount

		self.codon_freq = dict([(c,self.pseudocount) for c in translate.DNACodons()])
		self.nucleotide_freq = dict([(nt,self.pseudocount) for nt in 'ATGC'])
		# The only reason to do this, versus just _accumulateCodonFrequencies(''.join(seqs)...), is a hedge against
		# the possibility that some sequences have bad lengths, and we want to keep everything in frame; this method
		# resets the reading frame every gene.
		for seq in seqs:
			if len(seq) % 3 == 0:
				self._accumulateCodonFrequencies(seq, self.codon_freq, self.nucleotide_freq)

		self._generateProbabilitiesFromFrequencies()
		self._generateScoresFromProbabilities()


	def initializeFromFrequencies(self, codon_freqs, nucleotide_freqs):
		self.codon_freq = dict(codon_freqs.items())
		self.nucleotide_freq = dict(nucleotide_freq.items())
		self._generateProbabilitiesFromFrequencies()
		self._generateScoresFromProbabilities()

	def initializeFromDictionary(self, codon_syn_scores):
		self.codon_syn_scores = dict(codon_syn_scores.items())

	def _generateProbabilitiesFromFrequencies(self):
		# Turn frequencies into probabilities
		total_codons = float(sum(self.codon_freq.values()))
		self.codon_prob = dict([(codon, self.codon_freq[codon]/total_codons) for codon in self.codon_freq.keys()])
		total_nucleotides = float(sum(self.nucleotide_freq.values()))
		self.codon_prob_from_nucleotide = {}
		for codon in self.codon_freq.keys():
			cprob = math.exp(sum([math.log(self.nucleotide_freq[nt]/total_nucleotides) for nt in codon]))
			self.codon_prob_from_nucleotide[codon] = cprob

		self.codon_prob_given_aa = {}
		for aa in translate.AAsAndStop():
			codons = translate.getCodons(aa, rna=False)
			#alt_codons[aa] = codons
			marginal_prob = sum([self.codon_prob[c] for c in codons])
			for codon in codons:
				# Compute the conditional probability of a codon, given that the amino acid is specified
				if marginal_prob > 0.0:
					self.codon_prob_given_aa[codon] = self.codon_prob[codon]/marginal_prob
				else:
					self.codon_prob_given_aa[codon] = 0.0

	def _accumulateCodonFrequencies(self, seq, codon_freq, nucleotide_freq):
		codons = cai.split(seq)
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

	def _generateScoresFromProbabilities(self):
		gc = translate.geneticCode()
		self.codon_syn_scores = {}
		for to_codon in self.codon_prob.keys():
			aa = gc[to_codon]
			sum_sc_i = 0.0
			# Go over all alternative codons and compute the average selection coefficient for moving from that codon to this one
			for from_codon in translate.getSynonyms(to_codon, rna=False):
				s_from_to = math.log(self.codon_prob[to_codon]/self.codon_prob[from_codon]) - math.log(self.codon_prob_from_nucleotide[to_codon]/self.codon_prob_from_nucleotide[from_codon])
				sum_sc_i += self.codon_prob_given_aa[from_codon]*s_from_to
			self.codon_syn_scores[to_codon] = sum_sc_i

	def getCodonProbabilitiesForMultipleSequences(self, seqs, pseudocount):
		return (codon_prob, codon_prob_from_nucleotide, codon_freq, nucleotide_freq)

	def getCodonProbabilityGivenAA(self, codon):
		return self.codon_prob_given_aa[codon]

	def getCodonProbability(self, codon):
		return self.codon_prob[codon]

	def getCodonProbabilityFromNucleotideComposition(self, codon):
		return self.codon_prob_from_nucleotide[codon]

	def getCodonRawFrequency(self, codon):
		return self.codon_freq[codon]

	def getNucleotideRawFrequency(self, nt):
		return self.nucleotide_freq[nt]

	def getSYangNielsen(self, seq):
		scores = [self.codon_syn_scores[to_codon] for to_codon in cai.split(seq)]
		return sum(scores)/len(scores)

	def getSYN(self, seq):
		return self.getSYangNielsen(seq)

	def getCodonSYNScores(self):
		return self.codon_syn_scores

	def __str__(self):
		s = ""
		if not self.nucleotide_freq is None:
			s += '\nnt\tnt.freq\n'
			for nt in 'ATGC':
				s += '{0}\t{1:d}\n'.format(nt, int(self.nucleotide_freq[nt]))
			s += '\n'
		if not self.codon_freq is None:
			if self.codon_prob is None:
				print "# generating probs from freqs"
				self._generateProbabilitiesFromFrequencies()
			if self.codon_syn_scores is None:
				print "# generating scores from probs"
				self._generateScoresFromProbabilities()
			s += 'aa\tcodon\tcodon.freq\tcodon.prob\tcodon.prob.from.nt\tcodon.cond.prob.given.aa\tsyn\n'
			for aa in translate.AAsAndStop():
				codons = translate.getCodons(aa, rna=False)
				for codon in codons:
					s += '{0}\t{1}\t{2:d}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\n'.format(aa, codon, int(self.codon_freq[codon]), self.codon_prob[codon], self.codon_prob_from_nucleotide[codon], self.codon_prob_given_aa[codon], self.codon_syn_scores[codon])
			s += '\n'
		return s

def getSYN(seq, syn_scores):
	scores = []
	for to_codon in cai.split(seq):
		try:
			scores.append(syn_scores[to_codon])
		except KeyError:
			continue
	if len(scores) > 0:
		res = sum(scores)/len(scores)
	else:
		res = None
	return res

if __name__=='__main__':
	parser = OptionParser(usage="%prog <genome filename> <genome dir> <format> <alignment filename> <tree or tree filename> [options]")
	parser.add_option("-o", "--out", dest="out_fname", type="string", default=None, help="output filename")
	parser.add_option("-f", "--format", dest="format", type="string", default="vanilla", help="format of ID in FASTA entry")
	parser.add_option("-d", "--dict-out", dest="score_dict_fname", type="string", default=None, help="score dictionary output filename")
	parser.add_option("-s", "--scores-out", dest="score_fname", type="string", default="vanilla", help="format of ID in FASTA entry")
	parser.add_option("-p", "--pseudocount", dest="pseudocount", type="float", default=0.0, help="pseudocount to be added to all frequencies")
	(options, args) = parser.parse_args()
	in_fname = args[0]

	info_outs = util.OutStreams(sys.stdout)
	data_outs = util.OutStreams()

	# Start up output
	if not options.out_fname is None:
		outf = file(options.out_fname, 'w')
		data_outs.addStream(outf)
	else:
		data_outs.addStream(sys.stdout)
	formatFxn = biofile.getIDFunction(options.format)
	cdna_dict = biofile.readFASTADict(in_fname, formatFxn)
	calc = Calculator()
	calc.initializeFromSequences(cdna_dict.values(), options.pseudocount)
	syn_dict = calc.getCodonSYNScores()
	data_outs.write("# Read {0}\n{1:d} sequences, {2:d} codons, {3:d} nucleotides\n".format(in_fname, len(cdna_dict.keys()), int(sum(calc.codon_freq.values())), int(sum(calc.nucleotide_freq.values()))))
	data_outs.write("# syn_scores = {0!s}\n".format(syn_dict))
	data_outs.write("{0!s}".format(calc))

	if not options.score_dict_fname is None:
		pickle.dump(syn_dict, file(options.score_dict_fname,'w'))

	if not options.score_fname is None:
		outf = file(options.score_fname, 'w')
		outf.write("orf\tsyn\n")
		orfs = cdna_dict.keys()
		n_written = 0
		for orf in sorted(orfs):
			seq = cdna_dict[orf]
			score = getSYN(seq, syn_dict)
			outf.write("{0}\t{1}\n".format(orf, util.formatNA(score,'{0:.5f}')))
			n_written += 1
		info_outs.write("# Wrote {0} lines to {1}\n".format(n_written, options.score_fname))


