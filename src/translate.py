#!python
"""Module for performing various basic biology operations.

Seed version by Jesse Bloom, 2004.
Expanded by D. Allan Drummond, 2004-2014."""

import os, sys, string, math, random
#-----------------------------------------------------------------------------------
class BioUtilsError(Exception):
    """Error using one of the bio utils."""

#---------------------------------------------------------------------------
# The universal genetic code
_genetic_code = {
		'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L',
		'CTA':'L', 'CTG':'L', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'GTT':'V',
		'GTC':'V', 'GTA':'V', 'GTG':'V', 'TCT':'S', 'TCC':'S', 'TCA':'S',
		'TCG':'S', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'ACT':'T',
		'ACC':'T', 'ACA':'T', 'ACG':'T', 'GCT':'A', 'GCC':'A', 'GCA':'A',
		'GCG':'A', 'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
		'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'AAT':'N', 'AAC':'N',
		'AAA':'K', 'AAG':'K', 'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
		'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W', 'CGT':'R',
		'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGT':'S', 'AGC':'S', 'AGA':'R',
		'AGG':'R', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
		'UUU':'F', 'UUC':'F', 'UUA':'L', 'UUG':'L', 'CUU':'L', 'CUC':'L',
		'CUA':'L', 'CUG':'L', 'AUU':'I', 'AUC':'I', 'AUA':'I', 'AUG':'M', 'GUU':'V',
		'GUC':'V', 'GUA':'V', 'GUG':'V', 'UCU':'S', 'UCC':'S', 'UCA':'S',
		'UCG':'S', 'CCU':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'ACU':'T',
		'ACC':'T', 'ACA':'T', 'ACG':'T', 'GCU':'A', 'GCC':'A', 'GCA':'A',
		'GCG':'A', 'UAU':'Y', 'UAC':'Y', 'UAA':'*', 'UAG':'*',
		'CAU':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'AAU':'N', 'AAC':'N',
		'AAA':'K', 'AAG':'K', 'GAU':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
		'UGU':'C', 'UGC':'C', 'UGA':'*', 'UGG':'W', 'CGU':'R',
		'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGU':'S', 'AGC':'S', 'AGA':'R',
		'AGG':'R', 'GGU':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}

_rna_codons = [x for x in _genetic_code.keys() if not 'T' in x]
_dna_codons = [x for x in _genetic_code.keys() if not 'U' in x]
_aa_dna_codons = [x for x in _dna_codons if not _genetic_code[x] == '*']
_aa_rna_codons = [x for x in _rna_codons if not _genetic_code[x] == '*']

# Yeast mitochondrial code, from http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG3
#     MITO              NUCLEAR
# AUA    Met  M          Ile  I
# CUU    Thr  T          Leu  L
# CUC    Thr  T          Leu  L
# CUA    Thr  T          Leu  L
# CUG    Thr  T          Leu  L
# UGA    Trp  W          Ter  *
#
# CGA    absent          Arg  R
# CGC    absent          Arg  R
_scer_mito_code = dict(_genetic_code.items())
_scer_mito_code_alts = {'AUA':'M', 'CUU':'T', 'CUC':'T', 'CUA':'T', 'CUG':'T', 'UGA':'W'}
for (cod,aa) in _scer_mito_code_alts.items():
	_scer_mito_code[cod] = aa
	_scer_mito_code[cod.replace('U','T')] = aa

def geneticCode(rna, code=None):
	res = _genetic_code
	if not code is None:
		if code == 'scer-mito':
			res = _scer_mito_code
	if rna:
		res = dict([(k,v) for (k,v) in res.items() if not 'T' in k])
	else:
		res = dict([(k,v) for (k,v) in res.items() if not 'U' in k])
	return res


def codons(seq):
	n = len(seq)
	over_nts = n % 3
	max_aas = (n-over_nts)/3
	for ci in xrange(max_aas):
		yield seq[3*ci : 3*(ci+1)]

#--------------------------------------------------------------------------------
def translate(seq, genetic_code=_genetic_code):
	"""Translates a gene sequence to a protein sequence.

	'seq' is the gene sequence to be translated. It can begin with any codon
	(does not have to be ATG), but it must be of the proper length
	(a multiple of 3).  It can contain a trailing stop codon, but if it
	contains stop codons before the end of the sequence it is not
	translated.  If the translation is successful, returns a string
	corresponding to the protein sequence.  If the translation
	fails, returns 'None'."""
	seq = seq.upper()
	if len(seq) % 3 != 0:
		return None # gene length not a multiple of three
	prot_length = len(seq) / 3
	prot = ""
	stop_codon = False
	for codon in codons(seq):
		try:
			aa = genetic_code[codon]
			if aa == '*':
				stop_codon = True
				break
			prot += aa
		except KeyError:
			prot = None
			break
	if stop_codon:
		if len(prot)<(prot_length-1):
			prot = None
	else:
		if not prot is None and len(prot)<prot_length:
			prot = None
	return prot

def Translate(seq):
	return translate(seq)

def translateRaw(seq, genetic_code=_genetic_code, bad_aa = 'x'):
	"""Translates a nucleotide sequence to a protein sequence.

	'seq' is the gene sequence to be translated. It can begin with any codon
	(does not have to be ATG). Length need not be a multiple of 3 nucleotides.
	'bad_aa' is the character used to indicate any codon not in the standard
	code. Internal stop codons are translated as '*'.

	If the translation is successful, returns a string corresponding to the
	protein sequence plus the stop codon if present."""
	seq = seq.upper()
	prot = ''.join([genetic_code.get(codon,bad_aa) for codon in codons(seq)])
	return prot

def reverseTranslate(prot, rna=False, bad_codon ='---'):
	gene = ""
	rev_code = dict([(aa,codon) for (codon,aa) in _genetic_code.items() if not 'U' in codon])
	for aa in prot:
		try:
			gene += rev_code[aa]
		except KeyError, ke:
			gene += bad_codon
	assert(len(gene)==3*len(prot))
	return gene

def randomReverseTranslate(prot, rna=False, bad_codon='---'):
	"""
	Translates a protein into codons, choosing codons randomly instead of
	deterministically.
	"""
	letters = 'ACDEFGHIKLMNPQRSTVWXY-'
	codon_choices = {}
	for aa in letters:
		codon_choices[aa] = get_codons_for_aa(aa, rna=rna)
		if len(codon_choices[aa])==0:
			codon_choices[aa] = [bad_codon]
	gene = ''
	#print prot
	for aa in prot:
		codon = random.choice(codon_choices[aa])
		gene += codon
	#print sequenceDiffs(prot, translate.TranslateRaw(gene))
	return gene


_three_letter_codes = {
	'A':'Ala', 'C':'Cys', 'D':'Asp', 'E':'Glu', 'F':'Phe', 'G':'Gly', 'H':'His', 'I':'Ile', 'K':'Lys', 'L':'Leu', 'M':'Met', 'N':'Asn', 'P':'Pro', 'Q':'Gln', 'R':'Arg', 'S':'Ser', 'T':'Thr', 'V':'Val', 'W':'Trp', 'Y':'Tyr'}

def threeLetterCodes():
	return _three_letter_codes

#---------------------------------------------------------------------------
def codonToAA(codon):
	"""Returns one-letter amino acid for codon.

	Argument is three-letter string 'codon'.
	'codon' should be upper case.
	Raises error if invalid codon.."""
	try:
		return _genetic_code[codon]
	except KeyError:
		raise BioUtilsError, "Invalid codon '%s'." % codon

def Codon_to_AA(codon):
	return codonToAA(codon)

#---------------------------------------------------------------------------------
def sequenceIdentity(aligned_seq1, aligned_seq2):
	num_identical = 0
	num_aligned = 0
	for i in range(min(len(aligned_seq1), len(aligned_seq2))):
		aa1 = aligned_seq1[i]
		aa2 = aligned_seq2[i]
		if aa1 != '-' and aa2 != '-':
			num_aligned += 1
			if aa1 == aa2:
				num_identical += 1
	seq_identity = 0.0
	if num_aligned > 0:
		seq_identity = float(num_identical)/num_aligned
	return seq_identity, num_identical, num_aligned

class SiteConsensus(object):
	def __init__(self, gap_threshold=0.0, no_consensus_char='.'):
		self.amino_acid = None
		self.proportion = None
		self.frequency = None
		self.gap_frequency = None
		self.gap_proportion = None
		self.entropy = None
		self._no_consensus_char = no_consensus_char
		#self._proportion_threshold = threshold
		self._gap_threshold = gap_threshold # Proportion of gaps above which there is no consensus

	def computeFrom(self, aa_list):
		if len(aa_list)>0:
			counts = [(aa_list.count(aa),aa) for aa in AAs()]
			counts.sort(reverse=True)
			ungapped_len = float(sum(counts))
			props = [ct/ungapped_len for (ct,aa) in counts if ct>0]
			self.entropy = sum([-p*math.log(p,20.0) for p in props])
			self.frequency = counts[0][0]
			self.gap_frequency = len(aa_list)-ungapped_len
			self.gap_proportion = self.gap_frequency/len(aa_list)
			self.proportion = self.frequency/ungapped_len
			self.amino_acid = counts[0][1]
			if self.gap_proportion >= self._gap_threshold:
				self.amino_acid = self._no_consensus_char

class ConsensusSequence(object):
	def __init__(self, threshold=0.0, no_consensus_char='.'):
		self._seq = ''
		self._len = -1
		self._no_consensus_char = no_consensus_char
		self._threshold = threshold
		self._cons_list = None

	def computeFrom(self, seq_list):
		n = len(seq_list[0])
		self._len = n
		self._cons_list = []
		self._seq = ''
		for ai in range(n):
			cons = SiteConsensus(self._threshold, self._no_consensus_char)
			cons.computeFrom([s[ai] for s in seq_list])
			self._cons_list.append(cons)
			self._seq += cons.amino_acid

	def identityTo(self, seq):
		assert len(seq)==len(self._seq)
		num_id = 0
		num_tot = 0
		for (ai,s) in enumerate(seq):
			cons = self._cons_list[ai]
			if self._seq[ai] != self._no_consensus_char:
				num_tot += 1
				if s==cons.amino_acid:
					num_id += 1
		return float(num_id)/num_tot

	def __getitem__(self, ind):
		return self._cons_list[ind]

	@property
	def noconsensus(self):
		return self._no_consensus_char

	def __str__(self):
		return self._seq

#---------------------------------------------------------------------------------
_dna_complement = {'A':'T','T':'A','U':'A','G':'C','C':'G','a':'t','t':'a','u':'a','g':'c','c':'g'}
_rna_complement = {'A':'U','U':'A','G':'C','C':'G','a':'u','u':'a','g':'c','c':'g'}
def complement(a, rna=True):
	if rna:
		return _rna_complement[a]
	return _dna_complement[a]
#---------------------------------------------------------------------------------
def reverse_complement(seq):
	return reverseComplement(seq)

def reverseComplement(seq):
	rna = 'U' in seq.upper()
	rc = [complement(x,rna) for x in seq]
	rc.reverse()
	return ''.join(rc)

def get_codons_for_aa(aa, rna=True):
	return getCodonsForAA(aa, rna)

def getCodonsForAA(aa, rna=True):
	gc = _genetic_code
	aa_codons = []
	if rna:
		aa_codons = [c for c in gc.keys() if gc[c] == aa and not 'T' in c]
	else:
		aa_codons = [c for c in gc.keys() if gc[c] == aa and not 'U' in c]
	return aa_codons

def getCodons(aa_or_stop, rna=True):
	gc = _genetic_code
	codons = []
	if rna:
		codons = [c for c in gc.keys() if gc[c] == aa_or_stop and not 'T' in c]
	else:
		codons = [c for c in gc.keys() if gc[c] == aa_or_stop and not 'U' in c]
	return sorted(codons)

def getSynonyms(codon, rna=True):
	gc = _genetic_code
	if rna:
		synonyms = [c for c in gc.keys() if gc[c] == gc[codon] and not 'T' in c]
	else:
		synonyms = [c for c in gc.keys() if gc[c] == gc[codon] and not 'U' in c]
	return synonyms

def degenerateAAs():
	return "ACDEFGHIKLNPQRSTVY"

def AAs():
	return "ACDEFGHIKLMNPQRSTVWY"

def AAsAndStop():
	return AAs()+'*'

def DNACodons():
	return _dna_codons

def RNACodons():
	return _rna_codons

def AADNACodons():
	return _aa_dna_codons

def AARNACodons():
	return _aa_rna_codons


if __name__=='__main__':
	__test_reverseTranslate()
	__test_randomReverseTranslate()
