#!/usr/bin/python
# Begin my_bio_utils.py
"""Module for performing various basic biology operations.

Original version by Jesse Bloom, 2004.
Expanded by D. Allan Drummond, 2004-2010."""
#
import re, os, sys, string, math, random
#-----------------------------------------------------------------------------------
class BioUtilsError(Exception):
    """Error using one of the bio utils."""

#--------------------------------------------------------------------------------
def translate(seq):
	"""Translates a gene sequence to a protein sequence.

	'seq' is the gene sequence to be translated. It can begin with any codon
	(does not have to be ATG), but it must be of the proper length
	(a multiple of 3).  It can contain a trailing stop codon, but if it
	contains stop codons before the end of the sequence it is not
	translated.  If the translation is successful, returns a string
	corresponding to the protein sequence.  If the translation
	fails, returns 'None'."""
	if len(seq) % 3 != 0:
		return None # gene length not a multiple of three
	prot_length = len(seq) / 3
	prot = ""
	for i in range(prot_length - 1):
		codon = seq[3 * i : 3 * (i + 1)]
		try:
			aa = codonToAA(codon)
		except BioUtilsError: # unrecognized codon
			return None
		if aa == '*': # premature stop codon
			return None
		prot += aa
	# last codon, might be stop codon
	codon = seq[3 * (prot_length - 1) : 3 * prot_length]
	aa = codonToAA(codon)
	if aa != '*':
		prot += aa
	assert len(prot) in [prot_length, prot_length - 1]
	return prot

def Translate(seq):
	return translate(seq)

def translateRaw(seq, bad_aa = 'x'):
	"""Translates a nucleotide sequence to a protein sequence.

	'seq' is the gene sequence to be translated. It can begin with any codon
	(does not have to be ATG), and the length must be at least 3 nucleotides.
	'bad_aa' is the character used to indicate any codon not in the standard
	code.

	If the translation is successful, returns a string corresponding to the
	protein sequence plus stop codons and ."""
	prot = ""
	max_aas = int(math.floor(len(seq)/3))
	codon = seq[0:3]
	for i in range(max_aas):
		codon = seq[3 * i : 3 * (i + 1)]
		try:
			aa = codonToAA(codon)
		except BioUtilsError: # unrecognized codon
			aa = bad_aa
		prot += aa
	return prot

def TranslateRaw(seq, bad_aa = 'x'):
	return translateRaw(seq, bad_aa)

def reverseTranslate(prot, bad_codon ='xxx'):
	gene = ""
	rev_code = dict([(aa,codon) for (codon,aa) in _genetic_code.items() if not 'U' in codon])
	for aa in prot:
		try:
			gene += rev_code[aa]
		except KeyError, ke:
			gene += bad_codon
	assert(len(gene)==3*len(prot))
	return gene

def ReverseTranslate(prot, bad_codon ='xxx'):
	return reverseTranslate(prot, bad_codon)


# Test the reverse translator
def __test_reverseTranslate():
	N = 1000
	aas = 'ACDEFGHIKLMNPQRSTVWY'
	for i in range(N):
		prot = ''.join([random.choice(aas) for xi in range(100)])
		gene = reverseTranslate(prot)
		newprot = translate(gene)
		assert(prot == newprot)

def randomReverseTranslate(prot, rna=False):
	"""
	Translates a protein into codons, choosing codons randomly instead of
	deterministically.
	"""
	letters = 'ACDEFGHIKLMNPQRSTVWXY-'
	codon_choices = {}
	for aa in letters:
		codon_choices[aa] = get_codons_for_aa(aa, rna=rna)
		if len(codon_choices[aa])==0:
			codon_choices[aa] = ['---']
	gene = ''
	#print prot
	for aa in prot:
		codon = random.choice(codon_choices[aa])
		gene += codon
	#print sequenceDiffs(prot, translate.TranslateRaw(gene))
	return gene

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

def geneticCode():
	return _genetic_code

def DNACodons():
	return _dna_codons

def RNACodons():
	return _rna_codons

def AADNACodons():
	return _aa_dna_codons

def AARNACodons():
	return _aa_rna_codons


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
#---------------------------------------------------------------------------------
def complement(a):
	a = a.upper()
	if a == "A":
		return "T"
	elif a == "T" or a == "U":
		return "A"
	elif a == "C":
		return "G"
	elif a == "G":
		return "C"
	# If we don't know the complement, return the base
	return a
#---------------------------------------------------------------------------------
def reverse_complement(seq):
	return reverseComplement(seq)

def reverseComplement(seq):
	rc = [x for x in seq]
	rc.reverse()
	for i in range(len(rc)):
		rc[i] = complement(rc[i])
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
	return codons

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

# End translate.py
