#! python

import sys, os, math, random, stats
from optparse import OptionParser
import util, biofile

# pK values from http://helixweb.nih.gov/emboss/html/iep.html
# pI code from Peter Collingridge, http://www.petercollingridge.co.uk/sites/files/peter/predictPI.txt

class ProteinProperties(object):
	def __init__(self):
		self.pKa     = {'D':3.9, 'E':4.3, 'H':6.1, 'C':8.3, 'Y':10.1, 'K':10.5, 'R':12, 'N-term':8, 'C-term':3.1}
		self.charges = {'D':-1,  'E':-1,  'H':1,  'C':-1,  'Y':-1,   'K':1,    'R':1,  'N-term':1, 'C-term':-1}

	def _aminoAcidCharge(self, amino_acid, pH):
		proportion = 1 / (1 + 10**(pH - self.pKa[amino_acid]))
		res = None
		if self.charges[amino_acid] == 1:
			res = proportion
		else:
			res = proportion-1.0 # more clearly, -1 * (1-proportion)
		return res

	def getCharge(self, sequence, pH):
		protein_charge = self._aminoAcidCharge('N-term', pH) + self._aminoAcidCharge('C-term', pH)
		for amino_acid in self.pKa.keys():
			protein_charge += sequence.count(amino_acid) * self._aminoAcidCharge(amino_acid, pH)
		return protein_charge

	def getIsoelectricPoint(self, sequence, tolerance=1e-4):
		min_pH, max_pH = 3, 13 
		done = False
		pI = None
		n_iterations = 0
		# Binary search for the pH at which the net charge is within the specified tolerance around zero.
		while not done:
			mid_pH = 0.5 * (max_pH + min_pH)
			protein_charge = self.getCharge(sequence, mid_pH)
			if protein_charge > tolerance:
				min_pH = mid_pH
			elif protein_charge < -tolerance:
				max_pH = mid_pH
			else:
				pI = mid_pH
				done = True
			n_iterations += 1
		return pI

if __name__=='__main__':
	parser = OptionParser(usage="%prog [options] <sequence>")
	(options, args) = parser.parse_args()
	seq = args[0]
	pp = ProteinProperties()
	print "length =", len(seq)
	pI = pp.getIsoelectricPoint(seq, tolerance=1e-4)
	print "pI = {0}".format(pI) 
	pH = 7.2
	print "charge at pH {0:1.1f} = {1:1.2f}".format(pH, pp.getCharge(seq, pH))
	print "charge at pH=pI = {:1.2f}".format(pp.getCharge(seq, pI))

	