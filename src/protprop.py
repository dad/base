#! python

import sys, os, math, random, stats
import util, biofile

# pK values from http://helixweb.nih.gov/emboss/html/iep.html

class ProteinProperties(object):
	def __init__(self):
		self.pkas = {'C':8.5, 'D':3.9, 'E':4.1, 'H':6.5, 'K':10.8, 'R':12.5, 'Y':10.1, 'N-term':8.6, 'C-term':3.6}
	
	def getPKa(self, seq):
		pkas = [self.pkas['N-term'], self.pkas['C-term']]
		for aa in seq:
			try:
				pkas.append(self.pkas[aa])
			except KeyError:
				pass
		return stats.mean(pkas)

	def getCharge(self, seq, pH=7.2):
		q_sum = 0.0
		for aa in seq:
			try:
				pka = self.pkas[aa]
				r = 10**(pH-pka)
				q_sum += r/(1.0+r)
			except KeyError:
				pass
		return q_sum
			