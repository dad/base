#!/usr/bin/python
"""Module for testing codon bias module "cai".
"""

import re, os, sys, string, math, random, unittest
import translate

class test_reverse_translate(unittest.TestCase):
	def test_reverseTranslate(self):
		N = 1000
		aas = 'ACDEFGHIKLMNPQRSTVWY'
		for i in range(N):
			prot = ''.join([random.choice(aas) for xi in range(100)])
			gene = translate.reverseTranslate(prot)
			newprot = translate.translate(gene)
			assert(prot == newprot)

	def test_randomReverseTranslate(self):
		N = 1000
		aas = 'ACDEFGHIKLMNPQRSTVWY'
		for i in range(N):
			prot = ''.join([random.choice(aas) for xi in range(100)])
			gene = translate.randomReverseTranslate(prot)
			newprot = translate.translate(gene)
			assert(prot == newprot)


if __name__=='__main__':
	unittest.main(verbosity=2)

