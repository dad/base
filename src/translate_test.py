#! python

import os, sys, string, math, random, unittest
import translate

class Tester(unittest.TestCase):
	def test001(self):
		"""reverse translate"""
		N = 1000
		aas = 'ACDEFGHIKLMNPQRSTVWY'
		for i in range(N):
			prot = ''.join([random.choice(aas) for xi in range(100)])
			gene = translate.reverseTranslate(prot)
			newprot = translate.translate(gene)
			self.assertTrue(prot == newprot)

	def test002(self):
		"""reverse complement"""
		s = 'ATGCatgc'
		self.assertTrue(translate.reverseComplement(s)=='gcatGCAT')
		self.assertTrue(translate.reverseComplement(translate.reverseComplement(s))==s)

	def test003(self):
		"""translation"""
		s = 'ATGCatTCT'
		#print translate.translate(s)
		self.assertTrue(translate.translate(s)=='MHS')

	def test004(self):
		"""translation with problems"""
		s = 'ATGCatTCTNNNTAAAGA'
		#print translate.translate(s)
		self.assertTrue(translate.translate(s) is None)
		#print translate.translateRaw(s,bad_aa='@')
		self.assertTrue(translate.translateRaw(s,bad_aa='@') == 'MHS@*R')

	def test005(self):
		"""codons"""
		s = 'ATGCatTCTNNNTAAAGA'
		c = [cod for cod in translate.codons(s)]
		self.assertTrue(c[0]=='ATG')
		self.assertTrue(c[1]=='Cat')
		self.assertTrue(c[-1]=='AGA')

	def test006(self):
		"""zero-length codons"""
		s = 'AT'
		c = [cod for cod in translate.codons(s)]
		self.assertTrue(len(c)==0)

	def test007(self):
		"""zero-length protein"""
		s = 'AT'
		self.assertTrue(translate.translate(s) is None)
		self.assertTrue(translate.translateRaw(s) == '')

	def test008(self):
		"""zero-length nucleotides"""
		s = ''
		c = [cod for cod in translate.codons(s)]
		self.assertTrue(len(c)==0)
		self.assertTrue(translate.translate(s) == '')
		self.assertTrue(translate.translateRaw(s) == '')

	def test009(self):
		"""odd-length coding sequence"""
		s = 'TCTCGTAAGTACGCAGC'
		c = [cod for cod in translate.codons(s)]
		self.assertTrue(len(c)==5)
		self.assertTrue(c[-1]=='GCA')
		self.assertTrue(translate.translate(s) is None)
		self.assertTrue(translate.translateRaw(s) == 'SRKYA')

class TestCompare(unittest.TestCase):
	def test_comp(self):
		"""compare"""
		a = 'ACDEFGHIKLMNPQRSTVWY'
		b = 'ACDEFGHIKLMNPQRSTVWY'
		c = 'ACDEFGHIKLMNPQRSTVWq'
		d = 'ACDEFGHIKLMNPQRS-VWY'
		simab = translate.compare(a, b)
		self.assertTrue(simab.identity == 1.0)
		simac = translate.compare(a, c)
		self.assertTrue(simac.identity == 19/20.0)
		simad = translate.compare(a, d)
		print simad
		self.assertTrue(simad.identity == 1.0)
		self.assertTrue(simad.fraction_aligned_x == 19/20.0)
		self.assertTrue(simad.fraction_aligned_y == 1.0)


if __name__=='__main__':
	unittest.main(verbosity=2)

