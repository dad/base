#! python

import sys, os, math, random, stats, unittest
import util, biofile, protprop, translate

class test001(unittest.TestCase):
	def setUp(self):
		self.pp = protprop.ProteinProperties()
		
	def test_pI(self):
		"""Same isoelectric point if you add non-ionizable amino acids"""
		pp = self.pp
		self.assertTrue(pp.getIsoelectricPoint('KKKKKKKKKKLL')==pp.getIsoelectricPoint('KKKKKKKKKK'))
		# net charge at isoelectric point is zero
		tol = 1e-5
		pI = pp.getIsoelectricPoint('KKKKKKKKKKLL', tolerance=tol)
		self.assertAlmostEqual(pp.getCharge('KKKKKKKKKKLL',pI),0.0,places=4)
	
	def test_YDR114C(self):
		"""Maximum number of iterations prevents overlong searches"""
		pp = self.pp
		# This sequence takes an eternity to search for pI.
		pp = protprop.ProteinProperties()
		ydr = "MPLFARLCQPQSRRMFSSISSFSALSVLRPQTGMLLNSSPLKTPSFTPLGFGLIGQRRWKSRGNTYQPSTLKRKRTFGFLARAKSKQGSKILKRRKLKGRWFLSH"
		pI = pp.getIsoelectricPoint(ydr)
		charge = pp.getCharge(ydr, pH=7.2)
	
	def test_length(self):
		"""Lengths of proteins with gaps or stop codons"""
		pp = self.pp
		L = pp.getLength("KKKKKKKKKK")
		self.assertTrue(L==10)
		L = pp.getLength("KKKKKKKKKK*")
		self.assertTrue(L==10)
		L = pp.getLength("KKKKKKK*KKK")
		self.assertTrue(L==7)
		L = pp.getLength("KKKKK--KKK*")
		self.assertTrue(L==8)

class test002(unittest.TestCase):
	def test_run(self):
		"""Composition"""
		comp = protprop.Composition()
		fname = "tmp_composition.txt"
		inf = file(fname, 'w')
		inf.write("aa\tproportion\n")
		for aa in translate.AAs():
			inf.write("{}\t{}\n".format(aa, 1.0/20))
		inf.close()
		inf = file(fname, 'r')
		comp.read(inf)
		self.assertAlmostEqual(comp['A'], 1.0/20)
		inf.close()
		os.remove(fname)

def genAlternatingMotif(aas_list, num_samples_list):
	"""Generate a motif that has"""
	res = ''
	for (i,n) in num_samples_list:
		res += ''.join(stats.sample_wr(aas_list[i],n))
	return res

class test003(unittest.TestCase):
	def test_run(self):
		"""Composition and motifs"""
		comp = protprop.Composition()
		pp = protprop.ProteinProperties()
		for xi in range(5):
			seq = genAlternatingMotif(['FY','P','NQ'], [(0,2),(1,1),(2,2),(0,2),(1,1),(2,2)])
			self.assertTrue(pp.count(seq, 'FY')==4)
			#self.assertTrue(comp.motif(seq)=='aabccddeff')

if __name__=='__main__':
	unittest.main(verbosity=2)
