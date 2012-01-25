#! python

import sys, os, math, random, stats, unittest
import util, biofile, protprop

class ProtPropBasic(unittest.TestCase):
	def setUp(self):
		self.pp = protprop.ProteinProperties()
		
	def test_pI(self):
		pp = self.pp
		# Same isoelectric point if you add non-ionizable amino acids
		self.assertTrue(pp.getIsoelectricPoint('KKKKKKKKKKLL')==pp.getIsoelectricPoint('KKKKKKKKKK'))
		# net charge at isoelectric point is zero
		tol = 1e-5
		pI = pp.getIsoelectricPoint('KKKKKKKKKKLL', tolerance=tol)
		self.assertTrue(abs(pp.getCharge('KKKKKKKKKKLL',pI))<=tol)
	
	def test_YDR114C(self):
		pp = self.pp
		# This sequence takes an eternity to search for pI; test that maximum number of iterations prevents that.
		pp = protprop.ProteinProperties()
		ydr = "MPLFARLCQPQSRRMFSSISSFSALSVLRPQTGMLLNSSPLKTPSFTPLGFGLIGQRRWKSRGNTYQPSTLKRKRTFGFLARAKSKQGSKILKRRKLKGRWFLSH"
		pI = pp.getIsoelectricPoint(ydr)
		charge = pp.getCharge(ydr, pH=7.2)
	
	def test_length(self):
		# Lengths of proteins with gaps or stop codons
		pp = self.pp
		L = pp.getLength("KKKKKKKKKK")
		self.assertTrue(L==10)
		L = pp.getLength("KKKKKKKKKK*")
		self.assertTrue(L==10)
		L = pp.getLength("KKKKKKK*KKK")
		self.assertTrue(L==7)
		L = pp.getLength("KKKKK--KKK*")
		self.assertTrue(L==8)
	

if __name__=='__main__':
	unittest.main(verbosity=2)
