#! python

import sys, os, math, random, stats, unittest
import util, biofile, protprop

class ProtPropBasic(unittest.TestCase):
	def test_pI(self):
		pp = protprop.ProteinProperties()
		# Same isoelectric point if you add non-ionizable amino acids
		self.assertTrue(pp.getIsoelectricPoint('KKKKKKKKKKLL')==pp.getIsoelectricPoint('KKKKKKKKKK'))
		# net charge at isoelectric point is zero
		tol = 1e-5
		pI = pp.getIsoelectricPoint('KKKKKKKKKKLL', tolerance=tol)
		self.assertTrue(abs(pp.getCharge('KKKKKKKKKKLL',pI))<=tol)
		
		

if __name__=='__main__':
	unittest.main(verbosity=2)
