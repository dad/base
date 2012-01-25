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
	
	def test_YDR114C(self):
		pp = protprop.ProteinProperties()
		ydr = "MPLFARLCQPQSRRMFSSISSFSALSVLRPQTGMLLNSSPLKTPSFTPLGFGLIGQRRWKSRGNTYQPSTLKRKRTFGFLARAKSKQGSKILKRRKLKGRWFLSH"
		ydr = "MPLFARLCQPQSRRMFSSISSFSALSVLRPQTGMLLNSSPLKTPSFTPLGFGLIGQRRWKSRGNTYQPSTLKRKRTFGFLARAKSKQGSKILKRRKLKGRWFLSH"
		pI = pp.getIsoelectricPoint(ydr, tolerance=0.1)
		print pI
		charge = pp.getCharge(ydr, pH=7.2)
		print charge


if __name__=='__main__':
	unittest.main(verbosity=2)
