#! python

import sys, os, math, random, stats, unittest
import util, biofile, protprop

class ProtPropBasic(unittest.TestCase):
	def test_pka(self):
		pp = protprop.ProteinProperties()
		self.assertTrue(pp.getPKa('KKKKKKKKKKLL')==pp.getPKa('KKKKKKKKKK'))
		

if __name__=='__main__':
	unittest.main(verbosity=2)
