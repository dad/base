#! python

import os, sys, math, unittest
import orthodbutil

class test001(unittest.TestCase):
	def test_run(self):
		"""Handle quoted names"""
		header = "7244:002ce8 FBpp0236088 gene=FBgn0208790 orthodb8_OG=EOG8MGTH1 orthodb8_level=32281 organism_name=`Drosophila virilis` uniprot_de=`GJ21671`"
		d = orthodbutil.headerDict(header)
		self.assertTrue(d['organism_name'] == 'Drosophila virilis')

if __name__=="__main__":
	unittest.main(verbosity=2)
