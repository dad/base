import time, os, random, string, sys, math, traceback, unittest
import biofile

class test001(unittest.TestCase):
	def test_run(self):
		mfr = biofile.MultipleFASTAReader(file('./test-biofile/test-multiple-fasta-001.fa', 'r'), biofile.UCSCMultipleFASTAHeader)
		for ex_list in mfr.exons():
			print len(ex_list)

if __name__=="__main__":
	unittest.main(verbosity=2)
