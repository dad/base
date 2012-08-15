import time, os, random, string, sys, math, traceback, unittest
import biofile

class test001(unittest.TestCase):
	def test_run(self):
		mfr = biofile.MultipleFASTAReader(file('./test-biofile/test-multiple-fasta-001.fa', 'r'), biofile.UCSCMultipleFASTAHeader)
		n = None
		while not mfr.atEnd():
			for ex_list in mfr.exons():
				if n is None:
					n = len(ex_list)
				print ex_list[0][0].id, len(ex_list), n
				self.assertTrue(len(ex_list) == n)

if __name__=="__main__":
	unittest.main(verbosity=2)
