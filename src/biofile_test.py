import time, os, random, string, sys, math, traceback, unittest
import biofile

def writeGFF(stream):
	stream.write()

class testExons(unittest.TestCase):
	def test_reading(self):
		mfr = biofile.MultipleFASTAReader(file('./test-biofile/test-multiple-fasta-001.fa', 'r'), biofile.UCSCExonHeader)
		n = 15
		while not mfr.atEnd():
			for ex_list in mfr.exons():
				#print ex_list[0][0].id, len(ex_list), n
				self.assertTrue(len(ex_list) == n)
			#print ''
	
class testCDS(unittest.TestCase):
	def test_reading(self):
		mfr = biofile.MultipleFASTAReader(file('./test-biofile/test-multiple-fasta-001.fa', 'r'), biofile.UCSCExonHeader)
		for cds_alignment in mfr.CDSs():
			L = None
			for entry in cds_alignment:
				if L is None:
					L = len(entry.sequence)
				#print entry.header, entry.sequence
				# Check that every entry in the alignment has the same length
				self.assertTrue(len(entry.sequence)==L)
	
	def test_length(self):
		mfr = biofile.MultipleFASTAReader(file('./test-biofile/test-multiple-fasta-002.fa', 'r'), biofile.UCSCExonHeader)
		for cds_alignment in mfr.CDSs():
			#print len(cds_alignment[0].sequence)
			for entry in cds_alignment:
				# Check that the second gene has the right length
				if entry.header.id == "CG17540-RC" and entry.header.species == 'dm3':
					#print len(entry.sequence)
					self.assertTrue(len(entry.sequence)==1092)

class test002(unittest.TestCase):
	def test_run(self):
		"""secondOrFirstField"""
		x = 'FIRST SECOND'
		self.assertTrue(biofile.firstField(x)=='FIRST')
		self.assertTrue(biofile.secondField(x)=='SECOND')
		y = 'FIRST'
		self.assertTrue(biofile.secondOrFirstField(y)=='FIRST')

class test003(unittest.TestCase):
	def test_run(self):
		"""GFF basic"""
		pass


if __name__=="__main__":
	unittest.main(verbosity=2)
