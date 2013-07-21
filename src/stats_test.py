#! python

import unittest
import stats

class test001(unittest.TestCase):
	def test(self):
		adjps = stats.adjustPValue([0.1,0.2,0.3,0.4])
		self.assertTrue(round(adjps[0],1)==0.4)
		self.assertTrue(round(adjps[1],1)==0.4)
		self.assertTrue(round(adjps[2],1)==0.4)
		self.assertTrue(round(adjps[3],1)==0.4)

		adjps = stats.adjustPValue([0.1,0.2,0.003])
		self.assertTrue(round(adjps[0],2)==0.15)
		self.assertTrue(round(adjps[1],1)==0.2)
		self.assertTrue(round(adjps[2],3)==0.009)

if __name__=='__main__':
	unittest.main(verbosity=2)
