#! python

import unittest
import stats

class test001(unittest.TestCase):
	def test(self):
		adjps = stats.adjustPValue([0.1,0.2])
		self.assertTrue(adjps[0]==0.2)

if __name__=='__main__':
	unittest.main(verbosity=2)
