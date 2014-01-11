#! python

import unittest
import stats
import scipy as sp
import scipy.stats as st

class test001(unittest.TestCase):
	"""Test Benjamini-Hochberg P-value adjustment"""
	def test(self):
		adjps = stats.adjustPValues([0.1,0.2,0.3,0.4], method='FDR')
		self.assertTrue(round(adjps[0],1)==0.4)
		self.assertTrue(round(adjps[1],1)==0.4)
		self.assertTrue(round(adjps[2],1)==0.4)
		self.assertTrue(round(adjps[3],1)==0.4)

		adjps = stats.adjustPValues([0.1,0.2,0.003], method='fdr')
		self.assertTrue(round(adjps[0],2)==0.15)
		self.assertTrue(round(adjps[1],1)==0.2)
		self.assertTrue(round(adjps[2],3)==0.009)

class test002(unittest.TestCase):
	"""Test Benjamini-Hochberg P-value adjustment, ordering"""
	def test(self):
		pvals = sp.array(range(100))/1.0e5
		adjps = stats.adjustPValues(pvals, method='FDR')
		for i in range(1,10):
			self.assertTrue(adjps[i]>adjps[i-1])

class test003(unittest.TestCase):
	"""Test Benjamini-Hochberg P-value adjustment, ordering"""
	def test(self):
		pvals = sp.array(range(100,-1,-1))/1.0e5
		adjps = stats.adjustPValues(pvals, method='FDR')
		for i in range(1,10):
			self.assertTrue(adjps[i]<adjps[i-1])

class test004(unittest.TestCase):
	"""Test Benjamini-Hochberg P-value adjustment, ordering"""
	def test(self):
		sp.random.seed(111)
		n = 3
		#pvals = sp.array(range(n-1,-1,-1))/1.0e5
		#ranks = st.rankdata(sp.random.random(len(pvals)))
		#pvals = [pvals[x-1] for x in ranks]
		pvals = [0.0, 2e-5, 1e-5]
		adjps = stats.adjustPValues(pvals, method='FDR')
		#print ""
		#for i in range(n):
		#	print "{:1.6E}\t{:1.6E}".format(pvals[i],adjps[i])
		for i in range(n):
			self.assertTrue(adjps[i] >= pvals[i])

class test005(unittest.TestCase):
	"""Restore list order"""
	def test(self):
		a = [1,3,5,2,4,0,7]
		indexed = sorted(zip(a, range(len(a))), reverse=True)
		(b,inds) = zip(*indexed)
		a2 = [-1]*len(b)
		for i in range(len(a2)):
			a2[inds[i]] = b[i]
		for i in range(len(a)):
			self.assertTrue(a[i] == a2[i])

class test006(unittest.TestCase):
	"""Histogram"""
	def test(self):
		a = [1,3,5,2,4,0,7,7.01]
		hist = stats.Histogram()
		hist.init(min(a), max(a), 10)
		hist.add(a)
		self.assertTrue(len(hist.extras) == 0)


if __name__=='__main__':
	unittest.main(verbosity=2)
