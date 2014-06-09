#! python

import unittest
import stats
import scipy as sp
import scipy.stats as st

class test001(unittest.TestCase):
	def test(self):
		"""Test Benjamini-Hochberg P-value adjustment"""
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
	def test(self):
		"""Test Benjamini-Hochberg P-value adjustment, ordering"""
		pvals = sp.array(range(100))/1.0e5
		adjps = stats.adjustPValues(pvals, method='FDR')
		for i in range(1,10):
			self.assertTrue(adjps[i]>adjps[i-1])

class test003(unittest.TestCase):
	def test(self):
		"""Test Benjamini-Hochberg P-value adjustment, ordering"""
		pvals = sp.array(range(100,-1,-1))/1.0e5
		adjps = stats.adjustPValues(pvals, method='FDR')
		for i in range(1,10):
			self.assertTrue(adjps[i]<adjps[i-1])

class test004(unittest.TestCase):
	def test(self):
		"""Test Benjamini-Hochberg P-value adjustment, ordering"""
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
	def test(self):
		"""Restore list order"""
		a = [1,3,5,2,4,0,7]
		indexed = sorted(zip(a, range(len(a))), reverse=True)
		(b,inds) = zip(*indexed)
		a2 = [-1]*len(b)
		for i in range(len(a2)):
			a2[inds[i]] = b[i]
		for i in range(len(a)):
			self.assertTrue(a[i] == a2[i])

class test006(unittest.TestCase):
	def test(self):
		"""Histogram"""
		a = [1,3,5,2,4,0,7,7.01]
		hist = stats.Histogram()
		hist.init(min(a), max(a), 10)
		hist.add(a)
		self.assertTrue(len(hist.extras) == 0)
		hist.add(12.0)
		self.assertTrue(len(hist.extras) == 1)


	def testuneq(self):
		"""Histogram distance, equal"""
		a = [1,3,5,2,4,0,7,7.01]
		b = [1,3,5,2,4,0,7,7.01]
		d = stats.chiSquaredHistogramDistance(a,b)
		self.assertAlmostEqual(d, 0.0)

	def testeq(self):
		"""Histogram distance, unequal"""
		a = [1,3,5,2,4,0,7,7.01]
		b = [1,3,5,2,4,0,7,8]
		d = stats.chiSquaredHistogramDistance(a,b)
		self.assertTrue(d>0.0)

	def testknown(self):
		"""Histogram distance, known"""
		a = [1,3]
		b = [3,1]
		d = stats.chiSquaredHistogramDistance(a,b)
		self.assertAlmostEqual(d,1.0)

	def test_bin(self):
		"""Histogram bin value retrieval"""
		a = [1,3,5,2,4,0,7,7.01]
		hist = stats.Histogram(a, n_bins=8)
		self.assertTrue(hist[7].count == 2)
		self.assertTrue(hist[-1].count == 0)

	def test_bin(self):
		"""Histogram total retrieval"""
		a = [1,3,5,2,4,0,7,7.01]
		hist = stats.Histogram(a, n_bins=8)
		#print hist.total
		self.assertTrue(hist.total==len(a))
		hist.add([3,-1])
		self.assertTrue(hist.total==len(a)+2)



if __name__=='__main__':
	unittest.main(verbosity=2)
