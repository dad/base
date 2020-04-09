#!/usr/bin/python
# Begin stats.py
"""Module for statistics.

Originally written by Jesse Bloom, 2004.
Expanded and maintained by D. Allan Drummond, 2004-2018."""
#
import re, math, os, string, random
import listrank, na
import scipy as sp
import scipy.special

#---------------------------------------------------------------------------------
class StatsError(Exception):
	"""Statistics error."""

def chiSquaredHistogramDistance(h1, h2):
	# DAD: implement
	assert(len(h1) == len(h2))
	chisq_terms = [(q - p)*(q - p)/(q + p) for (q,p) in zip(h1,h2) if q+p>0.0]
	return 0.5 * sum(chisq_terms)

# From http://stackoverflow.com/questions/3679694/a-weighted-version-of-random-choice
def weighted_choice(choices):
	total = sum(w for c, w in choices)
	r = random.uniform(0, total)
	upto = 0
	for c, w in choices:
		if upto + w >= r:
			return c
		upto += w
	assert False, "Shouldn't get here"
  
def weighted_choice_index(weights):
	total = sum(weights)
	r = random.uniform(0, total)
	upto = 0
	for (i,w) in enumerate(weights):
		if upto + w >= r:
			return i
		upto += w
	assert False, "Shouldn't get here"

def choose_index_pair_weighted(n, indices, weights):
	# Choose n pairs with replacement
	found = False
	sum_w = sum(weights)
	pweights = [w/sum_w for w in weights]
	assert len(indices) == len(weights)
	res = []
	while not found:
		pairs = sp.random.choice(indices, size=(int(n*1.1),2), p=pweights, replace=True)
		res += [(i,j) for (i,j) in pairs if i != j]
		found = (len(res) >= n)
	return res[:n]

def choose_index_pair(n, indices):
	# Choose n distinct pairs of indices (no pairs with i=j) with replacement
	found = False
	res = []
	while not found:
		pairs = sp.random.choice(indices, size=(int(n*1.1),2), replace=True)
		res += [(i,j) for (i,j) in pairs if i != j]
		found = (len(res) >= n)
	return res[:n]

class HistogramBin(object):
	def __init__(self, mid, width, count, total=None, cumcount=None):
		self._mid = mid
		self._width = width
		self._count = count
		self._total = total
		self._cumcount = cumcount
	
	@property
	def lower(self):
		return self._mid - self._width/2.0
	
	@property
	def upper(self):
		return self._mid + self._width/2.0
	
	@property
	def count(self):
		return self._count
	
	@property
	def density(self):
		res = None
		if not self._total is None:
			res = self._count/float(self._total)
		return res
	
	@property
	def mid(self):
		return self._mid

	@property
	def width(self):
		return self._width

	@property
	def cum(self):
		return self._cumcount

	@property
	def cumdensity(self):
		res = None
		if not self._total is None and not self._cumcount is None:
			res = self._cumcount/float(self._total)
		return res
	

class HistogramFactory(object):
	def integerHistogram(self, min, max):
		"""Return histogram with bins centered at each integer between min and max, inclusive"""
		h = Histogram(vals=None, n_bins= int(max)-int(min)+1, min_val=min-0.5, max_val=max+0.5)
		return h

class Histogram:
	def __init__(self, vals=None, n_bins=1, min_val=None, max_val=None):
		if not vals is None:
			if min_val is None:
				min_val = min(vals)
			if max_val is None:
				max_val = max(vals)
			self.init(min_val, max_val, n_bins)
			self.add(vals)
		else:
			self._bins = [0]*n_bins
			self._min_val = 0.0
			if not min_val is None:
				self._min_val = min_val
			self._max_val = 1.0
			if not max_val is None:
				self._max_val = max_val
			self._bin_width = (self._max_val-self._min_val)/float(n_bins)
			self._extras = []
			self._total_count = 0

	def init(self, min_val, max_val, n_bins):
		assert n_bins > 0
		self._bins = [0]*n_bins
		self._min_val = min_val
		self._max_val = max_val
		self._bin_width = (max_val-min_val)/float(n_bins)
		self._total_count = 0
		self._extras = []

	def getBinIndex(self, x):
		"""Retrieve the numeric index of the bin corresponding to value x"""
		b = -1
		if x == self._max_val: # final bin is [low, high], where others are [low,high)
			b = len(self._bins)-1
		else:
			b = math.floor((x-self._min_val)/self._bin_width)
		return int(b)

	def getBinCount(self, binindex):
		res = 0
		if binindex > 0 and binindex < len(self._bins):
			res = self._bins[binindex]
		return res

	def validBin(self, b):
		return b >= 0 and b < len(self._bins)

	def _add(self, x):
		b = self.getBinIndex(x)
		if self.validBin(b):
			self._bins[b] += 1
		else:
			#print x, b, len(self.bins), ((self.max_val-self.min_val)/self.bin_width
			self._extras.append(x)
		self._total_count += 1

	def add(self, x):
		if isinstance(x,list):
			for y in x:
				self._add(y)
		else:
			self._add(x)

	def values(self):
		pass

	def __str__(self):
		s = 'mid\tcount\tdensity\n'
		for b in self.bins:
			line = '{b.mid:f}\t{b.count:d}\t{b.density:f}\n'.format(b=b)
			s += line
		return(s)
	
	def write(self, stream, header=None):
		stream.write(str(self))

	@property
	def total(self):
		return self._total_count
	
	@property
	def extras(self):
		return self._extras[:]

	def _bin_mid(self,binindex):
		return self._min_val + self._bin_width*(binindex+0.5)

	@property
	def size(self):
		return len(self._bins)
	
	@property
	def bins(self):
		cum = 0
		for (bi, b) in enumerate(self._bins):
			#bin_mid = self._min_val + width*(bi+0.5)
			cum += b
			bini = HistogramBin(self._bin_mid(bi), self._bin_width, b, self._total_count, cum)
			yield bini

	def __getitem__(self, x):
		# Get 
		bi = self.getBinIndex(x)
		return HistogramBin(self._bin_mid(bi), self._bin_width, self._bins[bi], self._total_count, None)

class Summary:
	def __init__(self):
		self.mean = None
		self.median = None
		self.sd = None
		self.se = None
		self.variance = None
		self.n = None
		self.na = None
		self.sum = None

	def __str__(self):
		if self.n == 0:
			return "no data"
		else:
			return "mean = {me.mean:1.2E}, sd = {me.sd:1.2E}, N = {me.n:d}".format(me=self)


class Accumulator(object):
	def __init__(self, x=None, store=True):
		self._sum = 0.0
		self._sum_sq = 0.0
		self._n = 0
		self._store = store
		self._na = 0
		self._min = 1e100
		self._max = -self._min
		if self._store:
			self._data = []
		if not x is None:
			self.addAll(x)

	def add(self, x):
		if not na.isNA(x):
			self._sum += x
			self._sum_sq += x*x
			self._n += 1
			if self._min > x:
				self._min = x
			if self._max < x:
				self._max = x
			if self._store:
				self._data.append(x)
		else:
			self._na += 1

	def addAll(self, x_list):
		for x in x_list:
			self.add(x)

	@property
	def mean(self):
		mean = 0.0
		if self._n > 0:
			mean = self._sum/self._n
		return mean

	@property
	def min(self):
		return self._min

	@property
	def max(self):
		return self._max

	@property
	def median(self):
		res = None
		if self._store and self._n>0:
			res = sp.median(self._data)
		return res

	@property
	def variance(self):
		if self._n == 1:
			return 0.0
		if self._n < 1:
			return None
		mu = self.mean
		# Sample variance
		return (1.0/(self._n-1.0))*(self._sum_sq - self._n*mu*mu)

	@property
	def sum(self):
		return self._sum

	@property
	def n(self):
		return self._n

	@property
	def na(self):
		return self._na
	
	@property
	def sd(self):
		res = None
		if self.n == 1:
			res = 0.0
		if self.n > 1:
			res = sp.sqrt(self.variance)
		return res

	@property
	def cv(self):
		# Coefficient of variation, sd/mean
		res = 0.0
		if self.n >= 1:
			res = self.sd/self.mean
		return res

	@property
	def se(self):
		if self.n == 1:
			return 0.0
		if self.n < 1:
			return None
		return self.sd/sp.sqrt(self.n)

	def getSEConfidenceInterval(self, alpha):
		assert False, "Not implemented yet"

	@property
	def data(self):
		res = None
		if self._store:
			res = self._data
		return res

	@property
	def summary(self):
		s = Summary()
		s.mean = self.mean
		s.n = self.n
		s.na = self.na
		s.variance = self.variance
		if not s.variance is None:
			s.sd = self.sd
			s.se = self.se
		s.sum = self.sum
		s.median = self.median
		return s

class LogAccumulator(Accumulator):
	def __init__(self, x=None, store=True):
		self._nolog_sum = 0.0
		super(LogAccumulator,self).__init__(x, store)
	
	def add(self, x):
		if not na.isNA(x) and x>0.0:
			super(LogAccumulator,self).add(math.log(x))
			self._nolog_sum += x
		else:
			self._na += 1

	@property
	def mean(self):
		log_mean = super(LogAccumulator,self).mean
		return math.exp(log_mean)

	@property
	def median(self):
		med = None
		log_median = super(LogAccumulator,self).median
		if not log_median is None:
			med = math.exp(log_median)
		return med

	@property
	def sum(self):
		return self._nolog_sum

	@property
	def variance(self):
		if self._n == 1:
			return 0.0
		if self._n < 1:
			return None
		mu = math.log(self.mean)
		# Sample variance
		return (1.0/(self._n-1.0))*(self._sum_sq - self._n*mu*mu)

# Sample with replacement
# http://code.activestate.com/recipes/273085-sample-with-replacement/
def sample_wr(population, k):
	"Chooses k random elements (with replacement) from a population"
	n = len(population)
	_random, _int = random.random, int  # speed hack
	return [population[x] for x in [_int(_random() * n) for i in range(k)]]

def adjustPValues(p_values, method="fdr"):
	"""Adjust P values for multiple testing.
	
	Methods (case-insensitive):
		'FDR' or 'BH': Benjamini-Hochberg false discovery rate correction
		'Bonferroni': Bonferroni correction
	"""
	adjusted_p_values = p_values[:]
	n = len(p_values)
	if method.lower() == "bh" or method.lower() == 'fdr':
		ni = range(n,0,-1) # from n to 1
		# Sort the P values and keep track of the indices
		indexed_pv = sorted(zip(p_values, range(n)), reverse=True)
		(pvals,inds) = zip(*indexed_pv)
		# adjust
		newp = [(float(n)/ni[xi])*pvals[xi] for xi in range(n)]
		cum_min_p = [min(newp[0:xi]) for xi in range(1,n+1)]
		adjp_sorted = [min(p,1.0) for p in cum_min_p]
		# re-sort
		adjusted_p_values = [-1]*n
		for xi in range(n):
			adjusted_p_values[inds[xi]] = adjp_sorted[xi]
	elif method.lower() == 'bonferroni':
		adjusted_p_values = [min(n*p,1.0) for p in p_values]
	return adjusted_p_values

#---------
def statsSummary(numlist):
	return summary(numlist)

def summary(numlist):
	acc = Accumulator(numlist, store=True)
	return acc.summary
#---------------------------------------------------------------------------------
def median(numlist):
	return sp.median(numlist)

def Median(numlist):
	return median(numlist)
#----------------------------------------------------------------------------------
def mean(numlist):
	return sp.mean(numlist)

def Mean(numlist):
	return mean(numlist)
#----------------------------------------------------------------------------------
def geometricMean(numlist):
	"""Returns the geometric mean of a list of numbers.

	If any entries of the list are 'None' or '-', they are removed
	first."""
	if len(numlist)==0:
		raise StatsError("Empty list.")
	mean = 0.0
	log_sum = 0.0
	n = 0
	for x in numlist:
		if x in [None, '-']:
			continue
		assert isinstance(x, (int, float))
		if x > 0:
			log_sum += math.log(x)
			n += 1
	if n == 0:
		# Can happen if only entry is zero.
		return 0.0
	return math.exp(log_sum / float(n))
#----------------------------------------------------------------------------------
def weightedMean(numlist, weights):
	"""Returns the weighted mean of a list of numbers.

	If any entries of the list are 'None' or '-', they are removed
	first."""
	wxsum = 0.0
	wsum = 0.0

	assert len(numlist) == len(weights)

	for (x,w) in zip(numlist, weights):
		wxsum += x*w
		wsum += w
	if wsum == 0.0:
		return 0.0
	return wxsum/wsum

#--------------------------------------------------------------------------------
def Variance(numlist):
	"""Returns the sample variance of a list of numbers.
	If length is <2, then variance = 0."""
	sum1 = sum2 = 0.0
	n = 0.0
	for x in numlist:
		assert isinstance(x, int) or isinstance(x, float)
		sum1 += x
		sum2 += x * x
		n += 1.0
	if n < 2.0:
		return 0.0
	var = (1.0/n)*(sum2 - (1/n)*sum1*sum1)
	if var < 0.0: # Due to numerical problems only!
		var = 0.0
	return var
#-------------------------------------------------------------------------------
def sampleVariance(numlist):
	"""Returns the sample variance of a list of numbers.
	If length is <2, then variance = 0."""
	sum1 = sum2 = 0.0
	n = 0.0
	for x in numlist:
		assert isinstance(x, int) or isinstance(x, float)
		sum1 += x
		sum2 += x * x
		n += 1.0
	if n < 2.0:
		return 0.0
	var = (1.0/(n+1.0))*(sum2 - (1/n)*sum1*sum1)
	if var < 0.0: # Due to numerical problems only!
		var = 0.0
	return var
#-------------------------------------------------------------------------------
def StandardDeviation(numlist):
	"""Returns the sample standard deviation of a list of numbers.

	If any entries of the list are 'None' or '-', they are removed first."""
	v = Variance(numlist)
	#print v
	return math.sqrt(v)
#-------------------------------------------------------------------------------
def sampleStandardDeviation(numlist):
	"""Returns the sample standard deviation of a list of numbers.

	If any entries of the list are 'None' or '-', they are removed first."""
	v = sampleVariance(numlist)
	#print v
	return math.sqrt(v)
#-------------------------------------------------------------------------------
def Kendalls_Tau(xlist, ylist):
	"""Calculates Kendall's tau non-parametric correlation between two variables.

	The input data is given in the two lists 'xdata' and 'ydata' which should be
	of the same length.  If entry i of either list is 'None', this entry is
	disregarded in both lists.
	Returns Kendall's partial tau, the one-tailed P-value, and the number of
	data points as a tuple: (tau, P, N).
	Includes a correction for ties.
	Based on Gibbons, JD, "Nonparametric measures of association",
	Sage University Papers, pg 15 (1983)."""
	if len(xlist) != len(ylist):
		raise StatsError("Data sets have different lengths.")
	xdata = []
	ydata = []
	for i in range(len(xlist)):
		if xlist[i] != None and ylist[i] != None:
			xdata.append(xlist[i])
			ydata.append(ylist[i])
	assert len(xdata) == len(ydata)
	assert len(xdata) <= len(xlist) - xlist.count(None)
	assert len(ydata) <= len(ylist) - ylist.count(None)
	assert len(ydata) >= len(ylist) - xlist.count(None) - ylist.count(None)
	if len(xdata) == 0:
		raise StatsError("No valid data entries.")
	n = len(xdata)
	# compute the number of concordant and discordant pairs
	conc = disc = 0.0 # concordant and discordant pairs
	for i in range(n): # loop over all pairs
		xi = xdata[i]
		yi = ydata[i]
		for j in range(i + 1, n):
			xd = xi - xdata[j]
			yd = yi - ydata[j]
			prod = xd * yd
			if prod == 0.0: # this is a tie
				continue
			elif prod > 0.0:
				conc += 1
			else:
				disc += 1
	# compute the tie correction: sum(t * t - t)
	xcopy = []
	ycopy = []
	for i in range(n):
		xcopy.append(xdata[i])
		ycopy.append(ydata[i])
	xties = yties = 0.0
	while xcopy:
		xi = xcopy[0]
		t = xcopy.count(xi)
		xties = xties + t * t - t
		while xcopy.count(xi) > 0:
			xcopy.remove(xi)
	while ycopy:
		yi = ycopy[0]
		t = ycopy.count(yi)
		yties = yties + t * t - t
		while ycopy.count(yi) > 0:
			ycopy.remove(yi)
	# Compute tau
	n = float(n)
	denom = math.sqrt((n * n - n - xties) * (n * n - n - yties))
	try:
		tau = 2.0 * (conc - disc) / denom
	except ZeroDivisionError:
		raise StatsError("Too few entries: {:d}.".format(n))
	# Compute P-value
	z = 3.0 * tau * math.sqrt(n * (n - 1.0)) / math.sqrt(2.0 * (2.0 * n + 5.0))
	prob = Prob_Z(z)
	return (tau, prob, int(n))

#-------------------------------------------------------------------------------
def Prob_Z(z, twosided=False):
	return probZ(z, twosided)

def probZ(z, twosided=False):
	p = scipy.special.erfc(z / math.sqrt(2.0))/2.0
	if twosided:
		if z < 0.0:
			p = 2*(1-p)
		else:
			p = 2*p
	return p

#-------------------------------------------------------------------------------
def Kendalls_Tau2(xlist, ylist):
	"""Calculates Kendall's tau non-parametric correlation between two variables.

	The input data is given in the two lists 'xdata' and 'ydata' which should be
	of the same length.  If entry i of either list is 'None', this entry is
	disregarded in both lists.
	Returns Kendall's partial tau, the one-tailed P-value, and the number of
	data points as a tuple: (tau, P, N).
	Includes a correction for ties.
	Based on Numerical Recipes in C."""
	if len(xlist) != len(ylist):
		raise StatsError("Data sets have different lengths.")
	xdata = xlist
	ydata = ylist
	#for i in range(len(xlist)):
	#	if xlist[i] != None and ylist[i] != None:
	#		xdata.append(xlist[i])
	#		ydata.append(ylist[i])
	assert len(xdata) == len(ydata)
	#assert len(xdata) <= len(xlist) - xlist.count(None)
	#assert len(ydata) <= len(ylist) - ylist.count(None)
	#assert len(ydata) >= len(ylist) - xlist.count(None) - ylist.count(None)
	if len(xdata) == 0:
		raise StatsError("No valid data entries.")
	n = len(xdata)
	# compute the number of concordant and discordant pairs
	conc = disc = 0.0 # concordant and discordant pairs
	nx = ny = 0.0
	updown = 0
	for i in range(n): # loop over all pairs
		xi = xdata[i]
		yi = ydata[i]
		if xi and yi:
			for j in range(i + 1, n):
				if xdata[j] and ydata[j]:
					xd = xi - xdata[j]
					yd = yi - ydata[j]
					prod = xd * yd
					if prod != 0:
						nx += 1
						ny += 1
						if prod > 0:
							updown += 1
						else:
							updown -= 1
					else:
						if xd != 0:
							nx += 1
						if yd != 0:
							ny += 1
	# Compute tau
	n = float(n)
	denom = math.sqrt(nx*ny)
	try:
		tau = float(updown) / denom
	except ZeroDivisionError:
		raise StatsError("Too few entries: {:d}".format(n))
	# Compute P-value
	z = 3.0 * tau * math.sqrt(n * (n - 1.0)) / math.sqrt(2.0 * (2.0 * n + 5.0))
	prob = Prob_Z(z)
	return (tau, prob, int(n))
#----------------------------------------------------------------------------------
def Kendalls_Partial_Tau(xdata, ydata, zdata):
	"""Computes Kendall's partial tau of two variables controlling for a third.

	The correlation is between 'xdata' and 'ydata' controlling for 'zdata'.
	The data is given in lists that must be of the same length.
	Returns partial tau as a scalar number.
	Based on Gibbons JD, "Nonparametric measures of associations",
	Sage University Papers, pg 49 (1983)."""
	if not len(xdata) == len(ydata) == len(zdata):
		raise StatsError("Data sets have different lengths.")
	txy = Kendalls_Tau(xdata, ydata)[0]
	tyz = Kendalls_Tau(ydata, zdata)[0]
	txz = Kendalls_Tau(xdata, zdata)[0]
	partial_tau = (txy - txz * tyz) / math.sqrt((1 - txz * txz) * (1 - tyz * tyz))
	return partial_tau
#------------------------------------------------------------------------------
def pearsonCorrelation(x, y):
	"""Computes the Pearson linear correlation between two data sets.

	Call is '(r, p, n) = PearsonCorrelation(xdata, ydata)'
	The input data is given in the two lists 'xdata' and 'ydata' which
	should be
	of the same length.  If entry i of either list is 'None', this
	entry is
	disregarded in both lists.
	Returns Pearson's correlation coefficient, the two-tailed P-value,
	and the number of data points as a tuple '(r, p, n)'."""
	sum_sq_x = 0
	sum_sq_y = 0
	sum_coproduct = 0
	mean_x = x[0]
	mean_y = y[0]
	if len(x) != len(y):
		raise StatsError("Data sets are of different lengths.")
	n = len(x)
	for i in range(1,n):
		sweep = i / (i+1.0)
		delta_x = x[i] - mean_x
		delta_y = y[i] - mean_y
		sum_sq_x += delta_x * delta_x * sweep
		sum_sq_y += delta_y * delta_y * sweep
		sum_coproduct += delta_x * delta_y * sweep
		mean_x += delta_x / (i+1.0)
		mean_y += delta_y / (i+1.0)
	pop_sd_x = math.sqrt( sum_sq_x / n )
	pop_sd_y = math.sqrt( sum_sq_y / n )
	cov_x_y = sum_coproduct / n
	r = cov_x_y / (pop_sd_x * pop_sd_y)
	z = math.fabs(r) * math.sqrt(n) / math.sqrt(2.0)
	p = Prob_Z(z)
	if not (0.0 <= p <= 1.0):
		raise StatsError("Invalid P-value of %r." % r)
	return (r, p, n)

def PearsonCorrelation(xdata, ydata):
	return pearsonCorrelation(xdata, ydata)

def test_pearsonCorrelation():
	eps = 1e-6
	nv = 100
	for i in range(100):
		nv = random.randint(10,1000)
		add = random.random()
		x = [random.random() for xi in range(nv)]
		y = [random.random()+add*x[xi] for xi in range(nv)]
		(r,p,n) = PearsonCorrelation(x,y)
		(r2,p2,n2) = pearsonCorrelation2(x,y)
		print(nv, r, r2)
	return True
#------------------------------------------------------------------------------
def PartialPearsonCorrelation(xdata, ydata, zdata):
	"""Computes the Pearson linear correlation between two data sets controlling for a third set.

	Call is '(r, p, n) = PartialPearsonCorrelation(xdata, ydata, zdata)'
	The input data is given in the two lists 'xdata' and 'ydata' which
	should be of the same length.
	Returns Pearson's partial correlation coefficient, the two-tailed P-value,
	and the number of data points as a tuple '(r, p, n)'."""
	try:
		(rxy, dummy, n) = PearsonCorrelation(xdata, ydata)
		(ryz, dummy, n) = PearsonCorrelation(ydata, zdata)
		(rxz, dummy, n) = PearsonCorrelation(xdata, zdata)
		r = (rxy - ryz*rxz)/math.sqrt((1-ryz**2)*(1-rxz**2))
	except ZeroDivisionError:
		raise StatsError("Standard deviation is zero.")
	if not (-1.0000000001 <= r <= 1.000000001):
		raise StatsError("Invalid correlation coefficient of %r." % r)
	t = r*math.sqrt((n-3)/(1-r*r))
	z = t
	p = Prob_Z(z)
	if not (0.0 <= p <= 1.0):
		raise StatsError("Invalid P-value of %r." % r)
	return (r, p, n)
#------------------------------------------------------------------------------
def SpearmanRankCorrelation(xdata, ydata, ties="average"):
	"""Computes the nonparametric Spearman rank correlation between two data sets."""
	xranks = listrank.rank(xdata, ties=ties)
	yranks = listrank.rank(ydata, ties=ties)
	return PearsonCorrelation(xranks, yranks)

def rank(xdata, ties="average"):
	return listrank.rank(xdata, ties)

#------------------------------------------------------------------------------

def factorial(n):
	if n <= 1:
		return 1
	else:
		return n*factorial(n-1)

__factorial_cache = []
__max_factorial_cache = 100
for n in range(__max_factorial_cache):
	__factorial_cache.append(factorial(n))

def logFactorial(n):
	if n < __max_factorial_cache:
		return math.log(__factorial_cache[n])

	# Gosper's approximation to the factorial; much more accurate than Stirling's.
	res_approx = n*math.log(n)-n+math.log(math.sqrt((2*n+1/3.0)*math.pi))
	#res_exact = sum([math.log(i+1) for i in range(0,n)])
	return res_approx

def logChoose(n,k):
	#num = sum([math.log(i+1) for i in range(n-k,n)])
	#denom = log_factorial(k)
	#return num - denom
	return logFactorial(n) - (logFactorial(k)+logFactorial(n-k))

def logBinom(n,k,p):
	if p <= 0.0 or p == 1.0:
		raise StatsError("Log factorial for p <= 0, %f" % p)
	if p > 1.0:
		raise StatsError("Log factorial for p > 1, %f" % p)

	lc = logChoose(n,k)
	log_p = math.log(p)
	return lc + k*log_p + (n-k)*math.log(1.0-p)


def binomialTest(k, n, p = 0.5, exact = False):
	"""Computes the exact probability for a binomial process to produce at least
	k successes in n tries given binomial proportion p.

	Call is p = binomialTest(k, n, p=0.5, exact = False)
	Conditions are 0 < k <= n and n > 0.  If n*p and n*(1-p) > 30, a normal
	approximation is used.  Use 'exact = True' to force an exact calculation."""
	assert(k <= n)
	assert(k >= 0 and n > 0)
	n = int(n)
	k = int(k)
	p_value = 1.0

	# Trivial cases where p = 0 or p = 1
	if p == 0.0:  # Must then have k = 0
		if k > 0:
			return 0.0
		else:
			return 1.0
	if p == 1.0:  # Must then have k = n
		if k <= n:
			return 1.0

	if k == 0:
		# Probability of at least zero successes is 1
		p_value = 1.0
	elif k == n:
		# Probability of all successes
		p_value = p**n
	else:
		if not exact and n*p > 30 and n*(1-p) > 30:
			# Use normal approximation
			mu = n*p
			sd = math.sqrt(n*p*(1-p))
			z = (k-mu)/sd
			if z < 0.0:
				p_value = 1-Prob_Z(z)
			else:
				p_value = Prob_Z(z)
		else:
			p_value = p**n # The last term in the sum
			for j in range(k,n):
				# Compute logarithm of (n choose j) p^j (1-p)^ (n-j), the
				# binomial probability.  Use logarithm to avoid overflow
				# problems with potentially enormous factorials.
				log_p = logChoose(n,j) + j*math.log(p) + (n-j)*math.log(1-p)
				p_value += math.exp(log_p)
			if p_value > 1.0:
				p_value = 1.0
	return p_value
#------------------------------------------------------------------------------
def WilcoxonTest(list1, list2):
	"""Computes the Wilcoxon signed-rank probability that list1's median differs
	from list2's by the observed amount."""
	p_value = 1.0
	assert(len(list1)==len(list2))
	diffs = [b-a for (a,b) in zip(list1, list2)]
	ranks = zip([math.fabs(d) for d in diffs], range(1,len(diffs)+1))
	ranks.sort()
	W = 0.0
	for i in range(len(diffs)):
		if diffs[i] > 0:
			W += ranks[i][1]
	#print "#", W
	n = len(list1)
	mean = n*(n+1)/4.0
	stdev = math.sqrt(n*(n+1)*(2*n+1)/6.0)
	if stdev > 0:
		Z = (W - mean)/stdev
	else:
		if Median(list1) < Median(list2):
			return 0.0
		else:
			return 1.0
	p_value = Prob_Z(Z)
	return p_value

#------------------------------------------------------------------------------
def Means_Differ(pop1, pop2):
	(n1, m1, sd1, se1) = StatsSummary(pop1)
	(n2, m2, sd2, se2) = StatsSummary(pop2)
	denom = math.sqrt(se1**2 + se2**2) #sd1**2/n1 + sd2**2/n2)
	p = 1.0
	if denom > 0.0:
		z = (m1-m2)/denom
		p = Prob_Z(z)
	else: # No variance -- so just check if means differ
		if m1 != m2:
			p = 0.0
		else:
			p = 1.0
	if not (-1e-6 <= p <= 1.0+1e-6):
		raise StatsError("Invalid P-value of %r." % p)
	return p

#------------------------------------------------------------------------------
def Complementary_Error_Function(z):
	"""Calculates the error function of z.

	The complementary error function of z is defined as:
	erfc(z) = 2 / sqrt(pi) * integral(e^(t^2) dt) where the integral
	is from z to infinity.
	Can be used to calculate cumulative normal probabilities: given a
	distribution with mean m and standard deviation s,
	the probability of observing x > m  when x > 0 is:
	P = 0.5 * erfc((x - m) / (s * sqrt(2)))
	Calculated according to Chebyshev fitting given by Numerical Recipes
	in C, page 220-221."""
	x = math.fabs(z)
	t = 1.0 / (1.0 + 0.5 * x)
	ans = t * math.exp(-x * x - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))))
	return ans
#--------------------------------------------------------------------------------
def Poisson(n, k):
	"""Returns the Poisson probability of observing a number.

	'Poisson(n, k)' takes as input an integer n >= 0 and a real number k >= 0.0.
	Returns p, the probability of getting n counts when the average outcome
	is k, according to the Possion distribution.  Returns 'None' if there is
	an error."""
	p = math.exp(-k) * math.pow(k, n) / float(Factorial(n))
	assert 0.0 <= p <= 1.0, "Error, value of p is invalid probability: " + str(p)
	return p
#---------------------------------------------------------------------------
def Factorial(n):
	"""Returns the factorial of an integer."""
	x = 1
	for i in range(1, n + 1):
		x *= i
	return x
#----------------------------------------------------------------------------------
def Choose(n, k):
	num = 1
	for i in range(n-k+1, n+1):
		num *= i
	den = Factorial(k)
	return num/den

def powerSet(s):
	d = dict(zip((1<<i for i in range(len(s))), (set([e]) for e in s) ))
	subset = set()
	yield subset
	for i in range(1, 1<<len(s)):
		subset = subset ^ d[i & -i]
		yield subset

def generateChoices(seqin,k):
	'''returns a generator which returns combinations of argument sequences without replacement
	for example generateChoices((1,2,3),2) returns a generator; calling the next()
	method on the generator will return [1,2], [1,3], [2,3] and the
	StopIteration exception.  This will not create the whole list of
	combinations in memory at once.'''
	def rloop(seqin, comb, k):
		'''recursive looping function'''
		if k>0:                   # any more sequences to process?
			for i in range(len(seqin)):
				item = seqin[i]
				newcomb=comb+[item]     # add next item to current combination
				# call rloop w/ remaining seqs, newcomb
				for item in rloop(seqin[(i+1):],newcomb,k-1): # remove this and previous items from consideration
					yield item          # seqs and newcomb
		else:                           # processing last sequence
			yield comb                  # comb finished, add to list
	return rloop(seqin,[],k)

def generateCombinations(seqin):
	'''returns a generator which returns combinations of argument sequences
	for example generateCombinations((1,2),(3,4)) returns a generator; calling the next()
	method on the generator will return [1,3], [1,4], [2,3], [2,4] and
	StopIteration exception.  This will not create the whole list of
	combinations in memory at once.'''
	def rloop(seqin,comb):
		'''recursive looping function'''
		if seqin:                   # any more sequences to process?
			for item in seqin[0]:
				newcomb=comb+[item]     # add next item to current combination
				# call rloop w/ remaining seqs, newcomb
				for item in rloop(seqin[1:],newcomb):
					yield item          # seqs and newcomb
		else:                           # processing last sequence
			yield comb                  # comb finished, add to list
	return rloop(seqin,[])


#----------------------------------------------------------------------------------
def __getAVEA(a, b, c, d):
	"""Returns the number of at-risk+outcome cases, expected number, and variance.

	   For i'th level, given risk factor ("exposure") E and outcome ("disease") D, the annotated
	   2x2 contingency table is:

		   |  E  | ~E  | Total
	   -----------------------
		 D | ai  |  bi |  m1i
	   -----------------------
		~D | ci  |  di |  m0i
	   -----------------------
		   | n1i | n0i |  ni

	   The Mantel-Haenszel statistic =
		 chi-squared_MH = (A - E(A))^2 / Var(A)
					  A = sum_i a_i         -- Number of disease cases associated with at-risk factor
				   E(A) = sum_i n1i*m1i/ni  -- Expected disease+at-risk cases if no association
				 Var(A) = sum_i n1i*n0i*m1i*m0i/((ni-1)*ni^2)
	"""
	assert(a>-1)
	assert(b>-1)
	assert(c>-1)
	assert(d>-1)
	m1i = a+b
	m0i = c+d
	n1i = a+c
	n0i = b+d
	ni = a+b+c+d
	v = 0.0
	ea = 0.0
	if ni<=1:
		# Avoid divide-by-zero
		ea = n1i*m1i
	else:
		v = n1i*n0i*m1i*m0i/float((ni-1)*ni*ni)
		ea = n1i*m1i/float(ni)
	return (a,v,ea)

def getOddsRatio(table, pseudocount=0):
	(a,b,c,d) = table
	# Add pseudocount if necessary
	if pseudocount>0:
		(a,b,c,d) = tuple(sp.array(table)+pseudocount)
	res = None
	if b*c>0:
		res = (a*d)/(b*c)
	return res
	
def getZScore(table, pseudocount=0):
	(a,b,c,d) = table
	# Add pseudocount if necessary
	if pseudocount>0:
		(a,b,c,d) = tuple(sp.array(table)+pseudocount)
	res = __getAVEA(a,b,c,d)
	return __ZStat(*res)
	

def getSummaryMHStats(tables):
	sum_a = 0
	sum_v = 0
	sum_ea = 0
	sum_odds_num = 0
	sum_odds_den = 0

	for table in tables:
		(a1,b1,c1,d1) = table
		ni = float(sum(table))
		# Eliminate tables with ni < 1.0, as these can
		# give undefined variances.
		if ni <= 1.0:
			continue
		(a,v,ea) = __getAVEA(a1,b1,c1,d1)
		sum_a += a
		sum_v += v
		sum_ea += ea
		if ni>0.0:
			sum_odds_num += a1*d1/ni
			sum_odds_den += b1*c1/ni
	return (sum_a, sum_v, sum_ea, sum_odds_num, sum_odds_den)

def MantelHaenszelZ(a1,b1,c1,d1):
	(a,v,ea) = __getAVEA(a1,b1,c1,d1)
	z = __ZStat(a,v,ea)
	p = Prob_Z(z)
	if not (0.0 <= p <= 1.0):
		raise StatsError("Invalid P-value of %r." % z)
	return z, p

def __ZStat(a,v,ea):
	z = 0
	if v>0:
		z = (a-ea)/math.sqrt(v)
	return z

def MantelHaenszelSummaryZ(tables):
	"""Returns the score and overall

	   For i'th level, given risk factor ("exposure") E and outcome ("disease") D, the annotated
	   2x2 contingency table is:

		   |  E  | ~E  | Total
	   -----------------------
		 D | ai  |  bi |  m1i
	   -----------------------
		~D | ci  |  di |  m0i
	   -----------------------
		   | n1i | n0i |  ni

	   The Mantel-Haenszel statistic =
		 chi-squared_MH = (A - E(A) -0.5)^2 / Var(A)
					  A = sum_i a_i         -- Number of disease cases associated with at-risk factor
				   E(A) = sum_i n1i*m1i/ni  -- Expected disease+at-risk cases if no association
				 Var(A) = sum_i n1i*n0i*m1i*m0i/((ni-1)*ni^2)
				 0.5 is a continuity correction.

	   Because the contingency table includes information on the direction of association, but
	   the chi-squared statistic destroys this information, the quantity

		   Z = (A-E(A)-0.5)/sqrt(Var(A))

	   is instead returned, along with a probability of this Z-score assuming normality.
	"""
	(sum_a, sum_v, sum_ea, sum_odds_num, sum_odds_den) = getSummaryMHStats(tables)
	if sum_v == 0.0:
		raise StatsError("Variance of summary M-H tables is zero; can't compute Z or P.")
	z = __ZStat(sum_a,sum_v,sum_ea)
	p = Prob_Z(z)
	if not (0.0 <= p <= 1.0):
		raise StatsError("Invalid P-value of %r (Z=%E)." % (p, z))
	return z, p

def MantelHaenszelOddsRatio(tables):
	"""Returns the score and overall

	   For i'th level, given risk factor ("exposure") E and outcome ("disease") D, the annotated
	   2x2 contingency table is:

		   |  E  | ~E  | Total
	   -----------------------
		 D | ai  |  bi |  m1i
	   -----------------------
		~D | ci  |  di |  m0i
	   -----------------------
		   | n1i | n0i |  ni

	   The Mantel-Haenszel odds ratio =
		 chi-squared_MH = X/Y
					  X = sum_i a_i*d_i/n_i
					  Y = sum_i b_i*c_i/n_i
	"""
	(sum_a, sum_v, sum_ea, sum_odds_num, sum_odds_den) = getSummaryMHStats(tables)
	#print (sum_a, sum_v, sum_ea, sum_odds_num, sum_odds_den)
	if sum_odds_den==0:
		raise StatsError("Denominator of odds ratio is zero; can't compute M-H odds ratio.")
	else:
		odds_ratio = sum_odds_num/sum_odds_den
	return odds_ratio


class MHVarResult:
	"""Class for storing results of Mantel-Haenszel variance computation"""

	def __init__(self, odds_ratio, var_odds_ratio, n_tables, n_counts):
		self.odds_ratio = odds_ratio
		self.var_odds_ratio = var_odds_ratio
		self.sd_odds_ratio = math.sqrt(var_odds_ratio)
		self.ln_odds_ratio = math.log(odds_ratio)
		# From Robins et al. 1986, bottom p.312
		self.var_ln_odds_ratio = self.var_odds_ratio/(self.odds_ratio**2.0)
		self.sd_ln_odds_ratio = math.sqrt(self.var_ln_odds_ratio)
		self.n_tables = n_tables
		self.n_counts = n_counts

	# Get confidence interval (CI)
	# Currently only returns 95% CI
	def getConfidenceInterval(self, log=False):
		fn = math.exp
		if log:
			fn = lambda x: x
		sd_ln95 = math.sqrt(self.var_ln_odds_ratio)*1.96
		mh_lower_95 = fn(self.ln_odds_ratio - sd_ln95)
		mh_upper_95 = fn(self.ln_odds_ratio + sd_ln95)
		return (mh_lower_95, mh_upper_95)


def MantelHaenszelOddsRatioVariance(tables):
	"""Returns the Robins et al. variance phi_US(W) for Mantel-Haenszel odds ratio W, per
	   Robins et al. Biometrics June 42:311-323 (1986).

		   |    E    |   ~E    | Total
	   -------------------------------
		 D |   X_k   |   Y_k   |  t_k
	   -------------------------------
		~D | n_k-X_k | m_k-Y_k |  N_k - t_k
	   -------------------------------
		   |   n_k   |   m_k   |  N_k

		Var_US(W) = [sum_k P_k R_k/2R_t^2 + sum_k (P_k S_k + Q_k R_k)/(2 R_t S_t) + sum_k Q_k S_k/2S_t^2] (W)^2

		where

			P_k = (X_k + m_k - Y_k)/N_k
			Q_k = (Y_k + n_k - X_k)/N_k
			R_k = X_k(m_k - Y_k)/N_k
			S_k = Y_k(n_k - X_k)/N_k
			R_t  = sum_k R_k
			S_t  = sum_k S_k
			odds ratio = R_t / S_t

		Raises StatsError if either R_t or S_t = 0.
	"""
	N_t = 0
	R_t = 0
	S_t = 0
	sum_PR = 0.0
	sum_PS = 0.0
	sum_QR = 0.0
	sum_QS = 0.0
	n_tables = 0
	for table in tables:
		(X_k, Y_k, n_k_X_k, m_k_Y_k) = table
		N_k = float(X_k + Y_k + n_k_X_k + m_k_Y_k)
		if N_k <= 0:
			continue
		n_tables += 1
		P_k = (X_k + m_k_Y_k)/N_k
		Q_k = (Y_k + n_k_X_k)/N_k
		R_k = X_k*m_k_Y_k/N_k
		S_k = Y_k*n_k_X_k/N_k
		# Accumulate
		sum_PR += P_k * R_k
		sum_PS += P_k * S_k
		sum_QR += Q_k * R_k
		sum_QS += Q_k * S_k
		N_t += N_k
		R_t += R_k
		S_t += S_k
	# Now compute M-H odds ratio and variance
	if S_t <= 0.0:
		raise StatsError("Denominator of odds ratio is zero; can't compute M-H odds ratio or variance.")
	odds_ratio = R_t / S_t
	if R_t <= 0.0:
		raise StatsError("Denominator R_t <= 0; can't compute M-H odds ratio variance.")
	var_odds_ratio = (sum_PR/(2*R_t**2.0) + (sum_PS + sum_QR)/(2*R_t*S_t) + sum_QS/(2*S_t**2.0)) * odds_ratio**2.0

	res = MHVarResult(odds_ratio, var_odds_ratio, n_tables, N_t)
	return res

if __name__ == "__main__":
	test_pearsonCorrelation()

# End stats.py
