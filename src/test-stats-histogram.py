#! /usr/local/bin/python

import sys, os, math, string, random, pickle
import stats

def main():
	h = stats.Histogram()
	bin_min = -1.3
	bin_max = 10
	h.init(bin_min, bin_max, 10)
	for i in range(10000):
		x = random.random()*(bin_max-bin_min-2)+bin_min
		h.add(x)
	print h

	h = stats.Histogram()
	h.init(bin_min, bin_max, 10)
	try:
		h.add("waah")
	except TypeError, te:
		print te
	print h

	h = stats.Histogram()
	print h

main()