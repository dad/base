import sys, os, math, string
import util, biofile, geneutil


def findSequencePositions(start_position, end_position, start_sequence, end_sequence, sentinel, headers, seqs):
	"""Find sequence positions. If start_position is a number, it is returned. Otherwise, start_sequence is searched for."""
	res_start_position = None
	res_end_position = None

	def found_sentinel(x):
		return sentinel is None or (sentinel in x)

	if not start_position is None:
		res_start_position = start_position
	else:
		if not start_sequence is None:
			# Find this sequence
			for (h,s) in zip(headers,seqs):
				if found_sentinel(h):
					ind = geneutil.gappedFind(s, start_sequence, start=True)
					if ind > 0:
						res_start_position = ind
						break
	if not end_position is None:
		end_position = end_position
	else:
		if not end_sequence is None:
			# Find this sequence
			for (h,s) in zip(headers,seqs):
				if found_sentinel(h):
					endpos = geneutil.gappedFind(s, end_sequence, start=False)
					if endpos > 0:
						res_end_position = endpos
						break

	return res_start_position, res_end_position
