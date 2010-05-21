import sys, os, math, string

class BioFileError:
	"""Class for packaging errors reading biological data from files"""
	def __init__(self, s):
		return

def writeFASTA(seq_list, outfile, headers=None):
	n_written = 0
	if headers is None:
		headers = ["%d"%(j+1,) for j in range(len(seq_list))]
	for (hdr,seq) in zip(headers,seq_list):
		line = ">%s\n%s\n" % (hdr, seq)
		outfile.write(line)
		n_written += 1
	return n_written

def writeFASTADict(seq_dict, outfile):
	headers = seq_dict.keys()
	seqs = [seq_dict[k] for k in headers]
	return writeFASTA(seqs, outfile, headers=headers)

#-----------------------------------------------------------------------------------
def readFASTA(infile_name):
	"""Reads the sequences and headers from a FASTA file.

	'infile' is a FASTA file.  Reads a list of all of the headers in the FASTA file
	and all of the sequences in the FASTA file and returns them as the tuple
	(headers, sequences) where headers[i] is the header corresponding to
	sequences[i].
	Removes the '>' from the headers.
	Returns None if there is a problem processing the file."""
	infile_name = os.path.expanduser(infile_name)
	if not os.path.isfile(infile_name):
		raise BioFileError, "Cannot find the FASTA file %s." % infile_name
	infile = file(infile_name, 'r')
	seq = []
	headers = []
	sequences = []
	lines = infile.readlines()
	if len(lines)>0:
		for line in lines:  # read the lines from the file
			if line[0] == '#':
				# Skip comment
				continue
			if line[0] == '>':  # a new header
				if seq:
					frag = ''.join(seq)
					sequences.append(frag.upper())
					seq = []
				# Strip off leading '>'
				headers.append(line[1:].rstrip())
			else:
				frag = line.rstrip().upper()
				seq.append(frag)
		frag = ''.join(seq)
		sequences.append(frag.upper())
	infile.close()
	assert len(headers) == len(sequences), "Error, headers and sequences have different lengths."
	return (headers, sequences)

def firstField(x):
	h = x.split()[0].strip()
	return h

def getPeptideID(header):
	try:
		h = header.split("pep:")[1].split()[0].strip()
	except:
		h = firstField(header)
	return h

def getGeneID(header):
	try:
		h = header.split("gene:")[1].split()[0].strip()
	except:
		h = firstField(header)
	return h

def getTranscriptID(header):
	try:
		h = header.split("trans:")[1].split()[0].strip()
	except:
		h = firstField(header)
	return h

def getIDFunction(s):
	if s == 'pep' or s == 'peptide':
		return getPeptideID
	elif s == 'trans' or s == 'transcript':
		return getTranscriptID
	elif s == 'gene':
		return getGeneID
	return firstField

def readFASTADict(infile_name, key_fxn = firstField):
	#(headers, sequences) = readFASTADict(infile)
	infile_name = os.path.expanduser(infile_name)
	fdict = {}
	if not os.path.isfile(infile_name):
		raise BioFileError, "Cannot find the FASTA file %s." % infile_name
	f = file(infile_name, 'r')
	seq = []
	currHeader = ''
	for line in f:  # read the lines from the file
		if line[0] == '#':
			# Skip comment
			continue
		if line[0] == '>':  # a new header
			if seq:
				frag = ''.join(seq)
				fdict[currHeader] = frag.upper()
				seq = []
			try:
				s = line[1:].strip()
				currHeader = key_fxn(s)
			except Exception, e:
				#print line[1:].rstrip()
				continue
		else:
			frag = line.rstrip()
			seq.append(frag.upper())
			frag = ''.join(seq)
	fdict[currHeader] = frag.upper()
	f.close()
	return fdict
