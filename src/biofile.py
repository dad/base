import sys, os, math, string

class BioFileError:
	"""Class for packaging errors reading biological data from files"""
	def __init__(self, s):
		return

'''
From http://genome.ucsc.edu/FAQ/FAQformat.html#format3
GFF (General Feature Format) lines are based on the GFF standard file format. GFF lines have nine required fields that must be tab-separated. If the fields are separated by spaces instead of tabs, the track will not display correctly. For more information on GFF format, refer to http://www.sanger.ac.uk/resources/software/gff.

Here is a brief description of the GFF fields:

seqname - The name of the sequence. Must be a chromosome or scaffold.
source - The program that generated this feature.
feature - The name of this type of feature. Some examples of standard feature types are "CDS", "start_codon", "stop_codon", and "exon".
start - The starting position of the feature in the sequence. The first base is numbered 1.
end - The ending position of the feature (inclusive).
score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). If there is no score value, enter ".".
strand - Valid entries include '+', '-', or '.' (for don't know/don't care).
frame - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be '.'.
group - All lines with the same group are linked together into a single item.
Example:
Here's an example of a GFF-based track. Click here for a copy of this example that can be pasted into the browser without editing. NOTE: Paste operations on some operating systems will replace tabs with spaces, which will result in an error when the GFF track is uploaded. You can circumvent this problem by pasting the URL of the above example (http://genome.ucsc.edu/goldenPath/help/regulatory.txt) instead of the text itselfinto the custom annotation track text box.

track name=regulatory description="TeleGene(tm) Regulatory Regions"
chr22  TeleGene enhancer  1000000  1001000  500 +  .  touch1
chr22  TeleGene promoter  1010000  1010100  900 +  .  touch1
chr22  TeleGene promoter  1020000  1020000  800 -  .  touch2
'''

class GFFRecord:
	def __init__( self ):
		self.seqname = ''
		self.source = ''
		self.feature = ''
		self.start = 0
		self.end = 0
		self.score = 0
		self.strand = '+'
		self.frame = '.'
		self.group = ''
		self._new_start = None
		self._new_end = None

	def read(self, line):
		flds = line.strip().split('\t')
		self.readFromFields(flds)

	def readFromFields(self, flds):
		self.seqname = flds[0]
		self.source = flds[1]
		self.feature = flds[2]
		self.start = int(flds[3])
		self.end = int(flds[4])
		self.score = flds[5]
		self.strand = flds[6]
		self.frame = flds[7]
		self.group = flds[8]

	def update(self, ref_site, net_difference):
		at_end = False
		# ref_site is 0-based index, start/end are 1-based indices
		if self.start-1 <= ref_site and self._new_start is None:
			self._new_start = self.start + net_difference
		elif self.end-1 <= ref_site and self._new_end is None:
			self._new_end = self.end + net_difference
			at_end = True
			self.end = self._new_end
			self.start = self._new_start
		return at_end

	def __str__(self):
		return "{x.seqname}\t{x.source}\t{x.feature}\t{x.start}\t{x.end}\t{x.score}\t{x.strand}\t{x.frame}\t{x.group}\n".format(x=self)

	def write(self, stream):
		stream.write(str(self))

class GFFRecordTracker:
	def __init__(self, gff_rec):
		self.rec = gff_rec
		self.new_start = gff_rec.start
		self.new_end = gff_rec.end
		self.cur_site = 0
	
	def startTracking(self, cur_site):
		self.cur_site = cur_site
	
	def stopTracking(self):
		self.rec.end = self.new_end
		self.rec.start = self.new_start
		return rec

	def update(self, ref_site, net_difference):
		if self.start >= ref_site and self.new_start is None:
			self._new_start = self.start + net_difference
		elif self.end <= ref_site and self.new_end is None:
			self.new_end = self.end + net_difference
			at_end = True
			self.end = self.new_end
			self.start = self.new_start
		return at_end
'''
Col	Field	Description
1	CHROM	Chromosome name
2	POS	1-based position. For an indel, this is the position preceding the indel.
3	ID	Variant identifier. Usually the dbSNP rsID.
4	REF	Reference sequence at POS involved in the variant. For a SNP, it is a single base.
5	ALT	Comma delimited list of alternative seuqence(s).
6	QUAL	Phred-scaled probability of all samples being homozygous reference.
7	FILTER	Semicolon delimited list of filters that the variant fails to pass.
8	INFO	Semicolon delimited list of variant information.
9	FORMAT	Colon delimited list of the format of individual genotypes in the following fields.
10+	Sample(s)	Individual genotype information defined by FORMAT.

##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	mdy748-raw-reads-cleaned-map.sorted.bam
ref|NC_001133|	19	.	CCCAC	CCACACCAC	99	.	INDEL;DP=1623;AF1=0.5;CI95=0.5,0.5;DP4=78,42,73,14;MQ=60;PV4=0.0026,1,1,0.011	PL:GT:GQ	255,0,255:0/1:99
ref|NC_001133|	24	.	ACACCACACCAC	ACCCACCACACCAC	11.8	.	INDEL;DP=1747;AF1=0.5;CI95=0.5,0.5;DP4=109,52,32,2;MQ=60;PV4=0.0012,1,1,1	PL:GT:GQ	49,0,255:0/1:52
ref|NC_001133|	25	.	CACCACA	CCA	99	.	INDEL;DP=1736;AF1=0.5;CI95=0.5,0.5;DP4=92,51,16,14;MQ=60;PV4=0.3,1,1,1	PL:GT:GQ	255,0,222:0/1:99
'''
class VCFRecord:
	def __init__( self ):
		self.chromosome = ''
		self.position = 0
		self.id = ''
		self.ref = 0
		self.alt = 0
		self.qual = 0
		self.filter = '+'
		self.info = '.'
		self.format = ''
		self.other = None

	def read(self, line):
		flds = line.strip().split('\t')
		self.readFromFields(flds)

	def readFromFields(self, flds):
		self.chromosome = flds[0]
		self.position = int(flds[1])
		self.id = flds[2]
		self.ref = flds[3]
		self.alt = flds[4]
		self.qual = float(flds[5])
		self.filter = flds[6]
		self.info = flds[7]
		self.format = flds[8]
		if len(flds)>9:
			self.other = '/'.join(flds[9:])

	def __str__(self):
		s = "{x.chromosome}\t{x.position}\t{x.id}\t{x.ref}\t{x.alt}\t{x.qual}\t{x.filter}\t{x.info}\t{x.format}".format(x=self)
		if not self.other is None:
			s += '\t{0}'.format(self.other)
		return s

	def write(self, stream):
		stream.write(str(self))

	def getReferenceStart(self):
		return self.position

	def getReferenceEnd(self):
		return self.position + len(self.ref)

	start = property(getReferenceStart, None, None, None)
	end = property(getReferenceEnd, None, None, None)

#-----------------------------------

def writeFASTA(seq_list, outfile, headers=None):
	n_written = 0
	if headers is None:
		headers = ["{0}".format(j+1) for j in range(len(seq_list))]
	for (hdr,seq) in zip(headers,seq_list):
		line = ">{0}\n{1}\n".format(hdr, seq)
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
	seq = None
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
				h = line[1:].strip()
				headers.append(h)
				#print h
				#print line
			else:
				frag = line.strip()
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

def getFlybaseGeneID(header):
	flds = header.split()
	id = flds[0]
	flybase_dict = dict([tuple(x[:-1].split('=')) for x in flds[1:]])
	#print flybase_dict
	gene_ids = [x for x in flybase_dict["parent"].split(',') if x.startswith('FBgn')]
	return gene_ids[0]

def getFlybaseTranscriptID(header):
	flds = header.split()
	id = flds[0]
	flybase_dict = dict([tuple(x[:-1].split('=')) for x in flds[1:]])
	gene_ids = [x for x in flybase_dict["parent"].split(',') if x.startswith('FBtr')]
	return gene_ids[0]

def getIDFunction(s):
	if s == 'pep' or s == 'peptide':
		return getPeptideID
	elif s == 'trans' or s == 'transcript':
		return getTranscriptID
	elif s == 'gene':
		return getGeneID
	elif s == 'flybase-gene':
		return getFlybaseGeneID
	elif s == 'flybase-transcript':
		return getFlybaseTranscriptID
	return firstField

def readFASTADict(infile_name, key_fxn = firstField):
	fdict = {}
	(headers, sequences) = readFASTA(infile_name)
	for (i,h) in enumerate(headers):
		h_key = key_fxn(h)
		s = sequences[i]
		fdict[h_key] = s
	return fdict
