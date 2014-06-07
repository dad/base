#! python

import sys, os, math, datetime, time
from argparse import ArgumentParser
import util, biofile, translate

# pK values from http://helixweb.nih.gov/emboss/html/iep.html
# pI code from Peter Collingridge, http://www.petercollingridge.co.uk/sites/files/peter/predictPI.txt

class ProteinProperties(object):
	def __init__(self):
		self.pKa     = {'D':3.9, 'E':4.3, 'H':6.1, 'C':8.3, 'Y':10.1, 'K':10.67, 'R':12, 'N-term':8, 'C-term':3.1}
		self.charges = {'D':-1,  'E':-1,  'H':1,  'C':-1,  'Y':-1,   'K':1,    'R':1,  'N-term':1, 'C-term':-1}
		#self.charges = {'D':-1, 'E':-1, 'H':1, 'K':1, 'R':1, 'N-term':1, 'C-term':-1}
		self.hydrophobicity_scales = {
			'Kyte-Doolittle':{'A':1.800,'R':-4.500,'N':-3.500,'D':-3.500,'C':2.500,'Q':-3.500,'E':-3.500,'G':-0.400,'H':-3.200,'I':4.500,'L':3.800,'K':-3.900,'M':1.900,'F':2.800,'P':-1.600,'S':-0.800,'T':-0.700,'W':-0.900,'Y':-1.300,'V':4.200}
		}

	def _aminoAcidCharge(self, amino_acid, pH):
		proportion = 1 / (1 + 10**(pH - self.pKa[amino_acid]))
		res = None
		if self.charges[amino_acid] == 1:
			res = proportion
		else:
			res = proportion-1.0 # more clearly, -1 * (1-proportion)
		return res

	def getCharge(self, sequence, pH, include_termini=True):
		protein_charge = 0.0
		if include_termini:
			protein_charge += self._aminoAcidCharge('N-term', pH) + self._aminoAcidCharge('C-term', pH)
		for aa in self.charges.keys():
			protein_charge += sequence.count(aa) * self._aminoAcidCharge(aa, pH)
		return protein_charge

	def getIsoelectricPoint(self, sequence, tolerance=1e-4):
		min_pH, max_pH = 3, 13
		done = False
		n_iterations = 0
		max_iterations = 1000
		# Binary search for the pH at which the net charge is within the specified tolerance around zero.
		while not done and n_iterations<max_iterations:
			mid_pH = 0.5 * (max_pH + min_pH)
			protein_charge = self.getCharge(sequence, mid_pH)
			if protein_charge > tolerance:
				min_pH = mid_pH
			elif protein_charge < -tolerance:
				max_pH = mid_pH
			else:
				done = True
			n_iterations += 1
		return mid_pH
	
	def getHydrophobicity(self, sequence, scale='Kyte-Doolittle'):
		hyd_scale = self.hydrophobicity_scales[scale]
		hyd = 0.0
		n = 0
		res = 0
		for aa in hyd_scale.keys():
			naa = sequence.count(aa)
			hyd += naa*hyd_scale[aa]
			n += naa
		if n > 0:
			res = hyd/n
		return res
	
	def getComposition(self, sequence, normalize=False, aas=translate.AAs()):
		#aas = translate.AAs()
		if aas is None:
			aas = ''
		#seq_aas = aas + ''.join(sorted(list(set([aa for aa in sequence if not aa in aas]))))
		aa_counts = [(aa,sequence.count(aa)) for aa in aas]
		res = aa_counts
		if normalize:
			tot = float(sum([c for (aa,c) in aa_counts]))
			res = [(aa,c/tot) for (aa,c) in aa_counts]
		return res
	
	def getLength(self, sequence, stopchr="*", gapchr="-"):
		# If there's a "*" -- usual way to put in a stop codon -- should we just report the length up to the stop?
		# Let's do that.
		if stopchr in sequence:
			seq = sequence[0:sequence.rfind(stopchr)]
		else:
			seq = sequence
		res = len(seq.replace(gapchr,''))
		return res

	def count(self, sequence, aas):
		the_count = 0
		for aa in sequence:
			if aa in aas:
				the_count += 1
		return the_count

	def motif(self, sequence, aa_classes, symbol_map=chr):
		aa_dict = {}
		for (xi,aas) in enumerate(aa_classes):
			for aa in aas:
				aa_dict[aa] = ord('a')+xi
		mot = ''
		for aa in sequence:
			try:
				mot += symbol_map(aa_dict[aa])
			except KeyError:
				pass
		return mot

class Composition(object):
	def __init__(self):
		aas = translate.AAs()
		self._comp_dict = dict([(aa,0.0) for aa in aas]) # list of (aa, frequency) tuples
	
	def initFromList(self, tuple_list):
		self._comp_dict = dict(tuple_list)
	
	def initFromSequence(self, seq):
		pp = ProteinProperties()
		comp = pp.getComposition(seq, aas=translate.AAs())
		self._comp_dict = dict(comp)
	
	def normalize(self):
		tot = float(sum([c for (aa,c) in self._comp_dict.items()]))
		self._comp_dict = dict([(aa,c/tot) for (aa,c) in self._comp_dict.items()])
	
	def write(self, stream, header=True):
		if header:
			stream.write("aa\tproportion\n")
		for aa in sorted(self._comp_dict.keys()):
			write("{:s}\t{:1.4f}\n".format(aa, self._comp_dict[aa]))
	
	def read(self, stream, header=True):
		tab = util.readTable(stream, header=header)
		for flds in tab.dictrows:
			self._comp_dict[flds['aa']] = flds['proportion']
	
	def __getitem__(self, aa):
		return self._comp_dict[aa]
			


if __name__=='__main__':
	parser = ArgumentParser() #usage="%prog [-i fasta] [-s sequence]")
	parser.add_argument("-o", "--out", dest="out_fname", default=None, help="output filename")
	parser.add_argument("-i", "--in", dest="in_fname", default=None, help="input FASTA filename")
	parser.add_argument("-s", "--seq", dest="sequence", default=None, help="input sequence")
	parser.add_argument("-t", "--translate", dest="translate", action="store_true", default=False, help="translate the input sequences?")
	parser.add_argument("-b", "--begin", dest="begin_aa", type=int, default=1, help="beginning amino acid (1-based)")
	parser.add_argument("-e", "--end", dest="end_aa", type=int, default=None, help="ending amino acid (1-based, inclusive)")
	parser.add_argument("-x", "--exclude", dest="exclude", action="store_true", default=False, help="exclude rather than include begin/end region?")
	parser.add_argument("-a", "--aas", dest="aas", default=None, help="amino acids for frequency counts")
	parser.add_argument("-g", "--degap", dest="degap", action="store_true", default=False, help="remove gaps before applying begin/end coordinates?")
	parser.add_argument("-q", "--query", dest="query", default=None, help="specific sequence identifier to query")
	parser.add_argument("--pH", dest="pH", type=float, default=7.2, help="pH for charge determination")
	#parser.add_argument("-r", "--report", dest="report", action="store_true", default=False, help="write long report per protein?")
	options = parser.parse_args()
	
	outs = util.OutStreams()
	if not options.out_fname is None:
		outf = file(os.path.expanduser(options.out_fname),'w')
		outs.addStream(outf)
	else:
		outs.addStream(sys.stdout)

	# Write parameters	
	outs.write("# Run {}\n".format(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')))
	outs.write("# Parameters:\n")
	optdict = vars(options)
	for (k,v) in optdict.items():
		outs.write("#\t{k}: {v}\n".format(k=k, v=v))

	
	pp = ProteinProperties()
	aas = None
	if not options.aas is None:
		if options.aas.lower() == 'all':
			aas = translate.AAs()
		else:
			aas = [aa for aa in options.aas]
	
	# Single sequence?
	if not options.sequence is None:
		headers = ['Input']
		seqs = [options.sequence]
	else:
		(headers,seqs) = biofile.readFASTA(file(options.in_fname, 'r'))
	
	'''
	if options.report: # Write a long report per protein
		for (hdr, seq) in zip(headers,seqs):
			if options.degap:
				seq = seq.replace('-','')
			if not options.end_aa is None and options.end_aa<= len(seq):
				seq = seq[0:options.end_aa]
			#print options.end_aa, options.begin_aa
			seq = seq[options.begin_aa:]
			outs.write("length = {:d}\n".format(pp.getLength(seq)))
			pI = pp.getIsoelectricPoint(seq, tolerance=1e-4)
			outs.write("pI = {0}\n".format(pI))
			outs.write("charge at pH {0:1.1f} = {1:1.2f}\n".format(options.pH, pp.getCharge(seq, options.pH)))
			outs.write("charge at pH=pI = {:1.2f}\n".format(pp.getCharge(seq, pI)))
			scale = 'Kyte-Doolittle'
			outs.write("hydrophobicity ({}) = {:1.2f}\n".format(scale, pp.getHydrophobicity(seq, scale)))
			outs.write("composition:\n\taa\tcount\tproportion\n")
			for (aa,ct) in pp.getComposition(seq):
				outs.write("\t{}\t{}\t{:1.2f}\n".format(aa,ct,float(ct)/len(seq)))
	else: # Write compact
	'''
	
	outs.write("orf\tlength\tcharge\tpI\thydrophobicity")
	if not aas is None:
		outs.write("\t"+"\t".join(["f.{}".format(aa) for aa in aas])) # fractions
		outs.write("\t"+"\t".join(["n.{}".format(aa) for aa in aas])) # numbers
	outs.write("\n")
	for (h,seq) in zip(headers,seqs):
		if options.query:
			if not h.strip().startswith(options.query):
				continue
		if options.translate:
			seq = translate.translateRaw(seq)
		if options.degap:
			seq = seq.replace('-','')
		#print seq
		#print options.begin_aa, options.end_aa, len(seq)
		test_seq = ''
		if not options.exclude:
			if not options.end_aa is None and options.end_aa <= len(seq):
				seq = seq[0:(options.end_aa)]
			#print seq, seq[(options.begin_aa-1):]
			seq = seq[(options.begin_aa-1):]
		else: # Exclude the sequence
			assert options.end_aa < len(seq)
			assert options.begin_aa < options.end_aa
			seq = seq[0:(options.begin_aa-1)] + seq[(options.end_aa):]
		degapped_seq = seq.replace("-","")
		print degapped_seq
		#print test_seq
		#print seq
		line = "#{}\n{}\t{:d}\t{:1.4f}\t{:1.4f}\t{:1.4f}".format(h, biofile.firstField(h), pp.getLength(degapped_seq), pp.getCharge(degapped_seq, options.pH), pp.getIsoelectricPoint(degapped_seq), pp.getHydrophobicity(degapped_seq))
		if not aas is None:
			counts = dict(pp.getComposition(degapped_seq, aas))
			line += '\t' + '\t'.join(["{:1.4f}".format(counts[aa]/float(len(degapped_seq))) for aa in aas]) + '\t' + '\t'.join(["{:d}".format(counts[aa]) for aa in aas])
		outs.write(line + '\n')
	#print seqs
	if not options.out_fname is None:
		outf.close()

	
