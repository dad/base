#! python

import sys, os, math, random, datetime
from optparse import OptionParser
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

	def getCharge(self, sequence, pH):
		protein_charge = self._aminoAcidCharge('N-term', pH) + self._aminoAcidCharge('C-term', pH)
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
	
	def getComposition(self, sequence):
		std_aas = translate.AAs()
		aas = std_aas + ''.join(sorted(list(set([aa for aa in sequence if not aa in std_aas]))))
		aa_counts = [(aa,sequence.count(aa)) for aa in aas]
		return aa_counts
	
	def getLength(self, sequence, stopchr="*", gapchr="-"):
		# If there's a "*" -- usual way to put in a stop codon -- should we just report the length up to the stop?
		# Let's do that.
		if stopchr in sequence:
			seq = sequence[0:sequence.rfind(stopchr)]
		else:
			seq = sequence
		res = len(seq.replace(gapchr,''))
		return res
		

if __name__=='__main__':
	parser = OptionParser(usage="%prog [-i fasta] [-s sequence]")
	parser.add_option("-o", "--out", dest="out_fname", type="string", default=None, help="output filename")
	parser.add_option("-i", "--in", dest="in_fname", type="string", default=None, help="input FASTA filename")
	parser.add_option("-s", "--seq", dest="sequence", type="string", default=None, help="input sequence")
	parser.add_option("-t", "--translate", dest="translate", action="store_true", default=False, help="translate the input sequences?")
	parser.add_option("--pH", dest="pH", type="float", default=7.2, help="pH for charge determination")
	(options, args) = parser.parse_args()
	
	outs = util.OutStreams()
	if not options.out_fname is None:
		outf = file(os.path.expanduser(options.out_fname),'w')
		outs.addStream(outf)
	else:
		outs.addStream(sys.stdout)
	
	pp = ProteinProperties()
	if not options.sequence is None:
		if options.translate:
			seq = translate.translateRaw(options.sequence)
		else:
			seq = options.sequence
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
	
	if (not options.in_fname is None):# and os.path.isfile(os.path.expanduser(options.in_fname)):
		(headers,seqs) = biofile.readFASTA(file(options.in_fname,'r'))
		outs.write("# {}\n# pH={}\n".format(datetime.datetime.now().strftime("%A, %d %B %Y %I:%M%p"), options.pH))
		outs.write("orf\tlength\tcharge\tpI\thydrophobicity\n")
		for (h,seq) in zip(headers,seqs):
			if options.translate:
				seq = translate.translateRaw(seq)
			line = "#{}\n{}\t{:d}\t{:1.4f}\t{:1.4f}\t{:1.4f}\n".format(h, biofile.firstField(h), pp.getLength(seq), pp.getCharge(seq, options.pH), pp.getIsoelectricPoint(seq), pp.getHydrophobicity(seq))
			outs.write(line)

	if not options.out_fname is None:
		outf.close()

	