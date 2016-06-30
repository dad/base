import sys, os, math, string, argparse
import util, datetime, time, biofile
import protprop

if __name__=='__main__':
	parser = argparse.ArgumentParser() #usage="%prog [-i fasta] [-s sequence]")
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
	parser.add_argument("-m", "--merge", dest="merge", action="store_true", default=False, help="merge all sequences together?")
	parser.add_argument("--hydrophobicity-scale", dest="hydrophobicity_scale", type=str, default='Kyte-Doolittle', help="which hydrophobicity scale to use (Kyte-Doolittle, Hopp-Wood, Cornette, Eisenberg, Rose, Janin, Engelman-GES)")
	parser.add_argument("--pH", dest="pH", type=float, default=7.2, help="pH for charge determination")
	#parser.add_argument("-r", "--report", dest="report", action="store_true", default=False, help="write long report per protein?")
	options = parser.parse_args()
	
	outs = util.OutStreams()
	if not options.out_fname is None:
		outf = open(os.path.expanduser(options.out_fname),'w')
		outs.addStream(outf)
	else:
		outs.addStream(sys.stdout)

	# Write parameters	
	outs.write("# Run {}\n".format(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')))
	outs.write("# Parameters:\n")
	optdict = vars(options)
	for (k,v) in optdict.items():
		outs.write("#\t{k}: {v}\n".format(k=k, v=v))

	
	pp = protprop.ProteinProperties()
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
		(headers,seqs) = biofile.readFASTA(open(options.in_fname, 'r'))
	
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
	if options.merge:
		outs.write("# Merging {:d} sequences into one\n".format(len(seqs)))
		seqs = [''.join(seqs)]
		headers = ["merged"]
	gap = '-'
	for (h,seq) in zip(headers,seqs):
		if options.query:
			if not h.strip().startswith(options.query):
				continue
		if options.translate:
			seq = translate.translateRaw(seq)
		if options.degap:
			seq = seq.replace(gap,'')
		if not options.exclude:
			if not options.end_aa is None and options.end_aa <= len(seq):
				seq = seq[0:(options.end_aa)]
			seq = seq[(options.begin_aa-1):]
		else: # Exclude the sequence
			assert options.end_aa < len(seq)
			assert options.begin_aa < options.end_aa
			seq = seq[0:(options.begin_aa-1)] + seq[(options.end_aa):]
		degapped_seq = seq.replace(gap,"")
		line = "#{}\n{}\t{:d}\t{:1.4f}\t{:1.4f}\t{:1.4f}".format(h, biofile.firstField(h), pp.getLength(degapped_seq), pp.getCharge(degapped_seq, options.pH), pp.getIsoelectricPoint(degapped_seq), pp.getHydrophobicity(degapped_seq))
		if not aas is None:
			counts = protprop.Composition()
			counts.initFromSequence(degapped_seq)
			freqs = protprop.Composition()
			freqs.initFromSequence(degapped_seq)
			freqs.normalize()
			line += '\t' + '\t'.join(["{:1.4f}".format(freqs[aa]) for aa in aas]) + '\t' + '\t'.join(["{:d}".format(counts[aa]) for aa in aas])
		outs.write(line + '\n')
	if not options.out_fname is None:
		outf.close()
