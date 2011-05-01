import sys, os, math, string, random
import cai, translate, util, biofile
from optparse import OptionParser

if __name__=='__main__':
	parser = OptionParser(usage="%prog [options] sequence")
	parser.add_option("-o", "--out", dest="out_fname", type="string", default=None, help="output filename")
	parser.add_option("-r", "--reverse-translate", dest="reverse_translate", action="store_true", help="interpret sequences as amino acids, and convert back to codons?")
	parser.add_option("--species", dest="species", type="string", default="scer", help="species for identifying optimal codons")
	parser.add_option("--min-ra", dest="min_rel_adapt", type="float", default=1.0, help="minimum relative adaptiveness for optimization")
	parser.add_option("-p", "--optimize", dest="optimize", action="store_true", default=False, help="optimize the codons?")

	(options, args) = parser.parse_args()
	seq = args[0]

	info_outs = util.OutStreams([sys.stdout])
	data_outs = util.OutStreams()
	multi_outs = util.OutStreams([info_outs, data_outs])
	if not options.out_fname is None:
		outf = file(options.out_fname,'w')
		data_outs.addStream(outf)
	else:
		data_outs.addStream(sys.stdout)

	fname = os.path.expanduser(seq)
	if os.path.isfile(fname):
		(headers, seq_list) = biofile.readFASTA(fname)
		seqs = zip(headers,seq_list)
		info_outs.write("# Read %d sequences from %s\n" % (len(seqs), fname))
	else:
		seqs = [("<input>",seq)]
	if options.reverse_translate:
		new_seqs = []
		for (h,s) in seqs:
			rev_trans_seq = translate.reverseTranslate(s)
			new_seqs.append((h,rev_trans_seq))
		seqs = new_seqs

	# Obtain gene sequence using only optimal codons
	info_outs.write("# Using optimal codons for %s\n" % options.species)
	opt_codons = cai.getOptimalCodons(options.species)
	relad_dict = cai.getRelativeAdaptivenessValues(options.species)
	ln_relad_dict = dict([(k,math.log(v+0.001)) for (k,v) in relad_dict.items()])
	for (id, seq) in seqs:
		info_outs.write("# %s Fop = %1.4f, CAI = %1.4f, GC = %1.2f\n" % (id, cai.getFop(seq, opt_codons), cai.getCAI(seq, ln_relad_dict), cai.getGC(seq)))

	if options.optimize:
		info_outs.write("# Optimizing sequences...\n")
		gc = translate.geneticCode()
		codons = {}
		opt_codon_dict = dict([(gc[c],c) for c in opt_codons])
		opt_codon_dict['W'] = 'TGG'
		opt_codon_dict['M'] = 'ATG'

		opt_headers = []
		opt_seqs = []
		# optimize the codon sequences
		for (id, seq) in seqs:
			prot_seq = translate.translate(seq)
			if not prot_seq is None:
				for aa in translate.AAs():
					codons[aa] = [c for c in translate.getCodonsForAA(aa, rna=False) if relad_dict[c] >= options.min_rel_adapt]
				opt_seq = ''
				for aa in prot_seq:
					#opt_seq += opt_codon_dict[aa] #random.choice(codons[aa])
					opt_seq += random.choice(codons[aa])
				assert translate.translate(opt_seq) == prot_seq
				header_line = "%s Fop = %1.4f, CAI = %1.4f, GC = %1.2f" % (id, cai.getFop(opt_seq, opt_codons), cai.getCAI(opt_seq, ln_relad_dict), cai.getGC(opt_seq))
				info_outs.write("# Optimized %s\n" % header_line)
				opt_headers.append(header_line)
				opt_seqs.append(opt_seq)

		biofile.writeFASTA(opt_seqs, data_outs, headers=opt_headers)
		if not options.out_fname is None:
			info_outs.write("# Wrote %d optimized sequences to %s\n" % (len(opt_seqs), options.out_fname))
			outf.close()

