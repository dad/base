#! python

import sys, os, math, random, datetime, argparse
import util, biofile, translate, protprop

# pK values from http://helixweb.nih.gov/emboss/html/iep.html
# pI code from Peter Collingridge, http://www.petercollingridge.co.uk/sites/files/peter/predictPI.txt


if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Print charge profile for protein sequence")
	parser.add_argument("-o", "--out", dest="out_fname", default=None, help="output filename")
	parser.add_argument("-i", "--in", dest="in_fname", default=None, help="input FASTA filename")
	parser.add_argument("-s", "--seq", dest="sequence", default=None, help="input sequence")
	parser.add_argument("-w", "--window", dest="window", type=int, default=7, help="window over which to profile")
	parser.add_argument("-t", "--translate", dest="translate", action="store_true", default=False, help="translate the input sequences?")
	parser.add_argument("--pH", dest="pH", type=float, default=7.2, help="pH for charge determination")
	options = parser.parse_args()
	
	outs = util.OutStreams()
	if not options.out_fname is None:
		outf = file(os.path.expanduser(options.out_fname),'w')
		outs.addStream(outf)
	else:
		outs.addStream(sys.stdout)
	
	pp = protprop.ProteinProperties()
	if not options.sequence is None:
		if options.translate:
			seq = translate.translateRaw(options.sequence)
		else:
			seq = options.sequence
		seq_dict = {"input":seq}
	else:
		# Load from FASTA
		seq_dict = biofile.readFASTADict(options.in_fname)
		if options.translate:
			for k in seq_dict.keys():
				seq_dict[k] = translate.translate(seq_dict[k])

	outs.write("# {}\n".format(options))
	outs.write("pos\taa\tcharge\n")
	n_seqs = len(seq_dict.keys())
	for (seqid, seq) in seq_dict.items():
		if n_seqs>1:
			outs.write("# {}\n".format(seqid))
			outs.write("# Total protein charge at pH {} = {}\n".format(options.pH, pp.getCharge(seq, options.pH)))
		window_width = (options.window - 1)/2
		# Run over
		start_pos = 0
		focal_pos = 0
		end_pos = min(start_pos+options.window, len(seq))
		while focal_pos < len(seq):
			start_pos = max(0,focal_pos-window_width)
			end_pos = min(len(seq),focal_pos+window_width+1)
			win_seq = seq[start_pos:end_pos]
			#print start_pos, end_pos, win_seq
			if len(win_seq)>0:
				win_charge = pp.getCharge(win_seq, options.pH)/float(len(win_seq))
				line = "{pos}\t{aa}\t{charge}\n".format(pos=focal_pos+1, aa=seq[focal_pos], charge=win_charge)
				outs.write(line)
			focal_pos += 1
			
		
