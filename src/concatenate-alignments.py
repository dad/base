#! python

import sys, os
from argparse import ArgumentParser
import util, biofile

if __name__=='__main__':
	parser = ArgumentParser(description="Concatenate two or more FASTA files entry-wise")
	parser.add_argument(dest="in_fname", default=None, help="input FASTA filename")

	# Optional arguments
	parser.add_argument("-o", "--out", dest="out_fname", default=None, help="output filename")
	parser.add_argument("-a", "--append", dest="append_fasta", action='append', default=[], help="input FASTA file to append entry-wise")	
	parser.add_argument("--check-headers", dest="check_headers", action="store_true", default=False, help="compare headers before merging?")
	parser.add_argument("--fasta-out", dest="fasta_out_fname", default=None, help="output FASTA filename")
	options = parser.parse_args()
	
	info_outs = util.OutStreams(sys.stdout)
	data_outs = util.OutStreams()
	fasta_outs = util.OutStreams()

	# Start up output
	if not options.out_fname is None:
		outf = open(options.out_fname,'w')
		data_outs.addStream(outf)
	else:
		# By default, write to stdout
		data_outs.addStream(sys.stdout)
	if not options.fasta_out_fname is None:
		outf = open(options.fasta_out_fname,'w')
		fasta_outs.addStream(outf)
	else:
		# By default, write to stdout
		fasta_outs.addStream(sys.stdout)

	# Write out parameters
	data_outs.write("# Run started {}\n".format(util.timestamp()))
	data_outs.write("# Command: {}\n".format(' '.join(sys.argv)))
	data_outs.write("# Parameters:\n")
	optdict = vars(options)
	for (k,v) in optdict.items():
		data_outs.write("#\t{k}: {v}\n".format(k=k, v=v))

	# Read input
	if not os.path.isfile(options.in_fname):
		raise IOError("# Error: file {} does not exist".format(options.in_fname))
	with open(options.in_fname,'r') as inf:
		# Read a FASTA file?
		(headers, seqs) = biofile.readFASTA(inf)
	info_outs.write("# Read {:d} sequences\n".format(len(seqs)))

	new_headers = headers
	new_seqs = seqs
	for append_fname in options.append_fasta:
		if not os.path.isfile(append_fname):
			raise IOError("# Error: file {} does not exist".format(append_fname))
		with open(append_fname,'r') as inf:
			# Read a FASTA file?
			(app_headers, app_seqs) = biofile.readFASTA(inf)
			info_outs.write("# Read {:d} sequences\n".format(len(app_seqs)))
			assert len(app_seqs) == len(new_seqs)
			if options.check_headers:
				for (h1, h2) in zip(headers, app_headers):
					assert h1==h2, "# Error: headers do not match:\n\t{}\n\t{}".format(h1,h2)
			for (i,s) in enumerate(new_seqs):
				new_seqs[i] += app_seqs[i]

	# Write output
	biofile.writeFASTA(new_seqs, fasta_outs, headers=new_headers)
	n_written = len(new_seqs)

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.fasta_out_fname is None:
		info_outs.write("# Wrote {} entries to {}\n".format(n_written, options.fasta_out_fname))
		outf.close()

	
 