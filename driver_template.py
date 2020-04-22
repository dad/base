#! python

import sys, os, math, random, argparse
import util, biofile

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Generic script template")
	# Required arguments
	parser.add_argument(dest="in_fname", default=None, help="input filename")
	# Optional arguments
	parser.add_argument("-o", "--out", dest="out_fname", default=None, help="output filename")
	options = parser.parse_args()

	info_outs = util.OutStreams(sys.stdout)
	data_outs = util.OutStreams()

	# Start up output
	if not options.out_fname is None:
		outf = open(options.out_fname,'w')
		data_outs.addStream(outf)
	else:
		# By default, write to stdout
		data_outs.addStream(sys.stdout)

	# Write out parameters
	data_outs.write("# Run started {}\n".format(util.timestamp()))
	data_outs.write("# Command: {}\n".format(' '.join(sys.argv)))
	data_outs.write("# Parameters:\n")
	optdict = vars(options)
	for (k,v) in optdict.items():
		data_outs.write("#\t{k}: {v}\n".format(k=k, v=v))

	'''
	# Read input
	if not os.path.isfile(options.in_fname):
		raise IOError("# Error: file {} does not exist".format(options.in_fname))
	with open(options.in_fname,'r') as inf:
		# Read a FASTA file?
		(headers, seqs) = biofile.readFASTA(inf)
		# Read a tab-delimited file?
		dlr = util.DelimitedLineReader(inf, header=True)
		for flds in dlr.dictentries:
			pass
	'''
	
	# Write output
	dout = util.DelimitedOutput()
	dout.addHeader('column1','A string','s')
	dout.addHeader('column2','An integer','d')
	# Write the header descriptions
	dout.describeHeader(data_outs)
	# Write the header fields
	dout.writeHeader(data_outs)
	#format = dout.getFormat(named=True)
	n_written = 0
	for entry in []:
		# A dictionary of results, one result per addHeader call above
		result = dout.createResult(default=None)
		result['column1'] = "the answer is"
		result['column2'] = 42
		# Parse the values, convert Nones to NA, etc.
		line = dout.formatLine(result)
		# A more manual approach:
		# line = format.format(column1="the answer is", column2=42)
		data_outs.write(line)
		n_written += 1

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.out_fname))
		outf.close()

