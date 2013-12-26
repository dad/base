#! python

import sys, os, math, random, argparse, time, datetime
import util

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Generic script template")
	# Required argument
	parser.add_argument(dest="in_fname", default=None, help="input filename")
	# Optional argument
	parser.add_argument("-o", "--out", dest="out_fname", default=None, help="output filename")
	options = parser.parse_args()

	# Assign required arguments here

	info_outs = util.OutStreams(sys.stdout)
	data_outs = util.OutStreams()

	# Start up output
	if not options.out_fname is None:
		outf = file(options.out_fname,'w')
		data_outs.addStream(outf)
	else:
		# By default, write to stdout
		data_outs.addStream(sys.stdout)

	# Write out parameters
	data_outs.write("# Run {}\n".format(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')))
	data_outs.write("# Parameters:\n")
	optdict = vars(options)
	for (k,v) in optdict.items():
		data_outs.write("#\t{k}: {v}\n".format(k=k, v=v))

	# Read input
	if not options.in_fname is None:
		assert os.path.isfile(options.in_fname), "# Error: file {} does not exist".format(options.in_fname)
		inf = file(options.in_fname, 'r')
		dlr = util.DelimitedLineReader(inf, header=True)
		for flds in dlr.dictentries:
			pass
		inf.close()

	# Shut down output
	if not options.out_fname is None:
		outf.close()

