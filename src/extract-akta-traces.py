#! python

import sys, os, math, random, argparse
from pycorn import pc_res3
import util, biofile



if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Extract tab-delimited data from AKTA .res files")
	# Required arguments
	parser.add_argument(dest="in_fname", default=None, help="input filename")
	# Optional arguments
	parser.add_argument("-o", "--out", dest="out_fname", default=None, help="output filename")
	options = parser.parse_args()

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
	data_outs.write("# Run started {}\n".format(util.timestamp()))
	data_outs.write("# Command: {}\n".format(' '.join(sys.argv)))
	data_outs.write("# Parameters:\n")
	optdict = vars(options)
	for (k,v) in optdict.items():
		data_outs.write("#\t{k}: {v}\n".format(k=k, v=v))

	# Create the instance using default options
	resfile = pc_res3(options.in_fname)

	# Parse the file
	resfile.load()

	'''
	# List columns
	print resfile.keys()
	for k in resfile.keys():
		print len(resfile[k]['data'])
	'''

	target_columns = ['UV'] #,'Cond','Pressure','Flow','Temp']
	target_column = 'UV'
	#print resfile[target_column].keys()
	#sys.exit()


	# Write output
	dout = util.DelimitedOutput()
	dout.addHeader('ml','Elution volume (mL)','f')
	dout.addHeader(target_column,'{name} ({unit})'.format(name=resfile[target_column]['data_name'], unit=resfile[target_column]['unit']),'f')
	dout.describeHeader(data_outs)

	dout.writeHeader(data_outs)
	format = dout.getFormat(named=True)
	data_length = len(resfile[target_column]['data'])
	n_written = 0
	for i in range(data_length):
		datdict = {}
		for x in target_columns:
			entry = resfile[x]['data'][i]
			datdict['ml'] = entry[0]
			datdict[x] = entry[1]
		#print datdict
		line = format.format(**datdict)
		data_outs.write(line)
		n_written += 1

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.out_fname))
		outf.close()

