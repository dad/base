#! python

import sys, os, math, random, argparse
from pycorn import pc_res3
import util, biofile, stats

def hasExtension(fname, target_ext):
	return fname.split('.')[-1].strip() == target_ext

def replaceExtension(fname, new_ext='txt'):
	new_fname = '.'.join(fname.split('.')[:-1]) + '.txt'
	return new_fname

def processFile(in_fname, out_fname, options, data_outs, info_outs):
	resfile = pc_res3(in_fname)

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

	# Write out parameters
	data_outs.write("# Run started {}\n".format(util.timestamp()))
	data_outs.write("# Command: {}\n".format(' '.join(sys.argv)))
	data_outs.write("# Parameters:\n")
	optdict = vars(options)
	for (k,v) in optdict.items():
		data_outs.write("#\t{k}: {v}\n".format(k=k, v=v))
	data_outs.write("# Input filename: {}\n".format(in_fname))
	data_outs.write("# Output filename: {}\n".format(out_fname))

	# Write output
	dout = util.DelimitedOutput()
	dout.addHeader('ml','Elution volume (mL)','f')
	dout.addHeader(target_column,'{name} ({unit})'.format(name=resfile[target_column]['data_name'], unit=resfile[target_column]['unit']),'f')
	dout.describeHeader(data_outs)

	dout.writeHeader(data_outs)
	format = dout.getFormat(named=True)
	data_length = len(resfile[target_column]['data'])
	n_written = 0
	last_entry_volume = resfile[target_column]['data'][0][0]
	last_entry = False
	cur_entry_volume = None
	# Tuples of (mL,UV)
	elution_window = []
	for i in range(data_length):
		entry = resfile[target_column]['data'][i]
		cur_entry_volume = entry[0]
		cur_entry_value = entry[1]
		#print last_entry_volume, cur_entry_volume
		if cur_entry_volume-last_entry_volume>=options.resolution_ml or i==data_length-1:
			# Compute averages
			datdict = {
				'ml':stats.mean([x[0] for x in elution_window]),
				target_column:stats.mean([x[1] for x in elution_window])
			}
			# Write out
			line = format.format(**datdict)
			data_outs.write(line)
			n_written += 1
			# Reset
			last_entry_volume = cur_entry_volume
			elution_window = []
		elution_window.append((cur_entry_volume, cur_entry_value))

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, out_fname))
		outf.close()

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Extract tab-delimited data from AKTA .res files")
	# Required arguments
	parser.add_argument(dest="in_fname", default=None, help="input filename or directory")
	# Optional arguments
	parser.add_argument("--resolution-ml", dest="resolution_ml", default=0.005, help="output resolution, in mL: average output over elution-volume ranges no larger than this number")
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

	# Create the instance using default options
	if os.path.isdir(options.in_fname):
		basedir = options.in_fname
		for in_fname in os.listdir(options.in_fname):
			full_fname = os.path.join(basedir, in_fname)
			if hasExtension(full_fname, 'res'):
				out_fname = replaceExtension(full_fname, 'txt')
				print full_fname, out_fname
				data_outs = util.OutStreams(file(out_fname,'w'))
				processFile(full_fname, out_fname, options, data_outs, info_outs)
	else: # Not a directory
		processFile(options.in_fname, options.out_fname, options, data_outs, info_outs)

	sys.exit()

