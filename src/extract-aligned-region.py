#! python

import sys, os
from argparse import ArgumentParser
import util, biofile, translate, geneutil

if __name__=='__main__':
	parser = ArgumentParser() #usage="%prog [-i fasta] [-s sequence]")
	parser.add_argument(dest="in_fname", default=None, help="input FASTA filename")

	# Optional arguments
	parser.add_argument("-o", "--out", dest="out_fname", default=None, help="output filename")
	parser.add_argument("-s", "--seq", dest="sequence", default=None, help="input sequence")
	parser.add_argument("--start-sequence", dest="start_sequence", type=str, default=None, help="sequence that starts the domain")
	parser.add_argument("--end-sequence", dest="end_sequence", type=str, default=None, help="sequence that ends the domain")
	parser.add_argument("--start-position", dest="start_position", type=int, default=None, help="starting sequence position for search (1-based)")
	parser.add_argument("--end-position", dest="end_position", type=int, default=None, help="ending sequence position for search (1-based, inclusive)")
	parser.add_argument("-x", "--exclude", dest="exclude", action="store_true", default=False, help="exclude rather than include start/end region?")
	parser.add_argument("-q", "--query", dest="query", type=str, default=None, help="specific sequence identifier to query")
	parser.add_argument("--fasta-out", dest="fasta_out_fname", default=None, help="output FASTA filename")
	#parser.add_argument("-r", "--report", dest="report", action="store_true", default=False, help="write long report per protein?")
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

	# Check constraints
	# If user provides start/end sequences, then must provide a query sequence
	if not options.start_sequence is None or not options.end_sequence is None:
		assert not options.query is None, "Must provide a query sequence identifier (--query) to use --start-sequence or --end-sequence."

	# Read input
	if not os.path.isfile(options.in_fname):
		raise IOError("# Error: file {} does not exist".format(options.in_fname))
	with open(options.in_fname,'r') as inf:
		# Read a FASTA file?
		(headers, seqs) = biofile.readFASTA(inf)
	info_outs.write("# Read {:d} sequences\n".format(len(seqs)))

	# Find query sequence(s)
	query_ids = []
	for (xi, h) in enumerate(headers):
		if options.query in h:
				query_ids.append(xi)
	if len(query_ids) == 0:
		info_outs.write("# Could not find sequences '{}'; exiting\n".format(options.query))
		sys.exit()
	if len(query_ids) > 1:
		info_outs.write("# Found more than one sequence matching '{}'; using the first one: \n#\t{}\n".format(options.query, headers[xi]))
	# Pick the first one
	query_id = query_ids[0]
	
	gap = '-'
	# Find subsequence
	(start_position,end_position) = (options.start_position, options.end_position)
	if not options.start_sequence is None:
		assert not options.end_sequence is None
		start_index = geneutil.gappedFind(seqs[query_id], options.start_sequence, gapless=False, gap=gap)
		end_index = geneutil.gappedFind(seqs[query_id], options.end_sequence, start=False, gapless=False, gap=gap)
		print(start_index, end_index)
		print(seqs[query_id][start_index:end_index])
	elif not options.start_position is None:
		assert not options.end_position is None
		start_index = geneutil.gappedIndex(seqs[query_id], options.start_position, gap=gap)
		end_index = geneutil.gappedIndex(seqs[query_id], options.end_position, gap=gap)
	else: # No start/end given; 
		info_outs.write("# No starting position or sequence given; nothing to do. Exiting\n")

	new_headers = []
	new_seqs = []
	for (h,seq) in zip(headers,seqs):
		if not options.exclude:
			ex_seq = seq[start_index:end_index]
		else: # Exclude the sequence
			#assert options.end_aa < len(seq)
			#assert options.begin_aa < options.end_aa
			ex_seq = seq[0:start_index] + seq[end_index:]
		#degapped_seq = seq.replace(gap,"")
		new_seqs.append(ex_seq)
		new_headers.append(h)
	seqs = new_seqs
	headers = new_headers

	# Write output
	biofile.writeFASTA(seqs, fasta_outs, headers=headers)
	n_written = len(seqs)

	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.fasta_out_fname is None:
		info_outs.write("# Wrote {} entries to {}\n".format(n_written, options.fasta_out_fname))
		outf.close()

	
