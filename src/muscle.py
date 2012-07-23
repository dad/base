#! python

import sys, os, random, string
from optparse import OptionParser
import biofile, util, translate, stats

class MuscleError(Exception):
	"""MUSCLE alignment error"""

def alignSequences(seq_list, max_iters=16, exepath="~/develop/muscle3.8.31/muscle"):
	tmp_fasta_file = "tmp-muscle-in-{}.txt".format(''.join(random.sample(string.ascii_letters, 20)))
	tmpfile = file(tmp_fasta_file, 'w')
	# Write out the sequences
	for si in range(len(seq_list)):
		tmpfile.write('>seq{:d}\n{}\n'.format(si, seq_list[si]))
	tmpfile.close()

	outfile_name = os.path.join(os.getcwd(),"tmp-muscle-out-{}.txt".format(''.join(random.sample(string.ascii_letters, 20))))

	cmd = "muscle -in {} -out {} -quiet -maxiters {:d}".format(tmp_fasta_file, outfile_name, max_iters)
	#print cmd
	error = os.spawnv(os.P_WAIT, os.path.expanduser(exepath), [x for x in cmd.split()])

	if not error:
		seq_dict = biofile.readFASTADict(outfile_name)
		seqs = [seq_dict["seq%d" % i] for i in range(len(seq_list))]
		os.remove(outfile_name)
		os.remove(tmp_fasta_file)
		return seqs
	else:
		if not os.path.isfile(os.path.expanduser(exepath)):
			raise MuscleError, "Can't find muscle executable at {}".format(os.path.expanduser(exepath))
		else:
			raise MuscleError, "Muscle error code %d" % error

def alignGeneFromProtein(gene, prot_align):
	j = 0
	gene_align = []
	for i in range(len(prot_align)):
		if prot_align[i] == '-':
			gene_align.append('---')
		else:
			gene_align.append(gene[j : 3 + j])
			j += 3
	return ''.join(gene_align)

def alignProteinFromProtein(prot, prot_align):
	j = 0
	out_align = ""
	for i in range(len(prot_align)):
		if prot_align[i] == '-':
			out_align += "-"
		else:
			out_align += prot[j]
			j += 1
	return out_align

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Muscle alignment")
	parser.add_argument("in_fname", help="input filename")
	parser.add_argument("-t", "--translate", dest="translate", action="store_true", default=False, help="translate the input sequences?")
	parser.add_argument("-o", "--out", dest="out_fname", type="string", default=None, help="output filename")
	options = parser.parse_args()
	
	outs = util.OutStreams()
	if not options.out_fname is None:
		outf = file(os.path.expanduser(options.out_fname),'w')
		outs.addStream(outf)
	else:
		outs.addStream(sys.stdout)
	
	(headers, seqs) = biofile.readFASTA(file(options.in_fname,'r'))
	seqs_to_align = seqs
	if options.translate:
		seqs_to_align = [translate.translate(s) for s in seqs]
	alseqs = alignSequences(seqs_to_align)
	if options.translate:
		alseqs = [alignGeneFromProtein(g, s) for (g,s) in zip(seqs,alseqs)]
	for (h,s) in zip(headers,alseqs):
		outs.write(">{}\n{}\n".format(h,s))
	
	if not options.out_fname is None:
		outf.close()

	
	
