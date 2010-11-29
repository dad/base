#! python

import sys, os, random, string
import biofile

class MuscleError(Exception):
	"""MUSCLE alignment error"""

def alignSequences(seq_list, max_iters=16, exepath="~/develop/muscle3.8.31/muscle"):
	tmp_fasta_file = "tmp-muscle-in-%s.txt" % ''.join(random.sample(string.ascii_letters, 20))
	tmpfile = file(tmp_fasta_file, 'w')
	# Write out the sequences
	for si in range(len(seq_list)):
		tmpfile.write('>seq%d\n%s\n' % (si, seq_list[si]))
	tmpfile.close()

	outfile_name = os.path.join(os.getcwd(),"tmp-muscle-out-%s.txt" % ''.join(random.sample(string.ascii_letters, 20)))

	cmd = "muscle -in %s -out %s -quiet -maxiters %d" % (tmp_fasta_file, outfile_name, max_iters)
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
			raise MuscleError, "Can't find muscle executable at %s" % os.path.expanduser(exepath)
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
	return string.join(gene_align,'')
