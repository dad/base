#! /usr/local/bin/python

import sys, os, math, string, random
import translate

class MuscleError(Exception):
	"""MUSCLE alignment error"""

def align_sequences_with_header(seq_list, max_iters=16, exepath="~/develop/muscle3.6_src/muscle"):
	tmp_fasta_file = "tmp%d.txt" % random.randint(0,100000)
	tmpfile = file(tmp_fasta_file, 'w')
	for (hdr, seq) in seq_list:
		tmpfile.write('>%s\n%s\n' % (hdr, seq))
	tmpfile.close()

	outfile_name = "tmp-muscle-out-%d.txt" % random.randint(0,100000)

	cmd = "muscle -in %s -out %s -quiet -stable -maxiters %d" % (tmp_fasta_file, outfile_name, max_iters)
	#print cmd
	error = os.spawnv(os.P_WAIT, os.path.expanduser(exepath), [x for x in cmd.split()])
	if not error:
		seq_dict = translate.Read_FASTA_Dict(outfile_name)
		os.remove(outfile_name)
		os.remove(tmp_fasta_file)
		return seq_dict.items()
	else:
		raise MuscleError, "Muscle error code %d" % error

def align_sequences(seq_list, max_iters=16, exepath="~/develop/muscle3.6_src/muscle"):
	return alignSequences(seq_list, max_iters, exepath)

def alignSequences(seq_list, max_iters=16, exepath="~/develop/muscle3.6_src/muscle"):
	tmp_fasta_file = "tmp%d.txt" % random.randint(0,100000)
	full_tmp_fasta_file = tmp_fasta_file #os.path.join(os.getcwd(),tmp_fasta_file)
	tmpfile = file(tmp_fasta_file, 'w')
	#tmpfile = file(full_tmp_fasta_file, 'w')
	i = 0
	for seq in seq_list:
		tmpfile.write('>%d\n%s\n' % (i, seq))
		i += 1
	tmpfile.close()

	outfile_name = os.path.join(os.getcwd(),"tmp-muscle-out-%d.txt" % random.randint(0,100000))

	cmd = "muscle -in %s -out %s -quiet -stable -maxiters %d" % (full_tmp_fasta_file, outfile_name, max_iters)
	#print cmd
	error = os.spawnv(os.P_WAIT, os.path.expanduser(exepath), [x for x in cmd.split()])

	if not error:
		(hdrs, seqs) = translate.Read_FASTA(outfile_name)
		os.remove(outfile_name)
		os.remove(full_tmp_fasta_file)
		return seqs
	else:
		raise MuscleError, "Muscle error code %d" % error

def alignProfiles(seq_file1, seq_file2, max_iters=16, exepath="~/develop/muscle3.6_src/muscle"):
	return align_profiles(seq_file1, seq_file2, max_iters, exepath)

def align_profiles(seq_file1, seq_file2, max_iters=16, exepath="~/develop/muscle3.6_src/muscle"):
	outfile_name = os.path.join(os.getcwd(),"tmp-muscle-out-%d.txt" % random.randint(0,100000))

	cmd = "muscle -in1 %s -in2 %s -out %s -quiet -maxiters %d -profile" % (seq_file1, seq_file2, outfile_name, max_iters)
	#print cmd
	error = os.spawnv(os.P_WAIT, os.path.expanduser(exepath), [x for x in cmd.split()])

	if not error:
		aldict = translate.Read_FASTA_Dict(outfile_name)
		os.remove(outfile_name)
		return aldict
	else:
		raise MuscleError, "Muscle error code %d" % error

def align_gene_from_protein(gene, prot_align):
	return alignGeneFromProtein(gene, prot_align)

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
