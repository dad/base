import sys, os, math, string, random, argparse
import stats, util, biofile, mq, na

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Extraction of evidence from MaxQuant evidence files")
	parser.add_argument(dest="database_fname", help="filename of FASTA database containing proteins to extract")
	parser.add_argument(dest="target_orf", help="target protein (ORF) to extract")
	parser.add_argument("-o", "--out", dest="out_fname", default=None, help="output filename")
	parser.add_argument("-e", "--evidence", dest="evidence_fnames", action="append", default=[], help="additional evidence files to merge")
	parser.add_argument("-x", "--experiment", dest="experiments", action="append", default=None, help="experiments to assay")
	parser.add_argument("-g", "--debug", dest="debugging", action="store_true", default=False, help="debug mode?")
	parser.add_argument("-m", "--merge", dest="merge", action="store_true", default=False, help="merge the indicated experiments?")
	parser.add_argument("-t", "--tag", dest="tags", action="append", default=[], help="tags to restrict the analysis to specific tagged experiments")
	parser.add_argument("-u", "--unique", dest="unique_matches", action="store_true", default=False, help="use unique peptides only?")
	options = parser.parse_args()

	# Set up some output
	info_outs = util.OutStreams(sys.stdout)

	orf_dict = None
	if not options.database_fname is None:
		orf_dict = biofile.readFASTADict(options.database_fname)

	# Pull out target protein
	target_prot = orf_dict[options.target_orf]
	if target_prot[-1] == '*':
		target_prot = target_prot[0:-1]

	evidences = []
	for fi in range(len(options.evidence_fnames)):
		fname = options.evidence_fnames[fi]
		if options.experiments is None:
			# If no experiments are specified, we assume invert refers to whole evidence files.
			ed = mq.EvidenceDescriptor()
			ed.filename = os.path.expanduser(fname)
			ed.tags = options.tags
			evidences.append(ed)
		else:
			for xi in range(len(options.experiments)):
				ed = mq.EvidenceDescriptor()
				ed.filename = os.path.expanduser(fname)
				ed.experiment = options.experiments[xi]
				ed.tags = options.tags
				evidences.append(ed)

	# Now, all experiments are set up.  Read in the data, then aggregate it.
	# First, read in.
	eef = mq.ExperimentEvidenceFactory()
	experiments = eef.load(evidences, options.tags, options.experiments, options.unique_matches, orf_dict)

	for ex in experiments:
		info_outs.write("# Read {0!s}\n".format(ex))

	# Get set of all peptide and protein identifiers.
	shared_peptide_keys = set([pep.key for pep in experiments[0].peptides])
	all_peptide_keys = set(list(shared_peptide_keys))
	all_protein_keys = set([prot.key for prot in experiments[0].proteins])
	if len(experiments) > 1:
		for ex in experiments[1:]:
			shared_peptide_keys = shared_peptide_keys.intersection(set([pep.key for pep in ex.peptides]))
			all_peptide_keys = all_peptide_keys.union(set([pep.key for pep in ex.peptides]))
			all_protein_keys = all_protein_keys.union(set([prot.key for prot in ex.proteins]))

	# Select subset of data that appears in all experiments
	# Get mean and variance of log-transformed H/L ratios for that subset
	# Normalize the data relative to that subset
	# DAD: disabled for now -- awaiting implementation of peptide-spike mixing normalization.
	'''
	for ex in experiments:
		ex.normalizePeptideRatios(...)
	'''
	if options.debugging:
		print "Before merging..."
		for ex in experiments:
			print ex
			pep = ex.peptide_data.items()[0]
			print pep[0], pep[1].heavy_light_ratio_list

	if not (options.merge or len(experiments) == 1):
		info_outs.write("# Multiple experiments without merging has not been implemented.  Did you mean to use --merge?\n")
		sys.exit()

	# Merge data
	merged_ex = experiments[0]
	if len(experiments) > 1:
		info_outs.write("# Merging {0} experiments into one.\n".format(len(experiments)))
		for ex in experiments[1:]:
			merged_ex = merged_ex.merge(ex, options.normalize_intensity)

	
	outs = util.OutStreams()
	if not options.out_fname is None:
		outf = file(os.path.expanduser(options.out_fname),'w')
		outs.addStream(outf)
	else:
		outs.addStream(sys.stdout)

	outs.write(">{}\n{}\n".format(options.target_orf, target_prot))

	prot_rec = merged_ex.getProtein(options.target_orf)
	pep_list = []
	for pep in prot_rec.peptides:
		pos = target_prot.find(pep.sequence)
		if pos>-1:
			pep_list.append((pos, pep.sequence))
	pep_list.sort()
	gap = '-'
	len_prot = len(target_prot)
	n_peps = 0
	for (pos, pep) in pep_list:
		n_peps += 1
		pepid = "{}-{}".format(options.target_orf, n_peps)
		line = gap*pos + pep + gap*(len_prot-(len(pep)+pos))
		outs.write(">{}\n{}\n".format(pepid, line))
			
	