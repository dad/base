import sys, os, math, string, random, argparse
import stats, util, biofile, mq, na

# Example:
# extract-triple-evidence.py -e evidence.txt -x experimentname --out output.txt
# extract-triple-evidence.py -e evidence.txt -x experimentname --unique --out unique-output.txt
# extract-triple-evidence.py -e evidence.txt -x experiment1 -x experiment2 --merge --out merged-output.txt

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Extraction of evidence from MaxQuant evidence files")
	parser.add_argument("-i", "--in", dest="in_fname", default=None, help="input filename")
	parser.add_argument("-o", "--out", dest="out_fname", default=None, help="output filename")
	parser.add_argument("--peptide-out", dest="peptide_out_fname", default=None, help="peptide output filename")
	#parser.add_argument("--protein-out", dest="protein_out_fname", default=None, help="protein output filename")
	#parser.add_argument("-v", "--invert", dest="invert", default="", \
	#				  help="string (Y,T,1 all meaning 'true') indicating, for each experiment, whether to invert heavy/light labels, e.g. " + \
	#				  "011 or NYY or FTT indicates that heavy/light ratios in experiment 1 are to be merged with light/heavy ratios in experiments 2 and 3.")
	parser.add_argument("-e", "--evidence", dest="evidence_fnames", action="append", default=[], help="additional evidence files to merge")
	parser.add_argument("-x", "--experiment", dest="experiments", action="append", default=None, help="experiments to assay")
	parser.add_argument("-d", "--database", dest="database_fname", default=None, help="filename of FASTA database containing proteins to extract")
	parser.add_argument("-g", "--debug", dest="debugging", action="store_true", default=False, help="debug mode?")
	parser.add_argument("-m", "--merge", dest="merge", action="store_true", default=False, help="merge the indicated experiments?")
	parser.add_argument("-t", "--tag", dest="tags", action="append", default=[], help="tags to restrict the analysis to specific tagged experiments")
	parser.add_argument("-u", "--unique", dest="unique_matches", action="store_true", default=False, help="use unique peptides only?")
	parser.add_argument("--normalize-intensity", dest="normalize_intensity", action="store_true", help="normalize intensity when merging?")
	parser.add_argument("--normalize-ratio-by", dest="normalize_ratio_by_orf", default=None, help="ORF to use for normalization across runs")
	parser.add_argument("--ratio-sig", dest="ratio_significance_field", default="ratio_hl_normalized", help="field to use for ratio significance calculations")
	parser.add_argument("--abundance", dest="abundance_field", default="intensity", help="field to use for abundance calculations")
	options = parser.parse_args()

	# Set up some output
	info_outs = util.OutStreams(sys.stdout)

	orf_dict = None
	if not options.database_fname is None:
		orf_dict = biofile.readFASTADict(options.database_fname)

	evidences = []
	if not options.in_fname is None:
		#print "# Loading..."
		# Read more experiments from master file
		inf = file(os.path.expanduser(options.in_fname), 'r')
		dlr = util.DelimitedLineReader(inf, header=True)
		while not dlr.atEnd():
			flds = dlr.nextDict()
			if os.path.isfile(os.path.expanduser(flds['filename'])):
				ed = mq.EvidenceDescriptor()
				ed.filename = os.path.expanduser(flds['filename'])
				ed.invert = flds['invert'][0].lower() in ['1','y','t']
				ed.tags = [x.strip() for x in flds['tags'].split(',')]
				ed.experiment = flds['experiment']
				evidences.append(ed)
			else:
				info_outs.write("# Can't find file {0}\n".format(flds['filename']))
		inf.close()

	for fi in range(len(options.evidence_fnames)):
		fname = options.evidence_fnames[fi]
		if options.experiments is None:
			# If no experiments are specified, we assume invert refers to whole evidence files.
			ed = mq.EvidenceDescriptor()
			ed.filename = os.path.expanduser(fname)
			#if not na.isNA(options.invert):
			#	ed.invert = (options.invert.lower()[fi] in ['1','y','t'])
			#else:
			#	ed.invert = False
			ed.tags = options.tags
			evidences.append(ed)
		else:
			for xi in range(len(options.experiments)):
				ed = mq.EvidenceDescriptor()
				ed.filename = os.path.expanduser(fname)
				ed.experiment = options.experiments[xi]
				ed.invert = False
				#if len(options.invert) > 0:
				#	assert len(options.invert) == len(options.experiments)
				#	ed.invert = (options.invert.lower()[xi] in ['1','y','t'])
				ed.tags = options.tags
				evidences.append(ed)

	# Now, all experiments are set up.  Read in the data, then aggregate it.
	# First, read in.
	eef = mq.ExperimentEvidenceFactory()
	experiments = eef.load(evidences, options.tags, options.experiments, options.unique_matches, [], orf_dict)

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

	if options.normalize_intensity:
		info_outs.write("# Normalizing to {0} shared peptide intensities\n".format(len(shared_peptide_keys)))
		for ex in experiments:
			ex.normalizePeptideIntensities(shared_peptide_keys)

	if options.debugging:
		print "Before merging..."
		for ex in experiments:
			print ex
			pep = ex.peptide_data.items()[0]
			print pep[0], pep[1].heavy_light_ratio_list

	if not (options.merge or len(experiments) == 1):
		info_outs.write("# Multiple experiments without merging has not been implemented.  Did you mean to use --merge?\n")
		sys.exit()

	# Normalize
	if not options.normalize_ratio_by_orf is None:
		for ex in experiments:
			try:
				norm_prot = ex.protein_data[options.normalize_ratio_by_orf]
				ex.normalizeRatiosBy(norm_prot)
			except KeyError:
				info_outs.write("# Error: could not find normalization protein {0}\n".format(options.normalize_ratio_by_orf))

	# Merge data
	merged_ex = experiments[0]
	if len(experiments) > 1:
		info_outs.write("# Merging {0} experiments into one.\n".format(len(experiments)))
		for ex in experiments[1:]:
			merged_ex = merged_ex.merge(ex, options.normalize_intensity)

	if options.debugging:
		print "After merging..."
		for ex in experiments:
			print ex
			pep = ex.peptide_data.items()[0]
			print pep[0], pep[1].heavy_light_ratio_list

	# Calculate ratio significances
	# DAD: rewrite to incorporate intelligent model.
	#merged_ex.calculateProteinRatioSignificance(options.ratio_significance_window, options.ratio_significance_field, options.abundance_field)

	outs = util.OutStreams()
	if not options.out_fname is None:
		outf = file(os.path.expanduser(options.out_fname),'w')
		outs.addStream(outf)
	else:
		outs.addStream(sys.stdout)

	# Print out peptide data for the merged experiment?
	pep_outs = util.OutStreams()
	if not options.peptide_out_fname is None:
		pep_outf = file(options.peptide_out_fname,'w')
		pep_outs.addStream(pep_outf)

		pep_outs.write("# {0!s}\n".format(merged_ex))
		pep_header = "seq\tproteins\tn.proteins"
		for rat in ['hl','ml','hm']:
			pep_header += "\tratio.{0}\tratio.{0}.mean\tratio.{0}.normalized\tratio.{0}.normalized.mean\tratio.{0}.normalized.lower.95\tratio.{0}.normalized.upper.95\tratio.{0}.count\tratio.{0}.sd\tratio.{0}.normalized.sd".format(rat)
			pep_header += "\tiratio.{0}\tiratio.{0}.mean\tiratio.{0}.count\tiratio.{0}.sd".format(rat)
		pep_header += "\tintensity\tintensity.h\tintensity.m\tintensity.l\tms.ms.count\tn.fractions\tfractions\n"
		pep_outs.write(pep_header)
		n_written = 0
		pep_keys = sorted(list(all_peptide_keys))
		for p in pep_keys:
			pep = merged_ex.peptide_data[p]
			# Get proteins for this peptide
			prot_ids = pep.parent_proteins
			prot_ids = list(prot_ids)
			prot_ids.sort()
			if prot_ids == []:
				line = '{0}\tNA\t0'.format(pep.key)
			else:
				line = '{0}\t{1}\t{2}'.format(p, ",".join(prot_ids), len(prot_ids))

			output_fields = []
			for rat in ['hl','ml','hm']:
				ratio_stats = pep.getRatioSummary(rat)
				ratio_norm_stats = pep.getNormalizedRatioSummary(rat)
				output_fields.append(util.FieldFormatter(ratio_stats.median,"{0:e}"))
				output_fields.append(util.FieldFormatter(ratio_stats.mean,"{0:e}"))
				output_fields.append(util.FieldFormatter(ratio_norm_stats.median,"{0:e}"))
				output_fields.append(util.FieldFormatter(ratio_norm_stats.mean,"{0:e}"))
				rn_lower_95 = None
				rn_upper_95 = None
				if not na.isNA(ratio_norm_stats.se):
					rn_lower_95 = math.exp(math.log(ratio_norm_stats.mean)-1.96*ratio_norm_stats.se)
					rn_upper_95 = math.exp(math.log(ratio_norm_stats.mean)+1.96*ratio_norm_stats.se)
				output_fields.append(util.FieldFormatter(rn_lower_95,"{0:e}"))
				output_fields.append(util.FieldFormatter(rn_upper_95,"{0:e}"))
				output_fields.append(util.FieldFormatter(ratio_stats.n,"{0:d}"))
				output_fields.append(util.FieldFormatter(ratio_stats.sd,"{0:e}"))
				output_fields.append(util.FieldFormatter(ratio_norm_stats.sd,"{0:e}"))
				# Intensity ratios -- no "normalized" ratios here.
				iratio_stats = pep.getIntensityRatioSummary(rat)
				output_fields.append(util.FieldFormatter(iratio_stats.median,"{0:e}"))
				output_fields.append(util.FieldFormatter(iratio_stats.mean,"{0:e}"))
				output_fields.append(util.FieldFormatter(iratio_stats.n,"{0:d}"))
				output_fields.append(util.FieldFormatter(iratio_stats.sd,"{0:e}"))

			output_fields.append(util.FieldFormatter(pep.intensity,"{0:e}"))
			output_fields.append(util.FieldFormatter(pep.heavy_intensity,"{0:e}"))
			output_fields.append(util.FieldFormatter(pep.medium_intensity,"{0:e}"))
			output_fields.append(util.FieldFormatter(pep.light_intensity,"{0:e}"))
			output_fields.append(util.FieldFormatter(pep.msms_count,"{0:d}"))
			fractions = [str(s) for s in pep.fractions]
			#output_fields.append(util.FieldFormatter(len(fractions),"{0:d}"))
			output_fields.append(util.FieldFormatter(','.join(fractions),"{0}"))

			line += "\t" + "\t".join([str(f) for f in output_fields])
			pep_outs.write("{0}\n".format(line))
			n_written += 1
		info_outs.write("# Wrote {0} peptide records to {1}\n".format(n_written, options.peptide_out_fname))

	outs.write("# Merging {} experiments\n".format(len(experiments)))
	for ex in experiments:
		outs.write("# {0!s}\n".format(ex))
	outs.write("# Merged into:\n")
	outs.write("# {0!s}\n".format(merged_ex))
	header = "orf\tn.proteins\tn.peptides"
	for rat in ['hl','ml','hm']:
		header += "\tratio.{0}\tratio.{0}.mean\tratio.{0}.normalized\tratio.{0}.normalized.mean\tratio.{0}.normalized.lower.95\tratio.{0}.normalized.upper.95\tratio.{0}.count\tratio.{0}.sd\tratio.{0}.normalized.sd".format(rat)
		header += "\tiratio.{0}\tiratio.{0}.mean\tiratio.{0}.count\tiratio.{0}.sd".format(rat)
	header += "\tintensity\tintensity.h\tintensity.m\tintensity.l\tms.ms.count\n"
	outs.write(header)
	n_written = 0
	prot_list = [(x.id, x) for x in merged_ex.proteins]
	for (prot_id,prot) in sorted(prot_list):
		line = '{}'.format(prot_id)
		if options.debugging:
			for pep in prot.peptides:
				for ratio in pep.heavy_light_normalized_ratio_list:
					if not na.isNA(ratio):
						print ratio, type(ratio)
		output_fields = []
		output_fields.append(util.FieldFormatter(prot.degeneracy,"{0:d}"))
		output_fields.append(util.FieldFormatter(prot.n_peptides,"{0:d}"))
		for rat in ['hl','ml','hm']:
			ratio_stats = prot.getRatioSummary(rat)
			ratio_norm_stats = prot.getNormalizedRatioSummary(rat)
			output_fields.append(util.FieldFormatter(ratio_stats.median,"{0:e}"))
			output_fields.append(util.FieldFormatter(ratio_stats.mean,"{0:e}"))
			output_fields.append(util.FieldFormatter(ratio_norm_stats.median,"{0:e}"))
			output_fields.append(util.FieldFormatter(ratio_norm_stats.mean,"{0:e}"))
			rn_lower_95 = None
			rn_upper_95 = None
			if not ratio_norm_stats.se is None:
				rn_lower_95 = math.exp(math.log(ratio_norm_stats.mean)-1.96*ratio_norm_stats.se)
				rn_upper_95 = math.exp(math.log(ratio_norm_stats.mean)+1.96*ratio_norm_stats.se)
			output_fields.append(util.FieldFormatter(rn_lower_95,"{0:e}"))
			output_fields.append(util.FieldFormatter(rn_upper_95,"{0:e}"))
			output_fields.append(util.FieldFormatter(ratio_stats.n,"{0:d}"))
			output_fields.append(util.FieldFormatter(ratio_stats.sd,"{0:e}"))
			output_fields.append(util.FieldFormatter(ratio_norm_stats.sd,"{0:e}"))
			
			# Intensity ratios -- no "normalized" ratios here.
			iratio_stats = prot.getIntensityRatioSummary(rat)
			output_fields.append(util.FieldFormatter(iratio_stats.median,"{0:e}"))
			output_fields.append(util.FieldFormatter(iratio_stats.mean,"{0:e}"))
			output_fields.append(util.FieldFormatter(iratio_stats.n,"{0:d}"))
			output_fields.append(util.FieldFormatter(iratio_stats.sd,"{0:e}"))
		output_fields.append(util.FieldFormatter(prot.intensity,"{0:e}"))
		output_fields.append(util.FieldFormatter(prot.heavy_intensity,"{0:e}"))
		output_fields.append(util.FieldFormatter(prot.medium_intensity,"{0:e}"))
		output_fields.append(util.FieldFormatter(prot.light_intensity,"{0:e}"))
		output_fields.append(util.FieldFormatter(prot.msms_count,"{0:d}"))
		line += "\t" + "\t".join([str(f) for f in output_fields])
		outs.write("{}\n".format(line))
		n_written += 1
	if not options.out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.out_fname))
		outf.close()
