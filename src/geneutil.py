#! python

import sys, os, math, string, random, pickle
import translate, muscle, biofile, util
import re

def primerMinMelting(seq):
	"""Return melting temperature of primer"""
	degs = {'T':2.0, 'A':2.0, 'G':2.0, 'C':2.0}
	temp = sum([degs[x] for x in seq])
	return temp

def gappedFind(seq, substring, start=True, gap='-'):
	"""Find substring in seq, permitting gaps in seq.
		e.g.:
		geneutil.gappedFind('AAASS--SAA','SSS')==3
		geneutil.gappedFind('AAASS--SAA','SSS',start=False)==8
	"""
	gap_pattern = ('[{}]*'.format(gap)).join([a for a in substring])
	pat = re.compile(gap_pattern)
	res = re.search(pat, seq)
	ind = -1
	if not res is None:
		if start:
			ind = res.start()
		else:
			ind = res.end()
	return ind

def longestRun(seq, character_list, max_interruptions=0):
	"""Find the longest run of character in seq, permitting no more than max_interruptions.
		E.g. longestRun('AAAAA','A') = 5
		longestRun('AAATAA','A',1) = 6
		longestRun('AAATTAA','A',1) = 3
	"""
	#print seq
	# Trivial case where characters do not occur in sequence
	if isinstance(character_list, str):
		character_list = [x for x in character_list]

	# No characters
	if len(character_list)==0:
		return 0

	in_seq = True
	for character in character_list:
		in_seq = in_seq or (character in seq)
	if not in_seq:
		return 0
			
	# Trivial case where sequence is 1 characters long
	if len(seq)<2:
		return len(seq)
	chars = list(set([x for x in seq]))
	other_chars = chars[:]
	for character in character_list:
		if character in other_chars:
			other_chars.remove(character)
	# No other characters? Also trivial.
	if len(other_chars)==0:
		return len(seq)
	# Mask the sequence, making all non-character elements into mask_char
	mask_char = '@'
	assert not mask_char in chars
	masked_seq = seq
	for c in other_chars:
		masked_seq = masked_seq.replace(c,mask_char)
		#print masked_seq

	# Mask sequence 
	target_char = '~'
	if target_char in seq:
		raise ValueError, "# Current implementation does not work if {} in sequence".format(target_char)
	for c in character_list:
		masked_seq = masked_seq.replace(c,target_char)

	character = target_char

	longest_run = -1
	if max_interruptions == 0:
		# Easy case
		runs = masked_seq.split(mask_char)
		longest_run = max([len(r) for r in runs])
	elif max_interruptions>0:
		char_match_pat = re.compile('{}+'.format(character))
		char_matches = re.findall(char_match_pat, masked_seq)
		other_match_pat = re.compile('{}+'.format(mask_char))
		other_matches = re.findall(other_match_pat, masked_seq)
		if masked_seq[0] == mask_char:
			other_matches = other_matches[1:]
		if masked_seq[-1] == character:
			other_matches.append('')
		run_lens = [len(r) for r in char_matches]
		irun_lens = [len(r) for r in other_matches]
		for xi in range(len(run_lens)):
			run_length = run_lens[xi]
			irun_length = irun_lens[xi]
			last_irun_length = 0
			yi = xi+1
			while irun_length<=max_interruptions and yi<len(run_lens):
				run_length += run_lens[yi]
				last_irun_length = irun_length
				irun_length += irun_lens[yi]
				yi += 1
			# Include the last interrupting run length as well,
			# to yield the total tract length
			# E.g. AAATAAATTT = 7, not 6
			if run_length+last_irun_length>longest_run:
				longest_run = run_length+last_irun_length
	return longest_run

def maxSlidingCount(seq, character, windowlen=5):
    """Find the maximum numbers of character in sliding window of length windowlen in seq.
		E.g. maxSlidingCount('AAAAA','A') = 5
		maxSlidingCount('AAATAA','A') = 4
		maxSlidingCount('AAATTAA','A') = 3
	"""
    if windowlen <= len(seq):
        max_count = max([seq[i:(i+windowlen)].count(character) for i in range(len(seq)-windowlen+1)])
    else:
        max_count = seq.count(character)
    return max_count


def default_alignment_print_fxn(num_alignments, prots, alignment, headers, orf):
	print num_alignments, orf, len(alignment), " ".join(["%s-%s"%(x,y) for (x,y) in headers])

def default_filter_fxn(orf, seqs, filter_data=None):
	return len(seqs)>1

def makeAlignments(ortho_dict, cdna_dicts, filter_fxn=default_filter_fxn, filter_data=None, alignment_print_fxn=default_alignment_print_fxn):
	alignment_dict = {}
	num_aligns = 0
	#print cdna_dicts.keys()

	for orf in ortho_dict.keys():
		ortho_orfs = ortho_dict[orf]
		#print orf, ortho_orfs
		seqs = {}
		for (spec, sorf) in ortho_orfs:
			try:
				genome = cdna_dicts[spec]
				seq = genome[sorf]
				# Translate and so on
				prot = translate.translate(seq)
				if prot:
					seqs[spec] = (sorf, prot)
				else:
					print "# protein", sorf, "did not translate"
					#print seq
					#print translate.translateRaw(seq)
			except KeyError, ke:
				print "#", ke, spec, sorf, orf
				pass

		species = seqs.keys()
		if filter_fxn(orf, seqs, filter_data): #len(species) == len(genome_dicts.keys()): # Found as many orthologs as genomes
			prots = [seqs[key][1] for key in species]
			try:
				protal = muscle.alignSequences(prots, 16)
				hdrs = [(spec, seqs[spec][0]) for spec in species]
				alignment_dict[orf] = (len(protal), hdrs, protal)
				num_aligns += 1
				alignment_print_fxn(num_aligns, prots, protal, hdrs, orf)
			except muscle.MuscleError, me:
				print "#", me

	return alignment_dict

def fracAligned(seqid, numid, numal, length):
	return numal/float(length)

def multipleFractionAligned(seqs, gap="-"):
	non_gap_counts = [len(s)-s.count(gap) for s in seqs]
	all_ungapped_count = 0
	for i in range(len(seqs[0])):
		if not gap in [x[i] for x in seqs]:
			all_ungapped_count += 1
	return [all_ungapped_count/float(c) for c in non_gap_counts]


def alignStats(seqs, queryfxn, statsfxn):
	vals = []
	for i in range(len(seqs)-1):
		for j in range(i+1,len(seqs)):
			(seqid, numid, numal) = translate.sequenceIdentity(seqs[i], seqs[j])
			vals.append(queryfxn(seqid, numid, numal, len(seqs[i])))
	return statsfxn(vals)

def getAlignments(alignment_dict, tree, alignment_threshold, outstream):
	# Assemble set of acceptable alignments
	species_orf_dict = {}
	alignment_map = {}
	tree_species = [n.name for n in tree.leaves]
	tree_species_set = set(tree_species)
	for spec in tree_species:
		species_orf_dict[spec] = []
	bad_aligns = 0
	missing_species = 0
	num_usable = 0

	line = "# Species = (%s)\n# Found %d total ORFs with alignments\n" % (', '.join(tree_species), len(alignment_dict.keys()))
	outstream.write(line)
	for orf in alignment_dict.keys():
		(al_len, spec_orf_pairs, aligned_prots) = alignment_dict[orf]
		#print al_len, spec_orf_pairs
		spec_list = [xspec for (xspec,xorf) in spec_orf_pairs]
		specs = set(spec_list)
		# Check to make sure all species in tree exist
		if tree_species_set.intersection(specs) != tree_species_set:
			missing_species += 1
			continue
		# Now get the corresponding proteins
		aldict = dict(zip(spec_list, aligned_prots))
		tree_aligned_prots = [aldict[s] for s in tree_species]
		mfal = alignStats(tree_aligned_prots, fracAligned, min)
		if mfal < alignment_threshold:
			bad_aligns += 1
			continue
		spec_dict = dict(spec_orf_pairs)
		num_usable += 1
		for spec in tree_species:
			species_orf_dict[spec].append(spec_dict[spec])
			alignment_map[spec_dict[spec]] = orf

	line = "# Rejected %d ORFs with missing aligned species\n# Rejected %d ORFs with insufficient fraction aligned\n# Found %d usable alignments\n" % \
		   (missing_species, bad_aligns, num_usable)
	outstream.write(line)
	return species_orf_dict, alignment_map

def readGenomesFromFile(multi_files_fname, genome_dir, genome_dicts, column_index=1, load_fxn=biofile.firstField, species=None, outstream=None):
	if outstream is None:
		outstream = util.OutStreams()
	# Format for
	species_map = {}
	for line in file(multi_files_fname,'r').readlines():
		if line[0] != '#' and not line.strip() == '':  # skip comments and blank lines
			flds = line.strip().split()
			#print flds, column_index
			species_map[flds[0]] = flds[column_index]
	if species is None:
		species = species_map.keys()
	else:
		assert set(species).intersection(set(species_map.keys())) == set(species), "Not all specified species found in mapping file"

	for spec in species:
		genome_file = os.path.join(os.path.expanduser(genome_dir), species_map[spec])
		if not os.path.isfile(genome_file):
			outstream.write("# Cannot find file %s\n" % genome_file)
		genome = biofile.readFASTADict(genome_file, load_fxn)
		genome_dicts[spec] = genome
		outstream.write("# species=%s, genome file=%s has %d entries, example ID=%s\n" % (spec, genome_file, len(genome.keys()), genome.keys()[0]))
	return species_map

_aa_molecular_weight_dict = {
	'A': 71.079, 'C': 103.145, 'E': 129.116, 'D': 115.089,
	'G': 57.052, 'F': 147.177, 'I': 113.16, 'H': 137.141,
	'K': 128.17, 'M': 131.2, 'L': 113.16, 'N': 114.104,
	'Q': 128.131, 'P': 97.117, 'S': 87.078, 'R': 156.188,
	'T': 101.105, 'W': 186.213, 'V': 99.133, 'Y': 163.176, '*':0.0}

def getMolecularWeight(prot_seq):
	mw = sum([_aa_molecular_weight_dict[aa] for aa in prot_seq])
	return mw

def getAAFraction(seq, aa, pseudocount=0.0):
	return (seq.count(aa) + pseudocount)/len(seq)

if __name__=='__main__':
	cdna_dict = biofile.readFASTADict(sys.argv[1])
	keys = cdna_dict.keys()
	for i in range(len(keys)-1):
		for j in range(i+1, len(keys)):
			s1 = cdna_dict[keys[i]]
			s2 = cdna_dict[keys[j]]
			print keys[i], keys[j]
			print translate.sequenceIdentity(s1,s2)
