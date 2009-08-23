#! python

import sys, os, math, string, random, pickle
import translate, muscle, paml, biofile

def default_alignment_print_fxn(num_alignments, prots, alignment, headers, orf):
	print num_alignments, orf, len(alignment), " ".join(["%s-%s"%(x,y) for (x,y) in headers])

def default_filter_fxn(orf, seqs, filter_data=None):
	return len(seqs)>1

def make_alignments(ortho_dict, cdna_dicts, filter_fxn=default_filter_fxn, alignment_print_fxn=default_alignment_print_fxn):
	return makeAlignments(ortho_dict, cdna_dicts, filter_fxn, alignment_print_fxn)

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

def alignStats(seqs, queryfxn, statsfxn):
	vals = []
	for i in range(len(seqs)-1):
		for j in range(i+1,len(seqs)):
			(seqid, numid, numal) = translate.sequenceIdentity(seqs[i], seqs[j])
			vals.append(queryfxn(seqid, numid, numal, len(seqs[i])))
	return statsfxn(vals)

def getAlignments(alignment_dict, tree, alignment_threshold, outstreams):
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
	for outstream in outstreams:
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
	for outstream in outstreams:
		outstream.write(line)
	return species_orf_dict, alignment_map


# Compute distances between pairs of genes.
# Doesn't include calculation of kappa
def compute_pairwise_stats_nokappa(pair_dict, alignment_dict, gene_dicts, dist_fxn, begin_index, end_index, paml_options=[]):
	distance_cache = {}
	i = 0
	genes = pair_dict.keys()
	genes.sort()
	num_genes = len(genes)
	for gene_id in genes[begin_index:end_index]:
		[(mspec, morf), (relspec,relorf)] = pair_dict[gene_id]
		(num_aligned, hdrs, protal) = alignment_dict[gene_id]
		aldict = dict(zip([xspec for (xspec,xorf) in hdrs], protal))
		base_prot = aldict[mspec]
		base_gene = gene_dicts[mspec][morf]
		base_gene = muscle.align_gene_from_protein(base_gene, base_prot)
		# Special-casing...yikes
		try:
			query_prot = aldict[relspec]
		except KeyError:
			query_prot = aldict[relspec.split('-')[0]]
		query_gene = gene_dicts[relspec][relorf]
		query_gene = muscle.align_gene_from_protein(query_gene, query_prot)
		genes = [base_gene, query_gene]
		assert(len(query_gene)==3*len(query_prot) or len(query_gene)==3*(len(query_prot)+1) )
		try:
			(dNML, dSML, numNonsynonymousSites, numSynonymousSites) = dist_fxn(genes, options=paml_options)
		except paml.PAMLError, mpe:
			print "#", base_gene
			print "#", query_gene
			continue
		(seqIdentity, numIdentical, numAligned) = translate.sequenceIdentity(query_prot, base_prot)
		fracAligned = numAligned/float(len(query_prot.replace('-','')))
		# Correct for issues with PAML?
		correct_paml = False
		if correct_paml:
			if fracAligned == 1.0 and seqIdentity == 1.0:
				dNML = 0.0 # Obviously no nonsyn. changes if proteins are identical
			elif seqIdentity < 1.0 and 1.0/numNonsynonymousSites > dNML:
				dNML = 1.0/numNonsynonymousSites # Minimum possible change
			(ntseqIdentity, ntnumIdentical, ntnumAligned) = translate.sequenceIdentity(query_gene, base_gene)
			ntfracAligned = ntnumAligned/float(len(base_gene.replace('-','')))
			if ntfracAligned == 1.0 and ntseqIdentity == 1.0:
				dSML = 0.0 # Obviously no syn. changes if genes are identical
			elif ntseqIdentity < 1.0 and 1.0/numSynonymousSites > dSML:
				dSML = 1.0/numSynonymousSites # Minimum possible change

		# Cache this distance and alignment
		distance_cache[gene_id] = (dSML, dNML, numSynonymousSites, numNonsynonymousSites, fracAligned, seqIdentity)
		print "# %d distances of %d: %s-%s (%s-%s) %1.6f %1.6f %1.2f %1.2f" % (i, num_genes, mspec, relspec, morf, relorf, dNML, dSML, numSynonymousSites, numNonsynonymousSites)
		sys.stdout.flush()
		i += 1
	return distance_cache

# Compute distances between pairs of genes.
def computePairwiseStats(pair_dict, alignment_dict, gene_dicts, distFxn, begin_index, end_index, paml_options):
	distance_cache = {}
	i = 0
	genes = pair_dict.keys()
	genes.sort()
	num_genes = len(genes)
	#print begin_index, end_index
	for gene_id in genes[begin_index:end_index]:
		[(mspec, morf), (relspec,relorf)] = pair_dict[gene_id]
		(num_aligned, hdrs, protal) = alignment_dict[gene_id]
		aldict = dict(zip([xspec for (xspec,xorf) in hdrs], protal))
		try:
			base_prot = aldict[mspec]
			base_gene = gene_dicts[mspec][morf]
			base_gene = muscle.alignGeneFromProtein(base_gene, base_prot)
			query_prot = aldict[relspec]
			query_gene = gene_dicts[relspec][relorf]
			query_gene = muscle.alignGeneFromProtein(query_gene, query_prot)
			genes = [base_gene, query_gene]
		except KeyError, ke:
			print "# Gene or species not found: %s" % ke
			continue
		assert(len(query_gene)==3*len(query_prot) or len(query_gene)==3*(len(query_prot)+1) )
		try:
			(dNML, dSML, numNonsynonymousSites, numSynonymousSites, ts_tv_kappa) = distFxn(base_gene, query_gene, options=paml_options)
		except paml.PAMLError, mpe:
			print "# PAMLError %s" % mpe
			#print "#", base_gene
			#print "#", query_gene
			continue
		(seqIdentity, numIdentical, numAligned) = translate.sequenceIdentity(query_prot, base_prot)
		fracAligned = numAligned/float(len(query_prot.replace('-','')))
		# Correct for issues with PAML?
		correct_paml = False
		if correct_paml:
			if fracAligned == 1.0 and seqIdentity == 1.0:
				dNML = 0.0 # Obviously no nonsyn. changes if proteins are identical
			elif seqIdentity < 1.0 and 1.0/numNonsynonymousSites > dNML:
				dNML = 1.0/numNonsynonymousSites # Minimum possible change
			(ntseqIdentity, ntnumIdentical, ntnumAligned) = translate.sequenceIdentity(query_gene, base_gene)
			ntfracAligned = ntnumAligned/float(len(base_gene.replace('-','')))
			if ntfracAligned == 1.0 and ntseqIdentity == 1.0:
				dSML = 0.0 # Obviously no syn. changes if genes are identical
			elif ntseqIdentity < 1.0 and 1.0/numSynonymousSites > dSML:
				dSML = 1.0/numSynonymousSites # Minimum possible change
		# Cache this distance and alignment
		distance_cache[gene_id] = (dSML, dNML, numSynonymousSites, numNonsynonymousSites, fracAligned, seqIdentity, ts_tv_kappa)
		print "# %d distances of %d: %s-%s (%s-%s) %1.6f %1.6f %1.2f %1.2f %1.2f" % (i, num_genes, mspec, relspec, morf, relorf, dNML, dSML, numSynonymousSites, numNonsynonymousSites, ts_tv_kappa)
		sys.stdout.flush()
		i += 1
	return distance_cache

def compute_pairwise_stats(pair_dict, alignment_dict, gene_dicts, dist_fxn, begin_index, end_index, paml_options=[]):
	return computePairwiseStats(pair_dict, alignment_dict, gene_dicts, dist_fxn, begin_index, end_index, paml_options)


# Compute distances between genes in a tree.
def computeStats(ortho_dict, alignment_dict, gene_dicts, tree, dist_fxn, begin_index, end_index, paml_options=[]):
	distance_cache = {}
	n_dist = 0
	genes = list(set(ortho_dict.keys()).intersection(alignment_dict.keys()))
	genes.sort()
	num_genes = len(genes)
	tree_species = [n.name for n in tree.leaves()]
	print "# Tree species:", tree_species
	for gene_id in genes[begin_index:end_index]:
		# Align the genes
		(num_aligned, hdrs, protal) = alignment_dict[gene_id]
		if num_aligned < len(tree_species):
			continue
		aldict = dict(zip([xspec for (xspec,xorf) in hdrs], protal))
		ordict = dict(ortho_dict[gene_id])
		genes = []
		orfs = []
		prots = []
		for spec in tree_species:
			orf = ordict[spec]
			prot = aldict[spec]
			seq = gene_dicts[spec][orf]
			al_seq = muscle.align_gene_from_protein(seq, prot)
			genes.append(al_seq)
			orfs.append(orf)
			prots.append(prot)
		try:
			(dNML, dSML, numNonsynonymousSites, numSynonymousSites, ts_tv_kappa) = \
				dist_fxn(genes, seq_labels=tree_species, tree_string='%s'%tree, options=paml_options)
		except paml.PAMLError, mpe:
			for (o, g) in zip(orfs,genes):
				print "#", o, g
			continue
		min_seq_id = 1.0
		min_frac_aligned = 1.0
		for i in range(0,len(prots)-1):
			for j in range(i+1,len(prots)):
				(seqIdentity, numIdentical, numAligned) = translate.sequenceIdentity(prots[i], prots[j])
				fracAligned = numAligned/float(len(prots[j]) - prots[j].count('-'))
				if seqIdentity < min_seq_id:
					min_seq_id = seqIdentity
				if fracAligned < min_frac_aligned:
					min_frac_aligned = fracAligned
		# Cache this distance and alignment
		distance_cache[gene_id] = (dSML, dNML, numSynonymousSites, numNonsynonymousSites, min_frac_aligned, min_seq_id, ts_tv_kappa)
		print "# %d distances of %d: %s (%s) %1.6f %1.6f %1.2f %1.2f %1.2f %1.2f %1.2f" % (n_dist, num_genes, '-'.join(tree_species), '-'.join(orfs), dNML, dSML, numSynonymousSites, numNonsynonymousSites, ts_tv_kappa, min_frac_aligned, min_seq_id)
		sys.stdout.flush()
		n_dist += 1
	return distance_cache

def read_genomes_from_file(multi_files_fname, genome_dir, genome_dicts, column_index=1, load_fxn=biofile.firstField):
	return readGenomesFromFile(multi_files_fname, genome_dir, genome_dicts, column_index, load_fxn)

def readGenomesFromFile(multi_files_fname, genome_dir, genome_dicts, column_index=1, load_fxn=biofile.firstField, species=None):
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
			print "# Cannot find file %s" % genome_file
		genome = biofile.readFASTADict(genome_file, load_fxn)
		genome_dicts[spec] = genome
		print "#", spec, genome_file, len(genome.keys()), genome.keys()[0]
	return species_map

