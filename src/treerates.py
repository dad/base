#! python

import sys, os, math, string, random, pickle
sys.path = [os.path.expanduser('~/research/lib/')] + sys.path
import stats, translate, my_paml, phylip, orthologs, muscle


# Set up amino-acid degeneracy table
codon_degeneracy = {}
for codon in [c for c in translate._genetic_code.keys() if not 'U' in c]:
	codon_degeneracy[codon] = 0
	aa = translate._genetic_code[codon]
	for base in 'ATGC':
		new_codon = codon[0:2]+base
		if translate._genetic_code[new_codon] == aa:
			codon_degeneracy[codon] += 1
codon_degeneracy['---'] = 0

#for (k,v) in codon_degeneracy.items():
#	print k,v

def frac_aligned(seqid, numid, numal, length):
	return numal/float(length)

def sequence_identity(seqid, numid, numal, length):
	return seqid

def align_stats(seqs, queryfxn, statsfxn):
	vals = []
	for i in range(len(seqs)-1):
		for j in range(i+1,len(seqs)):
			(seqid, numid, numal) = orthologs.sequenceIdentity(seqs[i], seqs[j])
			vals.append(queryfxn(seqid, numid, numal, len(seqs[i])))
	return statsfxn(vals)

def get_tree_species(tree):
	tree = tree.replace('(',' ')
	tree = tree.replace(',',' ')
	tree = tree.replace(')',' ')
	return tree.split()

def get_tree_rates(ortho_dict, alignment_dict, cdna_dicts, tree_string, begin_index, end_index):
	tree = phylip.parse_tree(tree_string)
	tree_species = [n.name for n in tree.leaves()]
	
	n_genes = 0

	# Individual genes
	rate_ancestor_cache = {}

	# Assemble gene list
	keys = []
	for gene in ortho_dict.keys():
		try:
			(length, corr_keys, aligned_prots) = alignment_dict[gene]
			alignment_species = [spec for (spec,orf) in corr_keys]
		except KeyError, ke:
			print "# Couldn't find alignments for", ke
			continue

		# Must have concordance between tree species and alignment species
		if set(tree_species).intersection(set(alignment_species)) == set(tree_species):
			keys.append(gene)
			
	keys.sort()
	print "# Found", len(keys), "genes"
	for gene in keys[begin_index:min(end_index,len(keys))]:
		try:
			(length, corr_keys, aligned_prots) = alignment_dict[gene]
			alignment_species = [spec for (spec,orf) in corr_keys]
		except KeyError, ke:
			print "# Couldn't find alignments for", ke
			continue

		if set(tree_species).intersection(set(alignment_species)) != set(tree_species):
			continue

		msid = align_stats(aligned_prots, sequence_identity, stats.Median)
		mfal = align_stats(aligned_prots, frac_aligned, min)
		#msid = min_seq_id(aligned_prots)
		#mfal = min_frac_aligned(aligned_prots)
		#print gene, mfal
		#continue
		if mfal < 0.5:
			print "# rejected orf (mfal,msid) %s (%1.2f,%1.2f)" % (gene, mfal, msid)
			continue

		all_prots = [xprot for ((xgenome,xorf), xprot) in zip(corr_keys,aligned_prots) if xgenome in tree_species]
		all_genes = []
		for (xgenome, xorf) in corr_keys:
			if xgenome in tree_species:
				try:
					seq = cdna_dicts[xgenome][xorf]
				except KeyError:
					seq = cdna_dicts[xgenome+'-mit'][xorf]
				all_genes.append(seq)
		all_species = tree_species #[xgenome for (xgenome, xorf) in corr_keys]
		assert len(all_genes) == len(all_prots)
		assert len(all_species) == len(all_genes)
		all_aligned_genes = [muscle.align_gene_from_protein(xgene,xprot) for (xgene,xprot) in zip(all_genes,all_prots)]

		# Put the genes in the right order
		sub_gene_dict = dict(zip(all_species, all_aligned_genes))
		recon_gene_list = [sub_gene_dict[spec] for spec in tree_species]

		try:
			#(rates,ancestor) = (1,2)
			(dns, dss) = my_paml.Get_dNdS_Per_Codon(recon_gene_list, tree_species, tree, 100)
			#(rates, ancestor) = my_paml.Get_Site_Rates_And_Ancestor(recon_gene_list, tree_species, tree)
			(r, p, n) = stats.PearsonCorrelation(dns, dss)
			print "# %d %s %1.3f" % (n_genes, gene, r)
		except my_paml.PAMLError, pe:
			print "#",pe
			continue
		#print len(rates), len(ancestor)
		rate_ancestor_cache[gene] = (dns, dss, tree_string)#(gene, rates, ancestor, tree)
		n_genes += 1
		#if n_genes > 1:
		#	break
	return rate_ancestor_cache

def get_rate_correlation_windowed(dns, dss, aligned_prots, cdna_dicts, spec_orf_list, xfold_only, xfold_degeneracy):
	window_size = 1
	ns = []
	ss = []
	# Examine only codons of a certain degeneracy?
	xfold_ending = 'CTAG'
	xfold_wrong_aas = ''
	aligned_cdnas = []
	for xi in range(len(spec_orf_list)):
		(spec,orf) = spec_orf_list[xi]
		try:
			aligned_cdnas.append(muscle.align_gene_from_protein(cdna_dicts[spec][orf], aligned_prots[xi]))
		except KeyError:
			aligned_cdnas.append(muscle.align_gene_from_protein(cdna_dicts[spec+'-mit'][orf], aligned_prots[xi]))
		assert len(aligned_prots[xi]) == len(dns)		

	for site in range(0,len(dns),window_size):
		if xfold_only and (window_size == 1):
			codons = [aligned_cdna[3*site:3*site+3] for aligned_cdna in aligned_cdnas]
			wrong_degeneracy = False
			wrong_ending = False
			wrong_aa = False
			degs = [codon_degeneracy[codon] for codon in codons]
			for codon in codons:
				if codon_degeneracy[codon] != xfold_degeneracy:
					wrong_degeneracy = True
				if not codon[2] in xfold_ending:
					wrong_ending = True
				if not codon=='---' and translate._genetic_code[codon] in xfold_wrong_aas:
					wrong_aa = True
			if wrong_degeneracy or wrong_ending or wrong_aa:
				continue
		# Add up substitutions in window
		syn = 0
		nsyn = 0
		valid_sites = False
		for nsite in range(site, min(len(dns),site+window_size)):
			(s,n) = (dss[nsite], dns[nsite])
			# don't consider cases with missing counts
			if not (s is None) and not (n is None):
				syn += s
				nsyn += n
				valid_sites = True
		if valid_sites:
			ns.append(nsyn)
			ss.append(syn)
	return stats.PearsonCorrelation(ns,ss), ns, ss

def write_tree_rates(rate_dict, alignment_dict, cdna_dicts, out_filename, xfold_only, xfold_degeneracy):
	outfile = file(out_filename,'w')
	outfile.write("orf\tnal\tcor\tpcor\tncor\tns\tss\tmsid\tmfal\tadn\n")
	ngenes = 0
	for orf in rate_dict.keys():
		(length, spec_orf_list, aligned_prots) = alignment_dict[orf]
		(dns, dss, tree) = rate_dict[orf]
		root = phylip.parse_tree(tree)
		species = root.leaves()
		msid = align_stats(aligned_prots, sequence_identity, stats.Mean)
		mfal = align_stats(aligned_prots, frac_aligned, min)
		#msid = min_seq_id(aligned_prots)
		#mfal = min_frac_aligned(aligned_prots)
		if mfal < 0.5:
			print "# rejected orf (mfal,msid) %s (%1.2f,%1.2f)" % (orf, mfal, msid)			
			continue
		num_genes = len(species)

		#   Compute correlation between nonsyn-syn changes
		try:
			#target_orf = dict(spec_orf_list)[spec_name]
			((r,p,n), ns, ss) = get_rate_correlation_windowed(dns, dss, aligned_prots, cdna_dicts, spec_orf_list, xfold_only, xfold_degeneracy)
			outfile.write("%s\t%d\t%1.3f\t%1.3f\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.4f\n" % \
				  (orf,num_genes,r,p,n,sum(ns),sum(ss),msid,mfal,sum(ns)/float(len(dns))))
			ngenes += 1
			outfile.flush()
		except stats.StatsError, se:
			continue
	outfile.close()
	print "# Wrote", ngenes, "tree rates to", out_filename

def read_tree(fname):
	f = file(fname,'r')
	for line in f.readlines():
		if line[0] != "#":
			return line.strip()
	return ""

def check_alignments_and_rates(alignment_fname, tree_rate_fname):
	rate_cache = pickle.load(file(tree_rate_fname, 'r'))
	alignment_cache = pickle.load(file(alignment_fname, 'r'))
	for orf in alignment_cache.keys():
		(nal, hdrs, aligns) = alignment_cache[orf]
		print ""
		for na in range(nal):
			print "%10s\n%s" % (hdrs[na], aligns[na])

def convert_yeast_alignment(alignment_fname):
	[alignment_dict, corr_dict] = pickle.load(file(alignment_fname,'r'))
	alignment_dict_all = {}
	for orf in alignment_dict.keys():
		prots = alignment_dict[orf]
		hdrs = [(xorg,xorf) for (xorf,xorg) in corr_dict[orf]]
		length = len(prots)
		alignment_dict_all[orf] = (length, hdrs, prots)
	pickle.dump(alignment_dict_all, file("cerev-ortholog-alignments.p",'w'))

