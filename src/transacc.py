#! /usr/local/bin/python

import sys, os, math, string, random
sys.path.append(os.path.expanduser("~/research/lib/src"))
import translate, cai, stats

def sequenceDiffs(s1, s2):
	diffs = 0
	for i in range(len(s1)):
		if s1[i] != s2[i]:
			diffs += 1
	return diffs

def isTransition(nt1, nt2):
	l = [nt1,nt2]
	l.sort()
	return ''.join(l) in ["AG","CT","CU"]

# According to Freeland and Hurst 1998
def probMistranslation(codon1, codon2):
	if codon1==codon2:
		# This ain't a mistranslation!
		return 0.0
	elif sequenceDiffs(codon1, codon2)>1:
		# Double or triple changes have 0 probability
		return 0.0
	nt1 = None
	nt2 = None
	for i in range(3):
		if codon1[i] != codon2[i]:
			nt1 = codon1[i]
			nt2 = codon2[i]
			break
	# at site i, change from nt1 to nt2
	# weight transversions double, as they're twice as likely
	total_prob = 1+2*0.5+0.5+2*0.1+1+2*1
	trans = isTransition(nt1, nt2)
	if i == 0: # first position
		if trans:
			return 1.0/total_prob
		else:
			return 0.5/total_prob
	elif i == 1: # second position
		if trans:
			return 0.5/total_prob
		else:
			return 0.1/total_prob
	elif i == 2: # third position
		return 1.0/total_prob
	return 0.0

def isTolerated(aa, i, alignment, scores, conservation_cutoff):
	aas = [x[i] for x in alignment]
	#return aa in aas #aas.count(aa)>=0.1*len(aas)
	return scores[i]>conservation_cutoff

def probToleratedVariability(aa, i, alignment, scores, conservation_cutoff):
	aas = [x[i] for x in alignment]
	return float(scores[i]>conservation_cutoff)

def probToleratedAlignment(aa, i, alignment, scores, conservation_cutoff):
	aas = [x[i] for x in alignment]
	return float(aa in aas)

def pf(fs, sep='\t'):
	str = ["%1.8f"]*len(fs)
	return sep.join(str) % tuple(fs)

def siteErrorWeight(codon, rel_adapt, accuracy_factor):
	# Want to have the probability difference between each
	# to be dependent on a factor of our specification
	# Site weights have a mean of 1.
	w = 1.0
	if rel_adapt.has_key(codon):
		ra = rel_adapt[codon]
		if ra == 1.0:
			# Optimal codon has the minimum error probability
			w = 1.0/accuracy_factor
		else:
			w = 1.0 - ra
		#pf([ra,w])
	return w

def getSiteErrorProbabilities(gene, error_rate):
	# Probability of an error at each site
	codons = cai.splitByFrame(gene,0)
	site_error_wts = [siteErrorWeight(c, cai._yeast_relative_adaptiveness, 10.0) for c in codons]
	# Site weights sum to the length of the gene.
	weight_sum = sum(site_error_wts)
	site_error_wts = [len(codons)*w/weight_sum for w in site_error_wts]
	assert abs(sum(site_error_wts) - float(len(codons))) < 1e-5
	if False:
		for i in range(len(site_error_wts)):
			if prot[i] != bad_aa:
				print "%d\t%s\t%1.4f" % (i,prot[i],site_error_wts[i])
	#sew_red = [site_error_wts[i] for i in range(len(site_error_wts)) if prot[i] != bad_aa]
	sum_sews = sum(site_error_wts)
	site_error_probs = [w*error_rate for w in site_error_wts]
	return site_error_probs

def sign(x):
	if x>0:
		return 1
	elif x<0:
		return -1
	return 0

def hillClimb(gene, direction, max_steps, alignment, regions, tolerance_fxn, scores, conservation_cutoff, error_rate, tracer):
	# While
	codons = cai.splitByFrame(gene,0)
	step = 0
	total_steps = 0
	gc = translate._genetic_code
	all_alt_codons = [c for c in gc.keys() if not 'U' in c]
	prot = translate.TranslateRaw(gene,bad_aa='X')
	letters = 'ACDEFGHIKLNPQRSTVY'
	site_error_probs = getSiteErrorProbabilities(''.join(codons), error_rate)
	(orig_score, prob_acc, prob_no_error) = getTranslationOutcomes(''.join(codons), alignment, tolerance_fxn, scores, conservation_cutoff, all_alt_codons, site_error_probs)
	score = orig_score
	tracer.write("score\tprev.score\tdiff\torig.score\tstep\ttotal.step\n")
	tracer.write("# Starting hillclimb\n")
	while step < max_steps:
		# Pick synonymous codons at random
		# Pick AA at random
		aa = random.choice(letters)
		# Pick part of protein to randomize
		region = random.choice(regions)
		# Get indices of all synonymous codons with that amino acid
		indices = [i for i in range(region[0],region[1]) if prot[i]==aa]
		if len(indices)<2:
			continue
		inds = random.sample(indices, 2)
		# Swap them
		tmp = codons[inds[0]]
		codons[inds[0]] = codons[inds[1]]
		codons[inds[1]] = tmp
		# Score the resulting gene
		site_error_probs = getSiteErrorProbabilities(''.join(codons), error_rate)
		(new_score, prob_acc, prob_no_error) = getTranslationOutcomes(''.join(codons), alignment, tolerance_fxn, scores, conservation_cutoff, all_alt_codons, site_error_probs)
		# If score goes in the right direction, set step=0 and continue
		diff = new_score - score
		eps = 1e-5
		if abs(diff)>eps and sign(diff) == direction and checkSequence(''.join(codons), gene, prot):
			line = "#%s\n%s\t%d\t%d\n" %(''.join(codons), pf([new_score, score, diff, orig_score]), step, total_steps)
			print line,
			tracer.write(line)
			tracer.flush()
			step = 0
			score = new_score
		else:
			# Otherwise reverse swap
			tmp = codons[inds[0]]
			codons[inds[0]] = codons[inds[1]]
			codons[inds[1]] = tmp
			step += 1
		total_steps += 1
	return ''.join(codons)

def getTranslationOutcomes(gene, alignment, tolerance_fxn, scores, conservation_cutoff, all_alt_codons, site_error_probs):
	gc = translate._genetic_code
	bad_aa = 'X'
	prot = translate.TranslateRaw(gene, bad_aa)
	# Translate all possible codon point mutants of the protein
	# Compute the probability of folding given an arbitrary error
	prob_fold = 1.0
	prob_acc = 1.0
	prob_no_error = 1.0
	n_alternatives = 0
	codons = cai.splitByFrame(gene,0)
	prob_error_folds_list = []
	prob_syn_list = []
	if False:
		relads = [cai._yeast_relative_adaptiveness[c] for c in codons]
		print "# corr:", stats.Spearman_Rank_Correlation(relads, site_error_probs)
		print "# sd:", stats.StatsSummary(site_error_probs), min(site_error_probs), max(site_error_probs)
		print "# sd:", stats.StatsSummary(site_error_probs), min(site_error_probs), max(site_error_probs)
		for i in range(len(site_error_probs)):
			print site_error_probs[i], site_error_probs[i]/min(site_error_probs)
	#pf(sew_red)
	#diffs = [s/max(sew_red) for s in sew_red]
	#pf(diffs)
	#print min(sew_red), max(sew_red), stats.StatsSummary(diffs)

	for i in range(len(prot)):
		if prot[i] == bad_aa:
			continue
		codon = codons[i] #gene[3*i:3*i+3]
		#print codons
		#print ""
		site_prob = 0.0
		alt_codons = [c for c in all_alt_codons if sequenceDiffs(codon,c)==1]
		p_error_folds = 0.0 # Probability of folding given an error at this site.
		p_synonymous = 0.0 # Probability of a synonymous error
		p_missense_error_folds = 0.0 # Probability of folding given a missense error
		for ac in alt_codons:
			prob_mis_error = probMistranslation(codon, ac)
			site_prob += prob_mis_error
			aa = gc[ac]
			if aa == prot[i]:
				p_synonymous += prob_mis_error
				tol = 1.0
			else:
				#tol = False
				tol = tolerance_fxn(aa, i, alignment, scores, conservation_cutoff)
				p_missense_error_folds += prob_mis_error * tol
			prob_tolerated = tol #int(tol)
			p_error_folds += prob_mis_error * prob_tolerated
		# Probability of folding is the product of
		# the probability of this particular error, given an error at this site,
		# the probability of an error at this site,
		# and the probability that this error is tolerated.
		p_missense_error_folds = p_missense_error_folds/(1-p_synonymous)
		prob_error_folds_list.append(p_missense_error_folds)
		prob_syn_list.append(p_synonymous)
		prob_fold *= (1-site_error_probs[i]*(1-p_error_folds))
		prob_acc *= (1-site_error_probs[i]*(1-p_synonymous))
		prob_no_error *= (1-site_error_probs[i])

	return prob_fold, prob_acc, prob_no_error
