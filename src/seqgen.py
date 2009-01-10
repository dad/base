import sys, os, string, math, random
import stats, consensus, folder, translate

aas = 'ACDEFGHIKLMNPQRSTVWXY-'

def makeIndependentSiteSampledSequences(alignment, n_prots):
	seqs = []
	for i in range(n_prots):
		p = [random.choice([x[j] for x in alignment]) for j in range(len(alignment[0]))]
		seqs.append(''.join(p))
	return seqs

def makeConsensusSequences(alignment, n_prots, top_n):
	seqs = []
	probs = consensus.build_site_abs_probabilities(alignment, weights = None, letters=aas)
	pdicts = [dict([(aa,p) for (p,aa) in pl]) for pl in probs]
	cons_aas = []
	for j in range(len(alignment[0])):
		p_list = probs[j]
		p_list.sort(reverse=True)
		cons_aa_j = [aa for (prob, aa) in p_list[0:top_n]]  # Select the top N aa's by probability
		cons_aas.append(cons_aa_j)

	for i in range(n_prots):
		p = [random.choice(cons_aas[j]) for j in range(len(alignment[0]))]
		seqs.append(''.join(p))
	return seqs

def countMissedPairs(seq, pair_dict, pairs_to_check=None):
	missed_pairs = []
	if pairs_to_check is None:
		pairs_to_check = pair_dict.keys()
	for (i,j) in pairs_to_check:
		if not pair_dict[(i,j)].has_key((seq[i],seq[j])):
			missed_pairs.append((i,j))
	return missed_pairs

def pickByCumulativeProbability(probs, p):
	i = 0
	(prob, val) = probs[i]
	while prob < p:
		i += 1
		(prob, val) = probs[i]
	return val

def test_pickByCumulativeProbability():
	for n in range(100): # 100 tests
		raw_freqs = [random.random() for i in range(50)]
		vals = random.sample(xrange(10000),len(cum_probs))
		dist = makeCumulativeProbabilities(zip(raw_freqs, vals), normalize=True)
		small_diff = 1e-8
		for i in range(len(cum_probs)):
			p = cum_probs[i]
			v = pickByCumulativeProbability(dist, p-small_diff)
			vi = vals.index(v)
			assert vi == i
	print "# test_pickByCumulativeProbability: All tests passed"

def makeCumulativeProbabilitiesOld(probs, normalize=True):
	if normalize:
		sum_vals = float(sum(probs))
		probs = [x/sum_vals for x in probs]
	cum_probs = []
	running_sum = 0.0
	for p in probs:
		running_sum += p
		cum_probs.append(running_sum)
	return cum_probs

def makeCumulativeProbabilities(prob_val_pairs, normalize=True):
	prob_val_pairs.sort(reverse=True)
	probs = [p for (p,v) in prob_val_pairs]
	if normalize:
		sum_vals = float(sum(probs))
		probs = [x/sum_vals for x in probs]
	cum_probs = []
	running_sum = 0.0
	for p in probs:
		running_sum += p
		cum_probs.append(running_sum)
	cum_list = zip(cum_probs, [v for (p,v) in prob_val_pairs])
	return cum_list

def makePairProbabilityTable(alignment, weights=None):
	# At each pair site, make a list of probabilities
	pair_prob_dict = {}
	if weights is None:
		weights = [1.0/len(alignment) for x in alignment]
	n_seqs = len(alignment)
	L = len(alignment[0])
	for i in range(L-1):
		for j in range(i+1,L):
			pair_dict_ij = {}
			for k in range(len(alignment)):
				seq = alignment[k]
				# Add weights for each sequence as each pair is encountered.
				(aai,aaj) = (seq[i],seq[j])
				try:
					pair_dict_ij[(aai,aaj)] += weights[k]
				except KeyError:
					pair_dict_ij[(aai,aaj)] = weights[k]
			pair_prob_dict[(i,j)] = [(pair, freq/n_seqs) for (pair,freq) in pair_dict_ij.items()]
	pair_dict = {}
	for (i,j) in pair_prob_dict.keys():
		pair_dict[(i,j)] = dict(pair_prob_dict[(i,j)])
	return pair_dict

def makeProbabilityTable(alignment):
	# At each site, make a list of probabilities
	prob_dict = {}
	n_seqs = len(alignment)
	for i in range(len(alignment[0])-1):
		prob_dict_i = {}
		for seq in alignment:
			# Add weights for each sequence as each aa is encountered.
			aai = seq[i]
			try:
				prob_dict_i[aai] += 1.0
			except KeyError:
				prob_dict_i[aai] = 1.0
			prob_dict[i] = dict([(aa, freq/n_seqs) for (aa,freq) in prob_dict_i.items()])
	return prob_dict

def makePairCumulativeProbabilityTable(alignment, normalize=True):
	# At each pair site, make a sorted list of cumulative probabilities
	pair_cum_prob_dict = {}
	for i in range(len(alignment[0])-1):
		for j in range(i+1,len(alignment[0])):
			pair_dict_ij = {}
			for seq in alignment:
				# Add weights for each sequence as each pair is encountered.
				(aai,aaj) = (seq[i],seq[j])
				try:
					pair_dict_ij[(aai,aaj)] += 1.0
				except KeyError:
					pair_dict_ij[(aai,aaj)] = 1.0
			freqs_ij = [(c, pair) for (pair, c) in pair_dict_ij.items()]
			freqs_ij.sort(reverse=True)
			pair_cum_prob_dict[(i,j)] = makeCumulativeProbabilities(freqs_ij, normalize)
	# Implicitly, we will only care about pairs (i,j) with an entry in pair_cum_prob_dict.
	pair_dict = {}
	for (i,j) in pair_cum_prob_dict.keys():
		pair_dict[(i,j)] = dict([(pair,prob) for (prob,pair) in pair_cum_prob_dict[(i,j)]])
	return pair_dict

def makePairwiseSequences(alignment, n_prots, pair_cum_prob_dict):
	seqs = []
	tries = []
	# Implicitly, we will only care about pairs (i,j) with an entry in pair_cum_prob_dict.
	pair_dict = {}
	for (i,j) in pair_cum_prob_dict.keys():
		pair_dict[(i,j)] = dict([(pair,prob) for (prob,pair) in pair_cum_prob_dict[(i,j)]])

	L = len(alignment[0])
	max_tries = 2.0 * (L*(L-1.0)/2.0) * math.log((L*(L-1.0)/2.0))

	n_total_tries = 0
	while len(seqs)<n_prots:
		seq = makeIndependentSiteSampledSequences(alignment, 1)[0]
		missed_pairs = countMissedPairs(seq, pair_dict)
		n_missed_pairs = len(missed_pairs)
		n_tries = 0
		n_sites = len(seq)
		while n_missed_pairs > 0 and n_tries < max_tries:
			# Prioritize sites i, j by the number of missed pairs involving i or j?
			# Right now, just pick at random
			(i,j) = random.choice(missed_pairs)
			conflicts_resolved = False
			# Only try a fixed number of times.  If we get stuck, best to move to another conflicted pair and let this one relax.
			n_ij_tries = 0
			max_ij_tries = 2*(L-1) # The maximum number of disrupted pairs
			while not conflicts_resolved and n_ij_tries < max_ij_tries:
				# Get set of pairs involving (i,j)
				involved_pairs = [(x,i) for x in range(i) if x != i] + [(i,x) for x in range(i+1,n_sites) if x != i] + \
					[(x,j) for x in range(j) if x != j] + [(j,x) for x in range(j+1,n_sites) if x != j]
				missed_involved_pairs = list(set(missed_pairs).intersection(set(involved_pairs)))
				#print pair_dict.keys()
				if len(missed_involved_pairs) == 0:
					conflicts_resolved = True # No conflicts!
					continue
				# Sample from alignment to get new (i,j) amino acids
				rand = random.random()
				pair_probs = pair_cum_prob_dict[(i,j)]
				(aai, aaj) = pickByCumulativeProbability(pair_probs, rand)
				new_seq = [x for x in seq]
				new_seq[i] = aai
				new_seq[j] = aaj
				# Check to see if the number of missed pairs is reduced by this choice
				# If so, make the change, otherwise reject
				new_missed_pairs = countMissedPairs(''.join(new_seq), pair_dict, missed_involved_pairs)
				#print "# %d, %d -> failed %d, %d" % (i,j, len(new_missed_pairs), len(missed_involved_pairs))
				if len(new_missed_pairs) < len(missed_involved_pairs):
					seq = ''.join(new_seq)
					diff_pairs = list(set(missed_involved_pairs).difference(set(new_missed_pairs)))
					for pair in diff_pairs:
						missed_pairs.remove(pair)
					#print i,j, len(new_missed_pairs), len(missed_involved_pairs)
				if len(new_missed_pairs) == 0:
					conflicts_resolved = True
				n_ij_tries += 1
			# Increment total number of tries
			n_tries += n_ij_tries
			missed_pairs = countMissedPairs(seq, pair_dict)
			new_n_missed_pairs = len(missed_pairs)
			if new_n_missed_pairs < n_missed_pairs:
				n_missed_pairs = new_n_missed_pairs
				print new_n_missed_pairs
			n_tries += 1
		if not seq in alignment and not seq in seqs:
			print "# Generated sequence: %s" % seq
			seqs.append(seq)
			tries.append(n_tries)
			if len(seqs) != len(set(seqs)):
				print "# Huh? len(seqs) != len(set(seqs)), %d != %d, %s" % (len(seqs), len(set(seqs)), seq)
			#print len(seqs), n_total_tries + n_tries
			sys.stdout.flush()
			n_total_tries = 0
		else:
			n_total_tries += n_tries
	return seqs

def relaxPairwiseSequences(orig_seq, alignment, n_max_steps, pair_cum_prob_dict):
	pair_dict = {}
	for k in pair_cum_prob_dict.keys():
		pair_dict[k] = dict([(pair,prob) for (prob,pair) in pair_cum_prob_dict[k]])

	L = len(alignment[0])
	max_tries = (L*(L-1.0)/2.0) * math.log((L*(L-1.0)/2.0))

	seq = orig_seq
	missed_pairs = countMissedPairs(seq, pair_dict)
	n_missed_pairs = len(missed_pairs)
	if n_missed_pairs > 0:
		print "# Error: relaxPairwiseSequences called with %d missed pairs already" % n_missed_pairs
		return orig_seq
	n_steps = 0
	n_tries = 0
	n_sites = len(seq)
	sites = range(n_sites)
	pairs = pair_dict.keys()
	while n_steps < n_max_steps:
		# Prioritize sites i, j by the number of missed pairs involving i or j?
		# Right now, just pick at random
		(i,j) = random.choice(pairs)
		# Get set of pairs involving (i,j)
		involved_pairs = [(x,i) for x in range(i) if x != i] + [(i,x) for x in range(i+1,n_sites) if x != i] + \
			[(x,j) for x in range(j) if x != j] + [(j,x) for x in range(j+1,n_sites) if x != j]
		pairs_to_check = list(set(involved_pairs).intersection(set(pairs)))
		# Sample from alignment to get new (i,j) amino acids
		rand = random.random()
		pair_probs = pair_cum_prob_dict[(i,j)]
		(aai, aaj) = pickByCumulativeProbability(pair_probs, rand)
		new_seq = [x for x in seq]
		new_seq[i] = aai
		new_seq[j] = aaj
		# Check to see if no missed pairs are introduced by this choice
		# If so, make the change, otherwise reject
		new_missed_pairs = countMissedPairs(''.join(new_seq), pair_dict, pairs_to_check)
		#print "# %d, %d -> failed %d, %d" % (i,j, len(new_missed_pairs), len(missed_involved_pairs))
		if len(new_missed_pairs) == 0:
			seq = ''.join(new_seq)
			n_steps += 1
			n_tries += 1
		else:
			n_tries += 1
	print n_steps, n_tries
	return seq

def makePairwiseSiteSampledSequencesSimple(alignment, n_prots):
	# Does not include any of the optimizations in the makePairwiseSiteSampledSequences version above
	pair_dict = {}
	for i in range(len(alignment[0])-1):
		for j in range(i+1,len(alignment[0])):
			col_ij = [(x[i],x[j]) for x in alignment]
			pair_dict[(i,j)] = dict([(pair, 0) for pair in col_ij])
			for (aai,aaj) in col_ij:
				pair_dict[(i,j)][(aai,aaj)] += 1
	seqs = []
	while len(seqs) < n_prots:
		seq = makeIndependentSiteSampledSequences(alignment, 1)[0]
		missed_pairs = countMissedPairs(seq, pair_dict)
		best_score = len(missed_pairs)
		n_tries = 0
		while best_score > 0:
			(i,j) = random.choice(missed_pairs)
			n_seq = random.randint(0,len(alignment)-1)
			new_seq = [x for x in seq]
			new_seq[i] = alignment[n_seq][i]
			new_seq[j] = alignment[n_seq][j]
			new_seq = ''.join(new_seq)
			missed_pairs = countMissedPairs(new_seq, pair_dict)
			new_score = len(missed_pairs)
			if new_score <= best_score:
				seq = new_seq
				best_score = new_score
			n_tries += 1
		if not seq in alignment and not seq in seqs:
			seqs.append(''.join(seq))
		seqs.append(''.join(seq))
	return seqs

def makeCoupledSequences(alignment, n_prots, pairs_to_check):
	pair_dict = {}
	for i in range(len(alignment[0])-1):
		for j in range(i+1,len(alignment[0])):
			col_ij = [(x[i],x[j]) for x in alignment]
			pair_dict[(i,j)] = dict([(pair, 0) for pair in col_ij])
			for (aai,aaj) in col_ij:
				pair_dict[(i,j)][(aai,aaj)] += 1
	seqs = []
	while len(seqs) < n_prots:
		seq = makeIndependentSiteSampledSequences(alignment, 1)[0]
		missed_pairs = countMissedPairs(seq, pair_dict, pairs_to_check)
		best_score = len(missed_pairs)
		n_tries = 0
		while best_score > 0:
			(i,j) = random.choice(missed_pairs)
			n_seq = random.randint(0,len(alignment)-1)
			new_seq = [x for x in seq]
			new_seq[i] = alignment[n_seq][i]
			new_seq[j] = alignment[n_seq][j]
			new_seq = ''.join(new_seq)
			missed_pairs = countMissedPairs(new_seq, pair_dict, pairs_to_check)
			new_score = len(missed_pairs)
			if new_score <= best_score:
				seq = new_seq
				best_score = new_score
			n_tries += 1
		if not seq in alignment and not seq in seqs:
			seqs.append(''.join(seq))
	return seqs

def test():
	test_pickByCumulativeProbability()
	return
