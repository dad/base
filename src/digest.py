import sys, os, math, string, re
import biofile

_enzymes = {"trypsin":"[KR][^P]", "lysc":"[K][*]"}

def split(seq, pat):
	frags = []
	s = seq
	prev_pos = -1
	pos = s.find(pat, prev_pos+1)
	while pos > -1:
		frags.append(s[prev_pos+1:pos+1])
		prev_pos = pos
		pos = s.find(pat, prev_pos+1)
	if prev_pos > -1:
		frags.append(s[prev_pos+1:])
	if prev_pos == -1 and pos == -1:
		# No matches at all -- add whole sequence
		frags.append(seq)
	return frags

def getCleavageSites(seq, compiled_pattern):
	"""Get indices of sites 
	"""
	cleavage_sites = set()
	for mat in compiled_pattern.finditer(seq):
		cleavage_sites.add(mat.start()+1)
		sub_seq = seq[mat.start()+1:mat.end()+1]
		#print mat.start()+1, mat.end()+1, seq, sub_seq
		sub_sites = getCleavageSites(sub_seq, compiled_pattern)
		# Move sub_sites forward by the placement of sub_seq in seq
		moved_sub_sites = set([s+mat.start()+1 for s in sub_sites])
		cleavage_sites = cleavage_sites.union(moved_sub_sites)
	return sorted(cleavage_sites)
		
def digest(seq, pattern, num_missed=0):
	pat = re.compile(pattern)
	frags = set()
	if num_missed == 0:
		frags = list()
	cleavage_sites = list(getCleavageSites(seq, pat))
	#print ','.join([str(x) for x in cleavage_sites])
	if len(cleavage_sites) > 0:
		# Now, generate set of peptides
		for miss in range(0,num_missed+1):
			if miss < len(cleavage_sites):
				# Loop over pairs of beginning/end
				beg_ind = -1
				end_ind = beg_ind + miss + 1
				beg = 0
				end = cleavage_sites[end_ind]
				done = False
				while not done: #beg < len(seq) and end <= len(seq):
					#print beg_ind, end_ind, beg, end, len(cleavage_sites), len(seq)
					frag = seq[beg:end]
					if num_missed == 0:
						frags.append(frag)
					else:
						frags.add(frag)
					beg_ind += 1
					end_ind += 1
					if end < len(seq):
						beg = cleavage_sites[beg_ind]
						if end_ind < len(cleavage_sites):
							end = cleavage_sites[end_ind]
						else:
							end = len(seq)
					else:
						done = True
		if num_missed == 0:
			# Need to remove missed cleavages due to pattern matching,
			# e.g. "KT" matches, and needs to be broken down to "K" and "T"
			frags = [f for f in frags if pat.match(f) is None]
	else:
		# No cleavages
		frags = [seq]
	return list(frags)

def digestWithEnzyme(seq, enzyme, num_missed=0):
	res = None
	try:
		pat = _enzymes[enzyme]
		res = digest(seq, pat, num_missed)
	except KeyError, ke:
		raise Exception, "Enzyme {0} not found".format(enzyme)
	return res

def digestOld(seq, patterns, complete=True):
	all_frags = []
	digest = [seq]
	if complete:
		for pat in patterns:
			new_digest = []
			for frag in digest:
				frags = split(frag, pat)
				new_digest += frags
			digest = new_digest
	return digest

def test001():
	seq = "MPIMLEDYQKNFLELAIECQALRFGSFKLKSGRESPYFFNLGLFNTGKLLSNLATAYAIAIIQSDLKFDVIFGPAYKGIPLAAIVCVKLAEIGGSKFQNIQYAFNRKEAKDHGEGGIIVGSALENKRILIIDDVMTAGTAINEAFEIISNAKGQVVGSIIALDRQEVVSTDDKEGLSATQTVSKKYGIPVLSIVSLIHIITYLEGRITAEEKSKIEQYLQTYGASA"
	seq = "MPIMLEDYQKNFLELAIECQALRFGSFKLKSGRESPYFFNLGLFNTGKLLSNLATAYAIAIIQSDLKFDVIFGPAYK"
	frags = digestWithEnzyme(seq, "trypsin", 0)
	#print ''.join(frags)
	#print seq
	assert ''.join(frags) == seq
	print "\ttest001 passed"

def test002():
	seq = "MPIMLEDYQKNFLELAIECQALRFGSFKLKSGRESPYFFNLGLFNTGKLLSNLATAYAIAIIQSDLKFDVIFGPAYK"
	frags = digestWithEnzyme(seq, "trypsin", 1)
	double_count = 0
	for frag in frags:
		k_count = frag.count("K")
		r_count = frag.count("R")
		both_count = k_count + r_count
		if both_count > 1:
			double_count += 1
		assert both_count < 3
		assert both_count >= 1
	assert double_count > 0
	print "\ttest002 passed"

def test003():
	seq = "KPAARPLLKPLLRPPK"
	frags = digestWithEnzyme(seq, "trypsin", 0)
	assert frags[0] == seq
	assert len(frags) == 1
	print "\ttest003 passed"

def test004():
	seq = "MQFSTVASVAFVALANFVAAESAAAVSQITDGQIQATTTATTEATTTAAPSSTVETVSPSSTETISQQTENGAAKAAVGMGAGALAAAAMLL"
	frags = digestWithEnzyme(seq, "trypsin", 1)
	assert len(frags)>0
	print "\ttest004 passed"

def test005():
	seq = "VRKT"
	frags = digestWithEnzyme(seq, "trypsin", 0)
	assert set(frags) == set(["VR","K","T"])
	print "\ttest005 passed"

def test006():
	seq = "VRKT"
	frags = digestWithEnzyme(seq, "trypsin", 1)
	assert set(frags) == set(["VR","VRK","KT","K","T"])
	print "\ttest006 passed"

if __name__ == '__main__':
	fname = sys.argv[1]
	if fname == "__test__":
		print "Running tests..."
		test001()
		test002()
		test003()
		test004()
		test005()
		test006()
		print "All tests passed"
		sys.exit()
	patterns = sys.argv[2].split("/")
	complete = True
	if os.path.isfile(os.path.expanduser(fname)):
		(headers, seqs) = biofile.readFASTA(os.path.expanduser(fname))
	else:
		seqs = [fname]
	for seq in seqs:
		frags = digest(seq, patterns, complete)
		#print "\n%s\n---" % seq
		for f in frags:
			print len(f), f
	if complete:
		assert ''.join(frags) == seq
	

