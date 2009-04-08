import sys, os, math, string
import biofile

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
		

def digest(seq, patterns, complete=True):
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
	pat = "K/R"
	frags = digest(seq, pat, True)
	assert ''.join(frags) == seq
	print "\ttest001 passed"

if __name__ == '__main__':
	fname = sys.argv[1]
	if fname == "__test__":
		print "Running tests..."
		test001()
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
	

