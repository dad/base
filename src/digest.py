import sys, os, math, string
import biofile

def split(seq, pat):
	frags = []
	s = seq
	prev_pos = -1
	pos = s.find(pat, prev_pos+1)
	while pos > -1:
		print pat, prev_pos, pos
		frags.append(s[prev_pos+1:pos+1])
		prev_pos = pos
		pos = s.find(pat, prev_pos+1)
	if prev_pos > 0:
		frags.append(s[prev_pos+1:])
	if prev_pos == -1 and pos == -1:
		# No matches at all
		frags.append(seq)
	return frags
		

def digest(seq, patterns, complete=True):
	all_frags = []
	digest = [seq]
	if complete:
		for pat in patterns:
			new_digest = []
			for frag in digest:
				frags = split(frag, pat) #["%s%s" % (p,pat) for p in frag.split(pat)]
				# if seq.startswith(pat):
				#	frags[0] = "%s%s" % (frags[0], pat)
				if frag.endswith(pat):
					frags = frags[:-1]
				new_digest += frags
			digest = new_digest
	return digest

if __name__ == '__main__':
	fname = sys.argv[1]
	patterns = sys.argv[2].split("/")
	(headers, seqs) = biofile.readFASTA(os.path.expanduser(fname))
	for seq in seqs:
		frags = digest(seq, patterns, True)
		print "\n%s\n---" % seq
		for f in frags:
			print f
	

