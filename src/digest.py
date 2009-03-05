import sys, os, math, string
import biofile

def digest(seq, patterns, complete=True):
	for pat in patterns:
		frags = seq.split(pat)
		print frags
		if seq.startswith(pat):
			pass

if __name__ == '__main__':
	fname = sys.argv[1]
	patterns = sys.argv[2].split("/")
	(headers, seqs) = biofile.readFASTA(os.path.expanduser(fname))
	for seq in seqs:
		digest(seq, patterns, True)
	

