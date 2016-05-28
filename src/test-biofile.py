import sys, os, math, string
import biofile

if __name__=='__main__':
	(h,s) = biofile.readFASTA('test-biofile/test-biofile-001.fa')
	assert len(h) == 143
	cd = biofile.readFASTADict(os.path.expanduser('test-biofile/test-biofile-001.fa'))
	assert len(cd.keys()) == len(h)
