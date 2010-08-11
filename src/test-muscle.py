import sys, os, math, string, random
import muscle, translate, util

def mut(x, p):
	res = x
	if random.random() < p:
		res = random.choice(translate.AAs())
	return res

def test001():
	s1 = ''.join(util.sample_wr(translate.AAs(), 100))
	others = [''.join([mut(x, 0.2) for x in s1]) for i in range(9)]
	seqs = [s1]+others
	als = muscle.alignSequences(seqs)
	assert len(als) == len(seqs)
	for i in range(len(als)):
		assert als[i].replace("-",'') == seqs[i]
	return True

if __name__=='__main__':
	res = True
	res1 = test001()
	res = res and res1
	if res:
		print "# All tests passed."
