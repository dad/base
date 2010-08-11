import sys, os, math, string, random
import muscle, translate, util

def mut(x, pmut, pgap):
	res = x
	if random.random() < pmut:
		res = random.choice(translate.AAs())
	if random.random() < pgap:
		res = res + '-'
	return res

def test001():
	s1 = ''.join(util.sample_wr(translate.AAs(), 100))
	others = [''.join([mut(x, 0.2,0.1) for x in s1]) for i in range(9)]
	seqs = [s1]+others
	als = muscle.alignSequences(seqs)
	res = len(als) == len(seqs)
	for i in range(len(als)):
		res = res and (als[i].replace("-",'') == seqs[i].replace("-",''))
	return True

def test002():
	s1 = ''.join(util.sample_wr(translate.AAs(), 100))
	others = [''.join([mut(x, 0.2,0.1) for x in s1]) for i in range(9)]
	seqs = [s1]+others
	res = False
	try:
		als = muscle.alignSequences(seqs, exepath="~/develop/muscle3.8.13/muscle")
	except muscle.MuscleError, me:
		res = True
	return res

if __name__=='__main__':
	res = True
	res1 = test001()
	res = res and res1
	res2 = test002()
	res = res and res2
	if res:
		print "# All tests passed."
