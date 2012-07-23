import sys, os, math, string, random, stats
import muscle, translate, status

def mut(x, pmut, pgap):
	res = x
	if random.random() < pmut:
		res = random.choice(translate.AAs())
	if random.random() < pgap:
		res = res + '-'
	return res

def test001():
	s1 = ''.join(stats.sample_wr(translate.AAs(), 100))
	others = [''.join([mut(x, 0.2,0.1) for x in s1]) for i in range(9)]
	seqs = [s1]+others
	als = muscle.alignSequences(seqs)
	res = len(als) == len(seqs)
	for i in range(len(als)):
		res = res and (als[i].replace("-",'') == seqs[i].replace("-",''))
	return True

def test002():
	s1 = ''.join(stats.sample_wr(translate.AAs(), 100))
	others = [''.join([mut(x, 0.2,0.1) for x in s1]) for i in range(9)]
	seqs = [s1]+others
	res = False
	try:
		als = muscle.alignSequences(seqs, exepath=os.path.expanduser("~/develop/muscle3.8.13/muscle"))
	except muscle.MuscleError, me:
		res = True
	return res

if __name__=='__main__':
	tests = [test001, test002]
	test_results = [t() for t in tests]
	all_passed = reduce(lambda x,y: x and y, test_results)
	if all_passed:
		print "# All tests passed."
	else:
		for i in range(len(tests)):
			if not test_results[i]:
				print "# Test {0} failed".format(tests[i])
