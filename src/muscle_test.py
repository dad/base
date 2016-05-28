import sys, os, math, string, random, unittest
import muscle, translate, stats

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
		#als = muscle.alignSequences(seqs, exepath=os.path.expanduser("/Users/dad/develop/muscle3.8.13/muscle"))
		als = muscle.alignSequences(seqs)
	except muscle.MuscleError as me:
		print(me)
		res = True
	return res

class test001(unittest.TestCase):
	def test_remove_gaps(self):
		s1 = ''.join(stats.sample_wr(translate.AAs(), 100))
		others = [''.join([mut(x, 0.2,0.1) for x in s1]) for i in range(9)]
		seqs = [s1]+others
		als = muscle.alignSequences(seqs)
		res = len(als) == len(seqs)
		for i in range(len(als)):
			self.assertTrue(als[i].replace("-",'') == seqs[i].replace("-",''))

class test002(unittest.TestCase):
	def test_gapped_index(self):
		s1 = ''.join(stats.sample_wr(translate.AAs(), 50))
		others = [''.join([mut(x, 0.2,0.1) for x in s1]) for i in range(9)]
		seqs = [s1]+others
		res = False
		try:
			#als = muscle.alignSequences(seqs, exepath=os.path.expanduser("/Users/dad/develop/muscle3.8.13/muscle"))
			als = muscle.alignSequences(seqs)
			#print als
			self.assertTrue(len(als) == len(seqs))
			for (i, s) in enumerate(seqs):
				#print(als[i])
				#print(s)
				self.assertTrue(s == als[i].replace("-",''))
		except muscle.MuscleError as me:
			self.assertTrue(False)

if __name__=="__main__":
	unittest.main(verbosity=2)