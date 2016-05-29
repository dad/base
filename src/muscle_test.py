import sys, os, math, string, random, unittest
import muscle, translate, stats

def mut(x, pmut, pgap):
	res = x
	if random.random() < pmut:
		res = random.choice(translate.AAs())
	if random.random() < pgap:
		res = res + '-'
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
		# No gaps: pgap = 0.0
		others = [''.join([mut(x, 0.2,0.0) for x in s1]) for i in range(9)]
		seqs = [s1]+others
		res = False
		try:
			als = muscle.alignSequences(seqs)
			#print als
			self.assertTrue(len(als) == len(seqs))
			for (i, s) in enumerate(seqs):
				self.assertTrue(s == als[i].replace("-",''))
		except muscle.MuscleError as me:
			self.assertTrue(False)

if __name__=="__main__":
	unittest.main(verbosity=2)