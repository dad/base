import time, os, random, string, sys, math, traceback, unittest, pickle, gzip
import buildgene, collections, biofile, translate

class test001(unittest.TestCase):
	def setUp(self):
		chr_dict = biofile.readFASTADict("../data/S288C_reference_sequence_R64-2-1_20150113.fsa", buildgene.getChromosomeNumber)
		# Load sequence
		self.chromcol = buildgene.ChromosomeCollection()
		for chr_name in chr_dict:
			chromosome = buildgene.Chromosome(chr_name, chr_dict[chr_name])
			self.chromcol.addChromosome(chromosome)
		# Load genes
		with open("../data/saccharomyces_cerevisiae_R64-2-1_20150113.gff3",'r') as inf:
			self.chromcol.loadRegions(inf, buildgene.nameMapperFactory({"chrmt":"chrMito"}))
		# Load UTRs
		with open("./test/Pelechano2013_SuppData3_medianTranscripts.txt", 'r') as inf:
			self.chromcol.loadUTRLengths(inf)
		# Load counts
		#plus_counts_in_name = "./test/GSM1137787_Pab1_1_plus.sgr"
		plus_counts_in_name = "./test/GSM_sgr_test.txt"
		with open(plus_counts_in_name,'rt') as inf:
			self.chromcol.loadStrandCounts("+", inf, buildgene.nameMapperFactory({"Mito":"chrMito"}))

	def test_reading(self):
		"""Reading GFF feature"""
		gene = self.chromcol.getGene('YAL005C')
		#chrI	SGD	CDS	139503	141431	.	-	0	Parent=YAL005C_mRNA;Name=YAL005C_CDS;orf_classification=Verified
		self.assertTrue(gene.strand=='-')

	def test_reading_intron(self):
		"""Reading GFF feature with intron"""
		gene = self.chromcol.getGene('YAL003W')
			# chrI	SGD	CDS	142174	142253	.	+	0	Parent=YAL003W_mRNA;Name=YAL003W_CDS;orf_classification=Verified
			# chrI	SGD	CDS	142620	143160	.	+	1	Parent=YAL003W_mRNA;Name=YAL003W_CDS;orf_classification=Verified
			# chrI	SGD	intron	142254	142619	.	+	.	Parent=YAL003W_mRNA;Name=YAL003W_intron;orf_classification=Verified
		self.assertTrue(gene.strand=='+')
		regions = gene.getRegions()
		self.assertTrue(len(regions) == 5)
		regions = gene.getRegions("intron")
		self.assertTrue(len(regions) == 1)
		regions = gene.getRegions(['CDS','intron'])
		self.assertTrue(len(regions) == 3)
		regions = gene.getRegions(['CDS'])
		self.assertTrue(len(regions) == 2)

	def test_reading_5p3p(self):
		"""Reading Pelechano data"""
		gene = self.chromcol.getGene('YAL005C')
		utr5 = gene.getRegion("UTR5")
		# 1	-	YAL005C	Verified	90	62	146	52.6088301878386	34.8383748121315
		self.assertTrue(len(utr5)==62)
		for gene in self.chromcol.genes:
			if gene.strand == '+':
				self.assertTrue(gene.get5PrimeEnd()<gene.get3PrimeEnd())
			else:
				self.assertTrue(gene.strand == '-')
				self.assertTrue(gene.get5PrimeEnd()>gene.get3PrimeEnd())

class test002(unittest.TestCase):
	def test_getChromosomeNumber(self):
		"""Read fasta header and get the chromosome number"""
		test_line1 = '>ref|NC_001133| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=I]'
		self.assertTrue(buildgene.getChromosomeNumber(test_line1) == 'chrI')

		test_line2 = '>ref|NC_001134| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=II]'
		self.assertTrue(buildgene.getChromosomeNumber(test_line2) == 'chrII')

		test_line3 = '>ref|NC_001145| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XIII]'
		self.assertTrue(buildgene.getChromosomeNumber(test_line3) == 'chrXIII')
		
		test_line4 = '>ref|NC_001224| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [location=mitochondrion] [top=circular]'
		self.assertTrue(buildgene.getChromosomeNumber(test_line4) == 'chrMito')

class test003(unittest.TestCase):
	def setUp(self):
		chr_dict = biofile.readFASTADict("../data/S288C_reference_sequence_R64-2-1_20150113.fsa", buildgene.getChromosomeNumber)
		# Load sequence
		self.chromcol = buildgene.ChromosomeCollection()
		for chr_name in chr_dict:
			chromosome = buildgene.Chromosome(chr_name, chr_dict[chr_name])
			self.chromcol.addChromosome(chromosome)
		# Load genes
		with open("../data/saccharomyces_cerevisiae_R64-2-1_20150113.gff3",'r') as inf:
			self.chromcol.loadRegions(inf, buildgene.nameMapperFactory({"chrmt":"chrMito"}))
		# Load UTRs
		with open("./test/Pelechano2013_SuppData3_medianTranscripts.txt", 'r') as inf:
			self.chromcol.loadUTRLengths(inf)
		# Load counts
		#plus_counts_in_name = "./test/GSM1137787_Pab1_1_plus.sgr"
		plus_counts_in_name = "./test/GSM_sgr_test.txt"
		with open(plus_counts_in_name,'rt') as inf:
			self.chromcol.loadStrandCounts("+", inf, buildgene.nameMapperFactory({"Mito":"chrMito"}))

	def test_start_and_stop(self):
		"""
		Testing start and stop codons genome-wide
		"""
		stop_codons = ['TAG','TAA','TGA']
		start_codons = ['ATG','ATA'] # Some mitochondrial genes like Q0032 have an ATA start
		for gene in self.chromcol.genes:
			cds = gene.getCodingSequence()
			if len(cds)>0:
				self.assertTrue(cds[0:3] in start_codons)
				self.assertTrue(cds[-3:] in stop_codons)

	def test_utr_retrieval(self):
		"""
		Testing retrieval of specific 5'UTRs
		"""
		utr5 = self.chromcol.getGene('YAL061W').getRegion("UTR5")
		self.assertTrue(utr5.nucleotides[-10:] == 'ATATTTCAGA')
		utr5 = self.chromcol.getGene('YAL067C').getRegion("UTR5")
		#print(utr5)
		self.assertTrue(utr5.nucleotides[-10:] == 'AATAACATAC')
		#Test a mitochondrial gene
		#Test UTR3 of a gene with intron

class test004(unittest.TestCase):
	def setUp(self):
		chr_dict = biofile.readFASTADict("../data/S288C_reference_sequence_R64-2-1_20150113.fsa", buildgene.getChromosomeNumber)
		# Load sequence
		self.chromcol = buildgene.ChromosomeCollection()
		for chr_name in chr_dict:
			chromosome = buildgene.Chromosome(chr_name, chr_dict[chr_name])
			self.chromcol.addChromosome(chromosome)
		# Load genes
		with open("../data/saccharomyces_cerevisiae_R64-2-1_20150113.gff3",'r') as inf:
			self.chromcol.loadRegions(inf, buildgene.nameMapperFactory({"chrmt":"chrMito"}))
		# Load UTRs
		with open("./test/Pelechano2013_SuppData3_medianTranscripts.txt", 'r') as inf:
			self.chromcol.loadUTRLengths(inf)
		# Load counts
		#plus_counts_in_name = "./test/GSM1137787_Pab1_1_plus.sgr"
		plus_counts_in_name = "./test/GSM_sgr_test.txt"
		with open(plus_counts_in_name,'rt') as inf:
			self.chromcol.loadStrandCounts("+", inf, buildgene.nameMapperFactory({"Mito":"chrMito"}))
				
	def test_loading_SGR(self):
		"""Count loading"""
		chrom = self.chromcol.getChromosome('chrXVI')
		self.assertTrue(chrom.getCounts('+', 948042, 948042) == [4])

class test005(unittest.TestCase):
	def setUp(self):
		chr_dict = biofile.readFASTADict("../data/saccharomyces_cerevisiae_R64-1-1_20110208.fsa")#, buildgene.getChromosomeNumber)
		# Load sequence
		self.chromcol = buildgene.ChromosomeCollection()
		for chr_name in chr_dict:
			chromosome = buildgene.Chromosome(chr_name, chr_dict[chr_name])
			self.chromcol.addChromosome(chromosome)
		# Load genes
		with open("../data/saccharomyces_cerevisiae_R64-1-1_20110208.gff",'r') as inf:
			self.chromcol.loadRegions(inf)
		# Load UTRs
		with open("./test/Pelechano2013_SuppData3_medianTranscripts.txt", 'r') as inf:
			self.chromcol.loadUTRLengths(inf)
		# Load counts
		#plus_counts_in_name = "./test/GSM1137787_Pab1_1_plus.sgr"
		plus_counts_in_name = "./test/GSM_sgr_test.txt"
		with open(plus_counts_in_name,'rt') as inf:
			self.chromcol.loadStrandCounts("+", inf, buildgene.nameMapperFactory({"Mito":"chrMito"}))

	def test_load_sequence(self):
		self.assertTrue(self.chromcol.getChromosome('chrIII').sequence[0:20] == 'CCCACACACCACACCCACAC')

	def test_load_regions(self):
		gene = self.chromcol.getGene('YAL067C')
		#print(gene.get5PrimeEnd())
		# 5'UTR is 33 nucleotides
		# 5'CDS end is 9016
		self.assertTrue(gene.get5PrimeEnd() == 9016+33)
		# 3'UTR is 222 nucleotides
		# 3'CDS end is 7235
		self.assertTrue(gene.get3PrimeEnd() == 7235-222)

	def test_load_counts(self):
		# Mito	41	14
		# Mito	42	15
		# Mito	43	16
		self.assertTrue(self.chromcol['chrMito'].strand_counts['+'][40:43] == [14,15,16])
		self.assertTrue(self.chromcol['chrMito'].strand_counts['-'][40:43] == [0,0,0])

class test006(unittest.TestCase):
	def test_load_object(self):
		obj_fname = "../output/transcriptome-nucleotide-counts-pab1.pickle.gz"
		chromcol = pickle.load(gzip.open(obj_fname,'rb'))
		ssa1 = chromcol.getGene('YAL005C')
		cds = ssa1.getRegion("CDS")
		utr5 = ssa1.getRegion("UTR5")
		utr3 = ssa1.getRegion("UTR3")
		# 
		# print(sum(cds.counts),len(cds))
		# print(sum(utr5.counts),len(utr5))
		# print(sum(utr3.counts),len(utr3))
		# print(sum(cds.counts)/len(cds))
		# print(sum(utr5.counts)/len(utr5))
		# print(sum(utr3.counts)/len(utr3))
		self.assertTrue(len(cds)==1929)
		# More counts in CDS...
		self.assertTrue(sum(cds.counts)>sum(utr5.counts))
		# ...but lower density
		self.assertTrue(sum(cds.counts)/len(cds)<sum(utr5.counts)/len(utr5))
		# chrI	SGD	CDS	139503	141431	.	-	0	Parent=YAL005C;Name=YAL005C;gene=SSA1;Alias=SSA1,YG100;Ontology_term=GO:0000060,GO:		
		#   count data:
		# chrI	139503	42
		# chrI	139504	40
		# chrI	139505	39
		# chrI	139506	26
		# chrI	139507	25
		self.assertTrue(cds.counts[-1]==42)

class test007(unittest.TestCase):
	def setUp(self):
		chr_dict = biofile.readFASTADict("../data/saccharomyces_cerevisiae_R64-1-1_20110208.fsa")#, buildgene.getChromosomeNumber)
		# Load sequence
		self.chromcol = buildgene.ChromosomeCollection()
		for chr_name in chr_dict:
			chromosome = buildgene.Chromosome(chr_name, chr_dict[chr_name])
			self.chromcol.addChromosome(chromosome)
		# Load genes
		with open("./test/intron-noncoding.gff3",'r') as inf:
			self.chromcol.loadRegions(inf)
		# Load UTRs
		with open("./test/Pelechano2013_SuppData3_medianTranscripts.txt", 'r') as inf:
			self.chromcol.loadUTRLengths(inf)
	def test_coding_transcript(self):
		"""Coding transcript"""
		# Gene with introns and UTRs
		gene = self.chromcol.getGene("YAL030W")
		utr5 = gene.get5PrimeUTRSequence() # 22 nt
		utr3 = gene.get3PrimeUTRSequence() # 133 nt
		cds = gene.getCodingSequence()   # 
		seq1 = utr5 + cds + utr3
		seq2 = gene.getSplicedTranscript()
		self.assertTrue(seq1 == seq2)
	def test_noncoding_transcript(self):
		"""Noncoding transcript"""
		gene = self.chromcol.getGene("RDN37-1")
		self.assertTrue(gene.getSplicedTranscript() != '')

if __name__=="__main__":
	unittest.main(verbosity=2)

