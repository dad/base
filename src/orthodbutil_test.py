#! python

import os, sys, math, unittest
import orthodbutil

class test001(unittest.TestCase):
	def test_run(self):
		"""Handle quoted names"""
		header = "7244:002ce8 FBpp0236088 gene=FBgn0208790 orthodb8_OG=EOG8MGTH1 orthodb8_level=32281 organism_name=`Drosophila virilis` uniprot_de=`GJ21671`"
		d = orthodbutil.headerDict(header)
		self.assertTrue(d['organism_name'] == 'Drosophila virilis')

	def test_prime(self):
		"""Handle prime quote"""
		header = "1111077:0001b9 M1WGG3_CLAP2 gene=M1WGG3 orthodb8_OG=EOG8JHB13 orthodb8_level=4751 organism_name=`Claviceps purpurea 20.1` uniprot_de=`Probable tRNA 2`-O-ribose methyltransferase`"
		d = orthodbutil.headerDict(header)
		self.assertTrue(d['organism_name'] == 'Claviceps purpurea 20.1')

	def test_taxa(self):
		"""Handle taxon name translation and updated OrthoDB format"""
		header = '>10224:0029f1 "pub_gene_id":"Sakowv30031477m", "pub_og_id":"EOG091G08IZ", "og_name":"guanine nucleotide binding protein-like 3 (nucleolar) ","level":33208'
		d = orthodbutil.translateHeader(header)
		self.assertTrue(d['taxon'] == 'Saccoglossus kowalevskii')

	def test_taxa_brace(self):
		"""Handle taxon name translation and updated OrthoDB format with braces"""
		header = '>9531:00482e {"pub_gene_id":"105598636", "pub_og_id":"EOG091G08IZ", "og_name":"guanine nucleotide binding protein-like 3 (nucleolar) ","level":33208, "description":"guanine nucleotide binding protein-like 3 (nucleolar)-like"}'
		d = orthodbutil.translateHeader(header)
		#print(d)
		self.assertTrue(d['taxon'] == 'Cercocebus atys')
		self.assertTrue(d['pub_gene_id'] == '105598636')

if __name__=="__main__":
	unittest.main(verbosity=2)
