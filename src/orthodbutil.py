import os, random, string, sys
import ete3

def headerDict(header_string):
	"""Return dictionary from header string"""
	# >7244:002ce8 FBpp0236088 gene=FBgn0208790 orthodb8_OG=EOG8MGTH1 orthodb8_level=32281 organism_name=`Drosophila virilis` uniprot_de=`GJ21671`
	# Handling awful cases like uniprot_de=`Probable tRNA 2`-O-ribose methyltransferase`
	header_string = header_string.replace("`-", "'-")
	quote_split = header_string.split("`")
	def garble(x):
		return x.replace(" ", "@*#/*")
	def degarble(x):
		return x.replace("@*#/*", " ")
	reform = quote_split[0]
	xi = 1
	while xi < len(quote_split):
		# string in quotes
		reform += garble(quote_split[xi])
		# next string
		reform += quote_split[xi+1]
		xi = xi+2
	# Split 

	d = {}
	for entry in reform.split():
		if '=' in entry:
			sp = entry.split('=')
			d[sp[0]] = degarble(sp[1])
	
	return d

def translateHeader(header_string):
	# >10224:0029f1 "pub_gene_id":"Sakowv30031477m", "pub_og_id":"EOG091G08IZ", "og_name":"guanine nucleotide binding protein-like 3 (nucleolar) ","level":33208
	# first number is an NCBI taxon id
	# second number is a unique hexadecimal id
	def unquote(x):
		x = x.replace('"','')
		x = x.replace("'","")
		return x
	# Find first space, indicating end of IDs
	id_end = header_string.find(' ')
	id_flds = header_string[1:id_end].split(":")
	ncbi_taxon_id = int(id_flds[0])
	taxon_name = ete3.NCBITaxa().get_taxid_translator([ncbi_taxon_id])[ncbi_taxon_id]
	rest = header_string[id_end:]
	brace_begin = rest.find('{')
	if brace_begin>0:
		rest = rest[(brace_begin+1):]
		rest = rest.replace("}",'')
	flds = [x.strip() for x in rest.split(',')]
	#print(header_string[id_end:])
	#print(flds)
	res_dict = dict([(unquote(x.strip().split(':')[0]),unquote(x.strip().split(':')[1])) for x in rest.split(',')])
	res_dict["taxon"] = taxon_name
	#print(res_dict)
	#print(res_dict)
	return res_dict
