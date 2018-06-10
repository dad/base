import os, random, string, sys, re
import ete3

MISSING_TAXON = "<unknown>"

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
	# y = re.compile('(\"(.+)\":\"?(.+)\"?)?,')
	def unquote(x):
		x = x.replace('"','')
		x = x.replace("'","")
		return x
	def get_key_colon_value(x):
		y = x.strip().split(":")
		if not len(y)==2:
			print(x,y)
		return (y[0],y[1])
	# Find first space, indicating end of IDs
	entry_start = header_string.find('>')
	if entry_start<0:
		entry_start=0
	else:
		entry_start = entry_start+1
	id_end = header_string.find(' ')
	id_flds = header_string[entry_start:id_end].split(":")
	ncbi_taxon_id = int(id_flds[0])
	taxon_dict = ete3.NCBITaxa().get_taxid_translator([ncbi_taxon_id])
	if ncbi_taxon_id in taxon_dict:
		taxon_name = taxon_dict[ncbi_taxon_id]
	else:
		#print(id_flds)
		taxon_name = MISSING_TAXON
	rest = header_string[id_end:]
	brace_begin = rest.find('{')
	if brace_begin>0:
		rest = rest[(brace_begin+1):]
		rest = rest.replace("}",'')
	#print(rest)
	y = re.compile('("([^":]*)":"([^"]|"")*")|("([^"]*)":(\d+))')
	def pickcolon(x):
		y = [e for e in x if e.find(":")>0]
		return y[0]
	flds = [pickcolon(x).split(":") for x in y.findall(rest)]
	#flds = [get_key_colon_value(x) for x in rest.split(',')]
	#print(header_string[id_end:])
	#print(flds)
	#print(flds)
	res_dict = dict([(unquote(x[0]),unquote(x[1])) for x in flds])
	res_dict["taxon"] = taxon_name
	#print(res_dict)
	#print(res_dict)
	return res_dict
