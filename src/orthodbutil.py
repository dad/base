import os, random, string, sys

def headerDict(header_string):
	"""Return dictionary from header string"""
	# >7244:002ce8 FBpp0236088 gene=FBgn0208790 orthodb8_OG=EOG8MGTH1 orthodb8_level=32281 organism_name=`Drosophila virilis` uniprot_de=`GJ21671`
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
