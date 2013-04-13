import msutil

# Computes mass and molecular formula for a given peptide
# DAD: implement.

if __name__=="__main__":
	# DAD: add "--minus-water" option.
	peps = sys.argv[1:]
	for pep in peps:
		components = peptideComponents(pep)
		component_string = ' '.join(['{}({})'.format(comp,components[comp]) for comp in 'HCNOPS' if comp in components.keys()])
		print "{} = {:.5f}; {:s}".format(pep, getPeptideMass(pep), component_string)
