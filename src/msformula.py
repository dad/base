import msutil, sys
from argparse import ArgumentParser

# Computes mass and molecular formula for a given peptide
# DAD: implement.

if __name__=="__main__":
	# DAD: add "--minus-water" option.
	parser = ArgumentParser(usage="%prog [-i fasta] [-s sequence]")
	parser.add_argument("-p", "--peptide", dest="peptides", action="append", default=[], help="peptides to analyze")
	parser.add_argument("-w", "--nowater", dest="nowater", action="store_true", default=False, help="subtract water from result?")
	options = parser.parse_args()

	water_adjustment = 0.0
	if options.nowater:
		print "# Removing H2O from each peptide"
		water_adjustment = msutil.formulaMass('H2O')
	
	for pep in options.peptides:
		components = msutil.peptideComponents(pep)
		if options.nowater:
			components['H'] = components['H']-2
			components['O'] = components['O']-1
		component_string = ' '.join(['{}({})'.format(comp,components[comp]) for comp in 'HCNOPS' if comp in components.keys()])
		print "{} = {:.5f}; {:s}".format(pep, msutil.getPeptideMass(pep)-water_adjustment, component_string)
