
## Cribbed from greylag, a collection of programs for MS/MS protein analysis
## Copyright (C) 2006-2008  Stowers Institute for Medical Research
## http://greylag.org/

# FIX: is this mono or avg?  (possible error here is ~0.0007 amu)
proton_mass =   1.007276
electron_mass = 0.000549                # ?

# Reference values from NIST (http://physics.nist.gov/PhysRefData/)
monoisotopic_atomic_mass = {
    'H' :  1.00782503214,
    'C' : 12.00000000,
    'N' : 14.00307400529,
    'O' : 15.994914622115,
    'P' : 30.9737615120,
    'S' : 31.9720706912,
    }

# most prevalent only (1 in 1000)
isotopic_atomic_mass = {                # prevalence (in %)
    'C13' : 13.003354837810,            # 1.078
    'N15' : 15.00010889849,             # 0.3687
    'O18' : 17.99916049,                # 0.20514
    'S33' : 32.9714585012,              # 0.762
    'S34' : 33.9678668311,              # 4.2928
    }

average_atomic_mass = {
    'H' :  1.007947,
    'C' : 12.01078,
    'N' : 14.00672,
    'O' : 15.99943,
    'P' : 30.9737612,
    'S' : 32.0655,
    }

heavy_atomic_mass = {
    'H' :  1.007947, # light
    'C' : 13.003354837810, # heavy
    'N' : 15.00010889849, # heavy
    'O' : 15.99943, # light
    'P' : 30.9737612, # light
    'S' : 32.0655 # light
    }


# The xtandem average residue masses are about 0.002 amu higher than those
# calculated directly from the above average atomic masses.  None of the
# chemists consulted knew of any reason why, aside from lack of precision in
# the average atomic mass estimates.  This shouldn't matter very much, as
# fragmentation calculations should all be monoisotopic, and we can always
# widen the parent tolerance window a bit.


def formula_mass(formula, atomic_mass=monoisotopic_atomic_mass):
    """Return the mass of formula, using the given mass regime (monoisotopic
    by default).

    >>> formula_mass('H2O', { 'H':1, 'O':16 })
    18
    >>> # monoisotopic mass of glycine
    >>> str(round(formula_mass('C2H3ON'), 4))
    '57.0215'

    """
    parts = [ p or '1' for p in re.split(r'([A-Z][a-z]*)', formula)[1:] ]
    # parts for glycine = ['C', '2', 'H', '3', 'O', '1', 'N', '1']
    return sum(atomic_mass[parts[i]] * int(parts[i+1]) for i in range(0, len(parts), 2))

# FIX: selenocysteine (U), etc
# residue -> formula
residue_formula = {
    'A' : "C3H5ON",
    'C' : "C3H5ONS",
    'D' : "C4H5O3N",
    'E' : "C5H7O3N",
    'F' : "C9H9ON",
    'G' : "C2H3ON",
    'H' : "C6H7ON3",
    'I' : "C6H11ON",
    'K' : "C6H12ON2",
    'L' : "C6H11ON",
    'M' : "C5H9ONS",
    'N' : "C4H6O2N2",
    'P' : "C5H7ON",
    'Q' : "C5H8O2N2",
    'R' : "C6H12ON4",
    'S' : "C3H5O2N",
    'T' : "C4H7O2N",
    'V' : "C5H9ON",
    'W' : "C11H10ON2",
    'Y' : "C9H9O2N",
    }

residues = sorted(residue_formula.keys())
residues_w_brackets = residues + ['[', ']']


# [0][1] -> 'H' -> fragment mass of H for regime 0
mass_regime_atomic_masses = []

## End greylag crib

mass_cache = {}

def getMass(peptide, label, masses):
	key = "%s@@%s" % (peptide, label)
	try:
		mass = mass_cache[key]
	except KeyError:
		## Add mass to mass_cache
		mass = sum([masses[x] for x in peptide])
		mass_cache[key] = mass
	return mass