import sys, os, math, string, re, unittest

## Cribbed from greylag, a collection of programs for MS/MS protein analysis
## Copyright (C) 2006-2008  Stowers Institute for Medical Research
## http://greylag.org/

# FIX: is this mono or avg?  (possible error here is ~0.0007 amu)
class MSConstants:
	proton_mass =   1.00727646677
	electron_mass = 0.00054857991

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
		'*' : ""
		}
	
	water = "H2O"

	residues = sorted(residue_formula.keys())
	residues_w_brackets = residues + ['[', ']']
	
	# Masses for the TMT 6-plex reporter ions
	TMT6plex_reporter_ion_masses = {
		126:126.127725,
		127:127.124760,
		128:128.134433,
		129:129.131468,
		130:130.141141,
		131:131.138176
	}

def ppm(actual, observed):
	return 1e6*(observed - actual)/float(actual)

class MS2Spectrum(object):
	"""A set of (centroid,intensity) pairs with other information."""
	def __init__(self):
		self._centroid_pair_list = []
		self._charge = None
		self._MS1 = None
		self._mz = None
		self._scan = None
	
	def init(self, pair_list, charge, ms1, mz, scan):
		self._centroid_pair_list = pair_list
		self._charge = charge
		self._MS1 = ms1
		self._mz = mz
		self._scan = scan
		
	def readFromDTA(self, stream):
		self._centroid_pair_list = []
		# Read MS1 and charge
		flds = stream.readline().strip().split(' ')
		self._charge = int(flds[1])
		self._MS1 = float(flds[0])
		for line in stream:
			if line.strip() != '':
				flds = line.strip().split(' ')
				self._centroid_pair_list.append((float(flds[0]), float(flds[1])))
		self._centroid_pair_list.sort()
	
	def writeToDTA(self, outs):
		outs.write("{me._MS1} {me._charge}\n".format(me=self))
		for peak in self.peaks:
			outs.write("{p.mz} {p.intensity}\n".format(p=peak))
	
	def closestPeak(self, mz, ppmresolution):
		"""Find the peak that is closest in m/z to mz. Return None if there is no peak within +/- resolution in ppm."""
		# First just do brutally slow but guaranteed-correct linear search.
		closest_index = -1
		smallest_diff = 1e6
		for (i,x) in enumerate(self._centroid_pair_list):
			(mass, intens) = x
			diff = mass - mz
			if abs(diff) < smallest_diff:
				closest_index = i
				smallest_diff = abs(diff)
			else:
				if diff > 0:
					# We have moved beyond where a closer match is possible, because mz is ordered.
					break
		# Determine whether this match is close enough
		closepair = self._centroid_pair_list[closest_index]
		res = None
		if abs(ppm(mz,closepair[0])) <= ppmresolution:
			res = Peak(closepair[0], closepair[1])
		return res
	
	@property
	def scan(self):
		return self._scan
	
	@scan.setter
	def scan(self, s):
		self._scan = s
	
	@property
	def charge(self):
		return self._charge
	
	@charge.setter
	def scan(self, chg):
		self._charge = chg
	
	@property
	def MS1(self):
		return self._MS1
	
	@MS1.setter
	def scan(self, MS1):
		self._MS1 = MS1
	
	@property
	def peaks(self):
		for (m,i) in self._centroid_pair_list:
			yield Peak(m,i)
	
	def __str__(self):
		return ",".join(["({m},{z})".format(m=m, z=z) for (m,z) in self._centroid_pair_list])
		
class Peak(object):
	"""A simple peak object."""
	def __init__(self, mz, intensity):
		self.mz = mz
		self.intensity = intensity
	
	def __str__(self):
		s = "{},{}".format(self.mz, self.intensity)
		return s

# The xtandem average residue masses are about 0.002 amu higher than those
# calculated directly from the above average atomic masses.  None of the
# chemists consulted knew of any reason why, aside from lack of precision in
# the average atomic mass estimates.  This shouldn't matter very much, as
# fragmentation calculations should all be monoisotopic, and we can always
# widen the parent tolerance window a bit.

def formulaComponents(formula):
    # parts for glycine = ['C', '2', 'H', '3', 'O', '1', 'N', '1']
    parts = [ p or '1' for p in re.split(r'([A-Z][a-z]*)', formula)[1:] ]
    return dict([(parts[i],int(parts[i+1])) for i in range(0,len(parts),2)])

def peptideComponents(aa_sequence):
	component_counts = dict([(component, 0) for component in MSConstants.monoisotopic_atomic_mass.keys()])
	if aa_sequence != '':
		for aa in aa_sequence:
			part_dict = formulaComponents(MSConstants.residue_formula[aa])
			for k in part_dict.keys():
				component_counts[k] += part_dict[k]
		# add water to account for N and C termini
		component_counts['H'] += 2
		component_counts['O'] += 1
	return component_counts
		

def formulaMass(formula, atomic_mass=MSConstants.monoisotopic_atomic_mass):
    """Return the mass of formula, using the given mass regime (monoisotopic
    by default).

    >>> formulaMass('H2O', { 'H':1, 'O':16 })
    18
    >>> # monoisotopic mass of glycine
    >>> str(round(formulaMass('C2H3ON'), 4))
    '57.0215'
    """
    #parts = [ p or '1' for p in re.split(r'([A-Z][a-z]*)', formula)[1:] ]
    parts = formulaComponents(formula)
    mass = sum([atomic_mass[comp]*n for (comp,n) in parts.items()])
    return mass #sum(atomic_mass[parts[i]] * int(parts[i+1]) for i in range(0, len(parts), 2))

## End greylag crib

def getPeptideMass(aa_sequence):
	mass = sum([formulaMass(MSConstants.residue_formula[x]) for x in aa_sequence]) + formulaMass(MSConstants.water)
	return mass


##########
# Test cases
##########

class test001(unittest.TestCase):
	def test_run(self):
		# Test closest peak
		ms2 = MS2Spectrum()
		ms2.init([(1,1), (2,2), (3,3), (4,4)], 2, 1000, 500, "the_scan")
		pk = ms2.closestPeak(2.0001, 100)
		self.assertTrue(pk.mz==2)
		
class test002(unittest.TestCase):
	def test_run(self):
		# Test closest peak, outside resolution
		ms2 = MS2Spectrum()
		ms2.init([(1,1), (2,2), (3,3), (4,4)], 2, 1000, 500, "the_scan")
		pk = ms2.closestPeak(2.3, 0.1)
		# No peak within 0.1 of 2.3.
		self.assertTrue(pk == None)
		
class test003(unittest.TestCase):
	def test_run(self):
		# Test formula mass
		m = formulaMass('H2O', {'H':1, 'O':16})
		self.assertTrue(m == 18)
		
class test004(unittest.TestCase):
	def test_run(self):
		# Test reading and writing from DTA
		ms2 = MS2Spectrum()
		ms2.init([(1,1), (2,2), (3,3), (4,4)], 2, 1000, 500, "the_scan")
		tmpfile = os.tmpfile()
		ms2.writeToDTA(tmpfile)
		ms2_copy = MS2Spectrum()
		tmpfile.flush()
		tmpfile.seek(0)
		ms2_copy.readFromDTA(tmpfile)
		peaks = [p for p in ms2.peaks]
		peaks_copy = [p for p in ms2_copy.peaks]
		self.assertTrue(ms2.charge == ms2_copy.charge)
		self.assertTrue(ms2.MS1 == ms2_copy.MS1)
		for (p1,p2) in zip(peaks, peaks_copy):
			self.assertTrue(p1.mz==p2.mz)
			self.assertTrue(p1.intensity==p2.intensity)
		
		
		

if __name__=="__main__":
	unittest.main(verbosity=2)
