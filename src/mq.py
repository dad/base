import math, os
import stats, util, na

class ProteinQuant(object):
	def __init__(self, id):
		self.id = id
		self.peptide_dict = {}
		self.significance = None

	def add(self, peptide_quant):
		# Confirm that this peptide is in my collection; if not, add it
		try:
			pep = self.peptide_dict[peptide_quant.key]
		except KeyError:
			self.peptide_dict[peptide_quant.key] = peptide_quant
		# Add this protein as a client of this peptide
		peptide_quant.addProtein(self.id)

	def merge(self, protein_quant):
		# Merge non-destructively
		id = self.id
		if protein_quant.id != self.id:
			id = "{}+{}".format(self.id, protein_quant.id)
		merge_pq = self.__class__(id)
		# Copy my peptides; merge others
		merge_pq.peptide_dict = dict(self.peptide_dict.items())
		for k in protein_quant.peptide_dict.keys():
			try:
				pep = self.peptide_dict[k]
				merged_pep = pep.merge(protein_quant.peptide_dict[k])
				merge_pq.peptide_dict[k] = merged_pep
			except KeyError:
				merge_pq.peptide_dict[k] = protein_quant.peptide_dict[k].copy()
		return merge_pq

	def copy(self):
		# Merge non-destructively
		pq = self.__class__(self.id)
		pq.peptide_dict = dict(self.peptide_dict.items())
		return pq

	def getNormalizedHeavyLightRatioSummary(self):
		acc = stats.LogAccumulator(store=True)
		for pep in self.peptide_dict.values():
			acc.addAll(pep.heavy_light_normalized_ratio_list)
		return acc.summary

	@property
	def normalized_ratio(self):
		med = self.getNormalizedHeavyLightRatioSummary().median
		return med

	def getHeavyLightRatioSummary(self):
		acc = stats.LogAccumulator(store=True)
		for pep in self.peptides:
			acc.addAll(pep.heavy_light_ratio_list)
		s = acc.summary
		return s

	@property
	def ratio(self):
		res = None
		med = self.getHeavyLightRatioSummary().median
		if not na.isNA(med):
			res = math.exp(med)
		return res

	def getIntensitySummary(self):
		acc = stats.Accumulator(store=True)
		acc.addAll([pep.intensity for pep in self.peptides])
		return acc.summary

	@property
	def intensity(self):
		return self.getIntensitySummary().sum

	def getHeavyIntensitySummary(self):
		acc = stats.Accumulator(store=True)
		acc.addAll([pep.heavy_intensity for pep in self.peptides])
		return acc.summary

	@property
	def heavy_intensity(self):
		return self.getHeavyIntensitySummary().sum

	def getLightIntensitySummary(self):
		acc = stats.Accumulator(store=True)
		acc.addAll([pep.light_intensity for pep in self.peptides])
		return acc.summary

	@property
	def light_intensity(self):
		return self.getLightIntensitySummary().sum

	def getMSMSCountSummary(self):
		acc = stats.Accumulator(store=True)
		acc.addAll([pep.msms_count for pep in self.peptide_dict.values()])
		return acc.summary

	@property
	def msms_count(self):
		return int(self.getMSMSCountSummary().sum)

	def getTopN(self, n=3, fetch=lambda x: x.intensity):
		pep_int_dict = {}
		for pep in self.peptides:
			pep_inte = fetch(pep)
			try:
				inte = pep_int_dict[pep.sequence]
				if pep_inte>inte:
					pep_int_dict[pep.sequence] = pep_inte
			except KeyError:
				pep_int_dict[pep.sequence] = pep_inte
		# DAD: assert different peptides?
		peps = sorted(pep_int_dict.values(),reverse=True)
		if len(peps)>=n:
			res = sum(peps[0:n])
		else:
			res = None
		return res

	@property
	def top3(self):
		return self.getTopN(3)
	
	@property
	def heavy_top3(self):
		return self.getTopN(3, lambda x: x.heavy_intensity)
	
	@property
	def medium_top3(self):
		return self.getTopN(3, lambda x: x.medium_intensity)
	
	@property
	def light_top3(self):
		return self.getTopN(3, lambda x: x.light_intensity)
	

	@property
	def degeneracy(self):
		min_proteins = None
		for pep in self.peptide_dict.values():
			n = len(pep._parent_proteins)
			if min_proteins is None or (n < min_proteins):
				min_proteins = n
		return min_proteins
	
	@property
	def peptides(self):
		for pep in self.peptide_dict.values():
			yield pep
	
	@property
	def n_peptides(self):
		return len(self.peptide_dict)
	
	@property
	def ratios(self):
		for pep in self.peptides:
			for r in pep.ratios:
				yield r
	
	@property
	def ratio_count(self):
		rc = 0
		for pep in self.peptides:
			rc += pep.ratio_count
		return rc

	@property
	def normalized_ratios(self):
		for pep in self.peptides:
			for r in pep.normalized_ratios:
				yield r
	
	@property
	def key(self):
		return self.id

class ProteinQuant3(ProteinQuant):
	def __init__(self, id):
		ProteinQuant.__init__(self, id)
		self._ratio_getters = {"hl":self.getHeavyLightRatioSummary, "hm":self.getHeavyMediumRatioSummary, "ml":self.getMediumLightRatioSummary}
		self._normalized_ratio_getters = {"hl":self.getNormalizedHeavyLightRatioSummary, "hm":self.getNormalizedHeavyMediumRatioSummary, "ml":self.getNormalizedMediumLightRatioSummary}

	def getRatioSummary(self, ratio_string):
		# Dispatch
		return self._ratio_getters[ratio_string]()

	def getIntensityRatioSummary(self, ratio_string):
		acc = stats.LogAccumulator(store=True)
		for pep in self.peptides:
			acc.addAll(pep.getIntensityRatios(ratio_string))
		return acc.summary

	def getNormalizedRatioSummary(self, ratio_string):
		# Dispatch
		return self._normalized_ratio_getters[ratio_string]()

	def getHeavyMediumRatioSummary(self):
		acc = stats.LogAccumulator(store=True)
		for pep in self.peptides:
			acc.addAll(pep.heavy_medium_ratio_list)
		return acc.summary

	@property
	def ratio_hm(self):
		res = None
		med = self.getHeavyMediumRatioSummary().median
		if not na.isNA(med):
			res = math.exp(med)
		return res

	def getMediumLightRatioSummary(self):
		acc = stats.LogAccumulator(store=True)
		for pep in self.peptides:
			acc.addAll(pep.medium_light_ratio_list)
		return acc.summary

	@property
	def ratio_ml(self):
		res = None
		med = self.getMediumLightRatioSummary().median
		if not na.isNA(med):
			res = math.exp(med)
		return res

	def getNormalizedHeavyMediumRatioSummary(self):
		acc = stats.LogAccumulator(store=True)
		for pep in self.peptides:
			acc.addAll(pep.heavy_medium_normalized_ratio_list)
		return acc.summary

	@property
	def normalized_ratio_hm(self):
		res = None
		med = self.getNormalizedHeavyMediumRatioSummary().median
		if not na.isNA(med):
			res = math.exp(med)
		return res

	def getNormalizedMediumLightRatioSummary(self):
		acc = stats.LogAccumulator(store=True)
		for pep in self.peptides:
			acc.addAll(pep.medium_light_normalized_ratio_list)
		return acc.summary

	@property
	def normalized_ratio_ml(self):
		res = None
		med = self.getNormalizedMediumLightRatioSummary().median
		if not na.isNA(med):
			res = math.exp(med)
		return res

	def getMediumIntensitySummary(self):
		acc = stats.Accumulator(store=True)
		acc.addAll([pep.medium_intensity for pep in self.peptides])
		return acc.summary

	@property
	def medium_intensity(self):
		return self.getMediumIntensitySummary().sum

		
class ProteinQuantTMT10:
	def __init__(self):
		pass

class PeptideData(object):
	def __init__(self):
		self.ratio_hl = None
		self.ratio_hl_normalized = None
		self.intensity_l = None
		self.intensity_h = None
		self._fraction = None
		self.sequence = None
		self.modified_sequence = None		
		self.mods = None
		self._msms_count = None

class PeptideData3(object):
	def __init__(self):
		self.ratio_hl = None
		self.ratio_hl_normalized = None
		self.ratio_hm = None
		self.ratio_hm_normalized = None
		self.ratio_ml = None
		self.ratio_ml_normalized = None
		self.intensity_l = None
		self.intensity_h = None
		self.intensity_m = None
		self._fraction = None
		self.sequence = None
		self.modified_sequence = None		
		self.mods = None
		self._msms_count = None

class PeptideDataTMT10(object):
	def __init__(self):
		self.precursor_ion_fraction = 
		self.reporter_intensities = []
		self.sequence = None
		self.modified_sequence = None		
		self.mods = None
		self._msms_count = None


class PeptideDetectionEvent(object):
	"""
	An object representing a single peptide detection event.
	DAD: track charge state.
	"""
	def __init__(self, prop_dict):
		# Takes a property dictionary
		self.prop_dict = prop_dict

# DAD: Tons of NA checking here. Is any of it really necessary? Any reason not to catch NA's as they're added, and exclude them?
class PeptideQuant(object):
	def __init__(self, key):
		self._key = key
		self.sequence = None
		self.modified_sequence = None
		self.heavy_light_ratio_list = []
		self.heavy_light_normalized_ratio_list = []
		self.intensity_l_list = []
		self.intensity_h_list = []
		self._msms_count = 0
		self._parent_proteins = set()
		self._fractions = set()

	@property
	def key(self):
		return self._key

	def add(self, pep_data):
		if self.sequence is None:
			self.sequence = pep_data.sequence
		else:
			assert(self.sequence == pep_data.sequence)
		if self.modified_sequence is None:
			self.modified_sequence = pep_data.modified_sequence
		# DAD: this fails when it encounters 
		# self.modified_sequence = _RFIVFHNEFSEHTFV(ph)ER_ 
		# pep_data.modified_sequence = _RFIVFHNEFSEHTFVER_
		# Not sure if this matters.
		#else:
		#	print self.modified_sequence, pep_data.modified_sequence
		#	assert(self.modified_sequence == pep_data.modified_sequence)
		self.heavy_light_ratio_list.append(pep_data.ratio_hl)
		self.heavy_light_normalized_ratio_list.append(pep_data.ratio_hl_normalized)
		self.intensity_l_list.append(pep_data.intensity_l)
		self.intensity_h_list.append(pep_data.intensity_h)
		self.msms_count += pep_data.msms_count
		self._fractions.add(pep_data._fraction)

	def addProtein(self, id):
		self._parent_proteins.add(id)

	def merge(self, pep_quant):
		# Merge non-destructively
		key = self.key
		if pep_quant.key != self.key:
			key = "{0}+{1}".format(self.key, pep_quant.key)
		pep = self.__class__(key)
		pep.heavy_light_ratio_list = self.heavy_light_ratio_list + pep_quant.heavy_light_ratio_list
		pep.heavy_light_normalized_ratio_list = self.heavy_light_normalized_ratio_list + pep_quant.heavy_light_normalized_ratio_list
		pep.intensity_l_list = self.intensity_l_list + pep_quant.intensity_l_list
		pep.intensity_h_list = self.intensity_h_list + pep_quant.intensity_h_list
		pep._parent_proteins = self._parent_proteins.union(pep_quant._parent_proteins)
		pep.msms_count = self.msms_count + pep_quant.msms_count
		pep._fractions = self._fractions.union(pep_quant._fractions)
		return pep

	def copy(self):
		pep = self.__class__(self.key)
		pep.heavy_light_ratio_list = self.heavy_light_ratio_list[:]
		pep.heavy_light_normalized_ratio_list = self.heavy_light_normalized_ratio_list[:]
		pep.intensity_l_list = self.intensity_l_list[:]
		pep.intensity_h_list = self.intensity_h_list[:]
		pep._parent_proteins = set(list(self._parent_proteins))
		pep.msms_count = self.msms_count
		pep._fractions = set(list(self._fractions))
		return pep

	def normalizeRatiosBy(self, ratio, norm_ratio):
		self.heavy_light_ratio_list = [x/ratio for x in self.heavy_light_ratio_list if not na.isNA(x)]
		self.heavy_light_normalized_ratio_list = [x/norm_ratio for x in self.heavy_light_normalized_ratio_list if not na.isNA(x)]

	def normalizeHeavyIntensity(self, weight):
		new_int = [intens/weight for intens in self.intensity_h_list if not na.isNA(intens)]
		self.intensity_h_list = new_int

	def normalizeLightIntensity(self, weight):
		new_int = [intens/weight for intens in self.intensity_l_list if not na.isNA(intens)]
		self.intensity_l_list = new_int

	def getHeavyLightRatioSummary(self):
		acc = stats.LogAccumulator(store=True)
		acc.addAll(self.heavy_light_ratio_list)
		return acc.summary

	@property
	def ratio_hl(self):
		res = None
		med = self.getHeavyLightRatioSummary().median
		if not na.isNA(med):
			res = math.exp(med)
		return res

	def getNormalizedHeavyLightRatioSummary(self):
		acc = stats.LogAccumulator(store=True)
		acc.addAll(self.heavy_light_normalized_ratio_list)
		return acc.summary

	@property
	def normalized_ratio_hl(self):
		res = None
		med = self.getNormalizedHeavyLightRatioSummary().median
		if not na.isNA(med):
			res = math.exp(med)
		return res
	
	@property
	def ratio_count_hl(self):
		return len(self.heavy_light_ratio_list)

	def getIntensitySummary(self):
		acc = stats.Accumulator(store=True)
		for (h_int, l_int) in zip(self.intensity_h_list,self.intensity_l_list):
			acc.add(h_int+l_int)
		return acc.summary

	@property
	def intensity(self):
		return self.light_intensity + self.heavy_intensity

	def getHeavyIntensitySummary(self):
		acc = stats.Accumulator(store=True)
		acc.addAll(self.intensity_h_list)
		return acc.summary

	@property
	def heavy_intensity(self):
		return self.getHeavyIntensitySummary().sum

	def getLightIntensitySummary(self):
		acc = stats.Accumulator(store=True)
		acc.addAll(self.intensity_l_list)
		return acc.summary

	@property
	def light_intensity(self):
		return self.getLightIntensitySummary().sum

	@property
	def msms_count(self):
		return int(self._msms_count)

	@msms_count.setter
	def msms_count(self, x):
		self._msms_count = x

	@property
	def parent_proteins(self):
		for p in self._parent_proteins:
			yield p

	@property
	def ratios(self):
		for x in self.heavy_light_ratio_list:
			yield x

	@property
	def normalized_ratios(self):
		for x in self.heavy_light_normalized_ratio_list:
			yield x

	@property
	def fractions(self):
		for s in self._fractions:
			yield s
	
	@property
	def n_fractions(self):
		return len(self._fractions)

	def __str__(self):
		line = self.key
		floatformat = "{:1.4f}"
		line += '\n\tH/L      \t' + '\t'.join([na.formatNA(r,format=floatformat) for r in self.heavy_light_ratio_list]) + '\n'
		line += '\tH/L Normal.\t' + '\t'.join([na.formatNA(r,format=floatformat) for r in self.heavy_light_normalized_ratio_list]) + '\n'
		line += '\tIntensity H\t' + '\t'.join([na.formatNA(r,format=floatformat) for r in self.intensity_h_list]) + '\n'
		line += '\tIntensity L\t' + '\t'.join([na.formatNA(r,format=floatformat) for r in self.intensity_l_list]) + '\n'
		line += '\tMS/MS Count\t{0:d}\n'.format(self.msms_count)
		return line

class PeptideQuant3(PeptideQuant):
	def __init__(self, key):
		PeptideQuant.__init__(self, key)
		self.heavy_medium_ratio_list = []
		self.heavy_medium_normalized_ratio_list = []
		self.medium_light_ratio_list = []
		self.medium_light_normalized_ratio_list = []
		self.intensity_m_list = []
		# Ratios straight from MaxQuant
		self._ratio_getters = {"hl":self.getHeavyLightRatioSummary, "hm":self.getHeavyMediumRatioSummary, "ml":self.getMediumLightRatioSummary}
		self._normalized_ratio_getters = {"hl":self.getNormalizedHeavyLightRatioSummary, "hm":self.getNormalizedHeavyMediumRatioSummary, "ml":self.getNormalizedMediumLightRatioSummary}		
		# Intensities
		self._intensities = {"h":self.intensity_h_list, "m":self.intensity_m_list, "l":self.intensity_l_list}

	def add(self, pep_data):
		PeptideQuant.add(self, pep_data)
		self.heavy_medium_ratio_list.append(pep_data.ratio_hm)
		self.heavy_medium_normalized_ratio_list.append(pep_data.ratio_hm_normalized)
		self.medium_light_ratio_list.append(pep_data.ratio_ml)
		self.medium_light_normalized_ratio_list.append(pep_data.ratio_ml_normalized)
		self.intensity_m_list.append(pep_data.intensity_m)
		
	def getRatioSummary(self, ratio_string):
		# Dispatch
		return self._ratio_getters[ratio_string]()
	
	# Intensity ratios: return at the peptide/event level
	def getIntensityRatios(self, ratio_string):
		rstr = ratio_string.lower()
		numerator = self._intensities[rstr[0]]
		denominator = self._intensities[rstr[1]]
		n = len(numerator)
		ratios = [None]*n
		assert(n==len(denominator))
		for i in range(n):
			if denominator[i] >0:
				ratios[i] = numerator[i]/float(denominator[i])
		return ratios

	def getIntensityRatioSummary(self, ratio_string):
		acc = stats.LogAccumulator(store=True)
		acc.addAll(self.getIntensityRatios(ratio_string))
		return acc.summary

	def getNormalizedRatioSummary(self, ratio_string):
		# Dispatch
		return self._normalized_ratio_getters[ratio_string]()

	def getHeavyMediumRatioSummary(self):
		acc = stats.LogAccumulator(store=True)
		acc.addAll(self.heavy_medium_ratio_list)
		return acc.summary

	@property
	def ratio_hm(self):
		res = None
		med = self.getHeavyMediumRatioSummary().median
		if not na.isNA(med):
			res = math.exp(med)
		return res

	@property
	def ratio_count_hm(self):
		return len(self.heavy_medium_ratio_list)

	def getMediumLightRatioSummary(self):
		acc = stats.LogAccumulator(store=True)
		acc.addAll(self.medium_light_ratio_list)
		return acc.summary

	@property
	def ratio_ml(self):
		res = None
		med = self.getMediumLightRatioSummary().median
		if not na.isNA(med):
			res = math.exp(med)
		return res

	@property
	def ratio_count_ml(self):
		return len(self.medium_light_ratio_list)

	def getNormalizedHeavyMediumRatioSummary(self):
		acc = stats.LogAccumulator(store=True)
		acc.addAll(self.heavy_medium_normalized_ratio_list)
		return acc.summary

	@property
	def normalized_ratio_hm(self):
		res = None
		med = self.getNormalizedHeavyMediumRatioSummary().median
		if not na.isNA(med):
			res = math.exp(med)
		return res

	def getNormalizedMediumLightRatioSummary(self):
		acc = stats.LogAccumulator(store=True)
		acc.addAll(self.medium_light_normalized_ratio_list)
		return acc.summary

	@property
	def normalized_ratio_ml(self):
		res = None
		med = self.getNormalizedMediumLightRatioSummary().median
		if not na.isNA(med):
			res = math.exp(med)
		return res

	def getIntensitySummary(self):
		acc = stats.Accumulator(store=True)
		for (h_int, m_int, l_int) in zip(self.intensity_h_list, self.intensity_m_list, self.intensity_l_list):
			acc.add(h_int+m_int+l_int)
		return acc.summary

	@property
	def intensity(self):
		return self.light_intensity + self.heavy_intensity + self.medium_intensity

	def getMediumIntensitySummary(self):
		acc = stats.Accumulator(store=True)
		acc.addAll(self.intensity_m_list)
		return acc.summary

	@property
	def medium_intensity(self):
		return self.getMediumIntensitySummary().sum

	def normalizeMediumIntensity(self, weight):
		new_int = [intens/weight for intens in self.intensity_m_list if not na.isNA(intens)]
		self.intensity_m_list = new_int

	def __str__(self):
		line = PeptideQuant.__str__(self)
		# DAD: implement
		return line


class EvidenceDescriptor(object):
	def __init__(self):
		self.filename = None
		self.invert = None
		self.experiment = None
		self.tags = None
		self.description = None
		self.tracked_modifications = None

	def __str__(self):
		line = "{}, {}, {}, {}, {}".format(self.filename, self.experiment, self.invert, self.tags, self.tracked_modifications)
		return line

class ExperimentEvidence(object):
	# Fields that are experiment-specific
	experiment_specific_header_flds = ["experiment","peptides.seq","razor.peptides.seq","unique.peptides.seq", \
									   "ratio.hl","ratio.hl.normalized","ratio.hl.significance.a","ratio.hl.significance.b",\
									   "ratio.hl.variability","ratio.hl.count","intensity","intensity.h","intensity.l"]
	# Fields that should be carried over
	carryover_flds = ["peptides.seq","razor.peptides.seq","unique.peptides.seq","intensity",\
					  "ratio.hl","ratio.hl.normalized","intensity.h","intensity.m","intensity.l"]
	# Fields that need to be renamed to indicate, e.g., the meaning of h and l (Heavy and Light...)
	rename_flds = ["ratio.hl.significance.a","ratio.hl.significance.b","ratio.hl.variability",\
				   "ratio.hl.count","intensity.h","intensity.l"]
	# Fields that may need to be inverted if invert = true
	invert_flds = ["ratio.hl","ratio.hl.normalized"]

	def __init__(self, tracked_mods=None, unique_matches=False):
		self.filename = None
		self.invert = False
		self.experiment = None
		self.protein_data = {}
		self.peptide_data = {}
		self.tracked_modifications = tracked_mods
		self.unique_matches_only = unique_matches
		self.rawfiles = set()

	def initFrom(self, ex_desc):
		self.filename = ex_desc.filename
		self.invert = ex_desc.invert
		self.experiment = ex_desc.experiment
		self.protein_data = {}
		self.peptide_data = {}
		self.tracked_modifications = ex_desc.tracked_modifications

	def readDescriptors(self, in_stream):
		line = in_stream.readline()
		while line and (util.isComment(line) or util.isBlank(line)):
			line = in_stream.readline()
		if line and not line.strip() == '':
			flds = line.strip().split("\t")
			self.filename = flds[0]
			self.invert = (flds[1].lower() in ['1','y','t'])
			self.experiment = flds[2]
			res = os.path.isfile(os.path.expanduser(self.filename))
			if not res:
				print "# Evidence file not found: {me.filename}".format(me=self)
		else:
			res = False
		return res

	def parseFields(self, flds, orf_dict=None):
		# Filter out all lines except those corresponding to the specified experiment
		# Return True if we decided to parse this line
		parsed = False
		if self.experiment is None or not flds.has_key('experiment') or self.experiment == "{0}".format(flds['experiment']):
			# Process this line
			# Skip lines corresponding to reverse or contaminant peptide matches.
			if flds['reverse'] == '+' or flds['contaminant'] == '+':
				# We parsed this line, even though we didn't like it.
				parsed = True
			else:
				# Process this line
				# First get list of protein accessions and split into id_flds variable
				all_id_flds = flds["proteins"].split(';')  # To get all proteins to which this peptide could belong, even if they are not uniquely supported by at least one peptide (e.g., they are paralogs which share this peptide)
				if self.unique_matches_only and len(all_id_flds)>1:
					return False
				# DAD: how to handle cases of multiple proteins identified?
				id_flds = flds["leading.razor.protein"].split(';')  # To get the "most likely" protein, in the sense of a protein whose identification is uniquely supported by at least one peptide
				id_flds = flds["proteins"].split(';')
				n_prots = len(id_flds)
				
				# Get raw file
				self.rawfiles.add(flds['raw.file'])

				# Get intensities and ratio
				pep_data = PeptideData3()
				for rat in ['hl','ml','hm']:
					rat_str = 'ratio.{}'.format(rat)
					rat_norm_str = 'ratio.{}.normalized'.format(rat)
					ratio = flds.get(rat_str)
					ratio_normalized = flds.get(rat_norm_str)
					if self.invert:
						raise "Inverting of triple-ratios not implemented yet"
					setattr(pep_data, 'ratio_{}'.format(rat), ratio)
					setattr(pep_data, 'ratio_{}_normalized'.format(rat), ratio_normalized)
				pep_data.intensity_h = flds.get('intensity.h')
				pep_data.intensity_m = flds.get('intensity.m')
				pep_data.intensity_l = flds.get('intensity.l')
				pep_data._fraction = flds.get('fraction')
				pep_data.sequence = flds['sequence']
				pep_data.modified_sequence = flds['modified.sequence']
				pep_data.msms_count = flds['ms.ms.count']
				pep_data.mods = flds['modifications']
				# Use sequence as unique identifier for this peptide
				pep_key = pep_data.sequence
				# Use modified peptide as a key unless we are not tracking these modifications
				if pep_data.mods != 'Unmodified': # and self.isTrackedModification(pep_data.mods):
					pep_key = flds['modified.sequence']
				# Retrieve the relevant PeptideQuant entry (the one with this sequence as its key), or make a new one
				try:
					pep_entry = self.peptide_data[pep_key]
				except KeyError:
					pep_entry = PeptideQuant3(pep_key)
					self.peptide_data[pep_key] = pep_entry
				pep_entry.add(pep_data)
				# Add this peptide to all identified proteins
				for orf in id_flds:
					# By default, extract everything
					extract = True
					# But if a set of ORFs has been provided, only extract peptides corresponding to those ORFs
					if not orf_dict is None:
						extract = orf_dict.has_key(orf)
					if extract:
						# Add data, or update it
						try:
							prot_entry = self.protein_data[orf]
						except KeyError:
							prot_entry = ProteinQuant3(orf)
							self.protein_data[orf] = prot_entry
						prot_entry.add(pep_entry)
				parsed = True
		return parsed

	def readData(self, orf_dict):
		inf = file(os.path.expanduser(self.filename),'r')
		dlr = util.DelimitedLineReader(inf, strip=False, header_name_processor=util.maxQuantHeader)
		# Read in the data
		max_lines = 1e9
		line = 0
		while not dlr.atEnd() and line < max_lines:
			line += 1
			flds = dlr.nextDict()
			self.parseFields(flds, orf_dict)
		if line == max_lines:
			print "# Warning: max_lines ({0}) exceeded in ExperimentEvidence.readData()".format(max_lines)

	@property
	def peptides(self):
		for p in self.peptide_data.values():
			yield p
	
	@property
	def proteins(self):
		for prot in self.protein_data.values():
			yield prot
	
	def getProtein(self, prot_id, default=None):
		return self.protein_data.get(prot_id, default)
	
	def getPeptide(self, pep_id, default=None):
		return self.peptide_data.get(pep_id, default)
	
	
	def normalizeRatiosBy(self, norm_prot):
		ratio_stats = norm_prot.getHeavyLightRatioSummary()
		ratio_norm_stats = norm_prot.getNormalizedHeavyLightRatioSummary()
		print ratio_stats.mean, ratio_norm_stats.mean
		#print ratio_stats.mean, ratio_norm_stats.mean
		for pep in self.peptide_data.values():
			pep.normalizeRatiosBy(ratio_stats.mean, ratio_norm_stats.mean)

	def normalizePeptideRatios(self, norm_prot, target_log_mean=0.0):
		# DAD: implement
		target_peps = [self.peptide_data[k] for k in peptide_keys]
		acc = stats.LogAccumulator(store=False)
		for p in target_peps:
			acc.addAll([r for r in p.ratios])
		mean = acc.mean
		sd = acc.sd
		log_mean = math.log(mean)
		for p_key in self.peptide_data.keys():
			p = self.peptide_data[p_key]
			ratlist = p.heavy_light_ratio_list
			for xi in range(len(ratlist)):
				#ratlist[xi] = math.exp((math.log(ratlist[xi]) - (log_mean - target_log_mean))/(log_sd/target_log_sd))
				try:
					norm_log_ratio = (math.log(ratlist[xi]) - (log_mean - target_log_mean))/(log_sd/target_log_sd)
					ratlist[xi] = math.exp(norm_log_ratio)
				except OverflowError:
					print "overflow! --", p_key, xi, ratlist[xi], math.log(ratlist[xi]), norm_log_ratio
		return

	def normalizePeptideIntensities(self, peptide_keys):
		# DAD: implement
		target_peps = [self.peptide_data[k] for k in peptide_keys]
		h_acc = stats.Accumulator(store=True)
		m_acc = stats.Accumulator(store=True)
		l_acc = stats.Accumulator(store=True)
		for p in target_peps:
			h_acc.addAll(p.intensity_h_list)
			m_acc.addAll(p.intensity_m_list)
			l_acc.addAll(p.intensity_l_list)
		median_h = h_acc.median
		median_m = m_acc.median
		median_l = l_acc.median
		# Normalize all peptides
		for p_key in self.peptide_data.keys():
			p = self.peptide_data[p_key]
			p.normalizeHeavyIntensity(median_h)
			p.normalizeMediumIntensity(median_m)
			p.normalizeLightIntensity(median_l)
	
	def isTrackedModification(self, mods):
		# DAD: implement
		return False

	def merge(self, exev, normalize):
		"""Merge the ExperimentEvidence exev with this one.  Create a new ExperimentEvidence object to hold the merged data."""
		# DAD: need to think about semantics...do we require that tracked modifications be the same?

		# Normalization: the idea is that, if all proteins were detected in both experiments, adding the intensities and counts would be fine.
		# However, if there is missing data, and two experiments differ, say, 10-fold in intensity, then proteins detected in only one
		# experiment will have intensities that have relative intensities off by 10-fold.
		# One solution is to normalize to the total intensity of the shared subset of peptides (even better, shared charge states...).
		# Another is to fit a model in which each protein has an abundance, intensity linearly turns into abundance, and there is a different
		# intercept term for each experiment.

		# What to do if mixing is an issue?  In other words, if one sample is mixed 1:1 and the other
		# is mixed 1:2, then how to merge?  Ignore for now: just add intensities.

		# Create new, merged experiment
		res = ExperimentEvidence()
		res.peptide_data = {}
		# Copy over self
		for k in self.peptide_data.keys():
			res.peptide_data[k] = self.peptide_data[k].copy()
		# Now the merge step will make copies
		for k in exev.peptide_data.keys():
			try:
				# If there's an existing matching peptide record, merge
				pep = res.peptide_data[k]
				res.peptide_data[k] = pep.merge(exev.peptide_data[k])
			except KeyError:
				# Otherwise, copy over new peptide record from exev
				res.peptide_data[k] = exev.peptide_data[k].copy()

		# Now merge protein data
		res.protein_data = {}
		# Copy over self
		for k in self.protein_data.keys():
			res.protein_data[k] = self.protein_data[k].copy()
		# Now the merge step will make copies
		for k in exev.protein_data.keys():
			try:
				prot = res.protein_data[k]
				res.protein_data[k] = prot.merge(exev.protein_data[k])
			except KeyError:
				res.protein_data[k] = exev.protein_data[k].copy()
		return res

	def calculateProteinRatioSignificance(self, num_nearest_proteins, ratio_field="ratio_hl_normalized", abundance_field="intensity"):
		# Limit significance calculations to proteins with ratios
		recs_with_ratios = [r for r in self.protein_data.values() if not na.isNA(getattr(r,ratio_field))]
		# Sort proteins by estimated abundance
		recs = sorted(recs_with_ratios, key=lambda x: getattr(x, abundance_field))
		rec_norm_hl = [getattr(r,ratio_field) for r in recs]
		n = num_nearest_proteins
		half_n = int(math.ceil(n/2.0))
		for ti in range(len(recs)):
			# fetch nearest N proteins by intensity
			if ti < half_n:
				beg_i = 0
				end_i = int(min(n, len(recs)))
			elif ti + half_n >= len(recs):
				beg_i = len(recs)-half_n
				end_i = len(recs)
			else:
				beg_i = ti - half_n
				end_i = ti + half_n
			log_ratios = [math.log(getattr(r,ratio_field)) for r in recs[beg_i:end_i]]
			my_log_ratio = math.log(getattr(recs[ti],ratio_field))
			(n,m,sd,se) = stats.StatsSummary(log_ratios)
			z = (my_log_ratio - m)/sd
			p_z = stats.Prob_Z(z)
			recs[ti].significance = p_z

	def __str__(self):
		line = ""
		if not self.filename is None:
			hl_str = 'H/L->H/L'
			if self.invert:
				hl_str = 'H/L->L/H'
			line += "{}\t{}\t{}, ".format(self.filename, hl_str, self.experiment)
		line += '{} peptides, {} proteins'.format(len(self.peptide_data), len(self.protein_data))
		return line
	
class ExperimentEvidenceFactory(object):
	"""Class that parses evidence files into multiple ExperimentEvidence objects.
	"""
	# A few possibilities:
	# 1) No experiments were specified to MQ.  Then there's no Experiment column in evidence.txt.
	#    - We should push all data into a single ExperimentEvidence object
	#    - experiment_names is None
	# 2) Experiments were specified to MQ, but no experiment names are provided to us.
	#    - We should learn the experiments on the fly, and report as many ExperimentEvidence objects
	#      as there are experiment names
	#    - experiment_names is [] (empty)
	# 3) Experiments were specified to MQ, and some/all of those names are provided to us.
	#    - We should extract only those named experiments
	#    - experiment_names is ['a','b',...]
	def __init__(self):
		self.experiments = None
		self.learn_experiments = False

	def load(self, evidence_descs, filter_tags, filter_experiments, unique_matches, tracked_modifications, orf_dict):
		# evidence_descs is a list of EvidenceDescriptor variables
		# Read experiments
		# Filter based on tags
		self.experiments = []
		evidence_fnames = {}
		tag_set = set(filter_tags)
		for ed in evidence_descs:
			# If no filter specified, or our experiment is in the filter, analyze it.
			if filter_experiments is None or filter_experiments == [] or ed.experiment in filter_experiments:
				shared_tags = list(tag_set.intersection(set(ed.tags)))
				if len(filter_tags) == 0 or len(shared_tags) > 0:
					# This experiment should be included
					exev = ExperimentEvidence(tracked_modifications, unique_matches)
					exev.initFrom(ed)
					self.experiments.append(exev)
					try:
						evidence_fnames[ed.filename].append(exev)
					except KeyError:
						evidence_fnames[ed.filename] = [exev]
		for fname in evidence_fnames.keys():
			exp_fname = os.path.expanduser(fname)
			assert os.path.isfile(exp_fname), "# No file found for {}\n".format(exp_fname)
			inf = file(exp_fname,'r')
			dlr = util.DelimitedLineReader(inf, strip=False, header_name_processor=util.maxQuantHeader)
			#print dlr.headers
			# Read in the data
			max_lines = 1e8
			line = 0
			while not dlr.atEnd() and line < max_lines:
				line += 1
				flds = dlr.nextDict()
				# Let all the experiments try to parse this line
				relevant_experiments = evidence_fnames[exp_fname]
				for ex in relevant_experiments:
					res = ex.parseFields(flds, orf_dict)
			inf.close()
		return self.experiments

def test001():
	def run(self):
		return True

if __name__=='__main__':
	# test cases
	harness = util.TestHarness()
	harness.add(test001())
	harness.run()
	
