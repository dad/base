import math, os
import stats, util

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
			id = "%s+%s" % (self.id, protein_quant.id)
		merge_pq = ProteinQuant(id)
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
		pq = ProteinQuant(self.id)
		pq.peptide_dict = dict(self.peptide_dict.items())
		return pq

	def getPeptides(self):
		return self.peptide_dict.values()

	def getNormalizedHeavyLightRatioSummary(self):
		acc = stats.Accumulator(store=True)
		for pep in self.peptide_dict.values():
			for ratio in pep.heavy_light_normalized_ratio_list:
				if not util.isNA(ratio):
					acc.add(math.log(ratio))
		s = acc.getSummary()
		return s

	def getNormalizedHeavyLightRatio(self):
		res = None
		med = self.getNormalizedHeavyLightRatioSummary().median
		if not util.isNA(med):
			res = math.exp(med)
		return res

	def getHeavyLightRatioSummary(self):
		acc = stats.Accumulator(store=True)
		for pep in self.peptide_dict.values():
			for ratio in pep.heavy_light_ratio_list:
				if not util.isNA(ratio):
					acc.add(math.log(ratio))
		s = acc.getSummary()
		return s

	def getHeavyLightRatio(self):
		res = None
		med = self.getHeavyLightRatioSummary().median
		if not util.isNA(med):
			res = math.exp(med)
		return res

	def getIntensitySummary(self):
		acc = stats.Accumulator(store=True)
		for pep in self.peptide_dict.values():
			intens = pep.getIntensity()
			if not util.isNA(intens):
				acc.add(intens)
		return acc.getSummary()

	def getIntensity(self):
		return self.getIntensitySummary().sum

	def getHeavyIntensitySummary(self):
		acc = stats.Accumulator(store=True)
		for pep in self.peptide_dict.values():
			intens = pep.getHeavyIntensity()
			if not util.isNA(intens):
				acc.add(intens)
		return acc.getSummary()

	def getHeavyIntensity(self):
		return self.getHeavyIntensitySummary().sum

	def getLightIntensitySummary(self):
		acc = stats.Accumulator(store=True)
		for pep in self.peptide_dict.values():
			intens = pep.getLightIntensity()
			if not util.isNA(intens):
				acc.add(intens)
		return acc.getSummary()

	def getLightIntensity(self):
		return self.getLightIntensitySummary().sum

	def getAbundance(self):
		return self.getIntensity()/len(self.getPeptides())

	def getHeavyAbundance(self):
		# DAD: should isolate heavy peptides.
		return self.getHeavyIntensity()/len(self.getPeptides())

	def getLightAbundance(self):
		# DAD: should isolate light peptides.
		return self.getLightIntensity()/len(self.getPeptides())

	def getMSMSCountSummary(self):
		acc = stats.Accumulator(store=True)
		for pep in self.peptide_dict.values():
			msms_count = pep.getMSMSCount()
			if not util.isNA(msms_count):
				acc.add(msms_count)
		return acc.getSummary()

	def getMSMSCount(self):
		return int(self.getMSMSCountSummary().sum)

	def getDegeneracy(self):
		min_proteins = None
		for pep in self.peptide_dict.values():
			n = len(pep.parent_proteins)
			if min_proteins is None or (n < min_proteins):
				min_proteins = n
		return min_proteins

	# Properties:
	# property(fget=None, fset=None, fdel=None, doc=None)
	intensity = property(getIntensity)
	intensity_h = property(getHeavyIntensity)
	intensity_l = property(getLightIntensity)
	abundance = property(getAbundance)
	abundance_h = property(getHeavyAbundance)
	abundance_l = property(getLightAbundance)
	msms_count = property(getMSMSCount)
	ratio_hl = property(getHeavyLightRatio)
	ratio_hl_normalized = property(getNormalizedHeavyLightRatio)
	degeneracy = property(getDegeneracy)

class PeptideData(object):
	def __init__(self):
		self.ratio_hl = None
		self.ratio_hl_normalized = None
		self.intensity_l = None
		self.intensity_h = None
		self._slice = None
		self.sequence = None
		self.mods = None
		self.msms_count = None

class PeptideDetectionEvent(object):
	"""
	An object representing a single peptide detection event.
	"""
	def __init__(self, prop_dict):
		# Takes a property dictionary
		self.prop_dict = prop_dict

class PeptideQuant(object):
	def __init__(self, key):
		self.key = key
		self.sequence = None
		self.mod_sequence = None
		self.heavy_light_ratio_list = []
		self.heavy_light_normalized_ratio_list = []
		self.intensity_l_list = []
		self.intensity_h_list = []
		self._msms_count = 0
		self.parent_proteins = set()
		self._slices = set()

	def add(self, pep_data):
		if self.sequence is None:
			self.sequence = pep_data.sequence
		else:
			assert(self.sequence == pep_data.sequence)
		self.heavy_light_ratio_list.append(pep_data.ratio_hl)
		self.heavy_light_normalized_ratio_list.append(pep_data.ratio_hl_normalized)
		self.intensity_l_list.append(pep_data.intensity_l)
		self.intensity_h_list.append(pep_data.intensity_h)
		self.msms_count += pep_data.msms_count
		self._slices.add(pep_data._slice)

	def addProtein(self, id):
		self.parent_proteins.add(id)

	def merge(self, pep_quant):
		# Merge non-destructively
		key = self.key
		if pep_quant.key != self.key:
			key = "{0}+{1}".format(self.key, pep_quant.key)
		pep = PeptideQuant(key)
		pep.heavy_light_ratio_list = self.heavy_light_ratio_list + pep_quant.heavy_light_ratio_list
		pep.heavy_light_normalized_ratio_list = self.heavy_light_normalized_ratio_list + pep_quant.heavy_light_normalized_ratio_list
		pep.intensity_l_list = self.intensity_l_list + pep_quant.intensity_l_list
		pep.intensity_h_list = self.intensity_h_list + pep_quant.intensity_h_list
		pep.parent_proteins = self.parent_proteins.union(pep_quant.parent_proteins)
		pep.msms_count = self.msms_count + pep_quant.msms_count
		pep._slices = self.slices.union(pep_quant._slices)
		return pep

	def copy(self):
		pep = PeptideQuant(self.key)
		pep.heavy_light_ratio_list = self.heavy_light_ratio_list[:]
		pep.heavy_light_normalized_ratio_list = self.heavy_light_normalized_ratio_list[:]
		pep.intensity_l_list = self.intensity_l_list[:]
		pep.intensity_h_list = self.intensity_h_list[:]
		pep.parent_proteins = set(list(self.parent_proteins))
		pep.msms_count = self.msms_count
		pep._slices = set(list(self._slices))
		return pep

	def normalizeRatiosBy(self, ratio, norm_ratio):
		self.heavy_light_ratio_list = [x/ratio for x in self.heavy_light_ratio_list if not util.isNA(x)]
		self.heavy_light_normalized_ratio_list = [x/norm_ratio for x in self.heavy_light_normalized_ratio_list if not util.isNA(x)]

	def normalizeHeavyIntensity(self, weight):
		new_int = [intens/weight for intens in self.intensity_h_list if not util.isNA(intens)]
		self.intensity_h_list = new_int

	def normalizeLightIntensity(self, weight):
		new_int = [intens/weight for intens in self.intensity_l_list if not util.isNA(intens)]
		self.intensity_l_list = new_int

	def getHeavyLightRatioSummary(self):
		acc = stats.Accumulator(store=True)
		for ratio in self.heavy_light_ratio_list:
			if not util.isNA(ratio):
				acc.add(math.log(ratio))
		return acc.getSummary()

	def getHeavyLightRatio(self):
		res = None
		med = self.getHeavyLightRatioSummary().median
		if not util.isNA(med):
			res = math.exp(med)
		return res

	def getNormalizedHeavyLightRatioSummary(self):
		acc = stats.Accumulator(store=True)
		for ratio in self.heavy_light_normalized_ratio_list:
			if not util.isNA(ratio):
				acc.add(math.log(ratio))
		return acc.getSummary()

	def getNormalizedHeavyLightRatio(self):
		res = None
		med = self.getNormalizedHeavyLightRatioSummary().median
		if not util.isNA(med):
			res = math.exp(med)
		return res

	def getIntensitySummary(self):
		acc = stats.Accumulator(store=True)
		for intens in self.intensity_l_list:
			if not util.isNA(intens):
				acc.add(intens)
		for intens in self.intensity_h_list:
			if not util.isNA(intens):
				acc.add(intens)
		return acc.getSummary()

	def getIntensity(self):
		return self.getHeavyIntensity() + self.getLightIntensity()

	def getHeavyIntensitySummary(self):
		acc = stats.Accumulator(store=True)
		for intens in self.intensity_h_list:
			if not util.isNA(intens):
				acc.add(intens)
		return acc.getSummary()

	def getHeavyIntensity(self):
		return self.getHeavyIntensitySummary().sum

	def getLightIntensitySummary(self):
		acc = stats.Accumulator(store=True)
		for intens in self.intensity_l_list:
			if not util.isNA(intens):
				acc.add(intens)
		return acc.getSummary()

	def getLightIntensity(self):
		return self.getLightIntensitySummary().sum

	def getMSMSCount(self):
		return int(self._msms_count)

	def getParentProteins(self):
		return list(self.parent_proteins)

	def getHeavyLightRatios(self):
		return self.heavy_light_ratio_list

	def getMSMSCount(self):
		return self._msms_count

	def setMSMSCount(self, x):
		self._msms_count = x

	def getSlices(self):
		return self._slices

	intensity = property(getIntensity)
	intensity_h = property(getHeavyIntensity)
	intensity_l = property(getLightIntensity)
	msms_count = property(getMSMSCount,setMSMSCount)
	ratio_hl = property(getHeavyLightRatio)
	ratio_hl_normalized = property(getNormalizedHeavyLightRatio)
	slices = property(getSlices)

	def __str__(self):
		line = self.key
		line += '\n\tH/L      \t' + '\t'.join(["%1.4f" % r for r in self.heavy_light_ratio_list]) + '\n'
		line += '\tH/L Normal.\t' + '\t'.join(["%1.4f" % r for r in self.heavy_light_normalized_ratio_list]) + '\n'
		line += '\tIntensity H\t' + '\t'.join(["%1.0f" % r for r in self.intensity_h_list]) + '\n'
		line += '\tIntensity L\t' + '\t'.join(["%1.0f" % r for r in self.intensity_l_list]) + '\n'
		line += '\tMS/MS Count\t{0:d}\n'.format(self.msms_count)
		return line

class EvidenceDescriptor(object):
	def __init__(self):
		self.filename = None
		self.invert = None
		self.experiment = None
		self.tags = None
		self.description = None

	def __str__(self):
		line = "%s, %s, %s, %s" % (self.filename, self.experiment, self.invert, self.tags)
		return line

class ExperimentEvidence(object):
	# Fields that are experiment-specific
	experiment_specific_header_flds = ["experiment","peptides.seq","razor.peptides.seq","unique.peptides.seq", \
									   "ratio.hl","ratio.hl.normalized","ratio.hl.significance.a","ratio.hl.significance.b",\
									   "ratio.hl.variability","ratio.hl.count","intensity","intensity.h","intensity.l"]
	# Fields that should be carried over
	carryover_flds = ["peptides.seq","razor.peptides.seq","unique.peptides.seq","intensity",\
					  "ratio.hl","ratio.hl.normalized","intensity.h","intensity.l"]
	# Fields that need to be renamed to indicate, e.g., the meaning of h and l (Heavy and Light...)
	rename_flds = ["ratio.hl.significance.a","ratio.hl.significance.b","ratio.hl.variability",\
				   "ratio.hl.count","intensity.h","intensity.l"]
	# Fields that may need to be inverted if invert = true
	invert_flds = ["ratio.hl","ratio.hl.normalized"]

	def __init__(self):
		self.filename = None
		self.invert = False
		self.experiment = None
		self.protein_data = {}
		self.peptide_data = {}
		self.tracked_modifications = None

	def initFrom(self, ex_desc):
		self.filename = ex_desc.filename
		self.invert = ex_desc.invert
		self.experiment = ex_desc.experiment
		self.protein_data = {}
		self.peptide_data = {}
		self.tracked_modifications = None

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
				print "# Evidence file not found: %s" % self.filename
		else:
			res = False
		return res

	def parseFields(self, flds, orf_dict):
		# Filter out all lines except those corresponding to the specified experiment
		# Return True if we decided to parse this line
		parsed = False
		if self.experiment is None or not flds.has_key('experiment') or self.experiment == "{0}".format(flds['experiment']):
		#if not util.isNA(self.experiment) and flds["experiment"] == self.experiment:
			# Process this line
			# Skip lines corresponding to reverse or contaminant peptide matches.
			if flds['reverse'] == '+' or flds['contaminant'] == '+':
				# We parsed this line, even though we didn't like it.
				parsed = True
			else:
				# Process this line
				# First get list of protein accessions and split into id_flds variable
				#id_flds = flds["proteins"].split(';')  # To get all proteins to which this peptide could belong, even if they are not uniquely supported by at least one peptide (e.g., they are paralogs which share this peptide)
				id_flds = flds["leading.razor.protein"].split(';')  # To get the "most likely" protein, in the sense of a protein whose identification is uniquely supported by at least one peptide
				n_prots = len(id_flds)

				ratio_hl = flds['ratio.hl']
				ratio_hl_normalized = flds['ratio.hl.normalized']
				# Get intensities and ratio
				pep_data = PeptideData()
				#intensity = flds['intensity']
				# Invert the H/L ratios if needed
				if self.invert:
					if not util.isNA(ratio_hl):
						pep_data.ratio_hl = 1.0/ratio_hl
					if not util.isNA(ratio_hl_normalized):
						pep_data.ratio_hl_normalized = 1.0/ratio_hl_normalized
					pep_data.intensity_h = flds['intensity.l']
					pep_data.intensity_l = flds['intensity.h']
				else:
					pep_data.ratio_hl = ratio_hl
					pep_data.ratio_hl_normalized = ratio_hl_normalized
					pep_data.intensity_h = flds['intensity.h']
					pep_data.intensity_l = flds['intensity.l']
				try:
					pep_data._slice = flds['gel.slice']
				except KeyError:
					pass
				pep_data.sequence = flds['sequence']
				pep_data.msms_count = flds['ms.ms.count']
				pep_data.mods = flds['modifications']
				# Use sequence as unique identifier for this peptide
				pep_key = pep_data.sequence
				# Use modified peptide as a key only if we are tracking these modifications
				if pep_data.mods != 'Unmodified' and self.isTrackedModification(pep_data.mods):
					pep_key = flds['modified.sequence']
				# Retrieve the relevant PeptideQuant entry (the one with this sequence as its key), or make a new one
				try:
					pep_entry = self.peptide_data[pep_key]
				except KeyError:
					pep_entry = PeptideQuant(pep_key)
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
							prot_entry = ProteinQuant(orf)
							self.protein_data[orf] = prot_entry
						prot_entry.add(pep_entry)
				parsed = True
		return parsed

	def readData(self, orf_dict):
		inf = file(os.path.expanduser(self.filename),'r')
		dlr = util.DelimitedLineReader(inf, strip=False, header_name_processor=util.maxQuantHeader)
		# Read in the data
		max_lines = 1e7
		line = 0
		while not dlr.atEnd() and line < max_lines:
			line += 1
			flds = dlr.nextDict()
			self.parseFields(flds, orf_dict)

	def getPeptideKeys(self):
		return self.peptide_data.keys()

	def getProteinKeys(self):
		return self.protein_data.keys()

	def normalizeRatiosBy(self, norm_prot):
		ratio_stats = norm_prot.getHeavyLightRatioSummary()
		ratio_norm_stats = norm_prot.getNormalizedHeavyLightRatioSummary()
		print math.exp(ratio_stats.mean), math.exp(ratio_norm_stats.mean)
		#print ratio_stats.mean, ratio_norm_stats.mean
		for pep in self.peptide_data.values():
			pep.normalizeRatiosBy(math.exp(ratio_stats.mean), math.exp(ratio_norm_stats.mean))

	def normalizePeptideRatios(self, norm_prot, target_log_mean=0.0):
		# DAD: implement
		target_peps = [self.peptide_data[k] for k in peptide_keys]
		acc = stats.Accumulator(store=False)
		for p in target_peps:
			for ratio in p.heavy_light_ratio_list:
				acc.add(math.log(ratio))
		log_mean = acc.getMean()
		#log_median = acc.getMedian()
		log_sd = 1.0 #acc.getSD()
		#print math.exp(log_mean), log_sd, acc.getVariance(), acc.getN()
		(mean, sd) = (math.exp(log_mean), math.exp(log_sd))
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
		l_acc = stats.Accumulator(store=True)
		for p in target_peps:
			for intens in p.intensity_h_list:
				if not util.isNA(intens):
					h_acc.add(intens)
			for intens in p.intensity_l_list:
				if not util.isNA(intens):
					l_acc.add(intens)
		median_h = h_acc.getMedian()
		median_l = l_acc.getMedian()
		# Normalize all peptides
		for p_key in self.peptide_data.keys():
			p = self.peptide_data[p_key]
			p.normalizeHeavyIntensity(median_h)
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
		recs_with_ratios = [r for r in self.protein_data.values() if not util.isNA(getattr(r,ratio_field))]
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
			line += "%s\t%s\t%s, " % (self.filename, hl_str, self.experiment)
		line += '%d peptides, %d proteins' % (len(self.peptide_data), len(self.protein_data))
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

	def load(self, evidence_descs, filter_tags, filter_experiments, orf_dict):
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
					exev = ExperimentEvidence()
					exev.initFrom(ed)
					self.experiments.append(exev)
					try:
						evidence_fnames[ed.filename].append(exev)
					except KeyError:
						evidence_fnames[ed.filename] = [exev]
		for fname in evidence_fnames.keys():
			exp_fname = os.path.expanduser(fname)
			assert os.path.isfile(exp_fname), "# No file found for %s\n" % exp_fname
			inf = file(exp_fname,'r')
			dlr = util.DelimitedLineReader(inf, strip=False, header_name_processor=util.maxQuantHeader)
			#print dlr.headers
			# Read in the data
			max_lines = 1e8
			line = 0
			while not dlr.atEnd() and line < max_lines:
				line += 1
				flds = dlr.nextDict()
				# Let all the experiments  try to parse this line
				relevant_experiments = evidence_fnames[exp_fname]
				for ex in relevant_experiments:
					res = ex.parseFields(flds, orf_dict)
			inf.close()
		return self.experiments

class FieldFormatter:
	def __init__(self, var, format, transform=lambda x: x):
		self.var = var
		self.format = format
		self.transform = transform

	def __str__(self):
		res = None
		if not util.isNA(self.var):
			try:
				trans_var = self.transform(self.var)
				res = self.format.format(trans_var)
			except ValueError:
				pass
			except TypeError:
				pass
		else:
			res = 'NA'
		return res

