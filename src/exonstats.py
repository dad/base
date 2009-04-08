import sys, os, math, string, random, pickle
import stats, translate, muscle, cai

# Goal: compute statistics about gene structure from chromosomal data

class CodingExonRecord:
	def __init__( self ):
		self.gene_ID = ''
		self.transcript_ID = ''
		self.peptide_ID = ''
		self.exon_ID = ''
		self.strand = ''
		self.exon_start = 0
		self.exon_end = 0
		self.coding_start = 0
		self.coding_end = 0
		self.chromosome = ''

	def copy( self, it ):
		self.gene_ID = it.gene_ID
		self.transcript_ID = it.transcript_ID
		self.peptide_ID = it.peptide_ID
		self.exon_ID = it.exon_ID
		self.strand = it.strand
		self.exon_start = it.exon_start
		self.exon_end = it.exon_end
		self.coding_start = it.coding_start
		self.coding_end = it.coding_end
		self.chromosome = it.chromosome

	def getKey( self ):
		return "%s-%s" % (self.exon_ID, self.peptide_ID)

	def printAttributes( self ):
		print self

	def printRecord( self ):
		print self

	def __repr__(self):
		return "%s %s %s %s %s %s %i %i %i %i" % ( self.chromosome, self.gene_ID, self.transcript_ID, self.peptide_ID,
												   self.exon_ID, self.strand, self.exon_start, self.exon_end, self.coding_start, self.coding_end )

def readChromosomeMapping(chr_map_filename):
	f = file(chr_map_filename)
	lines = f.readlines()
	chr_map = {}
	for line in lines:
		if line[0] == '#':
			continue
		flds = line.strip().split()
		chr_map[flds[0]] = flds[1]
	f.close()
	return chr_map

def readHeaderMapping(header_filename):
	header_dict = {}
	for line in file(header_filename,'r').readlines():
		if not line[0]=='#':
			flds = line.strip().split('\t')
			std_name = flds[0]
			specific_name = flds[1]
			header_dict[std_name] = specific_name.lower()
	return header_dict

def readExons(coord_filename, header_dict):
	f = file(coord_filename)
	lines = f.readlines()
	# Sort out header information
	header_line = lines[0].strip().lower().split('\t')
	gene_id_ind = header_line.index(header_dict["gene.id"])
	transcript_id_ind = header_line.index(header_dict["transcript.id"])
	peptide_id_ind = header_line.index(header_dict["peptide.id"])
	exon_id_ind = header_line.index(header_dict["exon.id"])
	chromosome_ind = header_line.index(header_dict["chromosome"])
	exon_start_ind = header_line.index(header_dict["exon.start"])
	exon_end_ind = header_line.index(header_dict["exon.end"])
	coding_start_ind = header_line.index(header_dict["coding.start"])
	coding_end_ind = header_line.index(header_dict["coding.end"])
	strand_ind = header_line.index(header_dict["strand"])

	peptide_records = {}
	# Read exon records
	for line in lines[1:]:
		flds = line.split('\t')
		if len(line.strip().split('\t')) < 8:
			continue
		try:
			xrec = CodingExonRecord()
			xrec.gene_ID = flds[gene_id_ind].strip()
			xrec.peptide_ID = flds[peptide_id_ind].strip()
			xrec.transcript_ID = flds[transcript_id_ind].strip()
			xrec.exon_ID = flds[exon_id_ind].strip()
			xrec.chromosome = flds[chromosome_ind].strip()
			xrec.exon_start = int(flds[exon_start_ind])
			xrec.exon_end = int(flds[exon_end_ind])
			xrec.strand = flds[strand_ind].strip()
			try:
				xrec.coding_start = int(flds[coding_start_ind])
				xrec.coding_end = int(flds[coding_end_ind])
				key = xrec.getKey()
				if peptide_records.has_key(xrec.peptide_ID):
					peptide_records[xrec.peptide_ID][key] = xrec
				else:
					peptide_records[xrec.peptide_ID] = {key:xrec}
			except ValueError:
				continue
		except Exception,e:
			print "#", e
			continue

	f.close()

	return peptide_records

def buildCodingSequence(recs, chrseq, exclude_nt_from_boundary):
	seq = ''
	# Coding_start indicates the position in the exon rather than the chromosomal position
	# in this case, r.coding_start == 1 for first exon
	exon_frags = []
	exon_seq = ""
	for r in recs:
		if r.coding_start: # if this is a coding exon
			# frags are in 0-based indexes.
			exon_frags.append((r.exon_start-1, r.exon_end))
			exon_seq += chrseq[r.exon_start-1:r.exon_end]
	# If more than one exon...
	'''
	Example (strand = -1):
	exon               e.start  e.end	 c.s c.e e.l c.l
	ENSDARE00000615257 40356040 40356181 505 576 142  72
	ENSDARE00000615257 40356470 40356704 270 504 235 235
	ENSDARE00000615257 40357179 40357350  98 269 172 172
	ENSDARE00000615257 40361610 40361840   1  97 231  97
	'''
	#if len(exon_frags) > 1:
	# take only the 3' end of the first exon and only the 5' end of the last exon
	first_fragment_length = recs[0].coding_end - recs[0].coding_start + 1
	last_fragment_length = recs[-1].coding_end - recs[-1].coding_start + 1
	(s,e) = exon_frags[0]
	exon_frags[0] = (e - first_fragment_length,e)
	(s,e) = exon_frags[-1]
	exon_frags[-1] = (s, s + last_fragment_length)

	coding_seq = ""
	for (s,e) in exon_frags:
		coding_seq += chrseq[s:e]
	return coding_seq, recs[0]

'''
	if recs[0].peptide_ID == 'ENSDARP00000010075':
		print "*** checking ***"
		for i in range(len(recs)):
			r = recs[i]
			(s,e) = exon_frags[i]
			print r.exon_ID, s, e, r.exon_start, r.exon_end, r.coding_start, r.coding_end, r.exon_end - r.exon_start + 1, r.coding_end - r.coding_start + 1, e-s
			if r.strand == '-1':
				print translate.reverseComplement(chrseq[s:e])
			else:
				print chrseq[s:e]
		print translate.translate(coding_seq)
'''

def buildIntronSequence(recs, chrseq):
	seq = ''
	if len(recs) == 1:
		return seq, recs[0]
	intron_boundaries = []
	for i in range(len(recs)-1):
		(intron_start, intron_end) = (recs[i].exon_end, recs[i+1].exon_start-1)
		intron_boundaries.append((intron_start, intron_end))
	#print intron_boundaries
	for (intron_start, intron_end) in intron_boundaries:
		seq += chrseq[intron_start:intron_end]
	return seq, recs[0]

# Return the amount of coding exon sequence that is outside the given
# number of nucleotides from an exon boundary, along with total coding amount
def codingOutsideBoundary(recs, nt_bound):
	coding_total = 0
	coding_outside = 0
	for r in recs:
		#assert r.exon_end >= r.coding_end
		#assert r.exon_start <= r.coding_start
		exon_coding = r.coding_end - r.coding_start + 1
		exon_total = r.exon_end - r.exon_start + 1
		if exon_total <= 2*nt_bound:
			# No possibility of coding sequence outside, since
			# all of exon is within boundary.
			continue
		else:
			# How much coding sequence is outside?
			# If all exon is coding, answer is exon_coding - 2*nt_bound
			exon_coding_out = exon_coding - 2*nt_bound
			# Only modulator is that some sequence can be noncoding
			# at the ends of the exon
			exon_coding_out += min(nt_bound,(r.exon_end-r.coding_end)) + min(nt_bound,(r.coding_start-r.exon_start))
		if False:
			print r
			print exon_coding_out, exon_coding
		assert exon_coding >= exon_coding_out
		if exon_coding_out>0:
			coding_outside += exon_coding_out
		coding_total += exon_coding
	return coding_outside, coding_total

def writeStats(chromosomes, peptide_records, bounds, chromosome_load_fxn, outfile):
	header = "gene\tpep\ttrans\tchr\tlen.intron\tn.exons\tfrac.coding\tgc.intron\tgcind.intron\tgt.intron\tlen.coding\tfrac.%s\n" % \
				  ('\tfrac.'.join([str(b) for b in bounds],))
	outfile.write(header)
	print header,
	total_genes = 0
	total_errors = 0
	# Reconstruct genes
	for chr in chromosomes:
		print "# Chromosome", chr
		# Load this chromosome's data
		chrseq = chromosome_load_fxn(chr)
		n_genes = 0
		n_errors = 0
		n_wrong_length = 0
		n_bad_translation = 0
		n_different_translation = 0
		for (pepid, recdict) in peptide_records.items():
			recs = recdict.values()
			if recs[0].chromosome == chr:
				recs.sort( key = lambda e: e.exon_start)
				(seq, sentinel_rec) = buildCodingSequence(recs, chrseq, 0)
				(intron_seq, sentinel_rec_intron) = buildIntronSequence(recs, chrseq)
				# Reverse-complement the sequence if it's on the negative strand
				if recs[0].strand == '-1':
					#print "reversing strand"
					seq = translate.reverseComplement(seq)
					intron_seq = translate.reverseComplement(intron_seq)
				n_genes += 1

				# Write out statistics
				if False:
					print 'g', seq
					print 'i', intron_seq
				if len(intron_seq)>0:
					gc_intron = '%1.4f' % cai.getGC(intron_seq)
					intron_length = len(intron_seq)
					frac_coding = '%1.4f' % (len(seq)/(len(seq)+float(intron_length)),)
					gcind_intron = '%1.4f' % cai.getDinucleotideIndex(intron_seq, 'GC')
					gt_intron = '%1.4f' % cai.getContent(intron_seq, 'GT')
				else:
					gc_intron = 'NA'
					intron_length = 0
					gcind_intron = 'NA'
					if len(seq)>0:
						frac_coding = '1.0'
					else:
						frac_coding = 'NA'
				num_coding_exons = len([xr for xr in recs if xr.coding_end > xr.coding_start])
				# Test to ensure agreement between codingOutsideBoundary and buildCodingSequence
				if False:
					(coding_outside, total_coding) = codingOutsideBoundary(recs, 0)
					print len(seq), total_coding, coding_outside
					assert total_coding == coding_outside
					assert coding_outside == len(seq)
				fracs_inside = []
				for bound in bounds:
					(coding_outside, total_coding) = codingOutsideBoundary(recs, bound)
					if total_coding > 0:
						frac_inside_bound = "%1.4f" % (1-float(coding_outside)/total_coding,)
					else:
						frac_inside_bound = "NA"
					fracs_inside.append(frac_inside_bound)
				line = "%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%s\n" % \
					   (sentinel_rec.gene_ID, sentinel_rec.peptide_ID, sentinel_rec.transcript_ID, sentinel_rec.chromosome, intron_length, num_coding_exons, frac_coding, gc_intron, gcind_intron, gt_intron, len(seq), '\t'.join(fracs_inside))
				outfile.write(line)
				#print line,

		outfile.flush()
		print "# Processed %d genes with %d errors" % (n_genes, n_errors)
		print "#     %d length errors" % (n_wrong_length,)
		print "#     %d bad translation errors" % (n_bad_translation,)
		print "#     %d different translation errors" % (n_different_translation,)
		total_errors += n_errors
		total_genes += n_genes
		chrseqlist = None
		chrseq = None
	print "# Processed %d genes with %d errors total" % (total_genes, total_errors)
	outfile.close()

def checkGenes(chromosomes, exon_records, chromosome_load_fxn, cdnas, id_fxn):
	total_genes = 0
	total_errors = 0
	# Reconstruct genes
	for chr in chromosomes:
		print "# Chromosome", chr
		# Load this chromosome's data
		chrseq = chromosome_load_fxn(chr)
		n_genes = 0
		n_errors = 0
		n_wrong_length = 0
		n_bad_translation = 0
		n_different_translation = 0
		n_no_comparison = 0
		for (pepid, recdict) in exon_records.items():
			recs = recdict.values()
			if recs[0].chromosome == chr:
				recs.sort( key = lambda e: e.exon_start)
				(seq, sentinel_rec) = buildCodingSequence(recs, chrseq, 0)
				# Reverse-complement the sequence if it's on the negative strand
				if recs[0].strand == '-1':
					#print "reversing strand"
					seq = translate.reverseComplement(seq)
				n_genes += 1
				try:
					orf = id_fxn(sentinel_rec)
					cdna = cdnas[orf]
					#print orf
					#print seq
					#print cdna
					if len(seq) == len(cdna):
						if seq != cdna:
							n_errors += 1
					elif len(seq)-3 == len(cdna):
						if seq[0:-3] != cdna:
							n_errors += 1
					else:
						n_errors += 1
						n_wrong_length += 1
				except KeyError:
					n_no_comparison += 1
					continue
		print "# Processed %d genes with %d errors" % (n_genes, n_errors)
		print "#     %d length errors" % (n_wrong_length,)
		print "#     %d no comparison errors" % (n_no_comparison,)
		print "#     %d bad translation errors" % (n_bad_translation,)
		print "#     %d different translation errors" % (n_different_translation,)
		total_errors += n_errors
		total_genes += n_genes
	print "# Processed %d genes with %d errors total" % (total_genes, total_errors)

def writeGenes(chromosomes, exon_records, chromosome_load_fxn, id_fxn, exclude_nt_from_boundary, outfile):
	total_genes = 0
	total_errors = 0
	# Reconstruct genes
	for chr in chromosomes:
		print "# Chromosome", chr
		# Load this chromosome's data
		chrseq = chromosome_load_fxn(chr)
		n_genes = 0
		n_errors = 0
		n_wrong_length = 0
		n_bad_translation = 0
		n_different_translation = 0
		n_no_comparison = 0
		for (pepid, recdict) in exon_records.items():
			recs = recdict.values()
			if recs[0].chromosome == chr:
				strand_sign = int(recs[0].strand)
				recs.sort( key = lambda e: e.exon_start)
				(seq, sentinel_rec) = buildCodingSequence(recs, chrseq, exclude_nt_from_boundary)
				rseq = seq
				# Reverse-complement the sequence if it's on the negative strand
				if recs[0].strand == '-1':
					#print "reversing strand"
					seq = translate.reverseComplement(seq)
				n_genes += 1
				if len(seq) % 3 != 0:
					n_errors += 1
					n_wrong_length += 1
					continue
				if exclude_nt_from_boundary == 0:
					#if sentinel_rec.peptide_ID == 'ENSDARP00000076309':
					#	print seq
					#	print rseq
					#	print translate.translateRaw(seq)
					#	print translate.translateRaw(rseq)
					#	sys.exit()
					try:
						prot = translate.translate(seq)
						if False: #not prot:
							print "****"
							print id_fxn(sentinel_rec)
							print seq
							print rseq
							for rec in recs:
								print rec
							#print "^%s^\t%s" % (recs[0].strand, seq[0:])
							n_errors += 1
							n_bad_translation += 1
							continue
					except translate.BioUtilsError:
						continue
				# Translation is good... write it.
				line = ">%s\n%s\n" % (id_fxn(sentinel_rec), seq)
				outfile.write(line)
				#outfile.write(">%s\n%s\n" % (id_fxn(sentinel_rec), seq))

		print "# Processed %d genes with %d errors" % (n_genes, n_errors)
		print "#     %d length errors" % (n_wrong_length,)
		print "#     %d bad translation errors" % (n_bad_translation,)
		total_errors += n_errors
		total_genes += n_genes
	print "# Wrote %d genes with %d errors total" % (total_genes, total_errors)
