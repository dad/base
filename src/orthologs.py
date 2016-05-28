#! /usr/bin/python

import sys, os, math, string, random
sys.path.append(os.path.abspath('../'))
import cai, translate, my_paml, muscle

def makeProteinFileFromGeneFile(binPath, geneFileName):
	geneDict = translate.Read_FASTA_Dict(geneFileName)
	protFileName = "tmp%d.txt" % random.randint(0,100000)
	protFile = file(protFileName, 'w')
	protDict = {}
	for (gene, seq) in geneDict.items():
		prot = translate.Translate(seq)
		if prot:
			protFile.write(">%s\n%s\n" % (gene, prot))
			protDict[gene] = prot
	protFile.close()
	return protFileName, protDict

def setUpBLASTDataFiles(binPath, queryGeneFileName, targetGeneFileName, queryProtFileName, targetProtFileName, queryProtDict, targetProtDict):
	# We want to BLAST protein sequences, but align gene sequences.  Given gene sequences,
	if not queryProtFileName:
		(queryProtFileName, queryProtDict) = makeProteinFileFromGeneFile(binPath, queryGeneFileName)
	else:
		queryProtDict = translate.Read_FASTA_Dict(queryProtFileName)
	#if not os.path.isfile("%s.psq" % queryProtFileName):
	#	# Now format the query database
	#	formatCmd = "%s/formatdb -i %s -p T -o T" % (binPath, queryProtFileName)
	#	os.popen(formatCmd, 'r')
	if not targetProtFileName:
		(targetProtFileName, targetProtDict) = makeProteinFileFromGeneFile(binPath, targetGeneFileName)
	else:
		targetProtDict = translate.Read_FASTA_Dict(targetProtFileName)
	if not os.path.isfile("%s.psq" % targetProtFileName):
		# Now format the target database
		formatCmd = "%s/formatdb -i %s -p T -o T" % (binPath, targetProtFileName)
		os.popen(formatCmd, 'r')
	return queryProtFileName, targetProtFileName, queryProtDict, targetProtDict

def getBLASTHits(binPath, queryGeneFileName, targetGeneFileName, evalue = 1e-20, queryProtFileName=None, \
	targetProtFileName=None, queryProtDict=None, targetProtDict=None, allowSelfHits = 0, blastOutputFileName=None, readExistingResults=True):
	# Get best protein-protein BLAST hits
	# First, generate protein files from gene files
	(queryProtFileName, targetProtFileName, queryProtDict, targetProtDict) = setUpBLASTDataFiles(\
		binPath, queryGeneFileName, targetGeneFileName, queryProtFileName, targetProtFileName, queryProtDict, targetProtDict)
	cmd = None
	if blastOutputFileName:
		if readExistingResults:
			if os.path.isfile(blastOutputFileName):
				f = file(blastOutputFileName, 'r')
				print "# Reading existing BLAST output from %s" % blastOutputFileName
			else:
				raise IOError, "Previous BLAST results unreadable from %s" % blastOutputFileName
		else:
			cmd = "blastall -p blastp -i %s -d %s -m 8 -e %E -o %s" % (queryProtFileName, targetProtFileName, evalue, blastOutputFileName)
			os.spawnv(os.P_WAIT, os.path.abspath(os.path.join(binPath,"blastall")), [x for x in cmd.split()])
	else:
		blastOutputFileName = "tmpBLAST%d.txt" % random.randint(0,100000)
		cmd = "blastall -p blastp -i %s -d %s -m 8 -e %E -o %s" % (queryProtFileName, targetProtFileName, evalue, blastOutputFileName)
		os.spawnv(os.P_WAIT, os.path.abspath(os.path.join(binPath,"blastall")), [x for x in cmd.split()])
	if cmd:
		print "# %s" % cmd
	hits = {}
	curGene = ""
	f = file(blastOutputFileName, 'r')
	for line in f.readlines():
		flds = line.split('\t')
		curGene = flds[0]
		hitGene = flds[1]
		if curGene != hitGene or allowSelfHits:
			if hits.has_key(curGene):
				hits[curGene].append(hitGene)
			else:
				hits[curGene] = [hitGene]
	f.close()
	return hits, queryProtFileName, targetProtFileName, queryProtDict, targetProtDict
	
def getReciprocalShortestDistancePairs(binPath, queryGeneFileName, targetGeneFileName,  evalue=1e-20, queryProtFileName=None,\
									   targetProtFileName=None, allowSelfHits=0, blastOutputFileName=None):
	queryGeneDict = translate.Read_FASTA_Dict(queryGeneFileName)
	targetGeneDict = translate.Read_FASTA_Dict(targetGeneFileName)

	queryProtDict = {}
	targetProtDict = {}
	distanceCache = {}
	alignmentCache = {}

	# Get query-to-target hits, then filter for shortest-distance hits
	(hits, queryProtFileName, targetProtFileName, queryProtDict, targetProtDict) = \
		   getBLASTHits(binPath, queryGeneFileName, targetGeneFileName, evalue, \
						queryProtFileName, targetProtFileName, queryProtDict, targetProtDict, allowSelfHits, blastOutputFileName)
	queryToTargetSDHits = getShortestDistanceHits(hits, queryGeneDict, targetGeneDict, \
												  queryProtDict, targetProtDict, distanceCache, alignmentCache)

	# Get target-to-query hits, then filter for shortest-distance hits	
	(hits, targetProtFileName, queryProtFileName, targetProtDict, queryProtDict) = \
		   getBLASTHits(binPath, targetGeneFileName, queryGeneFileName, evalue, \
						targetProtFileName, queryProtFileName, targetProtDict, queryProtDict, allowSelfHits, blastOutputFileName)
	targetToQuerySDHits = getShortestDistanceHits(hits, targetGeneDict, queryGeneDict, \
												  targetProtDict, queryProtDict, distanceCache, alignmentCache)

	pairs = []
	for gene in queryToTargetSDHits.keys():
		try:
			if gene == targetToQuerySDHits[queryToTargetSDHits[gene]]:
				# This is a reciprocal shortest distance hit
				print "%s\t%s" % (gene, queryToTargetSDHits[gene])
				pairs.append((gene, queryToTargetSDHits[gene]))
		except KeyError:
			continue
	return pairs, queryGeneDict, targetGeneDict, queryProtDict, targetProtDict, distanceCache, alignmentCache

def cacheKeySeparator():
	return "@@"

def cacheKey(geneA, geneB):
	# Make a key such that key(geneA, geneB) == key(geneB, geneA)
	key = ''
	if geneA < geneB:
		key = "%s%s%s" % (geneA, cacheKeySeparator(), geneB)
	else:
		key = "%s%s%s" % (geneB, cacheKeySeparator(), geneA)
	return key

def getShortestDistanceHits(hits, queryGeneDict, targetGeneDict, queryProtDict, targetProtDict, distanceCache, alignmentCache):
	queryToTargetSDHits = {}
	totalHits = len(hits.keys())
	nHits = 0
	for (queryGene, qHits) in hits.items():
		minDist = 100.0
		minHitGene = None
		for targetGene in qHits:
			key = cacheKey(queryGene, targetGene)
			dNML = 0.0
			try:
				(dDML, dSML, dNML, dDNG, dSNG, dNNG, numSynonymousSites, numNonsynonymousSites, fracAligned, seqIdentity) = distanceCache[key]
			except KeyError:
				queryGeneSeq = queryGeneDict[queryGene]
				targetGeneSeq = targetGeneDict[targetGene]
				[alignedQueryProt, alignedTargetProt] = muscle.align_sequences([queryProtDict[queryGene], targetProtDict[targetGene]])
				alignedQueryGene = muscle.align_gene_from_protein(queryGeneSeq, alignedQueryProt)
				alignedTargetGene = muscle.align_gene_from_protein(targetGeneSeq, alignedTargetProt)
				(dDML, dSML, dNML, dDNG, dSNG, dNNG, numSynonymousSites, numNonsynonymousSites) = my_paml.Get_Distance_NS(alignedQueryGene, alignedTargetGene, 'codon')
				(seqIdentity, numIdentical, numAligned) = sequenceIdentity(alignedQueryProt, alignedTargetProt)
				fracAligned = numAligned/float(len(queryProtDict[queryGene]))
				if fracAligned == 1.0 and seqIdentity == 1.0:
					dNML = 0.0 # Obviously no nonsyn. changes if proteins are identical
				elif 1.0/numNonsynonymousSites > dNML:
					dNML = 1.0/numNonsynonymousSites # Minimum possible change
				(ntseqIdentity, ntnumIdentical, ntnumAligned) = sequenceIdentity(alignedQueryGene, alignedTargetGene)
				ntfracAligned = ntnumAligned/float(len(queryGeneSeq))
				if ntfracAligned == 1.0 and ntseqIdentity == 1.0:
					dSML = 0.0 # Obviously no syn. changes if genes are identical
				elif 1.0/numSynonymousSites > dSML:
					dSML = 1.0/numSynonymousSites # Minimum possible change
				# Cache this distance and alignment
				distanceCache[key] = (dDML, dSML, dNML, dDNG, dSNG, dNNG, numSynonymousSites, numNonsynonymousSites, fracAligned, seqIdentity)
				alignmentCache[key] = ((queryGene, alignedQueryGene, alignedQueryProt), (targetGene, alignedTargetGene, alignedTargetProt))
			# Check to see if this is a shorter distance
			if dNML < minDist and targetGene != queryGene:
				minDist = dNML
				minHitGene = targetGene
		nHits += 1
		(dDML, dSML, dNML, dDNG, dSNG, dNNG, numSynonymousSites, numNonsynonymousSites, fracAligned, seqIdentity) = distanceCache[cacheKey(queryGene, minHitGene)]
		print "# %d of %d: %s was %s, dN = %1.6f, fracAlign = %1.6f" % (nHits, totalHits, queryGene, minHitGene, minDist, fracAligned)
		queryToTargetSDHits[queryGene] = minHitGene
	return queryToTargetSDHits

def sequenceIdentity(alignedSeq1, alignedSeq2):
	numIdentical = 0
	numAligned = 0
	for i in range(len(alignedSeq1)):
		aa1 = alignedSeq1[i]
		aa2 = alignedSeq2[i]
		if aa1 != '-' and aa2 != '-':
			numAligned += 1
			if aa1 == aa2:
				numIdentical += 1
	seqID = 0.0
	if numAligned > 0:
		seqID = float(numIdentical)/numAligned 
	return seqID, numIdentical, numAligned
	
