"""Annotate all windows with functional classes.

The annotation is saved regardless of the percentage of overlap 
between a given functional class and a window.

In order to annotate the windows, it first parses the TSV files 
downloaded from *UCSC Table Browser*.

- **Use**::
	
	python3 annotateWindows.py

- **Example of Usage**::

	python3 ~/code/utils/annotateWindows.py

Pre-requisites
--------------

Before using this script, make sure all the required files were pre-computed:

a) Download .TSV files from UCSC Table browser
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Make sure you have the TSV files containing the annotation for these 4 tracks: 

  1. **CCDS**: Transcript coding sequence (CDS) genomic regions;
  2. **GENCODE V44**;
  3. **RepeatMasker**;
  4. **RefSeq**: Functional elements from RefSeq.

Alternatively, these 4 files can also be downloaded (as ``.tar.xz`` compressed files due to their sizes) 
from this repository, in the directory scripts/assets/annotation.

b) Windows for each chromosome in the dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Make sure to run ``2_computeWindows.py`` for **every** chromosome included in 
the attribute ``chromLst`` declared in file ``dataset.py``. 


Time, Memory & Disk space
-------------------------

Running the script on **80 cores** takes **1 hour 32 minutes** and requires **15 GB** memory. 

=====================  ========
Step                   Time (s)
=====================  ========
Load all annotations      19.90
Annotate windows        5493.07
**Total time**          5512.96
=====================  ========

**Output files**:
	
	The script creates the file ``hg38.1000.annotations.pickle``, which takes **127 MB** of space.

"""
import os
import pickle
import collections
from collections import namedtuple
import sys
import numpy as np
from concurrent import futures
import argparse

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils.dataset import Dataset
from utils.basicTypes import Time

####################################
# Parsers.

# CCDS.tsv
# The Consensus Coding Sequence Project is a list of transcript 
# coding sequence (CDS) genomic regions that are identically 
# annotated by RefSeq and Ensembl/GENCODE. 
# CCDS undergoes extensive manual review and you can consider 
# these a subset of either gene track, filtered for high quality. 
# As the name implies, it does not cover UTR regions or non-coding transcripts. 
def loadCCDS(filename, qChromLst):
	source="CCDS"
	annotations = {}
	for qChrom in qChromLst:
		annotations[qChrom] = []
	with open(filename) as f:
		for line in f:
			if (not line.startswith("#")):
				elems = line.split("\t")
				# CCDS30745.1	chr1	-	67093004	67127240	67093004	67127240	6	67093004,67095234,67096251,67115351,67125751,67127165,	67093604,67095421,67096321,67115464,67125909,67127240,
				if(len(elems) == 10):
					elemId = elems[0]
					qChrom = elems[1]
					strand = elems[2]
					begPos = int(elems[3])
					endPos = int(elems[4])

					begPosCDS = int(elems[5])
					endPosCDS = int(elems[6])

					nbExons   = int(elems[7])

					# Sanity check.
					if((begPos != begPosCDS) or (endPos != endPosCDS)):
						print(f"WARNING! Inconsistent positions: ({begPos},{endPos}) and ({begPosCDS},{endPosCDS}).")

					if (qChrom in qChromLst):
						annotations[qChrom].append((begPos, endPos, strand))
						if(begPos > endPos):
							print(f"WARNING! {qChrom} : Inverted positions={(begPos, endPos, strand)}")
	return annotations

# (posBeg, posEnd, transcriptClass, transcriptBiotype, geneInfo, id, source)
# GENCODE-V44-all.tsv
def loadGencode(filename, qChromLst, annotationsCCDS):
	source="GENCODE-V44"
	annotations       = {}
	annotationsCoding = {}

	for qChrom in qChromLst:
		annotations[qChrom] = []
		annotationsCoding[qChrom] = []

	nbCoding = 0
	nbCodingNonCCDS = 0

	with open(filename) as f:
		for line in f:
			if (not line.startswith("#")):
				elems = line.split("\t")
				# ENST00000473358.1	chr1	+	29553	31097	29553	29553	3	29553,30563,30975,	30039,30667,31097,	MIR1302-2HG	lncRNA	lncRNA	nonCoding	n/a
				if(len(elems) == 15):
					elemId = elems[0]
					qChrom = elems[1]
					strand = elems[2]
					begPos = int(elems[3])
					endPos = int(elems[4])

					begPosCDS = int(elems[5])
					endPosCDS = int(elems[6])

					nbExons    = int(elems[7])

					exonsStart = (elems[8].split(","))[:-1]
					exonsStart = [int(exonPos) for exonPos in exonsStart]

					exonsEnd   = (elems[9].split(","))[:-1]
					exonsEnd   = [int(exonPos) for exonPos in exonsEnd]

					geneType          = elems[10]
					transcriptBiotype = elems[11]
					transcriptDetails = elems[12]
					transcriptClass   = elems[13]

					proteinClassPanther = elems[14].strip() # remove newline char

					if(transcriptClass == "problem"):
						transcriptBiotype = transcriptBiotype + "_" + transcriptDetails

					if (qChrom in qChromLst):
						if(transcriptClass == "coding"):

							# Check if it exists in CCDS list.
							if ((begPosCDS, endPosCDS, strand) in annotationsCCDS[qChrom]):
								# Compute UTRs (3' and 5').
								begPosUTR = begPos
								endPosUTR = begPosCDS-1
								if((endPosUTR-begPosUTR) > 0):
									annotationsCoding[qChrom].append((begPosUTR, endPosUTR, strand, transcriptClass, "utr", elemId))
							
								begPosUTR = endPosCDS
								endPosUTR = endPos
								if((endPosUTR-begPosUTR) > 0):
									annotationsCoding[qChrom].append((begPosUTR, endPosUTR, strand, transcriptClass, "utr", elemId))

								# Compute introns/exons.
								for idxExon, (exonStart, exonEnd) in enumerate(zip(exonsStart, exonsEnd)):

									# Add intron.
									if(idxExon > 0):
										exonEndPrev = exonsEnd[idxExon-1]
										annotationsCoding[qChrom].append((exonEndPrev, exonStart-1, strand, transcriptClass, "intron", elemId))

									# Add exon.
									annotationsCoding[qChrom].append((exonStart, exonEnd, strand, transcriptClass, geneType, elemId))

							# Add ignored coding region in the GENCODE mapping, but not in CCDS mapping.
							else:

								# Compute UTRs (3' and 5').
								begPosUTR = begPos
								endPosUTR = begPosCDS-1
								if((endPosUTR-begPosUTR) > 0):
									annotations[qChrom].append((begPosUTR, endPosUTR, strand, transcriptClass, "utr", elemId))
							
								begPosUTR = endPosCDS
								endPosUTR = endPos
								if((endPosUTR-begPosUTR) > 0):
									annotations[qChrom].append((begPosUTR, endPosUTR, strand, transcriptClass, "utr", elemId))

								# Compute introns/exons.
								for idxExon, (exonStart, exonEnd) in enumerate(zip(exonsStart, exonsEnd)):

									# Add intron.
									if(idxExon > 0):
										exonEndPrev = exonsEnd[idxExon-1]
										annotations[qChrom].append((exonEndPrev, exonStart-1, strand, transcriptClass, "intron", elemId))

									# Add exon.
									annotations[qChrom].append((exonStart, exonEnd, strand, transcriptClass, geneType, elemId))

								nbCodingNonCCDS += 1
							nbCoding += 1
							
						else:

							if(nbExons == 1):
								annotations[qChrom].append((begPos, endPos, strand, transcriptClass, transcriptBiotype, elemId))
							else:
								for idxExon, (exonStart, exonEnd) in enumerate(zip(exonsStart, exonsEnd)):
									# Add intron.
									if(idxExon > 0):
										exonEndPrev = exonsEnd[idxExon-1]
										annotations[qChrom].append((exonEndPrev, exonStart-1, strand, transcriptClass, "intron", elemId))
									# Add exon (probably from lncRNA).
									annotations[qChrom].append((exonStart, exonEnd, strand, transcriptClass,  transcriptBiotype, elemId))

						if(begPos > endPos):
							print(f"{qChrom} : {(begPos, endPos, strand, transcriptClass, transcriptBiotype, elemId)}")

	print(f" - Total coding = {nbCoding}; Total coding NOT IN CCDS = {nbCodingNonCCDS}")
	return annotations, annotationsCoding

# RepeatMasker
# https://genome.ucsc.edu/cgi-bin/hgTrackUi?g=rmsk
# (posBeg, posEnd, repClass, repFamily, extraDet (None), id, source)
# Repeats classes (repClass)
# In full display mode, this track displays up to ten different classes of repeats:
#     Short interspersed nuclear elements (SINE), which include ALUs
#     Long interspersed nuclear elements (LINE)
#     Long terminal repeat elements (LTR), which include retroposons
#     DNA repeat elements (DNA)
#     Simple repeats (micro-satellites)
#     Low complexity repeats
#     Satellite repeats
#     RNA repeats (including RNA, tRNA, rRNA, snRNA, scRNA, srpRNA)
#     Other repeats, which includes class RC (Rolling Circle)
#     Unknown
def loadRepeats(filename, qChromLst):
	source="RepeatMasker"
	annotations = {}
	for qChrom in qChromLst:
		annotations[qChrom] = []

	with open(filename) as f:
		for line in f:
			if (not line.startswith("#")):
				elems = line.split("\t")
				# genoName	genoStart	genoEnd	strand	repName	repClass	repFamily
				# chr1	67108753	67109046	+	L1P5	LINE	L1
				if(len(elems) == 7):
					qChrom = elems[0]
					begPos = int(elems[1])
					endPos = int(elems[2])
					strand = elems[3]
					elemId = elems[4]
					repClass  = elems[5]
					repFamily = elems[6].strip() # remove newline char

					if (repClass.endswith("?") or repFamily.endswith("?")):
						repFamily = repClass + "_" + repFamily
						repClass  = "unknown"

					if (qChrom in qChromLst):
						annotations[qChrom].append((begPos, endPos, strand, repClass, repFamily, elemId))
						if(begPos > endPos):
							print(f"{qChrom} : {(begPos, endPos, strand, repClass, repFamily, elemId)}")
	return annotations

# Functional elements from RefSeq
# https://genome.ucsc.edu/cgi-bin/hgTables
#  Sequence Ontology (SO) term [soTerm]
#    - Regulatory elements (items labeled by INSDC regulatory class)
#    - Protein binding sites (items labeled by bound moiety)
#    - Mobile elements
#    - Recombination features
#    - Sequence features 
#    - Other
def loadFunctionalElements(filename, qChromLst):
	source="RepeatMasker"
	annotations = {}
	for qChrom in qChromLst:
		annotations[qChrom] = []

	with open(filename) as f:
		for line in f:
			if (not line.startswith("#")):
				elems = line.split("\t")
				# chrom	chromStart	chromEnd	name	strand	soTerm
				# chr1	5742550	5743135	DNase_I_hypersensitive_site	+	region
				if(len(elems) == 6):
					qChrom = elems[0]
					begPos = int(elems[1])
					endPos = int(elems[2])
					elemId = elems[3]
					strand = elems[4]
					soTerm = elems[5].strip() # remove newline char

					if (qChrom in qChromLst):
						annotations[qChrom].append((begPos, endPos, strand, soTerm, elemId, elemId))
						if(begPos > endPos):
							print(f"{qChrom} : {(begPos, endPos, strand, soTerm, elemId, elemId)}")
	return annotations

# WARNING: It ended up not being actually used in the main text (maybe in a future paper?).
def loadProteinPantherAnnotation(filename):
	source="ProteinPanther"
	annotations = {}
	with open(filename) as f:
		for line in f:
			if (not line.startswith("#")):
				elems = line.split("\t")
				# ENST00000641515.2	chr1	+	65418	71585	65564	70008	OR4F5	protein_coding	protein_coding	coding	G-protein coupled receptor
				if(len(elems) == 12):
					elemId = elems[0]
					qChrom = elems[1]
					strand = elems[2]
					begPos = int(elems[3])
					endPos = int(elems[4])

					begPosCDS = int(elems[5])
					endPosCDS = int(elems[6])

					geneType = elems[7]
					proteinPantherClass = elems[11].strip() # remove newline char
					
					annotations[geneType] = proteinPantherClass
	return annotations

def loadAllAnnotations(my_dataset):
	""" Parses all annotation files download from the UCSC Table Browser.
		It returns a dictionary where all annotations are grouped by chromosome.
	"""

	# Load annotation files.
	filenameCCDS    = os.path.join(my_dataset.dirAnnot, f"CCDS-all.tsv")
	annotationsCCDS = loadCCDS(filenameCCDS, my_dataset.chromLst)
	# 63,85% of "coding" annotations from GENCODE (87,792 out of 137,517) are not present in CCDS.
	print(f" - CCDS={sum([len(annotationsCCDS[qChrom]) for qChrom in annotationsCCDS.keys()])} annotations.")

	filenameGencode    = os.path.join(my_dataset.dirAnnot, f"GENCODE-V44-all.tsv")
	annotationsGencode, annotationsCCDS = loadGencode(filenameGencode, my_dataset.chromLst, annotationsCCDS)

	filenameRepeatMasker    = os.path.join(my_dataset.dirAnnot, f"RepeatMasker-all.tsv")
	annotationsRepeatMasker = loadRepeats(filenameRepeatMasker, my_dataset.chromLst)

	filenameRefSeqFuncElems    = os.path.join(my_dataset.dirAnnot, f"RefSeqFuncElems-all.tsv")
	annotationsRefSeqFuncElems = loadFunctionalElements(filenameRefSeqFuncElems, my_dataset.chromLst)

	# Print info about annotations (CCDS).
	for qChrom in my_dataset.chromLst:
		groupedInfo = {}
		for elem in annotationsCCDS[qChrom]:
			begPos, endPos, strand, transcriptClass, transcriptBiotype, elemId = elem
			if(transcriptClass not in groupedInfo):
				groupedInfo[transcriptClass] = {}

			transcriptBiotype = transcriptBiotype if(transcriptBiotype in ["utr","intron"]) else "other"
			if(transcriptBiotype not in groupedInfo[transcriptClass]):
				groupedInfo[transcriptClass][transcriptBiotype] = []
			groupedInfo[transcriptClass][transcriptBiotype].append(elem)

		print(f"{qChrom}")
		for transcriptClass in groupedInfo:
			totalItems = sum([len(groupedInfo[transcriptClass][transcriptBiotype]) 
							for transcriptBiotype in groupedInfo[transcriptClass]])
			print(f"\t{transcriptClass} = {totalItems} annotations.")
			for transcriptBiotype in groupedInfo[transcriptClass]:
				print(f"\t\t {transcriptBiotype} = {len(groupedInfo[transcriptClass][transcriptBiotype])} items.")

	# Print info about annotations (Gencode).
	for qChrom in my_dataset.chromLst:
		
		groupedInfo = {}
		for elem in annotationsGencode[qChrom]:
			begPos, endPos, strand, transcriptClass, transcriptBiotype, elemId = elem
			if(transcriptClass not in groupedInfo):
				groupedInfo[transcriptClass] = {}
			if(transcriptBiotype not in groupedInfo[transcriptClass]):
				groupedInfo[transcriptClass][transcriptBiotype] = []
			groupedInfo[transcriptClass][transcriptBiotype].append(elem)

		print(f"{qChrom}")
		for transcriptClass in groupedInfo:
			totalItems = sum([len(groupedInfo[transcriptClass][transcriptBiotype]) 
							for transcriptBiotype in groupedInfo[transcriptClass]])
			print(f"\t{transcriptClass} = {totalItems} annotations.")
			for transcriptBiotype in groupedInfo[transcriptClass]:
				print(f"\t\t {transcriptBiotype} = {len(groupedInfo[transcriptClass][transcriptBiotype])} items.")

	# Print info about annotations (RepeatMasker).
	for qChrom in my_dataset.chromLst:
		
		groupedInfo = {}
		for elem in annotationsRepeatMasker[qChrom]:
			begPos, endPos, strand, repClass, repFamily, elemId = elem
			if(repClass not in groupedInfo):
				groupedInfo[repClass] = {}
			if(repFamily not in groupedInfo[repClass]):
				groupedInfo[repClass][repFamily] = []
			groupedInfo[repClass][repFamily].append(elem)

		print(f"{qChrom}")
		for repClass in groupedInfo:
			totalItems = sum([len(groupedInfo[repClass][repFamily]) 
							for repFamily in groupedInfo[repClass]])
			print(f"\t{repClass} = {totalItems} annotations.")
			for repFamily in groupedInfo[repClass]:
				print(f"\t\t {repFamily} = {len(groupedInfo[repClass][repFamily])} items.")

	# Print info about annotations (RefSeq Functional Elements).
	for qChrom in my_dataset.chromLst:
		
		groupedInfo = {}
		for elem in annotationsRefSeqFuncElems[qChrom]:
			begPos, endPos, strand, soTerm, elemId, elemId = elem
			if(soTerm not in groupedInfo):
				groupedInfo[soTerm] = {}
			if(elemId not in groupedInfo[soTerm]):
				groupedInfo[soTerm][elemId] = []
			groupedInfo[soTerm][elemId].append(elem)

		print(f"{qChrom}")
		for soTerm in groupedInfo:
			totalItems = sum([len(groupedInfo[soTerm][elemId]) 
							for elemId in groupedInfo[soTerm]])
			print(f"\t{soTerm} = {totalItems} annotations.")
			for elemId in groupedInfo[soTerm]:
				print(f"\t\t {elemId} = {len(groupedInfo[soTerm][elemId])} items.")

	allAnnotations = [annotationsGencode, annotationsCCDS, annotationsRepeatMasker, annotationsRefSeqFuncElems]
	return allAnnotations

###############################
# Methods to load info from our study.

def computeEmptyWindows(my_dataset, qChrom, windowsIds):
	# Compute the region in the chromosome covered by windows with PCSs.
	chromSize  = my_dataset.chromSizes[qChrom]
	windowsIds = sorted(windowsIds,key=lambda windowId: windowId[0])
	minBegPos, minEndPos = windowsIds[0]
	maxBegPos, maxEndPos = windowsIds[-1]
	# Add empty windows.
	emptyWindows = []
	windowSize   = my_dataset.windowSize
	# Empty windows between the start of the chromosome and the first occurrence of a PCS.
	endPos = minBegPos
	begPos = endPos-windowSize
	while(begPos > 0):
		emptyWindows.append((begPos,endPos))
		endPos = begPos
		begPos = endPos-windowSize
	# Empty windows between the last occurrence of a PCS and the end of the chromosome.
	begPos = maxEndPos
	endPos = begPos+windowSize
	while(endPos < chromSize):
		emptyWindows.append((begPos,endPos))
		begPos = endPos
		endPos = begPos+windowSize
	return emptyWindows

def loadWindowIds(my_dataset):
	""" This function loads all computed windows and returns the window ids grouped by chromosome.
	"""
	windowsIds = {}
	for qChrom in my_dataset.chromLst:
		print(f"\t[{qChrom}] Loading windows.")
		# Check if input file exists.
		windows_filename = my_dataset.getOutFilename_computeWindows(qChrom)
		if (not os.path.isfile(windows_filename)):
			print(f"ERROR! File not found: {windows_filename}")
			sys.exit()
		# Load windows.
		windows      = pickle.load(open(windows_filename, 'rb'))
		# Sanity check: checks if all species have the same window ID list.
		# If not: other scripts must be checked.
		UCSCnames    = list(windows.keys())
		rdmSpecies   = UCSCnames[0]
		nbWindowsRef = len(windows[rdmSpecies])
		for UCSCname in UCSCnames:
			windowsPerSpecies = windows[UCSCname]
			if(len(windowsPerSpecies) != nbWindowsRef):
				print(f"ERROR! Potential problem: species have different windows: {UCSCname} (Exp.: {nbWindowsRef}; Obs.: {len(windowsPerSpecies)})")
				sys.exit()
		# If it is ok, get windows from any species.
		windows = windows[rdmSpecies]
		windowsIdsPerChr = [windowId for (windowId, PCSsizeDistrib) in windows.items()]
		# Compute empty windows.
		windowsIdsPerChr.extend(computeEmptyWindows(my_dataset, qChrom, windowsIdsPerChr))
		windowsIds[qChrom] = sorted(windowsIdsPerChr,key=lambda windowId: windowId[0])
	return windowsIds

####################################################################
# Annotation related methods.

def getAnnotationsPerChr(allAnnotations, qChrom):
	annotationsGencode, annotationsCCDS, annotationsRepeatMasker, annotationsRefSeqFuncElems = allAnnotations
	return {"Gencode-v4"      : sorted(annotationsGencode[qChrom],         key=lambda a: a[0]),
			"CCDS"            : sorted(annotationsCCDS[qChrom],            key=lambda a: a[0]),
			"RepeatMasker"    : sorted(annotationsRepeatMasker[qChrom],    key=lambda a: a[0]), 
			"RefSeqFuncElems" : sorted(annotationsRefSeqFuncElems[qChrom], key=lambda a: a[0])}

# WARNING! Assume that positions from annotations are sorted by begPos.
def intersectingAnnot(begPosWin, endPosWin, allAnnotations):
	annot = {}
	# Compute intersecting annotations.
	for annotType in allAnnotations:
		annot[annotType] = []
		for (begPos, endPos, strand, mainCat, subCat, elemId) in allAnnotations[annotType]:
			if (((begPos <= begPosWin) and (begPosWin < endPos))
			or  ((begPosWin <= begPos) and (begPos < endPosWin))):
				annot[annotType].append((begPos, endPos, strand, mainCat, subCat, elemId))
			# WARNING: Assume annotations are sorted by begPos.
			if(begPos > endPosWin): break
	return annot

# Format annotation to follow the annotation hierarchy and order from Figure 5.
def getFormattedAnnot(annotDesc,proteinClasses=None,CCDSonly=False):
	source,type,subtype = annotDesc

	label	= source + "_" + type + "_" + subtype
	category = (source,type,subtype)
	order	 = (0,0,0)

	# Special labels.
	if(source == "Gencode-v4"):

		# Genes found in "GENCODE" but not in "CCDS" end up here.
		# coding - display protein coding transcripts, including polymorphic pseudogenes
		if(type == "coding"):
			if(not CCDSonly):
				if(subtype == "utr"):
					order	 = (2,2,2)
				elif(subtype == "intron"):
					order	 = (2,2,1)
				else:
					order	 = (1,1,1)
			else:
				order = (0,0,0)
				
			# https://www.genenames.org/data/genegroup/#!/group/348
			# Immunoglobulins (IG) or antibodies are antigen receptors of the B cells of the adaptive immune response. 
			# elif(subtype == "IG_V_gene"):
			#	 label	= f"Protein Coding - ImmunoGlobulin Variable (IG V) gene"
			#	 category = ("Coding","Protein","ImmunoGlobulin Variable (IG V) gene")
			#	 order	 = (1,1,2)
						
		# nonCoding - display non-protein coding transcripts
		elif(type == "nonCoding"):
			if(subtype == "lncRNA"):
				order	 = (1,2,1)
			elif(subtype == "miRNA"):
				order	 = (1,2,2)
			elif(subtype == "snRNA"):
				order	 = (1,2,3)
			elif(subtype == "snoRNA"):
				order	 = (1,2,4)
			elif(subtype == "intron"):
				order	 = (2,2,1) # "2.1.5"
			else:
				order	 = (1,2,5)

		# problem - display problem transcripts (Biotypes of retained_intron, TEC, or disrupted_domain) 
		elif(type == "problem"):
			order    = (0,0,0)
		# pseudo - display pseudogene transcript annotations
		elif(type == "pseudo"):
			order	 = (1,3,1)
	
	# https://en.wikipedia.org/wiki/List_of_gene_families
	elif(source == "CCDS"):
		if(type == "coding"):
			if(subtype == "utr"):
				order	 = (2,2,2)
			elif(subtype == "intron"):
				order	 = (2,2,1)
			else:
				# proteinDetails = proteinClasses[subtype] if(proteinClasses and (subtype in proteinClasses) and (len(proteinClasses[subtype]) > 5)) else ""
				# proteinDetails = "" if(proteinDetails == "n/a") else proteinDetails
				# label	= f"Protein Coding ({proteinDetails})"
				# category = ("Coding","Protein",proteinDetails)
				# order	 = "1.1." +  str(len(proteinClasses[subtype]))
				order	 = (1,1,1)
				
	# NCBI recently announced a new release of functional regulatory elements. 
	# NCBI is now providing RefSeq and Gene records for non-genic functional elements that 
	# have been described in the literature and are experimentally validated. 
	# Elements in scope include :
	# experimentally-verified gene regulatory regions (e.g., enhancers, silencers, locus control regions), 
	# known structural elements (e.g., insulators, DNase I hypersensitive sites, matrix/scaffold-associated regions), 
	# well-characterized DNA replication origins, and clinically-significant sites of DNA recombination and genomic instability.
	elif(source == "RefSeqFuncElems"):
		if(type == "enhancer"):
			order	 = (2,1,1)
		elif(type == "silencer"):
			order	 = (2,1,2)
		elif(type == "locus_control_region"):
			order	 = (2,1,5)
		# Insulator is a type of cis-regulatory element known as a long-range regulatory element.
		elif(subtype == "insulator"):
			order	 = (2,1,4)
		# Regulatory regions in general, and promoters in particular, tend to be DNase-sensitive.
		elif(subtype == "DNase_I_hypersensitive_site"):
			order	 = (2,1,3)
		else:
			order	 = (0,0,0)

	# https://genome.ucsc.edu/cgi-bin/hgTrackUi?g=rmsk
	elif(source == "RepeatMasker"):
		if(type == "DNA"):
			if(subtype.startswith("TcMar")):
				order	 = (3,3,2)
			elif(subtype.startswith("hAT")):
				order	 = (3,3,3)
			else:
				order	 = (3,3,1)
		elif(type == "LINE"):
			if(subtype == "CR1"):
				order	 = (3,1,2)
			elif(subtype == "L1"):
				order	 = (3,1,3)
			elif(subtype == "L2"):
				order	 = (3,1,4)
			else:
				order	 = (3,1,1)
		elif(type == "LTR"):
			if(subtype.startswith("ERV")):
				if(subtype == "ERV1"):
					order = (3,5,2)
				elif (subtype == "ERVK"):
					order = (3,5,3)
				elif (subtype == "ERVL_MaLR"):
					order = (3,5,4)
				elif (subtype == "ERVL"):
					order = (3,5,5)
			else:
				order	 = (3,5,1)
		elif(type == "Low_complexity"):
			order	 = (3,8,1)
		elif(type == "SINE"):
			if(subtype == "Alu"):
				order = (3,2,2)
			elif(subtype == "MIR"):
				order = (3,2,3)
			else:
				order = (3,2,1)
		elif(type == "Satellite"):
			order = (3,6,1)
		elif(type == "Simple_repeat"):
			order = (3,7,1)
		elif("RNA" in type):
			order = (3,4,1)
		else:
			order = (3,9,1)
	return order

def formatAnnotations(windowId, annotPerWin,CCDSonly=False):
	begPos,endPos = windowId
	annotPerWinFormat = {}
	# Annotation type (annotType) refers to one of the four 
	# UCSC Table Browser tracks used in this study:
	# (1) "Gencode-v4"; (2) "CCDS"; (3) "RepeatMasker"; (4) "RefSeqFuncElems".
	# A window can be mapped to none, one or many annotations with the same annotation type.
	for (annotType, annotLst) in annotPerWin.items():
		# Gather all annotations that overlap a window, regardless of their overlap.
		for annot in annotLst:
			begPosAnnot, endPosAnnot, strand, mainCat, subCat, elemId = annot
			# Get a formatted "nice" label.
			annotForm = getFormattedAnnot((annotType, mainCat, subCat),proteinClasses=None,CCDSonly=CCDSonly)
			if(annotForm[0] > 0):
				# Compute intersection between annotation and window in bps.
				begIntersec, endIntersec = max(begPos,begPosAnnot), min(endPos,endPosAnnot)
				if(annotForm not in annotPerWinFormat):
					annotPerWinFormat[annotForm] = []
				annotPerWinFormat[annotForm].extend(list(range(begIntersec,endIntersec)))
	# Compute % overlap that an annotation has with the window.
	# It makes sure that each position is covered only once per annotation.
	windowSize = endPos-begPos
	annotPerWinFormat = {annot: len(set(annotPerWinFormat[annot]))/windowSize for annot in annotPerWinFormat}
	return annotPerWinFormat

def getAnnotationsPerWin(parallelInput):
	windowLst, annotPerChr, stepDebug = parallelInput
	annotations = {}
	winCnt = 0
	for windowId in windowLst:
		begPos, endPos = windowId
		windowId = (begPos, endPos)
		# Find all annotations that intersect window.
		annot = intersectingAnnot(begPos, endPos, annotPerChr)
		# Format annotation categories.
		annotPerWinFormat     = formatAnnotations(windowId, annot, CCDSonly=False)
		annotations[windowId] = annotPerWinFormat
		# Print debug info.
		winCnt += 1
		if((stepDebug > 0) and ((winCnt % stepDebug)==0)): print(f"\t{winCnt+1} out of {len(windowLst)} ({int(100*winCnt/len(windowLst))} %)")
	return annotations

def chunks(lst, n):
	"""Yield successive n-sized chunks from lst."""
	for i in range(0, len(lst), n):
		yield lst[i:i + n]

def annotateWindows(my_dataset, allAnnotations, nbcores):

	# Save overall time.
	timeTrack = Time()
	timeTrack.start()

	timeTrack.startStep("Load window IDs")
	windowsIds = loadWindowIds(my_dataset)
	timeTrack.stopStep()

	# Compute chromosome sizes and whole-genome size.
	print(f"Size stats:")
	sizeGenome = 0
	for (qChrom, windowsPerChrom) in windowsIds.items():
		sizeChrom   = sum([endPos-begPos for (begPos,endPos) in windowsPerChrom])
		sizeGenome += sizeChrom
		print(f"{qChrom} : {sizeChrom} bps")
	print(f"Whole-genome : {sizeGenome} bps")

	# Annotate windows for all chromosomes.
	print(f"Annotating windows per chromosome")
	annotatedWindows = {}
	for qChrom in my_dataset.chromLst:
		annotatedWindows[qChrom] = {}
		# Get annotations per chromosome.
		annotPerChr = getAnnotationsPerChr(allAnnotations, qChrom)
		# Absolute coordinates.
		timeTrack.startStep(f"[{qChrom}] Annotate {len(windowsIds[qChrom])} windows")
		if(nbcores == 1):
			stepDebug = 10 #int(len(windowsIds[qChrom])*0.01)
			annotatedWindows[qChrom] = getAnnotationsPerWin((windowsIds[qChrom], annotPerChr, stepDebug))
		else:
			stepDebug      = int(len(windowsIds[qChrom])*0.05)
			stepCur        = 0
			windowIdLsts   = chunks(windowsIds[qChrom], 10000)
			parallelInputs = [(windowLst, annotPerChr, -1) for windowLst in windowIdLsts]
			with futures.ProcessPoolExecutor(nbcores) as pool:
				for result in pool.map(getAnnotationsPerWin, parallelInputs):
					annotatedWindowsPartial = result
					for (windowId, annotPerWinFormat) in annotatedWindowsPartial.items():
						annotatedWindows[qChrom][windowId] = annotPerWinFormat
						stepCur += 1
						if((stepDebug > 0) and ((stepCur % stepDebug)==0)): print(f"\t{stepCur+1} out of {len(windowsIds[qChrom])} ({int(100*stepCur/len(windowsIds[qChrom]))} %)")
		timeTrack.stopStep()

	# Gather all annotations found in the genome.
	timeTrack.startStep(f"Gather all annotations found in the genome")
	annotationsFound = set()
	for qChrom in my_dataset.chromLst:
		for (windowId, annotPerWinFormat) in annotatedWindows[qChrom].items():
			annotationsFound.update(list(annotPerWinFormat.keys()))
	annotationsFound = list(annotationsFound)
	timeTrack.stopStep()

	print(f"Save annotation file")
	annotFilename = my_dataset.getOutFilename_annotatedWindows()
	with open(annotFilename, 'wb') as pickleFile:
		pickle.dump((annotationsFound, annotatedWindows), pickleFile, protocol=pickle.HIGHEST_PROTOCOL)

	timeTrack.stop()
	timeTrack.print()

####################################
# MAIN.
####################################
# Use: python3 annotateWindows.py
if (__name__ == '__main__'):

	parser = argparse.ArgumentParser(description="Annotate windows.")
	parser.add_argument("-cores", help="Number of cores to be used during annotation.", type=int, required=True)

	args       = parser.parse_args()
	nbcores    = int(args.cores)

	my_dataset = Dataset()

	# Save overall time.
	timeTrack = Time()
	timeTrack.start()

	# Parse .tsv files w/annotations downloaded from UCSC Table Browser.
	timeTrack.startStep("Load all annotations")
	allAnnotations = loadAllAnnotations(my_dataset)
	timeTrack.stopStep() # It takes around 40 seconds.

	timeTrack.startStep("Annotate windows")
	annotateWindows(my_dataset, allAnnotations, nbcores)
	timeTrack.stopStep()

	timeTrack.stop()
	timeTrack.print()
