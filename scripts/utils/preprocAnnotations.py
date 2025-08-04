"""Parser for TSV files downloaded from *UCSC Table Browser*.
It saves a pickle file containing annotations grouped by chromosome.
Four tracks are considered: 
  1. **CCDS**: Transcript coding sequence (CDS) genomic regions;
  2. **GENCODE V44**;
  3. **RepeatMasker**;
  4. **RefSeq**: Functional elements from RefSeq.

**Use**::
	
	python3 preprocAnnotations.py

Time, Memory & Disk space
-------------------------

Running the script on a single core takes **21.14 seconds** and requires a small amount of memory. 

**Output files**:
	
	The script creates the file ``annotations.hg38.pickle``, which takes **312.3 MB** of space.

"""
import os
import pickle
import collections
import sys

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

####################################
# MAIN.
####################################
# Use: python3 preprocAnnotations.py
if (__name__ == '__main__'):

	my_dataset = Dataset()

	# Save overall time.
	timeTrack = Time()
	timeTrack.start()

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

	annotFilename  = my_dataset.getOutFilename_annotations()
	allAnnotations = [annotationsGencode, annotationsCCDS, annotationsRepeatMasker, annotationsRefSeqFuncElems]
	with open(annotFilename, 'wb') as pickleFile:
		pickle.dump(allAnnotations, pickleFile, protocol=pickle.DEFAULT_PROTOCOL) # HIGHEST_PROTOCOL

	timeTrack.stop()
	timeTrack.print()
