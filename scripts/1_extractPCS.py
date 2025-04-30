"""Reads a chain file (extension '.all.chain') 
downloaded from the UCSC website 
and outputs one pickle file per chromosome with 
all PCSs of a pair of species.

To make sure that PCSs are being properly extracted,
this script also needs the FASTA files of both species, 
which can be downloaded in the UCSC website.

- **Use**::

	python3 1_extractPCS.py -refsp_ucscname [UCSC name] -othsp_ucscname [UCSC name] -chain_file [/path/to/chain] -refsp_fastadir [/path/to/fasta] -othsp_fastadir [/path/to/fasta] -pcs_dir [/path/output] -pcs_minsize [number]

- **Example of Usage (human and mouse)**::

	python3 1_extractPCS.py -refsp_ucscname hg38 -othsp_ucscname mm39 -chain_file /path/to/ucsc/chains/mm39/mm39.hg38.all.chain -othsp_fastadir /path/to/fasta/mm39/chr -refsp_fastadir /path/to/fasta/hg38/chr -pcs_dir /path/to/outputdir -pcs_minsize 5

	python3 ~/code/1_extractPCS.py -refsp_ucscname hg38 -othsp_ucscname mm39 -chain_file /flash/MillerU/Priscila/paper-validation/chains/mm39.hg38.all.chain -othsp_fastadir /bucket/MillerU/Priscila/fasta/mm39/chr -refsp_fastadir /bucket/MillerU/Priscila/fasta/hg38/chr -pcs_dir /flash/MillerU/Priscila/paper-validation/pcs -pcs_minsize 5

- **Input Parameters (all mandatory)**:
:-refsp_ucscname: UCSC name of the reference species that is being aligned (e.g. *hg38* for human).
:-othsp_ucscname: UCSC name of the species that is being aligned (e.g. *mm39* for mouse).
:-chain_file: file with extension ``.all.chain`` downloaded from UCSC (e.g.: hg38.papAnu4.all.chain.gz)
:-refsp_fastadir: 	directory with FASTA files following the format ``genomeName.chrName.bXXXeXXX.fa`` (e.g. mm39.chrY.b51040000e51060000.fa).
					FASTA files can be downloaded in the UCSC website, and then split into smaller files for more efficient reading.
:-othsp_fastadir: 	directory with FASTA files following the format ``genomeName.chrName.bXXXeXXX.fa`` (e.g. mm39.chrY.b51040000e51060000.fa).
					FASTA files can be downloaded in the UCSC website, and then split into smaller files for more efficient reading.
:-pcs_dir: Directory where files with PCSs are going to be saved.
:-pcs_minsize: minimum size (nb. of base pairs) for the PCS to be considered.

- **Output**: 
	A separate ``.pickle`` file for each chromosome, each including all PCSs 
	positioned on the corresponding chromosome of the reference genome. 
	The ``pickle`` file saves a list of ``Pcs`` named tuples, where each tuple 
	has the size of a PCS and its coordinates in both genomes.

.. note::
	**Required Directory Structure for FASTA files**: a single FASTA file 
	from UCSC (e.g. ``mm39.fa``) must be split in many files, 
	following the pattern ``genomeName.chrName.bXXXeXXX.fa`` (see Pre-Requisites below).

Pre-requisites
--------------

	Before using this script, make sure all the required files were downloaded:

a) Chain files from UCSC
^^^^^^^^^^^^^^^^^^^^^^^^

1. Download the "chain" file for the analyzed pairwise species 
(chained lastz alignments) in the UCSC website. 
For example, for human and mouse (*hg38* and *mm39*), go to::

	https://hgdownload.soe.ucsc.edu/goldenPath/mm39/vsHg38/

and download the file ``mm39.hg38.all.chain.gz``.

2. Unzip the download file::

	gzip -d mm39.hg38.all.chain.gz

3. The full directory of the unziped chain file (e.g. ``mm39.hg38.all.chain``) 
will be given to this script through the parameter ``-chain_file``.

b) FASTA files from UCSC
^^^^^^^^^^^^^^^^^^^^^^^^
4. Download the FASTA files for the analyzed pairwise species in the UCSC website (e.g. ``mm39.fa.gz``).
The FASTA files are needed to validate the PCS coordinates obtained from the chain file.

5. Run the script ``breakFasta.py`` using the unziped FASTA file (e.g. ``mm39.fa``). 
This script will *break* the unziped FASTA file into smaller files, following the pattern::

	genomeName.chrName.bXXXeXXX.fa

6. These files must be saved in a directory structure organized as follows::

	genomeName/chr/chrName/genomeName.chrName.bXXXeXXX.fa

Cluster resources
-----------------
In case you want to run this Python script stand-alone::

	srun -p compute -t 2:00:00 --mem 20G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash

Otherwise you can use the script ``../cluster/1-extractPCS-runAll.py``
to run this script for all 40 species used in this study.

Time
----

Stats on time (human-mouse alignment, whole-genome): **~52 minutes**

==================  ===================
Step                Time (s)
==================  ===================
Reading chains      136.36537051200867
Compatible chains   182.21501278877258
PCSs from chr1      120.18683695793152
PCSs from chr2      115.76346206665039
PCSs from chr3      61.49129343032837
PCSs from chr4      182.26222944259644
PCSs from chr5      194.65999698638916
PCSs from chr6      181.5871503353119
PCSs from chr7      177.63377285003662
PCSs from chr8      175.2457103729248
PCSs from chr9      135.67597675323486
PCSs from chr10     154.27861976623535
PCSs from chr11     143.3001594543457
PCSs from chr12     143.43108081817627
PCSs from chr13     102.48199248313904
PCSs from chr14     79.90870404243469
PCSs from chr15     91.05354380607605
PCSs from chr16     76.45854377746582
PCSs from chr17     72.0199842453003
PCSs from chr18     70.00186824798584
PCSs from chr19     97.53436303138733
PCSs from chr20     59.93760871887207
PCSs from chr21     32.661128997802734
PCSs from chr22     39.76104974746704
PCSs from chrX      195.10439586639404
PCSs from chrY      19.989043951034546
**Total time**      3103.8352103233337
==================  ===================

Storage
-------

Size of output files with PCSs for the *human* and *mouse* 
alignment (whole-genome): **~1085 MB**.

========================== ======
Output file                Size
========================== ======
hg38.mm39.chr1.pcs.pickle  98M 
hg38.mm39.chr2.pcs.pickle  104M 
hg38.mm39.chr3.pcs.pickle  72M 
hg38.mm39.chr4.pcs.pickle  59M 
hg38.mm39.chr5.pcs.pickle  74M 
hg38.mm39.chr6.pcs.pickle  66M 
hg38.mm39.chr7.pcs.pickle  56M 
hg38.mm39.chr8.pcs.pickle  55M 
hg38.mm39.chr9.pcs.pickle  44M 
hg38.mm39.chr10.pcs.pickle 46M 
hg38.mm39.chr11.pcs.pickle 65M 
hg38.mm39.chr12.pcs.pickle 50M 
hg38.mm39.chr13.pcs.pickle 33M 
hg38.mm39.chr14.pcs.pickle 39M 
hg38.mm39.chr15.pcs.pickle 34M 
hg38.mm39.chr16.pcs.pickle 30M 
hg38.mm39.chr17.pcs.pickle 24M 
hg38.mm39.chr18.pcs.pickle 30M 
hg38.mm39.chr19.pcs.pickle 14M 
hg38.mm39.chr20.pcs.pickle 28M 
hg38.mm39.chr21.pcs.pickle 12M 
hg38.mm39.chr22.pcs.pickle 11M 
hg38.mm39.chrX.pcs.pickle  41M 
hg38.mm39.chrY.pcs.pickle  387K 
========================== ======

Function details
----------------

Only relevant functions have been documented below. 
For more details on any function, check the comments in the souce code.

"""
import os
import re
import sys
from collections import defaultdict
import pickle
import glob
import time
import argparse

from utils.basicTypes import Chain, Block, Pcs, Time

########################################################
# Class used to validate if PCSs are being correctly extracted.

class PcsSingleCheck:
	"""This class makes sure that PCSs are being correctly processed. 
	For every new PCS, it determines whether the PCS is a good candidate 
	to be checked into details. Good candidates are long PCSs,
	PCSs in the reverse strand of one of the genomes, among others.
	It keeps a history of how many PCSs have been checked so far.

	:param desc: Brief description of which type of check is being performed. 
	:type desc: str
	:param cntUB: Maximum number of PCSs that will be checked in details.
	:type cntUB: int
	:param pcsSingleCheckFunc: Function that receives a named tuple PCS 
	and returns True (PCS should be checked) or False.
	:type pcsSingleCheckFunc: function
	"""
	def __init__(self, desc, cntUB, pcsSingleCheckFunc):
		"""Docstring for class Foo."""
		self.desc  = desc
		self.cntUB = cntUB
		self.pcsSingleCheckFunc = pcsSingleCheckFunc

		self.reset()

	def reset(self):
		self.cntCur	   = 0
		self.lastPCSstate = False

	def print(self):
		return f"{self.desc} = {self.cntCur}"

	def checkPCS(self, pcs):
		if((self.cntCur < self.cntUB) and (self.pcsSingleCheckFunc(pcs))):
			self.cntCur += 1
			self.lastPCSstate = True
		else:
			self.lastPCSstate = False
		return self.lastPCSstate

def checkPCSsize(pcs):
	pcs_size_lb = 80 # At least 70 base pairs to be considered "big".
	return (pcs.size >= pcs_size_lb)

def checkPCSrevT(pcs):
	return (pcs.tStrand == "-")

def checkPCSrevQ(pcs):
	return (pcs.qStrand == "-")

def initializePCSchecks(cntUB=100,prefixTarget="tar",prefixQuery="qry"):
	return [PcsSingleCheck(f"PCSs with large size", cntUB,checkPCSsize),
			PcsSingleCheck(f"PCSs in reverse strand ({prefixTarget})", cntUB,checkPCSrevT),
			PcsSingleCheck(f"PCSs in reverse strand ({prefixQuery})",  cntUB,checkPCSrevQ)]

def performPCSchecks(pcs, checks):
	res = False
	for c in checks:
		res = c.checkPCS(pcs)
		if (res) : break
	return res

########################################################
# DNA processing methods.

def reverseComplement(seq):
	"""This function returns a reverse complement of a DNA."""
	# Complement strand		
	seq = seq.translate(str.maketrans('ATCGatcg','TAGCtagc'))
	# Reverse strand
	seq = seq[::-1]
	return seq

def readDNA(dirFASTA, prefix, chrName, strand, posBeg, posEnd, debug=False):
	"""This function reads a DNA segment in a FASTA file."""
	isRevComp = (strand == "-")
	dna = ""

	# FASTA files that contain fragment.
	# genomeName.chrName.bXXXeXXX.fa
	selectedFiles = []
	dnafileNamePattern = os.path.join(dirFASTA, f"{prefix}.{chrName}.b*.fa")
	dnafileNamesLst	= glob.glob(dnafileNamePattern)
	dnafileNamesLstUpd = []
	for fileName in dnafileNamesLst:
		m = re.match(r".*b(\d+)e(\d+).fa$", fileName)
		if m:
			begDNAFrag = int(m.group(1))
			endDNAFrag = int(m.group(2))
			dnafileNamesLstUpd.append((begDNAFrag, fileName))
	dnafileNamesLstUpd.sort(key=lambda x: x[0])
	dnafileNamesLst = [fileName for (pos,fileName) in dnafileNamesLstUpd]
	if (len(dnafileNamesLst) == 0):
		print(f"WARNING! Chromosome not found: {chrName} in {dnafileNamePattern}")
		return dna

	for fileName in dnafileNamesLst:
		m = re.match(r".*b(\d+)e(\d+).fa$", fileName)
		if m:
			begDNAFrag = int(m.group(1))
			endDNAFrag = int(m.group(2))
			if (endDNAFrag >= posBeg):
				selectedFiles.append((begDNAFrag, fileName))
			if (posEnd <= endDNAFrag):
				break
	selectedFiles.sort(key=lambda x: x[0])
	selectedFiles = [fileName for (pos,fileName) in selectedFiles]
	#print(*[f.split("/")[-1] for f in selectedFiles],sep=", ")

	# Open FASTA files and read.
	#print("Reading FASTA... posBeg=" + str(posBeg) + " posEnd=" + str(posEnd))
	for fastafileName in selectedFiles:
		m = re.match(r".*b(\d+)e(\d+).fa$", fastafileName)
		if m:
			begDNAFrag = int(m.group(1))
			endDNAFrag = int(m.group(2))
			
			with open(fastafileName) as f:
				dnaFragment = f.read()

			# First file.
			if(begDNAFrag <= posBeg):
				if(posEnd <= endDNAFrag):								
					dna = dnaFragment[(posBeg-begDNAFrag):(posEnd-begDNAFrag)]
				else:
					dna = dnaFragment[(posBeg-begDNAFrag):]
			# Last file.
			elif(posEnd <= endDNAFrag):
				dna += dnaFragment[:(posEnd-begDNAFrag)]
			# Intermediary files.
			else:
				dna += dnaFragment
	if isRevComp: dna = reverseComplement(dna)
	if (len(dna) != (posEnd-posBeg)): print("ERROR! Different sizes: " + str(posEnd-posBeg) + " " + str(len(dna)))
	return dna

########################################################
# + : [tStart, tEnd)
# - : (tSize-tEnd, tSize-tStart]
def computePos(strand,size,start,end):
	if(strand == "+"):
		return start, end
	elif (strand == "-"):
		return (size-end), (size-start)
	else:
		print("WARNING! Unidentified strand=" + strand)
		return start,end

#@profile
def readChains(chainsFilename):
	"""This function parses the chains from Chain file."""

	patternChain	 = re.compile(r"\s*chain\s*([\d]+)\s*([\w]+)\s*([\d]+)\s*(\S+)\s*([\d]+)\s*([\d]+)\s*([\w]+)\s*([\d]+)\s*(\S+)\s*([\d]+)\s*([\d]+)\s*([\d]+)\s*$")
	patternBlock	 = re.compile(r"\s*([\d]+)\s+([\d]+)\s+([\d]+)\s*$")
	patternLastBlock = re.compile(r"\s*([\d]+)\s*$")

	chains=[]
	blocks={}

	# File structure:
	#chain 3000 chr6 172126628 + 168948273 168948415 chr19 61420004 - 47155418 47155564 3258284
	#26	0	4
	#18	4	0
	#26	0	4
	#18	4	0
	#26	0	4
	#20
	#size dt dq
	#	size -- the size of the ungapped alignment
	#	dt -- the difference between the end of this block and the beginning of the next block (reference/target sequence)
	#	dq -- the difference between the end of this block and the beginning of the next block (query sequence)
	# NOTE: The last line of the alignment section contains only one number: the ungapped alignment size of the last block.
	#
	with open(chainsFilename) as f:
		# chain score		 tName tSize		 tStrand tStart	 tEnd			qName qSize		 qStrand qStart	 qEnd			id 
		# chain 619040156 chr14 101161492 +			 18886862 100137521 chr12 120092757 +			 44683557 113394468 1
		readChainDetails = False
		chain = None
		for line in f:
			m = patternChain.match(line)
			if m:
				score = int(m.group(1))
				chainId = int(m.group(12))

				# Target sequence.
				tName = m.group(2)
				tSize = int(m.group(3))
				tStrand = m.group(4)
				tStart = int(m.group(5))
				tEnd	 = int(m.group(6))
				tStart, tEnd = computePos(tStrand,tSize,tStart,tEnd)

				# Query sequence.
				qName = m.group(7)
				qSize = int(m.group(8))
				qStrand = m.group(9)
				qStart = int(m.group(10))
				qEnd	 = int(m.group(11))
				qStart, qEnd = computePos(qStrand,qSize,qStart,qEnd)						
				chain = Chain(score, chainId, tName, tSize, tStrand, tStart, tEnd, qName, qSize, qStrand, qStart, qEnd)
				chains.append(chain)

				readChainDetails = True
				blocks[chain]	= []

			elif readChainDetails:
				m = patternBlock.match(line)
				if m:
					igs	= int(m.group(1))
					tgap = int(m.group(2))
					qgap = int(m.group(3))

					block = Block(igs, tgap, qgap)
					blocks[chain].append(block)
				else:
					m = patternLastBlock.match(line)
					if m:
						igs	= int(m.group(1))
						
						block = Block(igs, -1, -1)
						blocks[chain].append(block)

						readChainDetails = False
						chain = None
	return chains, blocks

def computePCScoords(begAbs,posRel,isRevComp,sizePCS,sizeMis):
	pcsAbsPosBeg = begAbs-posRel+sizeMis if (isRevComp) else begAbs+posRel-(sizeMis+sizePCS)
	return pcsAbsPosBeg

def dir_path(path):
	if os.path.isdir(path):
		return path
	else:
		raise argparse.ArgumentTypeError(f"ERROR! {path} is not a valid path. Check if the directory exists.")

def file_path(path):
	if os.path.isfile(path):
		return path
	else:
		raise argparse.ArgumentTypeError(f"ERROR! {path} is not a valid path. Check if the file exists.")

####################################
# MAIN.
####################################
if (__name__ == '__main__'):

	parser = argparse.ArgumentParser(description="Extract PCSs (Perfectly Conserved Sequences) from Chains (ordered aligned blocks).")

	parser.add_argument("-refsp_ucscname", help="UCSC name of the species used as reference point. In our study, it is always 'hg38' (human genome).", type=str, required=True)
	parser.add_argument("-othsp_ucscname", help="UCSC name of the other species.", type=str, required=True)

	parser.add_argument("-chain_file", help="File with the extension '*.all.chain' downloaded from UCSC website.", type=file_path, required=True)

	parser.add_argument("-refsp_fastadir", help="Directory with FASTA files in the format: genomeName.chrName.bXXXeXXX.fa", type=dir_path, required=True)
	parser.add_argument("-othsp_fastadir", help="Directory with FASTA files in the format: genomeName.chrName.bXXXeXXX.fa", type=dir_path, required=True)

	parser.add_argument("-pcs_dir", help="Directory where PCSs will be saved (output directory).", type=dir_path, required=True)
	parser.add_argument("-pcs_minsize", help="Minimum size (in number of base pairs) for the PCS to be considered.", type=int, required=True)

	args = parser.parse_args()

	prefixTarget   = args.othsp_ucscname
	prefixQuery	   = args.refsp_ucscname # In our study, query is always the human genome (hg38)
	chainsFilename = args.chain_file	 # File *.all.chain downloaded from UCSC website. 
	tdirFASTA	  = args.othsp_fastadir # Directory with FASTA files in the format: genomeName.chrName.bXXXeXXX.fa
	qdirFASTA	  = args.refsp_fastadir # Directory with FASTA files in the format: genomeName.chrName.bXXXeXXX.fa
	outDir		 = args.pcs_dir
	minPCSsize	 = int(args.pcs_minsize)
		
	qChromLst=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]

	print("******************************************************")
	print("*		PCS (Perfectly conserved sequences)		 *")
	print("* Extract PCSs from chains (ordered aligned blocks). *")
	print("******************************************************")
	print(f"- Reference genome: {prefixTarget}")
	print(f"- Query genome: {prefixQuery}")
	print(f"- Chain file: {chainsFilename}")
	print(f"- Directory for FASTA file of {prefixTarget}: {tdirFASTA}")
	print(f"- Directory for FASTA file of {prefixQuery}: {qdirFASTA}")
	print(f"- PCS directory (output): {outDir}")
	print(f"- Minimum size for PCS to be considered: >={minPCSsize} base pairs.")

	##############################
	# Check parameters given by the user.
	if(not os.path.isfile(chainsFilename)):
		sys.exit(f"ERROR! Invalid filepath for chain file (*.all.chain): {chainsFilename}")

	# Save times for each step and overall time.
	start_all = time.time()
	times_lst = []

	##############################
	# Read header and details of all chains.
	# This step takes between 3 and 4 minutes.
	print("Reading chains...")
	start = time.time()
	chains, blocks = readChains(chainsFilename)
	end   = time.time()
	times_lst.append(Time("Reading chains",end-start))
	print(f"Reading chains : t={end-start}") # It takes around 3--4 minutes.
	
	# Sort chains by score.
	chains.sort(reverse=True, key=lambda x: x.score)
	
	##############################
	# Check if chains make sense.

	# Check if all intervals are increasing in position.
	for idx, chain in enumerate(chains):
		if ((chain.tStart >= chain.tEnd)
		 or (chain.qStart >= chain.qEnd)):
			print(f"WARNING! Problematic chain: target={chain.tStart},{chain.tEnd}; query={chain.qStart},{chain.qEnd}")

	##############################
	# Select chains that are compatible.
	# This step takes around 4 minutes.

	print("Selecting compatible chains...")
	start = time.time()
	# Check if a sequence is compatible (non-intersecting) to a list of sequences.
	def isCompatible(pos, pos_lst):
		b1, e1 = pos
		for (b2, e2) in pos_lst:
			if (((b1 >= b2) and (b1 <= e2))
			or	((e1 >= b2) and (e1 <= e2))
			or	((b2 >= b1) and (b2 <= e1))
			or	((e2 >= b1) and (e2 <= e1))):
				return False
		return True

	# Check if a chain is compatible to **both** genomes.
	tChains = defaultdict(list)
	qChains = defaultdict(list)
	compatibleChainsCount = 0
	for idx, chain in enumerate(chains):
		if (idx % 100000 == 0):
			print("Step " + str(idx) + " out of " + str(len(chains)) + "(" + str(int((idx/len(chains))*100)) + "%)")

		# Select chains grouped by target chromosomes.
		tChrom	 = chain.tChrom
		tSelChains = [(c.tStart, c.tEnd) for c in tChains[tChrom]] if tChrom in tChains else []
		tComp	   = isCompatible((chain.tStart, chain.tEnd), tSelChains)

		# Select chains grouped by query chromosomes.
		qChrom	 = chain.qChrom
		qSelChains = [(c.qStart, c.qEnd) for c in qChains[qChrom]] if qChrom in qChains else []
		qComp	  = isCompatible((chain.qStart, chain.qEnd), qSelChains)

		# Check if current chain is compatible with already selected chains.
		if tComp and qComp :
			tChains[tChrom].append(chain)
			qChains[qChrom].append(chain)
			compatibleChainsCount += 1

	end = time.time()
	times_lst.append(Time("Compatible chains",end-start))
	print(f"Selecting compatible chains : t={end-start}") # It takes around 4 minutes.
	print(f"Compatible chains: {compatibleChainsCount} out of {len(chains)} ({100*(compatibleChainsCount/len(chains)):.2f})")

	# Sort chains per starting position in each chromosome.
	for qChrom in qChains.keys():
		qChains[qChrom] = sorted(qChains[qChrom], key=lambda x: x.qStart)

	##############################
	# Compute PCSs from selected chains.
	# This step takes between 2 and 3 minutes per chromosome.
	print("Computing PCSs from selected chains...")
	for qChrom in qChromLst:
		print(f"\tChains in chromosome {qChrom} ({len(qChains[qChrom])} chains)")
		if(qChrom not in qChains):
			print(f"WARNING! Chains for chromosome {qChrom} NOT found.")
			continue

		outFilename = os.path.join(outDir, f"{prefixQuery}.{prefixTarget}.{qChrom}.pcs.pickle")

		##############################
		# Validation setup.
		# This step is important to make sure that PCSs are correctly extracted.
		pcsChecks = initializePCSchecks(cntUB=150,prefixTarget=prefixTarget,prefixQuery=prefixQuery)

		##############################
		# Read DNA from chains and extract PCSs.
		start   = time.time()
		pcs_lst = []
		for chain in qChains[qChrom]:

			# At the end of an iteration, the chain should have all PCSs computed.
			tDNA = readDNA(os.path.join(tdirFASTA,chain.tChrom), prefixTarget, chain.tChrom, chain.tStrand, chain.tStart, chain.tEnd)
			qDNA = readDNA(os.path.join(qdirFASTA,chain.qChrom), prefixQuery,  chain.qChrom, chain.qStrand, chain.qStart, chain.qEnd)

			# Sanity check before continue.
			validSeqs = ((chain.tEnd-chain.tStart) == len(tDNA)) and ((chain.qEnd-chain.qStart) == len(qDNA))

			# Check blocks.
			if (not validSeqs):
				print("WARNING! Problem to read DNA segments from chain. Chain will be ignored.")
				problematicChain = chain
			else:
				# try:
				qIsRevComp = (chain.qStrand == "-")
				qOp		   = -1 if qIsRevComp else 1
				qBegAbs	   = chain.qEnd if qIsRevComp else chain.qStart
				qPosRel	= 0

				tIsRevComp = (chain.tStrand == "-")
				tOp		   = -1 if tIsRevComp else 1
				tBegAbs	= chain.tEnd if tIsRevComp else chain.tStart
				tPosRel	= 0
				
				for block in blocks[chain]:
					# Size of IGS (Intergap sequence), which may contain matches and mismatches, but no gaps.
					igs	= block.igs 

					# Read DNA sequence from block.
					tBlockSeq = tDNA[tPosRel:(tPosRel+igs)].upper()
					qBlockSeq = qDNA[qPosRel:(qPosRel+igs)].upper()

					# Compute PCS (Perfectly Conserved Sequences).
					sizePCS = 0
					sizeMis = 0

					def addPCS():
						if (sizePCS >= minPCSsize):

							tPosBeg = computePCScoords(tBegAbs,tPosRel,tIsRevComp,sizePCS,sizeMis)
							qPosBeg = computePCScoords(qBegAbs,qPosRel,qIsRevComp,sizePCS,sizeMis)

							pcs = Pcs(sizePCS, chain.tChrom, chain.tStrand, tPosBeg, chain.qChrom, chain.qStrand, qPosBeg)

							#####################################
							# Validation.
							# Check if the DNA sequences associated with a PCS should be validated.
							if(performPCSchecks(pcs, pcsChecks)):

								# Read DNA sequence from PCS.
								tPCS = readDNA(os.path.join(tdirFASTA,chain.tChrom), prefixTarget, chain.tChrom, chain.tStrand, tPosBeg, tPosBeg+sizePCS)
								qPCS = readDNA(os.path.join(qdirFASTA,chain.qChrom), prefixQuery,  chain.qChrom, chain.qStrand, qPosBeg, qPosBeg+sizePCS)

								if(tPCS.upper() != qPCS.upper()):
									print(f"PCS validation: ERROR!")
									print(f"PCSs do not match: {tPCS} (coords={prefixTarget},{chain.tChrom},{chain.tStrand},{tPosBeg}) and {qPCS} (coords={prefixQuery},{chain.qChrom},{chain.qStrand},{qPosBeg})")
									print(f"Chain start : Qry={qBegAbs}\tTar={tBegAbs}")
									print(f"Block size: {igs}\n{prefixTarget}:\n{tDNA}\n{prefixQuery}:\n{qDNA}")
									print(f"Extra info to help debugging: {tIsRevComp=} {tBegAbs=} {tPosRel=} {sizeMis=} {sizePCS=}")
									print(f"Extra info to help debugging: {qIsRevComp=} {qBegAbs=} {qPosRel=} {sizeMis=} {sizePCS=}")
									sys.exit()
								# else:
								# 	print(f"PCS validation: OK! [Stats on PCSs checked so far: {';'.join([c.print() for c in pcsChecks])}]")

							#####################################
							# Save PCS.
							pcs_lst.append(pcs)

					# Find PCSs in a block.
					for idx in range(len(tBlockSeq)):
						if (tBlockSeq[idx] == qBlockSeq[idx]):
							if (sizeMis > 0):
								addPCS()
								sizePCS = 0
								sizeMis = 0
							sizePCS += 1
						else:
							sizeMis += 1
						tPosRel += 1
						qPosRel += 1
					# Save last PCS of the block.
					addPCS()

					# Update next positions.
					tPosRel += block.tGap
					qPosRel += block.qGap

				# except Exception as e:
				# 	print(f"Problematic Chain:\n{e}\n{tDNA}\n{qDNA}\n")
		end = time.time()
		times_lst.append(Time(f"PCSs from {qChrom}",end-start))
		print(f"\tComputing PCSs from {qChrom} (total found={len(pcs_lst)} PCSs): t={end-start}")
		print(f"\tPCS validation: OK! [Stats: {';'.join([c.print() for c in pcsChecks])}]")

		##############################
		# Save PCSs.
		print("Save PCSs per chromosome...")
		with open(outFilename, 'wb') as pickleFile:
			pickle.dump(pcs_lst, pickleFile, protocol=pickle.HIGHEST_PROTOCOL)
		print(f"Pickle saved! ({outFilename})")

	end_all = time.time()
	times_lst.append(Time(f"Total time:",end_all-start_all))
	print("Stats on time:\n")
	for t in times_lst: print(f"\t{t.desc} {t.t}")
	print("PCS: Done!")

