"""Reads a chain file (extension '.all.chain') 
downloaded from the UCSC website 
and outputs one pickle file per chromosome with 
all PCSs of a pair of species.

To make sure that PCSs are being properly extracted,
this script also needs the FASTA files of both species, 
which can be downloaded in the UCSC website.

- **Use**::

	python3 1_extractPCS.py -sp_ucsc_name [UCSC name]

- **Example of Usage (human (reference genome) and mouse)**::

	python3 ~/code/1_extractPCS.py -sp_ucsc_name mm39

- **Input Parameter (mandatory)**:
:-sp_ucsc_name: UCSC name of the species that is being aligned with the reference species (e.g. *mm39* for mouse).

- **Other Parameters taken from** ``dataset.py``:
:refsp_ucscName: UCSC name of the reference species that is being aligned (e.g. *hg38* for human).
:chain_file: file with extension ``.all.chain`` downloaded from UCSC (e.g.: hg38.papAnu4.all.chain.gz)
:refsp_fastadir: 	directory with FASTA files following the format ``genomeName.chrName.bXXXeXXX.fa`` (e.g. mm39.chrY.b51040000e51060000.fa).
					FASTA files can be downloaded in the UCSC website, and then split into smaller files for more efficient reading.
:othsp_fastadir: 	directory with FASTA files following the format ``genomeName.chrName.bXXXeXXX.fa`` (e.g. mm39.chrY.b51040000e51060000.fa).
					FASTA files can be downloaded in the UCSC website, and then split into smaller files for more efficient reading.
:pcs_dir: Directory where files with PCSs are going to be saved.
:pcs_minsize: minimum size (nb. of base pairs) for the PCS to be considered.

.. note::
	Make sure that the required parameters described above are correctly defined in the file ``utils/dataset.py``.

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
In case you want to run this Python script stand-alone in a cluster that uses Slurm to manage jobs::

	srun -p compute -t 2:00:00 --mem 20G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash

Otherwise you can use the script ``../cluster/1_extractPCS_runAll.py``
to run this script for all 40 species used in this study.

Time, Memory & Disk space
-------------------------

For reference, here we include an **upper limit** on runtime, memory usage, and disk 
space required for running this script on each pairwise alignment of the 40 
vertebrate dataset examined in our study.

===========  =========  =========  ======
UCSC name    Time       Memory     Disk  
===========  =========  =========  ======
panPan3      01:01:19   5GB        0.95GB
panTro6      00:59:09   6GB        0.98GB
gorGor6      00:50:02   5GB        1.14GB
ponAbe3      01:20:13   5GB        1.98GB
papAnu4      06:27:13   38GB       2.76GB
macFas5      02:42:13   20GB       2.71GB
rhiRox1      03:29:32   23GB       2.98GB
chlSab2      00:59:52   6GB        2.99GB
nasLar1      04:26:24   15GB       2.25GB
rheMac10     01:01:55   6GB        2.80GB
calJac4      02:50:19   36GB       3.14GB
tarSyr2      04:57:09   42GB       2.62GB
micMur2      02:24:55   17GB       2.25GB
galVar1      09:48:31   40GB       2.56GB
mm39         00:58:45   9GB        1.10GB
oryCun2      01:23:21   19GB       1.67GB
rn7          01:27:01   24GB       1.07GB
vicPac2      01:23:54   17GB       2.19GB
bisBis1      02:50:09   24GB       1.86GB
felCat9      01:57:19   28GB       2.19GB
manPen1      03:18:57   18GB       1.99GB
bosTau9      01:10:49   15GB       1.59GB
canFam6      01:00:37   18GB       2.07GB
musFur1      01:00:55   17GB       2.21GB
neoSch1      01:20:42   19GB       2.46GB
equCab3      01:30:38   17GB       2.60GB
myoLuc2      01:46:55   20GB       1.60GB
susScr11     01:20:55   18GB       1.93GB
enhLutNer1   01:27:55   18GB       2.26GB
triMan1      01:20:32   20GB       2.16GB
macEug2      05:58:51   18GB       0.22GB
ornAna2      01:20:04   13GB       0.17GB
aptMan1      00:45:15   3GB        0.12GB
galGal6      00:45:15   3GB        0.09GB
thaSir1      01:30:08   5GB        0.07GB
aquChr2      00:30:03   3GB        0.12GB
melGal5      01:15:00   3GB        0.09GB
xenLae2      01:15:46   6GB        0.06GB
xenTro10     02:40:30   40GB       0.06GB
danRer11     04:08:26   14GB       0.04GB
===========  =========  =========  ======

Time per Run: Details
^^^^^^^^^^^^^^^^^^^^^

Stats on time of a single run (human-mouse alignment, whole-genome): **~65 minutes**
Details on computational time are available in the log of the run.

==================  =============
Step                Time (s)     
==================  =============
Reading chains      164.45332742 
Compatible chains   236.70236850 
PCSs from chr1      246.49562716 
PCSs from chr2      283.03650451 
PCSs from chr3      218.50386620 
PCSs from chr4      198.70389938 
PCSs from chr5      206.89367223 
PCSs from chr6      199.71864271 
PCSs from chr7      204.75429463 
PCSs from chr8      199.65774012 
PCSs from chr9      150.57288957 
PCSs from chr10     173.72448373 
PCSs from chr11     156.16801929 
PCSs from chr12     168.79578829 
PCSs from chr13     113.74493265 
PCSs from chr14     93.78934121  
PCSs from chr15     99.19000816  
PCSs from chr16     91.97542357  
PCSs from chr17     85.33002925  
PCSs from chr18     83.99309540  
PCSs from chr19     102.53493857 
PCSs from chr20     72.90905333  
PCSs from chr21     37.42282581  
PCSs from chr22     44.24664927  
PCSs from chrX      223.61548066 
PCSs from chrY      21.73898244  
**Total time**      3949.23640776
==================  =============

Storage per Run: Details
^^^^^^^^^^^^^^^^^^^^^^^^

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

Genome Browser
--------------

Here there are instructions in case you wish to check the corresponding 
DNA sequences of a PCS in UCSC's Genome Browser. Take for example the 
following PCS found in the pairwise alignment between human (hg38) and bonobo (panPan3)::

	Pcs(size=45, tChrom='chr10', tStrand='+', tPosBeg=87682, qChrom='chr16', qStrand='-', qPosBeg=21852)

Here, ``qChrom``, ``qStrand``, and ``qPosBeg`` correspond to the chromosome, strand, 
and starting coordinate in the query genome. In our dataset, the query genome is 
always the human genome (hg38). Similarly, ``tChrom``, ``tStrand``, and ``tPosBeg`` 
refer to the same information in the target genome (panPan3 in this example).

To compute the coordinates, first add ``+1`` to ``tPosBeg`` and ``qPosBeg``, as the 
coordinates in *Genome Browser* start at 1, whereas here they start at 0.

Next, by summing the size of the PCS (size=45) to 
the initial coordinates (``tPosBeg+1``  and ``qPosBeg+1``), it is possible 
to retrieve the DNA sequence in Genome Browser for both genomes.
 
For this example, by accessing::
	
	https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38

and typing ``chr16:21853-21897``, we obtain::

	TAATGAGTGTACTGCTCCTAAACAAGGAGTATCTGCATTTATTTA

If it is indicated that sequence should be read in the negative strand (qStrand='-'),
the reverse complement should be computed::

	TAAATAAATGCAGATACTCCTTGTTTAGGAGCAGTACACTCATTA

Similarly for the bonobo genome, by accessing::

	https://genome.ucsc.edu/cgi-bin/hgTracks?db=panPan3

and typing ``chr10:87683-87727``, we obtain::

	TAAATAAATGCAGATACTCCTTGTTTAGGAGCAGTACACTCATTA

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
import argparse

from utils.basicTypes import Chain, Block, Pcs, Time
from utils.dataset import Dataset

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
	"""This function parses the chains from a chain file (check the UCSC 
	website for more information on the format).
	"""

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

def isCompatible(pos, pos_lst):
	"""
	Check if a sequence is compatible (non-intersecting) to a list of sequences.
	"""
	b1, e1 = pos
	for (b2, e2) in pos_lst:
		if (((b1 >= b2) and (b1 <= e2))
		or	((e1 >= b2) and (e1 <= e2))
		or	((b2 >= b1) and (b2 <= e1))
		or	((e2 >= b1) and (e2 <= e1))):
			return False
	return True

def checkSpecies(dataset, ucscName):
	chainFilename = dataset.getChainFilename(ucscName,compressed=False)
	fastaDir = dataset.getDirFasta(ucscName)
	if ((chainFilename == None) or (fastaDir == None)): sys.exit()

####################################
# MAIN.
####################################
if (__name__ == '__main__'):

	parser = argparse.ArgumentParser(description="Extract PCSs (Perfectly Conserved Sequences) from Chains (ordered aligned blocks).")
	parser.add_argument("-sp_ucsc_name", help="UCSC name of the species that is being compared to the reference species. In our study, the reference species is 'hg38' (human genome), whereas the other species is one of the 40 vertebrates.", type=str, required=True)

	args       = parser.parse_args()
	my_dataset = Dataset()

	prefixTarget   = args.sp_ucsc_name
	checkSpecies(my_dataset, prefixTarget)

	prefixQuery	   = my_dataset.refsp_ucscName              # In our study, query is always the human genome (hg38)
	qChromLst      = my_dataset.chromLst

	chainsFilename = my_dataset.getChainFilename(prefixTarget,compressed=False) # File *.all.chain downloaded from UCSC website. 
	tdirFASTA	   = my_dataset.getDirFasta(prefixTarget)   # Directory with FASTA files in the format: genomeName.chrName.bXXXeXXX.fa
	qdirFASTA	   = my_dataset.getDirFasta(prefixQuery)    # Directory with FASTA files in the format: genomeName.chrName.bXXXeXXX.fa
	outDir		   = my_dataset.dirPCSs                     # Directory where PCSs will be saved (output directory).
	minPCSsize	   = my_dataset.minPCSsize                  # Minimum size (in number of base pairs) for the PCS to be considered.

	print("******************************************************")
	print("*        PCS (Perfectly conserved sequences)         *")
	print("* Extract PCSs from chains (ordered aligned blocks). *")
	print("******************************************************")
	print(f"- Reference genome: {prefixTarget}")
	print(f"- Query genome: {prefixQuery}")
	print(f"- Chain file: {chainsFilename}")
	print(f"- Directory for FASTA file of {prefixTarget}: {tdirFASTA}")
	print(f"- Directory for FASTA file of {prefixQuery}: {qdirFASTA}")
	print(f"- PCS directory (output): {outDir}")
	print(f"- Minimum size for PCS to be considered: >={minPCSsize} base pairs.")

	# Save times for each step and overall time.
	timeTrack = Time()
	timeTrack.start()

	##############################
	# Read header and details of all chains.
	# This step takes between 3 and 4 minutes.
	timeTrack.startStep("Reading chains")
	chains, blocks = readChains(chainsFilename)
	timeTrack.stopStep()

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

	timeTrack.startStep("Compatible chains")

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

	timeTrack.stopStep() # It takes around 4 minutes.
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

		outFilename = my_dataset.getOutFilename_extractPCS(prefixTarget,qChrom)

		##############################
		# Validation setup.
		# This step is important to make sure that PCSs are correctly extracted.
		pcsChecks = initializePCSchecks(cntUB=150,prefixTarget=prefixTarget,prefixQuery=prefixQuery)

		##############################
		# Read DNA from chains and extract PCSs.
		timeTrack.startStep(f"PCSs from {qChrom}")
		pcs_lst = []
		for chain in qChains[qChrom]:

			# At the end of an iteration, the chain should have all PCSs computed.
			tDNA = readDNA(os.path.join(tdirFASTA,chain.tChrom), prefixTarget, chain.tChrom, chain.tStrand, chain.tStart, chain.tEnd)
			qDNA = readDNA(os.path.join(qdirFASTA,chain.qChrom), prefixQuery,  chain.qChrom, chain.qStrand, chain.qStart, chain.qEnd)

			# Sanity check before continue.
			validSeqs = ((chain.tEnd-chain.tStart) == len(tDNA)) and ((chain.qEnd-chain.qStart) == len(qDNA))

			# Check blocks.
			if (not validSeqs):
				print("ERROR! Problem to read DNA segments from chain.")
				sys.exit()
			else:
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

		timeTrack.stopStep()

		print(f"\tTotal PCSs found in {qChrom}={len(pcs_lst)}")
		print(f"\tPCS validation: OK! [Stats: {';'.join([c.print() for c in pcsChecks])}]")

		##############################
		# Save PCSs.
		print("Save PCSs per chromosome...")
		with open(outFilename, 'wb') as pickleFile:
			pickle.dump(pcs_lst, pickleFile, protocol=pickle.HIGHEST_PROTOCOL)
		print(f"Pickle saved! ({outFilename})")

	timeTrack.stop()
	timeTrack.print()
	print("PCS: Done!")

