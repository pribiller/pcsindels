import os
from utils.basicTypes import CompRes

class Dataset:
	"""This class contains all needed information about a genome dataset.
	"""
	def __init__(self):

		# Reference species: Humans.
		# WARNING: In the UCSC pairwise alignments, hg38 corresponds to the **query** genome (even though here we call it *reference*).
		self.refsp_ucscName   = "hg38" 
		self.chromLst         = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]

		# 40 vertebrate species information.
		self.speciesUCSCnames = ['panPan3', 'panTro6', 'gorGor6', 'ponAbe3', 'papAnu4', 'macFas5', 'rhiRox1', 'chlSab2', 'nasLar1', 'rheMac10', 'calJac4', 'tarSyr2', 'micMur2', 'galVar1', 'mm39', 'oryCun2', 'rn7', 'vicPac2', 'bisBis1', 'felCat9', 'manPen1', 'bosTau9', 'canFam6', 'musFur1', 'neoSch1', 'equCab3', 'myoLuc2', 'susScr11', 'enhLutNer1', 'triMan1', 'macEug2', 'ornAna2', 'aptMan1', 'galGal6', 'thaSir1', 'aquChr2', 'melGal5', 'xenLae2', 'xenTro10', 'danRer11']
		self.divergenceTimes  = {'enhLutNer1': 94.0, 'melGal5': 319.0, 'equCab3': 94.0, 'thaSir1': 319.0, 'canFam6': 94.0, 'xenLae2': 352.0, 'danRer11': 429.0, 'myoLuc2': 94.0, 'gorGor6': 15.1, 'vicPac2': 94.0, 'panTro6': 12.1, 'nasLar1': 28.8, 'bisBis1': 94.0, 'triMan1': 99.0, 'aptMan1': 319.0, 'micMur2': 74.0, 'manPen1': 94.0, 'musFur1': 94.0, 'galVar1': 79.0, 'rn7': 87.0, 'xenTro10': 352.0, 'ponAbe3': 15.2, 'rheMac10': 28.8, 'ornAna2': 180.0, 'galGal6': 319.0, 'tarSyr2': 69.0, 'bosTau9': 94.0, 'macEug2': 160.0, 'papAnu4': 28.8, 'panPan3': 12.1, 'rhiRox1': 28.8, 'neoSch1': 94.0, 'mm39': 87.0, 'oryCun2': 87.0, 'aquChr2': 319.0, 'felCat9': 94.0, 'macFas5': 28.8, 'chlSab2': 28.8, 'calJac4': 43.0, 'susScr11': 94.0}
		self.commonNames      = {'enhLutNer1': 'Southern sea otter', 'melGal5': 'Turkey', 'equCab3': 'Horse', 'thaSir1': 'Garter snake', 'canFam6': 'Dog', 'xenLae2': 'African clawed frog', 'danRer11': 'Zebrafish', 'myoLuc2': 'Little brown bat', 'gorGor6': 'Gorilla', 'vicPac2': 'Alpaca', 'panTro6': 'Chimp', 'nasLar1': 'Proboscis Monkey', 'bisBis1': 'Bison', 'triMan1': 'Manatee', 'aptMan1': 'Brown kiwi', 'micMur2': 'Mouse lemur', 'manPen1': 'Chinese pangolin', 'musFur1': 'Ferret', 'galVar1': 'Malayan flying lemur', 'rn7': 'Rat', 'xenTro10': 'X. tropicalis', 'ponAbe3': 'Orangutan', 'rheMac10': 'Rhesus', 'ornAna2': 'Platypus', 'galGal6': 'Chicken', 'tarSyr2': 'Tarsier', 'bosTau9': 'Cow', 'macEug2': 'Wallaby', 'papAnu4': 'Baboon', 'panPan3': 'Bonobo', 'rhiRox1': 'Golden snub-nosed monkey', 'neoSch1': 'Hawaiian monk seal', 'mm39': 'Mouse', 'oryCun2': 'Rabbit', 'aquChr2': 'Golden eagle', 'felCat9': 'Cat', 'macFas5': 'Crab-eating macaque', 'chlSab2': 'Green monkey', 'calJac4': 'Marmoset', 'susScr11': 'Pig'}

		# Directories.
		self.dirChains   = "/bucket/MillerU/Priscila/ucsc/chains" # Compressed files.
		self.dirFasta    = "/bucket/MillerU/Priscila/fasta"
		self.dirPCSs     = "/flash/MillerU/Priscila/paper-validation/pcs"
		self.dirWindows  = "/flash/MillerU/Priscila/paper-validation/windows"
		self.dirSetupEvolTimes = "/flash/MillerU/Priscila/paper-validation/taus-setup"

		self.dirTemp    = "/flash/MillerU/Priscila/paper-validation/tmp" # Temporary area.
		self.dirLog     = "/flash/MillerU/Priscila/paper-validation/log" # Temporary area.

		# Key parameters of the method used for estimating evolutionary times.
		self.minPCSsize  = 5    # bps
		self.windowSize  = 1000 # bps
		self.alphas      = [1.1, 10]

		# Parameters for pre-computing evolutionary times (3_preCompEvolTimes.py).
		self.maxDiffSampSize = 2   # maximum difference in number of PCSs between reference sample size and sample size (nb. of PCSs) observed in real datasets.
		self.maxDiffWinSize  = 50  # maximum difference in the nb. of bps found in the window and the reference to compute estimate.
		self.maxNbTaus       = 500 # maximum number of taus for brute force.
		self.nbSamplesPerTau = 100 # number of PCS size distribution sampled for each tau.
		self.minDiffPcsSizes = 5   # minimum number of points in PCS size distribution graph.

		# Computational resources required depending on the job input (species/chromosome/etc.).
		self.computationalResources_extractPCS     = {"panPan3": CompRes("01:01:19", 5, 0.95),"panTro6": CompRes("00:59:09", 6, 0.98),"gorGor6": CompRes("00:50:02", 5, 1.14),"ponAbe3": CompRes("01:20:13", 5, 1.98),"papAnu4": CompRes("06:27:13", 38, 2.76),"macFas5": CompRes("02:42:13", 20, 2.71),"rhiRox1": CompRes("03:29:32", 23, 2.98),"chlSab2": CompRes("00:59:52", 6, 2.99),"nasLar1": CompRes("04:26:24", 15, 2.25),"rheMac10": CompRes("01:01:55", 6, 2.80),"calJac4": CompRes("02:50:19", 36, 3.14),"tarSyr2": CompRes("04:57:09", 42, 2.62),"micMur2": CompRes("02:24:55", 17, 2.25),"galVar1": CompRes("09:48:31", 40, 2.56),"mm39": CompRes("00:58:45", 9, 1.10),"oryCun2": CompRes("01:23:21", 19, 1.67),"rn7": CompRes("01:27:01", 24, 1.07),"vicPac2": CompRes("01:23:54", 17, 2.19),"bisBis1": CompRes("02:50:09", 24, 1.86),"felCat9": CompRes("01:57:19", 28, 2.19),"manPen1": CompRes("03:18:57", 18, 1.99),"bosTau9": CompRes("01:10:49", 15, 1.59),"canFam6": CompRes("01:00:37", 18, 2.07),"musFur1": CompRes("01:00:55", 17, 2.21),"neoSch1": CompRes("01:20:42", 19, 2.46),"equCab3": CompRes("01:30:38", 17, 2.60),"myoLuc2": CompRes("01:46:55", 20, 1.60),"susScr11": CompRes("01:20:55", 18, 1.93),"enhLutNer1": CompRes("01:27:55", 18, 2.26),"triMan1": CompRes("01:20:32", 20, 2.16),"macEug2": CompRes("05:58:51", 18, 0.22),"ornAna2": CompRes("01:20:04", 13, 0.17),"aptMan1": CompRes("00:45:15", 3, 0.12),"galGal6": CompRes("00:45:15", 3, 0.09),"thaSir1": CompRes("01:30:08", 5, 0.07),"aquChr2": CompRes("00:30:03", 3, 0.12),"melGal5": CompRes("01:15:00", 3, 0.09),"xenLae2": CompRes("01:15:46", 6, 0.06),"xenTro10": CompRes("02:40:30", 40, 0.06),"danRer11": CompRes("04:08:26", 14, 0.04)}
		self.computationalResources_computeWindows = {"chr1": CompRes("01:10:00", 11, 0.34),"chr2": CompRes("01:30:00", 12, 0.36),"chr3": CompRes("01:00:00", 8, 0.28),"chr4": CompRes("01:00:00", 8, 0.28),"chr5": CompRes("01:00:00", 8, 0.26),"chr6": CompRes("01:00:00", 7, 0.23),"chr7": CompRes("01:00:00", 6, 0.19),"chr8": CompRes("01:00:00", 7, 0.2),"chr9": CompRes("00:30:34", 5, 0.15),"chr10": CompRes("00:45:43", 6, 0.17),"chr11": CompRes("00:45:24", 5, 0.19),"chr12": CompRes("00:45:18", 6, 0.19),"chr13": CompRes("00:45:01", 5, 0.14),"chr14": CompRes("00:45:10", 4, 0.14),"chr15": CompRes("00:45:46", 4, 0.12),"chr16": CompRes("00:45:06", 4, 0.12),"chr17": CompRes("00:30:16", 4, 0.1),"chr18": CompRes("00:30:13", 4, 0.11),"chr19": CompRes("00:20:52", 3, 0.08),"chr20": CompRes("00:30:15", 4, 0.1),"chr21": CompRes("00:20:09", 3, 0.06),"chr22": CompRes("00:07:34", 2, 0.05),"chrX": CompRes("00:29:55", 4, 0.2),"chrY": CompRes("00:10:10", 2, 0.03)}

	def getChainFilename(self,ucscName,compressed=True):
		fileExtension    = "chain.gz" if compressed else "chain"
		chainFilename    = f"{ucscName}.{self.refsp_ucscName}.all.{fileExtension}"
		chainFilename_full = os.path.join(self.dirChains, chainFilename)
		chainFilename_tmp  = os.path.join(self.dirTemp,   chainFilename)
		if (not os.path.isfile(chainFilename_full)):
			if (not os.path.isfile(chainFilename_tmp)):
				print(f"ERROR! Chain file not found for species {ucscName}. Check if the file exists: {chainFilename}.\nDirectories where file was searched:\n\t{self.dirChains}\n\t{self.dirTemp}")
				chainFilename = None
			else:
				chainFilename = chainFilename_tmp
		else:
			chainFilename = chainFilename_full
		return chainFilename

	def getDirFasta(self,ucscName):
		fastaDir = os.path.join(self.dirFasta, ucscName, "chr")
		if (not os.path.isdir(fastaDir)):
			print(f"ERROR! Directory not found for FASTA files of species {ucscName}. Check if directory exists: {fastaDir}.")
			fastaDir = None
		return fastaDir

	###########################
	# Extract PCS.
	def getCompRes_extractPCS(self,ucscName):
		return self.computationalResources_extractPCS[ucscName]

	def getLogFilename_extractPCS(self,ucscName):
		return os.path.join(self.dirLog, f"extractPCS.{ucscName}.log")

	def getOutFilename_extractPCS(self,ucscName,chrom):
		return os.path.join(self.dirPCSs, f"{self.refsp_ucscName}.{ucscName}.{chrom}.pcs.pickle")

	def getLogFilename_extractPCS(self,ucscName):
		return os.path.join(self.dirLog, f"extractPCS.{ucscName}.log")

	###########################
	# Compute windows.
	def getCompRes_computeWindows(self,chrom):
		return self.computationalResources_computeWindows[chrom]

	def getLogFilename_computeWindows(self,chrom):
		return os.path.join(self.dirLog, f"computeWindows.{chrom}.log")

	def getTmpFilename_computeWindows(self,chrom):
		return os.path.join(self.dirTemp, f"{self.refsp_ucscName}.{chrom}.mergedPCSs.pickle")

	def getOutFilename_computeWindows(self,chrom):
		return os.path.join(self.dirWindows, f"{self.refsp_ucscName}.{chrom}.{self.windowSize}.windows.pickle")

	###########################
	# Pre-compute evol. times.
	def getOutFilename_preCompEvolTimes(self,windowSize,alpha):
		return os.path.join(self.dirSetupEvolTimes, f"taus-setup.alpha{alpha}.windowSize{windowSize}.pickle")
		