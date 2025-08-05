import os
import re
import glob
from utils.basicTypes import CompRes

class Dataset:
	"""This class contains all needed information about a genome dataset.
	"""
	def __init__(self):

		# Reference species: Humans.
		# WARNING: In the UCSC pairwise alignments, hg38 corresponds to the **query** genome (even though here we call it *reference*).
		self.refsp_ucscName   = "hg38" 
		self.chromLst         = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
		# Source: https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes
		self.chromSizes       = {"chr1" : 248956422, "chr2" : 242193529, "chr3" : 198295559, "chr4" : 190214555, "chr5" : 181538259, "chr6" : 170805979, "chr7" : 159345973, "chr8" : 145138636, "chr9" : 138394717, "chr10" : 133797422, "chr11" : 135086622, "chr12" : 133275309, "chr13" : 114364328, "chr14" : 107043718, "chr15" : 101991189, "chr16" : 90338345, "chr17" : 83257441, "chr18" : 80373285, "chr19" : 58617616, "chr20" : 64444167, "chr21" : 46709983, "chr22" : 50818468, "chrX" : 156040895, "chrY" : 57227415}

		# 40 vertebrate species information.
		self.speciesUCSCnames = ['panPan3', 'panTro6', 'gorGor6', 'ponAbe3', 'papAnu4', 'macFas5', 'rhiRox1', 'chlSab2', 'nasLar1', 'rheMac10', 'calJac4', 'tarSyr2', 'micMur2', 'galVar1', 'mm39', 'oryCun2', 'rn7', 'vicPac2', 'bisBis1', 'felCat9', 'manPen1', 'bosTau9', 'canFam6', 'musFur1', 'neoSch1', 'equCab3', 'myoLuc2', 'susScr11', 'enhLutNer1', 'triMan1', 'macEug2', 'ornAna2', 'aptMan1', 'galGal6', 'thaSir1', 'aquChr2', 'melGal5', 'xenLae2', 'xenTro10', 'danRer11']
		self.divergenceTimes  = {'enhLutNer1': 94.0, 'melGal5': 319.0, 'equCab3': 94.0, 'thaSir1': 319.0, 'canFam6': 94.0, 'xenLae2': 352.0, 'danRer11': 429.0, 'myoLuc2': 94.0, 'gorGor6': 15.1, 'vicPac2': 94.0, 'panTro6': 12.1, 'nasLar1': 28.8, 'bisBis1': 94.0, 'triMan1': 99.0, 'aptMan1': 319.0, 'micMur2': 74.0, 'manPen1': 94.0, 'musFur1': 94.0, 'galVar1': 79.0, 'rn7': 87.0, 'xenTro10': 352.0, 'ponAbe3': 15.2, 'rheMac10': 28.8, 'ornAna2': 180.0, 'galGal6': 319.0, 'tarSyr2': 69.0, 'bosTau9': 94.0, 'macEug2': 160.0, 'papAnu4': 28.8, 'panPan3': 12.1, 'rhiRox1': 28.8, 'neoSch1': 94.0, 'mm39': 87.0, 'oryCun2': 87.0, 'aquChr2': 319.0, 'felCat9': 94.0, 'macFas5': 28.8, 'chlSab2': 28.8, 'calJac4': 43.0, 'susScr11': 94.0}
		self.commonNames      = {'enhLutNer1': 'Southern sea otter', 'melGal5': 'Turkey', 'equCab3': 'Horse', 'thaSir1': 'Garter snake', 'canFam6': 'Dog', 'xenLae2': 'African clawed frog', 'danRer11': 'Zebrafish', 'myoLuc2': 'Little brown bat', 'gorGor6': 'Gorilla', 'vicPac2': 'Alpaca', 'panTro6': 'Chimp', 'nasLar1': 'Proboscis Monkey', 'bisBis1': 'Bison', 'triMan1': 'Manatee', 'aptMan1': 'Brown kiwi', 'micMur2': 'Mouse lemur', 'manPen1': 'Chinese pangolin', 'musFur1': 'Ferret', 'galVar1': 'Malayan flying lemur', 'rn7': 'Rat', 'xenTro10': 'X. tropicalis', 'ponAbe3': 'Orangutan', 'rheMac10': 'Rhesus', 'ornAna2': 'Platypus', 'galGal6': 'Chicken', 'tarSyr2': 'Tarsier', 'bosTau9': 'Cow', 'macEug2': 'Wallaby', 'papAnu4': 'Baboon', 'panPan3': 'Bonobo', 'rhiRox1': 'Golden snub-nosed monkey', 'neoSch1': 'Hawaiian monk seal', 'mm39': 'Mouse', 'oryCun2': 'Rabbit', 'aquChr2': 'Golden eagle', 'felCat9': 'Cat', 'macFas5': 'Crab-eating macaque', 'chlSab2': 'Green monkey', 'calJac4': 'Marmoset', 'susScr11': 'Pig'}

		# Directories.
		self.dirChains   = "/bucket/MillerU/Priscila/ucsc/chains" # Compressed files.
		self.dirFasta    = "/bucket/MillerU/Priscila/fasta"
		self.dirPCSs     = "/bucket/MillerU/Priscila/Biller2025-MBE/pcs"
		self.dirWindows  = "/bucket/MillerU/Priscila/Biller2025-MBE/windows"
		self.dirSetupEvolTimes = "/bucket/MillerU/Priscila/Biller2025-MBE/taus-setup"
		self.dirEstEvolTimes   = "/bucket/MillerU/Priscila/Biller2025-MBE/taus-ests"
		self.dirSampEvolTimes  = "/bucket/MillerU/Priscila/Biller2025-MBE/taus-samps" #"/bucket/MillerU/Priscila/Biller2025-MBE/taus-samps"
		self.dirPlots    = "/flash/MillerU/Priscila/Biller2025-MBE/plots"
		self.dirIcons    = os.path.join(os.path.dirname(os.path.dirname(__file__)), "assets", "img",   "icons-species")
		self.dirTrees    = os.path.join(os.path.dirname(os.path.dirname(__file__)), "assets", "trees")
		self.dirAnnot    = os.path.join(os.path.dirname(os.path.dirname(__file__)), "assets", "annotation")

		self.dirTemp    = "/flash/MillerU/Priscila/Biller2025-MBE/tmp" # Temporary area.
		self.dirLog     = "/flash/MillerU/Priscila/Biller2025-MBE/log" # Temporary area.

		# Key parameters of the method used for estimating evolutionary times.
		self.minPCSsize  = 5    # bps
		self.windowSize  = 1000 # bps
		self.alphas      = [1.1, 10]

		# Parameters for pre-computing evolutionary times (3_setupEvolTimes.py).
		self.maxDiffSampSize = 2   # maximum difference in number of PCSs between reference sample size and sample size (nb. of PCSs) observed in real datasets.
		self.maxDiffWinSize  = 50  # maximum difference in the nb. of bps found in the window and the reference to compute estimate.
		self.maxNbTaus       = 500 # maximum number of taus for brute force.
		self.nbSamplesPerTau = 100 # number of PCS size distribution sampled for each tau.
		self.minDiffPcsSizes = 5   # minimum number of points in PCS size distribution graph.
		self.minSampRefSize  = 5   # minimum number of PCSs for reference sample size.
		self.maxSampRefSize  = 150 # maximum number of PCSs for reference sample size.
		
		# Parameters for sampling evolutionary times.
		self.nbSamplesPerWin = 10  # number of PCS size distribution sampled for each window.
		
		# Parameters for plotting results.
		self.annotationCategories={
			(1,1,1) : ("Coding","Protein",""),

			(1,2,1) : ("Coding","Long noncoding RNAs (LncRNA)",""),
			(1,2,2) : ("Coding","Non-protein","Micro RNAs (miRNA)"),
			(1,2,3) : ("Coding","Non-protein","Small nuclear RNAs (snRNA)"),
			(1,2,4) : ("Coding","Non-protein","Small nucleolar RNAs (snoRNAs)"),
			(1,2,5) : ("Coding","Non-protein","Other RNAs"),

			(1,3,1) : ("Coding","Pseudogene",""),

			(2,1,1) : ("Gene regulatory regions","Enhancers",""),
			(2,1,2) : ("Gene regulatory regions","Silencers",""),
			(2,1,3) : ("Gene regulatory regions","DNase I hypersensitive sites",""),
			(2,1,4) : ("Gene regulatory regions","Insulators",""),
			(2,1,5) : ("Gene regulatory regions","Locus control regions",""),

			(2,2,1) : ("Gene regulatory regions", "Introns", ""),
			(2,2,2) : ("Gene regulatory regions", "UTRs (untranslated regions)",""),

			(3,1,1) : ("Repeats","Long interspersed nuclear elements (LINE)",""),
			(3,1,2) : ("Repeats","Long interspersed nuclear elements (LINE)","CR1"),
			(3,1,3) : ("Repeats","Long interspersed nuclear elements (LINE)","L1"),
			(3,1,4) : ("Repeats","Long interspersed nuclear elements (LINE)","L2"),

			(3,2,1) : ("Repeats","Short interspersed nuclear elements (SINE)",""),
			(3,2,2) : ("Repeats","Short interspersed nuclear elements (SINE)","ALU"),
			(3,2,3) : ("Repeats","Short interspersed nuclear elements (SINE)","MIR"),

			(3,3,1) : ("Repeats","DNA repeat elements (DNA)",""),
			(3,3,2) : ("Repeats","DNA repeat elements (DNA)","TcMar"),
			(3,3,3) : ("Repeats","DNA repeat elements (DNA)","hAT"),

			(3,4,1) : ("Repeats","RNA repeats",""),

			(3,5,1) : ("Repeats","Long terminal repeat elements (LTR)",""),
			(3,5,2) : ("Repeats","Long terminal repeat elements (LTR)","ERV1"),
			(3,5,3) : ("Repeats","Long terminal repeat elements (LTR)","ERVK"),
			(3,5,4) : ("Repeats","Long terminal repeat elements (LTR)","ERVL_MaLR"),
			(3,5,5) : ("Repeats","Long terminal repeat elements (LTR)","ERVL"),

			(3,6,1) : ("Repeats","Satellite repeats",""),
			(3,7,1) : ("Repeats","Simple repeats (micro-satellites)",""),
			(3,8,1) : ("Repeats","Low complexity repeats",""),
			(3,9,1) : ("Repeats","Others",""),

			(10,1,1): ("Unannotated", "", "")
		}
		self.unannotatedCategory = (10,1,1)
		# A functional class must cover over 90% of the window for
		# the window to be annotated with that class.
		self.thresholdAnnot      = 0.9

		# Computational resources required depending on the job input (species/chromosome/etc.).
		# Memory and time needed to run each job in 1 core. One job per species (40 jobs in total).
		self.computationalResources_extractPCS        = {"panPan3": CompRes("01:01:19", 5, 0.95),"panTro6": CompRes("00:59:09", 6, 0.98),"gorGor6": CompRes("00:50:02", 5, 1.14),"ponAbe3": CompRes("01:20:13", 5, 1.98),"papAnu4": CompRes("06:27:13", 38, 2.76),"macFas5": CompRes("02:42:13", 20, 2.71),"rhiRox1": CompRes("03:29:32", 23, 2.98),"chlSab2": CompRes("00:59:52", 6, 2.99),"nasLar1": CompRes("04:26:24", 15, 2.25),"rheMac10": CompRes("01:01:55", 6, 2.80),"calJac4": CompRes("02:50:19", 36, 3.14),"tarSyr2": CompRes("04:57:09", 42, 2.62),"micMur2": CompRes("02:24:55", 17, 2.25),"galVar1": CompRes("09:48:31", 40, 2.56),"mm39": CompRes("00:58:45", 9, 1.10),"oryCun2": CompRes("01:23:21", 19, 1.67),"rn7": CompRes("01:27:01", 24, 1.07),"vicPac2": CompRes("01:23:54", 17, 2.19),"bisBis1": CompRes("02:50:09", 24, 1.86),"felCat9": CompRes("01:57:19", 28, 2.19),"manPen1": CompRes("03:18:57", 18, 1.99),"bosTau9": CompRes("01:10:49", 15, 1.59),"canFam6": CompRes("01:00:37", 18, 2.07),"musFur1": CompRes("01:00:55", 17, 2.21),"neoSch1": CompRes("01:20:42", 19, 2.46),"equCab3": CompRes("01:30:38", 17, 2.60),"myoLuc2": CompRes("01:46:55", 20, 1.60),"susScr11": CompRes("01:20:55", 18, 1.93),"enhLutNer1": CompRes("01:27:55", 18, 2.26),"triMan1": CompRes("01:20:32", 20, 2.16),"macEug2": CompRes("05:58:51", 18, 0.22),"ornAna2": CompRes("01:20:04", 13, 0.17),"aptMan1": CompRes("00:45:15", 3, 0.12),"galGal6": CompRes("00:45:15", 3, 0.09),"thaSir1": CompRes("01:30:08", 5, 0.07),"aquChr2": CompRes("00:30:03", 3, 0.12),"melGal5": CompRes("01:15:00", 3, 0.09),"xenLae2": CompRes("01:15:46", 6, 0.06),"xenTro10": CompRes("02:40:30", 40, 0.06),"danRer11": CompRes("04:08:26", 14, 0.04)}
		# Memory and time needed to run each job in 1 core. One job per chromosome (24 jobs in total).
		self.computationalResources_computeWindows    = {"chr1": CompRes("01:10:00", 11, 0.34),"chr2": CompRes("01:30:00", 12, 0.36),"chr3": CompRes("01:00:00", 8, 0.28),"chr4": CompRes("01:00:00", 8, 0.28),"chr5": CompRes("01:00:00", 8, 0.26),"chr6": CompRes("01:00:00", 7, 0.23),"chr7": CompRes("01:00:00", 6, 0.19),"chr8": CompRes("01:00:00", 7, 0.2),"chr9": CompRes("00:30:34", 5, 0.15),"chr10": CompRes("00:45:43", 6, 0.17),"chr11": CompRes("00:45:24", 5, 0.19),"chr12": CompRes("00:45:18", 6, 0.19),"chr13": CompRes("00:45:01", 5, 0.14),"chr14": CompRes("00:45:10", 4, 0.14),"chr15": CompRes("00:45:46", 4, 0.12),"chr16": CompRes("00:45:06", 4, 0.12),"chr17": CompRes("00:30:16", 4, 0.1),"chr18": CompRes("00:30:13", 4, 0.11),"chr19": CompRes("00:20:52", 3, 0.08),"chr20": CompRes("00:30:15", 4, 0.1),"chr21": CompRes("00:20:09", 3, 0.06),"chr22": CompRes("00:07:34", 2, 0.05),"chrX": CompRes("00:29:55", 4, 0.2),"chrY": CompRes("00:10:10", 2, 0.03)}
		# Memory and time needed to run each job in 80 cores. One job per alpha (2 jobs in total).
		self.computationalResources_setupEvolTimes    = CompRes("3:00:00", 35, 6.290)
		# Memory and time needed to run each job in 30 cores. One job per species/alpha (80 jobs in total).
		self.computationalResources_estimateEvolTimes = {"panPan3": CompRes("02:32:16", 27, 4.03),"panTro6": CompRes("02:37:07", 29, 4.16),"gorGor6": CompRes("02:32:49", 31, 4.15),"ponAbe3": CompRes("02:42:49", 33, 3.88),"papAnu4": CompRes("02:32:03", 33, 3.62),"macFas5": CompRes("01:55:55", 33, 3.59),"rhiRox1": CompRes("02:34:24", 31, 3.83),"chlSab2": CompRes("02:31:06", 33, 3.85),"nasLar1": CompRes("02:34:28", 31, 3.11),"rheMac10": CompRes("02:20:27", 32, 3.64),"calJac4": CompRes("02:32:40", 27, 3.31),"tarSyr2": CompRes("02:12:27", 25, 3.24),"micMur2": CompRes("02:12:30", 25, 2.69),"galVar1": CompRes("02:12:18", 25, 3.22),"mm39": CompRes("02:12:46", 24, 2.14),"oryCun2": CompRes("02:12:04", 25, 2.57),"rn7": CompRes("02:12:31", 24, 2.12),"vicPac2": CompRes("02:12:22", 25, 2.99),"bisBis1": CompRes("02:12:00", 23, 2.71),"felCat9": CompRes("01:58:46", 25, 3.0),"manPen1": CompRes("02:32:00", 24, 2.82),"bosTau9": CompRes("01:40:54", 24, 2.26),"canFam6": CompRes("02:06:19", 25, 2.87),"musFur1": CompRes("02:34:22", 25, 3.08),"neoSch1": CompRes("02:36:12", 26, 3.22),"equCab3": CompRes("02:36:21", 25, 3.25),"myoLuc2": CompRes("02:24:06", 24, 2.29),"susScr11": CompRes("02:31:27", 25, 2.72),"enhLutNer1": CompRes("02:33:34", 25, 3.08),"triMan1": CompRes("02:31:40", 25, 2.99),"macEug2": CompRes("01:55:15", 20, 0.73),"ornAna2": CompRes("01:55:16", 18, 0.54),"aptMan1": CompRes("01:35:15", 18, 0.41),"galGal6": CompRes("01:55:10", 18, 0.33),"thaSir1": CompRes("01:52:28", 18, 0.26),"aquChr2": CompRes("01:55:58", 16, 0.39),"melGal5": CompRes("01:55:43", 19, 0.34),"xenLae2": CompRes("01:55:43", 18, 0.24),"xenTro10": CompRes("01:56:51", 19, 0.23),"danRer11": CompRes("01:56:48", 20, 0.17)}
		# Memory and time needed to run each job in 80 cores. One job per species/alpha (80 jobs in total).
		self.computationalResources_sampleEvolTimes   = {"panPan3": CompRes("08:18:15", 26, 0.067),"panTro6": CompRes("08:22:54", 14, 0.069),"gorGor6": CompRes("08:36:19", 15, 0.067),"ponAbe3": CompRes("08:12:14", 12, 0.06),"papAnu4": CompRes("08:06:18", 12, 0.056),"macFas5": CompRes("07:42:12", 10, 0.056),"rhiRox1": CompRes("08:06:59", 10, 0.059),"chlSab2": CompRes("08:22:47", 10, 0.059),"nasLar1": CompRes("07:40:08", 12, 0.048),"rheMac10": CompRes("08:07:51", 10, 0.056),"calJac4": CompRes("06:59:01", 13, 0.05),"tarSyr2": CompRes("06:52:29", 16, 0.049),"micMur2": CompRes("05:31:45", 11, 0.041),"galVar1": CompRes("06:31:39", 11, 0.048),"mm39": CompRes("04:49:00", 9, 0.033),"oryCun2": CompRes("05:34:08", 9, 0.039),"rn7": CompRes("04:35:49", 8, 0.033),"vicPac2": CompRes("06:32:32", 10, 0.045),"bisBis1": CompRes("05:43:21", 9, 0.041),"felCat9": CompRes("06:25:28", 10, 0.045),"manPen1": CompRes("06:21:07", 13, 0.043),"bosTau9": CompRes("04:58:30", 9, 0.035),"canFam6": CompRes("06:07:25", 10, 0.043),"musFur1": CompRes("06:28:23", 10, 0.046),"neoSch1": CompRes("06:58:50", 11, 0.048),"equCab3": CompRes("06:57:51", 11, 0.049),"myoLuc2": CompRes("04:58:24", 9, 0.035),"susScr11": CompRes("05:43:36", 10, 0.041),"enhLutNer1": CompRes("06:20:01", 10, 0.046),"triMan1": CompRes("06:36:52", 10, 0.045),"macEug2": CompRes("01:59:04", 8, 0.012),"ornAna2": CompRes("01:30:12", 8, 0.009),"aptMan1": CompRes("01:15:55", 8, 0.007),"galGal6": CompRes("01:08:40", 7, 0.006),"thaSir1": CompRes("00:58:08", 9, 0.005),"aquChr2": CompRes("01:13:07", 10, 0.007),"melGal5": CompRes("01:09:33", 6, 0.006),"xenLae2": CompRes("00:56:29", 7, 0.005),"xenTro10": CompRes("00:58:07", 7, 0.005),"danRer11": CompRes("00:48:06", 7, 0.003)}

		# WARNING: Only disk space was checked for this script.
		self.computationalResources_breakFasta        = {"panPan3": CompRes("23:59:59", 50, 2.84), "panTro6": CompRes("23:59:59", 50, 2.84), "gorGor6": CompRes("23:59:59", 50, 2.84), "ponAbe3": CompRes("23:59:59", 50, 2.85), "papAnu4": CompRes("23:59:59", 50, 2.76), "macFas5": CompRes("23:59:59", 50, 2.74), "rhiRox1": CompRes("23:59:59", 50, 2.70), "chlSab2": CompRes("23:59:59", 50, 2.60), "nasLar1": CompRes("23:59:59", 50, 2.81), "rheMac10": CompRes("23:59:59", 50, 2.77), "calJac4": CompRes("23:59:59", 50, 2.70), "tarSyr2": CompRes("23:59:59", 50, 3.22), "micMur2": CompRes("23:59:59", 50, 2.27), "galVar1": CompRes("23:59:59", 50, 2.97), "mm39": CompRes("23:59:59", 50, 2.54), "oryCun2": CompRes("23:59:59", 50, 2.55), "rn7": CompRes("23:59:59", 50, 2.47), "vicPac2": CompRes("23:59:59", 50, 2.02), "bisBis1": CompRes("23:59:59", 50, 2.75), "felCat9": CompRes("23:59:59", 50, 2.35), "manPen1": CompRes("23:59:59", 50, 2.05), "bosTau9": CompRes("23:59:59", 50, 2.53), "canFam6": CompRes("23:59:59", 50, 2.15), "musFur1": CompRes("23:59:59", 50, 2.25), "neoSch1": CompRes("23:59:59", 50, 2.24), "equCab3": CompRes("23:59:59", 50, 2.33), "myoLuc2": CompRes("23:59:59", 50, 1.89), "susScr11": CompRes("23:59:59", 50, 2.33), "enhLutNer1": CompRes("23:59:59", 50, 2.26), "triMan1": CompRes("23:59:59", 50, 2.89), "macEug2": CompRes("23:59:59", 50, 2.86), "ornAna2": CompRes("23:59:59", 50, 1.86), "aptMan1": CompRes("23:59:59", 50, 1.42), "galGal6": CompRes("23:59:59", 50, 0.99), "thaSir1": CompRes("23:59:59", 50, 1.33), "aquChr2": CompRes("23:59:59", 50, 1.11), "melGal5": CompRes("23:59:59", 50, 1.05), "xenLae2": CompRes("23:59:59", 50, 2.53), "xenTro10": CompRes("23:59:59", 50, 1.35), "danRer11": CompRes("23:59:59", 50, 1.56)}

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
	def getCompRes_setupEvolTimes(self):
		return self.computationalResources_setupEvolTimes

	def getLogFilename_setupEvolTimes(self,alpha):
		return os.path.join(self.dirLog, f"setupEvolTimes.alpha{float(alpha)}.log")

	def getOutFilename_setupEvolTimes(self,windowSizeRef,alpha):
		return os.path.join(self.dirSetupEvolTimes, f"evolTimes-setup.alpha{float(alpha)}.winSize{windowSizeRef}.pickle")
	
	def getOutFilenamePattern_setupEvolTimes(self,alpha):
		return os.path.join(self.dirSetupEvolTimes, f"evolTimes-setup.alpha{float(alpha)}.winSize*.pickle")

	def getOutFilenames_setupEvolTimes(self,alpha):
		""" This function retrieves all setup files that 
		were computed so far for a given alpha.

		:returns: A dictionary where keys are reference window sizes, 
			and values are filepaths where the setup information
			on a reference window size can be found.
		"""
		setupFilenames = glob.glob(self.getOutFilenamePattern_setupEvolTimes(float(alpha)))
		winSizeRef_all = {}
		for fileName in setupFilenames:
			m = re.match(r".*winSize(\d+).pickle$", fileName)
			if m: 
				winRefSize = int(m.group(1))
				winSizeRef_all[winRefSize] = fileName
		return winSizeRef_all

	###########################
	# Estimate evol. times.
	def getCompRes_estimateEvolTimes(self,ucscName):
		return self.computationalResources_estimateEvolTimes[ucscName]

	def getLogFilename_estimateEvolTimes(self,ucscName,alpha):
		return os.path.join(self.dirLog, f"estimateEvolTimes.alpha{float(alpha)}.{ucscName}.log")

	def getOutFilenamePattern_estimateEvolTimes(self,ucscName,alpha):
		return os.path.join(self.dirEstEvolTimes, f"evolTimes-ests.{ucscName}.*.alpha{float(alpha)}.pickle")

	def getOutFilename_estimateEvolTimes(self, ucscName, chrom, alpha):
		return os.path.join(self.dirEstEvolTimes, f"evolTimes-ests.{ucscName}.{chrom}.alpha{float(alpha)}.pickle")

	###########################
	# Sample evol. times.
	def getCompRes_sampleEvolTimes(self,ucscName):
		return self.computationalResources_sampleEvolTimes[ucscName]

	def getLogFilename_sampleEvolTimes(self,ucscName,alpha):
		return os.path.join(self.dirLog, f"sampleEvolTimes.alpha{float(alpha)}.{ucscName}.log")

	def getOutFilenamePattern_sampleEvolTimes(self,ucscName,alpha):
		return self.getOutFilename_sampleEvolTimes(ucscName, alpha)

	def getOutFilename_sampleEvolTimes(self, ucscName, alpha):
		return os.path.join(self.dirSampEvolTimes, f"pcsDistrib-samp.{ucscName}.alpha{float(alpha)}.pickle")

	###########################
	# Plot figures.

	# Figure 2 + supp related figures.
	def getOutFilename_plot_PcsDistribComp(self, alpha, isSuppFig, fileExtension):
		return os.path.join(self.dirPlots, f"pcsDistrib-comp{('-supp' if isSuppFig else '')}.alpha{float(alpha)}.{fileExtension}")

	# Figure 3.
	def getOutFilename_plot_EvolTimesComp(self, paperRef, corr=False):
		return os.path.join(self.dirPlots, f"evolTimes-comp{paperRef}{('-corr' if corr else '')}.pdf")

	# Figure 5 + supp related figures.
	def getOutFilename_annotatedWindows(self):
		return os.path.join(self.dirAnnot, f"{self.refsp_ucscName}.{self.windowSize}.annotations.pickle")

	# Supplementary figure (Fig. S7 - Functional Classes).
	def getOutFilename_plot_FuncClasses(self, annotOverlapThresh,matrix_type):
		return os.path.join(self.dirPlots, f"funcClasses-stats.{matrix_type}.overlapThresh{annotOverlapThresh}.pdf")
