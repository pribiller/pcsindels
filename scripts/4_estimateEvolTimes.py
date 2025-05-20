"""Reads all PCS size distributions across the windows of a designated chromosome in 
a given species (stored in one of the outputs files of ``2_computeWindows.py``)
and, along with the setup files produced by the script ``3_setupEvolTimes.py``, 
estimates the evolutionary times for each window.

- **Use**::

	python3 4_estimateEvolTimes.py -sp_ucsc_name [UCSC name] -chr [chromosome name] -alpha [any number > 1] -cores [nb. of cores] [--overwrite, optional]

- **Example of Usage (human (reference genome) and mouse, all windows mapping to chromosome 16 in humans)**::

	python3 ~/code/4_estimateEvolTimes.py -sp_ucsc_name mm39 -chr chr16 -alpha 1.1 -cores 1 --overwrite

- **Input Parameter (mandatory)**:

:-sp_ucsc_name: UCSC name of the species that is being aligned with the reference species (e.g. *mm39* for mouse).
:-chr: 			Chromosome name in the reference species (e.g. chr1, chr2, ..., chrX, chrY).
:-alpha: 		It determines how frequent longer indels can occur. 
				The parameter alpha can take any value above 1 (1 is **not** included): (1, ∞).
				If alpha is near 1, larger indels are more likely to occur.
				If alpha is above 5 (a hard upper limit set internally with ``max_alpha``), 
				the model turns into a substitution only model.

:-cores: 		If ``cores=1``, the script will execute serially, which may take several hours to complete. It is advisable to utilize as many available cores as possible.
:--overwrite: 	Optional flag. If this flag is specified, any existing output files will be overwritten during the run.

- **Other Parameters taken from** ``dataset.py``:

:refsp_ucscName: UCSC name of the reference species that is being aligned (e.g. *hg38* for human).
:dirWindows: 		Directory where windows and their PCS size distributions are saved (input files).
:dirSetupEvolTimes: Directory where setup files are saved (input files).
:dirEstEvolTimes: Directory where estimates of evolutionary time will be saved (output files).

.. note::
	Make sure that the required parameters described above are correctly defined in the file ``utils/dataset.py``.

- **Output**: 

	One ``.pickle`` file is created as output, detailing the posterior 
	evolutionary times for each window. All windows mapped to the specified 
	chromosome of the reference genome (*human* in our study) are included. 
	The PCS size distributions from these windows, representing perfectly conserved 
	sequences between the reference genome (e.g. ``hg38``) and the species 
	specified in the input (e.g. ``mm39``), are used as input for the 
	evolutionary time estimation method.
	
	The output ``.pickle`` file is structured as follows:

	1. 	A dictionary where each key represents a window identifier, 
		defined as a tuple containing the start and end coordinates 
		within the chromosome of the reference genome. The corresponding 
		values are NumPy arrays of size ``n``, where each index ``i`` 
		corresponds to an evolutionary time, and the value at that index 
		indicates the likelihood of the window being at evolutionary time ``i``.
	2.	A dictionary that maps window identifiers to their related 
		reference window sizes and reference sample values.
	3.	A dictionary that maps reference window sizes to selected evolutionary times.

	.. note::
		If you are interested in finding the evolutionary time at position ``i`` for a window, 
		first check the reference window size in the second dictionary, then use it to get 
		the evolutionary time values from the third dictionary. This indirect way of accessing information
		avoids saving the same data multiple times.

Pre-requisites
--------------

	Before using this script, make sure all the required files were pre-computed:

a) Window file for the chromosome specified in the input
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Make sure to run ``2_computeWindows.py`` for the chromosome specified in the 
input parameter ``-chr``.

b) Setup files for the estimation method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Make sure to run ``3_setupEvolTimes.py`` for the alpha value specified in the 
input parameter ``-alpha``.


Cluster resources
-----------------
In case you want to run this Python script stand-alone in a cluster that uses Slurm to manage jobs::

	srun -p compute -t 2:00:00 --mem 20G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash

Otherwise you can use the script ``../cluster/4_estimateEvolTimes_runAll.py``
to run this script for all 40 species used in this study (whole genome).

Time, Memory & Disk space
-------------------------

For reference, here we include a run example, with runtime, memory usage, and disk 
space required for running this script on each pairwise alignment of the 40 
vertebrate dataset examined in our study. **30 cores** were used in these runs.

======================  =========  =========  ========  =========  =========  ========
Desc.                         Parameter α=1.1                 Parameter α=10          
----------------------  ------------------------------  ------------------------------
UCSC name               Time       Memory     Disk      Time       Memory     Disk    
======================  =========  =========  ========  =========  =========  ========
panPan3                 02:12:16   27GB       4.030GB   02:08:24   22GB       3.710GB 
panTro6                 02:07:07   29GB       4.160GB   02:08:34   23GB       3.830GB 
gorGor6                 02:12:49   31GB       4.150GB   02:14:49   26GB       3.810GB 
ponAbe3                 02:22:49   33GB       3.880GB   02:08:45   25GB       3.560GB 
papAnu4                 02:02:03   33GB       3.620GB   02:04:07   27GB       3.320GB 
macFas5                 01:50:55   33GB       3.590GB   01:44:34   26GB       3.300GB 
rhiRox1                 02:04:24   31GB       3.830GB   02:06:27   26GB       3.520GB 
chlSab2                 02:11:06   33GB       3.850GB   01:58:09   27GB       3.530GB 
nasLar1                 02:04:28   31GB       3.110GB   01:58:45   25GB       2.860GB 
rheMac10                01:59:27   32GB       3.640GB   02:08:37   27GB       3.340GB 
calJac4                 02:02:40   27GB       3.310GB   02:00:49   24GB       3.040GB 
tarSyr2                 01:52:27   25GB       3.240GB   01:57:15   21GB       2.970GB 
micMur2                 01:56:30   25GB       2.690GB   01:51:04   20GB       2.470GB 
galVar1                 01:59:18   25GB       3.220GB   02:00:21   21GB       2.950GB 
mm39                    01:38:46   24GB       2.140GB   01:47:34   18GB       1.960GB 
oryCun2                 01:53:04   25GB       2.570GB   01:48:55   19GB       2.360GB 
rn7                     01:56:31   24GB       2.120GB   01:45:28   18GB       1.940GB 
vicPac2                 01:58:22   25GB       2.990GB   01:58:28   20GB       2.740GB 
bisBis1                 01:43:00   23GB       2.710GB   01:56:04   20GB       2.480GB 
felCat9                 01:28:46   25GB       3.000GB   01:56:15   20GB       2.750GB 
manPen1                 02:02:00   24GB       2.820GB   01:53:13   20GB       2.580GB 
bosTau9                 01:40:54   24GB       2.260GB   01:46:43   19GB       2.070GB 
canFam6                 02:06:19   25GB       2.870GB   01:57:34   20GB       2.630GB 
musFur1                 02:04:22   25GB       3.080GB   01:57:49   20GB       2.830GB 
neoSch1                 02:06:12   26GB       3.220GB   01:59:31   20GB       2.960GB 
equCab3                 02:06:21   25GB       3.250GB   01:59:49   21GB       2.980GB 
myoLuc2                 01:54:06   24GB       2.290GB   01:29:52   19GB       2.100GB 
susScr11                02:01:27   25GB       2.720GB   01:41:52   19GB       2.500GB 
enhLutNer1              02:03:34   25GB       3.080GB   01:54:00   20GB       2.830GB 
triMan1                 02:01:40   25GB       2.990GB   01:46:48   20GB       2.740GB 
macEug2                 01:40:15   20GB       0.730GB   01:24:21   16GB       0.660GB 
ornAna2                 01:40:16   18GB       0.540GB   01:30:07   16GB       0.500GB 
aptMan1                 01:35:15   18GB       0.410GB   01:34:37   16GB       0.380GB 
galGal6                 01:40:10   18GB       0.330GB   01:37:57   13GB       0.300GB 
thaSir1                 01:42:28   18GB       0.260GB   01:27:48   15GB       0.240GB 
aquChr2                 01:35:58   16GB       0.390GB   01:27:22   15GB       0.350GB 
melGal5                 01:40:43   19GB       0.340GB   01:19:10   15GB       0.310GB 
xenLae2                 01:40:43   18GB       0.240GB   01:22:27   13GB       0.220GB 
xenTro10                01:36:51   19GB       0.230GB   01:30:52   15GB       0.210GB 
danRer11                01:36:48   20GB       0.170GB   01:30:24   15GB       0.160GB 
======================  =========  =========  ========  =========  =========  ========

Time per Run: Details
^^^^^^^^^^^^^^^^^^^^^

Stats on time of a single run (human-mouse alignment, chromosome 16, 1 core): **~5 minutes**
Details on computational time are available in the log of the run.

====================  =======
Step                  Time (s)
====================  =======
Read windows            14.10
Read setup              91.50
Compute evol. times    135.49
**Total time**         241.34
====================  =======

Storage per Run: Details
^^^^^^^^^^^^^^^^^^^^^^^^

Size of output files with PCSs for the *human* and *mouse* 
alignment, chromosome 16: **~61 MB**.

Function details
----------------

Only relevant functions have been documented below. 
For more details on any function, check the comments in the souce code.

"""

import os
import sys
import pickle
from collections import defaultdict
from concurrent import futures
import numpy as np
from mpmath import *
import argparse

from utils.indelsModel import IndelsModel, EvolTimesSolver, getSampleSize
from utils.basicTypes import Time
from utils.dataset import Dataset

mp.dps = 30
mp.pretty = True

def getSampleSizeRef(sampleSizeRef_lst, pcsSizeDistrib, maxDiffSampSize):
	sampleSizeObsAdj = min(getSampleSize(pcsSizeDistrib), max(sampleSizeRef_lst))
	sampleSizeRef    = min(sampleSizeRef_lst, key=lambda val: abs(val-sampleSizeObsAdj))
	if(abs(sampleSizeRef-sampleSizeObsAdj) > maxDiffSampSize):
		print(f"ERROR! Observed sample size without reference sample size [{sampleSizeObsAdj=}]. Available sample sizes: {', '.join([str(v) for v in sampleSizeRef_lst])}")
		sys.exit()
	return sampleSizeRef

def printIgnoredWindows(pcs_distrib_allwin, minDiffPcsSizes):
	winIgnored_lst  = [(winEndPos-winBegPos, PCSsizeDistrib) for ((winBegPos, winEndPos), PCSsizeDistrib) in pcs_distrib_allwin.items() if len(PCSsizeDistrib) < minDiffPcsSizes]
	maxPCSsize_dict = defaultdict(int)
	for (winSize, PCSsizeDistrib) in winIgnored_lst:
		maxPCSsize = max(list(PCSsizeDistrib.keys())) if (len(PCSsizeDistrib.keys()) > 0) else 0
		maxPCSsize_dict[maxPCSsize] += 1
	maxPCSsize_lst = sorted(list(maxPCSsize_dict.items()), key=lambda val: val[0])
	maxPCSsize_lst_str = [f"     Max.PCS={val[0]}\tCnt.={val[1]}\n" for val in maxPCSsize_lst]
	print(f" - Ignored Windows:\n{''.join(maxPCSsize_lst_str)}")
	nbwin_total   = len(pcs_distrib_allwin)
	nbwin_ignored = len(winIgnored_lst)
	print(f" - Ignored windows: {nbwin_ignored} out of {nbwin_total} ({(100*(nbwin_ignored/nbwin_total)):.2f} %).\n - Criteria to ignore windows: less than {minDiffPcsSizes} distinct PCS sizes above {my_dataset.minPCSsize} base pairs.")

####################################################################
# "Reading inputs" related methods.

def readWindows(my_dataset, ucscName, qChrom):
	# Check if input file exists.
	windows_filename = my_dataset.getOutFilename_computeWindows(qChrom)
	if (not os.path.isfile(windows_filename)):
		print(f"ERROR! File not found: {windows_filename}")
		sys.exit()
	# Load windows.
	pcs_distrib_all = pickle.load(open(windows_filename, 'rb'))
	if(ucscName not in pcs_distrib_all):
		print(f"ERROR! Species not found: {ucscName} ({qChrom}).")
		sys.exit()
	return pcs_distrib_all[ucscName]


def readSampleSizeRefs(winSizeRef_filepaths):
	winSizeRef_samplerefs = defaultdict(list)
	for winSizeRef, setupFilepath in winSizeRef_filepaths.items():
		# Load setup information.
		ts_sel, setupInfo = pickle.load(open(setupFilepath, 'rb'))
		winSizeRef_samplerefs[winSizeRef] = list(sorted(list(setupInfo.keys())))
	return winSizeRef_samplerefs

def readSetup(my_dataset, pcs_distrib_allwin, alpha):
	""" This function retrieves all reference window/sample sizes that 
		were pre-calculated using the script ``3_setupEvolTimes.py`` 
		and matches each window from the dataset with a corresponding 
		reference window/sample size from the collected data.

		:returns: The function returns:
			1. A dictionary mapping reference window/sample sizes to a list of windows identified by their IDs;
			2. A dictionary mapping reference window sizes to the file paths where their setup details are stored.
		
	"""

	winSizeRef_obslst = defaultdict(lambda: defaultdict(list))
	maxDiffWinSize    = my_dataset.maxDiffWinSize
	minDiffPcsSizes   = my_dataset.minDiffPcsSizes

	# Get all setup files.
	winSizeRef_filepaths = my_dataset.getOutFilenames_setupEvolTimes(alpha)

	# Read reference sample sizes for each reference window size.
	winSizeRef_samplerefs = readSampleSizeRefs(winSizeRef_filepaths)

	# Create a map between windows of the dataset 
	# and reference window size values.
	winSizeObs_lst = [(winEndPos-winBegPos, winBegPos, winEndPos) for ((winBegPos, winEndPos), PCSsizeDistrib) in pcs_distrib_allwin.items() if len(PCSsizeDistrib) >= minDiffPcsSizes]
	winSizeObs_lst = list(sorted(winSizeObs_lst, key=lambda win: win[0]))

	winSizeRef_lst = list(sorted(list(winSizeRef_filepaths.keys())))

	# Print some stats on ignored windows.
	printIgnoredWindows(pcs_distrib_allwin, minDiffPcsSizes)

	# Find a reference window size for each window.
	curIdx = 0
	for (winSizeObs, winBegPos, winEndPos) in winSizeObs_lst:
		windowId     = (winBegPos, winEndPos)
		winSizeRef   = None
		searchForRef = True
		while(searchForRef):
			winSizeRef    = winSizeRef_lst[curIdx]
			searchForRef  = (abs(winSizeObs-winSizeRef) > maxDiffWinSize)
			if(searchForRef):
				curIdx += 1
				searchProblem = ((winSizeObs < winSizeRef) or (curIdx == len(winSizeRef_lst)))
				if searchProblem: # Cases that should not happen.
					print(f"ERROR! Observed window size without reference window size [{winSizeObs=}]. {winSizeRef} {curIdx} {winSizeRef_lst}")
					sys.exit()

		if(searchForRef != None):
			# Find a reference sample size.
			sampleSizeRef = getSampleSizeRef(winSizeRef_samplerefs[winSizeRef], pcs_distrib_allwin[windowId], my_dataset.maxDiffSampSize)
			# Add to the reference tuple.
			winSizeRef_obslst[winSizeRef][sampleSizeRef].append(windowId)
		else:
			print(f"ERROR! Observed window size without reference window size [{winSizeObs=}].")
			sys.exit()
	return winSizeRef_obslst, winSizeRef_filepaths
	

####################################################################
# "Estimate evolutionary time" related methods.

def estimateEvolTimes_winRefSize(parallelInput):
	winSizeRef, my_dataset, setupFilepath, windows_info, debugMessage = parallelInput

	# Save times for each step and overall time.
	timeTrack = Time()
	timeTrack.start()

	# Load setup information.
	if(debugMessage): timeTrack.startStep(f"{debugMessage} Load setup file")
	ts_sel, setupInfo = pickle.load(open(setupFilepath, 'rb'))
	if(debugMessage): timeTrack.stopStep()

	# Calculate the posterior for each window.
	# Using the KDE objects associated with each evolutionary time, 
	# it computes the likelihood that a window is at a certain evolutionary time.
	win_posteriors = {}
	win_refs = {}
	if(debugMessage): timeTrack.startStep(f"{debugMessage} Compute posterior for {sum([len(windowId_lst) for sampleSizeRef, (windowId_lst, pcsSizeDistrib_lst) in windows_info.items()])} windows")
	for sampleSizeRef, (windowId_lst, pcsSizeDistrib_lst) in windows_info.items():
		solver = setupInfo[sampleSizeRef]
		# numpy matrix M_posterior: 
		# rows = evolutionary times / cols = window lists
		M_posterior = solver.estimate_ts(pcsSizeDistrib_lst)
		for colIdx, windowId in enumerate(windowId_lst):
			win_posteriors[windowId] = M_posterior[:, colIdx]
			win_refs[windowId] = (winSizeRef,sampleSizeRef)
	if(debugMessage): timeTrack.stopStep()

	# Return results for all sample sizes for a given window size.
	return (winSizeRef, timeTrack, ts_sel, win_refs, win_posteriors)


####################################################################
# Functions related to parallelizing computation.

def initParallelInputs(my_dataset, winSizeRef_obslst, 
					winSizeRef_filepaths, pcs_distrib_allwin, printStep):
	parallelInputs_lst = []
	winSizeRef_lst = winSizeRef_filepaths.keys()
	nbchars        = len(str(len(winSizeRef_lst)))
	for runIdx, winSizeRef in enumerate(winSizeRef_lst):
		setupFilepath  = winSizeRef_filepaths[winSizeRef]
		# Retrieve window IDs and PCS size distributions for all reference 
		# sample sizes related to the assessed reference window size.
		windows_info = {}
		for sampleSizeRef, windowId_lst in winSizeRef_obslst[winSizeRef].items():
			pcsSizeDistrib_lst = [pcs_distrib_allwin[windowId] for windowId in windowId_lst]
			windows_info[sampleSizeRef] = (windowId_lst, pcsSizeDistrib_lst)
		debugMessage = f"[{str(runIdx+1).rjust(nbchars)}/{len(winSizeRef_lst)}] " if ((runIdx % printStep) == 0) else ""
		parallelInputs_lst.append((winSizeRef, my_dataset, setupFilepath, windows_info, debugMessage))
	return parallelInputs_lst

def processIteration(result, winSizeRefs_ts, winIds_refs, winIds_ests):
	winSizeRef, it_timeTrack, winSizeRef_evolTimes, winIds_refs_new, winIds_ests_new = result

	# Update results.
	winSizeRefs_ts[winSizeRef] = winSizeRef_evolTimes
	winIds_refs = winIds_refs | winIds_refs_new
	winIds_ests = winIds_ests | winIds_ests_new
	
	return winSizeRefs_ts, winIds_refs, winIds_ests

####################################
# MAIN.
####################################
if (__name__ == '__main__'):

	parser = argparse.ArgumentParser(description="Estimate evolutionary time for each window of a given species/chromosome based on their PCS size distribution.")	
	parser.add_argument("-sp_ucsc_name", help="UCSC name of the species that is being compared to the reference species. In our study, the reference species is 'hg38' (human genome), whereas the other species is one of the 40 vertebrates.", type=str, required=True)
	parser.add_argument("-chr", help="Chromosome from the reference species for which the PCS size distributions will be used as input.", type=str, required=True)
	parser.add_argument("-alpha", help="It determines how frequent longer indels can occur. The parameter alpha can take any value above 1: (1, ∞). If alpha is near 1, larger indels are more likely to occur. If alpha is above 5 (a hard upper limit set internally with [max_alpha]), the model turns into a substitution only model.", type=float, required=True)
	parser.add_argument("-cores", help="Number of cores to be used during estimation.", type=int, required=True)
	parser.add_argument("--overwrite", action='store_true', help="If this flag is specified, any existing output files will be overwritten during the run.")

	args       = parser.parse_args()
	my_dataset = Dataset()

	prefixTarget   = args.sp_ucsc_name		 # prefix of species
	qChrom		   = args.chr		         # chromosome
	windowSize     = my_dataset.windowSize
	minPCSsize     = my_dataset.minPCSsize
	alpha          = float(args.alpha)       # how frequent longer indels can occur
	nbcores        = int(args.cores)
	overwriteFiles = args.overwrite

	print("***************************************************************")
	print("*                Estimate evolutionary times                  *")
	print("* Estimate a posterior of evolutionary times for each window. *")
	print("* All windows mapped to the specified chromosome of the       *")
	print("* reference genome are included. The PCS size distributions   *")
	print("* from these windows, representing perfectly conserved        *")
	print("* sequences between the reference genome (e.g. ``hg38``) and  *")
	print("* the species specified in the input, are used as input for   *")
	print("* the evolutionary time estimation method.                    *")
	print("***************************************************************")
	print(f"- Propensity for indels: α={alpha}")
	print(f"- Reference genome: {prefixTarget}")
	print(f"- Query genome: {my_dataset.refsp_ucscName}")
	print(f"- Chromosome: {qChrom}")
	print(f"- Window size: ~{windowSize} base pairs.")
	print(f"- Window directory (input): {my_dataset.dirWindows}")
	print(f"- Setup directory (input): {my_dataset.dirSetupEvolTimes}")
	print(f"- Evolutionary time estimates directory (output): {my_dataset.dirEstEvolTimes}")
	print(f"- Log directory with details of the run (output): {my_dataset.dirLog}")
	print(f"- Minimum size for PCS to be considered: >={minPCSsize} base pairs.")

	# Save times for each step and overall time.
	timeTrack = Time()
	timeTrack.start()

	# Read windows.
	timeTrack.startStep("Read windows")
	pcs_distrib_allwin = readWindows(my_dataset, prefixTarget, qChrom)
	timeTrack.stopStep()

	# Read setup files.
	timeTrack.startStep("Read setup")
	winSizeRef_obslst, winSizeRef_filepaths = readSetup(my_dataset, pcs_distrib_allwin, alpha)
	timeTrack.stopStep()

	# Compute evolutionary time estimates.
	# Compute evolutionary time estimates for all windows of a given species/chromosome.
	timeTrack.startStep("Compute evol. times")
	parallelInputs = initParallelInputs(my_dataset, winSizeRef_obslst, winSizeRef_filepaths, pcs_distrib_allwin, int(len(winSizeRef_filepaths)*0.1))
	print(f"Reference window sizes={len(parallelInputs)}; Cores to process window sizes={nbcores}.")
	winSizeRefs_ts, winIds_refs, winIds_ests = {}, {}, {}
	if(nbcores == 1):		
		for parallelInput in parallelInputs:
			result = estimateEvolTimes_winRefSize(parallelInput)
			winSizeRefs_ts, winIds_refs, winIds_ests = processIteration(result, winSizeRefs_ts, winIds_refs, winIds_ests)
	else:
		with futures.ProcessPoolExecutor(nbcores) as pool:
			for result in pool.map(estimateEvolTimes_winRefSize, parallelInputs):
				winSizeRefs_ts, winIds_refs, winIds_ests = processIteration(result, winSizeRefs_ts, winIds_refs, winIds_ests)
	timeTrack.stopStep()

	# Save results for chromosome.
	outputFilepath = my_dataset.getOutFilename_estimateEvolTimes(prefixTarget, qChrom, alpha)
	with open(outputFilepath, 'wb') as pickleFile:
		# It saves info about evolutionary times, window IDs, 
		# and the posteriors (likelihood values) of each window.
		pickle.dump((winIds_ests, winIds_refs, winSizeRefs_ts),pickleFile, protocol=pickle.HIGHEST_PROTOCOL)

	timeTrack.stop()
	timeTrack.print()
	print(f"Done!")
