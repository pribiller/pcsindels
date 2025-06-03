"""Process the posterior evolutionary times for each window 
of a given pairwise alignment, then sample evolutionary 
times from each window. The sampled evolutionary times 
are then used to derive distributions for evolutionary 
time and PCS size at the whole genome level as well 
as for each chromosome.

- **Use**::

	python3 5_sampleEvolTimes.py.py -sp_ucsc_name [UCSC name] -alpha [any number > 1] -cores [nb. of cores] [--overwrite, optional]

- **Example of Usage (human (reference genome) and mouse)**::

	python3 ~/code/5_sampleEvolTimes.py -sp_ucsc_name mm39 -alpha 1.1 -cores 1 --overwrite
	python3 ~/code/5_sampleEvolTimes.py -sp_ucsc_name mm39 -alpha 1.1 -cores 80 --overwrite

- **Input Parameter (mandatory)**:

:-sp_ucsc_name: UCSC name of the species that is being aligned with the reference species (e.g. *mm39* for mouse).
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
:dirEstEvolTimes: Directory where estimates of evolutionary time will be saved (input files).
:dirSampEvolTimes: Directory where (1) sampled PCS size distribution and (2) sampled evolutionary time distribution will be saved (output files).

.. note::
	Make sure that the required parameters described above are correctly defined in the file ``utils/dataset.py``.

Pre-requisites
--------------

Before using this script, make sure all the required files were pre-computed:

a) Files with estimated evolutionary times per window
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Make sure to run ``4_estimateEvolTimes.py`` for the alpha value specified in the 
input parameter ``-alpha``.

Cluster resources
-----------------

In case you want to run this Python script stand-alone in a cluster that uses Slurm to manage jobs::

	srun -p compute -t 2:00:00 --mem 20G --nodes=1 --ntasks=1 --cpus-per-task=80 --pty bash

Otherwise you can use the script ``../cluster/5_sampleEvolTimes_runAll.py``
to run this script for all 40 species used in this study (whole genome).

Time, Memory & Disk space
-------------------------

For reference, here we include a run example, with runtime, memory usage, and disk 
space required for running this script on each pairwise alignment of the 40 
vertebrate dataset examined in our study. **80 cores** were used in these runs.

===========  =========  =========  ========
Desc.        Time       Memory     Disk    
===========  =========  =========  ========
panPan3      17:18:51   28GB       0.004GB 
panTro6      18:04:45   23GB       0.004GB 
gorGor6      17:14:12   17GB       0.002GB 
ponAbe3      16:27:18   12GB       0.001GB 
papAnu4      15:34:33   13GB       0.002GB 
macFas5      16:00:06   12GB       0.001GB 
rhiRox1      16:16:28   12GB       0.002GB 
chlSab2      16:48:04   12GB       0.001GB 
nasLar1      13:54:34   12GB       0.001GB 
rheMac10     15:39:37   12GB       0.001GB 
calJac4      14:29:21   11GB       0.002GB 
tarSyr2      13:48:03   12GB       0.002GB 
micMur2      12:01:39   13GB       0.002GB 
galVar1      14:11:28   12GB       0.002GB 
mm39         09:48:15   8GB        0.002GB 
oryCun2      11:24:42   10GB       0.002GB 
rn7          09:23:12   8GB        0.002GB 
vicPac2      13:35:36   11GB       0.002GB 
bisBis1      12:01:16   10GB       0.002GB 
felCat9      13:22:03   11GB       0.002GB 
manPen1      12:05:07   10GB       0.002GB 
bosTau9      10:04:48   10GB       0.002GB 
canFam6      12:30:33   11GB       0.002GB 
musFur1      13:43:46   14GB       0.002GB 
neoSch1      14:44:31   11GB       0.002GB 
equCab3      14:12:16   12GB       0.002GB 
myoLuc2      09:59:07   10GB       0.002GB 
susScr11     11:51:39   12GB       0.002GB 
enhLutNer1   13:57:12   11GB       0.002GB 
triMan1      13:24:52   11GB       0.002GB 
macEug2      03:44:52   8GB        0.002GB 
ornAna2      02:52:52   8GB        0.002GB 
aptMan1      02:14:18   7GB        0.002GB 
galGal6      01:55:25   7GB        0.002GB 
thaSir1      01:35:42   8GB        0.002GB 
aquChr2      02:08:16   8GB        0.002GB 
melGal5      01:59:52   7GB        0.002GB 
xenLae2      01:29:09   8GB        0.002GB 
xenTro10     01:28:15   6GB        0.002GB 
danRer11     01:08:39   7GB        0.002GB 

Time per Run: Details
^^^^^^^^^^^^^^^^^^^^^

Stats on time of a single run (human-mouse alignment, chromosome 16, 80 cores): **~12 minutes**
Details on computational time are available in the log of the run.

==========================================  =======
Step                                        Time (s)
==========================================  =======
[chr16] Load data (obs. PCSs + estimates)     16.80
[chr16] Sample results                       689.99
**Total time**                               706.80
==========================================  =======

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
import random

from utils.indelsModel import IndelsModel, EvolTimesSolver, getSampleSize
from utils.basicTypes import Time
from utils.dataset import Dataset
from utils.outputs import splitList, readWindows, readEvolutionaryTimes

mp.dps = 30
mp.pretty = True

####################################################################
# Functions related to sampling of evolutionary times.

def sample_pcs(model, solver, ts_sampled, winPcsObs, nbSamplesPerWin):
	# Compute sampled PCS size distribution.
	nbPCSs_obs = getSampleSize(winPcsObs)
	maxPCS_obs = max(list(winPcsObs.keys()))
	winSamples = []
	for t in ts_sampled:
		PCSdistrib_samp = solver.computeSample_t(t, model, 1, nbPCSs_obs, 1)[0]
		maxPCS_est = max(list(PCSdistrib_samp.keys()))
		winSamples.append((t, PCSdistrib_samp, maxPCS_est))
	winSamples = sorted(winSamples, key=lambda samp: abs(samp[2]-maxPCS_obs))
	winSamples = winSamples[:nbSamplesPerWin]
	random.shuffle(winSamples)
	return winSamples

def sampleEvolTimes(parallelInput):
	my_dataset, alpha, winSizeRefs_ts, windows = parallelInput

	nbSamplesPerWin = my_dataset.nbSamplesPerWin
	nbSamplesTotal  = nbSamplesPerWin*3

	winSizeRef_prev = -1
	model  = None # Model depends on the reference window size.
	solver = EvolTimesSolver(None)

	PCSdistrib_obs_win = defaultdict(int)
	PCSdistrib_est_win = [defaultdict(float) for sampleIdx in range(nbSamplesPerWin)]
	taudistrib_est_win = [defaultdict(float) for sampleIdx in range(nbSamplesPerWin)]

	for stepIdx, (windowId, winSizeRef, winPcsObs, winEvolTimes) in enumerate(windows):
		# Initialize indel model with reference window size.
		if(winSizeRef != winSizeRef_prev):
			model = IndelsModel(winSizeRef, alpha, my_dataset.minPCSsize)
			winSizeRef_prev = winSizeRef

		# Sample evolutionary times from window posterior.
		ts = winSizeRefs_ts[winSizeRef]
		ts_sampled = None
		try:
			ts_sampled = np.random.choice(ts, size=nbSamplesTotal, p=np.array(winEvolTimes)/np.sum(winEvolTimes))
		except Exception as e:
			print(f"ERROR! {e}\n\t- Nan found: {windowId=} {winPcsObs=} {winEvolTimes=} {winSizeRef=}")
			sys.exit()

		# Compute sampled PCS size distribution.
		winSamples = sample_pcs(model, solver, ts_sampled, winPcsObs, nbSamplesPerWin)

		# Update results.
		for (PCSsize, PCScnt) in winPcsObs.items(): 
			PCSdistrib_obs_win[PCSsize] += PCScnt
		for sampleIdx, (t, PCSdistrib_samp, maxPCS_est) in enumerate(winSamples):
			taudistrib_est_win[sampleIdx][t] += 1
			for (PCSsize, PCScnt) in PCSdistrib_samp.items():
				PCSdistrib_est_win[sampleIdx][PCSsize] += PCScnt
	return (PCSdistrib_obs_win, PCSdistrib_est_win, taudistrib_est_win)

####################################################################
# Functions related to parallelizing computation.

def initParallelInputs(my_dataset, prefixTarget, qChrom, alpha, nbcores):

	timeTrack = Time()
	timeTrack.start()

	# Read windows.
	timeTrack.startStep("Load windows")
	pcs_distrib_allwin = readWindows(my_dataset, prefixTarget, qChrom)
	timeTrack.stopStep()

	# Read evolutionary times.
	timeTrack.startStep("Load evol. times")
	winIds_ests, winIds_refs, winSizeRefs_ts = readEvolutionaryTimes(my_dataset, prefixTarget, qChrom, alpha)
	timeTrack.stopStep()

	timeTrack.startStep("Create inputs")
	# For each window.
	parallelInputs=[]
	for (windowId, winEvolTimes) in winIds_ests.items():
		winSizeRef, sampleSizeRef = winIds_refs[windowId]
		winPcsObs = pcs_distrib_allwin[windowId]
		parallelInputs.append((windowId, winSizeRef, winPcsObs, winEvolTimes))

	# Sort windows by reference window size, to speed up sampling step later.
	parallelInputs = sorted(parallelInputs, key=lambda win: win[1])
	timeTrack.stopStep()

	# Split windows uniformly across available cores.
	print(f"Windows w/estimates={len(parallelInputs)}; Cores to process windows={nbcores}; {int(len(parallelInputs)/nbcores)} windows per core.")
	parallelInputs = splitList(parallelInputs, nbcores)
	parallelInputs = [(my_dataset, alpha, winSizeRefs_ts, parallelInput) for parallelInput in parallelInputs]

	timeTrack.stop()
	timeTrack.print()
	return parallelInputs

def processIteration(chrom, result, info_saved):
	PCSdistrib_obs_whl, PCSdistrib_est_whl, taudistrib_est_whl, PCSdistrib_obs_chr, PCSdistrib_est_chr, taudistrib_est_chr = info_saved
	PCSdistrib_obs_win, PCSdistrib_est_win, taudistrib_est_win = result
	# Update observed PCS size distribution.
	for (PCSsize, PCScnt) in PCSdistrib_obs_win.items(): 
		PCSdistrib_obs_whl[PCSsize] += PCScnt
		PCSdistrib_obs_chr[chrom][PCSsize] += PCScnt
	# Update estimates.
	for sampleIdx, (taudistrib_samp, PCSdistrib_samp) in enumerate(zip(taudistrib_est_win, PCSdistrib_est_win)):
		# Update sampled evolutionary time distribution.
		for (t, cnt) in taudistrib_samp.items():
			taudistrib_est_whl[sampleIdx][t] += cnt
			taudistrib_est_chr[chrom][sampleIdx][t] += cnt
		# Update sampled PCS size distribution.
		for (PCSsize, PCScnt) in PCSdistrib_samp.items():
			PCSdistrib_est_whl[sampleIdx][PCSsize] += PCScnt
			PCSdistrib_est_chr[chrom][sampleIdx][PCSsize] += PCScnt
	return info_saved

def initializeOutputs(chromLst, nbSamplesPerWin):
	PCSdistrib_obs_whl = defaultdict(int)
	PCSdistrib_est_whl = [defaultdict(float) for sampleIdx in range(nbSamplesPerWin)]
	taudistrib_est_whl = [defaultdict(float) for sampleIdx in range(nbSamplesPerWin)]
	PCSdistrib_obs_chr = {qChrom: defaultdict(int) for qChrom in chromLst}
	PCSdistrib_est_chr = {qChrom: [defaultdict(float) for sampleIdx in range(nbSamplesPerWin)] for qChrom in chromLst}
	taudistrib_est_chr = {qChrom: [defaultdict(float) for sampleIdx in range(nbSamplesPerWin)] for qChrom in chromLst}	
	return (PCSdistrib_obs_whl, PCSdistrib_est_whl, taudistrib_est_whl, PCSdistrib_obs_chr, PCSdistrib_est_chr, taudistrib_est_chr)

####################################
# MAIN.
####################################
if (__name__ == '__main__'):

	parser = argparse.ArgumentParser(description="Sample evolutionary times from each window and derive a distribution of evolutionary times and PCS sizes at whole-genome level.")
	parser.add_argument("-sp_ucsc_name", help="UCSC name of the species that is being compared to the reference species. In our study, the reference species is 'hg38' (human genome), whereas the other species is one of the 40 vertebrates.", type=str, required=True)
	parser.add_argument("-alpha", help="It determines how frequent longer indels can occur. The parameter alpha can take any value above 1: (1, ∞). If alpha is near 1, larger indels are more likely to occur. If alpha is above 5 (a hard upper limit set internally with [max_alpha]), the model turns into a substitution only model.", type=float, required=True)
	parser.add_argument("-cores", help="Number of cores to be used during estimation.", type=int, required=True)
	parser.add_argument("--overwrite", action='store_true', help="If this flag is specified, any existing output files will be overwritten during the run.")

	args       = parser.parse_args()
	my_dataset = Dataset()

	prefixTarget    = args.sp_ucsc_name		      # prefix of species
	chromLst        = my_dataset.chromLst
	windowSize      = my_dataset.windowSize
	minPCSsize      = my_dataset.minPCSsize
	alpha           = float(args.alpha)           # how frequent longer indels can occur
	nbSamplesPerWin = my_dataset.nbSamplesPerWin  # how many evolutionary times will be sampled from each window.
	nbcores         = int(args.cores)
	overwriteFiles  = args.overwrite

	print("************************************************************")
	print("*                Sample evolutionary times                 *")
	print("* Process the posterior evolutionary times for each window *")
	print("* of a given pairwise alignment, then sample evolutionary  *")
	print("* times from each window. The sampled evolutionary times   *")
	print("* are then used to derive distributions for evolutionary   *")
	print("* time and PCS size at the whole genome level as well as   *")
	print("* for each chromosome.                                     *")
	print("************************************************************")
	print(f"- Propensity for indels: α={alpha}")
	print(f"- Reference genome: {prefixTarget}")
	print(f"- Query genome: {my_dataset.refsp_ucscName}")
	print(f"- Window size: ~{windowSize} base pairs.")
	print(f"- Window directory (input): {my_dataset.dirWindows}")
	print(f"- Evolutionary time estimates directory (input): {my_dataset.dirEstEvolTimes}")
	print(f"- Evolutionary time samples directory (output): {my_dataset.dirSampEvolTimes}")
	print(f"- Log directory with details of the run (output): {my_dataset.dirLog}")
	print(f"- Minimum size for PCS to be considered: >={minPCSsize} base pairs.")

	# Save times for each step and overall time.
	timeTrack = Time()
	timeTrack.start()

	# Information saved:
	# - observed and estimated PCS size distribution per chromosome as well as whole-genome.
	# - sampled evolutionary time distribution per chromosome as well as whole-genome.
	info_saved = initializeOutputs(chromLst, nbSamplesPerWin)

	outputFilepath = my_dataset.getOutFilename_sampleEvolTimes(prefixTarget, alpha)
	if(os.path.isfile(outputFilepath) and (not overwriteFiles)):
		print(f"WARNING! Output file already exists: {outputFilepath}. Skipping computation [{overwriteFiles=}].")
		sys.exit()

	# Read estimated evolutionary times for each chromosome.
	for qChrom in chromLst:

		# Read windows and evolutionary times.
		timeTrack.startStep(f"[{qChrom}] Load data (obs. PCSs + estimates)")
		parallelInputs = initParallelInputs(my_dataset, prefixTarget, qChrom, alpha, nbcores)
		timeTrack.stopStep()

		# Process windows in a parallelized way.
		timeTrack.startStep(f"[{qChrom}] Sample results")
		if(nbcores == 1):		
			for parallelInput in parallelInputs:
				result     = sampleEvolTimes(parallelInput)
				info_saved = processIteration(qChrom, result, info_saved)
		else:
			with futures.ProcessPoolExecutor(nbcores) as pool:
				for result in pool.map(sampleEvolTimes, parallelInputs):
					print(f"\tProcessing results")
					info_saved = processIteration(qChrom, result, info_saved)
		timeTrack.stopStep()
	
	# Save results for pairwise alignment.
	with open(outputFilepath, 'wb') as pickleFile:
		# It saves sampled evolutionary times and PCSs for each chromosome and also whole-genome.
		pickle.dump(info_saved, pickleFile, protocol=pickle.HIGHEST_PROTOCOL)

	timeTrack.stop()
	timeTrack.print()
	print(f"Done!")
