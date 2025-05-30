"""Reads all PCS size distributions across the windows of a designated chromosome in 
a given species (stored in one of the outputs files of ``2_computeWindows.py``)
and, along with the setup files produced by the script ``3_setupEvolTimes.py``, 
estimates the evolutionary times for each window.

- **Use**::

	python3 5_1_fig_PCSsizeDistribComp.py -sp_ucsc_name [UCSC name] -alpha [any number > 1] -cores [nb. of cores] [--overwrite, optional]

- **Example of Usage (human (reference genome) and mouse)**::

	python3 ~/code/5_1_fig_PCSsizeDistribComp.py -sp_ucsc_name mm39 -alpha 1.1 -cores 1 --overwrite
	python3 ~/code/5_1_fig_PCSsizeDistribComp.py -sp_ucsc_name mm39 -alpha 1.1 -cores 70 --overwrite

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

max to ask: 62,68,62,60,57,52 cores.
It takes 10 minutes for whole genome, 10 sample, 70 cores.

"""

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

	stepPrint = max(int(len(windows)*0.01),1)
	for stepIdx, (windowId, winSizeRef, winPcsObs, winEvolTimes) in enumerate(windows):

		if(stepIdx % stepPrint == 0):
			print(f"\t{stepIdx} out of {len(windows)} windows.")

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

def processIteration(chrom, result, PCSdistrib_obs_whl, PCSdistrib_est_whl, taudistrib_est_whl, 
	PCSdistrib_obs_chr, PCSdistrib_est_chr, taudistrib_est_chr):
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
			taudistrib_est_chr[sampleIdx][t] += cnt
		# Update sampled PCS size distribution.
		for (PCSsize, PCScnt) in taudistrib_samp.items():
			PCSdistrib_est_whl[sampleIdx][PCSsize] += PCScnt
			PCSdistrib_est_chr[chrom][sampleIdx][PCSsize] += PCScnt
	return PCSdistrib_obs_whl, PCSdistrib_est_whl, taudistrib_est_whl, PCSdistrib_obs_chr, PCSdistrib_est_chr, taudistrib_est_chr

####################################
# MAIN.
####################################
if (__name__ == '__main__'):

	parser = argparse.ArgumentParser(description="Process data and plot figure with PCS size distribution comparison.")
	parser.add_argument("-sp_ucsc_name", help="UCSC name of the species that is being compared to the reference species. In our study, the reference species is 'hg38' (human genome), whereas the other species is one of the 40 vertebrates.", type=str, required=True)
	parser.add_argument("-alpha", help="It determines how frequent longer indels can occur. The parameter alpha can take any value above 1: (1, ∞). If alpha is near 1, larger indels are more likely to occur. If alpha is above 5 (a hard upper limit set internally with [max_alpha]), the model turns into a substitution only model.", type=float, required=True)
	parser.add_argument("-cores", help="Number of cores to be used during estimation.", type=int, required=True)
	parser.add_argument("--overwrite", action='store_true', help="If this flag is specified, any existing output files will be overwritten during the run.")

	args       = parser.parse_args()
	my_dataset = Dataset()

	prefixTarget    = args.sp_ucsc_name		 # prefix of species
	windowSize      = my_dataset.windowSize
	minPCSsize      = my_dataset.minPCSsize
	alpha           = float(args.alpha)       # how frequent longer indels can occur
	nbSamplesPerWin = my_dataset.nbSamplesPerWin
	nbcores         = int(args.cores)
	overwriteFiles  = args.overwrite

	print("************************************************************")
	print("*          Plot PCS size distribution comparison           *")
	print("* Process the posterior evolutionary times for each window *")
	print("* from a given pairwise alignment. Plot the comparisons of *")
	print("* estimated and observed PCS size distributions across the *")
	print("* whole genome and each chromosome, along with an inset    *")
	print("* displaying the sampled evolutionary time distribution in *")
	print("* each plot.                                               *")
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
	PCSdistrib_obs_whl = defaultdict(int)
	PCSdistrib_est_whl = [defaultdict(float) for sampleIdx in range(nbSamplesPerWin)]
	taudistrib_est_whl = [defaultdict(float) for sampleIdx in range(nbSamplesPerWin)]
	PCSdistrib_obs_chr = {}
	PCSdistrib_est_chr = {}
	taudistrib_est_chr = {}
	
	# Read estimated evolutionary times for each chromosome.
	chromLst = my_dataset.chromLst
	for qChrom in chromLst:

		# Read windows and evolutionary times.
		timeTrack.startStep(f"[{qChrom}] Load data (obs. PCSs + estimates)")
		parallelInputs = initParallelInputs(my_dataset, prefixTarget, qChrom, alpha, nbcores)
		timeTrack.stopStep()

		# Process windows in a parallelized way.
		timeTrack.startStep(f"[{qChrom}] Sample results")
		if(nbcores == 1):		
			for parallelInput in parallelInputs:
				result = sampleEvolTimes(parallelInput)
				PCSdistrib_obs_whl, PCSdistrib_est_whl, taudistrib_est_whl, PCSdistrib_obs_chr, PCSdistrib_est_chr, taudistrib_est_chr = processIteration(result, PCSdistrib_obs_whl, PCSdistrib_est_whl, taudistrib_est_whl, PCSdistrib_obs_chr, PCSdistrib_est_chr, taudistrib_est_chr)
		else:
			with futures.ProcessPoolExecutor(nbcores) as pool:
				for result in pool.map(sampleEvolTimes, parallelInputs):
					winSizeRefs_ts, winIds_refs, winIds_ests = processIteration(result, winSizeRefs_ts, winIds_refs, winIds_ests)
		timeTrack.stopStep()
	
	# Save results for pairwise alignment.
	outputFilepath = my_dataset.getOutFilename_fig_PCSsizeDistribComp(prefixTarget, alpha)
	with open(outputFilepath, 'wb') as pickleFile:
		# It saves info about evolutionary times, window IDs, 
		# and the posteriors (likelihood values) of each window.
		pickle.dump((PCSdistrib_obs_whl, PCSdistrib_est_whl, taudistrib_est_whl, 
			PCSdistrib_obs_chr, PCSdistrib_est_chr, taudistrib_est_chr), pickleFile, protocol=pickle.HIGHEST_PROTOCOL)

	# Plot information.
	# fig=plt.figure()

	timeTrack.stop()
	timeTrack.print()
	print(f"Done!")
