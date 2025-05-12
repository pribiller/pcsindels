# Use: python3 3_preCompEvolTimes.py -alpha 1.1 -cores 1
#	   python3 ~/code/3_preCompEvolTimes.py -alpha 1.1 -cores 1
# It takes around 25 minutes per alpha value using 80 cores.

import os
import sys
import re
import pickle
import numpy as np
from collections import defaultdict
import argparse

from mpmath import *

from concurrent import futures

from scipy.optimize import lsq_linear
from scipy.stats import gaussian_kde
from sklearn.decomposition import PCA

from utils.basicTypes import Pcs, Time
from utils.dataset import Dataset
from utils.indelsModel import IndelsSolver

mp.dps    = 30
mp.pretty = True

####################################################################
# Sampling related methods.
def samplePCSsNonPar(input):
	LocalProcRandGen = np.random.RandomState()
	PCSsizes, PCSsizes_probs, N = input
	sample = {}
	N_samp = 0
	while(N_samp < N):
		sampleSize = min(N-N_samp,10000)
		PCSsizes_sampled = LocalProcRandGen.choice(PCSsizes, sampleSize, p=PCSsizes_probs, replace=True)
		for PCSsize in PCSsizes_sampled:
			if (PCSsize not in sample.keys()): sample[PCSsize] = 0
			sample[PCSsize] += 1
		N_samp += sampleSize

	PCSsizes_sampled = list(sorted(sample.keys()))
	PCScnts_sampled  = [sample[PCSsize] for PCSsize in PCSsizes_sampled]
	return PCSsizes_sampled, PCScnts_sampled
	
def samplePCSs(PCSsizes, PCSsizes_probs, N, nbCores=1):
	if(nbCores == 1):
		return samplePCSsNonPar((PCSsizes, PCSsizes_probs, N))
	else:
		# Prepare parallel inputs.
		N_par = int(N/(nbCores))
		N_rem = N % (nbCores)
		parallelInputs = [(PCSsizes, PCSsizes_probs, N_par) for idxSample in range(nbCores-N_rem)]
		parallelInputs.extend([(PCSsizes, PCSsizes_probs, N_par+1) for idxSample in range(N_rem)])		
		# Parallel sampling.
		#print(f"\t Parallel sampling: nbCores={nbCores} / samples per core = {N_par}")
		sample = {}
		with futures.ProcessPoolExecutor(nbCores) as pool:
			for result in pool.map(samplePCSsNonPar, parallelInputs):
				PCSsizes_sampled, PCScnts_sampled = result
				# Aggregating results.
				for (PCSsize,PCScnt) in zip(PCSsizes_sampled,PCScnts_sampled):
					if (PCSsize not in sample.keys()): sample[PCSsize] = 0
					sample[PCSsize] += PCScnt
		# Return aggregated results.
		PCSsizes_sampled = list(sorted(sample.keys()))
		PCScnts_sampled  = [sample[PCSsize] for PCSsize in PCSsizes_sampled]
		return PCSsizes_sampled, PCScnts_sampled

####################################################################
# Slope / intercept estimation.
def linearize(PCSsizes_obs,PCScnts_obs):
	y = np.log(PCScnts_obs + [1e-1])       # "real" y-axis.
	x = PCSsizes_obs+[max(PCSsizes_obs)+1] # "real" x-axis.
	# Prepare the design matrix A (with a column of ones for the intercept)
	A = np.vstack([np.ones(len(x)), x]).T  # Shape (100, 2)
	# Target vector b
	b = y  # Shape (100,)
	# Define bounds for the coefficients
	# Ideally, the slope should be negative, but in low samples this might not be the case.
	bounds = ([-np.inf, -np.inf], [np.inf, -1e-10])  # Intercept can be any value, slope <= 0
	# Solve the least squares problem
	result = lsq_linear(A, b, bounds=bounds, tol=1e-20)
	# Extract the optimized parameters
	intercept_opt, slope_opt = result.x
	return (slope_opt,intercept_opt)
	
####################################################################
# Functions involved in setting up the estimation method, i.e.: 
# (1) Selection of evolutionary times;
# (2) Computation of expected PCS size distributions;
# (3) Sampling M times N points for each the PCS size distribution;
# (4) Initialization of a KDE object based on the sample.
def selectEvolTimes(solver, maxNbTaus):
	t_span  = (0.0, 1.0)
	seltaus = solver.find_ts(t_span,maxNbTaus)
	seltaus = seltaus[1:] # Ignore t=0
	return seltaus

# WARNING! Sample bias.
# Ignores that some taus might have considerably more PCSs under minPCSsize,
# thus creating a sample with PCSs above minPCSsize would be very unlikely.
def computeEvolTimes_PcsProbs(solver, seltaus):
	seltaus_PCSprobs = []
	for t in seltaus:
		# Computes the expected counts for each PCS size given an evolutionary time t.
		PCSprobs_t = solver.comp_c_kt(t) # idx 0 = min PCS size
		PCSprobs_t[-1] = 0.0 # Ignore max PCS size (=window size).
		# Normalize counts to get probabilities.
		totalCount = sum(PCSprobs_t)
		PCSprobs_t = PCSprobs_t/totalCount
		seltaus_PCSprobs.append(PCSprobs_t)
	return seltaus_PCSprobs

def computeEvolTimes_PcsSamp(my_dataset, t_PCSsizes, t_PCSprobs, nbPCSs):
	nbSamplesPerTau = my_dataset.nbSamplesPerTau
	minDiffPcsSizes = my_dataset.minDiffPcsSizes

	slopes    =[]
	intercepts=[]
	for sampleIdx in range(nbSamplesPerTau):
		
		# Sample PCS sizes from the PCS size distribution. 
		# Use same conditions as the ones found in observed data (same count of PCSs).
		# Consider only PCS sizes above min. PCS size.
		PCSsizes_est_samp = []
		nb_attempts	      = 0
		max_nb_attempts   = 100
		# Discard PCS size distribution if number of distinct PCS sizes is below [minDiffPcsSizes].
		while (len(PCSsizes_est_samp) < minDiffPcsSizes): 
			PCSsizes_est_samp, PCScnt_est_samp = samplePCSs(t_PCSsizes,t_PCSprobs, nbPCSs, 1)
			# Makes the sampling process easier if it is taking too much time.
			nb_attempts += 1
			if(nb_attempts > max_nb_attempts): minDiffPcsSizes=minDiffPcsSizes-1

		# Exponential behavior is expected to be linear in (constrained) log-y.
		slope, intercept = linearize(PCSsizes_est_samp, PCScnt_est_samp)
		slopes.append(slope)
		intercepts.append(intercept)
	return (slopes,intercepts)

def computeEvolTimes_PcsKde(samp_info):

	# If variables are correlated, gaussian_kde gets stuck.
	# gaussian_kde does not currently support data that lies in a 
	# lower-dimensional subspace of the space in which it is expressed. 
	slopes, intercepts = samp_info

	# Combine x and y into a 2D array
	x = np.array(slopes)
	y = np.array(intercepts)
	# Combine x and y into a 2D array
	data     = np.vstack([x, y]) # KDE: shape (# of dims, # of data).
	data_pca = data.T # PCA: shape (n_samples, n_features) (flipped data).
	kde  = None
	pca  = None
	try:
		kde  = gaussian_kde(data, bw_method='scott') # 2-D array with shape (# of dims, # of data).
	# This exception was happening for the slope and nb. bps; they were not 2 distinct dimensions in few cases.
	# My guess is that the slope becomes random/uninformative when there are only few points.
	# In our case, this exception happens usually for the smallest sample size allowed (sampleSizeRef=5),
	# and rarely for other sample sizes (sampleSizeRef=8,11).
	except np.linalg.LinAlgError as e:
		print(f"[{t=} {nbSamplesPerTau=} {sampleSize=}] Gaussian KDE exception caught: {e}\n{slopes}\n{intercepts}")
		# Perform PCA to reduce dimensions
		pca = PCA(n_components=1)  # Reduce to 1 dimension
		pca.fit(data_pca)
		reduced_data = pca.transform(data_pca)
		# Now try Gaussian KDE on the reduced data.
		kde = gaussian_kde(reduced_data.T) # Gaussian KDE: shape (# of dims, # of data)
	pdf_values = kde.evaluate(data) if (pca == None) else kde.evaluate(pca.transform(data_pca).T)
	return (kde, pca, pdf_values.max())

def computeSetupInfo_sampRefSize(my_dataset, ts, t_PCSsizes, ts_PCSprobs, sampleSizeRef):
	ts_setup = []
	for t_idx, (t, t_PCSprobs) in enumerate(zip(ts, ts_PCSprobs)):
		samp_info = computeEvolTimes_PcsSamp(my_dataset, t_PCSsizes, t_PCSprobs, sampleSizeRef)
		kde_t, pca_t, max_likelihood_t = computeEvolTimes_PcsKde(samp_info)
		ts_setup.append((kde_t, pca_t, max_likelihood_t))
	return ts_setup

def computeSetupInfo_winRefSize(parallelInput):
	"""
	This function handles the setup of the estimation 
	method for a given window size.
	During the setup, the function selects evolutionary times 
	given an alpha and a window size. 
	
	A given window size also has a list of associated 
	PCS counts found in observed data. 

	Given this list, the method calculates additional 
	information pertinent to the chosen evolutionary times.
	"""
	winSizeRef, sampleSizeRef_lst, alpha, my_dataset = parallelInput

	# Save times for each step and overall time.
	timeTrack = Time()
	timeTrack.start()

	solver  = IndelsSolver(winSizeRef, alpha, my_dataset.minPCSsize)

	# Select evolutionary times.
	timeTrack.startStep("Select evol. times")
	ts_sel = selectEvolTimes(solver, my_dataset.maxNbTaus)
	timeTrack.stopStep()

	# Compute PCS size distribution for each selected evolutionary time.
	timeTrack.startStep("Compute PCS size distribution for selected evol. times")
	ts_sel_PCSprobs = computeEvolTimes_PcsProbs(solver, ts_sel)
	timeTrack.stopStep()

	# Sample PCS size distribution and create a 'KDE object' 
	# for each combination of [winRefSize][sampleRefSize].
	# KDE = Kernel Density Estimation. 
	# KDE computes the probability that a point was sampled 
	# from the same distribution as a group of points.
	t_PCSsizes = solver.sizes
	setupInfo  = {}
	timeTrack.startStep("Sample PCSs")
	for sampleIdx, sampleSizeRef in enumerate(sampleSizeRef_lst):
		print(f"[{alpha} - {winSizeRef}] [{sampleIdx+1} out of {len(sampleSizeRef_lst)}] Starting sampling... ({sampleSizeRef=})")
		# Setup for each reference sample size.
		setupInfo[sampleSizeRef] = computeSetupInfo_sampRefSize(my_dataset, ts_sel, solver.sizes, ts_sel_PCSprobs, sampleSizeRef)
	timeTrack.stopStep()

	# Save results for all sample sizes for a given window size.
	outFilename = my_dataset.getOutFilename_preCompEvolTimes(winSizeRef,alpha)
	with open(outFilename, 'wb') as pickleFile:
		pickle.dump((ts_sel, setupInfo),pickleFile, protocol=pickle.HIGHEST_PROTOCOL)

	timeTrack.stop()
	return (winSizeRef, timeTrack)

########################################
# Functions related to reference values 
# for window size and sample size, based on real data.

def loadObservedValues(dirWindows, windowSize, prefixQuery, qChromLst, speciesUCSCnames):
	""" This function loads all computed windows and returns the observed
	window sizes and sample sizes (i.e. the total count of PCSs in a window).
	"""
	obsVals = defaultdict(int)
	# Get all observed window sizes and sample sizes (= total count of PCSs).
	for qChrom in qChromLst:
		print(f"\t[{qChrom}] Loading observed values.")
		# Check if input file exists.
		pcs_distrib_filename = os.path.join(dirWindows, f"{prefixQuery}.{qChrom}.{windowSize}.windows.pickle")
		if (not os.path.isfile(pcs_distrib_filename)):
			print(f"ERROR! File not found: {pcs_distrib_filename}")
			sys.exit()
		# Load windows.
		pcs_distrib_all = pickle.load(open(pcs_distrib_filename, 'rb'))
		for ucscName_other in speciesUCSCnames:
			if(ucscName_other not in pcs_distrib_all):
				print(f"ERROR! Species not found: {ucscName_other} ({qChrom}).")
				sys.exit()
			# Get the PCS size distribution of the species/chromosome.
			pcs_distrib_allwin = pcs_distrib_all[ucscName_other]
			for (windowId, PCSsizeDistrib) in pcs_distrib_allwin.items():
				winBegPos, winEndPos = windowId
				winSizeObs  = winEndPos-winBegPos
				sampSizeObs = sum(list(PCSsizeDistrib.values()))
				if(sampSizeObs > 0): obsVals[(winSizeObs,sampSizeObs)] += 1 
	return obsVals

def getReferenceValues(obsVals, maxDiffWinSize, maxDiffSampSize):
	""" This function computes the reference values for window size and 
	sample size, based on the real data from the dataset.

	For each window in the chromosome of the reference species, there 
	exists a PCS size distribution linked to that window, which varies 
	based on the other species included in the analysis.

	Every window has two associated values: ``window size`` and ``sample size``.
	``Window size`` is defined as the difference between the start and 
	end coordinates of the window, while ``sample size`` is the total 
	number of PCSs observed in that window.

	Due to the high computational cost of estimating evolutionary time 
	for all possible combinations of window sizes and sample sizes in 
	the dataset, we use this clustering function.

	The clustering function is implemented as follows. 

	Each window is assigned to a cluster, which is defined by 
	a ``reference window size`` and a ``reference sample size``.
	Clusters are defined such that the window sizes within a cluster 
	cannot differ from the reference window size by more than ``maxDiffWinSize``. 
	There is also an upper limit on the difference between the 
	reference and observed sample sizes, set by ``maxDiffSampSize``.
	"""
	refVals     = {}
	refVals_lst = []
	obsVals_lst = sorted(list(obsVals.keys()), key=lambda val: (val[0], val[1]))
	for (winSizeObs,sampSizeObs) in obsVals_lst:
		# Find reference value for each observed value.
		refKey  = None
		for (w,s) in reversed(refVals_lst):
			if ((winSizeObs-w) <= maxDiffWinSize):
				if ((sampSizeObs-s) <= maxDiffSampSize):
					refKey = (winSizeObs,sampSizeObs)
					break
			else:
				break
		# Add a new reference value in case the existing ones
		# cannot be linked to the current observed value.
		if(refKey == None):
			refKey = (winSizeObs,sampSizeObs)
			refVals_lst.append(refKey)
		refVals[(winSizeObs,sampSizeObs)] = refKey
	return refVals

def checkPCSsizeDistribFiles(dirWindows, windowSize, prefixQuery, qChromLst):
	""" This function checks whether the PCS size distributions 
	have been computed for every window across all chromosomes in the dataset.
	"""
	for qChrom in qChromLst:
		# Check if input files exist.
		pcs_distrib_filename = os.path.join(dirWindows, f"{prefixQuery}.{qChrom}.{windowSize}.windows.pickle")
		if (not os.path.isfile(pcs_distrib_filename)):
			print(f"ERROR! PCS size distribution not found for chromosome '{qChrom}' (filepath: {pcs_distrib_filename}).")
			sys.exit()

####################################################################
# Functions related to parallelizing computation.
def initParallelInputs(refVals, alpha, my_dataset):
	refVals_lst = sorted(list(set(refVals.values())), key=lambda val: (val[0], val[1]))
	print(f"- Total reference values: {len(refVals_lst)}")
	# Group values by reference window size.
	parallelInputs = defaultdict(list)
	for (winSizeRef,sampSizeRef) in refVals_lst:
		parallelInputs[winSizeRef].append(sampSizeRef)
	# Create list of parallel inputs.
	parallelInputs_lst = []
	winSizeRef_lst = sorted(list(parallelInputs))
	for winSizeRef in winSizeRef_lst:
		parallelInputs_lst.append((winSizeRef, parallelInputs[winSizeRef], alpha, my_dataset))
	return parallelInputs_lst

def processIteration(cntRuns, stepRun, totRuns, result, it_timeTrack_max):
	winSizeRef, it_timeTrack = result
	if((cntRuns % stepRun) == 0): 
		print(f"{cntRuns} out of {totRuns}.\n\tHighest running time in an iteration so far:")
		it_timeTrack_max.print()
	if((it_timeTrack_max == None) or (it_timeTrack_max.t_global < it_timeTrack.t_global)):
		it_timeTrack_max = it_timeTrack
	return cntRuns+1, it_timeTrack_max

####################################
# MAIN.
####################################
if (__name__ == '__main__'):

	parser = argparse.ArgumentParser(description="Extract PCSs (Perfectly Conserved Sequences) from Chains (ordered aligned blocks).")
	parser.add_argument("-alpha", help="It determines how frequent longer indels can occur. The parameter alpha can take any value above 1: (1, ∞). If alpha is near 1, larger indels are more likely to occur. If alpha is above 5 (a hard upper limit set internally with [max_alpha]), the model turns into a substitution only model.", type=float, required=True)
	parser.add_argument("-cores", help="Number of cores to be used during setup.", type=int, required=True)
							
	args       = parser.parse_args()
	my_dataset = Dataset()

	alpha   = float(args.alpha)
	nbcores = int(args.cores)

	windowSize = my_dataset.windowSize   # Window size (in number of base pairs).
	minPCSsize = my_dataset.minPCSsize	 # Minimum size (in number of base pairs) for the PCS to be considered.

	prefixQuery	     = my_dataset.refsp_ucscName     # In our study, query is always the human genome (hg38)
	qChromLst        = my_dataset.chromLst           # All chromosomes from query genome (chr1, chr2, ..., chrY)
	qChromLst        = ["chr16"]
	speciesUCSCnames = my_dataset.speciesUCSCnames

	maxDiffWinSize   = my_dataset.maxDiffWinSize
	maxDiffSampSize  = my_dataset.maxDiffSampSize
	nbSamplesPerTau  = my_dataset.nbSamplesPerTau

	dirOut		     = my_dataset.dirSetupEvolTimes  # Directory where pre-computed info will be saved (output directory).
	dirWin           = my_dataset.dirWindows

	print("***************************************************")
	print("*      Evolutionary Time Inference: Setup         *")
	print("* Pre-compute the Kernel Density Estimation (KDE) *")
	print("* for various window sizes and PCS counts (sample *")
	print("* sizes) derived from the dataset.                *")
	print("***************************************************")
	print(f"Core parameter values:")
	print(f"----------------------")
	print(f" - Query genome: {prefixQuery}")
	print(f" - Window size: {windowSize}")
	print(f" - Minimum PCS size: {minPCSsize}")
	print(f" - Number of cores: {nbcores}")
	print(f"Setup-specific parameter values:")
	print(f"--------------------------------")
	print(f" - Input directory with windows and their PCS size distribution: {dirWin}")
	print(f" - Propensity for indels: α={alpha}")
	print(f" - Number of sampled PCS size distributions per evolutionary time: {nbSamplesPerTau}")
	print(f" - Constraint on reference window sizes: max(|[obs. window size]-[ref. window size]|) ≤ {maxDiffWinSize}")
	print(f" - Constraint on reference sample sizes: max(|[obs. sample size]-[ref. sample size]|) ≤ {maxDiffSampSize}")
	print(f" - Output directory for pre-computed data: {dirOut}\n\n")

	timeTrack = Time()
	timeTrack.start()

	#######################################
	# Check pre-requisites.

	# Check whether the PCS size distributions have been computed for every 
	# window across all chromosomes and genomes in the dataset.
	checkPCSsizeDistribFiles(dirWin, windowSize, prefixQuery, qChromLst)

	# Load observed window sizes and sample sizes (total count of PCSs).
	timeTrack.startStep("Load observed window sizes and PCS counts.")
	obsVals = loadObservedValues(dirWin, windowSize, prefixQuery, qChromLst, speciesUCSCnames)
	timeTrack.stopStep()
	# Compute reference values for window and sample sizes.
	timeTrack.startStep("Compute reference size values.")
	refVals = getReferenceValues(obsVals, maxDiffWinSize, maxDiffSampSize)
	timeTrack.stopStep()

	#######################################
	# Compute evolutionary times for reference values.
	timeTrack.startStep("Setup method.")
	parallelInputs  = initParallelInputs(refVals, alpha, my_dataset)
	cntRuns, stepRun, totRuns = 0, 1, len(parallelInputs)
	w_timeTrack_max = None
	if(nbcores == 1):		
		for parallelInput in parallelInputs:
			result = computeSetupInfo_winRefSize(parallelInput)
			cntRuns, it_timeTrack_max = processIteration(cntRuns, stepRun, totRuns, result, it_timeTrack_max)
	else:
		with futures.ProcessPoolExecutor(nbcores) as pool:
			for (w, w_timeTrack) in pool.map(computeSetupInfo_winRefSize, parallelInputs):
				cntRuns, it_timeTrack_max = processIteration(cntRuns, stepRun, totRuns, result, it_timeTrack_max)
	timeTrack.stopStep()

	timeTrack.stop()
	timeTrack.print()
	print("Done!")
