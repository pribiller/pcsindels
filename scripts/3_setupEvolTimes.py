"""Setup procedure for computing evolutionary times. The setup process includes the following steps:

1. 	Read all windows in the dataset (computed via ``2_computeWindows.py``) to obtain 
	their **window sizes** and the number of PCSs in each window, here referred 
	to as **sample sizes**. Based on these observed values, the initial step of 
	the setup method will compute **reference window sizes** and **reference sample sizes** 
	for the dataset, which are detailed further in the documentation for the 
	method ``getReferenceValues``. The result of this process will be a list of tuples 
	in the format of (**ref. window size**, **ref. sample size**);
2. 	For each distinct *reference window size* calculated in the previous step, 
	a corresponding set of evolutionary times is chosen. This selection is based 
	on the *distinctiveness* of their expected PCS size distributions, 
	as detailed in the method ``find_ts`` implemented in ``utils.indelsModel``;
3. 	After associating each reference window size with a set of evolutionary times, 
	the reference sample sizes are considered. For each reference sample size, 
	the expected PCS size distribution for a given evolutionary time is sampled 
	``N`` times to initialize a KDE object. A KDE object uses *Kernel Density Estimation* 
	to assess the likelihood of a point to be sampled from the same distribution. 
	The sample size is taken into account in our method because for smaller samples 
	the observed PCS size distribution can significantly deviate from the 
	expected size distribution, and this can be used to derive an uncertainty 
	measure for the estimate.

- **Use**::

	python3 3_setupEvolTimes.py -alpha [any number > 1] -cores [nb. of cores] --overwrite

- **Example of Usage (high probability of indels; no parallelization; no overwrite)**::

	python3 ~/code/3_setupEvolTimes.py -alpha 1.1 -cores 1

- **Input Parameters**:

:-alpha: 	It determines how frequent longer indels can occur. 
			The parameter alpha can take any value above 1 (1 is **not** included): (1, ∞).
			If alpha is near 1, larger indels are more likely to occur.
			If alpha is above 5 (a hard upper limit set internally with ``max_alpha``), 
			the model turns into a substitution only model.

:-cores: If ``cores=1``, the script will execute serially, which may take several hours to complete. It is advisable to utilize as many available cores as possible.
:--overwrite: Optional flag. If this flag is specified, any existing output files will be overwritten during the run.

- **Other Parameters taken from** ``dataset.py``:

:maxDiffWinSize: Constraint on reference window sizes: ``max(|obs. window size-ref. window size|) ≤ maxDiffWinSize``
:maxDiffSampSize: Constraint on reference sample sizes: ``max(|obs. sample size-ref. sample size|) ≤ maxDiffSampSize``
:minSampRefSize: Constraint on reference sample sizes: ``minSampRefSize ≤ ref. sample size ≤ maxSampRefSize``
:maxSampRefSize: Constraint on reference sample sizes: ``minSampRefSize ≤ ref. sample size ≤ maxSampRefSize``
:nbSamplesPerTau: Number of sampled PCS size distributions per evolutionary time.
:dirWindows: Directory where windows computed with the script ``2_computeWindows.py`` were saved (input files).
:dirSetupEvolTimes: Directory where setup files will be saved (output files).

The choice of ``maxDiffWinSize`` and ``maxDiffSampSize`` directly affects the runtime for this step.
These two parameters control the maximum difference between an observed window (or sample) size in the dataset
and the reference window size (or sample) size computed in this script. 
Higher allowed differences result in fewer reference values required for the dataset, 
thereby decreasing the computational time needed. As reference, in our analysis we set 
``maxDiffWinSize=50`` and ``maxDiffSampSize=2``.

It is also important to note that our analysis did not explore how increasing the 
values of  ``maxDiffWinSize`` and ``maxDiffSampSize`` affect the estimated evolutionary times.

.. note::
	Make sure that the required parameters described above are correctly defined in the file ``utils/dataset.py``.

- **Output**: 
	One ``.pickle`` file for each reference window size.
	The pickle file contains the following objects:

	1. 	A list of selected evolutionary times for the reference window size;
	2. 	A dictionary where the keys are reference sample sizes associated with the window size,
		and the values are lists of KDE objects, one for each evolutionary time.

	All files are saved in the directory specified in the parameter ``dirSetupEvolTimes``.

.. note::
	Importantly, **the output files can be reused for other datasets**. 
	The dataset serves as a basis for calculating the reference values, 
	but all output files produced in this step are independent of the dataset. 
	Therefore, if another dataset has similar parameter values 
	(such as window size, minimum PCS size, alpha, etc.), the time needed 
	can be significantly decreased by reusing the setup files from previous computations.

Pre-requisites
--------------

	Before using this script, make sure all the required files were computed:

a) Windows for each chromosome in the dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Make sure to run ``2_computeWindows.py`` for **every** chromosome included in 
the attribute ``chromLst`` declared in file ``dataset.py``. 

Cluster resources
-----------------

In case you want to run this Python script stand-alone::

	srun -p compute -t 2:00:00 --mem 20G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash

Otherwise you can use the script ``../cluster/3_setupEvolTimes_runAll.py``
to run this script for the two α parameter values used in our analysis: 

* **α=1.1**, for a model with indels, 
* **α=10**, for substitutions only.

Time, Memory & Disk space
-------------------------

As a point of reference, we provide statistics on runtime, 
memory consumption, and disk space needed to execute this script 
on the 40 vertebrate dataset analyzed in our study.

The run below used **80 cores**.

======  =========  =========  ========
Desc.   Time       Memory     Disk    
======  =========  =========  ========
1.1     02:38:48   29GB       6.290GB 
10      01:50:30   29GB       5.700GB 
======  =========  =========  ========

Time per Run: Details
^^^^^^^^^^^^^^^^^^^^^

For the 40 vertebrate dataset used in our study, and by setting ``maxDiffWinSize=50`` and ``maxDiffSampSize=2``, 
**4,295** pairs of (**reference window size**, **reference sample size**) were computed.

Regarding the reference window size, there were **117** values in total.
Each reference window size is processed on a separate core, so the more 
cores specified in the input, the faster the method will execute.

* Stats on time for all dataset: **~ 1 hour 45 minutes**

**80 cores, 4,295 reference pairs (all pairs)**

===========================================  ========
Step                                         Time (s)
===========================================  ========
Load observed window sizes and PCS counts.     560.36
Compute reference size values.                   1.41
Setup method.                                 5773.33
**Total time**                                6337.82
===========================================  ========

Next, as reference, we also report the time needed for 
**one** reference window size with **336** associated reference sample sizes, 
resulting in 336 pairs of (**reference window size**, **reference sample size**) 
processed and one output file saved.

* Stats on time of a single reference window size: **~9 minutes**

**1 core, 336 reference pairs**

==============================================  =============
Step                                            Time (s)     
==============================================  =============
Select evolutionary times                       43.15036654  
Compute PCS size distribution for evol. times   11.44942665  
Sample PCSs                                     432.17269206 
**Total time**                                  487.30657291 
==============================================  =============

More details on computational time can be found in the log of the run.

Storage per Run: Details
^^^^^^^^^^^^^^^^^^^^^^^^

Total size of output files (**117 files**, one for each reference window size): **~6.3 GB**.

Details of each output file, including the reference window size, file size, and filename:

========  =====    ============================================= 
Win.Ref.  Size     Filename
========  =====    ============================================= 
1000      21M      evolTimes-setup.alpha1.1.winSize1000.pickle
1051      22M      evolTimes-setup.alpha1.1.winSize1051.pickle
1102      23M      evolTimes-setup.alpha1.1.winSize1102.pickle
1153      24M      evolTimes-setup.alpha1.1.winSize1153.pickle
1204      25M      evolTimes-setup.alpha1.1.winSize1204.pickle
1255      25M      evolTimes-setup.alpha1.1.winSize1255.pickle
1306      28M      evolTimes-setup.alpha1.1.winSize1306.pickle
1357      28M      evolTimes-setup.alpha1.1.winSize1357.pickle
1408      29M      evolTimes-setup.alpha1.1.winSize1408.pickle
1459      33M      evolTimes-setup.alpha1.1.winSize1459.pickle
1510      33M      evolTimes-setup.alpha1.1.winSize1510.pickle
1561      38M      evolTimes-setup.alpha1.1.winSize1561.pickle
1612      40M      evolTimes-setup.alpha1.1.winSize1612.pickle
1663      41M      evolTimes-setup.alpha1.1.winSize1663.pickle
1714      43M      evolTimes-setup.alpha1.1.winSize1714.pickle
1765      46M      evolTimes-setup.alpha1.1.winSize1765.pickle
1816      48M      evolTimes-setup.alpha1.1.winSize1816.pickle
1867      50M      evolTimes-setup.alpha1.1.winSize1867.pickle
1918      52M      evolTimes-setup.alpha1.1.winSize1918.pickle
1969      52M      evolTimes-setup.alpha1.1.winSize1969.pickle
2020      52M      evolTimes-setup.alpha1.1.winSize2020.pickle
2071      53M      evolTimes-setup.alpha1.1.winSize2071.pickle
2122      56M      evolTimes-setup.alpha1.1.winSize2122.pickle
2173      57M      evolTimes-setup.alpha1.1.winSize2173.pickle
2224      58M      evolTimes-setup.alpha1.1.winSize2224.pickle
2275      60M      evolTimes-setup.alpha1.1.winSize2275.pickle
2326      61M      evolTimes-setup.alpha1.1.winSize2326.pickle
2377      62M      evolTimes-setup.alpha1.1.winSize2377.pickle
2428      63M      evolTimes-setup.alpha1.1.winSize2428.pickle
2479      64M      evolTimes-setup.alpha1.1.winSize2479.pickle
2530      64M      evolTimes-setup.alpha1.1.winSize2530.pickle
2581      65M      evolTimes-setup.alpha1.1.winSize2581.pickle
2632      65M      evolTimes-setup.alpha1.1.winSize2632.pickle
2683      66M      evolTimes-setup.alpha1.1.winSize2683.pickle
2734      66M      evolTimes-setup.alpha1.1.winSize2734.pickle
2785      67M      evolTimes-setup.alpha1.1.winSize2785.pickle
2836      67M      evolTimes-setup.alpha1.1.winSize2836.pickle
2887      69M      evolTimes-setup.alpha1.1.winSize2887.pickle
2938      70M      evolTimes-setup.alpha1.1.winSize2938.pickle
2989      71M      evolTimes-setup.alpha1.1.winSize2989.pickle
3040      72M      evolTimes-setup.alpha1.1.winSize3040.pickle
3091      75M      evolTimes-setup.alpha1.1.winSize3091.pickle
3142      76M      evolTimes-setup.alpha1.1.winSize3142.pickle
3193      74M      evolTimes-setup.alpha1.1.winSize3193.pickle
3244      79M      evolTimes-setup.alpha1.1.winSize3244.pickle
3295      76M      evolTimes-setup.alpha1.1.winSize3295.pickle
3346      76M      evolTimes-setup.alpha1.1.winSize3346.pickle
3397      75M      evolTimes-setup.alpha1.1.winSize3397.pickle
3448      70M      evolTimes-setup.alpha1.1.winSize3448.pickle
3499      71M      evolTimes-setup.alpha1.1.winSize3499.pickle
3550      68M      evolTimes-setup.alpha1.1.winSize3550.pickle
3601      78M      evolTimes-setup.alpha1.1.winSize3601.pickle
3652      70M      evolTimes-setup.alpha1.1.winSize3652.pickle
3703      70M      evolTimes-setup.alpha1.1.winSize3703.pickle
3754      70M      evolTimes-setup.alpha1.1.winSize3754.pickle
3805      71M      evolTimes-setup.alpha1.1.winSize3805.pickle
3856      73M      evolTimes-setup.alpha1.1.winSize3856.pickle
3907      73M      evolTimes-setup.alpha1.1.winSize3907.pickle
3958      70M      evolTimes-setup.alpha1.1.winSize3958.pickle
4009      68M      evolTimes-setup.alpha1.1.winSize4009.pickle
4060      70M      evolTimes-setup.alpha1.1.winSize4060.pickle
4111      68M      evolTimes-setup.alpha1.1.winSize4111.pickle
4162      71M      evolTimes-setup.alpha1.1.winSize4162.pickle
4213      67M      evolTimes-setup.alpha1.1.winSize4213.pickle
4264      68M      evolTimes-setup.alpha1.1.winSize4264.pickle
4315      67M      evolTimes-setup.alpha1.1.winSize4315.pickle
4366      70M      evolTimes-setup.alpha1.1.winSize4366.pickle
4417      68M      evolTimes-setup.alpha1.1.winSize4417.pickle
4468      63M      evolTimes-setup.alpha1.1.winSize4468.pickle
4520      67M      evolTimes-setup.alpha1.1.winSize4520.pickle
4572      67M      evolTimes-setup.alpha1.1.winSize4572.pickle
4623      65M      evolTimes-setup.alpha1.1.winSize4623.pickle
4675      62M      evolTimes-setup.alpha1.1.winSize4675.pickle
4726      63M      evolTimes-setup.alpha1.1.winSize4726.pickle
4777      62M      evolTimes-setup.alpha1.1.winSize4777.pickle
4828      70M      evolTimes-setup.alpha1.1.winSize4828.pickle
4880      71M      evolTimes-setup.alpha1.1.winSize4880.pickle
4931      65M      evolTimes-setup.alpha1.1.winSize4931.pickle
4982      62M      evolTimes-setup.alpha1.1.winSize4982.pickle
5033      67M      evolTimes-setup.alpha1.1.winSize5033.pickle
5084      62M      evolTimes-setup.alpha1.1.winSize5084.pickle
5138      65M      evolTimes-setup.alpha1.1.winSize5138.pickle
5192      62M      evolTimes-setup.alpha1.1.winSize5192.pickle
5244      62M      evolTimes-setup.alpha1.1.winSize5244.pickle
5310      63M      evolTimes-setup.alpha1.1.winSize5310.pickle
5366      60M      evolTimes-setup.alpha1.1.winSize5366.pickle
5437      52M      evolTimes-setup.alpha1.1.winSize5437.pickle
5488      58M      evolTimes-setup.alpha1.1.winSize5488.pickle
5546      58M      evolTimes-setup.alpha1.1.winSize5546.pickle
5606      48M      evolTimes-setup.alpha1.1.winSize5606.pickle
5670      53M      evolTimes-setup.alpha1.1.winSize5670.pickle
5735      60M      evolTimes-setup.alpha1.1.winSize5735.pickle
5794      52M      evolTimes-setup.alpha1.1.winSize5794.pickle
5851      50M      evolTimes-setup.alpha1.1.winSize5851.pickle
5914      57M      evolTimes-setup.alpha1.1.winSize5914.pickle
5987      38M      evolTimes-setup.alpha1.1.winSize5987.pickle
6041      45M      evolTimes-setup.alpha1.1.winSize6041.pickle
6124      43M      evolTimes-setup.alpha1.1.winSize6124.pickle
6229      25M      evolTimes-setup.alpha1.1.winSize6229.pickle
6280      35M      evolTimes-setup.alpha1.1.winSize6280.pickle
6375      27M      evolTimes-setup.alpha1.1.winSize6375.pickle
6476      27M      evolTimes-setup.alpha1.1.winSize6476.pickle
6553      42M      evolTimes-setup.alpha1.1.winSize6553.pickle
6609      30M      evolTimes-setup.alpha1.1.winSize6609.pickle
6661      17M      evolTimes-setup.alpha1.1.winSize6661.pickle
6819      43M      evolTimes-setup.alpha1.1.winSize6819.pickle
6942      35M      evolTimes-setup.alpha1.1.winSize6942.pickle
7010      14M      evolTimes-setup.alpha1.1.winSize7010.pickle
7068      30M      evolTimes-setup.alpha1.1.winSize7068.pickle
7192      15M      evolTimes-setup.alpha1.1.winSize7192.pickle
7271      22M      evolTimes-setup.alpha1.1.winSize7271.pickle
7338      19M      evolTimes-setup.alpha1.1.winSize7338.pickle
7443      19M      evolTimes-setup.alpha1.1.winSize7443.pickle
7974       7M      evolTimes-setup.alpha1.1.winSize7974.pickle
8032      15M      evolTimes-setup.alpha1.1.winSize8032.pickle
8119      12M      evolTimes-setup.alpha1.1.winSize8119.pickle
8361      22M      evolTimes-setup.alpha1.1.winSize8361.pickle
========  =====    ============================================= 

Function details
----------------

Only relevant functions have been documented below. 
For more details on any function, check the comments in the souce code.

"""

import os
import sys
import re
import pickle
import numpy as np
from collections import defaultdict
import argparse

from mpmath import *

from concurrent import futures

from utils.basicTypes import Pcs, Time
from utils.dataset import Dataset
from utils.indelsModel import IndelsModel, EvolTimesSelector, EvolTimesSolver, getSampleSize

mp.dps    = 30
mp.pretty = True

####################################################################
# Functions involved in setting up the estimation method, i.e.: 
# (1) Selection of evolutionary times;
# (2) Computation of expected PCS size distributions;
# (3) Sampling M times N points for each the PCS size distribution;
# (4) Initialization of a KDE object based on the sample.
def selectEvolTimes(model, maxNbTaus):
	selector = EvolTimesSelector(model)
	t_span   = (0.0, 1.0)
	seltaus  = selector.find_ts(t_span,maxNbTaus)
	seltaus  = seltaus[1:] # Ignore t=0
	return seltaus

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
	winSizeRef, sampleSizeRef_lst, alpha, overwriteFiles, debugMessage, my_dataset = parallelInput

	# Save times for each step and overall time.
	timeTrack = Time()
	timeTrack.start()

	outFilename = my_dataset.getOutFilename_setupEvolTimes(winSizeRef,alpha)
	if ((not overwriteFiles) and os.path.isfile(outFilename)):
		print(f"WARNING! Setup file for α={alpha} and reference window size={winSizeRef} was already computed. Skipping computation...")
	else:
		model  = IndelsModel(winSizeRef, alpha, my_dataset.minPCSsize)

		# Select evolutionary times.
		if(debugMessage): timeTrack.startStep(f"{debugMessage} Select evolutionary times")
		ts_sel = selectEvolTimes(model, my_dataset.maxNbTaus)
		if(debugMessage): timeTrack.stopStep()

		# Sample PCS size distribution and create a 'KDE object' 
		# for each combination of [winRefSize][sampleRefSize].
		# KDE = Kernel Density Estimation. 
		# KDE computes the probability that a point was sampled 
		# from the same distribution as a group of points.
		t_PCSsizes = model.sizes
		setupInfo  = {}
		if(debugMessage): timeTrack.startStep(f"{debugMessage} Sample PCSs")
		# Setup for each reference sample size.
		for sampleIdx, sampleSizeRef in enumerate(sampleSizeRef_lst):
			solver = EvolTimesSolver(ts_sel)
			solver.initialize_ts(model, my_dataset.nbSamplesPerTau, sampleSizeRef, my_dataset.minDiffPcsSizes)
			setupInfo[sampleSizeRef] = solver
		if(debugMessage): 
			timeTrack.stopStep()
			print(f" Average time per sample size = {timeTrack.t_local/len(sampleSizeRef_lst)}")


		# Save results for all sample sizes for a given window size.
		with open(outFilename, 'wb') as pickleFile:
			pickle.dump((ts_sel, setupInfo),pickleFile, protocol=pickle.HIGHEST_PROTOCOL)

	timeTrack.stop()
	return (winSizeRef, timeTrack)

########################################
# Functions related to reference values 
# for window size and sample size, based on real data.

def loadObservedValues(my_dataset, qChromLst):
	""" This function loads all computed windows and returns the observed
	window sizes and sample sizes (i.e. the total count of PCSs in a window).
	"""
	obsVals = defaultdict(int)
	speciesUCSCnames = my_dataset.speciesUCSCnames
	# Get all observed window sizes and sample sizes (= total count of PCSs).
	for qChrom in qChromLst:
		print(f"\t[{qChrom}] Loading observed values.")
		# Check if input file exists.
		windows_filename = my_dataset.getOutFilename_computeWindows(qChrom)
		if (not os.path.isfile(windows_filename)):
			print(f"ERROR! File not found: {windows_filename}")
			sys.exit()
		# Load windows.
		pcs_distrib_all = pickle.load(open(windows_filename, 'rb'))
		for ucscName_other in speciesUCSCnames:
			if(ucscName_other not in pcs_distrib_all):
				print(f"ERROR! Species not found: {ucscName_other} ({qChrom}).")
				sys.exit()
			# Get the PCS size distribution of the species/chromosome.
			pcs_distrib_allwin = pcs_distrib_all[ucscName_other]
			for (windowId, PCSsizeDistrib) in pcs_distrib_allwin.items():
				winBegPos, winEndPos = windowId
				winSizeObs  = winEndPos-winBegPos
				sampSizeObs = getSampleSize(PCSsizeDistrib)
				if(sampSizeObs > 0): obsVals[(winSizeObs,sampSizeObs)] += 1 
	return obsVals

def getReferenceValues(my_dataset, obsVals):

	""" This function computes the reference values for window size and 
	sample size, based on the real data from the dataset.

	For each window in the chromosome of the reference species, there 
	exists a PCS size distribution linked to that window, which varies 
	based on the other species included in the analysis.

	Every window has two associated values: **window size** and **sample size**.
	**Window size** is defined as the difference between the start and 
	end coordinates of the window, while **sample size** is the total 
	number of PCSs observed in that window.

	Due to the high computational cost of estimating evolutionary time 
	for all possible combinations of window sizes and sample sizes in 
	the dataset, we use this clustering function.

	The clustering function is implemented as follows. 

	Each window is assigned to a cluster, which is defined by 
	a **reference window size** and a **reference sample size**.
	Clusters are defined such that the window sizes within a cluster 
	cannot differ from the reference window size by more than ``maxDiffWinSize``. 
	There is also an upper limit on the difference between the 
	reference and observed sample sizes, set by ``maxDiffSampSize``.
	"""

	# Parameters to determine reference values.
	maxDiffWinSize  = my_dataset.maxDiffWinSize
	maxDiffSampSize = my_dataset.maxDiffSampSize
	minSampRefSize  = my_dataset.minSampRefSize
	maxSampRefSize  = my_dataset.maxSampRefSize

	refVals     = {}
	refVals_lst = []
	obsVals_lst = sorted(list(obsVals.keys()), key=lambda val: (val[0], val[1]))
	for (winSizeObs,sampSizeObs) in obsVals_lst:
		# Find reference value for each observed value.
		refKey  = None
		# Ignore window if the number of PCSs in the window is too low.
		if(sampSizeObs < minSampRefSize): continue
		# Adjust sample size.
		sampSizeObsAdj = min(max(sampSizeObs, minSampRefSize), maxSampRefSize)
		# The list is sorted by window size, so the list of 
		# reference values will also be sorted by window sizes.
		for (w,s) in reversed(refVals_lst):
			if (abs(winSizeObs-w) <= maxDiffWinSize):
				if (abs(sampSizeObsAdj-s) <= maxDiffSampSize):
					refKey = (w,s)
					break
				else:
					refKey = (w, None)
			else:
				break
		# Add a new reference value in case the existing ones
		# cannot be linked to the current observed value.
		if((refKey == None) or (refKey[1] == None)):
			refKey = (winSizeObs,sampSizeObsAdj) if (refKey == None) else (refKey[0],sampSizeObsAdj)
			refVals_lst.append(refKey)
		refVals[(winSizeObs,sampSizeObs)] = refKey
	return refVals

def checkPCSsizeDistribFiles(dirWindows, windowSize, prefixQuery, qChromLst):
	""" This function checks whether the PCS size distributions 
	have been computed for every window across all chromosomes in the dataset.
	"""
	for qChrom in qChromLst:
		# Check if input files exist.
		windows_filename = os.path.join(dirWindows, f"{prefixQuery}.{qChrom}.{windowSize}.windows.pickle")
		if (not os.path.isfile(windows_filename)):
			print(f"ERROR! PCS size distribution not found for chromosome '{qChrom}' (filepath: {windows_filename}).")
			sys.exit()

####################################################################
# Functions related to parallelizing computation.
def initParallelInputs(refVals, alpha, overwriteFiles, nbcores, my_dataset):
	refVals_lst = sorted(list(set(refVals.values())), key=lambda val: (val[0], val[1]))
	print(f"- Total reference values: {len(refVals_lst)}")
	# Group values by reference window size.
	parallelInputs = defaultdict(list)
	for (winSizeRef,sampSizeRef) in refVals_lst:
		parallelInputs[winSizeRef].append(sampSizeRef)
	# Create list of parallel inputs.
	parallelInputs_lst = []
	winSizeRef_lst = sorted(list(parallelInputs))
	nbchars = len(str(len(winSizeRef_lst)))
	for runIdx, winSizeRef in enumerate(winSizeRef_lst):
		debugMessage = f"[{str(runIdx+1).rjust(nbchars)}/{len(winSizeRef_lst)}] " if ((runIdx % nbcores) == 0) else ""
		parallelInputs_lst.append((winSizeRef, parallelInputs[winSizeRef], alpha, overwriteFiles, debugMessage, my_dataset))
	return parallelInputs_lst

def processIteration(cntRuns, stepRun, totRuns, result, it_timeTrack_max):
	winSizeRef, it_timeTrack = result
	if((it_timeTrack_max == None) or (it_timeTrack_max.t_global < it_timeTrack.t_global)):
		it_timeTrack_max = it_timeTrack
	if((cntRuns % stepRun) == 0): 
		print(f"[Max runtime in the last {stepRun} iterations]")
		it_timeTrack_max.print()
		it_timeTrack_max = None
	return cntRuns+1, it_timeTrack_max

####################################
# MAIN.
####################################
if (__name__ == '__main__'):

	parser = argparse.ArgumentParser(description="Extract PCSs (Perfectly Conserved Sequences) from Chains (ordered aligned blocks).")
	parser.add_argument("-alpha", help="It determines how frequent longer indels can occur. The parameter alpha can take any value above 1: (1, ∞). If alpha is near 1, larger indels are more likely to occur. If alpha is above 5 (a hard upper limit set internally with [max_alpha]), the model turns into a substitution only model.", type=float, required=True)
	parser.add_argument("-cores", help="Number of cores to be used during setup.", type=int, required=True)
	parser.add_argument("--overwrite", action='store_true', help="If this flag is specified, any existing output files will be overwritten during the run.")

	args       = parser.parse_args()
	my_dataset = Dataset()

	alpha   = float(args.alpha)
	nbcores = int(args.cores)
	overwriteFiles = args.overwrite

	windowSize = my_dataset.windowSize   # Window size (in number of base pairs).
	minPCSsize = my_dataset.minPCSsize	 # Minimum size (in number of base pairs) for the PCS to be considered.

	prefixQuery	     = my_dataset.refsp_ucscName     # In our study, query is always the human genome (hg38)
	qChromLst        = my_dataset.chromLst           # All chromosomes from query genome (chr1, chr2, ..., chrY)
	speciesUCSCnames = my_dataset.speciesUCSCnames

	maxDiffWinSize   = my_dataset.maxDiffWinSize
	maxDiffSampSize  = my_dataset.maxDiffSampSize
	minSampRefSize   = my_dataset.minSampRefSize
	maxSampRefSize   = my_dataset.maxSampRefSize
	nbSamplesPerTau  = my_dataset.nbSamplesPerTau

	dirOut		     = my_dataset.dirSetupEvolTimes  # Directory where pre-computed info will be saved (output directory).
	dirWin           = my_dataset.dirWindows

	print("***************************************************")
	print("*      Evolutionary Time Inference: Setup         *")
	print("* Pre-compute the Kernel Density Estimation (KDE) *")
	print("* for various window sizes and PCS counts (sample *")
	print("* sizes) derived from the dataset.                *")
	print("***************************************************")
	print(f" Core parameter values:")
	print(f" ----------------------")
	print(f" - Query genome: {prefixQuery}")
	print(f" - Window size: {windowSize}")
	print(f" - Minimum PCS size: {minPCSsize}")
	print(f" - Number of cores: {nbcores}")
	print(f" Setup-specific parameter values:")
	print(f" --------------------------------")
	print(f" - Input directory with windows and their PCS size distribution: {dirWin}")
	print(f" - Propensity for indels: α={alpha}")
	print(f" - Number of sampled PCS size distributions per evolutionary time: {nbSamplesPerTau}")
	print(f" - Constraint on reference window sizes: max(|[obs. window size]-[ref. window size]|) ≤ {maxDiffWinSize}")
	print(f" - Constraint on reference sample sizes: max(|[obs. sample size]-[ref. sample size]|) ≤ {maxDiffSampSize}")
	print(f" - Constraint on reference sample sizes: {minSampRefSize} ≤ [ref. sample size]| ≤ {maxSampRefSize}")
	print(f" - Output directory for pre-computed data: {dirOut}")
	print(f" - Log directory with details of the run (output): {my_dataset.dirLog}")
	print(f" - Overwrite output files? {overwriteFiles}\n\n")

	timeTrack = Time()
	timeTrack.start()

	#######################################
	# Check pre-requisites.

	# Check whether the PCS size distributions have been computed for every 
	# window across all chromosomes and genomes in the dataset.
	checkPCSsizeDistribFiles(dirWin, windowSize, prefixQuery, qChromLst)

	# Load observed window sizes and sample sizes (total count of PCSs).
	timeTrack.startStep("Load observed window sizes and PCS counts.")
	obsVals = loadObservedValues(my_dataset, qChromLst) # ~ 10 minutes
	timeTrack.stopStep()
	# Compute reference values for window and sample sizes.
	timeTrack.startStep("Compute reference size values.")
	refVals = getReferenceValues(my_dataset, obsVals)
	timeTrack.stopStep()

	#######################################
	# Compute evolutionary times for reference values.
	timeTrack.startStep("Setup method.")
	parallelInputs  = initParallelInputs(refVals, alpha, overwriteFiles, nbcores, my_dataset)
	print(f"Reference window sizes={len(parallelInputs)}; Cores to process window sizes={nbcores}.")
	cntRuns, stepRun, totRuns = 0, int(len(parallelInputs)*0.1), len(parallelInputs)
	it_timeTrack_max = None
	if(nbcores == 1):		
		for parallelInput in parallelInputs:
			result = computeSetupInfo_winRefSize(parallelInput)
			cntRuns, it_timeTrack_max = processIteration(cntRuns, stepRun, totRuns, result, it_timeTrack_max)
	else:
		with futures.ProcessPoolExecutor(nbcores) as pool:
			for result in pool.map(computeSetupInfo_winRefSize, parallelInputs):
				cntRuns, it_timeTrack_max = processIteration(cntRuns, stepRun, totRuns, result, it_timeTrack_max)
	timeTrack.stopStep()

	timeTrack.stop()
	timeTrack.print()
	print("Done!")
