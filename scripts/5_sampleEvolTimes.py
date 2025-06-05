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

a) Window file for all chromosomes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Make sure to run ``2_computeWindows.py`` for all chromosomes of the reference genome 
(human in our case: *hg38*).

b) Files with estimated evolutionary times per window
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

======================  =========  =========  ========  =========  =========  ========
Desc.                         Parameter α=1.1                 Parameter α=10          
----------------------  ------------------------------  ------------------------------
UCSC name               Time       Memory     Disk      Time       Memory     Disk    
======================  =========  =========  ========  =========  =========  ========
panPan3                 10:38:06   24GB       0.009GB   07:05:55   26GB       0.009GB 
panTro6                 10:26:52   15GB       0.009GB   07:19:46   13GB       0.009GB 
gorGor6                 10:19:39   14GB       0.007GB   07:10:25   14GB       0.007GB 
ponAbe3                 08:46:34   10GB       0.004GB   06:56:24   10GB       0.004GB 
papAnu4                 08:13:30   9GB        0.004GB   07:26:07   9GB        0.004GB 
macFas5                 08:21:27   9GB        0.004GB   06:00:43   9GB        0.004GB 
rhiRox1                 08:03:42   9GB        0.004GB   06:37:15   9GB        0.004GB 
chlSab2                 08:46:12   10GB       0.003GB   03:27:19   9GB        0.000GB 
nasLar1                 06:54:55   9GB        0.003GB   05:35:06   9GB        0.004GB 
rheMac10                08:05:15   9GB        0.003GB   06:23:03   9GB        0.003GB 
calJac4                 07:15:16   9GB        0.003GB   06:27:27   9GB        0.003GB 
tarSyr2                 07:39:09   15GB       0.003GB   06:03:52   10GB       0.003GB 
micMur2                 06:10:46   10GB       0.003GB   05:20:00   10GB       0.003GB 
galVar1                 06:29:39   13GB       0.003GB   05:32:54   10GB       0.003GB 
mm39                    04:57:12   8GB        0.003GB   03:59:24   9GB        0.003GB 
oryCun2                 06:19:43   9GB        0.003GB   04:41:00   8GB        0.003GB 
rn7                     05:18:09   7GB        0.003GB   03:58:33   7GB        0.003GB 
vicPac2                 06:45:18   14GB       0.003GB   05:22:36   9GB        0.003GB 
bisBis1                 06:04:58   9GB        0.003GB   04:57:54   10GB       0.003GB 
felCat9                 06:48:25   10GB       0.003GB   06:04:21   9GB        0.003GB 
manPen1                 06:49:04   9GB        0.003GB   05:01:36   9GB        0.003GB 
bosTau9                 05:33:34   9GB        0.003GB   04:01:31   9GB        0.003GB 
canFam6                 06:38:34   12GB       0.003GB   05:13:01   9GB        0.003GB 
musFur1                 07:07:16   9GB        0.003GB   06:25:33   9GB        0.003GB 
neoSch1                 07:40:54   10GB       0.003GB   05:37:21   10GB       0.003GB 
equCab3                 07:23:28   10GB       0.003GB   05:45:48   10GB       0.003GB 
myoLuc2                 05:56:58   14GB       0.003GB   04:13:25   8GB        0.003GB 
susScr11                06:14:24   9GB        0.003GB   04:37:46   9GB        0.003GB 
enhLutNer1              08:25:34   9GB        0.003GB   06:03:54   9GB        0.003GB 
triMan1                 06:48:24   9GB        0.003GB   05:12:03   9GB        0.003GB 
macEug2                 02:04:03   7GB        0.003GB   01:38:45   7GB        0.003GB 
ornAna2                 01:39:16   7GB        0.002GB   01:20:48   7GB        0.002GB 
aptMan1                 01:25:54   6GB        0.002GB   01:06:15   7GB        0.002GB 
galGal6                 01:09:00   6GB        0.002GB   00:58:13   7GB        0.002GB 
thaSir1                 01:06:57   8GB        0.002GB   00:53:36   8GB        0.002GB 
aquChr2                 01:19:43   8GB        0.002GB   01:05:13   6GB        0.002GB 
melGal5                 01:19:22   8GB        0.002GB   01:04:42   6GB        0.002GB 
xenLae2                 00:58:37   7GB        0.002GB   00:48:54   7GB        0.002GB 
xenTro10                01:00:18   7GB        0.002GB   00:49:39   7GB        0.002GB 
danRer11                00:50:30   6GB        0.002GB   00:44:39   8GB        0.002GB 
======================  =========  =========  ========  =========  =========  ========

Time per Run: Details
^^^^^^^^^^^^^^^^^^^^^

Stats on time of a single run (human-mouse alignment, all chromosomes, 80 cores): **~3 hours 18 minutes**
Details on computational time are available in the log of the run.

==========================================  =========
Step                                         Time (s)
==========================================  =========
[chr1] Load data (obs. PCSs + estimates)        45.22
[chr1] Sample results                          967.04
[chr2] Load data (obs. PCSs + estimates)        49.74
[chr2] Sample results                         1050.78
[chr3] Load data (obs. PCSs + estimates)        33.68
[chr3] Sample results                          717.64
[chr4] Load data (obs. PCSs + estimates)        40.90
[chr4] Sample results                          614.35
[chr5] Load data (obs. PCSs + estimates)        34.67
[chr5] Sample results                          756.51
[chr6] Load data (obs. PCSs + estimates)        33.89
[chr6] Sample results                          722.64
[chr7] Load data (obs. PCSs + estimates)        26.56
[chr7] Sample results                          566.90
[chr8] Load data (obs. PCSs + estimates)        29.50
[chr8] Sample results                          599.79
[chr9] Load data (obs. PCSs + estimates)        22.56
[chr9] Sample results                          452.59
[chr10] Load data (obs. PCSs + estimates)       24.31
[chr10] Sample results                         491.95
[chr11] Load data (obs. PCSs + estimates)       28.77
[chr11] Sample results                         599.64
[chr12] Load data (obs. PCSs + estimates)       24.69
[chr12] Sample results                         519.79
[chr13] Load data (obs. PCSs + estimates)       21.11
[chr13] Sample results                         362.63
[chr14] Load data (obs. PCSs + estimates)       19.26
[chr14] Sample results                         414.56
[chr15] Load data (obs. PCSs + estimates)       17.94
[chr15] Sample results                         349.07
[chr16] Load data (obs. PCSs + estimates)       19.27
[chr16] Sample results                         323.72
[chr17] Load data (obs. PCSs + estimates)       14.70
[chr17] Sample results                         256.59
[chr18] Load data (obs. PCSs + estimates)       16.56
[chr18] Sample results                         325.18
[chr19] Load data (obs. PCSs + estimates)       12.01
[chr19] Sample results                         166.99
[chr20] Load data (obs. PCSs + estimates)       14.11
[chr20] Sample results                         310.32
[chr21] Load data (obs. PCSs + estimates)        8.08
[chr21] Sample results                         156.26
[chr22] Load data (obs. PCSs + estimates)        8.82
[chr22] Sample results                         142.33
[chrX] Load data (obs. PCSs + estimates)        28.60
[chrX] Sample results                          403.82
[chrY] Load data (obs. PCSs + estimates)         7.52
[chrY] Sample results                           30.64
**Total time**                               11884.55
==========================================  =========

Storage per Run: Details
^^^^^^^^^^^^^^^^^^^^^^^^

Total size of output files (**80 files**, one for each pairwise alignment, given α=1.1 and α=10.0): **~261 MB**.

Details of each output file (α=1.1), including the file size and filename:

===========  =======    ============================================= 
UCSC name    Size       Filename
===========  =======    ============================================= 
panPan3      9.2 MB     pcsDistrib-samp.panPan3.alpha1.1.pickle
panTro6      9.3 MB     pcsDistrib-samp.panTro6.alpha1.1.pickle
gorGor6      6.9 MB     pcsDistrib-samp.gorGor6.alpha1.1.pickle
ponAbe3      4.2 MB     pcsDistrib-samp.ponAbe3.alpha1.1.pickle
papAnu4      3.8 MB     pcsDistrib-samp.papAnu4.alpha1.1.pickle
macFas5      3.7 MB     pcsDistrib-samp.macFas5.alpha1.1.pickle
rhiRox1      3.7 MB     pcsDistrib-samp.rhiRox1.alpha1.1.pickle
chlSab2      3.4 MB     pcsDistrib-samp.chlSab2.alpha1.1.pickle
nasLar1      3.5 MB     pcsDistrib-samp.nasLar1.alpha1.1.pickle
rheMac10     3.4 MB     pcsDistrib-samp.rheMac10.alpha1.1.pickle
calJac4      3.2 MB     pcsDistrib-samp.calJac4.alpha1.1.pickle
tarSyr2      3.0 MB     pcsDistrib-samp.tarSyr2.alpha1.1.pickle
micMur2      3.0 MB     pcsDistrib-samp.micMur2.alpha1.1.pickle
galVar1      2.9 MB     pcsDistrib-samp.galVar1.alpha1.1.pickle
mm39         2.8 MB     pcsDistrib-samp.mm39.alpha1.1.pickle
oryCun2      2.9 MB     pcsDistrib-samp.oryCun2.alpha1.1.pickle
rn7          2.8 MB     pcsDistrib-samp.rn7.alpha1.1.pickle
vicPac2      2.9 MB     pcsDistrib-samp.vicPac2.alpha1.1.pickle
bisBis1      2.9 MB     pcsDistrib-samp.bisBis1.alpha1.1.pickle
felCat9      2.9 MB     pcsDistrib-samp.felCat9.alpha1.1.pickle
manPen1      2.8 MB     pcsDistrib-samp.manPen1.alpha1.1.pickle
bosTau9      2.8 MB     pcsDistrib-samp.bosTau9.alpha1.1.pickle
canFam6      2.9 MB     pcsDistrib-samp.canFam6.alpha1.1.pickle
musFur1      2.9 MB     pcsDistrib-samp.musFur1.alpha1.1.pickle
neoSch1      2.9 MB     pcsDistrib-samp.neoSch1.alpha1.1.pickle
equCab3      2.9 MB     pcsDistrib-samp.equCab3.alpha1.1.pickle
myoLuc2      2.7 MB     pcsDistrib-samp.myoLuc2.alpha1.1.pickle
susScr11     2.9 MB     pcsDistrib-samp.susScr11.alpha1.1.pickle
enhLutNer1   2.9 MB     pcsDistrib-samp.enhLutNer1.alpha1.1.pickle
triMan1      2.8 MB     pcsDistrib-samp.triMan1.alpha1.1.pickle
macEug2      2.6 MB     pcsDistrib-samp.macEug2.alpha1.1.pickle
ornAna2      2.6 MB     pcsDistrib-samp.ornAna2.alpha1.1.pickle
aptMan1      2.6 MB     pcsDistrib-samp.aptMan1.alpha1.1.pickle
galGal6      2.5 MB     pcsDistrib-samp.galGal6.alpha1.1.pickle
thaSir1      2.3 MB     pcsDistrib-samp.thaSir1.alpha1.1.pickle
aquChr2      2.6 MB     pcsDistrib-samp.aquChr2.alpha1.1.pickle
melGal5      2.5 MB     pcsDistrib-samp.melGal5.alpha1.1.pickle
xenLae2      2.1 MB     pcsDistrib-samp.xenLae2.alpha1.1.pickle
xenTro10     2.1 MB     pcsDistrib-samp.xenTro10.alpha1.1.pickle
danRer11     1.7 MB     pcsDistrib-samp.danRer11.alpha1.1.pickle
===========  =======    ============================================= 

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
	winSizeRef, my_dataset, alpha, ts, windows = parallelInput

	nbSamplesPerWin = my_dataset.nbSamplesPerWin
	nbSamplesTotal  = nbSamplesPerWin*3

	winSizeRef_prev = -1
	model  = None # Model depends on the reference window size.
	solver = EvolTimesSolver(None)

	PCSdistrib_obs_win = defaultdict(int)
	PCSdistrib_est_win = [defaultdict(float) for sampleIdx in range(nbSamplesPerWin)]
	taudistrib_est_win = [defaultdict(float) for sampleIdx in range(nbSamplesPerWin)]

	# Initialize indel model with reference window size.
	model = IndelsModel(winSizeRef, alpha, my_dataset.minPCSsize)

	for stepIdx, (windowId, winPcsObs, winEvolTimes) in enumerate(windows):

		# Sample evolutionary times from window posterior.
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
	return (winSizeRef, PCSdistrib_obs_win, PCSdistrib_est_win, taudistrib_est_win)

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
	parallelInputs=defaultdict(list)
	for (windowId, winEvolTimes) in winIds_ests.items():
		winSizeRef, sampleSizeRef = winIds_refs[windowId]
		winPcsObs = pcs_distrib_allwin[windowId]
		parallelInputs[(winSizeRef,0)].append((windowId, winPcsObs, winEvolTimes))

	# Re-balance windows per core.
	# Windows are grouped by their reference window size, and smaller 
	# reference window sizes tend to have many more windows associated than 
	# bigger reference window sizes.
	windowsPerJob  = max(int(np.quantile([len(parallelInput) for parallelInput in parallelInputs.values()], 0.50)),30)
	winSizeRef_lst = list(sorted(list(parallelInputs.keys()), key=lambda dictKey: dictKey[0]))
	print(f"Ideal number of windows per core = {windowsPerJob} windows. ")
	for (winSizeRef, lstcount) in winSizeRef_lst:
		windows = parallelInputs[(winSizeRef,lstcount)]
		win_nb  = len(windows)
		print(f"\t{winSizeRef} : {win_nb} windows.")
		if(win_nb > windowsPerJob):
			windows_lst = splitList(windows, int(1+win_nb/windowsPerJob))
			print(f"\t\t {winSizeRef} : Split into {len(windows_lst)} lists.")
			for windows_idx, windows_new in enumerate(windows_lst):
				parallelInputs[(winSizeRef,windows_idx)]=windows_new
	timeTrack.stopStep()

	# Split windows according to their reference window size across available cores.
	print(f"Windows w/estimates={len(winIds_ests)}; Cores to process windows={nbcores}.")
	parallelInputs = [(winSizeRef, my_dataset, alpha, winSizeRefs_ts[winSizeRef], windows) for ((winSizeRef,lstcount), windows) in parallelInputs.items()]

	timeTrack.stop()
	timeTrack.print()
	return parallelInputs

def processIteration(chrom, result, info_saved):
	PCSdistrib_obs_whl, PCSdistrib_est_whl, taudistrib_est_whl, PCSdistrib_obs_chr, PCSdistrib_est_chr, taudistrib_est_chr = info_saved
	winSizeRef, PCSdistrib_obs_win, PCSdistrib_est_win, taudistrib_est_win = result

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
	return (PCSdistrib_obs_whl, PCSdistrib_est_whl, taudistrib_est_whl, PCSdistrib_obs_chr, PCSdistrib_est_chr, taudistrib_est_chr)

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
					winSizeRef = result[0]
					info_saved = processIteration(qChrom, result, info_saved)
		timeTrack.stopStep()
	
	# Save results for pairwise alignment.
	with open(outputFilepath, 'wb') as pickleFile:
		# It saves sampled evolutionary times and PCSs for each chromosome and also whole-genome.
		pickle.dump(info_saved, pickleFile, protocol=pickle.HIGHEST_PROTOCOL)

	timeTrack.stop()
	timeTrack.print()
	print(f"Done!")
