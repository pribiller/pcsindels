"""Plot ``Figure 3 (panel B)`` from the paper, i.e., a comparison
of evolutionary estimates based on pairwise alignments (ours) 
versus those based on phylogenetic analysis 
(Kuderna et al. [2023] and Upham et al. [2019]).

The comparison depends on how many species from our 40-vertebrate 
dataset overlap with those from other studies. Specifically, 
there are 15 species shared with the dataset from Kuderna et al. (2023), 
and 31 species shared with the dataset from Upham et al. (2019).

In the two generated plots, our estimates are shown in the y-axis,
whereas those from other studies are shown in the x-axis.
Each dot corresponds to a species shared by both datasets, 
with the color of the dot indicating its divergence time from humans, 
following the color code specified in the paper in Figure 3, panel A 
(green for primates, yellow for mammals, and red for vertebrates).

For more information about the two studies used in this comparison, please check:

  * A global catalog of whole-genome diversity from 233 primate species (Kuderna et al., Science, 2023);
  * Inferring the mammal tree: species-level sets of phylogenies for questions in ecology, evolution, and conservation (Upham et al., PLoS Biol, 2019).

- **Use**::
	
	python3 6_2_plotFig_EvolTimesComp.py

- **Example of Usage**::

	python3 ~/code/6_2_plotFig_EvolTimesComp.py

- **Input Parameter**:

To ensure the graphs match those used in the paper,
the parameters are hard-coded in the script and 
cannot be modified via command line.

Pre-requisites
--------------

Before using this script, make sure all the required files were pre-computed:

a) Files with sampled evolutionary times
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Make sure to run ``5_sampleEvolTimes.py`` for **α=10.0**.

b) Logs from evolutionary time estimates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Make sure to keep the logs from ``4_estimateEvolTimes.py`` for **α=10.0**. It contains the information regarding windows without estimates.

Time, Memory & Disk space
-------------------------

Running the script on a single core takes **22.44 seconds** and requires a small amount of memory. All output files take very few space (< 25 kb).

**Output files**:
	1. ``evolTimes-compKuderna2023.pdf``: The PDF output file containing the comparison with the estimates from Kuderna et al. (2023);
	2. ``evolTimes-compUpham2019-corr.pdf``: The PDF output file containing the (adjusted) comparison with the estimates from Upham et al. (2019);
	3. ``evolTimes-compUpham2019.pdf``: The PDF output file containing the (non-adjusted) comparison with the estimates from Upham et al. (2019) (inset plot);

Function details
----------------

Only relevant functions have been documented below. 
For more details on any function, check the comments in the souce code.

"""

import os
import pickle
import random
from collections import defaultdict
import re
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.backends.backend_pdf import PdfPages
from colorir import Grad
from Bio import Phylo

from utils.dataset import Dataset
from utils.basicTypes import Time

###############################
# Methods to load info from other studies.

def parseMapping(mappingFilename, UCSConly=True):
	mappingSpecies = {}
	with open(mappingFilename) as f:
		for line in f:
			m = re.match(r'(.*);(.*);(.*);(.*);(.*);(.*);(.*);', line)
			if m:
				prefix	        = m.group(1)
				speciesName     = m.group(2)
				labelUpham2019	= m.group(3)
				labelIrisarri2017_raxml = m.group(4)
				labelIrisarri2017_phylobayes = m.group(5)
				labelKuderna2023  = m.group(6) # Science paper
				labelKuderna2024  = m.group(7) # Nature paper
				if(UCSConly and (prefix == "NONE")):
					continue
				prefix = prefix if (prefix != "NONE") else labelKuderna2024
				mappingSpecies[prefix] = (speciesName, {"Upham2019": labelUpham2019, 
														"Irisarri2017_raxml": labelIrisarri2017_raxml, 
														"Irisarri2017_phylobayes": labelIrisarri2017_phylobayes,
														"Kuderna2023": labelKuderna2023,
														"Kuderna2024": labelKuderna2024})
	return mappingSpecies

def getDistances(mappingSpecies, tree, refName, paperRef):
	""" Get distance (i.e., the sum of branch lengths in the path joining two species),
	between humans (reference species) and all other species that also appear in our dataset.
	"""
	alldists = {}
	for prefix, (speciesName, labels) in mappingSpecies.items():
		label = labels[paperRef]
		if(label != "NONE"):
			#print(f"{prefix} {speciesName} {label} {refName}")
			mrca   = tree.common_ancestor({"name": refName}, {"name": label})
			distAB = tree.distance({"name": refName},{"name": label})
			distA  = tree.distance(mrca,{"name": label})
			distB  = tree.distance(mrca,{"name": refName})
			#print(f"{distAB} {distA+distB} {distA} {distB}")
			alldists[prefix] = distAB # distA,distB,distAB,abs((distAB/2.0)-distA),abs((distAB/2.0)-distB)
	return alldists

def loadData(my_dataset, mappingSpecies, paperRef):

	# Reference names.	
	refNames = {"Upham2019"   : "Homo_sapiens_HOMINIDAE_PRIMATES",
				"Kuderna2024" : "Homo_sapiens",
				"Kuderna2023" : "Homo_sapiens",
				"Irisarri2017_raxml" : "Homo_sapiens",
				"Irisarri2017_phylobayes" : "Homo_sapiens"}

	# Newick trees.	
	treeFilenames = {"Upham2019"  : "Upham2019_RAxML_bipartitions.result_FIN4_raw_rooted_wBoots_4098mam1out_OK.newick",
					"Kuderna2024" : "Kuderna2024_447-mammalian-2022v1.nh",
					"Kuderna2023" : "Kuderna2023_science.abn7829_data_s3.nw.tree",
					"Irisarri2017_raxml" : "Irisarri2017_misgen_50-RAXML-PROTGAMMAGTR-100xRAPIDBP.tre",
					"Irisarri2017_phylobayes" : "Irisarri2017_jack50000-0DP-all-CATG4.con.ann.tre"}

	treeFilename = os.path.join(my_dataset.dirTrees,treeFilenames[paperRef])
	tree    = Phylo.read(treeFilename, "newick")
	refName = refNames[paperRef]
	return getDistances(mappingSpecies, tree,  refName,  paperRef)

###############################
# Methods to load info from our study.

def meanEvolTimes(taudistrib_est_onespecies, empty_win_info=(-1,-1), empty_win_tau=-1):
	mean_vals       = []
	empty_win_cnt, total_win_cnt = empty_win_info
	empty_windows   = [empty_win_tau for i in range(empty_win_cnt)] if (empty_win_cnt > 0) else []
	# For each sampled collection of evolutionary times (whole-genome level, every 
	# window with "enough" information has a sampled evolutionary time associated with).
	for taudistrib_est in taudistrib_est_onespecies:
		# Gather taus from all windows with estimates.
		all_taus  = [tauval for (tauval, taucnt) in taudistrib_est.items() for idx in range(int(taucnt))]
		if ((empty_win_cnt > 0) and (len(all_taus) != (total_win_cnt-empty_win_cnt))):
			print(f"WARNING! Expected number of windows with estimates: {(total_win_cnt-empty_win_cnt)}; Found: {len(all_taus)}.")
		all_taus.extend(empty_windows)
		# Bootstrap for one sample of evolutionary times.
		bootstrap_per_sample = 5
		bootstrap_size       = 400
		for bootstrap_idx in range(bootstrap_per_sample):
			# Select random evolutionary times from whole-genome sample.
			evoltimes_rdm = random.sample(all_taus, bootstrap_size)
			# Compute centrality measure.
			mean_vals.append(np.mean(evoltimes_rdm))
	return mean_vals

def loadOurData(alpha, my_dataset, empty_windows=None, empty_mu=-1):
	dists = {}
	for UCSCname in my_dataset.speciesUCSCnames:

		# Load sampled distributions.
		sampEvolTimesFilepath = my_dataset.getOutFilename_sampleEvolTimes(UCSCname, alpha)
		if(not os.path.isfile(sampEvolTimesFilepath)):
			print(f"[{UCSCname}] WARNING! File not found: {sampEvolTimesFilepath}. Skipping computation!")
			continue
		PCSdistrib_obs_whl, PCSdistrib_est_whl, taudistrib_est_whl, PCSdistrib_obs_chr, PCSdistrib_est_chr, taudistrib_est_chr = pickle.load(open(sampEvolTimesFilepath, 'rb'))
		# Compute mean evolutionary time for a species.
		if ((not empty_windows) or (empty_mu < 0)):
			dists[UCSCname] = meanEvolTimes(taudistrib_est_whl)

		# Compute mean evolutionary time for a species taking into account empty windows.
		# Apply correction to empty windows if a value for a high mutation rate is provided (empty_mu).
		else:
			# Evolutionary time estimate for an empty window.
			empty_tau       = empty_mu*my_dataset.divergenceTimes[UCSCname] 
			dists[UCSCname] = meanEvolTimes(taudistrib_est_whl, empty_windows[UCSCname], empty_tau)
	return dists

def getEmptyWindows(alpha, my_dataset):
	empty_dict = {}
	for UCSCname in my_dataset.speciesUCSCnames:
		logFilename = my_dataset.getLogFilename_estimateEvolTimes(UCSCname,alpha)
		# Parse log file of 4_estimateEvolTimes.py to get the number of discarded windows.
		winIgnor_cnt = 0
		winTotal_cnt = 0
		with open(logFilename, 'r') as file:
			for line in file:
				# "- Ignored windows: 11966 out of 116941 (10.23 %)." -> Information is per chromosome.
				m = re.match(r'.*Ignored windows: (\d+) out of (\d+) \(([\d\.]+) .\).*', line)
				if m:
					winIgnor_cnt += int(m.group(1))
					winTotal_cnt += int(m.group(2))
		empty_dict[UCSCname] = (winIgnor_cnt, winTotal_cnt)
	return empty_dict

###############################
# Auxiliary methods for plotting data.

def rgba_to_hex(rgba):
	# Convert each component from float [0.0, 1.0] to int [0, 255]
	r = int(rgba[0] * 255)
	g = int(rgba[1] * 255)
	b = int(rgba[2] * 255)
	a = int(rgba[3] * 255)  # Alpha is optional for hex representation

	# Format the hex string
	hex_color = '#{:02x}{:02x}{:02x}'.format(r, g, b)
	return hex_color
	
def makeColormap(my_dataset):
	# Create colormap.
	gradSpecies = Grad(["#1ca470","#fac784","#fac784","#fac784","#f77a7d","#f77a7d"])
	cmapSpecies = gradSpecies.to_cmap()
	# Choose colors for species.
	colormapSpecies={}
	min_t = min(list(my_dataset.divergenceTimes.values()))
	max_t = max(list(my_dataset.divergenceTimes.values()))
	for UCSCname in my_dataset.speciesUCSCnames:
		t_norm = (my_dataset.divergenceTimes[UCSCname]-min_t)/(max_t-min_t)
		color  = rgba_to_hex(cmapSpecies(t_norm))
		colormapSpecies[UCSCname] = color
	return colormapSpecies

###############################
# Plot methods.

def plotComparisonEvolTimes(paperRef, my_dataset, alldists_ours, alldists_other, outliers, sameAxis, hasEmptyWindows):

	colormapSpecies = makeColormap(my_dataset)

	# Prepare points to plot.
	yvals  = []
	yvals_err = []
	xvals  = []
	colors = []
	for UCSCname in my_dataset.speciesUCSCnames:
		# Check if there are outliers to be removed.
		if ((UCSCname in outliers) or (UCSCname not in alldists_ours) or (UCSCname not in alldists_other)): continue
		xvals.append(alldists_other[UCSCname])
		yvals.append(np.mean(alldists_ours[UCSCname]))
		yvals_err.append(np.std(alldists_ours[UCSCname]))
		colors.append(colormapSpecies[UCSCname])
	
	# Create plot.
	plotFilename = my_dataset.getOutFilename_plot_EvolTimesComp(paperRef, hasEmptyWindows)
	pp = PdfPages(plotFilename)
	fig, ax = plt.subplots(1, 1, figsize=(9, 7))

	# Labels for axis.
	author = (paperRef.split("_")[0])[:-4]
	year   = (paperRef.split("_")[0])[-4:]
	ylabel = "Evolutionary time (τ; α ➞ ∞)"
	xlabel = f"Substitution distance (tree from {author} et al. {year})"
	ax.set_ylabel(ylabel)
	ax.set_xlabel(xlabel)

	# Plot diagonal line.
	if(sameAxis):
		min_val   = 0.0
		max_val   = np.max(np.max(np.array([xvals,yvals,list(np.array(yvals)+np.array(yvals_err))])),0)
		extra_gap = (max_val-min_val)*0.05
		ax.set_xlim((min_val,max_val+extra_gap))
		ax.set_ylim((min_val,max_val+extra_gap))
		diag=np.linspace(min_val, max_val+extra_gap, 1000)
		ax.plot(diag,diag,linestyle='--',linewidth=5,color="darkgray")

	# Plot error bars.
	for (x,y,y_err,c) in zip(xvals, yvals, yvals_err, colors):
		ax.errorbar(x, y, yerr=[[y_err],[y_err]], fmt='o', markersize=12,capsize=10, color=c, 
					ecolor="gray", markeredgecolor="gray", elinewidth=2, markeredgewidth=2) #, alpha=0.5

	ax.xaxis.get_label().set_fontsize(18)
	ax.yaxis.get_label().set_fontsize(18)
	ax.set_facecolor("#F5F5F5") # "whitesmoke"
	ax.tick_params(axis='both', labelsize=15)
	for spine in ax.spines.values():
		spine.set_visible(False)

	print(f"Saving plot {plotFilename}...")
	pp.savefig(fig)
	pp.close()

def makeFigure3(alpha, my_dataset):

	# Mapping between UCSC names (our dataset) and species names from other datasets.
	mappingFilename = os.path.join(my_dataset.dirTrees, "mapspecies.csv")
	mappingSpecies  = parseMapping(mappingFilename)

	# Load evolutionary time estimates from other studies.
	print(f"Loading evolutionary time estimates from other studies...")
	alldists_Upham2019   = loadData(my_dataset, mappingSpecies, "Upham2019")
	alldists_Kuderna2023 = loadData(my_dataset, mappingSpecies, "Kuderna2023")

	# Load info on how many windows do not have an estimate in our study (not enough PCSs).
	# - Criteria to ignore windows: less than 5 distinct PCS sizes above 5 base pairs (more details on the log of 4_estimateEvolTimes.py).
	print(f"Computing empty windows...")
	empty_windows = getEmptyWindows(alpha, my_dataset)

	# Load evolutionary time estimates from our study taking into account amount of empty windows.
	print(f"Loading evolutionary time estimates from our study...")
	alldists_without_empty_ours = loadOurData(alpha, my_dataset)
	alldists_with_empty_ours    = loadOurData(alpha, my_dataset, empty_windows, 0.018)

	# Plot data.
	plotComparisonEvolTimes("Kuderna2023", my_dataset, alldists_without_empty_ours, alldists_Kuderna2023, ["mm39"],    True,  False)
	plotComparisonEvolTimes("Upham2019",   my_dataset, alldists_without_empty_ours, alldists_Upham2019,   ["ornAna2"], False, False)
	plotComparisonEvolTimes("Upham2019",   my_dataset, alldists_with_empty_ours,    alldists_Upham2019,   ["ornAna2"], True,  True)
	

####################################
# MAIN.
####################################
if (__name__ == '__main__'):
	alpha	   = 10.0 # Substitutions
	my_dataset = Dataset()

	# Save overall time.
	timeTrack = Time()
	timeTrack.start()

	makeFigure3(alpha, my_dataset)

	timeTrack.stop()
	timeTrack.print()

	print(f"Done!")

