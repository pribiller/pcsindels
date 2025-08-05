"""Plot ``Figure 5`` from the paper and related plots 
from Supplementary Figures, i.e., plots related to the 
distribution of slow-evolving regions across functional classes.

All annotations for human (hg38) were obtained from 
**UCSC Table Browser**: 
  * Coding sequence coordinates were taken from the tracks **GENCODE-V44** and **CCDS**; 
  * Repeats from the track **RepeatMasker**; 
  * Regulation elements from the **RefSeq** track.

In the generated plots, rows indicate different functional classes, 
columns indicate different pairwise alignments between human and another
vertebrate (left: closest species to human; right: farthest species
from human). For each pairwise alignment, the set of 
*slow-evolving* windows was defined by selecting the windows with the
lowest evolutionary times, summing up to 250 MB (around 8% of
the human genome). Variations on this threshold were tested
and generated as Supplementary Figures.

For each functional class, the y-axis shows the percentage 
of windows identified as *slow-evolving* that fall
under that class, compared to the baseline (solid black line in the
middle of each graph), which denotes the percentage of windows
in the entire genome under the same functional class.

- **Use**::
	
	python3 6_4_plotFig_FuncClasses.py

- **Example of Usage**::

	python3 ~/code/6_4_plotFig_FuncClasses.py

- **Input Parameter**:

To ensure the graphs match those used in the paper,
the parameters are hard-coded in the script and 
cannot be modified via command line.

Pre-requisites
--------------

Before using this script, make sure all the required files were pre-computed:

a) Files with sampled evolutionary times
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Make sure to run ``5_sampleEvolTimes.py`` for **α=1.1**.

b) Annotated windows 
^^^^^^^^^^^^^^^^^^^^^
Make sure to run ``annotateWindows.py``, to annotate all windows used in the analysis.

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
from matplotlib.colors import LogNorm
from matplotlib.transforms import Bbox
from colorir import Grad
import seaborn as sns

from utils.dataset import Dataset
from utils.basicTypes import Time


####################################################################
# Annotation related methods.

def loadAnnotations(my_dataset, threshold_annotation):

	# Load protein panther classification.
	# It complements the gene annotation by providing a more specific identification of the protein family.
	# WARNING: It ended up not being actually used in the main text (maybe in a future paper?).
	# pantherFilename = my_dataset.getOutFilename_annotationsPanther()
	# proteinClasses  = loadProteinPantherAnnotation(pantherFilename)

	# Load annotations for windows.
	annotFilename = my_dataset.getOutFilename_annotatedWindows()
	print(f"Loading file {annotFilename}")
	if (not os.path.isfile(annotFilename)):
		print(f"ERROR! Annotation file not found: {annotFilename}")
		sys.exit(1)
	annotationsFound, annotatedWindows = pickle.load(open(annotFilename, 'rb'))

	print(f"Computing chromosome sizes")
	totalBps = 0
	totalWin = 0
	windowSizes = defaultdict(int)
	for qChrom in annotatedWindows.keys():
		chromSize = 0
		for ((begPos,endPos), windowAnnot) in annotatedWindows[qChrom].items():
			winSize   = endPos-begPos
			chromSize += winSize
			totalWin  += 1
			windowSizes[winSize] += 1
		print(f"{qChrom} : {chromSize} bps")
		totalBps += chromSize
	print(f"- Total Bps={totalBps} / Total Win={totalWin}")
	windowSizesSorted = sorted(list(windowSizes.keys()))
	for winSize in windowSizesSorted:
		print(f"{winSize} : {windowSizes[winSize]} windows.")


	# Filter annotations based on the % of overlap with a window (threshold_annotation).
	if(threshold_annotation > 0.0):
		annotatedWindowsAboveThresh = {}
		for qChrom in my_dataset.chromLst:
			if (qChrom not in annotatedWindows):
				print(f"WARNING! {qChrom} not annotated. Skipping it.")
				continue
			annotatedWindowsAboveThresh[qChrom] = {windowId : {annotType : annotOverlap for (annotType, annotOverlap) in windowAnnot.items() if annotOverlap > threshold_annotation} for (windowId, windowAnnot) in annotatedWindows[qChrom].items()}
		annotatedWindows = annotatedWindowsAboveThresh
	return (annotationsFound, annotatedWindows)

###############################
# Methods to load info from our study.

def loadEvolEstimates(my_dataset, UCSCname, alpha, threshold_conserved=-1.0):
	# Load sampled distributions.
	sampEvolTimesFilepath = my_dataset.getOutFilename_sampleEvolTimes(UCSCname, alpha)
	if(not os.path.isfile(sampEvolTimesFilepath)):
		print(f"[{UCSCname}] WARNING! File not found: {sampEvolTimesFilepath}. Aborting computation!")
		sys.exit(1)
	PCSdistrib_obs_whl, PCSdistrib_est_whl, taudistrib_est_whl, PCSdistrib_obs_chr, PCSdistrib_est_chr, taudistrib_est_chr, taudistrib_est_det_chr = pickle.load(open(sampEvolTimesFilepath, 'rb'))
	# Filter most conserved windows.
	if(threshold_conserved > 0):
		# Gather all windows with estimates and sort them by their evolutionary times (from the most conserved ones to the fast-evolving ones).
		taudistrib_est_det_chr_sorted = sorted([(qChrom, windowId, taudistrib_est_det_chr[qChrom][windowId]) 
								for qChrom in taudistrib_est_det_chr for windowId in taudistrib_est_det_chr[qChrom]], key=lambda winInfo: winInfo[2])
		size_so_far = 0
		taudistrib_est_det_chr = {qChrom : {} for qChrom in taudistrib_est_det_chr}
		for (qChrom, windowId, evolEstAvg) in taudistrib_est_det_chr_sorted:
			begPos, endPos = windowId
			size_so_far   += (endPos-begPos)
			taudistrib_est_det_chr[qChrom][windowId] = evolEstAvg
			if(size_so_far > threshold_conserved): break

	return taudistrib_est_det_chr

###############################
# Auxiliary methods for plotting data.

def groupWindowsByAnnotType(my_dataset, annotatedWindows):
	# Group windows per annotation type.
	windowsPerAnnot = defaultdict(list)
	for qChrom in annotatedWindows:
		for (windowId, windowAnnot) in annotatedWindows[qChrom].items():
			if (len(windowAnnot) == 0):
				begPos, endPos = windowId
				windowsPerAnnot[my_dataset.unannotatedCategory].append(((qChrom,windowId),endPos-begPos))
			else:
				for (annotType, annotOverlapBps) in windowAnnot.items():
					windowsPerAnnot[annotType].append(((qChrom,windowId),annotOverlapBps))
	return windowsPerAnnot

###############################
# Plot methods.

def sortAnnotations(annotations):
	""" Sort annotations. Each item in the list is a different annotation.
		An annotation is a tuple of three numbers: (Main category: X, Sub-category: Y, Sub-sub-category: Z).
		These numbers specify the hierarchical order (X.Y.Z) used in the plots.
	"""
	return sorted(annotations, key=lambda order: (order[0],order[1],order[2]))

def makeAnnotationMatrix(alpha, my_dataset, matrix_type, threshold_annotation, threshold_conserved=-1.0):
	timeTrack = Time()
	timeTrack.start()

	# A functional class can have an overlap of any size with a window.
	timeTrack.startStep("Load annotations")
	annotations, annotatedWindows = loadAnnotations(my_dataset, threshold_annotation)	
	timeTrack.stopStep()
	
	# Group windows by annotation type.
	timeTrack.startStep("Group annotations by type")
	windowsPerAnnot = groupWindowsByAnnotType(my_dataset, annotatedWindows)
	timeTrack.stopStep()

	# Create a matrix where rows = species and cols = annotations.
	# M[species_i][annot_j] = % of windows annotated as ``annot_j`` with estimates in ``species_i``.
	colsLabels = sortAnnotations(list(annotations))
	rowsLabels = my_dataset.speciesUCSCnames
	M = np.zeros((len(rowsLabels), len(colsLabels)))
	if (matrix_type in ["all","annot"]): M = np.zeros((len(rowsLabels), len(colsLabels)), dtype=int)

	# Total windows annotated, regardless if they have evol. estimates or not.
	colsBars = np.array([len(windowsPerAnnot[annotType]) for annotType in colsLabels])
	print(f"Sum={sum(colsBars)} / List={colsBars}")

	# Make M.
	timeTrack.startStep(f"Create matrix")
	for (rowIdx, UCSCname) in enumerate(rowsLabels):
		# timeTrack.startStep(f"Create matrix: row {rowIdx+1} out of {len(rowsLabels)} [{UCSCname}]")
		evolEstPerWindows = loadEvolEstimates(my_dataset, UCSCname, alpha, threshold_conserved)
		for (colIdx, annotType) in enumerate(colsLabels):
			if(annotType not in windowsPerAnnot): continue
			# Row content: count of all windows with estimates.
			if (matrix_type == "all"):
				M[rowIdx][colIdx] = len([1 for ((qChrom,windowId),annotOverlapBps) in windowsPerAnnot[annotType]])
			elif (matrix_type == "annot"):
				M[rowIdx][colIdx] = len([1 for ((qChrom,windowId),annotOverlapBps) in windowsPerAnnot[annotType] if ((qChrom in evolEstPerWindows) and (windowId in evolEstPerWindows[qChrom]))])
			else:
				M[rowIdx][colIdx] = len([1 for ((qChrom,windowId),annotOverlapBps) in windowsPerAnnot[annotType] if ((qChrom in evolEstPerWindows) and (windowId in evolEstPerWindows[qChrom]))])
				M[rowIdx][colIdx] = 100*M[rowIdx][colIdx]/colsBars[colIdx]
		# timeTrack.stopStep()
	timeTrack.stopStep()
	
	timeTrack.stop()
	timeTrack.print()
	return (M,rowsLabels,colsLabels,colsBars)

def makeFigureSuppS7(alpha, my_dataset, matrix_type):

	threshold_annotation = 0.0
	M,rowsLabels,colsLabels,colsBars = makeAnnotationMatrix(alpha, my_dataset,matrix_type,threshold_annotation)

	# Plot heatmap of annotations with bars on top.
	xlabels = [my_dataset.annotationCategories[annotOrder] for annotOrder in colsLabels]
	xlabels = [f"{mainCat} - {subCat} {subsubCat}" for (mainCat, subCat, subsubCat) in xlabels]

	ylabels = [f"{my_dataset.commonNames[UCSCname]}" for UCSCname in rowsLabels]

	# Options used for debugging.
	if(matrix_type == "all") or (matrix_type == "annot"):
		sns.set(font_scale=1)
		fig, ax = plt.subplots(figsize=(32,28))
		titleStr = "All annotated windows --- with and without estimates associated with" if(matrix_type == "all") else "Annotated windows only with estimates associated with"
		plt.title(titleStr)
		ax = sns.heatmap(M, ax=ax, annot=True, cmap="Blues", mask=(M==0), fmt="d", norm=LogNorm(), 
						 xticklabels=xlabels, yticklabels=ylabels) # fmt="d", mask=(matrix==0) | up_triang
		# Annotate the heatmap with median values on top
		max_all  = np.max(M,axis=0)
		for i, max_val in enumerate(max_all):
			ax.text(i + 0.5, -0.5, f"{max_val}", ha='center', va='center', color='black', fontsize=10, fontweight='bold')
		# Save plot.
		plotFilename = my_dataset.getOutFilename_plot_FuncClasses(threshold_annotation, matrix_type)
		pp = PdfPages(plotFilename)
		pp.savefig(fig,bbox_inches = 'tight')
		pp.close()
		return True

	sns.set(font_scale=2,rc={'axes.facecolor':'#ececec', 'figure.facecolor':(0,0,0,0)}) # 'figure.figsize':(18,24),
	fig, (ax_top, ax) = plt.subplots(2, 1, figsize=(32,28), gridspec_kw={'height_ratios': [1, 5]})
	plt.rcParams["axes.grid"] = False
	plt.tick_params(axis='both', which='major', labelsize=10, labelbottom = False, bottom=False, top = False, labeltop=True)	

	# Barplot
	ytopvals = colsBars
	xtopvals = list(range(ytopvals.size))
	ax_top.bar(xtopvals,ytopvals,width=1.0) #,width=1.0
	ax_top.margins(x=0, tight=True)
	ax_top.tick_params(
		axis='x',		 # changes apply to the x-axis
		which='both',	 # both major and minor ticks are affected
		bottom=False,	 # ticks along the bottom edge are off
		top=True,		 # ticks along the top edge are off
		labelbottom=False,
		labeltop=True)
	ax_top.set_yscale('log')
	ax_top.set_ylabel('Frequency (log)')
	ax_top.set_xticks(xtopvals)
	ax_top.set_xticklabels(xlabels, fontsize=20, rotation=90)

	# Heatmap
	ax = sns.heatmap(M, ax=ax, annot=True, cmap="Blues", mask=(M==0), vmin=0, vmax=100, fmt='', yticklabels=ylabels) # fmt="d", mask=(matrix==0) | up_triang, norm=LogNorm(), 

	# Apply the custom formatting function to the annotations
	for t in ax.texts:
		val = float(t.get_text())
		t.set_text('{:.2g}'.format(val) if (val < 100) else f"{int(val)}")
		t.set_fontsize(10 if (val < 1) else 20)

	# Annotate the heatmap with median values on top
	max_all  = np.max(M,axis=0)
	for i, max_val in enumerate(max_all): ax.text(i + 0.5, -0.5, f"{max_val}", ha='center', va='center', color='black', fontsize=10, fontweight='bold')

	# Set the font size of x and y tick labels
	ax.set_yticklabels(ax.get_yticklabels(), fontsize=20)

	# Adjust top bar plot size. (code needs to be added at the end, just before plt.show())
	(x0m, y0m), (x1m, y1m) = ax.get_position().get_points()  # main heatmap
	(x0h, y0h), (x1h, y1h) = ax_top.get_position().get_points()  # horizontal histogram
	ax_top.set_position(Bbox([[x0m, y0h], [x1m, y1h]]))

	# Save plot.
	plotFilename = my_dataset.getOutFilename_plot_FuncClasses(threshold_annotation, matrix_type)
	pp = PdfPages(plotFilename)
	pp.savefig(fig,bbox_inches = 'tight')
	pp.close()

def makeFigure5(alpha, my_dataset, threshold_conserved=0.08):

	# A functional class must cover over 90% of the window for
	# the window to be annotated with that class.
	annotations, annotatedWindows = loadAnnotations(my_dataset, my_dataset.thresholdAnnot)
	
	# Sort windows by their evolutionary time estimates and 
	# take the slowest windows that sum up to 250 MB (~8% of the genome).
	threshold_conserved_Rands2014 = 0.08

####################################
# MAIN.
####################################
if (__name__ == '__main__'):
	alpha	     = 1.1
	my_dataset = Dataset()

	# Save overall time.
	timeTrack = Time()
	timeTrack.start()

	timeTrack.startStep("Figure Supp S7")
	makeFigureSuppS7(alpha, my_dataset, "annot_pc")   # Supp Fig: Functional classes stats.
	timeTrack.stopStep()

	timeTrack.startStep("Figure Supp S7---Extra")
	makeFigureSuppS7(alpha, my_dataset, "all")
	makeFigureSuppS7(alpha, my_dataset, "annot")
	timeTrack.stopStep()

	# makeFigure5(alpha, my_dataset)        # Figure 5.
	# makeFigure5(alpha, my_dataset)        # Supp Fig: Same as Figure 5, but a different threshold for conservation (1%).
	# makeFigure5(alpha, my_dataset)        # Supp Fig: Same as Figure 5, but a different threshold for conservation (5%).
	# makeFigure5(alpha, my_dataset)        # Supp Fig: Same as Figure 5, but a different threshold for conservation (10%).

	timeTrack.stop()
	timeTrack.print()

	print(f"Done!")
