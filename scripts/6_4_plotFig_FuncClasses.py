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
lowest evolutionary times, summing up to 308 MB (around 10% of
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

Running the script on a single core takes **16 minutes** and requires a small amount of memory. 
In total, the output files (svg+pdf plots) require 3.9 MB of disk space.

======================================  =======
Step                                    Time (s)
======================================  =======
Figure 5                                 248.59
Figure 5 - SI (conserv. thresh.: ~1%)    214.66
Figure 5 - SI (conserv. thresh.: ~5%)    231.99
Figure Supp S8                           243.96
**Total time**                           939.20
======================================  =======

**Output files**:

A SVG version of the files listed below is also saved in the same directory:

	1. ``funcClasses-conservDistrib.conservedThresh308824475.pdf``: The PDF output file containing the plots for Figure 5;
	2. ``funcClasses-conservDistrib.conservedThresh30882447.pdf``: The PDF output file for Supp. Figure with similar plots to Figure 5, but using a different conservation threshold (1%);
	3. ``funcClasses-conservDistrib.conservedThresh154412238.pdf``: The PDF output file for Supp. Figure with similar plots to Figure 5, but using a different conservation threshold (5%);
	4. ``funcClasses-stats.overlapThresh0.0.pdf``: The PDF output file for Supp. Figure S8;

Function details
----------------

Only relevant functions have been documented below. 
For more details on any function, check the comments in the souce code.

"""

import os
import pickle
import random
import sys
from collections import defaultdict
import re
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
from matplotlib.transforms import Bbox
from matplotlib import gridspec
import textwrap

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
	# windowSizesSorted  = sorted(list(windowSizes.keys()))
	# windowSizesDetails = ' '.join([f"{winSize} bps = {windowSizes[winSize]} windows; " for winSize in windowSizesSorted])
	# print(f"Window sizes details : {windowSizesDetails}")
	
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
			begPos, endPos    = windowId
			annotOverlapTotBp = 0
			for (annotType, annotOverlapPc) in windowAnnot.items():
				windowsPerAnnot[annotType].append(((qChrom,windowId),annotOverlapPc*(endPos-begPos)))
				annotOverlapTotBp += annotOverlapPc*(endPos-begPos)
			# Empty windows.
			if(annotOverlapTotBp == 0):
				windowsPerAnnot[my_dataset.unannotatedCategory].append(((qChrom,windowId),endPos-begPos))
	return windowsPerAnnot

###############################
# Plot methods.

def sortAnnotations(annotations):
	""" Sort annotations. Each item in the list is a different annotation.
		An annotation is a tuple of three numbers: (Main category: X, Sub-category: Y, Sub-sub-category: Z).
		These numbers specify the hierarchical order (X.Y.Z) used in the plots.
	"""
	return sorted(annotations, key=lambda order: (order[0],order[1],order[2]))

def get_aspect(ax):
	# Total figure size
	figW, figH = ax.get_figure().get_size_inches()
	# Axis size on figure
	_, _, w, h = ax.get_position().bounds
	# Ratio of display units
	disp_ratio = (figH * h) / (figW * w)
	# Ratio of data units
	# Negative over negative because of the order of subtraction
	data_ratio = sub(*ax.get_ylim()) / sub(*ax.get_xlim())

	return disp_ratio / data_ratio

def makeAnnotationMatrix(alpha, my_dataset, unit, threshold_annotation, threshold_conserved=-1.0):
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
	# M[species_i][annot_j] = number of windows annotated as ``annot_j`` with estimates in ``species_i``.
	colsLabels = sortAnnotations(list(windowsPerAnnot.keys()))
	rowsLabels = my_dataset.speciesUCSCnames
	M = np.zeros((len(rowsLabels), len(colsLabels)))

	# Cols headers.
	# Total windows annotated, regardless if they have evol. estimates or not.
	if (unit == "win"):
		colsBars = np.array([len(windowsPerAnnot[annotType]) for annotType in colsLabels])
	# Total bps annotated, regardless if they have evol. estimates or not.
	elif(unit == "bps"):
		colsBars = np.array([sum([annotOverlapBps for ((qChrom,windowId),annotOverlapBps) in windowsPerAnnot[annotType]]) for annotType in colsLabels])
	else:
		print(f"ERROR! Unit not recognized: {unit}. Valid units:\n\t'win': Stats will be based on the number of windows;\n\t'bps': Stats will be based on the number of base pairs.\nAborting execution.")
		sys.exit(1)

	# Compute matrix M.
	timeTrack.startStep(f"Create matrix")
	rowsBars = []
	for (rowIdx, UCSCname) in enumerate(rowsLabels):
		timeTrack.startStep(f"Create matrix: row {rowIdx+1} out of {len(rowsLabels)} [{UCSCname}]")

		# Load evolutionary estimates.
		evolEstPerWindows = loadEvolEstimates(my_dataset, UCSCname, alpha, threshold_conserved)
		
		# Rows headers.
		# Total windows with evol. estimates, regardless if they have annotation or not.
		if (unit == "win"):
			rowsBars.append(sum([len(windowsPerChr) for (qChrom, windowsPerChr) in evolEstPerWindows.items()]))
		# Total bps with evol. estimates, regardless if they have annotation or not.
		else: 
			rowsBars.append(sum([endPos-begPos for (qChrom, windowsPerChr) in evolEstPerWindows.items() for (begPos, endPos) in windowsPerChr]))

		# Cells M[i][j].
		for (colIdx, annotType) in enumerate(colsLabels):
			if(annotType not in windowsPerAnnot): continue
			# Row content: Count of all windows (or base pairs) with a given annotation that have estimates.
			M[rowIdx][colIdx] = sum([1 if (unit == "win") else annotOverlapBps for ((qChrom,windowId),annotOverlapBps) in windowsPerAnnot[annotType] if ((qChrom in evolEstPerWindows) and (windowId in evolEstPerWindows[qChrom]))])
		timeTrack.stopStep()
	rowsBars = np.array(rowsBars)
	timeTrack.stopStep()

	timeTrack.stop()
	timeTrack.print()
	return (M,rowsLabels,rowsBars,colsLabels,colsBars)

def makeFigureSuppS7(alpha, my_dataset):

	threshold_annotation = 0.0
	M,rowsLabels,rowsBars,colsLabels,colsBars = makeAnnotationMatrix(alpha, my_dataset, "win", threshold_annotation)

	# Compute % of windows annotated with estimates.
	M=100*M/colsBars

	# Plot heatmap of annotations with bars on top.
	xlabels = [my_dataset.annotationCategories[annotOrder] for annotOrder in colsLabels]
	xlabels = [f"{mainCat} - {subCat} {subsubCat}" for (mainCat, subCat, subsubCat) in xlabels]
	ylabels = [f"{my_dataset.commonNames[UCSCname]}" for UCSCname in rowsLabels]

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
	plotFilename = my_dataset.getOutFilename_plot_FuncClasses(threshold_annotation,"svg")
	plt.savefig(plotFilename, format="svg")
	plotFilename = my_dataset.getOutFilename_plot_FuncClasses(threshold_annotation,"pdf")
	pp = PdfPages(plotFilename)
	pp.savefig(fig,bbox_inches = 'tight')
	pp.close()

def getColors(annotations):
	colors = {(1,1,1) : "#2a7fff",
			(1,2,1) : "#3771c8",
			(1,3,1) : "#00ccff", 

			(1,2,2) : "#37abc8", 
			(1,2,3) : "#88cede", 
			(1,2,4) : "#37abc8", 
			(1,2,5) : "#88a9dd", 

			(2,2,1) : "#d4aa00",
			(2,2,2) : "#ffcc00",
			
			(2,1,1) : "#ff6600",
			(2,1,2) : "#deaa87",
			(2,1,3) : "#ff8838",
			(2,1,4) : "#ff9955",
			(2,1,5) : "#ffb380", 
			
			(3,3,1) : "#37c871",
			(3,3,2) : "#37c871",
			(3,3,3) : "#37c871", 
			(3,4,1) : "#37c871",

			(3,1,1) : "#447821",
			(3,1,2) : "#447821",
			(3,1,3) : "#447821",
			(3,1,4) : "#447821", 

			(3,5,1) : "#88aa00",
			(3,5,3) : "#88aa00",
			(3,5,4) : "#88aa00",
			(3,5,5) : "#88aa00",
			(3,5,2) : "#88aa00", 
			
			(3,2,2) : "#56e0aa",
			(3,8,1) : "#56e0aa",
			(3,6,1) : "#56e0aa",
			(3,7,1) : "#56e0aa",

			(3,9,1) : "#56e0aa",
			
			(10,1,1) : "#6f7c91"}
	return [colors[annot] if annot in colors else "#ff0000" for annot in annotations]

def plotFuncClassDistrib(my_dataset, ax, label, desc, color, x, y_totbps, y_selbps):
	""" Plot the distribution of a given annotation in conserved windows for different species.
	"""

	# Background.
	ax.patch.set_facecolor('none')
	# Borders.
	for spine in ax.spines.values():
		spine.set_color('white')
			
	# Y limits.
	y_base_avg = np.mean(y_totbps)
	y_diff	   = np.max(np.abs(y_selbps-y_base_avg))
	y_diff	   = y_diff + y_diff*0.15
	y_lim_lb   = y_base_avg-y_diff
	y_lim_ub   = y_base_avg+y_diff
	ax.set_ylim((y_lim_lb, y_lim_ub))
	
	# Remove x-axis.
	ax.set_xticklabels([])

	# y-ticks
	tick_positions    = [y_base_avg-y_diff, y_base_avg, y_base_avg+y_diff] 
	if (y_base_avg-y_diff < 0):
		if((y_base_avg-np.min(y_selbps)) > 0.75*y_diff):
			tick_positions    =	[np.min(y_selbps), y_base_avg, y_base_avg+y_diff]
		else:
			tick_positions    =	[y_base_avg, y_base_avg+y_diff]

	tick_labels       = [f"{yval*100:.1g}" if(yval*100 < 1) else f"{int(round(yval*100))}" for yval in tick_positions] #{yval*100:.4f}

	# In case the labels end up the same because of lack of precision: adds more precision to the middle label.
	if ((tick_labels[0] == tick_labels[1]) or ((len(tick_labels) == 3) and (tick_labels[1] == tick_labels[2]))): tick_labels[1] = f"{tick_positions[1]*100:.2g}"

	label_middle      = tick_labels[0] if (len(tick_labels) == 2) else tick_labels[1]
	tiny_font         = 8 
	small_font        = 10 # 11
	medium_font       = 12
	big_font          = 14 # 15
	label_middle_size = tiny_font if (len(label_middle) > 4) else big_font
	if (len(label_middle) == 4): label_middle_size = medium_font
	label_exts_size   = tiny_font if(sum([1 if (len(label) > 4) else 0 for label in tick_labels]) > 0) else small_font

	label_sizes       = [label_exts_size, label_middle_size, label_exts_size] if (len(tick_positions)==3) else [label_middle_size, label_exts_size]
	label_aligns      = ["bottom", "center", "top"] if (len(tick_positions)==3) else ["center", "top"]
	ax.set_yticks(tick_positions)
	ax.set_yticklabels(tick_labels)
	for tick, size, align in zip(ax.get_yticklabels(), label_sizes, label_aligns):
		tick.set_fontsize(size)
		tick.set_verticalalignment(align)
	ax.tick_params(axis="y",direction="in", pad=-6)
	ax.set_axisbelow(False)

	# Hide tick marks.
	ax.tick_params(axis='both', which='both', length=0)

	# Zebra background.
	highlighted_species  = ["panTro6","mm39","myoLuc2","danRer11"]
	highlighted_color    = "#d3cccc"
	highlighted_opacity  = 1.00
	zebra_colors         = ['#fdfdfd','#f7f4f4']
	zebra_opacity        = [1.0,1.0]
	zebra_barsize        = 0.98
	for idx, UCSCname in enumerate(my_dataset.speciesUCSCnames):
		bg_color = highlighted_color if (UCSCname in highlighted_species) else zebra_colors[int(idx % 2 == 0)]
		bg_alpha = highlighted_opacity if (UCSCname in highlighted_species) else zebra_opacity[int(idx % 2 == 0)]
		ax.axvspan(idx-zebra_barsize/2, idx+zebra_barsize/2, color=bg_color, alpha=bg_alpha, lw=0)

	# Ratio to check if a class is under-represented or over-represented in slow evolving.
	under_color_hatch = "#f08080"
	under_color_bg    = "#f07d7d"
	over_color_hatch  = "#66cdaa"
	over_color_bg     = "#64cdaa"
	ax.fill_between(x, [min(yval,y_totbps) for yval in y_selbps], y_totbps, color=under_color_bg, alpha=0.2)
	ax.fill_between(x, [max(yval,y_totbps) for yval in y_selbps], y_totbps, color=over_color_bg,  alpha=0.2)

	ax.fill_between(x, [min(yval,y_totbps) for yval in y_selbps], y_totbps, edgecolor=under_color_hatch, facecolor="none", hatch='////', alpha=1.0)
	ax.fill_between(x, [max(yval,y_totbps) for yval in y_selbps], y_totbps, edgecolor=over_color_hatch,  facecolor="none", hatch='////', alpha=1.0)

	# Data.
	ax.hlines(y=y_totbps, xmin=min(x), xmax=max(x), linewidth=1, color="black")  # Baseline  (% annotation in the whole-genome)
	ax.scatter(x, y_selbps, color=color) # Conserved (% annotation in the most conserved windows of a given pairwise species (human and another vertebrate))

	# title
	#ax.set_title(f"{desc}",fontsize=5)
	desc = "\n".join(textwrap.wrap(desc, 15, break_long_words=True))
	ax.text(-0.25, 0.5, f"{desc}", rotation=90, fontsize=8, verticalalignment='center', horizontalalignment='center', transform=ax.transAxes)

	# set aspect ratio
	ax_aspect = 0.77/2.33 # Height: 0.77 inches; Width: 2.33 inches.
	ax.set_box_aspect(ax_aspect)

	# Information used in the main text.
	labelStr = " ".join(my_dataset.annotationCategories[label])
	print(f"Annotation {labelStr}\n\t - Percentage of the number of bps in conserved regions (avg. for all species) : {np.mean(y_selbps)}\n\t - Percentage of the number of bps in the whole genome (avg. for all species) : {np.mean(y_totbps)}\n\t - Ratio conserved/whole-genome : {np.mean(y_selbps)/np.mean(y_totbps)}")

def getAnnotOrder(annotGrid_paper, rowsLabels):
	rowsLabels_sel    = [annotSel for (annotSel, labelSel) in annotGrid_paper]
	rowsLabels_notsel = [annot for annot in rowsLabels if annot not in rowsLabels_sel]
	rowsOrder         = rowsLabels_sel + rowsLabels_notsel
	rowsOrder         = [rowsLabels.index(annot) for annot in rowsOrder if annot in rowsLabels]
	return rowsOrder

def getAnnotTitles(annotGrid_paper,rowsLabels,my_dataset):
	titles_dict = {annotSel : labelSel for (annotSel, labelSel) in annotGrid_paper}
	for annot in rowsLabels:
		if (annot not in titles_dict): titles_dict[annot] = " ".join(my_dataset.annotationCategories[annot])
	return titles_dict

def makeFigure5(alpha, my_dataset, threshold_conserved=250e6):

	# TODO: Include un-annotated.

	# A functional class must cover over 90% of the window for
	# the window to be annotated with that class.
	threshold_annotation =  my_dataset.thresholdAnnot
	threshold_annotation = 0.0

	# Sort windows by their evolutionary time estimates and 
	# take the slowest windows that sum up to 250 MB (~8% of the genome).
	threshold_conserved  = threshold_conserved

	# Compute matrix of annotations taking both thresholds into account.
	# Flip cols and rows.
	M,colsLabels,colsBars,rowsLabels,rowsBars = makeAnnotationMatrix(alpha, my_dataset, "bps", threshold_annotation, threshold_conserved=threshold_conserved)
	M = M.T
	rowColors = getColors(rowsLabels)

	totalBps  = sum([chromSize for (chrom, chromSize) in my_dataset.chromSizes.items()])
	selecBps  = colsBars.T # threshold_conserved

	# Normalize values: 
	# - Total values are normalized by the total number of base pairs in the human genome; 
	# - Selected values are normalized by the total number of base pairs in 8% of the genome.
	rowsBpsTotal = rowsBars / totalBps
	rowsBpsCons  = M / selecBps
	
	# Plots info.
	y = np.arange(len(rowsLabels))
	x = np.arange(len(colsLabels))

	fig = plt.figure(figsize=(3,40))
	gs  = gridspec.GridSpec(len(rowsLabels), 1, wspace=0.0, hspace=0.08)

	# fig = plt.figure(figsize=(3,62)) # Ideal size: (3,62) (other sizes will distort scale)
	# gs  = gridspec.GridSpec(len(rowsLabels), 1, wspace=0.0, hspace=0.90)

	ymin   = 0.1
	ymax   = 2

	annotGrid_paper   = [((1,1,1),"Proteins"), ((1,2,1), "LncRNAs"), ((1,2,2), "micro RNAs"), ((1,2,4), "snoRNAs"), ((1,2,3),"snRNAs"), ((1,2,5),"Other RNAs"), ((1,3,1),"Pseudogenes"), ((3,1,1),"LINE"), ((3,1,2),"LINE (CR1)"), ((3,1,3),"LINE (L1)"), ((3,1,4),"LINE (L2)"),((3,5,1),"LTR"),((3,5,5),"LTR (ERVL)"),((3,5,4),"LTR (MaLR)"),((3,5,2),"LTR (ERV1)"),((3,5,3),"LTR (ERVK)"),((10,1,1),"Un-annotated"),((2,2,1),"Introns"),((2,2,2),"UTRs"),((2,1,1),"Enhancers"),((2,1,3),"DNAase HSs"),((2,1,5),"Locus ctrl."),((2,1,4),"Insulators"),((2,1,2),"Silencers"),((3,3,1),"DNA"),((3,3,3),"DNA (hAT)"),((3,3,2),"DNA (TcMar)"),((3,4,1),"RNA"),((3,2,2),"SINE (ALU)"),((3,6,1),"Satellite"),((3,7,1),"Micro-satellite"),((3,8,1),"Low comp."),((3,9,1),"Others")]

	rowsIdxOrder = getAnnotOrder(annotGrid_paper, rowsLabels)
	rowsTitles   = getAnnotTitles(annotGrid_paper,rowsLabels,my_dataset)
	
	for idx, rowIdx in enumerate(rowsIdxOrder):

		label    = rowsLabels[rowIdx]
		color    = rowColors[rowIdx]
		desc     = rowsTitles[label]
		y_totbps = rowsBpsTotal[rowIdx] # whole-genome bps (Normalized by the total number of base pairs).
		y_selbps = rowsBpsCons[rowIdx]  # conserved bps (Normalized by the total number of base pairs in conserved regions).

		ax = plt.subplot(gs[idx])
		plotFuncClassDistrib(my_dataset, ax, label, desc, color, x, y_totbps, y_selbps)

	# Save plot.
	plotFilename = my_dataset.getOutFilename_plot_FuncClassesDistribInConserv(threshold_conserved,"svg")
	plt.savefig(plotFilename, format="svg")
	plotFilename = my_dataset.getOutFilename_plot_FuncClassesDistribInConserv(threshold_conserved,"pdf")
	pp = PdfPages(plotFilename)
	pp.savefig(fig,bbox_inches = 'tight')
	pp.close()

####################################
# MAIN.
####################################
if (__name__ == '__main__'):
	alpha	     = 1.1
	my_dataset = Dataset()

	# Save overall time.
	timeTrack = Time()
	timeTrack.start()

	timeTrack.startStep("Figure 5")
	threshold_conserved_10pc = 308824475 # 308.82 MB (~10%)
	makeFigure5(alpha, my_dataset, threshold_conserved_10pc)       # Figure 5.
	timeTrack.stopStep()

	timeTrack.startStep("Figure 5 - SI (conserv. thresh.: ~1%)")
	threshold_conserved_1pc = 30882447 # 30.88 MB (~1%)
	makeFigure5(alpha, my_dataset, threshold_conserved_1pc)         # Supp Fig: Same as Figure 5, but a different threshold for conservation (1%).
	timeTrack.stopStep()

	timeTrack.startStep("Figure 5 - SI (conserv. thresh.: ~5%)")
	threshold_conserved_5pc = 154412238 # 154.41 MB (~5%)
	makeFigure5(alpha, my_dataset, threshold_conserved_5pc)         # Supp Fig: Same as Figure 5, but a different threshold for conservation (5%).
	timeTrack.stopStep()

	timeTrack.startStep("Figure Supp S8")
	makeFigureSuppS7(alpha, my_dataset)                             # Supp Fig: Functional classes stats.
	timeTrack.stopStep()

	timeTrack.stop()
	timeTrack.print()

	print(f"Done!")
