"""Plot ``Figure 4`` from the paper, i.e., a comparison
between the new estimated indel rates in the lineage relationships 
between human and 40 other vertebrates. It also includes, 
as a reference, **Direct estimates** and 
**Indirect estimates** from previous studies (listed below).

Direct estimates
-----------------

**Direct estimates** refers to methods that count mutations between generations 
in present-day individuals. The following studies were used as reference:

=================== ===== ============   =================================  ===========================  ================ ==============================  =========================================  ================================================
Authors             Year  Indel sizes    Indel rate estimation (original)   Indel rate estimation (CI)   Generation Time  Generation Time Interval (CI)   Indel rate PPPY (Per Position Per Year)    Indel rate PPPY (Per Position Per Year) (CI) 
=================== ===== ============   =================================  ===========================  ================ ==============================  =========================================  ================================================
Kloosterman et al.  2015  1-20           0.68*(10**(-9))                    -                            29.27            (24.385, 34.155)                2.3231886903593173e-11                     (1.990918179344908e-11, 2.788584149604477e-11)
Besenbacher et al.  2016  1-35           0.929*(10**(-9))                   -                            30.26            -                               3.07*(10**(-11))                           (2.91*10**(-11), 3.25*(10**(-11)))
Maretty et al.      2017  1-10           1.3*(10**(-9))                     -                            27.7             -                               4.70e-11                                   -
Besenbacher et al.  2015  1-50           1.5e-9                             (1.2e-9, 1.9e-9)             28.4             -                               5.28169014084507e-11                       (4.225352112676056e-11, 6.690140845070423e-11)
Kondrashov (del.)   2002  1-             0.526*(10**(-9))                   (0.216e-9,0.836e-9)          20               -                               2.63e-11                                   (1.58e-11, 4.18e-11)
Kondrashov (ins.)   2002  1-             0.182*(10**(-9))                   (0.072e-9,0.292e-9)          20               -                               0.91e-11                                   (0.36e-11, 1.46e-11)
Palamara et al.     2015  1-20           1.26*(10**(-9))                    (1.2e-09, 1.32e-09)          29               -                               4.3448275862068967e-11                     (4.137931034482759e-11, 4.5517241379310344e-11)
=================== ===== ============   =================================  ===========================  ================ ==============================  =========================================  ================================================

.. note::
	Note that the indel sizes vary in each study (see values in the ``Indel sizes'' column). As lower the maximum indel size is, higher the mutation rate estimation.

Indirect estimates
-------------------

**Indirect estimates** refers to estimates based on the evolutionary distance separating
two species divided by (twice) their divergence time. The following studies were used as reference:

=================== ===== ============ ============   =================================  ===========================  ================ ==============================  =========================================  ================================================
Authors             Year  Species      Indel sizes    Indel rate estimation (original)   Indel rate estimation (CI)   Generation Time  Generation Time Interval (CI)   Indel rate PPPY (Per Position Per Year)    Indel rate PPPY (Per Position Per Year) (CI) 
=================== ===== ============ ============   =================================  ===========================  ================ ==============================  =========================================  ================================================
Nachman and Crowell 2000  Chimp        1-4            2.3*(10**(-9))                     -                            20               -                               4.95049504950495e-11                       (3.712871287128713e-11, 6.188118811881188e-11)
Lunter              2007  Mouse        1-             0.053                              -                            2*87*(10**6)     (2*81.3*(10**6), 2*91*(10**6))  30.46e-11                                  (29.12087912087912e-11, 32.595325953259533e-11)
=================== ===== ============ ============   =================================  ===========================  ================ ==============================  =========================================  ================================================


Plots
------

In the generated plot, the new estimates are shown in orange, 
with the rectangle borders representing the standard deviation, 
and the middle point showing the mean indel rate. 
Indirect and direct estimates from previous studies are
indicated in green and blue, respectively. All indel rates were
adjusted to "per position per year" (PPPY) in order to make
them comparable. The orange dashed line indicates the average
indel rate across all species. If evolution were uniform across all
lineages, values should be concentrated around this point.

- **Use**::
	
	python3 6_3_plotFig_MutRatesComp.py

- **Example of Usage**::

	python3 ~/code/6_3_plotFig_MutRatesComp.py

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

b) Logs from evolutionary time estimates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Make sure to keep the logs from ``4_estimateEvolTimes.py`` for **α=1.1**. It contains the information regarding windows without estimates.

Time, Memory & Disk space
-------------------------

Running the script on a single core takes **3 minutes** (143.07 seconds) and requires a small amount of memory. 
In total, the output file requires 2.4 MB of disk space.

**Output files**:

The file ``mutRates-comp.highEvolTimeQuant0.99.svg'' contains the plot for Figure 4.

Function details
----------------

Only relevant functions have been documented below. 
For more details on any function, check the comments in the souce code.

"""

import os
import pickle
import random
import sys
from collections import namedtuple
import re
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle
from matplotlib.offsetbox import AnnotationBbox
from matplotlib.backends.backend_pdf import PdfPages

from utils.dataset import Dataset
from utils.basicTypes import Time

import textwrap
import cairosvg
import skunk

# Create a namedtuple to save information about previous studies on mutation rates.
MutRateStudy = namedtuple("MutRateStudy", "author publYear methodType ucscName mutRatePPPY mutRatePPPY_lb, mutRatePPPY_ub")

###############################
# Methods to load info from our study.

def meanEvolTimes(taudistrib_est_onespecies):
	mean_vals = []
	std_vals  = []
	# For each sampled collection of evolutionary times (whole-genome level, every 
	# window with "enough" information has a sampled evolutionary time associated with).
	for taudistrib_est in taudistrib_est_onespecies:
		# Gather taus from all windows with estimates.
		all_taus  = [tauval for (tauval, taucnt) in taudistrib_est.items() for idx in range(int(taucnt))]
		mean_vals.append(np.mean(all_taus))
		std_vals.append(np.std(all_taus))
	return (np.mean(mean_vals), np.mean(std_vals))

def computeDistribQuantile(taudistrib_est_onespecies, quantile):
	tau_quantile = []
	for taudistrib_est in taudistrib_est_onespecies:
		# Gather taus from all windows with estimates.
		all_taus  = [tauval for (tauval, taucnt) in taudistrib_est.items() for idx in range(int(taucnt))]	
		tau_quantile.append(np.quantile(all_taus,quantile))
	return np.mean(tau_quantile)

def loadOurData(alpha, my_dataset):
	dists = {}
	for UCSCname in my_dataset.speciesUCSCnames:
		# Load sampled distributions.
		sampEvolTimesFilepath = my_dataset.getOutFilename_sampleEvolTimes(UCSCname, alpha)
		if(not os.path.isfile(sampEvolTimesFilepath)):
			print(f"[{UCSCname}] WARNING! File not found: {sampEvolTimesFilepath}. Skipping computation!")
			continue
		PCSdistrib_obs_whl, PCSdistrib_est_whl, taudistrib_est_whl, PCSdistrib_obs_chr, PCSdistrib_est_chr, taudistrib_est_chr, taudistrib_est_det_chr = pickle.load(open(sampEvolTimesFilepath, 'rb'))
		# Compute mean evolutionary time for a species.
		dists[UCSCname] = meanEvolTimes(taudistrib_est_whl)
		print(f" {UCSCname} : mean indel rate={dists[UCSCname]}")
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

def computeEvolTimeEmpty(my_dataset,UCSCname, alpha, empty_evoltime_quantile):
	# Load sampled distributions.
	sampEvolTimesFilepath = my_dataset.getOutFilename_sampleEvolTimes(UCSCname, alpha)
	if(not os.path.isfile(sampEvolTimesFilepath)):
		print(f"[{UCSCname}] WARNING! File not found: {sampEvolTimesFilepath}. Skipping computation!")
		sys.exit()
	PCSdistrib_obs_whl, PCSdistrib_est_whl, taudistrib_est_whl, PCSdistrib_obs_chr, PCSdistrib_est_chr, taudistrib_est_chr, taudistrib_est_det_chr = pickle.load(open(sampEvolTimesFilepath, 'rb'))
	# Evolutionary time estimate for an empty window, 
	# based on the evolutionary time of a high tail of the distribution of evolutionary times.
	return computeDistribQuantile(taudistrib_est_whl, empty_evoltime_quantile)

def getNewEstimates(my_dataset, alpha, empty_evoltime_quantile=-1):
	indel_estimates_new  = []

	# Load info on how many windows do not have an estimate in our study (not enough PCSs).
	# - Criteria to ignore windows: less than 5 distinct PCS sizes above 5 base pairs (more details on the log of 4_estimateEvolTimes.py).
	print(f"Computing empty windows...")
	empty_windows = getEmptyWindows(alpha, my_dataset)

	# Load evolutionary time estimates from our study taking into account amount of empty windows.
	print(f"Loading evolutionary time estimates from our study...")
	alldists = loadOurData(alpha, my_dataset)

	print(f"Computing new indel rates...")
	for UCSCname in my_dataset.speciesUCSCnames:
		# Convert divergence time from millions of years to years.
		t = my_dataset.divergenceTimes[UCSCname]*(10**6)*2
				
		# Estimate indel mutation rate.
		evoltime_avg, evoltime_std = alldists[UCSCname]

		# Evolutionary time estimate for an empty window.
		winIgnor_cnt, winTotal_cnt = empty_windows[UCSCname]
		empty_pc        = winIgnor_cnt/winTotal_cnt
		empty_evoltime  = computeEvolTimeEmpty(my_dataset, UCSCname, alpha, empty_evoltime_quantile)

		mutRate_avg     = ((1.0-empty_pc)*evoltime_avg/t) + (empty_pc*(empty_evoltime/t))
		mutRate_std_lb  = ((1.0-empty_pc)*(evoltime_avg-evoltime_std)/t) + (empty_pc*(empty_evoltime/t))
		mutRate_std_ub  = ((1.0-empty_pc)*(evoltime_avg+evoltime_std)/t) + (empty_pc*(empty_evoltime/t))

		indel_estimates_new.append(MutRateStudy(f"New", 2025, "Indirect", UCSCname, mutRate_avg, mutRate_std_lb, mutRate_std_ub))
		print(f"\t{UCSCname} : {mutRate_avg}\t{mutRate_std_lb}\t{mutRate_std_ub}")

	return indel_estimates_new

###############################
# Methods to load info from other studies.
# Literature review on indel estimates.

def getDirectEstimates():
	""" It returns a list where each entry corresponds to a previous study.
	All studies returned by this method are *Direct Estimates*, i.e., they count mutations that occur between generations in present-day individuals.

	.. note::
		The indel sizes vary in each study (see entry *Indel sizes* in each tuple). As lower the maximum indel size is, higher the mutation rate estimation.

	"""
	indel_estimates_direct = [	MutRateStudy(f"Kloosterman et al.", 2015, "Direct", "hg38", 2.3231886903593173e-11, 1.990918179344908e-11, 2.788584149604477e-11),
								MutRateStudy(f"Besenbacher et al.", 2016, "Direct", "hg38", 3.07*(10**(-11)), 2.91*(10**(-11)), 3.25*(10**(-11))),
								MutRateStudy(f"Maretty et al.",     2017, "Direct", "hg38", 4.70e-11, None, None),
								MutRateStudy(f"Besenbacher et al.", 2015, "Direct", "hg38", 5.28169014084507e-11, 4.225352112676056e-11, 6.690140845070423e-11),
								#MutRateStudy(f"Montgomery et al.", 2013, "Direct", "hg38", None, None, None), # indel sizes: (1, 50)
								MutRateStudy(f"Kondrashov (del.)",  2002, "Direct", "hg38", 2.63e-11, 1.58e-11, 4.18e-11),
								MutRateStudy(f"Kondrashov (ins.)",  2002, "Direct", "hg38", 0.91e-11, 0.36e-11, 1.46e-11),
								MutRateStudy(f"Palamara (ins.)",    2015, "Direct", "hg38", 4.3448275862068967e-11, 4.137931034482759e-11, 4.5517241379310344e-11)
							]
	return indel_estimates_direct

def getIndirectEstimates():
	""" It returns a list where each entry corresponds to a previous study.
	All studies returned by this method are *Indirect Estimates*, i.e., they compute their estimate based on the evolutionary distance and the divergence time separating two species.
	"""
	indel_estimates_indir  = [ 	MutRateStudy(f"Nachman \n and Crowell", 2000, "Indirect", "panTro6", 4.95049504950495e-11, 3.712871287128713e-11, 6.188118811881188e-11),
								MutRateStudy(f"Lunter",                 2007, "Indirect", "mm39",    30.46e-11, 29.12087912087912e-11, 32.595325953259533e-11)
							]
	return indel_estimates_indir

def getExtrapolatedEstimates():
	""" It returns a list where each entry corresponds to a previous study.
	All studies returned by this method are *Extrapolated Estimates*, i.e., they estimate the indel rate based on the substitution rate.

	These studies have the generation time unclear and, therefore, were **left out of the analysis**.
	"""
	indel_estimates_extrap = [ (f"Lynch (del.)",		 2010, "Extrap.",  (1,  50),	0.58*(10**(-9)),	 (None, None),		25,				   (20, 30),			   2.32e-11, (1.933333333333333e-11, 2.8999999999999997e-11), "hg38"),
							   (f"Lynch (ins.)",		 2010, "Extrap.",  (1,  50),	0.20*(10**(-9)),	 (None, None),		25,				   (20, 30),			   0.80e-11, (0.66e-11,1e-11), "hg38"),
							   (f"Ramu et al.",		     2013, "Extrap.",  (1, None),   1.06*(10**(-9)), (0.235e-09, 2.75e-09),   25,				   (20, 30), 4.2400000000000004e-11, (3.5333333333333334e-11, 5.3000000000000004e-11), "hg38"),
							 ]
	return indel_estimates_extrap

###############################
# Plot methods.

def plotErrorRect(ax, study, ycoord, color):
	rect = (None,None,None,None)
	if(study.mutRatePPPY_lb and study.mutRatePPPY_ub):
		rec_height = 0.6 # 0.6
		rec_width  = study.mutRatePPPY_ub-study.mutRatePPPY_lb
		xrecbeg    = study.mutRatePPPY_lb
		yrecbeg    = ycoord - rec_height/2
		ax.add_patch(Rectangle( (xrecbeg, yrecbeg), rec_width, rec_height,color=color, alpha=0.5))
		rect       = (xrecbeg,yrecbeg,rec_width,rec_height)
	return rect

def plotIcon(ax, iconFilename, xval, yval, rect):
	xrecbeg,yrecbeg,rec_width,rec_height = rect
	space   = 2e-11
	xposbeg = xval+space if (yval % 2 == 0) else xval-space
	if (xrecbeg):
		xposbeg  += (rec_width/2 if (yval % 2 == 0) else (-rec_width/2))
	skunkName = f"{os.path.basename(iconFilename)}-{random.randint(1, 1000)}"
	skunkIm   = skunk.Box(30, 30, skunkName)
	ab = AnnotationBbox(skunkIm, (xposbeg, yval), xycoords='data', frameon=False) # 'axes fraction'
	ax.add_artist(ab)
	return skunkName

def createLines(ax, my_dataset, rows, corr=0):
	""" Creates horizontal lines separating direct and indirect estimates.
	"""
	# Information of previous row (first row).
	previous_row_method, previous_row_ucscName, previous_row_author  = (rows[0]).methodType, (rows[0]).ucscName, (rows[0]).author
	# Check rows, starting at the second row.
	for rowIdx, study in zip(range(1,len(rows)),rows[1:]):
		cur_row_method, cur_row_ucscName, cur_row_author  = study.methodType, study.ucscName, study.author
		# Check if the type of method (direct, indirect) changed.
		if(previous_row_method != cur_row_method):
			# Plot a gray thick line separating Direct/Indirect methods.
			ax.axhline(rowIdx+corr, color='#999999', lw=4)	
			ax.hlines(rowIdx+corr, 0, -1.00, color='#999999', lw=4, clip_on=False, transform=ax.get_yaxis_transform())
		# Update previous info.
		previous_row_method, previous_row_ucscName, previous_row_author = cur_row_method, cur_row_ucscName, cur_row_author

def yLabels(ax, my_dataset, rows, corr=0):
	""" Creates y-labels.

	The y-label consists of the name of the species that was compared with human 
	(or "Human" if it is a direct comparison), and from which study (author + publication year)
	the estimate comes from.

	"""

	# Remove current labels.
	ax.set_ylabel('')
	
	# Labels.
	fontSizeLbl  = 16  # "xx-small"
	charsperline = 20 # 20
	
	# How rows are organized:
	# (0) Estimates are sorted by direct / indirect methods;
	# (1) Inside each method type, they are sorted by how far they are from humans;
	# (2) Grouped by species;
	# (3) If there is more than one publication per species, then they are sorted by publication year.
	ucscName_prev,   ucscName_prev_rowIdx   = "", 0
	methodType_prev, methodType_prev_rowIdx = "", 0
	
	methodType_xpos = -0.99
	ucscName_xpos   = -0.70
	author_xpos     = -0.30
	for rowIdx, study in enumerate(rows):

		methodType_cur = study.methodType
		ucscName_cur   = my_dataset.commonNames[study.ucscName] if (study.ucscName != my_dataset.refsp_ucscName) else "Human"
		author_cur     = f"{study.author} ({study.publYear})"

		ax.text(author_xpos, rowIdx+0.5+corr, author_cur, fontsize=fontSizeLbl, ha='center', va='center', clip_on=False, transform=ax.get_yaxis_transform())
	
		# Time to print previous species.
		# The name of the species is printed only once, vertically aligned in the center of all estimates regarding the species.
		if(ucscName_prev != ucscName_cur):
			if(ucscName_prev):
				ax.text(ucscName_xpos, rowIdx-((rowIdx-ucscName_prev_rowIdx)/2)+corr, textwrap.fill(ucscName_prev, charsperline), fontsize=fontSizeLbl, ha='center', va='center',
						clip_on=False, transform=ax.get_yaxis_transform())
			ucscName_prev, ucscName_prev_rowIdx = ucscName_cur, rowIdx

		# Time to print previous type of the method (Direct or Indirect).
		if(methodType_prev != methodType_cur):
			if(methodType_prev):
				ax.text(methodType_xpos, rowIdx-((rowIdx-methodType_prev_rowIdx)/2)+corr, textwrap.fill(methodType_prev, charsperline), fontsize=fontSizeLbl, ha='center', va='center',
						clip_on=False, transform=ax.get_yaxis_transform())
			methodType_prev, methodType_prev_rowIdx = methodType_cur, rowIdx
			
	# Last labels.
	rowIdx = len(rows)
	if(ucscName_prev):
		ax.text(ucscName_xpos, rowIdx-((rowIdx-ucscName_prev_rowIdx)/2)+corr, textwrap.fill(ucscName_prev, charsperline), fontsize=fontSizeLbl, ha='center', va='center',
				clip_on=False, transform=ax.get_yaxis_transform())
	if(methodType_prev):
		ax.text(methodType_xpos, rowIdx-((rowIdx-methodType_prev_rowIdx)/2)+corr, textwrap.fill(methodType_prev, charsperline), fontsize=fontSizeLbl, ha='center', va='center',
				clip_on=False, transform=ax.get_yaxis_transform())

def plotMutationRatePerType(ax, mutRateEstsAll, mutRateEsts, markerStyle):
	skunkMap = {}
	markerColor, markerType, markerSize = markerStyle

	# Plot data.
	xvals = [study.mutRatePPPY for study in mutRateEsts]
	yvals = [mutRateEstsAll.index(study) for study in mutRateEsts]
	ax.scatter(xvals, yvals, marker=markerType, color=markerColor, s=markerSize)

	# Plot additional info.
	for (study, xval, yval) in zip(mutRateEsts,xvals,yvals):

		# Plot rectangles (error boxes).
		rect = plotErrorRect(ax, study, yval, markerColor)

		# Plot icons.
		if (study.ucscName != my_dataset.refsp_ucscName):
			iconFilename = my_dataset.getIconFilename(study.ucscName)
			if(iconFilename): 
				skunkName = plotIcon(ax, iconFilename, xval, yval, rect)
				skunkMap[skunkName] = iconFilename
	return skunkMap

def plotMutationRateComparison(my_dataset,new_ests,direct_ests,indirect_ests, empty_evoltime_quantile):

	data = [("New", new_ests), ("Direct", direct_ests), ("Indirect", indirect_ests)]
	markerStyles  = {"New"     : ('darkorange',  'o', 100), 
					"Indirect" : ('yellowgreen', 'o', 100), 
					"Direct"   : ('royalblue',   '*', 200)}

	# Gather all estimates, regardless of their type (direct, indirect or new), and sort them by the following criteria: 
	# (1) divergence time, (2) name of species and (3) publication year.
	# This order determines their position in the y-axis.
	mutRates_all = sorted(new_ests+direct_ests+indirect_ests, key=lambda study: (0 if (study.ucscName == my_dataset.refsp_ucscName) else my_dataset.divergenceTimes[study.ucscName], study.ucscName, int(study.publYear)))

	# Plot estimates.
	fig, ax = plt.subplots(figsize=(18,22),tight_layout=True)
	plt.rcParams['font.sans-serif'] = ['DejaVu Sans']

	# Find image size.
	skunkMap = {}

	# Group estimates per type (new, direct, or indirect).
	for (groupType, mutRateEsts) in data:
		skunkMap.update(plotMutationRatePerType(ax, mutRates_all, mutRateEsts, markerStyles[groupType]))

	ax.set(yticklabels=[])  # remove the tick labels
	ax.tick_params(left=False)
	
	# x-axis.
	ax.ticklabel_format(axis='x', scilimits=[-11, -11])
	ax.tick_params(axis='x', which='major', labelsize=16)
	ax.tick_params(axis='x', which='minor', labelsize=16)
	ax.xaxis.get_offset_text().set_fontsize(16)
	ax.xaxis.get_offset_text().set_verticalalignment("baseline")
	
	# y-labels.
	yLabels(ax, my_dataset, mutRates_all, corr=0.0) #-0.5
	# Line separating "Direct" and "Indirect" methods.
	#createLines(ax, my_dataset, mutRates_all, corr=-0.5)
	
	plt.title("Indel rates comparison", fontsize=18)
	posxlabel = -0.00000000008
	ax.set_xlabel("Estimated indel rate (PPPY - per position per year)", fontsize=16, ha='left', x=ax.transAxes.inverted().transform(ax.transData.transform((posxlabel, 0)))[0])
	#ax.set_xlim((0,max([study.mutRatePPPY_ub if (study.mutRatePPPY_ub) else study.mutRatePPPY for study in mutRates_all])))

	# Avg. line
	ax.axvline(x=np.mean([new_est.mutRatePPPY for new_est in new_ests]),color='orange', linestyle='--')
	
	# Add SVG images.
	plt.subplots_adjust(left=0.1,right=0.8)
	svgFile = skunk.insert(skunkMap)

	# Save file.
	plotFilename = my_dataset.getOutFilename_plot_MutRatesComp(empty_evoltime_quantile, "svg")
	cairosvg.svg2svg(bytestring=svgFile, write_to=plotFilename)

def makeFigure4(alpha, my_dataset):
	
	# Evolutionary time for empty windows: 
	# High tail of the evolutionary time distribution.
	empty_evoltime_quantile = 0.99

	# Load new and previous estimated indel mutation rates.
	new_ests      = getNewEstimates(my_dataset, alpha, empty_evoltime_quantile=empty_evoltime_quantile)
	direct_ests   = getDirectEstimates()
	indirect_ests = getIndirectEstimates()

	# Plot data.
	plotMutationRateComparison(my_dataset,new_ests,direct_ests,indirect_ests, empty_evoltime_quantile)

####################################
# MAIN.
####################################
if (__name__ == '__main__'):
	alpha	   = 1.1
	my_dataset = Dataset()

	# Save overall time.
	timeTrack = Time()
	timeTrack.start()

	timeTrack.startStep("Figure 4")
	makeFigure4(alpha, my_dataset)
	timeTrack.stopStep()

	timeTrack.stop()
	timeTrack.print()

	print(f"Done!")
