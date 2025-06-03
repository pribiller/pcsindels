import os
import pickle
import random
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.backends.backend_pdf import PdfPages
import textwrap
import cairosvg
import skunk
from colorir import Grad

from utils.dataset import Dataset
from utils.basicTypes import Time

###############################
# Methods to prepare data for plotting (averaging, binning, etc.).

def findBin(val,bins):
	""" Given a list of bins and a value, find in which bin the value should go.

		:param bins: a **sorted** list of values. A bin is defined as (bins[i-1], bins[i]].
		:type bins: list of floats.

	"""
	for binIdx, b in enumerate(bins):
		if(val <= b):
			return (binIdx, b)
	return (len(bins)-1,bins[-1])

def prepareDataForPlotting(dictVals):
	""" Input is provided as a dictionary containing binned data, 
		with keys representing bin values (like PCS sizes) and 
		values indicating the counts of PCS sizes within those bins. 
		The process computes the average and standard deviation of 
		the observed counts for each bin.
	"""
	xvals = list(sorted(dictVals.keys()))
	yvals = np.array([np.mean(dictVals[xval]) for xval in xvals])
	yvals_err_ub = np.array([np.std(dictVals[xval]) for xval in xvals])
	yvals_err_lb = np.array([np.std(dictVals[xval]) for xval in xvals])
	return (xvals,yvals,yvals_err_lb,yvals_err_ub)

###############################
# Plot methods.

def plotPCSsizeDistribBig(ax,PCSdist_obs_onespecies,PCSdist_est_onespecies, min_PCSsize, PCSsize_nbbins):
	""" Plot the observed and estimated distribution of PCS sizes for a pair of species.
	Due to the limited space, in order to improve clarity and avoid overlapping dots, 
	the x-axis values are organized into bins.
	"""	

	# Create bins for plots.
	max_PCSsize = max(max([max(PCSdistrib.keys()) for PCSdistrib in PCSdist_est_onespecies]),max(PCSdist_obs_onespecies.keys()))
	# Log bins.
	PCSsize_bins = np.logspace(np.log10(min_PCSsize),np.log10(max_PCSsize), PCSsize_nbbins)
	PCSsize_bins = list(sorted(set([int(val) for val in PCSsize_bins])))

	N_obs = sum(list(PCSdist_obs_onespecies.values()))
	N_est = [sum(list(PCSdist_est.values())) for PCSdist_est in PCSdist_est_onespecies]

	print(f"{N_obs=}")
	print(f"{N_est=}")

	# Estimated data.
	PCSdistrib_est_all_binned = defaultdict(list)
	for PCSdist_est in PCSdist_est_onespecies:
		for (PCSsize, PCScnt) in PCSdist_est.items():
			# First bin PCS cnts, then compute average and standard deviation later.
			PCSsize_binIdx, PCSsize_binVal = findBin(PCSsize,PCSsize_bins)
			PCSdistrib_est_all_binned[PCSsize_binVal].append(PCScnt)
	xvals_est,yvals_est,yvals_est_err_lb,yvals_est_err_ub = prepareDataForPlotting(PCSdistrib_est_all_binned)

	# Observed data.
	PCSdistrib_obs_all_binned = defaultdict(list)
	for (PCSsize,PCScnt) in PCSdist_obs_onespecies.items():
		PCSsize_binIdx, PCSsize_binVal = findBin(PCSsize,PCSsize_bins)
		PCSdistrib_obs_all_binned[PCSsize_binVal].append(PCScnt)
	xvals_obs,yvals_obs,yvals_obs_err_lb,yvals_obs_err_ub = prepareDataForPlotting(PCSdistrib_obs_all_binned)
	
	# Plot data.
	ax.errorbar(xvals_obs, yvals_obs, yerr=[yvals_obs_err_lb,yvals_obs_err_ub], markersize=5, fmt='o', color='royalblue',  ecolor='royalblue',  elinewidth=8, capsize=15, markeredgewidth=10)
	ax.errorbar(xvals_est, yvals_est, yerr=[yvals_est_err_lb,yvals_est_err_ub], markersize=5, fmt='o', color='darkorange', ecolor='darkorange', elinewidth=8, capsize=15, markeredgewidth=10)

	ax.fill_between(xvals_obs,yvals_obs-yvals_obs_err_lb,yvals_obs+yvals_obs_err_ub,facecolor='royalblue',edgecolor="none",alpha=0.2)
	ax.fill_between(xvals_est,yvals_est-yvals_est_err_lb,yvals_est+yvals_est_err_ub,facecolor='darkorange',edgecolor="none",alpha=0.2)
	
	ax.set_xscale("log")
	ax.set_yscale("log")

	ax.set_xlabel("PCS size (log)")
	ax.set_ylabel("Count (log)")
	# change font size for all axis
	ax.xaxis.get_label().set_fontsize(60)
	ax.yaxis.get_label().set_fontsize(60)
	ax.set_facecolor("#F5F5F5") # "whitesmoke"
	ax.tick_params(axis='both', labelsize=45)
	for spine in ax.spines.values():
		spine.set_visible(False)

def plotTauDistribSampled_sum(ax_fg,ax_bg,taudistrib_est_onespecies,tau_bins,min_tau,max_tau,colors,cmap):
	""" Plot the estimated distribution of evolutionary times for a pair of species.
	"""	

	# Prepare data.
	taudistrib_est_all_binned = defaultdict(list)
	for taudistrib_est in taudistrib_est_onespecies:
		# Bin tau distribution of one sample.
		taudistrib_est_one_binned = {tauBin: 0 for tauBin in tau_bins}
		for (tauVal, tauCnt) in taudistrib_est.items():
			tau_binIdx, tau_binVal = findBin(tauVal,tau_bins)
			taudistrib_est_one_binned[tau_binVal] += tauCnt
		# Aggregate binned results.
		for (tauBin, tauCnt) in taudistrib_est_one_binned.items():
			taudistrib_est_all_binned[tauBin].append(tauCnt)
			
	xvals,yvals,yvals_err_lb,yvals_err_ub = prepareDataForPlotting(taudistrib_est_all_binned)
	
	# Plot data.
	y_max	  = np.max(yvals)
	bar_height = y_max*0.10
	bar_xcoord = -bar_height
	
	filling	= ax_bg.imshow(  [[0.,1.], [0.,1.]], # Image data.
						  origin='lower', aspect='auto', cmap=cmap,
						  extent=[tau_bins[0], tau_bins[-1], bar_xcoord, 0.0], # The bounding box in data coordinates that the image will fill.
						  interpolation = 'bicubic',zorder=-1)
	
	for idx, (x,y,y_err_lb,y_err_ub,c) in enumerate(zip(xvals, yvals, yvals_err_lb, yvals_err_ub, colors)):
		ax_fg.errorbar(x, y, yerr=y_err_lb, fmt='o', markersize=6,capsize=10, color=c, ecolor=c)
		x0, y0 = (xvals[idx-1], yvals[idx-1]) if (idx > 0) else (0,0)
		ax_fg.fill_between([x0,x],[y0,y],color=c,alpha=0.25)
	
	ax_fg.set_xlim((min_tau,max_tau))
	ax_fg.set_xscale('log')
	ax_fg.set_ylim(bar_xcoord,y_max*1.05)

	ax_bg.axis('off')
	ax_bg.grid(False)
	ax_bg.set_xticks([]) # for major ticks
	ax_bg.set_xticks([], minor=True) # for minor ticks
	
	ax_fg.axis('off')
	ax_fg.grid(False)
	plt.text(0.50, -0.30, "Evol.times",
		 horizontalalignment='center',
		 fontsize=55,
		 transform = ax_fg.transAxes)
	
	# Remove the spines (lines around the plot).
	for spine in ax_fg.spines.values():
		spine.set_visible(False)
	for spine in ax_bg.spines.values():
		spine.set_visible(False)
		
	ax_fg.set_facecolor('#F5F5F5')

def makeFigure2(alpha, my_dataset):

	print("Plotting graphs to be used in papers...")

	# Colorbar for plots.
	gradEvolTimes = Grad(["#1a9850","#ffcc00","#d73027"])
	cmapEvolTimes = gradEvolTimes.to_cmap()

	# Parameters for plots.
	PCSsize_nbbins = 20
	PCSsize_min	= my_dataset.minPCSsize

	# Log bins.
	t_min	 = 5e-5
	t_max	 = 0.50
	t_bins_nb = 30
	t_bins	  = np.logspace(np.log10(t_min),np.log10(t_max), t_bins_nb)

	# Colors for tau distribution.
	color_min = 0.001 # min(t_bins) 
	color_max = 0.02  # max(t_bins)
	colors = [max(0.0,min(1.0,(t-t_min)/(t_max-t_min))) for t in t_bins]
	colors = [t_idx/len(t_bins) for t_idx, t in enumerate(t_bins)]
	colors = [cmapEvolTimes(normtau) for normtau in colors]

	plotFilename = my_dataset.getOutFilename_plot_PcsDistribComp(alpha)
	# pp = PdfPages(plotFilename)
	
	# Create grid for species.
	nb_rows = 6
	nb_cols = 7
	speciesGrid = [  ['panPan3', 'panTro6', 'gorGor6', 'ponAbe3',  'papAnu4'],
					 ['macFas5', 'rhiRox1', 'chlSab2', 'nasLar1',  'rheMac10', 'calJac4'],
					 ['tarSyr2', 'micMur2', 'galVar1', 'mm39',	 'rn7',	  'oryCun2'],
					 ['vicPac2', 'bisBis1', 'felCat9', 'manPen1',  'bosTau9',  'canFam6'],
					 ['musFur1', 'neoSch1', 'equCab3', 'myoLuc2',  'susScr11', 'enhLutNer1'],
					 ['triMan1', 'macEug2', 'ornAna2', 'aptMan1',  'galGal6',  'thaSir1'],
					 ['aquChr2', 'melGal5', 'xenLae2', 'xenTro10', 'danRer11']]
	# For test.
	# speciesGrid = [  ['panPan3'],
	#				  ['macFas5', 'rhiRox1'],
	#				  ['tarSyr2'],
	#				  ['vicPac2', 'bisBis1'],
	#				  ['musFur1'],
	#				  ['triMan1', 'macEug2'],
	#				  ['aquChr2']]
	
	# One graph for the PCS size distribution, and another 
	# graph for evolutionary time distribution (inset).
	print(f"\t Plot data...")
	fig=plt.figure(figsize=(12*(nb_cols+1), 10*(nb_rows+1)))

	subfigs  = fig.subfigures(nrows=nb_rows, ncols=nb_cols, wspace=0.05, hspace=0.10)
	skunkMap = {}
	for colIdx, groupSpecies in enumerate(speciesGrid):
		for rowIdx, UCSCname in enumerate(groupSpecies):

			divTime	   = my_dataset.divergenceTimes[UCSCname]
			commonName = my_dataset.commonNames[UCSCname]
			
			print(f"> {commonName} ({UCSCname}); divergenge-time estimate: {divTime} mya.")
			
			# Load sampled distributions.
			sampEvolTimesFilepath = my_dataset.getOutFilename_sampleEvolTimes(UCSCname, alpha)
			if(not os.path.isfile(sampEvolTimesFilepath)):
				print(f"[{UCSCname}] WARNING! File not found: {sampEvolTimesFilepath}. Skipping computation!")
				continue
			PCSdistrib_obs_whl, PCSdistrib_est_whl, taudistrib_est_whl, PCSdistrib_obs_chr, PCSdistrib_est_chr, taudistrib_est_chr = pickle.load(open(sampEvolTimesFilepath, 'rb'))

			# Create graphs for species.
			subfig = subfigs[rowIdx][colIdx]
			#subfig.suptitle(f"{commonName} ({UCSCname}); divergenge-time estimate: {divTime} mya.",fontsize=32)
			ax=subfig.subplots(1, 1)

			# Plot PCS size distribution.
			plotPCSsizeDistribBig(ax,PCSdistrib_obs_whl,PCSdistrib_est_whl, PCSsize_min, PCSsize_nbbins)
			
			# Plot tau distribution.

			# Create an inset of width 30% and height 40% of the parent 
			# axes' bounding box at the lower left corner (loc=3)
			ax_inset_fg = inset_axes(ax, width="100%", height="100%", loc=3, bbox_to_anchor=(.05, .12, .55, .30), bbox_transform=ax.transAxes) # (start_x,start_y,width,height)
			ax_inset_bg = ax_inset_fg.twiny()

			plotTauDistribSampled_sum(ax_inset_fg,ax_inset_bg,taudistrib_est_whl,t_bins,t_min,t_max,colors,cmap=cmapEvolTimes)
			
			# Adjust z-order of plots.
			ax.grid(False)
			ax.set_zorder(ax_inset_fg.get_zorder()-1)
			# ax.patch.set_visible(False)
			ax.set_facecolor("#F5F5F5")
			
			# Plot icons.
			iconFilename = [f for f in os.listdir(my_dataset.dirIcons) if ((UCSCname in f) and (f.endswith("svg")))]
			iconFilename = os.path.join(my_dataset.dirIcons, iconFilename[0]) if (len(iconFilename) > 0) else ""
			if(iconFilename):
				skunkName = f"{os.path.basename(iconFilename)}-{random.randint(1, 1000)}"
				newax = subfig.add_axes([0.62, 0.60, 0.30, 0.50], anchor='C') # zorder=-1 / [0:left / 1:right, 0:bottom / 1:top, w, h]
				charsperline  = 12
				
				breakSmallText = ((len(commonName) <= charsperline) and (" " in commonName))
				breakBigText   = (len(commonName) > charsperline)
				
				formattedName = textwrap.fill(commonName, charsperline) if (not breakSmallText) else "\n".join(commonName.split(" "))
								
				ypos  = -0.20 if (breakSmallText or breakBigText) else -0.10
				newax.set_title(f"{formattedName}", y=ypos, fontsize=55)
				newax.axis('off')
				newax.set_zorder(ax.get_zorder()+1)
				skunk.connect(newax, skunkName)
				skunkMap[skunkName] = iconFilename

	# Add SVG images.
	# svg = skunk.insert(skunkMap)
	# skunk.display(svg)

	# cairosvg.svg2pdf(bytestring=svg, write_to=plotFilename)
	
	# plt.show()
	# pp.savefig(fig)
	# pp.close()

####################################
# MAIN.
####################################
# Use: python3 ~/code/6_1_plotFig_PcsDistribComp.py
if (__name__ == '__main__'):
	alpha	   = 1.1
	my_dataset = Dataset()

	# Save overall time.
	timeTrack = Time()
	timeTrack.start()

	makeFigure2(alpha, my_dataset)

	timeTrack.stop()
	timeTrack.print()

	print(f"Done!")

