"""Computes the PCS size distribution for every window 
in a chromosome, considering all species present in the dataset.

- **Use**::

	python3 2_computeWindows.py -chr [chromosome name]

- **Example of Usage (vertebrate dataset with 40 species)**::

	python3 ~/code/2_computeWindows.py -chr chr16

- **Input Parameter (mandatory)**:

:-chr: Chromosome name in the reference species (e.g. chr1, chr2, ..., chrX, chrY).

- **Other Parameters taken from** ``dataset.py``:
:-refsp_ucscname: UCSC name of the reference species that is being aligned (e.g. *hg38* for human).
:-win_size: Size of windows in base pairs.
:-species_UCSC_names: list of UCSC names of all species included in the dataset (40 vertebrates in our study).
:-pcs_dir: Directory where the PCSs for each pairwise alignment included in the dataset are saved.
:-win_dir: Directory where windows and their PCS size distributions will be saved.

.. note::
	Make sure that the required parameters described above are correctly defined in the file ``utils/dataset.py``.

- **Output**: 
	One ``.pickle`` file for the specified chromosome with a dictionary.
	In this dictionary, the keys are the UCSC names of the species, 
	while the values are dictionaries that map window coordinates 
	on the chromosome to their PCS size distributions.
	These distributions are also organized in a dictionary format, where 
	the keys indicate size and the values denote counts.
	
Pre-requisites
--------------

	Before using this script, make sure all the required files were computed:

a) PCSs from each pairwise alignment in the dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Make sure to run ``1_extractPCS.py`` for **every** species included in 
the attribute ``speciesUCSCnames`` declared in file ``dataset.py``. 
In our study, this list includes the UCSC names of the 40 vertebrate species 
used in our analysis, whose data (fasta files and chain files) was retrieved 
from the UCSC website.

Cluster resources
-----------------

In case you want to run this Python script stand-alone::

	srun -p compute -t 2:00:00 --mem 20G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash

Otherwise you can use the script ``../cluster/2_computeWindows_runAll.py``
to run this script for all 24 chromosomes (chr1, chr2, ..., chrX, chrY).

Time, Memory & Disk space
-------------------------

For reference, here we include an **upper limit** on runtime, memory usage, and disk 
space required for running this script on the 40 vertebrate dataset examined in our study.

======  =========  =========  ========
Desc.   Time       Memory     Disk    
======  =========  =========  ========
chr1    01:09:13   9GB        0.340GB 
chr2    01:29:36   11GB       0.360GB 
chr3    00:59:03   7GB        0.280GB 
chr4    01:01:40   7GB        0.280GB 
chr5    00:54:51   7GB        0.260GB 
chr6    01:00:43   5GB        0.230GB 
chr7    00:51:09   6GB        0.190GB 
chr8    00:47:49   5GB        0.200GB 
chr9    00:32:51   3GB        0.150GB 
chr10   00:36:37   4GB        0.170GB 
chr11   00:29:28   5GB        0.190GB 
chr12   00:30:16   5GB        0.190GB 
chr13   00:20:57   3GB        0.140GB 
chr14   00:19:24   3GB        0.140GB 
chr15   00:17:15   3GB        0.120GB 
chr16   00:17:42   3GB        0.120GB 
chr17   00:18:00   3GB        0.100GB 
chr18   00:21:40   3GB        0.110GB 
chr19   00:12:39   2GB        0.080GB 
chr20   00:18:12   3GB        0.100GB 
chr21   00:10:03   2GB        0.060GB 
chr22   00:07:28   1GB        0.050GB 
chrX    00:35:28   5GB        0.200GB 
chrY    00:04:18   1GB        0.030GB 
======  =========  =========  ========


Time per Run: Details
^^^^^^^^^^^^^^^^^^^^^

Stats on time of a single run (chr16, 40 vertebrate species): **~15 minutes**

More details on computational time can be found in the log of the run.

================================  =============
Step                              Time (s)     
================================  =============
Merging PCSs                      715.09147525 
Computing windows                 0.57280803   
Computing PCS size distribution   172.84002709 
**Total time**                    902.55868173 
================================  =============

Storage per Run: Details
^^^^^^^^^^^^^^^^^^^^^^^^

Size of output files with all windows of one chromosome (chr16): **~134 MB**.

==================================  ======
Output files                        Size
==================================  ======
hg38.chr16.1000.windows.pickle      119M
hg38.chr16.mergedPCSs.pickle         15M
==================================  ======

.. note::
	The temporary file ``hg38.chr16.mergedPCSs.pickle`` can be removed after the output file ``hg38.chr16.1000.windows.pickle`` is succesfully computed.

Function details
----------------

Only relevant functions have been documented below. 
For more details on any function, check the comments in the souce code.

"""
import os
import pickle
from collections import Counter, defaultdict
import sys
import re
import random
import itertools
import argparse

from utils.basicTypes import Pcs, Time
from utils.dataset import Dataset

####################################
# "Merge PCSs" related code

def isOverlap(pcs1,pcs2):
	b1, e1, b2, e2 = (pcs1.qPosBeg, pcs1.qPosBeg+pcs1.size, pcs2.qPosBeg, pcs2.qPosBeg+pcs2.size)
	return (((b1 <= b2) and (e1 > b2)) or ((b2 <= b1) and (e2 > b1)))

def isContiguous(pcs1,pcs2):
	b1, e1, b2, e2 = (pcs1.qPosBeg, pcs1.qPosBeg+pcs1.size, pcs2.qPosBeg, pcs2.qPosBeg+pcs2.size)
	return ((b1 == e2) or (b2 == e1))

def printPCS(p):
	return f"[{p.qPosBeg}\t{p.qPosBeg+p.size}); {p.size}"

def mergePCS_pairwise_check(pcs_lst_cur, idx_to_add):
	"""This function checks if the new PCS was properly added to the list. Two constrains are checked:

		1. The new PCS should not overlap the previous or the next PCS;
		2. The list must keep its property of being sorted by position.

	"""
	p_prev = pcs_lst_cur[idx_to_add-1] if (idx_to_add > 0) else None
	p_new  = pcs_lst_cur[idx_to_add]
	p_next = pcs_lst_cur[idx_to_add+1] if (idx_to_add < (len(pcs_lst_cur)-1)) else None

	listOK = (((p_prev == None) or ((p_prev.qPosBeg + p_prev.size) <= p_new.qPosBeg ))
		  and ((p_next == None) or ((p_new.qPosBeg  + p_new.size)  <= p_next.qPosBeg)))

	if (not listOK):
		print(f"ERROR! Problem to add a PCS to the list of merged PCSs:")
		print(f"PCS before: {printPCS(p_prev)}")
		print(f"PCS new:    {printPCS(p_new)}")
		print(f"PCS after:  {printPCS(p_next)}")
		sys.exit()

def mergePCS_pairwise_findPos(p_new, pcs_lst_cur, idx_to_add):
	"""This function determines the index to insert the new PCS 
	while keeping the list sorted and the PCSs non-overlapping.
	"""
	searchForPosition = True
	nb_pcs_merged = 0
	posBeg_merge  = -1
	while (searchForPosition):
		p_cur = pcs_lst_cur[idx_to_add]
		if(isOverlap(p_cur,p_new) or isContiguous(p_cur,p_new)): 
		#if isOverlap(p_cur,p_new): 
			nb_pcs_merged += 1
			if(posBeg_merge < 0): posBeg_merge = idx_to_add
		searchForPosition = ((p_new.qPosBeg + p_new.size) > p_cur.qPosBeg)
		if(searchForPosition): 
			idx_to_add += 1
			searchForPosition = (idx_to_add < len(pcs_lst_cur))
	# Adjust index if needed.
	idx_to_add = posBeg_merge if(nb_pcs_merged > 0) else idx_to_add
	return idx_to_add, nb_pcs_merged

def mergePCS_pairwise_updLst(p_new, pcs_lst_cur, idx_to_add, nb_pcs_merged):
	"""This function updates the PCS list with a new PCS (if needed).
	"""
	if(nb_pcs_merged > 0):
		posBeg_merge  = idx_to_add
		posEnd_merge  = posBeg_merge+nb_pcs_merged
		pcs_lst_merge = pcs_lst_cur[posBeg_merge:posEnd_merge] + [p_new]
		
		qPosBeg = min([p.qPosBeg for p in pcs_lst_merge])
		qPosEnd = max([p.qPosBeg+p.size for p in pcs_lst_merge])
		sizePCS = qPosEnd-qPosBeg

		pcs_largest = max(pcs_lst_merge, key=lambda pcs: pcs.size)
		pcs_merged  = Pcs(sizePCS, pcs_largest.tChrom, pcs_largest.tStrand, pcs_largest.tPosBeg, pcs_largest.qChrom, pcs_largest.qStrand, qPosBeg)

		# Remove PCSs from list in case more than one PCS is merged. 
		# Otherwise, just replace the PCS by the merged PCS.
		if(nb_pcs_merged > 1): del pcs_lst_cur[posBeg_merge:posEnd_merge]
	# Add pcs.
	p_add = pcs_merged if (nb_pcs_merged > 0) else p_new
	if(nb_pcs_merged == 1): # Replace previous PCS by merged PCS.
		pcs_lst_cur[idx_to_add] = p_add
	elif (idx_to_add < len(pcs_lst_cur)):
		pcs_lst_cur.insert(idx_to_add, p_add)
	else:
		pcs_lst_cur.append(p_add)
	return pcs_lst_cur

def mergePCS_pairwise(pcs_lst_cur,pcs_lst_new):
	"""This function adds a list of new PCSs to the current 
	list of already-merged PCSs. 

	A new PCS can be:

		* *Discarded*, if it is encompassed by another PCS in the list;
		* *Added*, if its position does not overlap any PCS position in the list;
		* *Merged*, if its position overlaps one or more PCSs in the list.

	"""

	# Make sure new PCSs are sorted.
	pcs_lst_new = sorted(pcs_lst_new, key=lambda pcs: pcs.qPosBeg)

	# Check if current list is empty.
	if(not pcs_lst_cur): return pcs_lst_new

	# New PCS can be: 
	# - Discarded, if it is encompassed by another PCS in the list;
	# - Added, if its position does not overlap any PCS position in the list;
	# - Merged, if its position overlaps one or more PCSs in the list.
	idx_to_add = 0
	pcs_lst_cur_len = len(pcs_lst_cur)
	print(f"\t{len(pcs_lst_new)} new PCSs will be checked...")

	nb_pcss_merged = defaultdict(int)
	debug=False
	for p_idx, p_new in enumerate(pcs_lst_new):
		nb_pcs_merged = 0

		# Find in which index the new PCS should be included.
		# List is sorted by position.
		idx_to_add, nb_pcs_merged = mergePCS_pairwise_findPos(p_new, pcs_lst_cur, idx_to_add)

		# Update list (if needed).
		pcs_lst_cur = mergePCS_pairwise_updLst(p_new, pcs_lst_cur, idx_to_add, nb_pcs_merged)

		# Check if list is still sorted after modifications.
		mergePCS_pairwise_check(pcs_lst_cur, idx_to_add)

		# Save some stats on merged PCSs.
		nb_pcss_merged[nb_pcs_merged] += 1

	print(f"\tStats on merged PCSs")
	print("\n".join("\t{}\t{}".format(k, v) for k, v in sorted(nb_pcss_merged.items(), key=lambda t: t[0])))
	print(f"\tMerged list before: {pcs_lst_cur_len} items; Merged list after: {len(pcs_lst_cur)} items.")
	return pcs_lst_cur

def consistencyCheck(pcs_lst):
	totPCSs     = len(pcs_lst)
	nb_sample   = 1000
	sample_size = min(100,totPCSs)
	for idx in range(nb_sample):
		begIdx    = random.randint(0, totPCSs-sample_size)
		pcs_pairs = itertools.combinations(pcs_lst[begIdx:begIdx+sample_size], 2)
		for idx, (pcs1, pcs2) in enumerate(pcs_pairs):
			if isOverlap(pcs1,pcs2):
				print(f"ERROR! While merging PCSs of all species, it is expected that two PCSs will never overlap. However, two PCSs were found to be overlapping: {pcs1} and {pcs2}.")
				sys.exit()

def mergePCS_all(qChrom, my_dataset):
	timeTrack = Time()
	timeTrack.start()
	
	pcs_lst_merged = []
	mergedPcsFilename = my_dataset.getTmpFilename_computeWindows(qChrom)
	if os.path.isfile(mergedPcsFilename):
		print(f"Loading merged PCSs from {mergedPcsFilename}...")
		pcs_lst_merged = pickle.load(open(mergedPcsFilename, 'rb'))
	else:	
		speciesUCSCnames = my_dataset.speciesUCSCnames
		for ucscName_other in speciesUCSCnames:
			# Load PCSs.
			timeTrack.startStep(f"[{ucscName_other}] Load PCSs")
			pcsFilename = my_dataset.getOutFilename_extractPCS(ucscName_other,qChrom)
			pcs_lst     = pickle.load(open(pcsFilename, 'rb'))
			timeTrack.stopStep()
			# Merge PCSs.
			timeTrack.startStep(f"[{ucscName_other}] Merge PCSs")
			pcs_lst_merged = mergePCS_pairwise(pcs_lst_merged,pcs_lst)
			timeTrack.stopStep()
			# Check consistency (make sure PCS positions are non-overlapping).
			timeTrack.startStep(f"[{ucscName_other}] Check overlap")
			consistencyCheck(pcs_lst_merged)
			timeTrack.stopStep()

		print(f"Save temporary file with merged PCSs in {mergedPcsFilename}.")
		with open(mergedPcsFilename, 'wb') as pickleFile:
			pickle.dump(pcs_lst_merged, pickleFile, protocol=pickle.HIGHEST_PROTOCOL)

	timeTrack.stop()
	timeTrack.print()
	return pcs_lst_merged


####################################
# "Window" related code

# TODO: Maybe instead of returning a list of tuples,
#       it is better to return a list of numbers,
#       where i and i+1 define the windows coordinates?
def computeWindows(pcs_lst, windowSize):
	""" This method computes coordinates for windows 
    in the human/reference genome. This method 
    ensures that the window coordinates are defined 
    in such a way as to prevent the disruption of 
    any PCS found  in the species included in the dataset.
	
	:returns: A list of window coordinates, each approximately equal 
		to the specified ``windowSize`` (either exactly that size or 
		slightly larger). The coordinates are formatted as ``[begPos, endPos)``, 
		where ``endPos`` is excluded from the interval.
	
	"""

	# Make sure PCSs are sorted.
	pcs_lst = sorted(pcs_lst, key=lambda pcs: pcs.qPosBeg)
	win_lst = []
	pcs_sizes_stats=defaultdict(int)
	for p in pcs_lst: pcs_sizes_stats[p.size]+=1
	print(f"\tStats on merged PCSs (N={len(pcs_lst)})")
	print("\n".join("\t{}\t{}".format(k, v) for k, v in sorted(pcs_sizes_stats.items(), key=lambda t: t[0])))
	# Initialize window position.
	curPcsIdx = 0
	curPcs    = pcs_lst[curPcsIdx]
	winEndPos = curPcs.qPosBeg
	windowSize_stats=defaultdict(int)
	while (curPcsIdx <= len(pcs_lst)):

		winBegPos = winEndPos
		winEndPos = winBegPos+windowSize
		# Check if window coordinates conflict with a PCS.
		checkPCSs = True
		while checkPCSs:
			# Case 1: PCS is fully encompassed in the window. Check next PCS.
			qPosEnd = curPcs.qPosBeg + curPcs.size
			if (qPosEnd < winEndPos):
				curPcsIdx += 1
			else:
				# Case 2.1: PCS is after the interval. Stop looking at PCSs.
				if(curPcs.qPosBeg >= winEndPos):
					checkPCSs = False
				# Case 2.2: PCS conflicts with window. Solve conflict.
				else:
					winEndPos = qPosEnd
					curPcsIdx += 1
			if(curPcsIdx < len(pcs_lst)):
				curPcs    = pcs_lst[curPcsIdx]
			else:
				checkPCSs = False
		# Add window to the list.
		win_lst.append((winBegPos,winEndPos))
		windowSize_stats[winEndPos-winBegPos] += 1
	print(f"\tStats on windows")
	print("\n".join("\t{}\t{}".format(k, v) for k, v in sorted(windowSize_stats.items(), key=lambda t: t[0])))
	return win_lst

####################################
# "PCS distribution" related code

def distribPCS_single(win_lst, pcs_lst):
	""" This method computes the PCS size distribution for each window given in a list.

	:param win_lst: List of positions representing non-overlapping consecutive windows. 
		Each position is a tuple in the format ``(begPos, endPos)``, where ``begPos`` and ``endPos`` 
		are positive integers indicating absolute positions on the chromosome. 
		The list is sorted in ascending order by these positions.

	:type desc: list of tuples of integers
	:param pcs_lst: list of PCSs, containing their positions in the chromosome. List should be sorted by position.
	:type pcs_lst: list of named tuples ``Pcs``

	:returns: a dictionary that maps tuples of window coordinates ``(begPos, endPos)``
		on the chromosome to another dictionary. This inner dictionary contains the 
		distribution of PCS sizes, with PCS sizes as keys and their respective counts as values.
	
	"""
	
	# Make sure PCSs are sorted.
	pcs_lst = sorted(pcs_lst, key=lambda pcs: pcs.qPosBeg)

	pcs_distrib_allwin = {}
	curPcsIdx = 0
	curPcs    = pcs_lst[curPcsIdx]
	for windowId in win_lst:
		winBegPos, winEndPos = windowId
		# Find all PCSs encompassed by the window.
		winPCSs_lst = []
		checkPCSs   = (curPcsIdx < len(pcs_lst))
		while checkPCSs:
			# Case 1: This condition should always be true.
			if(curPcs.qPosBeg >= winBegPos):
				qPosEnd = curPcs.qPosBeg + curPcs.size
				# Case 1.1: PCS is inside the window.
				if(qPosEnd <= winEndPos): # The last position is not included neither in the window nor in the PCS.
					winPCSs_lst.append(curPcs.size)
					# Check next PCS.
					curPcsIdx += 1
					if(curPcsIdx < len(pcs_lst)):
						curPcs    = pcs_lst[curPcsIdx]
					else:
						checkPCSs = False
				# Case 1.2: PCS appears after window.
				elif(curPcs.qPosBeg >= winEndPos):
					checkPCSs = False
				# Case 1.3: Unexpected: PCS overlaps window.
				else:
					print(f"ERROR! Unexpected scenario: PCS and window positions partially overlap.")
					sys.exit()
			# Case 2: Unexpected: PCS without window.
			else:
				print(f"ERROR! Unexpected scenario: PCS without window. {winBegPos=} {printPCS(curPcs)}")
				sys.exit()
		# Add PCS size distribution.
		pcs_distrib_allwin[windowId] = Counter(winPCSs_lst)
	return pcs_distrib_allwin

def distribPCS_all(qChrom, win_lst, my_dataset):
	""" This method computes the PCS size distribution for every window in a chromosome, 
	considering all species present in the dataset.
	"""
	pcs_distrib_all  = {}
	speciesUCSCnames = my_dataset.speciesUCSCnames
	for ucscName_other in speciesUCSCnames:
		print(f"\t[{ucscName_other}] Computing PCS size distribution")
		pcsFilename = my_dataset.getOutFilename_extractPCS(ucscName_other,qChrom)
		pcs_lst     = pickle.load(open(pcsFilename, 'rb'))

		# Computes the PCS distribution of each window of a species/chromosome.
		pcs_distrib_allwin = distribPCS_single(win_lst, pcs_lst)
		
		# Saves PCS size distribution of the species/chromosome.
		pcs_distrib_all[ucscName_other] = pcs_distrib_allwin
	return pcs_distrib_all

####################################
# "Check parameters" related code

def checkInputFiles(qChrom, my_dataset):
	speciesMissing   = []
	speciesUCSCnames = my_dataset.speciesUCSCnames
	for ucscName_other in speciesUCSCnames:
		pcsFilename = my_dataset.getOutFilename_extractPCS(ucscName_other,qChrom)
		if (not os.path.isfile(pcsFilename)): speciesMissing.append(ucscName_other)
	if (len(speciesMissing) > 0):
		print(f" WARNING! PCS file of {qChrom} not found for {len(speciesMissing)} out of {len(speciesUCSCnames)} species. These species will be ignored: {', '.join(speciesMissing)}")

####################################
# MAIN.
####################################
if (__name__ == '__main__'):

	parser = argparse.ArgumentParser(description="Compute windows and their PCS size distributions.")
	parser.add_argument("-chr", help="Chromosome from the reference species for which the PCS size distributions will be calculated.", type=str, required=True)
	
	args         = parser.parse_args()
	my_dataset   = Dataset()

	prefixQuery	 = my_dataset.refsp_ucscName # In our study, query is always the human genome (hg38)
	qChrom		 = args.chr
	windowSize   = my_dataset.windowSize
	minPCSsize   = my_dataset.minPCSsize

	print("***********************************************************")
	print("*                    Compute windows                      *")
	print("* Compute window coordinates in the reference genome of a *")
	print("* chromosome. It also computes the PCS size distribution  *") 
	print("* of each window for each genome in the dataset.          *")
	print("***********************************************************")
	print(f"- Reference genome: {prefixQuery}")
	print(f"- Chromosome: {qChrom}")
	print(f"- Window size: ~{windowSize} base pairs.")
	print(f"- PCS directory (input): {my_dataset.dirPCSs}")
	print(f"- Window directory (output): {my_dataset.dirWindows}")
	print(f"- Temporary directory (output): {my_dataset.dirTemp}")
	print(f"- Minimum size for PCS to be considered: >={minPCSsize} base pairs.")

	checkInputFiles(qChrom, my_dataset)

	# Save times for each step and overall time.
	timeTrack = Time()
	timeTrack.start()

	# Read PCSs from each species, and merge them into a single list of PCSs.
	timeTrack.startStep("Merging PCSs")
	pcs_lst_merged = mergePCS_all(qChrom, my_dataset)
	timeTrack.stopStep()

	# Compute windows based on a PCS list.
	timeTrack.startStep("Computing windows")
	win_lst = computeWindows(pcs_lst_merged,windowSize)
	timeTrack.stopStep()

	# Compute PCS size distribution for every window in a chromosome, 
	# considering all species present in the dataset.
	timeTrack.startStep("Computing PCS size distribution")
	pcs_distrib_all = distribPCS_all(qChrom, win_lst, my_dataset)
	timeTrack.stopStep()

	# Save PCS size distribution (pickle format).
	pcs_distrib_filename = my_dataset.getOutFilename_computeWindows(qChrom)
	with open(pcs_distrib_filename, 'wb') as pickleFile:
		pickle.dump(pcs_distrib_all, pickleFile, protocol=pickle.HIGHEST_PROTOCOL)
		print("Pickle saved!")

	timeTrack.stop()
	timeTrack.print()
	print("Done!")
