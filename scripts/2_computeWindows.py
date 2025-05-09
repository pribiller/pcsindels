"""Computes the PCS size distribution for every window 
in a chromosome, considering all species present in the dataset.

- **Use**::

	python3 2_computeWindows.py -chr [chromosome name]

- **Example of Usage (vertebrate dataset with 40 species)**::

	python3 2_computeWindows.py -chr chr16

- **Input Parameter (mandatory)**:

:-chr: Chromosome name of the reference species (e.g. chr1, chr2, ..., chrX, chrY).

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
the metadata file given in ``-metadata_file``. In our study,
this file corresponds to ``speciesMetadata.pickleProtocol4.pickle``,
and contains the information (UCSC name, divergence time, etc.) of 
40 vertebrate species, whose data was downloaded from the UCSC website.

Cluster resources
-----------------

In case you want to run this Python script stand-alone::

	srun -p compute -t 2:00:00 --mem 20G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash

Otherwise you can use the script ``../cluster/2_computeWindows_runAll.py``
to run this script for all 24 chromosomes (chr1, chr2, ..., chrX, chrY).

Time
----

To be added.

Storage
-------

To be added.

Function details
----------------

Only relevant functions have been documented below. 
For more details on any function, check the comments in the souce code.

"""
import os
import pickle
from collections import Counter
import sys
import re
import itertools
import argparse

from utils.basicTypes import Pcs, Time
from utils.dataset import Dataset

####################################
# "Merge PCSs" related code

def isOverlap(pcs1,pcs2):
	b1, e1, b2, e2 = (pcs1.qPosBeg, pcs1.qPosBeg+pcs1.size, pcs2.qPosBeg, pcs2.qPosBeg+pcs2.size)
	return (((b1 <= b2) and (e1 >= b2)) or ((b2 <= b1) and (e2 >= b1)))

# Merge PCSs.
def mergePCS_pairwise(pcs_lst_cur,pcs_lst_new):

	# Make sure new PCSs are sorted.
	pcs_lst_new = sorted(pcs_lst_new, key=lambda pcs: pcs.qPosBeg)

	# Check if current list is empty.
	if(not pcs_lst_cur): return pcs_lst_new

	# New PCS can be: 
	# - Discarded, if it is encompassed by another PCS in the list;
	# - Added, if its position does not overlap any PCS position in the list;
	# - Merged, if its position overlaps one or more PCSs in the list.
	idx_to_add = 0
	for p_new in pcs_lst_new:
		nb_pcs_merged = 0

		# Find in which index the new PCS could be included.
		# List is sorted by position.
		searchForPosition = True
		while (searchForPosition):
			p_cur = pcs_lst_cur[idx_to_add]
			if(isOverlap(p_cur,p_new)): nb_pcs_merged += 1
			searchForPosition = (p_new.qPosBeg <= (p_cur.qPosBeg + p_cur.size))
			if(searchForPosition): idx_to_add += 1
			searchForPosition = (idx_to_add < len(pcs_lst_cur))

		# Update list (if needed).
		if(nb_pcs_merged > 0):
			pcs_lst_merge = pcs_lst_cur[idx_to_add:(idx_to_add+nb_pcs_merged)] + [p_new]

			qPosBeg = min([p.qPosBeg for p in pcs_lst_merge])
			qPosEnd = max([p.qPosBeg+p.size for p in pcs_lst_merge])
			sizePCS = qPosEnd-qPosBeg

			pcs_largest = max(pcs_lst_merge, key=lambda pcs: pcs.size)
			pcs_merged  = Pcs(sizePCS, pcs_largest.tChrom, pcs_largest.tStrand, pcs_largest.tPosBeg, pcs_largest.qChrom, pcs_largest.qStrand, qPosBeg)

			del pcs_lst_cur[idx_to_add:(idx_to_add+nb_pcs_merged)]

		# Add pcs.
		p_add = pcs_merged if (nb_pcs_merged > 0) else p_new
		if (idx_to_add < len(pcs_lst_cur)):
			pcs_lst_cur.insert(idx_to_add, p_add)
		else:
			pcs_lst_cur.append(p_add)
	return pcs_lst_cur


def mergePCS_all(qChrom, my_dataset):
	pcs_lst_merged   = []
	speciesUCSCnames = my_dataset.speciesUCSCnames
	for ucscName_other in speciesUCSCnames:
		pcsFilename = my_dataset.getOutFilename_extractPCS(ucscName_other,qChrom)
		pcs_lst     = pickle.load(open(pcsFilename, 'rb'))
		pcs_lst_merged = mergePCS_pairwise(pcs_lst_merged,pcs_lst)
		# Check consistency (make sure PCS positions are non-overlapping).
		for pcs1, pcs2 in itertools.combinations(pcs_lst_merged, 2):
			if isOverlap(pcs1,pcs2):
				print(f"ERROR! While merging PCSs of all species, it is expected that two PCSs will never overlap. However, two PCSs were found to be overlapping after processing {ucscName_other}: {pcs1} and {pcs2}.")
				sys.exit()
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

	# Initialize window position.
	curPcsIdx = 0
	curPcs    = pcs_lst[curPcsIdx]
	winEndPos = curPcs.qPosBeg 
	while (curPcsIdx <= len(pcs_lst)):

		winBegPos = winEndPos
		winEndPos = winBegPos+windowSize

		# Check if window coordinates conflict with a PCS.
		checkPCSs = True
		while checkPCSs:
			# Case 1: PCS is fully encompassed in the window. Check next PCS.
			if (curPcs.endPos < winEndPos):
				curPcsIdx += 1
				curPcs    = pcs_lst[curPcsIdx]
			else:
				# Case 2.1: PCS is after the interval. Stop looking at PCSs.
				if(curPcs.begPos > winEndPos):
					checkPCSs = False
				# Case 2.2: PCS conflicts with window. Solve conflict.
				else:
					winEndPos = curPcs.endPos+1
					curPcsIdx += 1
					curPcs    = pcs_lst[curPcsIdx]
		
		# Add window to the list.
		win_lst.append((winBegPos,winEndPos))

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
	
	pcs_distrib_allwin = {}
	curPcsIdx = 0
	curPcs    = pcs_lst[curPcsIdx]
	for windowId in win_lst:
		winBegPos, winEndPos = windowId
		# Find all PCSs encompassed by the window.
		winPCSs_lst = []
		checkPCSs   = True
		while checkPCSs:
			# Case 1: This condition should always be true.
			if(curPcs.begPos >= winBegPos):
				# Case 1.1: PCS is inside the window.
				if(curPcs.endPos < winEndPos): # The last position (winEndPos) is not encompassed by the window.
					winPCSs_lst.append(curPcs.size)
				# Case 1.2: PCS appears after window.
				elif(curPcs.begPos >= winEndPos):
					checkPCSs = False
				# Case 1.3: Unexpected: PCS overlaps window.
				else:
					print(f"ERROR! Unexpected scenario: PCS and window positions partially overlap.")
					sys.exit()
			# Case 2: Unexpected: PCS without window.
			else:
				print(f"ERROR! Unexpected scenario: PCS without window.")
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
		pcsFilename = my_dataset.getOutFilename_extractPCS(ucscName_other,qChrom)
		pcs_lst     = pickle.load(open(pcsFilename, 'rb'))

		# Computes the PCS distribution of each window of a species/chromosome.
		pcs_distrib_allwin = distribPCS_single(win_lst, pcs_lst)
		
		# Saves PCS size distribution of the species/chromosome.
		pcs_distrib_all[ucscName_other] = pcs_distrib_allwin
	return pcs_distrib_all

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
