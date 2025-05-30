"""  Functions for loading and handling data produced during the pipeline.
"""
import os
import sys
import pickle

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils.dataset import Dataset
 
def splitList(input_list, m):
	# Calculate the size of each chunk
	avg_chunk_size = len(input_list) // m
	remainder = len(input_list) % m
	# Create the sublists using a list comprehension
	return [input_list[i*avg_chunk_size + min(i,remainder):(i+1)*avg_chunk_size + min(i+1,remainder)] for i in range(m)]

def readWindows(my_dataset, ucscName, qChrom):
	# Check if input file exists.
	windows_filename = my_dataset.getOutFilename_computeWindows(qChrom)
	if (not os.path.isfile(windows_filename)):
		print(f"ERROR! File not found: {windows_filename}")
		sys.exit()
	# Load windows.
	pcs_distrib_all = pickle.load(open(windows_filename, 'rb'))
	if(ucscName not in pcs_distrib_all):
		print(f"ERROR! Species not found: {ucscName} ({qChrom}).")
		sys.exit()
	return pcs_distrib_all[ucscName]

def readEvolutionaryTimes(my_dataset, prefixTarget, qChrom, alpha):
	evoltimes_filename = my_dataset.getOutFilename_estimateEvolTimes(prefixTarget, qChrom, alpha)
	if (not os.path.isfile(evoltimes_filename)):
		print(f"ERROR! File not found: {evoltimes_filename}")
		sys.exit()
	# Load estimated evolutionary times.
	winIds_ests, winIds_refs, winSizeRefs_ts = pickle.load(open(evoltimes_filename, 'rb'))
	return (winIds_ests, winIds_refs, winSizeRefs_ts)
