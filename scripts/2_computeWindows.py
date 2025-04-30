import os
import pickle
import collections
import glob
import sys
import re
import numpy as np

from utils.basicTypes import Pcs, Time


verbose=False

####################################
# "Window" related code

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
		if(nb_pcs_merged == 0):
			if (idx_to_add < len(pcs_lst_cur)):
				pcs_lst_cur.insert(idx_to_add, p_new)
			else:
				pcs_lst_cur.append(p_new)
		else:
			pcs_lst_merge = pcs_lst_cur[idx_to_add:(idx_to_add+nb_pcs_merged)] + [p_new]

			qPosBeg = p_new.qPosBeg
			qPosEnd = p_new.qPosBeg + p_new.size
			pcs_largest = p_new

			[p.qPosBeg for p in pcs_lst_merge]

			pcs_merged = Pcs(sizePCS, p_new.tChrom, p_new.tStrand, p_new.tPosBeg, p_new.qChrom, p_new.qStrand, qPosBeg)









def mergePCS(qChrom, speciesUCSCnames, ucscName_human, dirPCSs):
	pcs_lst_merged = []
	for ucscName_other in speciesUCSCnames:
		pcsFilename = os.path.join(dirPCSs, f"{ucscName_human}.{ucscName_other}.{qChrom}.pcs.pickle")
		pcs_lst     = pickle.load(open(pcsFilename, 'rb'))
		pcs_lst_merged = mergePCS_pairwise(pcs_lst_merged,pcs_lst)
		# Check consistency (make sure PCS positions are non-overlapping).
	return pcs_lst_merged

	


	
	minPos, maxPos, minSize, maxSize = getMinMaxPosSize(allPCSregions)
	# Matrix: (PCS size, (binned) position in the chromosome)
	Z = np.zeros((1, maxPos-minPos+1)) #, dtype=np.uint8
	for prefixTarget, PCSregions in allPCSregions.items():
		#print(prefixTarget)
		for region in PCSregions:
			size = region[1]-region[0]
			for pos in range(region[0],region[1]):
				Z[0][pos-minPos] += 1
	return Z[0]

def computeWindows(Z, windowSize):
	Z_wdw = np.copy(Z)
	Z_wdw[Z_wdw == 0] = -1
	totSize   = len(Z)
	
	windowLastPos = 0		
	while (windowLastPos < totSize):

		# Update region.
		if(Z_wdw[windowLastPos] < 0): Z_wdw[windowLastPos] = 0

		# Find next position.
		windowNextPos = min(windowLastPos + windowSize, totSize)
		while((windowNextPos < totSize) and (Z[windowNextPos] > 0)):
			windowNextPos += 1
		windowLastPos = windowNextPos
	return Z_wdw

# 0 : position without PCS (gaps)
# >0: position with at least 1 PCS
# <0: position without PCS, but converted into "useful" region
def countRegions(Z):
	indicesNoPCS  = np.where(Z == 0)[0]
	
	prevIdx = 0
	qty0s   = 1
	
	regions = []
	
	for idxNoPCS in indicesNoPCS:
		intervalSize = idxNoPCS-prevIdx
		
		# Non consecutive Gap positions
		if(intervalSize > 1):

			begRegion  = prevIdx+1
			endRegion  = idxNoPCS
			gapL	   = qty0s

			# Compute PCS distribution in the region.
			region = Z[begRegion:endRegion]
			indicesPCS = np.where(region > 0)[0]
			PCSdistrib = []
			prevIdxPCS = -1
			pcsSize	= 0
			for idxPCS in indicesPCS:					
				intervalSizePCS = idxPCS-prevIdxPCS
				if((prevIdxPCS >= 0) and (intervalSizePCS > 1)):
					if(region[prevIdxPCS] > 0): 
						PCSdistrib.append(pcsSize)
					pcsSize = 1
				else:
					pcsSize += 1
				prevIdxPCS=idxPCS
			# Last PCS.
			if (pcsSize > 1): PCSdistrib.append(pcsSize)
			if(sum(PCSdistrib) != len(indicesPCS)):
				print(f"ERROR! Diff amount of PCS= {sum(PCSdistrib)} {len(indicesPCS)}")
				
			# Add info about the region	
			regions.append((begRegion, endRegion, gapL, PCSdistrib))
			
			qty0s   = 1
			
		# Consecutive 1s
		elif(intervalSize == 1):
			qty0s += 1
		
		# Update previous idx
		prevIdx=idxNoPCS

	# Last region
	# 0001222200 -> info about the region was already added
	# 000122220  -> info about the region was already added 
	# 0000112112112
	if (Z[-1] != 0):
		begRegion  = (indicesNoPCS[-1])+1
		endRegion  = len(Z)
		gapL	   = qty0s
		
		# Compute PCS distribution in the region.
		region = Z[begRegion:endRegion]
		indicesPCS = np.where(region > 0)[0]
		PCSdistrib = []
		prevIdxPCS = -1
		pcsSize	= 0
		for idxPCS in indicesPCS:					
			intervalSizePCS = idxPCS-prevIdxPCS
			if((prevIdxPCS >= 0) and (intervalSizePCS > 1)):
				if(region[prevIdxPCS] > 0): 
					PCSdistrib.append(pcsSize)
				pcsSize = 1
			else:
				pcsSize += 1
			prevIdxPCS=idxPCS
		# Last PCS.
		if (pcsSize > 1): PCSdistrib.append(pcsSize)
		if(sum(PCSdistrib) != len(indicesPCS)):
			print(f"ERROR! Diff amount of PCS= {sum(PCSdistrib)} {len(indicesPCS)}")
				
		regions.append((begRegion, endRegion, gapL, PCSdistrib))
		
	return regions

def countPCSs(regions,nbPCSlb=10):
	cntNbPCS	 = []
	goodRegions  = []
	emptyRegions = []
	for idx, region in enumerate(regions):
		begRegion, endRegion, gapL, PCSdistrib = region
		size  = endRegion-begRegion
		nbPCS = len(PCSdistrib)
		if (nbPCS == 0):
			#print(f"ERROR! Empty PCS distribution! b:{begRegion} e:{endRegion}")
			emptyRegions.append(region)
		else:
			cntNbPCS.append(nbPCS)
			if (nbPCS > nbPCSlb):
				goodRegions.append(region)
	return collections.Counter(cntNbPCS), goodRegions, emptyRegions

def findRegion(pos, regions):
	for idx, region in enumerate(regions):
		begRegion, endRegion, gapL, PCSdistrib = region
		if ((begRegion <= pos) and (pos <= endRegion)):
			return idx
	return -1

def findRegionsLongPCSs(allPCSregions, regions, longPCSlb=500):
	selRegions = []
	minPos, maxPos, minSize, maxSize = getMinMaxPosSize(allPCSregions)
	for prefixTarget, PCSregions in allPCSregions.items():
		#print(f"\t{prefixTarget}")
		for region in PCSregions:
			begPos = region[0]-minPos
			endPos = region[1]-minPos
			size   = endPos-begPos
			if (size >= longPCSlb):
				selRegionIdx = findRegion(begPos, regions)
				if (selRegionIdx >= 0): selRegions.append(selRegionIdx)
	selRegionsUnique = list(set(selRegions))
	return [regions[idx] for idx in selRegionsUnique]

def findRegionsSharedPCS(Z_cnt, regions, cntlb=10):
	selRegions = []
	posSharedPCS = np.where(Z_cnt > cntlb)[0]
	
	for idx, region in enumerate(regions):
		begRegion, endRegion, gapL, PCSdistrib = region
		posSharedPCS = np.where(Z_cnt[begRegion:endRegion] > cntlb)[0]
		if (len(posSharedPCS) > 0):
			selRegions.append(region)
			
	return selRegions

def computePCSperRegion(Z, allPCSregions):

	minPos, maxPos, minSize, maxSize = getMinMaxPosSize(allPCSregions)
	indicesNoPCS  = np.where(Z == 0)[0]
	distribPCSperSpecies = {}
	
	for prefixTarget in allPCSregions.keys():

		#print(prefixTarget)
		distribPCSperSpecies[prefixTarget] = []
		PCSregions = allPCSregions[prefixTarget]

		# Set -1s
		Z_sp = np.empty(maxPos-minPos+1)
		Z_sp.fill(-1)
		
		# Set 1s
		for region in PCSregions:
			begPos_i = region[0]
			endPos_i = region[1]
			size_i   = endPos_i-begPos_i
			Z_sp[(begPos_i-minPos):(endPos_i-minPos)] = 1

		# Set 0s
		for posNoPCS in indicesNoPCS: 
			Z_sp[posNoPCS] = 0
			
		# Compute distribution PCS size per region
		regions = countRegions(Z_sp)
		distribPCSperSpecies[prefixTarget] = regions
		
	return distribPCSperSpecies

def getMetadataSpecies(prefixTargetLst, speciesMetadataFilename):
	timeDiv={}
	commonNames={}
	if (os.path.isfile(speciesMetadataFilename)):
		speciesMetadataTmp=pickle.load(open(speciesMetadataFilename, 'rb'))
		for prefixTarget in prefixTargetLst:
			m = re.match(r'(\D+)\d+\D*', prefixTarget)
			speciesAbbrev = m.group(1) if m else prefixTarget
			for s in speciesMetadataTmp.keys():
				if((speciesMetadataTmp[s]["Abbrev"]).startswith(speciesAbbrev)):
					timeDiv[prefixTarget] = speciesMetadataTmp[s]["DivTime"]
					commonNames[prefixTarget] = (speciesMetadataTmp[s]["CommonName"])
					#print("\tTime="+str(speciesMetadataTmp[s]["DivTime"]))
					break
		if (len(timeDiv.keys()) != len(prefixTargetLst)):
			sys.exit(f"ERROR! Some metadata info could not be retrieved: Found={len(timeDiv.keys())} Expected={len(allPCSregions.keys())}.")	
	else:
		sys.exit(f"ERROR! Metadata file not found in the directory {dirMetadata}.")
	return timeDiv, commonNames

####################################
# MAIN.
####################################
# Use: python3 2-computeWindows.py chr1 1000
if (__name__ == '__main__'):
	isCluster = True if os.path.isdir("/bucket/MillerU/Priscila") else False

	prefixQuery	 = "hg38"
	qChrom		 = sys.argv[1]
	windowSize   = int(sys.argv[2])

	dirMetadata  = "/bucket/MillerU/Priscila/paper/data/inputs" # Directory "data" is available for download at Zenodo.
	dirPCSs      = "/flash/MillerU/Priscila/paper-validation/pcs"
	dirWindows   = "/flash/MillerU/Priscila/paper-validation/windows"

	####################################
	# Load species information.
	speciesMetadataFilename = os.path.join(dirMetadata, "speciesMetadata.pickleProtocol4.pickle")
	if(not os.path.isfile(speciesMetadataFilename)):
		print(f"ERROR! File with metadata for species not found ({speciesMetadataFilename}).")
		sys.exit()
	speciesUCSCnames, divergenceTimes, commonNames = pickle.load(open(speciesMetadataFilename, 'rb'))
	# Sort by divergence time.
	speciesUCSCnames = sorted(speciesUCSCnames, key=lambda species: (divergenceTimes[species], commonNames[species]))


	####################################
	# Read PCSs from each species, and merge them in a single list of PCSs.
	pcs_lst_merged = mergePCS(speciesUCSCnames, dirPCSs)

	for ucscName_other in speciesUCSCnames:


	####################################
	# Read PCSs.

	allPCSs = None
	if (os.path.isfile(coordsFilename)):
		allPCSregions = pickle.load(open(coordsFilename, 'rb'))
	else:
		sys.exit("ERROR! PCS coords file not found.")

	####################################
	# Get divergence time from human (Mya)

	prefixTargetLst = list(allPCSregions.keys())
	speciesMetadataFilename = os.path.join(dirMetadata, "speciesMetadata-v1.pickle")
	timeDiv, commonNames    = getMetadataSpecies(prefixTargetLst, speciesMetadataFilename)

	speciesSortedByDivTime = [(prefixTarget, timeDiv[prefixTarget], commonNames[prefixTarget]) for prefixTarget in prefixTargetLst]
	speciesSortedByDivTime = sorted(speciesSortedByDivTime, key=lambda tup: (tup[1],tup[2]))

	####################################
	# Compute windows.

	minPos, maxPos, minSize, maxSize = getMinMaxPosSize(allPCSregions)

	Z_cnt = mergePCS(allPCSregions)
	Z_wdw = computeWindows(Z_cnt, windowSize)
	regions_wdw = countRegions(Z_wdw)
	nbWindows   = len(regions_wdw)

	distribPCSperSpecies = computePCSperRegion(Z_wdw, allPCSregions)

	####################################
	# Save windows.

	for idxSpecies, (prefixTarget, divTime, commonName) in enumerate(speciesSortedByDivTime):
		print(prefixTarget)
		regions = distribPCSperSpecies[prefixTarget]
		if (len(regions) != nbWindows):
			print(f"ERROR! {prefixTarget} has diff. nb. of windows: found={len(regions)} expected={nbWindows}")
			
		else:
			winFilename = prefixQuery + "." + qChrom + "." + prefixTarget + "." + str(windowSize) + ".win"
			with open(os.path.join(dirOut, winFilename), 'w') as winFile:
				# Write headers.
				winFile.write(f"#Human ({prefixQuery}), {qChrom}; {commonName} ({prefixTarget}); Divergence Time: {divTime}; Window Size: {windowSize}; Window Count: {nbWindows} \n")
				winFile.write(f"#WindowID WindowStartPos WindowEndPos WindowSize PCSsumSize PCScount PCSdistribution\n")
				# Write content.		   
				for idx, region in enumerate(regions):
					begPos, endPos, gapL, PCSdistrib = region
					PCSdistrib = sorted(PCSdistrib)
					PCSdistribStr = ",".join([str(pcs) for pcs in PCSdistrib])
					winFile.write(f"{idx} {begPos+minPos} {endPos+minPos} {endPos-begPos} {sum(PCSdistrib)} {len(PCSdistrib)} {PCSdistribStr}\n")
				
	####################################
	# Save pickle.
	with open(os.path.join(dirOut, prefixQuery + "." + qChrom + "." + str(windowSize) + '.windows.pickle'), 'wb') as pickleFile:
		pickle.dump(distribPCSperSpecies, pickleFile, protocol=pickle.HIGHEST_PROTOCOL)
		print("Pickle saved!")

	print("Done!")
