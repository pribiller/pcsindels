import os
import re
import sys
import subprocess
import pickle
import argparse

from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from utils.dataset import Dataset
from utils.basicTypes import CompRes

def createSlurmScript(	slurmFilename, memory, cores, duration, 
						dirCode, outFilename, logFilename, tmpFilename,
						chrom, overwriteFiles):

	python_module = "python/3.11.4"
	python_env    = "~/python-env-2024"

	jobName = f"{chrom}.pcs"
		
	with open(slurmFilename, "w") as slurm_file:
		codeLines  = f"module load {python_module}\nsource {python_env}/bin/activate\nstart_time=$(date +%s)\n\n"
		if(overwriteFiles): codeLines += f"rm -f -- {tmpFilename}\n"
		codeLines += f"python3 {dirCode}/2_computeWindows.py -chr {chrom}\n"
		codeLines += "\nend_time=$(date +%s)\ntotal_time=$(( end_time - start_time ))\n\ndeactivate\n"
		codeLines += 'echo "Job ID: $SLURM_JOB_ID"\n'
		codeLines += f"ls -l {outFilename}" 
		codeLines += " | awk '{sum += $5} END {printf \"Total disk usage: %.3f GB\\n\", sum / 1073741824}'\n"
		codeLines += f"ls -l {tmpFilename}" 
		codeLines += " | awk '{sum += $5} END {printf \"Temporary disk usage: %.3f GB\\n\", sum / 1073741824}'\n"
		codeLines += 'echo "\nTotal time to run: $total_time seconds"\nsstat --allsteps -j $SLURM_JOB_ID --format=JobID,MaxRSS\n'
		slurm_file.write(f"#!/bin/bash\n#SBATCH -p compute\n#SBATCH -t {duration}\n#SBATCH --output={logFilename}\n#SBATCH --job-name={jobName}\n#SBATCH --mem={memory}G\n#SBATCH --nodes=1\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task={cores}\n\n{codeLines}\n")
	return slurmFilename

####################################
# MAIN.
####################################
# Usage:   python3 2_computeWindows_runAll.py [--overwrite, optional]
# Example: python3 ~/code/cluster/2_computeWindows_runAll.py --overwrite

if (__name__ == '__main__'):

	parser = argparse.ArgumentParser(description="[Cluster Version] Compute windows.")
	parser.add_argument("--overwrite", action='store_true', help="If this flag is specified, any existing output files will be overwritten during the run.")
	
	args = parser.parse_args()
	overwriteFiles = args.overwrite

	# This data structure contains all needed information about the dataset.
	my_dataset  = Dataset()

	# Directory organization.
	dirCode     = "/home/p/priscila-do/code" 
	dirSlurm    = "/home/p/priscila-do/code/slurm"

	dirPCSs     = my_dataset.dirPCSs
	dirWin      = my_dataset.dirWindows
	dirTemp     = my_dataset.dirTemp

	print("***********************************************************")
	print("*              Compute windows - CLUSTER                  *")
	print("* Compute window coordinates in the reference genome of a *")
	print("* chromosome. It also computes the PCS size distribution  *") 
	print("* of each window for each genome in the dataset.          *")
	print("***********************************************************")
	print(f"- Directory for Python code: {dirCode}")
	print(f"- Directory for Slurm scripts: {dirSlurm}")
	print(f"- Directory for PCSs (input): {dirPCSs}")
	print(f"- Directory for windows (output): {dirWin}")
	print(f"- Directory for temporary files: {dirTemp}")
	print(f"- Overwrite output files? {overwriteFiles}")
	
	# Check if all directories exist.
	for dirPath in [dirCode, dirSlurm, dirPCSs, dirWin, dirTemp]:
		if (not os.path.exists(dirPath)):
			print(f"ERROR! Directory not found ({dirPath}).")
			sys.exit()		

	# Species information.
	ucscName_human   = my_dataset.refsp_ucscName
	chromLst         = my_dataset.chromLst
	speciesUCSCnames = my_dataset.speciesUCSCnames

	# For each chromosome.
	for qChrom in chromLst:

		# Check if input files exist.
		pcsFilenames = [my_dataset.getOutFilename_extractPCS(ucscName_other,qChrom) for ucscName_other in speciesUCSCnames]

		logFilename  = my_dataset.getLogFilename_computeWindows(qChrom)
		outFilename  = my_dataset.getOutFilename_computeWindows(qChrom)
		tmpFilename  = my_dataset.getTmpFilename_computeWindows(qChrom)

		skipRun = ((not overwriteFiles) and os.path.isfile(outFilename))
		if(skipRun): 
			print(f"WARNING! All files for {qChrom} were already computed. Skipping run...")
			continue

		skipRun = (not all(os.path.isfile(filename) for filename in pcsFilenames))
		if (skipRun):
			print(f"ERROR! {qChrom} : Check if PCS files exist for all species.")
			continue

		# Slurm Parameters.
		comp_res = my_dataset.getCompRes_computeWindows(qChrom)
		duration = comp_res.time # hh:mm:ss
		memory	 = comp_res.mem  # GB
		cores    = 1             # number of cores

		# Create slurm script.
		print(f"[{qChrom}] Create slurm script...")
		slurmFilename = os.path.join(dirSlurm,f"{qChrom}.win.slurm")
		createSlurmScript(	slurmFilename, memory, cores, duration, 
							dirCode, outFilename, logFilename, tmpFilename,
							qChrom, overwriteFiles)

		# Run slurm script.
		print(f"[{qChrom}] Running slurm script... {slurmFilename}")
		process = subprocess.Popen(['sbatch', slurmFilename]) 
