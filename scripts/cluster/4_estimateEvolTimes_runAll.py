import os
import re
import sys
import subprocess
import pickle
import argparse
import itertools

from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from utils.dataset import Dataset
from utils.basicTypes import CompRes

def createSlurmScript(	slurmFilename, memory, cores, duration, 
						dirCode, outFilenamePattern, logFilename, 
						ucscName_other, chromLst, alpha, overwriteFiles):

	python_module = "python/3.11.4"
	python_env    = "~/python-env-2024"

	jobName = f"alpha{alpha}.taus"
	
	# python3 4_estimateEvolTimes.py -sp_ucsc_name [UCSC name] -chr [chromosome name] -alpha [any number > 1] -cores [nb. of cores] [--overwrite, optional]
	# evolTimes-ests.mm39.chr16.alpha1.1.winSize3805.pickle

	with open(slurmFilename, "w") as slurm_file:
		codeLines  = f"module load {python_module}\nsource {python_env}/bin/activate\nstart_time=$(date +%s)\n\n"
		for chrom in chromLst:
			codeLines += f"python3 {dirCode}/4_estimateEvolTimes.py -sp_ucsc_name {ucscName_other} -chr {chrom} -alpha {alpha} -cores {cores} {'--overwrite' if (overwriteFiles) else ''}\n"
		codeLines += "\nend_time=$(date +%s)\ntotal_time=$(( end_time - start_time ))\n\ndeactivate\n"
		codeLines += 'echo "Job ID: $SLURM_JOB_ID"\n'
		codeLines += f"ls -l {outFilenamePattern}" 
		codeLines += " | awk '{sum += $5} END {printf \"Total disk usage: %.3f GB\\n\", sum / 1073741824}'\n"
		codeLines += 'echo "\nTotal time to run: $total_time seconds"\nsstat --allsteps -j $SLURM_JOB_ID --format=JobID,MaxRSS\n'
		slurm_file.write(f"#!/bin/bash\n#SBATCH -p compute\n#SBATCH -t {duration}\n#SBATCH --output={logFilename}\n#SBATCH --job-name={jobName}\n#SBATCH --mem={memory}G\n#SBATCH --nodes=1\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task={cores}\n\n{codeLines}\n")
	return slurmFilename

####################################
# MAIN.
####################################
# Usage:   python3 4_estimateEvolTimes_runAll.py --cores [nb. of cores] [--overwrite, optional]
# Example: python3 ~/code/cluster/4_estimateEvolTimes_runAll.py --cores 10 --overwrite
#          python3 ~/code/cluster/4_estimateEvolTimes_runAll.py --cores 30 --overwrite

if (__name__ == '__main__'):

	parser = argparse.ArgumentParser(description="[Cluster Version] Evolutionary Time Inference: Estimation.")
	parser.add_argument("--cores", help="Number of cores to be used during setup.", type=int, required=True)
	parser.add_argument("--overwrite", action='store_true', help="If this flag is specified, any existing output files will be overwritten during the run.")
	
	args = parser.parse_args()
	overwriteFiles = args.overwrite
	nbcores = int(args.cores)

	# This data structure contains all needed information about the dataset.
	my_dataset  = Dataset()

	# Directory organization.
	dirCode  = "/home/p/priscila-do/code" 
	dirSlurm = "/home/p/priscila-do/code/slurm"

	dirWin   = my_dataset.dirWindows
	dirTemp  = my_dataset.dirTemp
	dirSetupEvolTimes = my_dataset.dirSetupEvolTimes
	dirEstEvolTimes   = my_dataset.dirEstEvolTimes

	print("******************************************************************")
	print("*           Estimate evolutionary times - CLUSTER                *")
	print("* Estimate a posterior distribution of evolutionary times for    *")
	print("* each window across all chromosomes. Each window can yield      *")
	print("* various estimates based on the alpha parameter and the species *")
	print("* in the dataset, which affects the PCS size distribution.       *")
	print("* For instance, in the 40-vertebrate dataset, using two alpha    *")
	print("* values (1.1 and 10) results in 80 evolutionary time posteriors *")
	print("* for a single window.                                           *")
	print("******************************************************************")
	print(f"- Number of cores: {nbcores}")
	print(f"- Directory for Python code: {dirCode}")
	print(f"- Directory for Slurm scripts: {dirSlurm}")
	print(f"- Directory for windows (input): {dirWin}")
	print(f"- Directory for setup files (input): {dirSetupEvolTimes}")
	print(f"- Directory for estimates (output): {dirEstEvolTimes}")
	print(f"- Directory for temporary files: {dirTemp}")
	print(f"- Overwrite output files? {overwriteFiles}")

	# Check if all directories exist.
	for dirPath in [dirCode, dirSlurm, dirWin, dirSetupEvolTimes, dirEstEvolTimes, dirTemp]:
		if (not os.path.exists(dirPath)):
			print(f"ERROR! Directory not found ({dirPath}).")
			sys.exit()		

	# Check if input files exist.
	winFilenames = [my_dataset.getOutFilename_computeWindows(qChrom) for qChrom in my_dataset.chromLst]
	abortRun     = (not all(os.path.isfile(filename) for filename in winFilenames))
	if (abortRun):
		print(f"ERROR! Check if window files exist for all chromosomes.")
		sys.exit()

	# Method information.
	alphas    = my_dataset.alphas
	UCSCnames = my_dataset.speciesUCSCnames

	# For each alpha and each species being compared to the reference genome.
	for (alpha, ucscName_other) in itertools.product(alphas,UCSCnames):

		# Slurm Parameters.
		comp_res = my_dataset.getCompRes_estimateEvolTimes(ucscName_other)
		duration = comp_res.time # hh:mm:ss
		memory	 = comp_res.mem  # GB
		
		logFilename = my_dataset.getLogFilename_estimateEvolTimes(ucscName_other,alpha)
		outFilenamePattern = my_dataset.getOutFilenamePattern_estimateEvolTimes(ucscName_other,alpha)

		# Create slurm script.
		print(f"[α={alpha}, {ucscName_other}] Create slurm script...")
		slurmFilename = os.path.join(dirSlurm,f"alpha{alpha}.{ucscName_other}.ests.slurm")
		createSlurmScript(	slurmFilename, memory, nbcores, duration, 
							dirCode, outFilenamePattern, logFilename, 
							ucscName_other, my_dataset.chromLst, alpha, overwriteFiles)
		# Run slurm script.
		print(f"[α={alpha}] Running slurm script... {slurmFilename}")
		process = subprocess.Popen(['sbatch', slurmFilename]) 
